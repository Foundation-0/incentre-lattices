#!/usr/bin/env python3
"""Second method for the rigidity counts of results/rigidity-nd, n = 2..8: arithmetic modulo a prime.

Standard library only. No floating point. Imports nothing from this repository, and nothing from rigidity_exact.py,
so the two scripts share no code.

How it differs from rigidity_exact.py:
  * rigidity_exact.py divides the irrational factor c = 2 - sqrt2 out of every row and works over Q(i).
  * This script keeps the TRUE coordinates p(i,s) = s c e_i with c = 2 - sqrt2 and the true cross bar vector
    (2 - 2c) e_i, and does all arithmetic in the field GF(p). For a prime p = 1 (mod 8) both sqrt2 and i exist in GF(p),
    so c is an ordinary number there. The Bloch phases z_i are random nonzero elements of GF(p) (fixed seed).
  * It also checks that the intra and cross bars have equal length, which the other script does not: 2 c^2 and
    (2 - 2c)^2 must agree in GF(p).

What a rank mod p shows: the rank of an integer-valued matrix mod p is at most its rank over the rationals, so a zero-mode
count mod p is at least the true count. Three different primes agreeing with the formulas is evidence, not proof.
Evidence label: VERIFIED (finite check).

Usage:
  python3 -I -O -B modp_check.py            run n = 2..8 for three primes, exit 0 only if every count and control passes
  python3 -I -O -B modp_check.py --corrupt  expect one zero mode too many: must report MISMATCH and exit 1
"""
import argparse
import itertools
import random

PRIMES = (1000000009, 998244353, 2013265921)  # all = 1 (mod 8)


def require(cond, msg):
    """Check that survives python -O (never a bare assert)."""
    if not cond:
        raise SystemExit("MISMATCH: " + msg)


def must_fail(name, accepted, log):
    require(not accepted, f"control '{name}' was wrongly accepted")
    log.append(name)


def sqrt_mod(a, p):
    """Square root of a modulo an odd prime p (Tonelli-Shanks); the caller checks the result."""
    a %= p
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1:
            tt = tt * tt % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


def rank_mod(rows, p):
    m = [r[:] for r in rows]
    rk, nc = 0, len(m[0])
    for col in range(nc):
        piv = next((r for r in range(rk, len(m)) if m[r][col] % p), None)
        if piv is None:
            continue
        m[rk], m[piv] = m[piv], m[rk]
        iv = pow(m[rk][col], p - 2, p)
        m[rk] = [x * iv % p for x in m[rk]]
        for r in range(len(m)):
            if r != rk and m[r][col] % p:
                f = m[r][col]
                m[r] = [(a - f * b) % p for a, b in zip(m[r], m[rk])]
        rk += 1
        if rk == len(m):
            break
    return rk


def joints(n):
    return [(i, s) for i in range(n) for s in (1, -1)]


def system(n, p, c, z=None, drop=None, off_axis=False):
    """Rows of the rigidity matrix with true coordinates (positions s c e_i), mod p. Intra bars always; cross bars with phases
    z when z is given. drop: index of one intra bar to delete. off_axis: move joint (0,+1) by c along axis 1."""
    js = joints(n)
    idx = {v: k for k, v in enumerate(js)}
    pos = {(i, s): [(s * c) % p if k == i else 0 for k in range(n)] for i, s in js}
    if off_axis:
        pos[(0, 1)][1] = c % p
    bars = [(a, b) for a, b in itertools.combinations(js, 2) if a[0] != b[0]]
    if drop is not None:
        bars = [b for k, b in enumerate(bars) if k != drop]
    rows, sq = [], set()
    for a, b in bars:
        d = [(pos[b][k] - pos[a][k]) % p for k in range(n)]
        sq.add(sum(x * x for x in d) % p)
        r = [0] * (n * len(js))
        for k in range(n):
            r[n * idx[a] + k] = (-d[k]) % p
            r[n * idx[b] + k] = d[k]
        rows.append(r)
    if z is not None:
        for i in range(n):
            a, b = (i, 1), (i, -1)
            # b in the next cell sits at pos[b] + 2 e_i; the bar vector from a is (2 - 2c) e_i
            d = [(pos[b][k] + (2 if k == i else 0) - pos[a][k]) % p for k in range(n)]
            sq.add(sum(x * x for x in d) % p)
            r = [0] * (n * len(js))
            for k in range(n):
                r[n * idx[a] + k] = (-d[k]) % p
                r[n * idx[b] + k] = d[k] * z[i] % p
            rows.append(r)
    return rows, sq


def counts(rows, n, p):
    nc = 2 * n * n
    rk = rank_mod(rows, p)
    return (nc - rk, len(rows) - rk)  # (zero modes, self-stresses)


def expected(n):
    """Formulas stated in results/rigidity-nd/PROOF.md (n >= 3) and the n = 2 values it reports."""
    if n >= 3:
        return {"cp": (n * (n + 1) // 2, n * (n - 3) // 2), "q0": (n * (n + 1) // 2, n * (n - 1) // 2),
                "zm1": (n * (n - 1) // 2, n * (n - 3) // 2), "gen": (n * (n - 1) // 2, n * (n - 3) // 2)}
    return {"cp": (4, 0), "q0": (3, 1), "zm1": (2, 0), "gen": (2, 0)}


def run(corrupt=False):
    out, log = {}, []
    for p in PRIMES:
        require(p % 8 == 1, f"prime {p} is not 1 mod 8")
        s2, ii = sqrt_mod(2, p), sqrt_mod(p - 1, p)
        require(s2 * s2 % p == 2 and ii * ii % p == p - 1, f"square roots mod {p} are wrong")
        c = (2 - s2) % p
        require(2 * c * c % p == (2 - 2 * c) ** 2 % p, f"mod {p}: intra and cross bar lengths differ")
        rnd = random.Random(p)
        res_p = {}
        for n in range(2, 9):
            zg = [rnd.randrange(2, p - 1) for _ in range(n)]
            exp = expected(n)
            got = {}
            for key, z in (("cp", None), ("q0", [1] * n), ("zm1", [p - 1] * n), ("gen", zg)):
                rows, sq = system(n, p, c, z)
                require(len(sq) == 1, f"mod {p} n={n} {key}: bar lengths not all equal")
                got[key] = counts(rows, n, p)
            planes = []
            for k in range(n + 1):
                rows, _ = system(n, p, c, [1] * k + zg[k:])
                planes.append(counts(rows, n, p))
            for key in exp:
                want = exp[key]
                if corrupt and key == "cp" and n >= 3:
                    want = (want[0] + 1, want[1])
                require(got[key] == want, f"mod {p} n={n} {key}: got {got[key]}, expected {want}")
            if n >= 3:
                require(planes == [(n * (n - 1) // 2 + k, n * (n - 3) // 2 + k) for k in range(n + 1)],
                        f"mod {p} n={n}: planes {planes}")
            res_p[str(n)] = {"counts": {k: list(v) for k, v in got.items()}, "planes": [list(x) for x in planes]}
        out[str(p)] = res_p

        # controls, run for each prime at n = 3 and n = 4
        for n in (3, 4):
            exp = expected(n)
            base_rows, _ = system(n, p, c)
            base = counts(base_rows, n, p)
            drop_rows, _ = system(n, p, c, drop=1)
            after = counts(drop_rows, n, p)
            must_fail(f"mod {p}: counts unchanged after deleting one intra bar (n={n})", after == base, log)
            wrong_c = 1  # not 2 - sqrt2: bars of different length
            wrong_rows, wrong_sq = system(n, p, wrong_c, [1] * n)
            must_fail(f"mod {p}: bar lengths still equal with c = 1 (n={n})", len(wrong_sq) == 1, log)
            off_rows, off_sq = system(n, p, c, off_axis=True)
            must_fail(f"mod {p}: intra bar lengths still equal after moving a joint (n={n})", len(off_sq) == 1, log)
            rows, _ = system(n, p, c, [rnd.randrange(2, p - 1) for _ in range(n)])
            must_fail(f"mod {p}: generic counts equal the q = 0 counts (n={n})",
                      counts(rows, n, p) == exp["q0"], log)
        # n = 2 (the square) is not rigid: claiming the n >= 3 count for it must be rejected
        rows, _ = system(2, p, c)
        must_fail(f"mod {p}: the square (n=2) has the rigid count n(n+1)/2 = 3 zero modes",
                  counts(rows, 2, p)[0] == 3, log)
    return out, log


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--corrupt", action="store_true",
                    help="expect one zero mode too many for the cross-polytope: must report MISMATCH and exit 1")
    args = ap.parse_args()
    out, log = run(args.corrupt)
    if args.corrupt:
        raise SystemExit("MISMATCH: --corrupt run reached the end without a mismatch")
    print(f"modp_check: OK (n = 2..8; primes {', '.join(map(str, PRIMES))}; {len(log)} controls rejected)")
    for p, res in out.items():
        print(p, " ".join(f"n={n}:{tuple(v['counts']['gen'])}" for n, v in res.items()))


if __name__ == "__main__":
    main()
