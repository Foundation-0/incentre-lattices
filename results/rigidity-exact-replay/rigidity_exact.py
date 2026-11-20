#!/usr/bin/env python3
"""Exact-arithmetic replay of the rigidity counts of results/rigidity-nd, n = 2..8.

Standard library only (fractions). No floating point, no tolerance. Imports nothing from this repository.
Evidence label: VERIFIED (a finite check at the stated n and wavevectors). It is not a proof for all n or all q.

Why exact arithmetic is possible although c = 2 - sqrt2 is irrational:
  every bar vector of the framework is a nonzero multiple of an integer direction
    intra bar  p(i,s) - p(j,t) = c (s e_i - t e_j)         -> direction s e_i - t e_j
    cross bar  (2 - 2c) e_i                                   -> direction e_i
  Scaling a row of a matrix by a nonzero number does not change the rank or the dimension of the left kernel
  (the self-stresses are the row dependencies). So the matrix built from the integer directions has the same rank,
  zero-mode count and self-stress count as the matrix with the true coordinates. The factor c is divided out exactly.

What this does NOT check:
  * that c = 2 - sqrt2 gives every bar the same length (arithmetic in Q(sqrt2); the float run of rigidity_nd.py covers it
    and the elastic check of results/independent-checks in pull request 102 uses Q(sqrt2));
  * elastic constants;
  * any statement for all wavevectors q: the Bloch phase z_i = exp(2 i q_i) is tried at exact points of Q(i) only,
    namely z = 1, z = -1 and Pythagorean unit points (a+bi)/m. Two generic wavevector families are tried (all z_i = -1,
    and one Pythagorean point per axis) plus the planes with k phases equal to 1;
  * counts only: nothing here identifies the zero modes as rotations.

Controls (must fail): each broken input must be rejected by the same code path that accepts the real framework.
  * delete one bar: rank and stress count must change (rigidity is lost at n = 3, a stress is lost at n >= 4);
  * one joint reduced to degree n-1: the framework must stop being rigid;
  * one joint moved off its axis: the exact squared bar lengths stop being all equal. This control tests LENGTHS ONLY.
    The counts do not change under this move, so the counts cannot reject it;
  * wrong tensions: the explicit stress basis with a sign twist must not pass the exact stress check;
  * periodic: delete the cross bar X_n at generic phases: the zero-mode count must change.
Recorded but not a control: a joint moved along its own axis keeps rank and counts (as the theorem says for any
  0 < t < 1) while the bar lengths change.

Usage:
  python3 -I -O -B rigidity_exact.py                run n = 2..8, write exact_counts.json, exit 0 only if every check passes
  python3 -I -O -B rigidity_exact.py --check        recompute and compare with the stored exact_counts.json
  python3 -I -O -B rigidity_exact.py --nmax 5       smaller range
  python3 -I -O -B rigidity_exact.py --corrupt      expect one count too high: must report MISMATCH and exit 1
"""
import argparse
import itertools
import json
import pathlib
import sys
from fractions import Fraction as F

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "exact_counts.json"


def require(cond, msg):
    """Check that survives python -O (never a bare assert)."""
    if not cond:
        raise SystemExit("MISMATCH: " + msg)


def must_fail(name, verdict_accepts, log):
    """A control passes only if the checker REJECTED the broken input."""
    require(not verdict_accepts, f"control '{name}' was wrongly accepted")
    log.append(name)


class QI:
    """Exact element of Q(i)."""
    __slots__ = ("re", "im")

    def __init__(self, re=0, im=0):
        self.re, self.im = F(re), F(im)

    def __add__(s, o): return QI(s.re + o.re, s.im + o.im)
    def __sub__(s, o): return QI(s.re - o.re, s.im - o.im)
    def __mul__(s, o): return QI(s.re * o.re - s.im * o.im, s.re * o.im + s.im * o.re)
    def iszero(s): return s.re == 0 and s.im == 0

    def inv(s):
        n = s.re * s.re + s.im * s.im
        return QI(s.re / n, -s.im / n)


ONE = QI(1)


def rank(rows):
    """Exact rank by Gauss-Jordan elimination over Q(i)."""
    m = [r[:] for r in rows]
    if not m:
        return 0
    nc = len(m[0])
    rk = 0
    for col in range(nc):
        p = next((r for r in range(rk, len(m)) if not m[r][col].iszero()), None)
        if p is None:
            continue
        m[rk], m[p] = m[p], m[rk]
        inv = m[rk][col].inv()
        m[rk] = [x * inv for x in m[rk]]
        for r in range(len(m)):
            if r != rk and not m[r][col].iszero():
                f = m[r][col]
                m[r] = [a - f * b for a, b in zip(m[r], m[rk])]
        rk += 1
        if rk == len(m):
            break
    return rk


def joints(n):
    return [(i, s) for i in range(n) for s in (1, -1)]


def positions(n):
    """Coordinates divided by the common factor c (exact integers): p(i,s) = s e_i."""
    return {(i, s): [s if k == i else 0 for k in range(n)] for i, s in joints(n)}


def intra_bars(n):
    return [((i, s), (j, t)) for (i, s), (j, t) in itertools.combinations(joints(n), 2) if i != j]


def rigidity_rows(n, pos, bars, z=None, cross=False):
    """Rows of the (periodic if cross) rigidity matrix. Row of bar a-b: -d at a, +d at b, with d the integer direction
    (position difference over c for intra bars, e_i times the phase z_i on the far end for cross bars)."""
    js = joints(n)
    idx = {v: k for k, v in enumerate(js)}
    rows = []
    for a, b in bars:
        d = [pos[b][k] - pos[a][k] for k in range(n)]
        r = [QI(0) for _ in range(n * len(js))]
        for k in range(n):
            r[n * idx[a] + k] = QI(-d[k])
            r[n * idx[b] + k] = QI(d[k])
        rows.append(r)
    if cross:
        for i in range(n):
            a, b = (i, 1), (i, -1)
            r = [QI(0) for _ in range(n * len(js))]
            r[n * idx[a] + i] = QI(-1)
            r[n * idx[b] + i] = z[i]
            rows.append(r)
    return rows


def counts(n, pos, bars, z=None, cross=False):
    rows = rigidity_rows(n, pos, bars, z, cross)
    rk = rank(rows)
    ncols = n * len(joints(n))
    return {"bars": len(rows), "rank": rk, "zero_modes": ncols - rk, "self_stresses": len(rows) - rk}


UNIT_POINTS = [(3, 4, 5), (5, 12, 13), (8, 15, 17), (7, 24, 25), (20, 21, 29), (9, 40, 41), (12, 35, 37), (11, 60, 61)]


def generic_z(n):
    out = []
    for a, b, m in UNIT_POINTS[:n]:
        require(a * a + b * b == m * m, f"({a},{b},{m}) is not a Pythagorean triple")
        out.append(QI(F(a, m), F(b, m)))
    return out


def exact_stress_basis_check(n, twist=False):
    """Explicit stresses of PROOF.md (a): tension w_ij on the four bars of the pair (i,j), sum_j w_ij = 0 at every i.
    Exact kernel of the unsigned K_n incidence matrix, then exact check that y^T R = 0 on the intra rows.
    With twist=True the tension on each bar is multiplied by the sign of its first end: a deliberately wrong stress."""
    pairs = list(itertools.combinations(range(n), 2))
    m = [[F(1) if i in p else F(0) for p in pairs] for i in range(n)]
    piv, rk = [], 0
    for col in range(len(pairs)):
        p = next((r for r in range(rk, n) if m[r][col] != 0), None)
        if p is None:
            continue
        m[rk], m[p] = m[p], m[rk]
        iv = 1 / m[rk][col]
        m[rk] = [x * iv for x in m[rk]]
        for r in range(n):
            if r != rk and m[r][col] != 0:
                f = m[r][col]
                m[r] = [a - f * b for a, b in zip(m[r], m[rk])]
        piv.append(col)
        rk += 1
    free = [c for c in range(len(pairs)) if c not in piv]
    basis = []
    for fcol in free:
        w = [F(0)] * len(pairs)
        w[fcol] = F(1)
        for r, pc in enumerate(piv):
            w[pc] = -m[r][fcol]
        basis.append(w)
    pos = positions(n)
    bars = intra_bars(n)
    rows = rigidity_rows(n, pos, bars)
    for w in basis:
        y = [w[pairs.index((a[0], b[0]))] * (a[1] if twist else 1) for a, b in bars]
        for col in range(len(rows[0])):
            tot = QI(0)
            for yr, row in zip(y, rows):
                tot = tot + QI(yr) * row[col]
            if not tot.iszero():
                return len(basis), False
    return len(basis), True


def squared_lengths(pos, bars):
    """Sorted set of exact squared lengths (units of c^2)."""
    return sorted({sum((pos[a][k] - pos[b][k]) ** 2 for k in range(len(pos[a]))) for a, b in bars})


def build_certificate(nmax=8, corrupt=False):
    """Run every count and control for n = 2..nmax and return the certificate dict (raises SystemExit on a mismatch)."""
    require(2 <= nmax <= 8, "--nmax must be between 2 and 8")
    ns = range(2, nmax + 1)

    rows_out = []
    for n in ns:
        pos = positions(n)
        ib = intra_bars(n)
        triv = n * (n + 1) // 2
        o = {"n": n, "dofs_per_cell": n * 2 * n, "bars_per_cell": len(ib) + n, "intra_bars": len(ib)}
        sqlen = squared_lengths(pos, ib)
        require(sqlen == [2], f"n={n}: intra squared lengths {sqlen}")
        o["intra_squared_length_over_c2"] = sqlen
        cp = counts(n, pos, ib)
        o["crosspolytope"] = cp
        o["crosspolytope_nontrivial_motions"] = cp["zero_modes"] - triv
        o["crosspolytope_rigid"] = cp["zero_modes"] == triv
        if n >= 3:
            expect_zero = triv + 1 if corrupt else triv
            require(cp["zero_modes"] == expect_zero, f"n={n}: cross-polytope zero modes {cp['zero_modes']}, expected {expect_zero}")
            require(cp["self_stresses"] == n * (n - 3) // 2, f"n={n}: stress dim {cp['self_stresses']}")
            dim, ok = exact_stress_basis_check(n)
            require(dim == n * (n - 3) // 2 and ok, f"n={n}: explicit stress basis")
            o["explicit_stress_basis_dimension"] = dim
            o["explicit_stress_basis_is_stress"] = ok
        else:
            require(cp["zero_modes"] == triv + 1 and cp["self_stresses"] == 0, f"n=2 {cp}")
            o["note"] = "n = 2: the square has one internal flex (outside the theorem)"
        q0 = counts(n, pos, ib, [ONE] * n, True)
        o["q0"] = q0
        pred_q0 = (triv, n * (n - 1) // 2) if n >= 3 else (3, 1)
        require((q0["zero_modes"], q0["self_stresses"]) == pred_q0, f"n={n}: q0 {q0}")
        gen = []
        for zs in ([QI(-1)] * n, generic_z(n)):
            c = counts(n, pos, ib, zs, True)
            gen.append((c["zero_modes"], c["self_stresses"]))
        pred = (n * (n - 1) // 2, n * (n - 3) // 2) if n >= 3 else (2, 0)
        require(all(g == pred for g in gen), f"n={n}: generic {gen} != {pred}")
        o["generic_z_minus1_and_pythagorean_zero_modes_selfstresses"] = gen
        planes = []
        for k in range(n + 1):
            zs = generic_z(n)
            for i in range(k):
                zs[i] = ONE
            c = counts(n, pos, ib, zs, True)
            planes.append((k, c["zero_modes"], c["self_stresses"]))
            if n >= 3:
                require((c["zero_modes"], c["self_stresses"]) == (pred[0] + k, pred[1] + k), f"n={n} k={k}: {c}")
        o["planes_k_zeromodes_selfstresses"] = planes
        rows_out.append(o)

    log, detail = [], {}
    for n in ns:
        if n < 3:
            continue
        pos = positions(n)
        ib = intra_bars(n)
        triv = n * (n + 1) // 2
        base = counts(n, pos, ib)
        # delete one bar
        c2 = counts(n, pos, ib[1:])
        detail[f"delete_one_bar_n{n}"] = {"before": base, "after": c2}
        must_fail(f"rank and stress count unchanged after deleting one bar (n={n})",
                  (c2["rank"], c2["self_stresses"]) == (base["rank"], base["self_stresses"]), log)
        # joint (0,+) has 2n-2 intra bars; keep n-1 of them
        at = [b for b in ib if (0, 1) in b]
        drop = set(at[n - 1:])
        c3 = counts(n, pos, [b for b in ib if b not in drop])
        detail[f"joint_degree_n-1_n{n}"] = c3
        must_fail(f"still rigid with one joint of degree n-1 (n={n})", c3["zero_modes"] == triv, log)
        # move joint (0,+) off its axis (position in units of c): lengths only
        pos4 = {k: v[:] for k, v in pos.items()}
        pos4[(0, 1)][1] = 1
        sq = squared_lengths(pos4, ib)
        c4 = counts(n, pos4, ib)
        detail[f"move_joint_offaxis_n{n}"] = {"squared_lengths": [str(x) for x in sq], "counts": c4}
        must_fail(f"all intra bars still equal length after moving a joint off axis (n={n}; lengths only)", sq == [2], log)
        # record: move the same joint along its axis; counts must not change, lengths do
        pos5 = {k: v[:] for k, v in pos.items()}
        pos5[(0, 1)][0] = F(1, 3)
        c5 = counts(n, pos5, ib)
        sq5 = squared_lengths(pos5, ib)
        require(c5["rank"] == base["rank"] and c5["self_stresses"] == base["self_stresses"],
                f"n={n}: axis move changed the counts")
        require(sq5 != [2], f"n={n}: axis move did not change the lengths")
        detail[f"move_joint_along_axis_n{n}"] = {"squared_lengths": [str(x) for x in sq5], "counts": c5}
        # a wrong set of tensions must not pass the explicit stress check
        # (for n = 3 the basis is empty, so there is nothing to twist)
        if n >= 4:
            _, twisted_ok = exact_stress_basis_check(n, twist=True)
            must_fail(f"sign-twisted tensions accepted as a self-stress (n={n})", twisted_ok, log)
        # periodic: delete the last cross bar at generic phases
        rows_p = rigidity_rows(n, pos, ib, generic_z(n), True)
        full = len(rows_p[0]) - rank(rows_p)
        less = len(rows_p[0]) - rank(rows_p[:-1])
        must_fail(f"generic-q zero-mode count unchanged after deleting cross bar X_n (n={n})", less == full, log)

    cert = {
        "result": "rigidity-exact-replay",
        "evidence_label": "VERIFIED (finite check): n = 2..%d at chosen exact points; not a proof for all n or all q" % nmax,
        "claims_replayed": {
            "cross_polytope": "infinitesimally rigid for n >= 3 with n(n-3)/2 self-stresses (n = 2: one internal flex)",
            "generic_z": "n(n-1)/2 zero modes, n(n-3)/2 self-stresses per cell (two generic phase families)",
            "q0": "n(n+1)/2 zero modes, n(n-1)/2 self-stresses (n = 2: 3 and 1)",
            "planes": "k phases equal to 1 add k to each count",
        },
        "arithmetic": "exact Q(i) with Fraction; the common factor c is divided out of every row; no floats, no tolerance",
        "not_covered": ["equal bar length from c = 2 - sqrt2", "elastic constants", "all wavevectors", "identification of zero modes as rotations"],
        "counts_by_n": rows_out,
        "must_fail_controls_rejected": log,
        "control_details": detail,
    }
    return cert


def dump(cert):
    return json.dumps(cert, indent=2, sort_keys=True, default=str) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--nmax", type=int, default=8, help="largest n (default 8; at most 8 because of the Pythagorean points)")
    ap.add_argument("--corrupt", action="store_true",
                    help="expect n(n+1)/2 + 1 zero modes for the rigid cross-polytope: must report MISMATCH and exit 1")
    ap.add_argument("--check", action="store_true",
                    help="recompute and compare with the stored exact_counts.json instead of writing it (exit 1 on any difference)")
    args = ap.parse_args()
    cert = build_certificate(args.nmax, args.corrupt)
    if args.corrupt:
        raise SystemExit("MISMATCH: --corrupt run reached the end without a mismatch")
    n_ctl = len(cert["must_fail_controls_rejected"])
    if args.check:
        require(args.nmax == 8, "--check compares the full n = 2..8 file")
        require(OUT.read_text() == dump(cert), "recomputed counts differ from the stored exact_counts.json")
        print(f"rigidity-exact-replay: stored file matches (n = 2..{args.nmax}; {n_ctl} controls rejected)")
        return
    OUT.write_text(dump(cert))
    print(f"rigidity-exact-replay: OK (n = 2..{args.nmax}; {n_ctl} controls rejected)")


if __name__ == "__main__":
    main()
