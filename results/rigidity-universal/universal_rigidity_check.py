#!/usr/bin/env python3
"""Exact check of the universal-rigidity note (standard library only, Fractions).

Parts
  A  Maxwell's identity on random frameworks: sum(omega |l|^2) = 0 and trace(Omega) = 2 sum(omega).
  B  The finite cross-polytope: stress space dimension, explicit stress family, Omega = -(W kron J2),
     every nonzero stress is indefinite.
  C  The simplex map: sending vertex k to e_k keeps every bar length of an equal-bar framework.
  D  The periodic incentre framework with the lattice fixed: equilibrium, S = T, Omega(0) positive
     semidefinite with kernel the constants, the sum-of-squares form, the fixed-lattice energy identity.
  E  Controls that must fail, each run through the same check functions as the real inputs.

Usage: python3 universal_rigidity_check.py            run everything (exit 0 only if all pass)
       python3 universal_rigidity_check.py --corrupt  run part D with a wrong bridging tension (must exit 1)
No bare assert is used: every check raises CheckFailed, so -O does not change the result.
"""
import itertools
import random
import sys
from fractions import Fraction as F

SEED = 20261121
CROSS_N = range(2, 8)          # finite cross-polytope dimensions
PERIODIC_N = range(2, 7)       # periodic framework dimensions
RATIONAL_C = [F(1, 2), F(1, 3), F(2, 5), F(3, 4)]


class CheckFailed(Exception):
    pass


def require(cond, msg):
    if not cond:
        raise CheckFailed(msg)


# ---------------------------------------------------------------- exact linear algebra

def rref(rows, ncols):
    m = [list(map(F, r)) for r in rows]
    piv, r = [], 0
    for c in range(ncols):
        p = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if p is None:
            continue
        m[r], m[p] = m[p], m[r]
        inv = 1 / m[r][c]
        m[r] = [x * inv for x in m[r]]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[r])]
        piv.append(c)
        r += 1
        if r == len(m):
            break
    return m[:r], piv


def rank(rows, ncols):
    return len(rref(rows, ncols)[1]) if rows else 0


def nullspace(rows, ncols):
    """Basis of {x : rows x = 0}."""
    red, piv = rref(rows, ncols) if rows else ([], [])
    free = [c for c in range(ncols) if c not in piv]
    basis = []
    for f in free:
        x = [F(0)] * ncols
        x[f] = F(1)
        for r, c in enumerate(piv):
            x[c] = -red[r][f]
        basis.append(x)
    return basis


def psd_rank(M):
    """(is_psd, rank) of a symmetric rational matrix by diagonal-pivot elimination. Exact."""
    n = len(M)
    a = [list(map(F, r)) for r in M]
    rk = 0
    for k in range(n):
        if a[k][k] < 0:
            return False, rk
        if a[k][k] == 0:
            if any(a[k][j] != 0 for j in range(n)):
                return False, rk
            continue
        rk += 1
        for i in range(k + 1, n):
            if a[i][k] != 0:
                f = a[i][k] / a[k][k]
                for j in range(k, n):
                    a[i][j] -= f * a[k][j]
        # keep symmetry exact: clear row k to the right
        for j in range(k + 1, n):
            a[k][j] = F(0)
    return True, rk


def quad(M, x):
    return sum(M[i][j] * x[i] * x[j] for i in range(len(x)) for j in range(len(x)))


def rand_frac(rng, lo=-9, hi=9):
    return F(rng.randint(lo, hi), rng.randint(1, 5))


# ---------------------------------------------------------------- finite frameworks (parts A to C)

def sqdist(p, q):
    return sum((a - b) ** 2 for a, b in zip(p, q))


def stress_matrix(V, bars, omega):
    Om = [[F(0)] * V for _ in range(V)]
    for (a, b), w in zip(bars, omega):
        Om[a][a] += w
        Om[b][b] += w
        Om[a][b] -= w
        Om[b][a] -= w
    return Om


def equilibrium_defects(points, bars, omega):
    """Vertices where sum_j omega_ij (p_i - p_j) is not zero."""
    d = len(points[0])
    res = [[F(0)] * d for _ in points]
    for (a, b), w in zip(bars, omega):
        for k in range(d):
            diff = points[a][k] - points[b][k]
            res[a][k] += w * diff
            res[b][k] -= w * diff
    return [i for i, r in enumerate(res) if any(x != 0 for x in r)]


def stress_space(points, bars):
    d, V = len(points[0]), len(points)
    rows = []
    for i in range(V):
        for k in range(d):
            row = [F(0)] * len(bars)
            for e, (a, b) in enumerate(bars):
                if a == i:
                    row[e] += points[a][k] - points[b][k]
                if b == i:
                    row[e] += points[b][k] - points[a][k]
            rows.append(row)
    return nullspace(rows, len(bars))


def maxwell_defects(points, bars, omega):
    """Problems with Maxwell's identity for an equilibrium stress; empty list means all hold."""
    out = []
    if equilibrium_defects(points, bars, omega):
        out.append("not an equilibrium")
    s = sum(w * sqdist(points[a], points[b]) for (a, b), w in zip(bars, omega))
    if s != 0:
        out.append("sum omega |l|^2 = %s" % s)
    Om = stress_matrix(len(points), bars, omega)
    if sum(Om[i][i] for i in range(len(points))) != 2 * sum(omega):
        out.append("trace != 2 sum omega")
    return out


def random_framework(rng, d):
    V = d + rng.randint(2, 4)
    while True:
        pts = [tuple(F(rng.randint(-6, 6)) for _ in range(d)) for _ in range(V)]
        if len(set(pts)) == V:
            break
    pairs = list(itertools.combinations(range(V), 2))
    bars = [pr for pr in pairs if rng.random() < 0.8]
    return pts, bars


def check_maxwell(rng):
    total = 0
    for d in (1, 2, 3, 4):
        for _ in range(8):
            pts, bars = random_framework(rng, d)
            basis = stress_space(pts, bars)
            combos = list(basis)
            if basis:
                coef = [rng.randint(-3, 3) for _ in basis]
                combos.append([sum(k * b[e] for k, b in zip(coef, basis)) for e in range(len(bars))])
            for om in combos:
                bad = maxwell_defects(pts, bars, om)
                require(not bad, "Maxwell identity fails: %s" % bad)
                total += 1
    require(total > 50, "too few stresses sampled (%d)" % total)
    return total


def cross_polytope(n):
    pts = []
    for i in range(n):
        for s in (1, -1):
            pts.append(tuple(F(s if k == i else 0) for k in range(n)))
    bars = [(a, b) for a, b in itertools.combinations(range(2 * n), 2) if a // 2 != b // 2]
    return pts, bars


def axis_pair_weights(n, w):
    """Stress on all bars of the cross-polytope: w[(i,j)] on the four bars between axes i and j."""
    pts, bars = cross_polytope(n)
    return [F(w[(a // 2, b // 2)]) for a, b in bars]


def zero_row_sum_weights(n):
    """Basis of {w on pairs i<j : sum_j w_ij = 0 for every i}, as dicts."""
    pairs = list(itertools.combinations(range(n), 2))
    rows = [[F(1) if i in pr else F(0) for pr in pairs] for i in range(n)]
    return [dict(zip(pairs, v)) for v in nullspace(rows, len(pairs))]


def indefinite_witnesses(n, w):
    """Exact values (x^T Omega x, y^T Omega y) at two vectors built from one nonzero pair weight w_ij.
    They are -8 w_ij and +8 w_ij, so they have opposite signs."""
    pts, bars = cross_polytope(n)
    om = axis_pair_weights(n, w)
    Om = stress_matrix(2 * n, bars, om)
    pair = next((pr for pr in sorted(w) if w[pr] != 0), None)
    if pair is None:
        return None
    i, j = pair
    x = [F(0)] * (2 * n)
    y = [F(0)] * (2 * n)
    for s in range(2):
        x[2 * i + s] = F(1)
        x[2 * j + s] = F(1)
        y[2 * i + s] = F(1)
        y[2 * j + s] = F(-1)
    return quad(Om, x), quad(Om, y)


def check_cross_polytope(rng, n):
    pts, bars = cross_polytope(n)
    require(all(sqdist(pts[a], pts[b]) == 2 for a, b in bars), "cross-polytope bars not equal")
    require(len(bars) == 2 * n * (n - 1), "bar count")
    dim = len(stress_space(pts, bars))
    expect = max(0, n * (n - 3) // 2)
    require(dim == expect, "n=%d stress dimension %d, expected %d" % (n, dim, expect))
    fam = zero_row_sum_weights(n)
    require(len(fam) == expect, "n=%d family size %d" % (n, len(fam)))
    stresses = [axis_pair_weights(n, w) for w in fam]
    for om in stresses:
        require(not maxwell_defects(pts, bars, om), "family member fails Maxwell or equilibrium")
    require(rank(stresses, len(bars)) == len(stresses), "family dependent")
    tested = 0
    for _ in range(4 if fam else 0):
        coef = [rng.randint(-3, 3) for _ in fam]
        w = {pr: sum(k * f[pr] for k, f in zip(coef, fam)) for pr in fam[0]}
        if not any(w.values()):
            continue
        om = axis_pair_weights(n, w)
        Om = stress_matrix(2 * n, bars, om)
        # Omega = -(W kron J2), W symmetric with zero diagonal
        for a in range(2 * n):
            for b in range(2 * n):
                want = F(0) if a // 2 == b // 2 else -F(w[tuple(sorted((a // 2, b // 2)))])
                require(Om[a][b] == want, "Omega is not -(W kron J2)")
        require(sum(Om[i][i] for i in range(2 * n)) == 0, "trace not zero")
        tr2 = sum(Om[i][j] * Om[j][i] for i in range(2 * n) for j in range(2 * n))
        require(tr2 == 8 * sum(v * v for v in w.values()), "trace Omega^2 != 8 sum w^2")
        lo, hi = indefinite_witnesses(n, w)
        require(lo * hi < 0, "no indefinite witnesses")
        ok, _ = psd_rank(Om)
        require(not ok, "a nonzero cross-polytope stress passed the PSD test")
        tested += 1
    return dim, tested


def simplex_image(V):
    return [tuple(F(1) if k == v else F(0) for k in range(V)) for v in range(V)]


def simplex_map_defects(points, bars, image, length_sq):
    """Problems with the claim: image keeps every bar at squared length length_sq and is not a congruence."""
    out = []
    for a, b in bars:
        if sqdist(points[a], points[b]) != length_sq:
            out.append("input bar not of the stated length")
            break
    for a, b in bars:
        if sqdist(image[a], image[b]) != length_sq:
            out.append("image changes a bar length")
            break
    changed = any(sqdist(points[a], points[b]) != sqdist(image[a], image[b])
                  for a, b in itertools.combinations(range(len(points)), 2))
    if not changed:
        out.append("image is a congruence")
    return out


def check_simplex(n):
    pts, bars = cross_polytope(n)
    bad = simplex_map_defects(pts, bars, simplex_image(2 * n), 2)
    require(not bad, "n=%d simplex map: %s" % (n, bad))
    # the antipodal pair goes from squared distance 4 to 2
    return sqdist(pts[0], pts[1]), sqdist(simplex_image(2 * n)[0], simplex_image(2 * n)[1])


# ---------------------------------------------------------------- periodic incentre framework (part D)

def incentre(n, c):
    """Dots (i,s) at s c e_i, index 2i (s=+) and 2i+1 (s=-). Bars (tail, head, shift), lattice 2 Z^n."""
    pts = [tuple(F(s) * c if k == i else F(0) for k in range(n)) for i in range(n) for s in (1, -1)]
    bars = []
    for a, b in itertools.combinations(range(2 * n), 2):
        if a // 2 != b // 2:
            bars.append((a, b, (0,) * n))
    for i in range(n):
        bars.append((2 * i, 2 * i + 1, tuple(2 if k == i else 0 for k in range(n))))
    return pts, bars


def bar_vec(pts, bar):
    a, b, L = bar
    return [pts[b][k] + L[k] - pts[a][k] for k in range(len(pts[0]))]


def periodic_stress(n, c, w=None, tau_shift=F(0)):
    """Tension w_ij on the four bars of axes i, j and c r_i/(1-c) (+ tau_shift) on cross bar i."""
    pts, bars = incentre(n, c)
    if w is None:
        w = {pr: F(1) for pr in itertools.combinations(range(n), 2)}
    r = [sum(w[tuple(sorted((i, j)))] for j in range(n) if j != i) for i in range(n)]
    om = []
    for k, (a, b, L) in enumerate(bars):
        if any(L):
            i = a // 2
            om.append(c * r[i] / (1 - c) + tau_shift)
        else:
            om.append(F(w[tuple(sorted((a // 2, b // 2)))]))
    return om


def periodic_equilibrium_defects(pts, bars, om):
    d = len(pts[0])
    res = [[F(0)] * d for _ in pts]
    for (a, b, L), w in zip(bars, om):
        l = bar_vec(pts, (a, b, L))
        for k in range(d):
            res[a][k] += w * l[k]
            res[b][k] -= w * l[k]
    return [i for i, r in enumerate(res) if any(x != 0 for x in r)]


def omega0(pts, bars, om):
    return stress_matrix(len(pts), [(a, b) for a, b, _ in bars], om)


def virial_defects(pts, bars, om, c, n, r):
    """S = T, S = diag(4 c r_i) and the lattice trace term; r_i is the row sum of the intra weights."""
    out = []
    S = [[F(0)] * n for _ in range(n)]
    T = [[F(0)] * n for _ in range(n)]
    for bar, w in zip(bars, om):
        l = bar_vec(pts, bar)
        L = bar[2]
        for i in range(n):
            for j in range(n):
                S[i][j] += w * l[i] * l[j]
                T[i][j] += w * l[i] * L[j]
    if S != T:
        out.append("S != T")
    want = [[4 * c * r[i] if i == j else F(0) for j in range(n)] for i in range(n)]
    if S != want:
        out.append("S != diag(4 c r_i)")
    if sum(T[i][i] for i in range(n)) != 4 * c * sum(r):
        out.append("lattice trace term != 8 c W (W = sum of w_ij)")
    return out


def sos_matrix(n, c, w):
    """Matrix of 4 [ m^T L_w m + sum r_i d_i^2/(1-c) ] in the 2n dot variables."""
    r = [sum(w[tuple(sorted((i, j)))] for j in range(n) if j != i) for i in range(n)]
    Mm = [[F(0)] * (2 * n) for _ in range(n)]
    Dd = [[F(0)] * (2 * n) for _ in range(n)]
    for i in range(n):
        Mm[i][2 * i] = Mm[i][2 * i + 1] = F(1, 2)
        Dd[i][2 * i], Dd[i][2 * i + 1] = F(1, 2), F(-1, 2)
    Lw = [[F(0)] * n for _ in range(n)]
    for (i, j), v in w.items():
        Lw[i][i] += v
        Lw[j][j] += v
        Lw[i][j] -= v
        Lw[j][i] -= v
    out = [[F(0)] * (2 * n) for _ in range(2 * n)]
    for a in range(2 * n):
        for b in range(2 * n):
            s = sum(Mm[i][a] * Lw[i][j] * Mm[j][b] for i in range(n) for j in range(n))
            s += sum(Dd[i][a] * r[i] / (1 - c) * Dd[i][b] for i in range(n))
            out[a][b] = 4 * s
    return out


def energy(pts, bars, om, disp, D):
    """sum omega |l'|^2 for dots moved by disp (each in R^D), lattice 2Z^n embedded in R^D."""
    total = F(0)
    n = len(pts[0])
    for (a, b, L), w in zip(bars, om):
        l = [pts[b][k] + L[k] - pts[a][k] if k < n else F(0) for k in range(D)]
        l = [l[k] + disp[b][k] - disp[a][k] for k in range(D)]
        total += w * sum(x * x for x in l)
    return total


def energy_identity_defects(rng, n, c, pts, bars, om, Om):
    D = n + 2
    out = []
    for _ in range(3):
        disp = [[rand_frac(rng) for _ in range(D)] for _ in pts]
        e0 = energy(pts, bars, om, [[F(0)] * D for _ in pts], D)
        e1 = energy(pts, bars, om, disp, D)
        quadsum = sum(quad(Om, [disp[v][k] for v in range(len(pts))]) for k in range(D))
        if e1 - e0 != quadsum:
            out.append("energy identity fails")
            break
    return out


def kernel_defects(Om, expect_dim):
    """Omega PSD, kernel dimension expect_dim, constants in the kernel."""
    out = []
    V = len(Om)
    ok, rk = psd_rank(Om)
    if not ok:
        out.append("not positive semidefinite")
    elif V - rk != expect_dim:
        out.append("kernel dimension %d, expected %d" % (V - rk, expect_dim))
    if any(sum(row) != 0 for row in Om):
        out.append("constants not in kernel")
    return out


def periodic_defects(rng, n, c, om):
    """All part D checks on one stress; returns the list of problems."""
    pts, bars = incentre(n, c)
    out = list(map(lambda s: "equilibrium at vertex %d" % s, periodic_equilibrium_defects(pts, bars, om)))
    Om = omega0(pts, bars, om)
    out += kernel_defects(Om, 1)
    return out


def check_periodic(rng, n, c, tau_shift=F(0)):
    pts, bars = incentre(n, c)
    require(len(bars) == 2 * n * (n - 1) + n, "bar count")
    w1 = {pr: F(1) for pr in itertools.combinations(range(n), 2)}
    om = periodic_stress(n, c, tau_shift=tau_shift)
    bad = periodic_defects(rng, n, c, om)
    require(not bad, "n=%d c=%s: %s" % (n, c, bad))
    r1 = [F(n - 1)] * n
    bad = virial_defects(pts, bars, om, c, n, r1)
    require(not bad, "n=%d c=%s virial: %s" % (n, c, bad))
    # the same identities for a random intra weight w_ij (any real numbers balance at q = 0 with cross tension c r_i/(1-c))
    wr = {pr: F(rng.randint(-3, 3)) for pr in itertools.combinations(range(n), 2)}
    rr = [sum(wr[tuple(sorted((i, j)))] for j in range(n) if j != i) for i in range(n)]
    omr = periodic_stress(n, c, w=wr)
    require(not periodic_equilibrium_defects(pts, bars, omr), "n=%d c=%s: random w not in equilibrium" % (n, c))
    bad = virial_defects(pts, bars, omr, c, n, rr)
    require(not bad, "n=%d c=%s virial for random w: %s" % (n, c, bad))
    Om = omega0(pts, bars, om)
    require(Om == sos_matrix(n, c, w1), "n=%d c=%s: sum-of-squares matrix differs from Omega(0)" % (n, c))
    bad = energy_identity_defects(rng, n, c, pts, bars, om, Om)
    require(not bad, "n=%d c=%s: %s" % (n, c, bad))
    # zero modes at q = 0: periodic displacement u, bar stretch l . (u_head - u_tail)
    rows = []
    for bar in bars:
        l = bar_vec(pts, bar)
        row = [F(0)] * (n * len(pts))
        for k in range(n):
            row[n * bar[0] + k] -= l[k]
            row[n * bar[1] + k] += l[k]
        rows.append(row)
    nullity = n * len(pts) - rank(rows, n * len(pts))
    require(nullity == n * (n + 1) // 2, "n=%d zero modes at q=0: %d" % (n, nullity))
    return nullity - n  # nontrivial modes beyond the n translations


# ---------------------------------------------------------------- controls (part E)

def must_reject(name, problems):
    if not problems:
        raise CheckFailed("control '%s' was accepted" % name)
    print("  control rejected: %-58s (%s)" % (name, problems[0]))


def run_controls(rng):
    # 1. Wrong bridging tension: the same equilibrium and kernel check must reject it.
    for n, c in ((3, F(1, 2)), (4, F(1, 3))):
        pts, bars = incentre(n, c)
        for shift in (F(1), F(-1, 3)):
            om = periodic_stress(n, c, tau_shift=shift)
            must_reject("wrong bridging tension n=%d c=%s shift %s" % (n, c, shift), periodic_defects(rng, n, c, om))
            must_reject("wrong bridging tension breaks S = T, n=%d c=%s shift %s" % (n, c, shift),
                        virial_defects(pts, bars, om, c, n, [F(n - 1)] * n))
        om = periodic_stress(n, c)
        om[0] = -om[0]
        must_reject("one intra bar sign flipped n=%d c=%s" % (n, c), periodic_defects(rng, n, c, om))
        neg = [-x for x in periodic_stress(n, c)]
        must_reject("whole stress negated n=%d c=%s" % (n, c), periodic_defects(rng, n, c, neg))
    # 2. A stress in equilibrium that is not positive semidefinite (one weight made very negative).
    n, c = 3, F(1, 2)
    w = {(0, 1): F(1), (0, 2): F(1), (1, 2): F(-5)}
    om = periodic_stress(n, c, w=w)
    pts, bars = incentre(n, c)
    require(not periodic_equilibrium_defects(pts, bars, om), "control setup: equilibrium expected")
    must_reject("equilibrium stress with a negative weight n=3", periodic_defects(rng, n, c, om))
    # 3. Kernel too big: weights on a disconnected graph give a two-dimensional kernel.
    n, c = 4, F(1, 2)
    w = {pr: F(0) for pr in itertools.combinations(range(n), 2)}
    w[(0, 1)] = F(1)
    w[(2, 3)] = F(1)
    om = periodic_stress(n, c, w=w)
    pts, bars = incentre(n, c)
    require(not periodic_equilibrium_defects(pts, bars, om), "control setup: equilibrium expected")
    must_reject("disconnected weights: kernel larger than constants n=4", periodic_defects(rng, n, c, om))
    # 4. Corrupt Omega: one symmetric entry changed must break the sum-of-squares identity and the energy identity.
    n, c = 3, F(1, 3)
    pts, bars = incentre(n, c)
    om = periodic_stress(n, c)
    Om = omega0(pts, bars, om)
    bad = [row[:] for row in Om]
    bad[0][2] += 1
    bad[2][0] += 1
    w1 = {pr: F(1) for pr in itertools.combinations(range(n), 2)}
    must_reject("corrupt Omega vs sum-of-squares matrix", [] if bad == sos_matrix(n, c, w1) else ["differs"])
    must_reject("corrupt Omega vs energy identity", energy_identity_defects(rng, n, c, pts, bars, om, bad))
    # 5. Maxwell identity detector on a stress that is not in equilibrium.
    pts, bars = random_framework(rng, 2)
    while not stress_space(pts, bars):
        pts, bars = random_framework(rng, 2)
    om = stress_space(pts, bars)[0]
    k = next(e for e, x in enumerate(om) if x != 0)
    om2 = list(om)
    om2[k] *= 2
    must_reject("Maxwell identity on a non-equilibrium stress", maxwell_defects(pts, bars, om2))
    must_reject("Maxwell sum alone on a non-equilibrium stress",
                [p for p in maxwell_defects(pts, bars, om2) if p.startswith("sum omega")])
    # 6. The indefiniteness detector must not flag a PSD matrix (the periodic Omega(0) and the square below).
    n, c = 3, F(1, 2)
    pts, bars = incentre(n, c)
    ok, _ = psd_rank(omega0(pts, bars, periodic_stress(n, c)))
    require(ok, "the PSD test rejected the periodic Omega(0)")
    print("  control accepted as required: the PSD test says PSD for the periodic Omega(0)")
    # 7. Simplex map: a corrupted image changes a bar length, the identity image is a congruence,
    #    and V = d + 1 (regular tetrahedron) admits no non-congruent image.
    n = 4
    pts, bars = cross_polytope(n)
    img = simplex_image(2 * n)
    img[0] = tuple(x + (1 if k == 0 else 0) for k, x in enumerate(img[0]))
    must_reject("simplex image with one vertex moved", simplex_map_defects(pts, bars, img, 2))
    must_reject("identity image is a congruence", simplex_map_defects(pts, bars, pts, 2))
    tet = [tuple(map(F, p)) for p in ((0, 0, 0), (1, 1, 0), (1, 0, 1), (0, 1, 1))]
    tbars = list(itertools.combinations(range(4), 2))
    must_reject("regular tetrahedron, V = d + 1: simplex map is a congruence",
                simplex_map_defects(tet, tbars, simplex_image(4), 2))
    # 8. Unequal bars: a nonzero positive semidefinite stress exists. The unit square with both diagonals.
    sq = [tuple(map(F, p)) for p in ((0, 0), (1, 0), (1, 1), (0, 1))]
    sbars = [(0, 1), (1, 2), (2, 3), (0, 3), (0, 2), (1, 3)]
    space = stress_space(sq, sbars)
    require(len(space) == 1, "square with diagonals: stress dimension %d" % len(space))
    om = space[0]
    if om[0] < 0:
        om = [-x for x in om]
    require(om == [1, 1, 1, 1, -1, -1], "square stress %s" % om)
    require(not maxwell_defects(sq, sbars, om), "square stress fails Maxwell")
    Om = stress_matrix(4, sbars, om)
    ok, rk = psd_rank(Om)
    require(ok and rk == 1, "square stress should be PSD of rank 1")
    require(sum(Om[i][i] for i in range(4)) == 2 * sum(om) == 4, "square trace")
    lens = {sqdist(sq[a], sq[b]) for a, b in sbars}
    require(lens == {1, 2}, "square bar lengths should be unequal")
    must_reject("simplex map on unequal bars (square with diagonals)", simplex_map_defects(sq, sbars, simplex_image(4), 1))
    print("  unequal bars: square with diagonals has the PSD stress %s of rank %d (the equal-bar obstruction does not apply)"
          % (list(map(str, om)), rk))


# ---------------------------------------------------------------- driver

def main(argv):
    corrupt = "--corrupt" in argv
    rng = random.Random(SEED)
    print("seed %d" % SEED)
    if corrupt:
        print("part D with a wrong bridging tension (must fail)")
        for n in (3,):
            check_periodic(rng, n, F(1, 2), tau_shift=F(1))
        print("UNEXPECTED: corrupted input accepted")
        return 0
    cnt = check_maxwell(rng)
    print("A  Maxwell identity: %d stresses on random frameworks (d = 1..4), no violation" % cnt)
    print("B  finite cross-polytope")
    for n in CROSS_N:
        dim, tested = check_cross_polytope(rng, n)
        note = "%d random stresses, all indefinite" % tested if dim else "no nonzero stress to test"
        print("   n=%d stress dimension %d (= max(0, n(n-3)/2)), %s" % (n, dim, note))
    print("C  simplex map, vertex k to e_k in R^(2n)")
    for n in CROSS_N:
        a, b = check_simplex(n)
        print("   n=%d every bar kept, antipodal squared distance %s -> %s" % (n, a, b))
    print("D  periodic incentre framework, lattice fixed")
    for n in PERIODIC_N:
        for c in RATIONAL_C:
            extra = check_periodic(rng, n, c)
        print("   n=%d c in %s: equilibrium, S = T, SOS = Omega(0), PSD, kernel = constants, energy identity; "
              "%d nontrivial zero modes at q = 0" % (n, [str(c) for c in RATIONAL_C], extra))
    print("E  controls")
    run_controls(rng)
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1:]))
    except CheckFailed as exc:
        print("CHECK FAILED: %s" % exc)
        sys.exit(1)
