"""Finite counter-rotation ('tilt') mechanism of the ball/bar contact framework.
Clusters (squares in 2D, octahedra in 3D) rotate rigidly, checkerboard sign; cross bars keep length L=2-2t;
lattice contracts. Find the range before two non-bonded balls (radius R = L/2 at t*) overlap.
Also verify bar lengths stay exact (finite mechanism) and report the macroscopic strain (Poisson ratio)."""
import numpy as np
import itertools
from scipy.spatial.transform import Rotation as Rot

T = 2 - np.sqrt(2); L = 2 - 2 * T; RB = L / 2


def rot2(th):
    c, s = np.cos(th), np.sin(th)
    return np.array([[c, -s], [s, c]])


def build_state(n, Rmat_even, t=T):
    """clusters on a 2x..x2 checkerboard; even clusters rotated by Rmat_even, odd by its inverse.
    Cross bar: from +t e_i of cluster c to -t e_i of cluster c+e_i. With R_odd=R_even^T,
    the bar vector d_i = a_i - t (R + R^T) e_i is identical for every cluster; choose a_i (lattice step)
    = s_i + L * u_i with s_i = t (R+R^T) e_i, u_i unit chosen to keep a_i along the axis direction if
    possible (orthogonal lattice), i.e. a_i = alpha_i e_i + perpendicular part fixed by s_i."""
    R = Rmat_even
    S = t * (R + R.T)
    steps = []
    for i in range(n):
        s = S[:, i]
        perp = s.copy(); perp[i] = 0
        # bar vector b = a_i - s must have length L. Choose a_i = alpha e_i  => b = alpha e_i - s
        need = L ** 2 - perp @ perp
        if need < 0:
            return None
        alpha = s[i] + np.sqrt(need)
        a = np.zeros(n); a[i] = alpha
        steps.append(a)
    return np.array(steps).T, R


def positions(n, steps, R, t=T, ncell=3):
    pts = []; bonded = set()
    cells = list(itertools.product(range(-ncell, ncell + 1), repeat=n))
    index = {}
    for c in cells:
        par = sum(c) % 2
        Rc = R if par == 0 else R.T
        centre = steps @ np.array(c, float)
        for i in range(n):
            for s in (1, -1):
                v = np.zeros(n); v[i] = s * t
                index[(c, i, s)] = len(pts)
                pts.append(centre + Rc @ v)
    pts = np.array(pts)
    # bonds: within cluster (non-antipodal) and cross bars
    for c in cells:
        for (i, s), (j, s2) in itertools.combinations([(i, s) for i in range(n) for s in (1, -1)], 2):
            if i != j:
                bonded.add(frozenset((index[(c, i, s)], index[(c, j, s2)])))
        for i in range(n):
            c2 = list(c); c2[i] += 1; c2 = tuple(c2)
            if c2 in dict.fromkeys(cells):
                bonded.add(frozenset((index[(c, i, 1)], index[(c2, i, -1)])))
    return pts, bonded, index, cells


def check(n, R, ncell=2):
    st = build_state(n, R)
    if st is None:
        return None
    steps, R = st
    pts, bonded, index, cells = positions(n, steps, R, ncell=ncell)
    # bar length check
    blen = [np.linalg.norm(pts[a] - pts[b]) for a, b in map(tuple, bonded)]
    # min non-bonded distance among balls of the central cell's cluster vs all
    centre_ids = [index[((0,) * n, i, s)] for i in range(n) for s in (1, -1)]
    dmin = np.inf
    for a in centre_ids:
        d = np.linalg.norm(pts - pts[a], axis=1)
        for b in np.argsort(d)[1:40]:
            if frozenset((a, b)) not in bonded:
                dmin = min(dmin, d[b]); break
    strain = np.diag(steps) / 2 - 1
    return np.array(blen), dmin, strain, steps


if __name__ == "__main__":
    print("2D counter-rotating squares: theta(deg), bar-length error, min non-bonded distance / 2R, strain xx, yy")
    for deg in [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]:
        bl, dmin, eps, _ = check(2, rot2(np.radians(deg)))
        print(deg, f"{np.abs(bl - L).max():.1e}  {dmin / (2 * RB):.4f}  {eps[0]:.5f} {eps[1]:.5f}")
    # find collision angle by bisection
    def margin2(th):
        return check(2, rot2(th))[1] - 2 * RB
    lo, hi = 0.01, np.radians(90)
    ths = np.linspace(lo, hi, 400)
    m = [margin2(x) for x in ths]
    k = next((i for i in range(len(m)) if m[i] < 0), None)
    if k:
        a, b = ths[k - 1], ths[k]
        for _ in range(60):
            c = (a + b) / 2
            a, b = (c, b) if margin2(c) > 0 else (a, c)
        bl, dmin, eps, _ = check(2, rot2(a))
        print(f"2D collision at theta = {np.degrees(a):.6f} deg, strain = {eps}, area change = {(1+eps[0])*(1+eps[1])-1:.5f}")

    print("\n3D tilt systems (even octahedra rotated by R, odd by R^T):")
    for name, axis in [("[001] (a0a0c-like)", [0, 0, 1]), ("[111] (a-a-a--like)", [1, 1, 1]), ("[011]", [0, 1, 1])]:
        ax = np.array(axis, float) / np.linalg.norm(axis)
        def res(th):
            return check(3, Rot.from_rotvec(th * ax).as_matrix(), ncell=2)
        ths = np.linspace(0.005, np.radians(90), 300)
        prev = None; coll = None
        for th in ths:
            r = res(th)
            if r is None:
                break
            if r[1] < 2 * RB - 1e-12:
                coll = (prev, th); break
            prev = th
        if coll:
            a, b = coll
            for _ in range(50):
                c = (a + b) / 2
                r = res(c)
                a, b = (c, b) if (r is not None and r[1] >= 2 * RB) else (a, c)
            r = res(a)
            print(f"  {name}: collision at {np.degrees(a):.4f} deg; max bar err {np.abs(r[0]-L).max():.1e}; strains {np.round(r[2],5)}; vol change {np.prod(1+r[2])-1:.5f}")
        r = res(np.radians(5))
        e = r[2]
        print(f"     small angle 5deg strains {np.round(e,6)}; ratio e_zz/e_xx = {e[2]/e[0]:.4f}")
