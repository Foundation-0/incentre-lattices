"""Thread A: the 'pyramid insphere' rule on every polytope / kaleidoscope we can build.

Rule (user's): cut polytope P from its centre into pyramids over its facets; in each pyramid put the
point equidistant from all its walls (base facet + lateral walls through the centre). Ball radius R = that
distance. Report: number of balls, whether all edge-neighbour balls touch, and hole/ball = (t - R)/R where t =
distance of dots from centre (the largest central ball touching all balls).

Checked numerically: the pyramid insphere is the largest ball in conv(0, facet), found by linear programming
over centre and radius (scipy linprog), and tangency to every wall is then checked; not by formula. Then
compare with the classical radius-ratio value for the dual (coordination) polyhedron:
circumradius/(half edge) - 1.
"""
import itertools as it
import numpy as np
from scipy.spatial import ConvexHull
from scipy.optimize import linprog

phi = (1 + 5 ** 0.5) / 2


def perms_signs(v, even_only=False):
    out = set()
    idx = list(it.permutations(range(len(v))))
    if even_only:
        def parity(p):
            p = list(p); s = 0
            for i in range(len(p)):
                for j in range(i + 1, len(p)):
                    s += p[i] > p[j]
            return s % 2
        idx = [p for p in idx if parity(p) == 0]
    for p in idx:
        w = [v[i] for i in p]
        for signs in it.product([1, -1], repeat=len(v)):
            out.add(tuple(round(s * x, 12) for s, x in zip(signs, w)))
    return np.array(sorted(out), float)


def facets(V):
    """Return list of (unit normal n, offset h, vertex indices) for convex hull of V (merged coplanar)."""
    hull = ConvexHull(V)
    fac = {}
    for eq in hull.equations:
        n, h = eq[:-1], -eq[-1]
        key = tuple(np.round(np.r_[n, h], 8))
        fac[key] = (n / np.linalg.norm(n), h / np.linalg.norm(n))
    out = []
    for n, h in fac.values():
        idx = np.where(np.abs(V @ n - h) < 1e-7)[0]
        out.append((n, h, idx))
    return out


def ridges_of_facet(V, n, h, idx):
    """Ridges (codim-2 faces) of a facet: facets of the facet polytope, via hull in facet's affine frame."""
    P = V[idx]
    c = P.mean(0)
    # orthonormal basis of the facet hyperplane
    d = V.shape[1]
    B = np.linalg.svd(np.eye(d) - np.outer(n, n))[0][:, : d - 1]
    Q = (P - c) @ B
    if d - 1 == 1:
        return [[idx[np.argmin(Q[:, 0])]], [idx[np.argmax(Q[:, 0])]]]
    sub = facets(Q)
    return [[idx[i] for i in s[2]] for s in sub]


def pyramid_insphere(V, n, h, idx):
    """Max ball inside pyramid conv(0, facet) : LP over centre x, radius r. Returns centre, r."""
    d = V.shape[1]
    A, b = [], []
    # base: n.x + r <= h
    A.append(np.r_[n, 1.0]); b.append(h)
    for ridge in ridges_of_facet(V, n, h, idx):
        R = V[ridge]
        # wall = hyperplane through origin and ridge; normal m orthogonal to ridge points (span dim d-1)
        M = R
        _, s, vt = np.linalg.svd(M)
        m = vt[-1]
        # orient m so that facet centroid is on the inside (m.x >= 0 inside)
        if m @ V[idx].mean(0) < 0:
            m = -m
        A.append(np.r_[-m, 1.0]); b.append(0.0)  # -m.x + r <= 0
    res = linprog(c=np.r_[np.zeros(d), -1.0], A_ub=np.array(A), b_ub=np.array(b),
                  bounds=[(None, None)] * d + [(0, None)], method="highs")
    x, r = res.x[:d], res.x[d]
    # check equidistance from ALL walls (insphere exists) ?
    dists = [h - n @ x] + [abs(a[:d] @ x) for a in A[1:]]
    tangent_all = np.allclose(dists, r, atol=1e-7)
    return x, r, tangent_all


def analyse(name, V):
    V = np.asarray(V, float)
    V = V - V.mean(0)
    F = facets(V)
    dots, radii, tang = [], [], []
    for n, h, idx in F:
        x, r, ta = pyramid_insphere(V, n, h, idx)
        dots.append(x); radii.append(r); tang.append(ta)
    dots = np.array(dots); radii = np.array(radii)
    t = np.linalg.norm(dots, axis=1)
    D = np.linalg.norm(dots[:, None] - dots[None], axis=2)
    np.fill_diagonal(D, np.inf)
    # neighbour = facets sharing a ridge
    nb_touch, nb_total, overlap = 0, 0, False
    ridge_sets = [set(map(frozenset, [frozenset(r) for r in ridges_of_facet(V, *f)])) for f in F]
    for i in range(len(F)):
        for j in range(i + 1, len(F)):
            if ridge_sets[i] & ridge_sets[j]:
                nb_total += 1
                if abs(D[i, j] - radii[i] - radii[j]) < 1e-6:
                    nb_touch += 1
            if D[i, j] < radii[i] + radii[j] - 1e-6:
                overlap = True
    hole = np.min(t - radii)
    ratio = hole / radii.mean() if np.ptp(radii) < 1e-9 else hole / radii.max()
    return dict(name=name, n_balls=len(F), equal_R=np.ptp(radii) < 1e-9, ins_all=all(tang),
                touch=f"{nb_touch}/{nb_total}", overlap=overlap, ratio=ratio, R=radii, t=t)


def classic_ratio(Q):
    """radius ratio for coordination polyhedron with vertices Q (equal edges assumed): Rc/(e/2)-1."""
    Q = np.asarray(Q, float); Q = Q - Q.mean(0)
    Rc = np.linalg.norm(Q, axis=1).mean()
    D = np.linalg.norm(Q[:, None] - Q[None], axis=2); np.fill_diagonal(D, np.inf)
    return Rc / (D.min() / 2) - 1


def prism(n, h=None):
    a = 2 * np.sin(np.pi / n)  # edge for unit circumradius polygon
    h = a if h is None else h
    return [(np.cos(2 * np.pi * k / n), np.sin(2 * np.pi * k / n), s * h / 2) for k in range(n) for s in (1, -1)]


def antiprism(n):
    a = 2 * np.sin(np.pi / n)
    # height for equal edges
    hh = np.sqrt(a ** 2 - (2 * np.sin(np.pi / (2 * n))) ** 2)
    return [(np.cos(np.pi * k / n), np.sin(np.pi * k / n), (hh / 2) * (1 if k % 2 == 0 else -1)) for k in range(2 * n)]


def polygon(n):
    return [(np.cos(2 * np.pi * k / n), np.sin(2 * np.pi * k / n)) for k in range(n)]


SOLIDS = {
    "tetrahedron": [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)],
    "cube": perms_signs((1, 1, 1)),
    "octahedron": perms_signs((1, 0, 0)),
    "dodecahedron": np.vstack([perms_signs((1, 1, 1)), perms_signs((0, phi, 1 / phi), even_only=True)]),
    "icosahedron": perms_signs((0, 1, phi), even_only=True),
    "cuboctahedron": perms_signs((1, 1, 0)),
    "rhombic dodecahedron": np.vstack([perms_signs((1, 1, 1)), perms_signs((2, 0, 0))]),
    "truncated octahedron": perms_signs((0, 1, 2)),
    "triakis? (skip)": None,
    "trigonal prism (uniform)": prism(3),
    "hexagonal prism (uniform)": prism(6),
    "square antiprism (uniform)": antiprism(4),
    "triangular bipyramid": [(1, 0, 0), (-0.5, 3 ** 0.5 / 2, 0), (-0.5, -3 ** 0.5 / 2, 0), (0, 0, 1.2), (0, 0, -1.2)],
    "tesseract (4D)": perms_signs((1, 1, 1, 1)),
    "16-cell (4D)": perms_signs((1, 0, 0, 0)),
    "24-cell (4D)": perms_signs((1, 1, 0, 0)),
    "5-cell (4D)": [(1, 1, 1, -1 / 5 ** 0.5), (1, -1, -1, -1 / 5 ** 0.5), (-1, 1, -1, -1 / 5 ** 0.5),
                    (-1, -1, 1, -1 / 5 ** 0.5), (0, 0, 0, 4 / 5 ** 0.5)],
    "triangle (2D)": polygon(3), "square (2D)": polygon(4), "pentagon (2D)": polygon(5),
    "hexagon (2D)": polygon(6),
}

if __name__ == "__main__":
    print(f"{'polytope P':30s} {'balls':>5s} {'eqR':>4s} {'insph':>5s} {'touch':>7s} {'ovl':>4s} {'hole/ball':>10s}")
    rows = []
    for k, V in SOLIDS.items():
        if V is None:
            continue
        a = analyse(k, V)
        rows.append(a)
        print(f"{k:30s} {a['n_balls']:5d} {str(a['equal_R']):>4s} {str(a['ins_all']):>5s} {a['touch']:>7s} "
              f"{str(a['overlap']):>4s} {a['ratio']:10.4f}" + ("" if a['equal_R'] else f"   R range {a['R'].min():.3f}-{a['R'].max():.3f}"))
    print()
    print("Classical radius-ratio limits of coordination polyhedra (equal edges), for comparison:")
    CQ = {
        "triangle CN3": polygon(3), "square planar CN4": polygon(4),
        "tetrahedron CN4": SOLIDS["tetrahedron"], "octahedron CN6": SOLIDS["octahedron"],
        "trigonal prism CN6": prism(3), "cube CN8": SOLIDS["cube"],
        "square antiprism CN8": antiprism(4), "icosahedron CN12": SOLIDS["icosahedron"],
        "cuboctahedron CN12": SOLIDS["cuboctahedron"], "hexagonal prism CN12": prism(6),
        "truncated tetrahedron CN12": [v for v in perms_signs((3, 1, 1)) if np.prod(np.sign(v)) > 0],
        "dodecahedron CN20": SOLIDS["dodecahedron"], "truncated octahedron CN24": SOLIDS["truncated octahedron"],
        "truncated cube CN24": perms_signs((2 ** 0.5 - 1, 1, 1)),
        "16-cell (4D) CN8": SOLIDS["16-cell (4D)"], "24-cell (4D) CN24": SOLIDS["24-cell (4D)"],
        "tesseract (4D) CN16": SOLIDS["tesseract (4D)"],
    }
    for k, Q in CQ.items():
        print(f"  {k:30s} {classic_ratio(Q):.4f}")
