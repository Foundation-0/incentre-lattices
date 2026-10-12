"""Check of Theorem 2 (hole/ball ratio of the pyramid-incentre ("dot") construction).

For a convex polytope P containing O, cut P into the pyramids conv(O, F) over its facets F and
find, by linear programming, the largest ball in each pyramid. We report:
  - whether that ball touches every facet of its pyramid (an inscribed ball exists),
  - the radii r_F and the distances D_F = |O d_F| (equal?),
  - the distance between dots of ridge-adjacent facets divided by 2r (contact?),
  - hole/ball = (D - r)/r and the classical radius-ratio value of the dot polytope Q,
  - the equal-dihedral formula 1/sin(phi/2) - 1 where it applies.
Cases: Platonic solids, rhombic dodecahedron, two Catalan solids, bipyramids of several heights
(mirror walls, hypothesis H1), the pentagonal icositetrahedron (dual of the snub cube) and the
tetragonal trapezohedron (equal dihedral angles, hypothesis H2, only at one shape), the
triangular prism (not isohedral), and the regular 4-polytopes.
Also: inscribability tests for the equilateral CN7/CN8/CN9 deltahedra.
Run: python radius_ratio.py
"""
import itertools
import numpy as np
from scipy.optimize import linprog, least_squares
from scipy.spatial import ConvexHull, HalfspaceIntersection

phi_ = (1 + 5 ** 0.5) / 2


def facets_of(V):
    """group hull simplices into facets; return list of (unit normal, offset, vertex index set)."""
    H = ConvexHull(V)
    fac = []
    for eq, simp in zip(H.equations, H.simplices):
        nrm, off = eq[:-1], -eq[-1]
        for f in fac:
            if np.allclose(f[0], nrm, atol=1e-9) and abs(f[1] - off) < 1e-9:
                f[2].update(simp.tolist())
                break
        else:
            fac.append([nrm, off, set(simp.tolist())])
    return fac


RESULTS = []   # one info dict per analysed polytope (read by verify.py)


def analyse(V, name, Qclassical=None):
    V = np.asarray(V, float)
    n = V.shape[1]
    fac = facets_of(V)
    dots, rad, touches_all = [], [], []
    adj = []
    walls = []   # unit normals of the walls (planes through O and a ridge), one per ridge
    for a, (na, ha, Sa) in enumerate(fac):
        A = [na]; b = [ha]
        cen = V[list(Sa)].mean(0)
        for g, (ng, hg, Sg) in enumerate(fac):
            if g == a:
                continue
            common = sorted(Sa & Sg)
            if len(common) < n - 1:
                continue
            M = V[common]
            if np.linalg.matrix_rank(M - M[0]) < n - 2:
                continue
            # wall through O and the ridge
            _, s, vh = np.linalg.svd(M)
            w = vh[-1]
            if w @ cen > 0:
                w = -w
            A.append(w); b.append(0.0)
            if a < g:
                adj.append((a, g))
                walls.append(w / np.linalg.norm(w))
        A = np.array(A); b = np.array(b)
        norms = np.linalg.norm(A, axis=1)
        res = linprog(np.r_[np.zeros(n), -1], A_ub=np.c_[A, norms], b_ub=b,
                      bounds=[(None, None)] * n + [(0, None)], method="highs")
        x, r = res.x[:n], res.x[n]
        slack = b - A @ x - r * norms
        touches_all.append(bool(np.all(np.abs(slack) < 1e-7)))
        dots.append(x); rad.append(r)
    dots = np.array(dots); rad = np.array(rad)
    D = np.linalg.norm(dots, axis=1)
    r = rad.mean()
    adjd = np.array([np.linalg.norm(dots[a] - dots[g]) for a, g in adj])
    # global non-overlap check
    allD = np.linalg.norm(dots[:, None] - dots[None], axis=2) + np.eye(len(dots)) * 9
    ratio = (D.mean() - r) / r
    line = (f"{name:38s} facets={len(fac):3d} insphere-in-every-pyramid={all(touches_all)!s:5s} "
            f"r spread={np.ptp(rad):.1e} D spread={np.ptp(D):.1e} adjacent/2r in [{adjd.min()/(2*r):.4f},{adjd.max()/(2*r):.4f}] "
            f"min pair/2r={allD.min()/(2*r):.4f} hole/ball={ratio:.4f}")
    if Qclassical is not None:
        line += f"  classical(Q)={Qclassical:.4f}"
    print(line)
    phis = np.array([np.arccos(np.clip(fac[a][0] @ fac[g][0], -1, 1)) for a, g in adj])  # angle between normals
    offs = np.array([f[1] for f in fac])
    # H1: reflect the vertex set in every wall; each wall is a mirror of P iff the image is the vertex set again
    h1_err = 0.0
    for w in walls:
        Vr = V - 2 * np.outer(V @ w, w)
        h1_err = max(h1_err, float(np.linalg.norm(Vr[:, None] - V[None], axis=2).min(axis=1).max()))
    info = dict(name=name, h1_max_mirror_error=h1_err, h1_walls_checked=len(walls),
                r_min=float(rad.min()), r_max=float(rad.max()), normal_angle_spread=float(np.ptp(phis)), facet_offset_spread=float(np.ptp(offs)),
                equal_dihedral_formula=float(1 / np.sin(phis.mean() / 2) - 1), dim=n, facets=len(fac), insphere_in_every_pyramid=all(touches_all),
                r_spread=float(np.ptp(rad)), D_spread=float(np.ptp(D)),
                adjacent_min_over_2r=float(adjd.min() / (2 * r)), adjacent_max_over_2r=float(adjd.max() / (2 * r)),
                min_pair_over_2r=float(allD.min() / (2 * r)), hole_over_ball=float(ratio),
                classical=None if Qclassical is None else float(Qclassical))
    RESULTS.append(info)
    return ratio, dots, r


def classical(Qverts):
    """circumradius/(edge/2) - 1 for an inscribed equal-edged polytope (edge = shortest distance)."""
    Q = np.asarray(Qverts, float)
    Rc = np.linalg.norm(Q, axis=1)
    d = np.linalg.norm(Q[:, None] - Q[None], axis=2) + np.eye(len(Q)) * 99
    return Rc.mean() / (d.min() / 2) - 1


def perms_signs(v, even_only=False):
    out = set()
    for p in itertools.permutations(range(3)):
        sign = np.linalg.det(np.eye(3)[list(p)])
        if even_only and sign < 0:
            continue
        for s in itertools.product((1, -1), repeat=3):
            out.add(tuple(s[k] * v[p[k]] for k in range(3)))
    return np.array(sorted(out))


tet = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
cube = np.array(list(itertools.product((1, -1), repeat=3)), float)
octa = np.vstack([np.eye(3), -np.eye(3)])
ico = np.array([p for p in perms_signs((0, 1, phi_), even_only=True)])
dod = np.vstack([cube, perms_signs((0, 1 / phi_, phi_), even_only=True)])
rhd = np.vstack([cube, 2 * octa])
cubocta = perms_signs((1, 1, 0))
trunc_oct = perms_signs((0, 1, 2))


def polar(Q):
    """polar dual {x : x.q <= 1} of a polytope with O inside (vertices of the dual)."""
    hs = np.c_[Q, -np.ones(len(Q))]
    return HalfspaceIntersection(hs, np.zeros(Q.shape[1])).intersections


def run_cases():
    """Re-run every case of Theorem 2; fills RESULTS and returns the deltahedra data."""
    print("== H1: every wall is a mirror (Wythoff case) ==")
    analyse(tet, "tetrahedron -> Q tetrahedron (CN4)", classical(tet))
    analyse(cube, "cube -> Q octahedron (CN6)", classical(octa))
    analyse(octa, "octahedron -> Q cube (CN8)", classical(cube))
    analyse(dod, "dodecahedron -> Q icosahedron (CN12)", classical(ico))
    analyse(ico, "icosahedron -> Q dodecahedron (CN20)", classical(dod))
    analyse(rhd, "rhombic dodecahedron -> Q cuboctahedron", classical(cubocta))
    analyse(polar(trunc_oct), "tetrakis hexahedron -> Q trunc. octahedron", classical(trunc_oct))
    analyse(polar(perms_signs((1, 1, 1 + 2 ** 0.5))), "triakis octahedron -> Q truncated cube",
            classical(perms_signs((1, 1, 1 + 2 ** 0.5))))
    for k, hgt in [(3, 0.4), (3, 1.0), (3, 2.5), (4, 0.6), (4, 1.7), (6, 1.3)]:
        eq = [(np.cos(2 * np.pi * m / k), np.sin(2 * np.pi * m / k), 0) for m in range(k)]
        bip = np.array(eq + [(0, 0, hgt), (0, 0, -hgt)])
        side = 2 * np.sin(np.pi / k)  # uniform prism: circumradius^2 = (side/2/sin(pi/k))^2 + (side/2)^2
        prism_cl = np.sqrt((side / 2 / np.sin(np.pi / k)) ** 2 + (side / 2) ** 2) / (side / 2) - 1
        analyse(bip, f"{k}-bipyramid h={hgt} -> uniform {k}-prism", prism_cl)
    print("\n== H2: tangential with all dihedral angles equal (P = polar of inscribed equal-edged Q) ==")
    t = 1.839286755214161  # tribonacci constant
    ev = [p for p in perms_signs((1, 1 / t, t), even_only=True) if np.sum(np.array(p) > 0) % 2 == 0]
    odd = []
    for p in itertools.permutations(range(3)):
        if np.linalg.det(np.eye(3)[list(p)]) > 0:
            continue
        for s in itertools.product((1, -1), repeat=3):
            if sum(x > 0 for x in s) % 2 == 1:
                v = (1, 1 / t, t)
                odd.append(tuple(s[k] * v[p[k]] for k in range(3)))
    snub = np.array(ev + odd, float)
    print(f"  snub cube vertices: {len(snub)}, norms spread {np.ptp(np.linalg.norm(snub,axis=1)):.1e}")
    analyse(polar(snub), "pentagonal icositetrahedron -> snub cube", classical(snub))

    def antiprism(z):
        return np.array([(np.cos(m * np.pi / 4), np.sin(m * np.pi / 4), z * (-1) ** m) for m in range(8)])
    # uniform square antiprism: lateral edge = base edge
    s8 = 2 * np.sin(np.pi / 8)
    zu = 0.5 * np.sqrt(s8 ** 2 - (2 * np.sin(np.pi / 8) * 0 + np.linalg.norm(antiprism(0)[0] - antiprism(0)[1])) ** 2 + s8 ** 2)
    # solve lateral edge = side for z
    f = lambda z: np.linalg.norm(antiprism(z)[0] - antiprism(z)[1]) - np.linalg.norm(antiprism(z)[0] - antiprism(z)[2])
    from scipy.optimize import brentq
    zu = brentq(f, 0.01, 2)
    for z in (0.7 * zu, zu, 1.3 * zu):
        Q = antiprism(z)
        tag = "uniform" if z == zu else f"stretched x{z/zu:.1f}"
        analyse(polar(Q), f"trapezohedron (dual of {tag} antiprism)", classical(Q) if z == zu else None)
    print("\n== outside the hypotheses ==")
    prism = np.array([(np.cos(2 * np.pi * m / 3), np.sin(2 * np.pi * m / 3), z) for m in range(3) for z in (np.sin(np.pi / 3), -np.sin(np.pi / 3))])
    analyse(prism, "triangular prism (square sides)")
    analyse(cubocta, "cuboctahedron (not facet-transitive)")
    print("\n== 4D, regular polytopes (H1) ==")
    tess = np.array(list(itertools.product((1, -1), repeat=4)), float)
    c16 = np.vstack([np.eye(4), -np.eye(4)])
    c24 = np.array([v for v in itertools.product((1, -1, 0), repeat=4) if sum(abs(x) for x in v) == 2], float)
    s5 = np.vstack([np.eye(5) - 0.2, ])  # 5-cell in the hyperplane sum=0 -> project to R^4
    B = np.linalg.svd(np.ones((1, 5)))[2][1:]
    s5 = (np.eye(5) - 0.2) @ B.T
    analyse(s5, "5-cell -> 5-cell", classical(s5))
    analyse(tess, "tesseract -> 16-cell", classical(c16))
    analyse(c16, "16-cell -> tesseract", classical(tess))
    analyse(c24, "24-cell -> 24-cell (dual position)", classical(polar(c24)))

    print("\n== equilateral deltahedra for CN7, CN8 (dodecahedral), CN9: are they inscribed? ==")
    # pentagonal bipyramid J13
    rq = 1 / (2 * np.sin(np.pi / 5)); zap = np.sqrt(1 - rq ** 2)
    print(f"  J13 pentagonal bipyramid (CN7): equator radius {rq:.4f}, apex distance {zap:.4f} -> not inscribed")
    # tricapped trigonal prism J51 (unit edges); centre fixed by D3h symmetry at the origin
    rp = 1 / np.sqrt(3); capd = rp * np.cos(np.pi / 3) + 1 / np.sqrt(2)
    print(f"  J51 tricapped trigonal prism (CN9): prism-vertex radius {np.hypot(rp,0.5):.4f}, cap radius {capd:.4f} -> not inscribed")
    # snub disphenoid J84 with D2d symmetry: orbit A (+-p,0,z1),(0,+-p,-z1); orbit B (0,+-q,z2),(+-q,0,-z2)
    def pts(x):
        p, z1, q, z2 = x
        return np.array([(p, 0, z1), (-p, 0, z1), (0, p, -z1), (0, -p, -z1),
                         (0, q, z2), (0, -q, z2), (q, 0, -z2), (-q, 0, -z2)])
    def res(x):
        P = pts(x); H = ConvexHull(P)
        E = set()
        for s in H.simplices:
            for a, b in itertools.combinations(sorted(s), 2):
                E.add((a, b))
        return np.array([np.linalg.norm(P[a] - P[b]) - 1 for a, b in sorted(E)] + [0] * (18 - len(E)))
    sol = least_squares(res, [0.6, 0.6, 0.9, 0.2])
    P = pts(sol.x); H = ConvexHull(P)
    print(f"  J84 snub disphenoid (CN8 dodecahedral): {len(H.simplices)} faces, max |edge-1| = {np.abs(res(sol.x)).max():.1e}, "
          f"orbit radii {np.linalg.norm(P[0]):.4f} and {np.linalg.norm(P[4]):.4f} -> {'inscribed' if abs(np.linalg.norm(P[0])-np.linalg.norm(P[4]))<1e-6 else 'not inscribed'}")
    return {"J13": (float(rq), float(zap)), "J51": (float(np.hypot(rp, 0.5)), float(capd)),
            "J84": (float(np.linalg.norm(P[0])), float(np.linalg.norm(P[4])), float(np.abs(res(sol.x)).max()), len(H.simplices))}


if __name__ == "__main__":
    run_cases()
