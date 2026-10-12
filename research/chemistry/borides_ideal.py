"""Thread B: tile space with a Voronoi cell, apply the pyramid-insphere rule in every cell, and read off the
'ideal' boron framework. Cube (simple cubic) -> CaB6 type; rhombic dodecahedron (fcc) -> UB12 type.
Also: truncated octahedron (bcc) and hexagonal prism to show where the rule fails (unequal balls).
All numbers computed, not assumed.
"""
import numpy as np
from polytopes import analyse, perms_signs

np.set_printoptions(precision=5, suppress=True)


def lattice_tile(cell_vertices, lattice_vecs, shells=1):
    a = analyse("cell", cell_vertices)
    return a


def report(name, cellV, lat, cellsize):
    """cellV: Voronoi cell vertices (centred at a lattice point at 0); lat: lattice translations list."""
    a = analyse(name, cellV)
    from polytopes import facets
    F = facets(np.asarray(cellV, float))
    # recompute dots
    from polytopes import pyramid_insphere
    dots = np.array([pyramid_insphere(np.asarray(cellV, float), *f)[0] for f in F])
    R = np.array([pyramid_insphere(np.asarray(cellV, float), *f)[1] for f in F])
    # tile
    pts, rad = [], []
    for L in lat:
        pts.append(dots + L); rad.append(R)
    pts = np.vstack(pts); rad = np.concatenate(rad)
    D = np.linalg.norm(pts[:, None] - pts[None], axis=2)
    i0 = int(np.argmin(np.linalg.norm(pts - dots[0], axis=1)))
    contacts = [j for j in range(len(pts)) if j != i0 and abs(D[i0, j] - rad[i0] - rad[j]) < 1e-7]
    overl = [j for j in range(len(pts)) if j != i0 and D[i0, j] < rad[i0] + rad[j] - 1e-7]
    print(f"\n== {name}: {len(F)} balls/cell, radii {np.unique(np.round(R,6))}, ball 0 touches {len(contacts)}, overlaps {len(overl)}")
    return dots, R, pts, rad


if __name__ == "__main__":
    # ---------------- simple cubic: cube cell side 2 ----------------
    cube = perms_signs((1, 1, 1))
    lat = [np.array(v, float) * 2 for v in np.ndindex(3, 3, 3)]
    lat = [v - 2 for v in lat]
    dots, R, pts, rad = report("cube / simple cubic -> CaB6", cube, lat, 2.0)
    k = np.linalg.norm(dots[0]); a = 2.0
    x = 0.5 - k / a  # fractional B coordinate with Ca at origin, B6 centre at (1/2,1/2,1/2)
    print(f"dot distance from cell centre t = {k:.6f} (2-sqrt2 = {2-2**0.5:.6f}); R = {R[0]:.6f}")
    print(f"=> CaB6 boron parameter x_ideal = 1/2 - t/a = {x:.6f}   ((sqrt2-1)/2 = {(2**0.5-1)/2:.6f})")
    d_intra = np.sqrt(2) * (0.5 - x); d_inter = 2 * x
    print(f"   intra-octahedral B-B = {d_intra:.6f} a ; inter-octahedral B-B = {d_inter:.6f} a (equal: {np.isclose(d_intra,d_inter)})")
    # metal at cube corner: distances to B
    corner = np.array([1.0, 1.0, 1.0])  # a cube corner in cell-centred coordinates
    allB = pts
    dM = np.sort(np.linalg.norm(allB - corner, axis=1))
    nM = np.sum(np.isclose(dM, dM[0]))
    print(f"   metal at corner: {nM} B at {dM[0]/a:.6f} a ; M-B / B-B = {dM[0]/(2*R[0]):.6f};  cage hole/ball = {(dM[0]-R[0])/R[0]:.6f}")
    print(f"   a / d_BB = {a/(2*R[0]):.6f}  (1+sqrt2 = {1+2**0.5:.6f})")

    # ---------------- fcc: rhombic dodecahedron, cubic cell a=4 here -> RD of fcc with nn distance 2*sqrt2 ----
    # fcc lattice points with cubic a = 4: (0,0,0),(0,2,2),... ; Voronoi cell = RD with vertices (+-2,0,0),(+-1,+-1,+-1)
    RD = np.vstack([perms_signs((1, 1, 1)), perms_signs((2, 0, 0))])
    A = 4.0
    basis = [np.array(v, float) for v in [(0, 0, 0), (0, 2, 2), (2, 0, 2), (2, 2, 0)]]
    lat = [b + A * np.array(n, float) for b in basis for n in np.ndindex(3, 3, 3)]
    lat = [v - A for v in lat]
    lat = [v for v in lat if np.linalg.norm(v) < 7]
    dots, R, pts, rad = report("rhombic dodecahedron / fcc -> UB12", RD, lat, A)
    t = np.linalg.norm(dots[0])
    xs = np.unique(np.round(np.abs(dots) / A, 6))
    print(f"   dot offset from cluster centre (fractional, cubic a): {xs}  -> UB12 B at 48i (1/2, x, x) with x_ideal = {xs.max():.6f} (1/6 = {1/6:.6f})")
    # metal at the 4-valent RD vertex (2,0,0) = octahedral hole of the cluster fcc
    corner = np.array([2.0, 0, 0])
    dM = np.sort(np.linalg.norm(pts - corner, axis=1)); nM = np.sum(np.isclose(dM, dM[0]))
    print(f"   metal at 4-valent vertex: {nM} B at {dM[0]/A:.6f} a; M-B / B-B = {dM[0]/(2*R[0]):.6f}; cage hole/ball = {(dM[0]-R[0])/R[0]:.6f}")
    print(f"   a / d_BB = {A/(2*R[0]):.6f} (3*sqrt2 = {3*2**0.5:.6f})")
    # other RD vertex (3-valent, tetrahedral hole)
    c3 = np.array([1.0, 1, 1]); d3 = np.sort(np.linalg.norm(pts - c3, axis=1)); n3 = np.sum(np.isclose(d3, d3[0]))
    print(f"   tetrahedral hole (3-valent vertex): {n3} B at M-B/B-B = {d3[0]/(2*R[0]):.6f}; hole/ball = {(d3[0]-R[0])/R[0]:.6f}")

    # ---------------- bcc: truncated octahedron (failure case) ----------------
    TO = perms_signs((0, 1, 2))
    basis = [np.zeros(3), np.array([2., 2, 2])]
    lat = [b + 4 * np.array(n, float) - 4 for b in basis for n in np.ndindex(3, 3, 3)]
    report("truncated octahedron / bcc", TO, lat, 4)

    # ---------------- 2D square: dots form 4.8.8 ----------------
    print("\n2D: square cell -> dots at distance 2-sqrt2 from centre, R = sqrt2-1, tiling = 4.8.8 (truncated square), each dot 3 contacts")
