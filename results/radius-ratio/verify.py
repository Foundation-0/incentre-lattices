#!/usr/bin/env python3
"""Replay for results/radius-ratio: the radius-ratio theorem (Theorem 3 of preprint v2; Theorem 2 of v1),
hole/ball = classical radius ratio of the dot polytope.

Re-runs every case in radius_ratio.py (numpy + scipy: ConvexHull, HalfspaceIntersection, linprog),
checks the theorem's conclusions on all cases that satisfy hypothesis H1 or H2, and checks that
the documented counterexamples really violate the conclusion (must-fail controls), and asserts
the premises and the specific numbers that PROOF.md quotes:
  * tetragonal trapezohedra polar to square antiprisms of height 0.7x and 1.3x the uniform one:
    facet-transitive (the S8 rotoreflection (x, y, z) -> (R_{pi/4}(x, y), -z) and the sigma_d mirror
    y -> -y map the vertex set onto itself to 1e-12, and the group they generate is transitive on
    the 8 kites; pyramid inradii and dot distances are equal to within 1e-6), an insphere centred at O
    (all facets at one distance from O to within 1e-6), every pyramid has an inscribed ball, but adjacent
    dots are NOT all at distance 2r: exactly 8 of the 16 contacts stay closed and the others open,
  * the triangular prism and the cuboctahedron (not facet-transitive: pyramid inradii differ by more than 1e-3);
  * each counterexample lies outside H1 and outside H2, each by more than 1e-3 (a control that met a
    hypothesis and failed the conclusion would refute the theorem, so it must not be logged as a
    rejected control): reflecting the vertex set in some wall leaves a vertex more than 1e-3 off the
    set, and the angle between adjacent facet normals or the distance of the facets from O varies by
    more than 1e-3;
  * H1 is checked numerically: the vertex set of P is reflected in every wall (plane through O and
    a ridge) and must map to itself, for every case that is claimed under H1; the H2-only cases
    are asserted not to be mirror-symmetric in this sense (a vertex more than 1e-3 off the set);
  * each counterexample has an inscribed ball in every pyramid (the premise of the control);
  * the inscribed-ball test itself must reject a 2x4x6 cuboid, none of whose pyramids has an
    inscribed ball (must-fail control; the cuboid is neither a theorem case nor a counterexample);
  * prism inradii, trapezohedron contact ranges and deltahedra orbit radii equal PROOF.md's numbers.
Writes certificate.json. Usage: python verify.py
"""
import contextlib
import io
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "tools"))
sys.path.insert(0, str(HERE))
import certlib  # noqa: E402
import radius_ratio as rr  # noqa: E402

TOL = 1e-6
H1_TOL = 1e-9          # a reflection of the vertex set must land on the vertex set to this distance
SYM_TOL = 1e-12        # a D4d generator must map a stretched trapezohedron's vertex set (and facet normals) onto itself to this distance
QUOTE_TOL = 6e-4       # PROOF.md quotes these numbers to three or four decimals
# Numbers quoted in PROOF.md (prism inradii for unit edge; adjacent-dot ranges in units of 2r; orbit radii for unit edge)
PRISM_INRADII_UNIT_EDGE = (0.1340, 0.1667)
TRAPEZOHEDRON_MAX_OVER_2R = {"x0.7": 1.017, "x1.3": 1.006}
DELTAHEDRA = {"J13": (0.851, 0.526), "J84": (0.930, 0.677), "J51": (0.764, 0.996)}
H2_FORMULA = ("pentagonal icositetrahedron -> snub cube", "trapezohedron (dual of uniform antiprism)")
COUNTER_PREFIX = ("trapezohedron (dual of stretched", "triangular prism", "cuboctahedron (not")
DETECTOR_CONTROL = "cuboid 2x4x6 (no inscribed ball)"   # must-fail input for the inscribed-ball test
# D4d generators of the stretched trapezohedra (axis z; an antiprism vertex on the x-axis): the S8 rotoreflection
# (x, y, z) -> (R_{pi/4}(x, y), -z) and the sigma_d mirror y -> -y
_C, _S = np.cos(np.pi / 4), np.sin(np.pi / 4)
D4D_GENERATORS = {"S8": np.array([[_C, -_S, 0], [_S, _C, 0], [0, 0, -1]]), "sigma_d": np.diag([1.0, -1.0, 1.0])}


def facet_orbit(V, name):
    """Require that each D4d generator maps the vertex set V onto itself and permutes the facets; return the size of
    the orbit of one facet under the group they generate, and the number of facets."""
    fac = rr.facets_of(V)
    normals = np.array([f[0] for f in fac])
    perms = []
    for gname, g in D4D_GENERATORS.items():
        err = float(np.linalg.norm((V @ g.T)[:, None] - V[None], axis=2).min(axis=1).max())
        certlib.require(err < SYM_TOL, f"{gname} does not map {name} onto itself: a vertex lands {err:.2e} off the set")
        d = np.linalg.norm((normals @ g.T)[:, None] - normals[None], axis=2)
        perm = d.argmin(axis=1)
        certlib.require(d.min(axis=1).max() < SYM_TOL and len(set(perm.tolist())) == len(fac), f"{gname} does not permute the facets of {name}")
        perms.append(perm)
    orbit, todo = {0}, [0]
    while todo:
        a = todo.pop()
        for perm in perms:
            b = int(perm[a])
            if b not in orbit:
                orbit.add(b)
                todo.append(b)
    return len(orbit), len(fac)


def theorem_holds(i):
    """The conclusion of the radius-ratio theorem (Theorem 3 of preprint v2; Theorem 2 of v1) for one analysed polytope."""
    return (i["insphere_in_every_pyramid"] and i["r_spread"] < TOL and i["D_spread"] < TOL
            and abs(i["adjacent_min_over_2r"] - 1) < TOL and abs(i["adjacent_max_over_2r"] - 1) < TOL
            and i["min_pair_over_2r"] > 1 - TOL
            and (i["classical"] is None or abs(i["hole_over_ball"] - i["classical"]) < TOL))


def main():
    with contextlib.redirect_stdout(io.StringIO()):
        delta = rr.run_cases()
    cases, counter, detector, log = [], [], [], []
    for i in rr.RESULTS:
        name = i["name"]
        if name == DETECTOR_CONTROL:
            # must-fail control for the inscribed-ball test itself; neither a theorem case nor a counterexample
            certlib.must_fail("inscribed-ball test: " + name, i["insphere_in_every_pyramid"], log)
            detector.append({"name": name, "all_pyramids_have_inscribed_ball": i["insphere_in_every_pyramid"]})
            continue
        if name.startswith(COUNTER_PREFIX):
            # must-fail controls: the theorem's conclusion must be violated
            certlib.must_fail("counterexample: " + name, theorem_holds(i), log)
            # premises of the control, so a geometry bug that dropped the inscribed balls cannot pass
            certlib.require(i["insphere_in_every_pyramid"], f"control premise: every pyramid has an inscribed ball ({name})")
            certlib.require(i["h1_max_mirror_error"] > 1e-3, f"control must lie outside H1 ({name})")
            certlib.require(i["normal_angle_spread"] > 1e-3 or i["facet_offset_spread"] > 1e-3, f"control must lie outside H2 ({name})")
            extra = {}
            if name.startswith("trapezohedron"):
                key = "x0.7" if "x0.7" in name else "x1.3"
                orbit, nfac = facet_orbit(i["vertices"], name)
                certlib.require(orbit == nfac == i["facets"] == 8, f"stretched trapezohedron is facet-transitive ({name}): orbit {orbit} of {nfac} facets")
                certlib.require(i["r_spread"] < TOL and i["D_spread"] < TOL, f"equal pyramid inradii and dot distances ({name})")
                certlib.require(i["facet_offset_spread"] < TOL, f"stretched trapezohedron has an insphere centred at O ({name})")
                closed = int(np.sum(np.abs(i["adjacent_over_2r"] - 1) < TOL))
                certlib.require(closed == 8 and len(i["adjacent_over_2r"]) == 16,
                                f"8 of the 16 contacts stay closed ({name}): {closed} of {len(i['adjacent_over_2r'])}")
                certlib.require(abs(i["adjacent_min_over_2r"] - 1) < TOL, f"some contact stays closed ({name})")
                certlib.require(abs(i["adjacent_max_over_2r"] - TRAPEZOHEDRON_MAX_OVER_2R[key]) < QUOTE_TOL,
                                f"PROOF.md contact range for {name}: {i['adjacent_max_over_2r']}")
                extra = {"facet_transitive_under_D4d_generators": True, "insphere_centred_at_O": True,
                         "ridges": len(i["adjacent_over_2r"]), "contacts_closed": closed}
            else:
                certlib.require(i["r_spread"] > 1e-3, f"control is not facet-transitive ({name})")
                if name.startswith("triangular prism"):
                    edge = 3 ** 0.5   # the prism in radius_ratio.py has edge sqrt(3); PROOF.md quotes unit edge
                    lo, hi = i["r_min"] / edge, i["r_max"] / edge
                    certlib.require(abs(lo - PRISM_INRADII_UNIT_EDGE[0]) < QUOTE_TOL and abs(hi - PRISM_INRADII_UNIT_EDGE[1]) < QUOTE_TOL,
                                    f"PROOF.md prism inradii: {lo}, {hi}")
                    extra = {"pyramid_inradii_unit_edge": [certlib.sig(lo, 4), certlib.sig(hi, 4)]}
            counter.append({"name": name, "adjacent_min_over_2r": certlib.sig(i["adjacent_min_over_2r"], 5),
                            "adjacent_max_over_2r": certlib.sig(i["adjacent_max_over_2r"], 5),
                            "all_pyramids_have_inscribed_ball": i["insphere_in_every_pyramid"],
                            "r_spread_above_1e-3": i["r_spread"] > 1e-3,
                            "h1_walls_are_mirrors": False, "h2_holds": False, **extra})
            continue
        certlib.require(theorem_holds(i), f"theorem fails on {name}: {i}")
        row = {"name": name, "dim": i["dim"], "facets": i["facets"],
               "hole_over_ball": certlib.sig(i["hole_over_ball"], 8),
               "classical_radius_ratio": None if i["classical"] is None else certlib.sig(i["classical"], 8),
               "all_adjacent_dots_at_distance_2r": True}
        h1 = i["h1_max_mirror_error"] < H1_TOL
        if name in H2_FORMULA:
            certlib.require(i["h1_max_mirror_error"] > 1e-3, f"{name} is claimed under H2 only, but all its walls are mirrors to within 1e-3")
            row["hypothesis"] = "H2"
        else:
            certlib.require(h1, f"H1 fails on {name}: reflecting the vertices in a wall leaves a vertex off the set by {i['h1_max_mirror_error']:.2e}")
            row["hypothesis"] = "H1"
            row["h1_walls_checked"] = i["h1_walls_checked"]
        row["h1_walls_are_mirrors"] = h1
        if name in H2_FORMULA:
            certlib.require(abs(i["equal_dihedral_formula"] - i["hole_over_ball"]) < TOL
                            and i["normal_angle_spread"] < TOL and i["facet_offset_spread"] < TOL,
                            f"equal-dihedral formula 1/sin(phi/2)-1 fails on {name}")
            row["equal_dihedral_formula_1_over_sin_half_phi_minus_1"] = certlib.sig(i["equal_dihedral_formula"], 8)
        cases.append(row)
    certlib.require(len(detector) == 1, "the inscribed-ball control (cuboid) did not run")
    # inscribability of the equilateral deltahedras (CN7, CN8 dodecahedral, CN9)
    j13, j51, j84 = delta["J13"], delta["J51"], delta["J84"]
    certlib.require(abs(j13[0] - j13[1]) > 0.1 and abs(j51[0] - j51[1]) > 0.1 and abs(j84[0] - j84[1]) > 0.1 and j84[2] < 1e-9
                    and j84[3] == 12, "deltahedra inscribability check changed")
    for key, got in (("J13", j13[:2]), ("J84", j84[:2]), ("J51", j51[:2])):
        certlib.require(all(abs(g - q) < QUOTE_TOL for g, q in zip(got, DELTAHEDRA[key])),
                        f"PROOF.md orbit radii for {key}: {got}")
    cert = {
        "result": "radius-ratio",
        "hypotheses": {"H1": "every wall (plane through O and a ridge) is a mirror of P, and one pyramid has an inscribed ball (mirror property checked numerically by reflecting the vertex set in every wall)",
                       "H2": "P has an insphere centred at O and all dihedral angles are equal"},
        "cases_satisfying_theorem": cases,
        "counterexamples_theorem_conclusion_fails": counter,
        "inscribed_ball_test_control": detector[0],
        "deltahedra_not_inscribed": {
            "J13_equator_radius_vs_apex_distance": [certlib.sig(x, 4) for x in j13],
            "J51_prism_vertex_radius_vs_cap_radius": [certlib.sig(x, 4) for x in j51],
            "J84_orbit_radii": [certlib.sig(x, 4) for x in j84[:2]],
            "J84_max_edge_error_below_1e-9": j84[2] < 1e-9,
            "J84_faces": j84[3]},
        "must_fail_controls_rejected": log,
        "tolerance": TOL,
        "inputs": certlib.inputs(HERE / "radius_ratio.py", HERE / "verify.py", HERE.parents[1] / "tools" / "certlib.py"),
    }
    certlib.write_json(HERE / "certificate.json", cert)
    print(f"radius-ratio: OK ({len(cases)} cases satisfy the theorem; {len(counter)} counterexample controls and "
          f"{len(detector)} inscribed-ball control rejected)")


if __name__ == "__main__":
    main()
