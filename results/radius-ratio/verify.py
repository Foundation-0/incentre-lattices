#!/usr/bin/env python3
"""Replay for results/radius-ratio (Theorem 2, hole/ball = classical radius ratio of the dot polytope).

Re-runs every case in radius_ratio.py (numpy + scipy: ConvexHull, HalfspaceIntersection, linprog),
checks the theorem's conclusions on all cases that satisfy hypothesis H1 or H2, and checks that
the documented counterexamples really violate the conclusion (must-fail controls):
  * tetragonal trapezohedra polar to square antiprisms of height 0.7x and 1.3x the uniform one
    (facet-transitive, an insphere exists, every pyramid has an inscribed ball, but adjacent
    dots are NOT at distance 2r: the contacts open),
  * the triangular prism and the cuboctahedron (not facet-transitive).
Writes certificate.json. Usage: python verify.py
"""
import contextlib
import io
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "tools"))
sys.path.insert(0, str(HERE))
import certlib  # noqa: E402
import radius_ratio as rr  # noqa: E402

TOL = 1e-6
H2_FORMULA = ("pentagonal icositetrahedron -> snub cube", "trapezohedron (dual of uniform antiprism)")
COUNTER_PREFIX = ("trapezohedron (dual of stretched", "triangular prism", "cuboctahedron (not")


def theorem_holds(i):
    """The conclusion of Theorem 2 for one analysed polytope."""
    return (i["insphere_in_every_pyramid"] and i["r_spread"] < TOL and i["D_spread"] < TOL
            and abs(i["adjacent_min_over_2r"] - 1) < TOL and abs(i["adjacent_max_over_2r"] - 1) < TOL
            and i["min_pair_over_2r"] > 1 - TOL
            and (i["classical"] is None or abs(i["hole_over_ball"] - i["classical"]) < TOL))


def main():
    with contextlib.redirect_stdout(io.StringIO()):
        delta = rr.run_cases()
    cases, counter, log = [], [], []
    for i in rr.RESULTS:
        name = i["name"]
        if name.startswith(COUNTER_PREFIX):
            # must-fail controls: the theorem's conclusion must be violated
            certlib.must_fail("counterexample: " + name, theorem_holds(i), log)
            counter.append({"name": name, "adjacent_min_over_2r": certlib.sig(i["adjacent_min_over_2r"], 5),
                            "adjacent_max_over_2r": certlib.sig(i["adjacent_max_over_2r"], 5),
                            "all_pyramids_have_inscribed_ball": i["insphere_in_every_pyramid"],
                            "r_spread_above_1e-3": i["r_spread"] > 1e-3})
            continue
        certlib.require(theorem_holds(i), f"theorem fails on {name}: {i}")
        row = {"name": name, "dim": i["dim"], "facets": i["facets"],
               "hole_over_ball": certlib.sig(i["hole_over_ball"], 8),
               "classical_radius_ratio": None if i["classical"] is None else certlib.sig(i["classical"], 8),
               "all_adjacent_dots_at_distance_2r": True}
        if name in H2_FORMULA:
            certlib.require(abs(i["equal_dihedral_formula"] - i["hole_over_ball"]) < TOL
                            and i["normal_angle_spread"] < TOL and i["facet_offset_spread"] < TOL,
                            f"equal-dihedral formula 1/sin(phi/2)-1 fails on {name}")
            row["equal_dihedral_formula_1_over_sin_half_phi_minus_1"] = certlib.sig(i["equal_dihedral_formula"], 8)
        cases.append(row)
    # inscribability of the equilateral deltahedras (CN7, CN8 dodecahedral, CN9)
    j13, j51, j84 = delta["J13"], delta["J51"], delta["J84"]
    certlib.require(abs(j13[0] - j13[1]) > 0.1 and abs(j51[0] - j51[1]) > 0.1 and abs(j84[0] - j84[1]) > 0.1 and j84[2] < 1e-9
                    and j84[3] == 12, "deltahedra inscribability check changed")
    cert = {
        "result": "radius-ratio",
        "hypotheses": {"H1": "every wall (plane through O and a ridge) is a mirror of P, and one pyramid has an inscribed ball",
                       "H2": "P has an insphere centred at O and all dihedral angles are equal"},
        "cases_satisfying_theorem": cases,
        "counterexamples_theorem_conclusion_fails": counter,
        "deltahedra_not_inscribed": {
            "J13_equator_radius_vs_apex_distance": [certlib.sig(x, 4) for x in j13],
            "J51_prism_vertex_radius_vs_cap_radius": [certlib.sig(x, 4) for x in j51],
            "J84_orbit_radii": [certlib.sig(x, 4) for x in j84[:2]],
            "J84_max_edge_error_below_1e-9": j84[2] < 1e-9,
            "J84_faces": j84[3]},
        "must_fail_controls_rejected": log,
        "tolerance": TOL,
        "inputs": certlib.inputs(HERE / "radius_ratio.py"),
    }
    certlib.write_json(HERE / "certificate.json", cert)
    print(f"radius-ratio: OK ({len(cases)} cases satisfy the theorem; {len(log)} counterexample controls rejected)")


if __name__ == "__main__":
    main()
