#!/usr/bin/env python3
"""Replay for results/fold-3d-cubic (SLOW target, about 1 minute).

The cubic-restricted [111] fold of the 3D contact framework (octahedra rotated alternately by
+-theta about [111], lattice kept cubic): exact values at 60 degrees, the coordination sequence
of the 60-degree contact net, and a grid scan showing that the cubic lock is not a jam
(trigonal lattices keep every bar and reach higher density).  Finite checks on a grid; the
grid connectivity is NOT a proof of a continuous path (see results/fold-close-packing for the
only path that is checked for continuity at a stated step size).
Must-fail controls: claiming 60 degrees as a jam, a wrong lattice constant, a wrong sequence.
Writes certificate.json.  Usage: python verify.py
"""
import pathlib
import sys

import numpy as np
import sympy as sp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "tools"))
sys.path.insert(0, str(HERE))
import certlib  # noqa: E402
import fold3d as eng  # noqa: E402

SMS = [6, 17, 42, 68, 110, 153, 212, 272, 348, 425]   # RCSR sms coordination sequence (from rcsr.net/data/3dall.txt)


def main():
    first, dev = eng.cubic_path()
    certlib.require(abs(first - 60.0) < 0.05 and dev < 1e-12, "cubic path: first contact at 60 degrees, bars preserved")
    ex = eng.exact_60()
    certlib.require(ex["a60_minus_8c_over_3"] == 0 and ex["all_contacts_exact"] and ex["contacts"] == [6] * 6
                    and abs(ex["density_minus_closed_form"]) < 1e-25, "exact values at 60 degrees")
    cs = eng.coordination_sequence(ex["a60"])
    certlib.require(cs == SMS, "coordination sequence")
    tri = eng.trigonal_scan()
    certlib.require(tri["connected"] == tri["feasible"] and tri["best_density"] > ex["density"] + 0.1 and tri["bar_dev"] < 1e-12,
                    "trigonal scan")

    log = []
    certlib.must_fail("60-degree cubic state claimed to be a jam (no denser feasible state)",
                      tri["best_density"] <= ex["density"] + 1e-9, log)
    P, lab = eng.cell3(np.pi / 3, (ex["a60"] + 0.01) * np.eye(3))
    Dm = np.linalg.norm(P - P[lab.index(((0, 0, 0), 0))], axis=1)
    certlib.must_fail("lattice constant a(60) + 0.01 claimed to keep six contacts at 2R for a dot",
                      int(((Dm > 1e-9) & (Dm < 2 * eng.R + 1e-6)).sum()) == 6, log)
    certlib.must_fail("coordination sequence claimed to match fcc (12, 42, 92, ...)", cs[:3] == [12, 42, 92], log)

    cert = {
        "result": "fold-3d-cubic",
        "model": "octahedra of 6 dots (c = 2 - sqrt2) rotated alternately by +-theta about [111]; lattice cubic; bars of length 2R kept",
        "cubic_path": {"first_nonbonded_contact_deg_grid_3000": certlib.sig(first, 5), "bar_length_deviation_below_1e-12": True},
        "exact_at_60_degrees": {
            "a": "8c/3 = 8(2 - sqrt(2))/3", "a_value": certlib.sig(ex["a60"], 10),
            "contacts_per_ball": ex["contacts"], "all_contacts_exactly_2R": True,
            "next_distance_over_2R": [certlib.sig(x, 5) for x in ex["next_over_2R"]],
            "density": "27 pi sqrt(2)/256", "density_value": certlib.sig(ex["density"], 10),
            "start_density_pi_R3": certlib.sig(ex["start_density"], 10)},
        "coordination_sequence_first_10": cs,
        "coordination_sequence_equals_RCSR_sms_first_10_terms": True,
        "coordination_sequence_note": "first 10 terms only; this is evidence for the net sms, not an isomorphism proof",
        "trigonal_grid_scan": {"grid": "theta 0..90 deg step 1; u -0.6..0.6 step 0.01",
                               "feasible_states": tri["feasible"], "grid_connected_to_start": tri["connected"],
                               "densest_density": certlib.sig(tri["best_density"], 6),
                               "densest_theta_deg": certlib.sig(tri["best_theta_deg"], 4), "densest_u": certlib.sig(tri["best_u"], 3),
                               "bar_deviation_below_1e-12": True},
        "must_fail_controls_rejected": log,
        "inputs": certlib.inputs(HERE / "fold3d.py"),
    }
    certlib.write_json(HERE / "certificate.json", cert)
    print(f"fold-3d-cubic: OK (a(60) = 8c/3 exact; trigonal scan reaches {tri['best_density']:.4f} > {ex['density']:.4f})")


if __name__ == "__main__":
    main()
