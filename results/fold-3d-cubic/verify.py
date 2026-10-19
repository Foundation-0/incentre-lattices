#!/usr/bin/env python3
"""Replay for results/fold-3d-cubic (SLOW target, about 1 minute).

The cubic-restricted [111] fold of the 3D contact framework (octahedra rotated alternately by
+-theta about [111], lattice kept cubic): exact values at 60 degrees, the cubic-path density at
30 degrees, the coordination sequence of the 60-degree contact net, a direct check that the
60-degree cubic state is not a jam (a trigonal shear of that state keeps every bar, overlaps
nothing and raises the density, at 100 sampled shears), and a grid scan of trigonal lattices
(every bar kept, density up to 0.6477).  Finite checks on samples and a grid; the grid
connectivity is NOT a proof of a continuous path (see results/fold-close-packing and
results/fold-m1 for sampled paths checked at step <= 0.01, and results/fold-continuous-fcc for
a continuous motion labelled PROVEN, internal review only).
Must-fail controls: claiming 60 degrees as a jam, the opposite shear claimed to raise the density,
a wrong lattice constant, the theta = 0 contact net claimed to have the sms sequence.
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
    d30 = 6 * 4 / 3 * np.pi * eng.R ** 3 / eng.a_cubic(np.pi / 6) ** 3
    certlib.require(abs(d30 - 0.26370) < 1e-5, "cubic path: density 0.2637 at 30 degrees")
    ex = eng.exact_60()
    certlib.require(ex["a60_minus_8c_over_3"] == 0 and ex["all_contacts_exact"] and ex["contacts"] == [6] * 6
                    and abs(ex["density_minus_closed_form"]) < 1e-25, "exact values at 60 degrees")
    certlib.require(all(abs(x - 1.3744) < 1e-4 for x in ex["next_over_2R"]), "next distance 1.3744 * 2R at 60 degrees")
    cs = eng.coordination_sequence(ex["a60"])
    certlib.require(cs == SMS, "coordination sequence")
    dus = np.linspace(-0.01, 0, 101)   # trigonal shear u = -K2(60) + du; du = 0 is the 60-degree cubic state
    shear = eng.trigonal_shear_60(dus)
    dens = [s[0] for s in shear]
    certlib.require(abs(dens[-1] - ex["density"]) < 1e-12, "shear at 60 degrees starts at the cubic state")
    certlib.require(all(s[2] < 1e-12 and s[1] >= 2 * eng.R - 1e-9 for s in shear), "shear at 60 degrees: bars kept, no overlap")
    certlib.require(all(d > ex["density"] + 1e-6 for d in dens[:-1]) and all(np.diff(dens) < 0),
                    "shear at 60 degrees: density strictly increasing as du decreases")
    tri = eng.trigonal_scan()
    certlib.require(tri["connected"] == tri["feasible"] and tri["best_density"] > ex["density"] + 0.1 and tri["bar_dev"] < 1e-12,
                    "trigonal scan")
    certlib.require(tri["feasible"] == 5008 and abs(tri["best_density"] - 0.647679) < 5e-6
                    and abs(tri["best_theta_deg"] - 31.0) < 1e-9 and abs(tri["best_u"] - 0.47) < 1e-9,
                    "trigonal scan: 5008 states, densest 0.647679 at 31 degrees, u = 0.47")

    log = []
    certlib.must_fail("60-degree cubic state claimed to be a jam (no bar-keeping, non-overlapping, denser trigonal shear)",
                      not any(s[2] < 1e-12 and s[1] >= 2 * eng.R - 1e-9 and s[0] > ex["density"] + 1e-9 for s in shear), log)
    certlib.must_fail("trigonal shear of the 60-degree cubic state the other way (du = +0.0001..+0.01) claimed to raise the density",
                      any(s[0] > ex["density"] for s in eng.trigonal_shear_60(-dus[:-1])), log)
    P, lab = eng.cell3(np.pi / 3, (ex["a60"] + 0.01) * np.eye(3))
    Dm = np.linalg.norm(P - P[lab.index(((0, 0, 0), 0))], axis=1)
    certlib.must_fail("lattice constant a(60) + 0.01 claimed to keep six contacts at 2R for a dot",
                      int(((Dm > 1e-9) & (Dm < 2 * eng.R + 1e-6)).sum()) == 6, log)
    certlib.must_fail("contact net at theta = 0 (five contacts per ball) claimed to have the sms coordination sequence",
                      eng.coordination_sequence(eng.a_cubic(0.0), t=0.0) == SMS, log)

    cert = {
        "result": "fold-3d-cubic",
        "model": "octahedra of 6 dots (c = 2 - sqrt2) rotated alternately by +-theta about [111]; lattice cubic; bars of length 2R kept",
        "cubic_path": {"first_nonbonded_contact_deg_grid_3000": certlib.sig(first, 5), "bar_length_deviation_below_1e-12": True},
        "cubic_path_density_at_30_deg": certlib.sig(d30, 7),
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
        "trigonal_motion_from_60_deg": {"family": "theta = 60 deg, u = -K2(60) + du (so b = du); du = 0 is the cubic state",
                                        "samples": "du = 0, -0.0001, ..., -0.01 (101 values)",
                                        "bar_deviation_below_1e-12": True, "min_distance_at_least_2R_minus_1e-9": True,
                                        "density_strictly_increasing_as_du_decreases": True,
                                        "du": -0.01, "density": certlib.sig(dens[0], 6)},
        "must_fail_controls_rejected": log,
        "inputs": certlib.inputs(HERE / "fold3d.py", HERE / "verify.py", HERE.parents[1] / "tools" / "certlib.py"),
    }
    certlib.write_json(HERE / "certificate.json", cert)
    print(f"fold-3d-cubic: OK (a(60) = 8c/3 exact; shear at 60 degrees reaches {dens[0]:.4f}; "
          f"trigonal scan reaches {tri['best_density']:.4f} > {ex['density']:.4f}; {len(log)} controls rejected)")


if __name__ == "__main__":
    main()
