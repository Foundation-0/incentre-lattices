#!/usr/bin/env python3
"""Replay for results/fold-2d (KNOWN linkage; incentre numbers are corollaries).

Exact (sympy): gap(theta) = 2 sqrt2 c cos(theta + pi/4), first contact at pi/4, L(pi/4) = 4R,
strain -R^2, density pi(sqrt2-1)^2 -> pi/4.  Numerical (numpy): lock angle and nearest
distances for dots placed at +-c e_i for several c; only c = 2 - sqrt2 ends in a square lattice
with a single contact distance (checked on a scan, not proved for all c).
Must-fail controls: wrong gap formula, wrong lock angle, and a non-incentre c claimed to end in a
single-distance square lattice are all rejected.  Writes certificate.json.  Usage: python verify.py
"""
import pathlib
import sys

import numpy as np
import sympy as sp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "tools"))
sys.path.insert(0, str(HERE))
import certlib  # noqa: E402
import fold2d as eng  # noqa: E402

SCAN = [0.45, 0.5, 0.55, eng.c, 0.62, 0.7]


def single_distance_square(nn):
    """Four nearest neighbours all at 2r (a square lattice with one contact distance)."""
    return bool(np.all(np.abs(nn[:4] - 1) < 1e-3))


def main():
    f = eng.exact_facts()
    claim = 2 * sp.sqrt(2) * eng.c_ * sp.cos(eng.th + sp.pi / 4)
    certlib.require(sp.simplify(f["gap"] - claim) == 0, "gap(theta) formula")
    certlib.require(f["gap_at_pi_over_4"] == 0 and f["L_minus_4R"] == 0 and f["strain_plus_R2"] == 0, "exact values at pi/4")
    certlib.require(sp.simplify(f["density_end"] - sp.pi / 4) == 0, "end density pi/4")
    # gap is positive for theta in [0, pi/4): cos(theta + pi/4) > 0 there
    scan = []
    for cc in SCAN:
        t, nn = eng.lock2(cc)
        scan.append({"c": certlib.sig(cc, 6), "lock_angle_deg_0.01_grid": round(t, 2),
                     "nearest_over_2r": [certlib.sig(x, 4) for x in nn],
                     "four_contacts_at_2r": single_distance_square(nn)})
    inc = [s for s in scan if abs(s["c"] - eng.c) < 1e-5][0]
    certlib.require(inc["lock_angle_deg_0.01_grid"] == 45.0 and inc["four_contacts_at_2r"], "incentre lock")
    certlib.require(all(not s["four_contacts_at_2r"] for s in scan if s is not inc), "other c must not be a single-distance square lattice")

    log = []
    wrong = 2 * sp.sqrt(2) * eng.c_ * sp.cos(eng.th - sp.pi / 4)
    certlib.must_fail("wrong gap formula cos(theta - pi/4)", sp.simplify(f["gap"] - wrong) == 0, log)
    t05, _ = eng.lock2(0.5)
    certlib.must_fail("lock at 45 degrees claimed for c = 0.5", abs(t05 - 45.0) < 0.05, log)
    _, nn62 = eng.lock2(0.62)
    certlib.must_fail("c = 0.62 claimed to end in a single-distance square lattice", single_distance_square(nn62), log)
    certlib.must_fail("final density claimed to be pi/3 instead of pi/4", sp.simplify(f["density_end"] - sp.pi / 3) == 0, log)

    cert = {
        "result": "fold-2d",
        "model": "squares of 4 dots at +-c e_i in cells of side 2; cross bars along the axes; balls radius R = sqrt2 - 1 (c = 2 - sqrt2)",
        "exact": {
            "L_of_theta": "2R + 2c cos(theta)",
            "gap_of_theta": "2 sqrt(2) c cos(theta + pi/4)",
            "gap_at_pi_over_4": 0,
            "L_at_pi_over_4_minus_4R": 0,
            "linear_strain_at_pi_over_4": "-(sqrt(2)-1)^2",
            "density_start": certlib.sig(sp.N(f["density_start"], 15), 10),
            "density_end": "pi/4",
            "density_end_value": certlib.sig(np.pi / 4, 10),
        },
        "scan_over_c": scan,
        "must_fail_controls_rejected": log,
        "inputs": certlib.inputs(HERE / "fold2d.py"),
    }
    certlib.write_json(HERE / "certificate.json", cert)
    print(f"fold-2d: OK ({len(scan)} values of c scanned; {len(log)} controls rejected)")


if __name__ == "__main__":
    main()
