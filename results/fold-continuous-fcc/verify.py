#!/usr/bin/env python3
"""Replay the exact finite certificate in the continuous fcc proof (standard library only)."""
import pathlib
import sys
from fractions import Fraction as F

import curve
import independent

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import certlib  # noqa: E402


def main():
    # Near-cancelling convergents detect incorrect floating sign decisions.
    for x, expected in [((0, 0), 0), ((1, 1), 1), ((-1, -1), -1),
                        ((1, -1), -1), ((-1, 1), 1),
                        ((665857, -470832), 1), ((-665857, 470832), -1),
                        ((886731088897, -627013566048), 1)]:
        certlib.require(curve.sign(tuple(map(F, x))) == expected, f"Q(sqrt2) sign {x}")

    proof = curve.exact_nonoverlap()
    second = independent.check()
    certlib.require(second["distance_polynomials"] == proof["distance_polynomials"] and
                    second["zero_coefficients"] == proof["zero_coefficients"] and
                    second["positive_coefficients"] + second["zero_coefficients"] == proof["bernstein_coefficients"] and
                    second["identically_zero_polynomials"] == proof["identically_zero_polynomials"],
                    "the two exact implementations disagree")
    certlib.require(proof["integer_shifts"] == 171 and proof["motif_differences"] == 19,
                    "enumeration coverage")
    certlib.require(proof["distance_polynomials"] == 3248 and
                    proof["identically_zero_polynomials"] == 18, "distance counts")
    # Each bridge has |a e_i+t C_i|^2 = a^2+2 t^2 = 2.
    certlib.require(all(curve.C[i][i] == 0 and sum(x*x for x in curve.C[i]) == 2
                        for i in range(3)), "bridge identity")
    # Bernstein conversion: reconstruct power coefficients independently by expansion.
    for cs in [curve.coefficients((1, -2, 1), (1, 0, -1)),
               [(F(i-2), F(3-i)) for i in range(5)]]:
        bs = curve.bernstein(cs)
        from math import comb
        rebuilt = [curve.ZERO for _ in range(5)]
        for k, b in enumerate(bs):
            for j in range(k, 5):
                rebuilt[j] = curve.add(rebuilt[j], curve.scale(
                    b, comb(4, k)*comb(4-k, j-k)*(-1)**(j-k)))
        certlib.require(rebuilt == cs, "Bernstein basis conversion")

    end = curve.endpoint()
    certlib.require(end["determinant"] == 12 and end["lattice_rows_have_even_sum"] and
                    end["motif_points_have_odd_sum"] and end["distinct_cosets"] and
                    end["contacts_per_ball"] == [12]*6, "fcc parity/index endpoint")
    controls = []
    # Same row lengths and starting lattice, but a sign change creates endpoint overlaps.
    wrong = ((0, 1, 1), (1, 0, 1), (1, 1, 0))
    accepted = True
    try:
        curve.exact_nonoverlap(wrong)
    except ValueError:
        accepted = False
    certlib.must_fail("sign-changed motion with overlapping endpoint", accepted, controls)
    interior_negative = [(F(3, 4), F(0)), (F(-4), F(0)), (F(4), F(0)), curve.ZERO, curve.ZERO]
    certlib.must_fail("positive endpoints with a negative interior",
                     all(curve.sign(b) >= 0 for b in curve.bernstein(interior_negative)), controls)
    certlib.must_fail("negative exact coefficient rounded toward zero",
                     curve.sign((F(-665857), F(470832))) >= 0, controls)
    cert = {"result": "fold-continuous-fcc",
            "evidence_label": "PROVEN by explicit curve and exact finite Bernstein certificate; internal review only",
            "arithmetic": "rational pairs a+b*sqrt(2); no floating arithmetic in this replay",
            "exact_nonoverlap": proof, "independent_exact_check": second, "endpoint": end,
            "bar_lengths_exact": True, "density_strictly_increasing": "written derivative argument, PROOF.md",
            "infinite_image_bound": "written symmetric-part bound: |n|^2 >= 12 strictly safe",
            "end_density": "pi/sqrt(18)", "must_fail_controls_rejected": controls,
            "inputs": certlib.inputs(HERE / "curve.py", HERE / "independent.py", HERE / "verify.py", HERE / "PROOF.md",
                                     ROOT / "tools" / "certlib.py")}
    certlib.write_json(HERE / "certificate.json", cert)
    print("fold-continuous-fcc: OK (3248 quartics, 16240 exact Bernstein coefficients; 3 controls rejected)")


if __name__ == "__main__":
    main()
