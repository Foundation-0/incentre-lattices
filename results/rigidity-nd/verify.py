#!/usr/bin/env python3
"""Replay for results/rigidity-nd (contact framework of the truncated n-cubic honeycomb, n = 2..8).

Dependencies: numpy, sympy.  Runtime: a few seconds.

What is replayed (all finite checks, for n = 2..8):
  * 2n^2 degrees of freedom and 2n^2 - n bars per cell, all of one length;
  * zero-mode and self-stress counts of the periodic rigidity matrix at random generic q, at q = 0
    and on the planes q_1 = ... = q_k = 0, against n(n-1)/2, n(n-3)/2 (and the q = 0 and plane
    counts of the theorem) for n >= 3, and against fixed n = 2 values (n = 2 is outside the theorem);
  * the cross-polytope alone is infinitesimally rigid with n(n-3)/2 self-stresses for n >= 3;
  * the explicit stress basis spans the generic self-stresses; the explicit q = 0 stresses;
    the generic zero modes are exactly the cluster rotations;
  * elastic constants with full internal relaxation: C11, C12, vanishing shear, nu_100 = 1/(n^2-n-1),
    for the incentre c = 2 - sqrt2 and for c = 0.3, 0.75;
  * the three steps of the written proof, symbolically (proof_steps.py).
Must-fail controls. The engine is rerun on a broken input and the result must be rejected: the elastic tensor
without internal relaxation (P = I) fails the C11 check at n = 4; moving one dot off its axis changes the q = 0
and the plane q_1 = 0 counts at n = 4; deleting one cross bar changes the generic zero-mode count at n = 5.
Outside the theorem's range: the n >= 3 generic zero-mode count n(n-1)/2 is rejected at n = 2.
Writes certificate.json.  Usage: python verify.py
"""
import contextlib
import io
import pathlib
import sys
from fractions import Fraction

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "tools"))
sys.path.insert(0, str(HERE))
import certlib  # noqa: E402
import proof_steps  # noqa: E402
import rigidity_nd as eng  # noqa: E402

NS = range(2, 9)
CS = (2 - np.sqrt(2), 0.3, 0.75)
TOL_ELASTIC = 1e-12  # absolute tolerance on C11, C12, shear, coupling, cubic form


def predicted(n):
    return (n * (n - 1) // 2, n * (n - 3) // 2)


def q0_ok(n, q0):
    return q0 == (n * (n + 1) // 2, n * (n - 1) // 2)


def plane_ok(n, k, z, s):
    """On the plane q_1 = ... = q_k = 0 there are k more zero modes and k more self-stresses than at generic q."""
    return z == predicted(n)[0] + k and s == predicted(n)[1] + k


def c11_exact(n):
    return Fraction(4 * (n - 1) ** 2, n ** 2 * 2 ** n)


def main():
    rows, elastic = [], []
    for n in NS:
        o = eng.check_n(n)
        pred = predicted(n)
        g = [tuple(int(v) for v in z) for z in o["generic"]]
        q0 = tuple(int(v) for v in o["q0"])
        planes = [(int(k), int(z), int(s)) for k, z, s in o["planes"]]
        cp = tuple(int(v) for v in o["crosspolytope"])
        sb_dim, sb_resid, sb_span = o["stress_basis"]
        row = {"n": n, "dofs_per_cell": int(o["dofs"]), "bars_per_cell": int(o["bars"]),
               "all_bar_lengths_equal": bool(abs(o["lengths"][1] - o["lengths"][0]) < 1e-12),
               "generic_q_zero_modes_selfstresses": g[0], "generic_q_all_three_samples_equal": len(set(g)) == 1,
               "q0_zero_modes_selfstresses": q0, "planes_k_zeromodes_selfstresses": planes,
               "crosspolytope_bars_rank_zeromodes_selfstresses": cp}
        certlib.require(row["dofs_per_cell"] == 2 * n * n and row["bars_per_cell"] == 2 * n * n - n,
                        f"n={n}: {row['dofs_per_cell']} dofs and {row['bars_per_cell']} bars per cell")
        certlib.require(row["all_bar_lengths_equal"], f"n={n}: bar lengths {o['lengths']}")
        if n >= 3:
            certlib.require(len(set(g)) == 1 and g[0] == pred, f"n={n}: generic counts {g} != {pred}")
            certlib.require(q0_ok(n, q0), f"n={n}: q=0 counts {q0}")
            certlib.require(all(plane_ok(n, k, z, s) for k, z, s in planes), f"n={n}: plane counts {planes}")
            certlib.require(cp == (2 * n * (n - 1), 2 * n * n - n * (n + 1) // 2, n * (n + 1) // 2, n * (n - 3) // 2),
                            f"n={n}: cross-polytope is not rigid with n(n-3)/2 stresses (bars, rank, zero modes, stresses {cp})")
            certlib.require(sb_dim == n * (n - 3) // 2 and sb_resid < 1e-9 and bool(sb_span), f"n={n}: stress basis")
            certlib.require(o["q0_stress"][0] < 1e-9 and o["q0_stress"][1] < 1e-9, f"n={n}: q=0 stress formula")
            certlib.require(o["modes_are_rotations"][0] < 1e-9 and bool(o["modes_are_rotations"][1]), f"n={n}: zero modes are rotations")
            row.update({"explicit_stress_basis_dimension": int(sb_dim), "explicit_stress_basis_spans_all_stresses": True,
                        "q0_stress_formula_tau_i_equals_sqrt2_sum_w_ij": True, "generic_zero_modes_equal_cluster_rotations": True})
        else:
            certlib.require(all(z == (2, 0) for z in g) and q0 == (3, 1), "n=2 counts (2,0) generic and (3,1) at q=0")
            certlib.require(planes == [(0, 2, 0), (1, 2, 0), (2, 3, 1)], f"n=2: plane counts {planes}")
            certlib.require(cp == (4, 4, 4, 0), f"n=2: cross-polytope (bars, rank, zero modes, self-stresses) {cp}")
            row["note"] = "n = 2 is outside the theorem: the K_2 kappa-block has one internal flex"
        rows.append(row)
    for n in NS:
        pn = Fraction(1, n * n - n - 1)
        c11 = c11_exact(n)
        c12 = Fraction(4, n ** 2 * 2 ** n)
        for c in CS:
            C11, C12, sh, cp_, nu, Cn = eng.elastic(n, c)
            # The exact values are proven in PROOF.md; this replay checks that the numerical relaxation
            # (rigidity_nd.elastic, threshold-free up to a checked spectral gap) agrees with them to 1e-12.
            # Residuals are printed on failure.
            tol = TOL_ELASTIC
            cubic_exp = float(c12) * np.ones((n, n)) + float(c11 - c12) * np.eye(n)
            resid = {"C11": abs(C11 - float(c11)), "C12": abs(C12 - float(c12)), "shear": sh, "coupling": cp_,
                     "cubic_form": float(np.abs(Cn - cubic_exp).max())}
            ok = all(v < tol for v in resid.values()) and abs(nu - float(pn)) < 1e-10
            detail = ", ".join(f"{k}={v:.3e}" for k, v in resid.items())
            certlib.require(ok, f"n={n}, c={c:.4f}: elastic constants (tol {tol:.3e}; {detail}; "
                                f"nu100 residual {abs(nu - float(pn)):.3e} vs 1e-10)")
        elastic.append({"n": n, "C11_exact": str(c11), "C12_exact": str(c12), "C11_over_C12": str(c11 / c12),
                        "nu100_exact": str(pn), "nu100_value": certlib.sig(float(pn), 10),
                        "matches_numerics_for_c_values": [certlib.sig(c, 6) for c in CS],
                        "numerical_agreement": "C11, C12 to 1e-12; shear and normal-shear coupling below 1e-12; nu100 to 1e-10",
                        "shear_moduli_zero": True})

    with contextlib.redirect_stdout(io.StringIO()):
        steps_ok = proof_steps.run()
    certlib.require(steps_ok, "symbolic proof steps")

    # must-fail controls: the engine recomputes from a broken input (the first: outside the theorem's range)
    # and the result must be rejected
    log = []
    n2 = eng.check_n(2)
    certlib.must_fail("n >= 3 generic zero-mode count n(n-1)/2 = 1 claimed at n = 2",
                      all(int(z[0]) == predicted(2)[0] for z in n2["generic"]), log)
    C11 = eng.elastic(4, 2 - np.sqrt(2), relax=False)[0]
    certlib.must_fail("C11 = 4(n-1)^2/(n^2 2^n) without internal relaxation (P = I) at n = 4",
                      abs(C11 - float(c11_exact(4))) < TOL_ELASTIC, log)
    n = 4
    shift = np.zeros((2 * n, n))
    shift[0, 1] = 0.05  # dot p(1,+) = c e_1 moved 0.05 along e_2, off its axis
    moved = eng.check_n(n, shift=shift)
    certlib.must_fail("q = 0 counts unchanged after moving dot p(1,+) by 0.05 along e_2 (n = 4)",
                      q0_ok(n, tuple(int(v) for v in moved["q0"])), log)
    z1, s1 = next((int(z), int(s)) for k, z, s in moved["planes"] if k == 1)
    certlib.must_fail("plane q_1 = 0 counts unchanged after moving dot p(1,+) by 0.05 along e_2 (n = 4)",
                      plane_ok(n, 1, z1, s1), log)
    n = 5
    Rq = eng.rigidity(n, 2 - np.sqrt(2), np.array([0.7, -1.3, 2.1, 0.4, -2.2]))
    V = 2 * n
    full = n * V - eng.rank(Rq)
    removed = n * V - eng.rank(np.delete(Rq, Rq.shape[0] - 1, axis=0))
    certlib.must_fail("generic zero-mode count unchanged after deleting a cross bar (n = 5; cross bars carry no generic-q stress)", removed == full, log)

    cert = {
        "result": "rigidity-nd",
        "claim": {"scope": "n >= 3; n = 2 values in counts_by_n and elastic_by_n",
                  "generic_q_zero_modes_per_cell": "n(n-1)/2", "generic_q_self_stresses_per_cell": "n(n-3)/2",
                  "q0_zero_modes": "n(n+1)/2", "q0_self_stresses": "n(n-1)/2",
                  "C11": "4(n-1)^2/(n^2 2^n)", "C12": "4/(n^2 2^n)", "C44": 0, "nu100": "1/(n^2-n-1)"},
        "evidence": {"proven": "all n >= 3 (written proof, PROOF.md)",
                     "verified": "counts n = 3..8; n = 2 counts separately; elastic n = 2..8 (finite numerical checks)"},
        "counts_by_n": rows,
        "elastic_by_n": elastic,
        "symbolic_proof_steps_n2_to_n8": True,
        "numerical_tolerances": {"rank_threshold": eng.TOL, "elastic": 1e-12},
        "must_fail_controls_rejected": log,
        "inputs": certlib.inputs(HERE / "rigidity_nd.py", HERE / "proof_steps.py", HERE / "verify.py",
                                 HERE.parents[1] / "tools" / "certlib.py"),
    }
    certlib.write_json(HERE / "certificate.json", cert)
    print(f"rigidity-nd: OK (n = 2..8; {len(log)} controls rejected)")


if __name__ == "__main__":
    main()
