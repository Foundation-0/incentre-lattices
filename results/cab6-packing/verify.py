#!/usr/bin/env python3
"""Replay for results/cab6-packing (maximal hard-sphere packing fraction of the CaB6 type).

Exact part (sympy): closed form, optimum position, Ca-B gap, fcc window endpoints.
Numerical part (numpy + scipy): brute-force maximisation over (h, scale), with the six B sites of one
cell checked against every Ca and B site in a 5x5x5 block of cells (Ca-Ca taken as 1), compared with
the closed form on a grid of size ratios.
Must-fail controls: wrong closed forms, "Ca touches B", a size ratio outside the window,
and the wrong Wyckoff convention for x are all rejected.
Writes certificate.json. Usage: python verify.py
"""
import contextlib
import io
import pathlib
import sys

import numpy as np
import sympy as sp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "tools"))
sys.path.insert(0, str(HERE))
import certlib  # noqa: E402

with contextlib.redirect_stdout(io.StringIO()):   # the engine prints its own report on import
    import cab6_packing as eng  # noqa: E402

s2 = sp.sqrt(2)
G0 = float(s2 - 1)
GRID = [0.1, 0.2, 0.3, 0.38, 0.40, 0.41, 0.4142, G0, 0.415, 0.42, 0.45, 0.5, 0.7, 1.0, 1.5]


def exact_part():
    rS = (s2 - 1) / 2
    phi = sp.pi / 6 * (30 * s2 - 41)
    out = {}
    certlib.require(sp.simplify(sp.Rational(4, 3) * sp.pi * (sp.Rational(1, 8) + 6 * rS ** 3) - phi) == 0, "closed form")
    h = sp.symbols("h", positive=True)
    hs = sp.solve(sp.Eq(h * s2, 1 - 2 * h), h)[0]
    x = sp.Rational(1, 2) - hs
    certlib.require(sp.simplify(hs - (1 - s2 / 2)) == 0 and sp.simplify(x - (s2 - 1) / 2) == 0, "optimum position")
    dCaB = sp.sqrt(x ** 2 + sp.Rational(1, 2))
    gap = sp.N(dCaB - sp.Rational(1, 2) - rS, 12)
    certlib.require(gap > 0, "Ca-B gap must be positive")
    lo = ((s2 - 1) / 6) ** sp.Rational(1, 3)
    hi = (4 + 7 * s2) ** sp.Rational(-1, 3)
    fcc = sp.pi / (3 * s2)
    certlib.require(sp.simplify(sp.pi / 6 * (1 + 6 * lo ** 3) - fcc) == 0, "left window end")
    certlib.require(sp.simplify(sp.pi / 6 * (s2 - 1) ** 3 * (6 + hi ** -3) - fcc) == 0, "right window end")
    nacl = sp.pi * (5 - 3 * s2) / 3
    certlib.require(sp.simplify(nacl - sp.pi / (3 * s2) * (1 + (s2 - 1) ** 3)) == 0, "NaCl value")
    out["phi_max_closed_form"] = "(pi/6)*(30*sqrt(2) - 41)"
    out["phi_max_value"] = certlib.sig(sp.N(phi, 15), 12)
    out["gamma_at_max"] = "sqrt(2) - 1"
    out["h_star"] = "1 - sqrt(2)/2"
    out["wyckoff_x_in_(x,1/2,1/2)_convention"] = "(sqrt(2) - 1)/2"
    out["ca_b_distance_over_a"] = certlib.sig(sp.N(dCaB, 12), 10)
    out["rCa_plus_rB_over_a"] = certlib.sig(sp.N(sp.Rational(1, 2) + rS, 12), 10)
    out["ca_b_gap_over_a"] = certlib.sig(gap, 8)
    out["ca_b_in_contact"] = False
    out["fcc_window_exact"] = {"lower": "((sqrt(2) - 1)/6)^(1/3)", "upper": "(4 + 7*sqrt(2))^(-1/3)",
                               "lower_value": certlib.sig(sp.N(lo, 12), 8), "upper_value": certlib.sig(sp.N(hi, 12), 8),
                               "phi_equals_fcc_at_both_ends": True}
    out["fcc_density"] = certlib.sig(sp.N(fcc, 12), 10)
    out["nacl_at_same_gamma"] = certlib.sig(sp.N(nacl, 12), 10)
    return out, float(phi), float(fcc), float(lo), float(hi)


def main():
    exact, phi_star, fcc, lo, hi = exact_part()
    # brute force against closed form
    rows, worst = [], 0.0
    for g in GRID:
        pb, hb = eng.phi_bf(g)
        pc = eng.phi_closed(g)
        worst = max(worst, abs(pb - pc))
        rows.append({"gamma": certlib.sig(g, 6), "closed_form": certlib.sig(pc, 8), "beats_fcc": bool(pc > fcc)})
        certlib.require(abs(pb - pc) < 1e-9, f"brute force disagrees with the closed form at gamma={g}")
    certlib.require(all(r["beats_fcc"] == (lo < r["gamma"] < hi) for r in rows if abs(r["gamma"] - lo) > 1e-5 and abs(r["gamma"] - hi) > 1e-5),
                    "beats-fcc pattern does not match the window")
    gg = np.linspace(0.05, 1.5, 3000)
    ph = np.array([eng.phi_closed(g) for g in gg])
    certlib.require(abs(ph.max() - phi_star) < 1e-3 and abs(gg[ph.argmax()] - G0) < 1e-3, "grid maximum")

    # must-fail controls
    log = []
    pb, _ = eng.phi_bf(G0)
    certlib.must_fail("wrong constant: (pi/6)(30 sqrt2 - 40) as the maximum", abs(np.pi / 6 * (30 * 2 ** 0.5 - 40) - pb) < 1e-9, log)
    certlib.must_fail("wrong branch: phi = (pi/6)(1+6 g^3) used for g > sqrt2-1",
                      abs(np.pi / 6 * (1 + 6 * 0.45 ** 3) - eng.phi_bf(0.45)[0]) < 1e-9, log)
    certlib.must_fail("gamma = 0.40 claimed to beat fcc", eng.phi_closed(0.40) > fcc, log)
    certlib.must_fail("gamma = 0.42 claimed to beat fcc", eng.phi_closed(0.42) > fcc, log)
    rS = (G0) / 2
    # nearest Ca-B distance at h* = 1 - sqrt2/2, measured in the 5^3 block: sqrt(x*^2 + 1/2) with x* = (sqrt2 - 1)/2
    dCaB = float(np.linalg.norm(eng.sites(1 - 2 ** 0.5 / 2)[:, None] - eng.shifts[None], axis=2).min())
    certlib.must_fail("Ca-B claimed to touch (distance equals rCa + rB)", abs(dCaB - (0.5 + rS)) < 1e-6, log)
    x_wrong = (2 ** 0.5 - 1) / 2   # used as the offset h from the octahedron centre instead of as Wyckoff x
    certlib.must_fail("x = (sqrt2-1)/2 read as the offset from the centre: intra and inter B-B distances differ",
                      abs(x_wrong * 2 ** 0.5 - (1 - 2 * x_wrong)) < 1e-9, log)

    cert = {
        "result": "cab6-packing",
        "model": "Ca radius rL on 1a (0,0,0); B radius rS on 6f (x,1/2,1/2) and permutations; x and lattice parameter a free; gamma = rS/rL",
        "exact": exact,
        "brute_force_vs_closed_form": {"cell_block": "5x5x5; the six B sites of one cell against every Ca and B site in the block; Ca-Ca taken as 1",
                                       "gamma_grid": rows, "max_abs_difference_below_1e-9": True},
        "must_fail_controls_rejected": log,
        "inputs": certlib.inputs(HERE / "cab6_packing.py", HERE / "verify.py", HERE.parents[1] / "tools" / "certlib.py"),
    }
    certlib.write_json(HERE / "certificate.json", cert)
    print(f"cab6-packing: OK (max |brute - closed| < 1e-9 on {len(GRID)} ratios; {len(log)} controls rejected)")


if __name__ == "__main__":
    main()
