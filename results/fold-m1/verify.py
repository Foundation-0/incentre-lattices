#!/usr/bin/env python3
"""Replay for results/fold-m1: the fold reaches close packing in the 1 x 1 x 1 cell (6 balls, 15 bars).

For each stored path (path_hcp_m1.npz: free density descent, no steering; path_fcc_m1.npz: steered to
an explicit fcc end state, density never decreasing):
  1. the locked checker autoresearch-fold/verify.py (bars, overlaps, midpoint, step <= 0.01, jamming);
  2. the independent checker results/fold-close-packing/chk_independent.py (with 10 sub-states per step);
  3. the exact end-state check results/fold-end-exact/exact_end.py (integer arithmetic);
  4. for the fcc path, det L strictly decreases at every step (density increases at every sample);
then must-fail controls, and certificate.json.

Dependencies: numpy, scipy (scipy only through autoresearch-fold/verify.py). About 20 seconds.
Usage: python verify.py
"""
import importlib.util
import pathlib
import sys
import tempfile

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FOLD = ROOT / "autoresearch-fold"
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "results" / "fold-close-packing"))
sys.path.insert(0, str(ROOT / "results" / "fold-end-exact"))
import certlib  # noqa: E402
import chk_independent as chk  # noqa: E402
import exact_end  # noqa: E402

_spec = importlib.util.spec_from_file_location("fold_verify", FOLD / "verify.py")
locked = importlib.util.module_from_spec(_spec)  # the hash-locked checker, loaded by path
_spec.loader.exec_module(locked)

PATHS = {"hcp": HERE / "path_hcp_m1.npz", "fcc": HERE / "path_fcc_m1.npz"}
CNA = {"hcp": {"(4, 2, 1)": 36, "(4, 2, 2)": 36}, "fcc": {"(4, 2, 1)": 72}}


def locked_summary(rep):
    keep = ["valid", "states", "m", "contacts_per_ball_min", "contacts_per_ball_max", "contact_pairs",
            "flex_bars_only", "flex_bars_plus_contacts", "first_order_jammed"]
    out = {k: rep[k] for k in keep}
    out["reasons"] = list(rep["reasons"])
    out["final_density"] = certlib.sig(rep["final_density"], 10)
    out["contact_degree_histogram"] = rep["contact_degree_histogram"]
    return out


def independent_summary(rep):
    return {
        "states": rep["states"], "bars": rep["bars"], "dots": rep["dots"],
        "bar_error_below_1e-9_relative": rep["max_bar_rel_err"] < 1e-9,
        "overlap_violations": rep["overlap_violations"],
        "min_nonbonded_over_2R_at_least_1_minus_1e-9": rep["min_nonbonded_over_2R"] > 1 - 1e-9,
        "max_dot_step": certlib.sig(rep["max_dot_step"], 4),
        "max_lattice_step": certlib.sig(rep["max_lattice_step"], 4),
        "interp_min_nonbonded_over_2R": certlib.sig(rep["interp_min_nonbonded_over_2R"], 6),
        "interp_max_bar_dev_over_2R": certlib.sig(rep["interp_max_bar_dev_over_2R"], 4),
        "start_density": certlib.sig(rep["start_density"], 10),
        "end_density": certlib.sig(rep["end_density"], 10),
        "end_density_equals_kepler_to_1e-12": abs(rep["end_density_minus_kepler"]) < 1e-12,
        "end_contacts_per_ball": rep["end_contacts_per_ball"],
        "end_cna": rep["end_cna"],
        "accepted": chk.accepts(rep),
    }


def exact_summary(rep):
    keep = ["unit", "lattice_rows", "det_L", "bars_exact", "no_overlap_exact", "contacts_per_ball",
            "density_is_pi_over_sqrt18_exactly", "cna_pairs", "structure", "exact_ok"]
    return {k: rep[k] for k in keep}


def main():
    cert = {"result": "fold-m1", "model": "m = 1 (6 balls, 15 bars)", "paths": {}}
    for name, path in PATHS.items():
        lk = locked.verify_file(str(path))
        certlib.require(lk["valid"] and lk["m"] == 1, f"locked checker rejected the {name} path: {lk['reasons']}")
        ind = chk.check(path)
        certlib.require(chk.accepts(ind), f"independent checker rejected the {name} path")
        certlib.require(ind["end_cna"] == CNA[name], f"{name}: end state CNA is {ind['end_cna']}")
        ex = exact_end.check(path)
        certlib.require(ex["exact_ok"] and ex["structure"] == name, f"{name}: exact end-state check failed")
        d = np.load(path)
        det = np.linalg.det(d["Ls"])
        entry = {"locked_checker": locked_summary(lk), "independent_checker": independent_summary(ind),
                 "exact_end_state": exact_summary(ex),
                 "det_L_strictly_decreasing_at_every_step": bool(np.all(np.diff(det) < 0))}
        if name == "fcc":
            certlib.require(entry["det_L_strictly_decreasing_at_every_step"], "fcc path: density decreases somewhere")
        cert["paths"][name] = entry

    # must-fail controls
    log = []
    d = np.load(PATHS["fcc"])
    for cname, edit in (
        ("one dot displaced by 0.05 (bar lengths broken)", lambda L, X: (L, _set(X, (100, 3, 0), X[100, 3, 0] + 0.05))),
        ("lattice scaled by 0.97 in one state (bars and overlaps broken)", lambda L, X: (_set(L, (150,), L[150] * 0.97), X)),
        ("states 100 to 199 deleted (step too large)", lambda L, X: (np.delete(L, range(100, 200), 0), np.delete(X, range(100, 200), 0))),
    ):
        L2, X2 = edit(d["Ls"].copy(), d["Xs"].copy())
        bad = locked.verify_arrays(1, L2, X2, want_rigidity=False)
        certlib.must_fail("locked checker: " + cname, bad["valid"], log)
        with tempfile.TemporaryDirectory() as td:
            p = pathlib.Path(td) / "bad.npz"
            np.savez(p, m=1, Ls=L2, Xs=X2)
            rep = chk.check(p, interp=False)
        certlib.must_fail("independent checker: " + cname, chk.accepts(rep, interp=False), log)
    h = np.load(PATHS["hcp"])
    certlib.must_fail("monotone test on the hcp path (its density dips)", bool(np.all(np.diff(np.linalg.det(h["Ls"])) < 0)), log)
    cert["must_fail_controls_rejected"] = log

    cert["inputs"] = certlib.inputs(PATHS["hcp"], PATHS["fcc"], FOLD / "verify.py",
                                    ROOT / "results" / "fold-close-packing" / "chk_independent.py",
                                    ROOT / "results" / "fold-end-exact" / "exact_end.py",
                                    HERE / "search_hcp_m1.py", HERE / "search_fcc_m1.py")
    certlib.write_json(HERE / "certificate.json", cert)
    print("fold-m1: OK (both checkers and the exact check accept the m = 1 hcp and fcc paths; controls rejected: %d)" % len(log))


def _set(a, idx, val):
    a = a.copy()
    a[idx] = val
    return a


if __name__ == "__main__":
    main()
