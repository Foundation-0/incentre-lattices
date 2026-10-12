#!/usr/bin/env python3
"""Replay for results/fold-close-packing.

Runs BOTH checkers on the hcp fold path and on the fcc extension of it:
  1. autoresearch-fold/verify.py          (locked checker: bars, overlaps, continuity, jamming)
  2. results/fold-close-packing/chk_independent.py   (independent re-implementation + CNA)
Regenerates the fcc extension from the hcp path and compares it with the stored file,
runs must-fail controls, and writes certificate.json.

Dependencies: numpy, scipy (scipy only through autoresearch-fold/verify.py).
Slow target (about 3 minutes): the independent checker interpolates 10 sub-states between
every pair of states.   Usage: python verify.py
"""
import pathlib
import sys
import tempfile

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FOLD = ROOT / "autoresearch-fold"
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(HERE))
import certlib  # noqa: E402
import chk_independent as chk  # noqa: E402
import make_fcc_path  # noqa: E402
locked = make_fcc_path.verify  # autoresearch-fold/verify.py (hash-locked), loaded by path


def _check_locked_checker_pin():
    """The locked checker is loaded by path; refuse to run unless it matches its pin."""
    import hashlib
    recorded = (FOLD / "verify.sha256").read_text().split()[0]
    actual = hashlib.sha256((FOLD / "verify.py").read_bytes()).hexdigest()
    if actual != recorded:
        raise SystemExit("REFUSING TO RUN: autoresearch-fold/verify.py does not match verify.sha256")


_check_locked_checker_pin()

HCP = FOLD / "best" / "path.npz"
FCC = HERE / "path_fcc.npz"
N_HCP_STATES = 793


def locked_summary(rep):
    keep = ["valid", "states", "m", "contacts_per_ball_min", "contacts_per_ball_max", "contact_pairs",
            "flex_bars_only", "flex_bars_plus_contacts", "first_order_jammed"]
    out = {k: rep[k] for k in keep}
    out["reasons"] = list(rep["reasons"])
    out["final_density"] = certlib.sig(rep["final_density"], 10)
    out["contact_degree_histogram"] = rep["contact_degree_histogram"]
    if "end_lattice_cos_angles" in rep:
        out["end_lattice_cos_angles"] = [round(float(v), 9) + 0.0 for v in rep["end_lattice_cos_angles"]]  # absolute rounding: cos = 0 is roundoff noise that differs by platform
    return out


def independent_summary(rep, extra=None):
    out = {
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
        "end_next_distance_over_2R": certlib.sig(rep["end_next_distance_over_2R"], 8),
        "end_cna": rep["end_cna"],
        "accepted": chk.accepts(rep),
    }
    if extra:
        out.update(extra)
    return out


def main():
    # 1. regenerate the fcc extension and compare with the stored file
    m, Ls, Xs = make_fcc_path.build()
    stored = np.load(FCC)
    dev = max(abs(stored["Ls"] - Ls).max(), abs(stored["Xs"] - Xs).max())
    certlib.require(stored["Xs"].shape == Xs.shape and dev < 1e-9,
                    f"stored path_fcc.npz differs from the regenerated path (max dev {dev:.1e})")
    hcp = np.load(HCP)
    certlib.require(np.array_equal(stored["Xs"][:N_HCP_STATES], hcp["Xs"]), "fcc path does not start with the hcp path")

    # 2. both checkers on both paths
    cert = {"result": "fold-close-packing", "model": "m = 2 (48 balls, 120 bars)", "paths": {}}
    for name, path in (("hcp", HCP), ("fcc", FCC)):
        lk = locked.verify_file(str(path))
        certlib.require(lk["valid"], f"locked checker rejected the {name} path: {lk['reasons']}")
        ind = chk.check(path)
        certlib.require(chk.accepts(ind), f"independent checker rejected the {name} path")
        extra = {}
        if name == "fcc":
            d = np.load(path)
            dens = d["Xs"].shape[1] * 4 / 3 * np.pi * locked.R ** 3 / np.abs(np.linalg.det(d["Ls"]))
            extra["min_density_during_slides"] = certlib.sig(dens[N_HCP_STATES:].min(), 6)
            extra["slide_states"] = int(len(dens) - N_HCP_STATES)
            certlib.require(ind["end_cna"] == {"(4, 2, 1)": 576}, "fcc end state is not pure fcc by CNA")
        else:
            certlib.require(ind["end_cna"] == {"(4, 2, 1)": 288, "(4, 2, 2)": 288}, "hcp end state is not hcp by CNA")
        cert["paths"][name] = {"locked_checker": locked_summary(lk), "independent_checker": independent_summary(ind, extra)}
    cert["fcc_regeneration_matches_stored_file_to_1e-9"] = True
    cert["checker_agreement"] = all(
        cert["paths"][k]["locked_checker"]["valid"] and cert["paths"][k]["independent_checker"]["accepted"]
        for k in ("hcp", "fcc"))

    # 3. must-fail controls (quick independent pass without interpolation)
    log = []
    d = np.load(HCP)
    for cname, edit in (
        ("one dot displaced by 0.05 (bar lengths broken)",
         lambda L, X: (L, _set(X, (400, 0, 0), X[400, 0, 0] + 0.05))),
        ("lattice scaled by 0.97 in one state (overlaps and bar lengths broken)",
         lambda L, X: (_set(L, (500,), L[500] * 0.97), X)),
    ):
        L2, X2 = edit(d["Ls"].copy(), d["Xs"].copy())
        bad = locked.verify_arrays(int(d["m"]), L2, X2, want_rigidity=False)
        certlib.must_fail("locked checker: " + cname, bad["valid"], log)
        with tempfile.TemporaryDirectory() as td:
            p = pathlib.Path(td) / "bad.npz"
            np.savez(p, m=int(d["m"]), Ls=L2, Xs=X2)
            rep = chk.check(p, interp=False)
        certlib.must_fail("independent checker: " + cname, chk.accepts(rep, interp=False), log)
    # unlifted slide: slide the upper layers without lifting them over the bridge site
    f = np.load(FCC)
    Ls2, Xs2 = f["Ls"].copy(), f["Xs"].copy()
    t = N_HCP_STATES + 100
    Xs2[t] = Xs2[N_HCP_STATES - 1] + (Xs2[t] - Xs2[N_HCP_STATES - 1]) * 0.5
    bad = locked.verify_arrays(int(f["m"]), Ls2, Xs2, want_rigidity=False)
    certlib.must_fail("locked checker: state replaced by a half-way blend (bars not preserved)", bad["valid"], log)
    with tempfile.TemporaryDirectory() as td:
        p = pathlib.Path(td) / "bad.npz"
        np.savez(p, m=int(f["m"]), Ls=Ls2, Xs=Xs2)
        rep = chk.check(p, interp=False)
    certlib.must_fail("independent checker: state replaced by a half-way blend (bars not preserved)",
                      chk.accepts(rep, interp=False), log)
    cert["must_fail_controls_rejected"] = log

    cert["inputs"] = certlib.inputs(HCP, FCC, FOLD / "verify.py", HERE / "chk_independent.py", HERE / "make_fcc_path.py")
    certlib.write_json(HERE / "certificate.json", cert)
    print("fold-close-packing: OK (both checkers accept the hcp and fcc paths; controls rejected: %d)" % len(log))


def _set(a, idx, val):
    a = a.copy()
    a[idx] = val
    return a


if __name__ == "__main__":
    main()
