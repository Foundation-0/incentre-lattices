#!/usr/bin/env python3
"""Replay for results/fold-m1: the fold reaches close packing in the 1 x 1 x 1 cell (6 balls, 15 bars).

For each stored path (path_hcp_m1.npz: free density descent, no steering; path_fcc_m1.npz: steered to
an explicit fcc end state, density never decreasing):
  1. the locked checker autoresearch-fold/verify.py (bars, overlaps, midpoint, step <= 0.01, jamming);
     its end state must have 12 contacts per ball and be first-order jammed (no flex with the contacts
     as bars, no first-order density ascent);
  2. the independent checker results/fold-close-packing/chk_independent.py (with 9 interior sub-states
     u = 0.1, ..., 0.9 per step);
  3. the exact end-state check results/fold-end-exact/exact_end.py (integer arithmetic); its lattice rows
     must be the end-state table of PROOF.md;
  4. for the fcc path, det L strictly decreases at every step (density increases at every sample);
then must-fail controls, and certificate.json. Both checkers are checked against their sha256 pins
(autoresearch-fold/verify.sha256, results/fold-close-packing/checkers.sha256) before they are imported.

Dependencies: numpy, scipy (scipy only through autoresearch-fold/verify.py). About 10 seconds.
Usage: python verify.py
"""
import hashlib
import importlib.util
import pathlib
import sys
import tempfile

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FOLD = ROOT / "autoresearch-fold"
FCP = ROOT / "results" / "fold-close-packing"


def _check_pins(folder, sums, needed):
    """Refuse to run unless every file pinned in folder/sums matches its sha256. Standard library only,
    and run before either checker is imported, so a modified checker cannot patch hashlib first."""
    rel, pinned = folder.relative_to(ROOT).as_posix(), set()
    for line in (folder / sums).read_text().splitlines():
        if not line.strip():
            continue
        digest, _, name = line.partition("  ")
        p = folder / name
        if not (p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest() == digest):
            raise SystemExit(f"REFUSING TO RUN: {rel}/{name} does not match {rel}/{sums}")
        pinned.add(name)
    if not needed <= pinned:
        raise SystemExit(f"REFUSING TO RUN: {rel}/{sums} does not pin {sorted(needed - pinned)}")


_check_pins(FOLD, "verify.sha256", {"verify.py"})
_check_pins(FCP, "checkers.sha256", {"chk_independent.py"})

sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(FCP))
sys.path.insert(0, str(ROOT / "results" / "fold-end-exact"))
import certlib  # noqa: E402
import chk_independent as chk  # noqa: E402
import exact_end  # noqa: E402

_spec = importlib.util.spec_from_file_location("fold_verify", FOLD / "verify.py")
locked = importlib.util.module_from_spec(_spec)  # the hash-locked checker, loaded by path
_spec.loader.exec_module(locked)

PATHS = {"hcp": HERE / "path_hcp_m1.npz", "fcc": HERE / "path_fcc_m1.npz"}
CNA = {"hcp": {"(4, 2, 1)": 36, "(4, 2, 2)": 36}, "fcc": {"(4, 2, 1)": 72}}
# the end-state table of PROOF.md: lattice rows in the cube frame of the octahedron, in the given unit
END = {"fcc": ([[2, -1, -1], [1, 2, -1], [1, -1, 2]], "K/1, K = 2 - sqrt2"),
       "hcp": ([[7, -1, 4], [-4, 7, -1], [-3, 3, 6]], "K/3, K = 2 - sqrt2")}


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
        certlib.require(lk.get("first_order_jammed") is True and lk.get("flex_bars_plus_contacts") == 0,
                        f"{name}: end state is not first-order jammed")
        certlib.require(lk.get("contacts_per_ball_min") == lk.get("contacts_per_ball_max") == 12,
                        f"{name}: locked checker does not count 12 contacts per ball at the end state")
        ind = chk.check(path)
        certlib.require(chk.accepts(ind), f"independent checker rejected the {name} path")
        certlib.require(ind["end_cna"] == CNA[name], f"{name}: end state CNA is {ind['end_cna']}")
        ex = exact_end.check(path)
        certlib.require(ex["exact_ok"] and ex["structure"] == name, f"{name}: exact end-state check failed")
        certlib.require((ex["lattice_rows"], ex["unit"]) == END[name], f"{name}: end lattice is not the PROOF.md table")
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
        rep = _independent(L2, X2)
        certlib.must_fail("independent checker: " + cname, chk.accepts(rep, interp=False), log)
    h = np.load(PATHS["hcp"])
    certlib.must_fail("monotone test on the hcp path (its density dips)", bool(np.all(np.diff(np.linalg.det(h["Ls"])) < 0)), log)

    # Controls that each break one property; the clause that rejects them is asserted, not only the rejection.
    # Overlap with bars intact and steps <= 0.01: the curve of results/fold-continuous-fcc, L = K((2 + a) I + t C)
    # with the dots fixed at the start state, run past its fcc end point s = 1 to s = 1.1. Its bars stay exactly
    # 2R for every real s (a^2 + 2 t^2 = 2), but past s = 1 balls overlap.
    s = np.linspace(0.0, 1.1, 441)
    t, a = 2 * s / (1 + s * s), np.sqrt(2) * (1 - s * s) / (1 + s * s)
    C = np.array([[0, -1, -1], [1, 0, -1], [1, -1, 0]])
    Lc = locked.K * ((2 + a)[:, None, None] * np.eye(3) + t[:, None, None] * C)
    Xc = np.repeat(d["Xs"][:1], len(s), axis=0)
    cname = "bars intact, steps <= 0.01, balls overlap (fold-continuous-fcc curve run past fcc to s = 1.1)"
    bad = locked.verify_arrays(1, Lc, Xc, want_rigidity=False)
    kinds = {r.split(":")[0] for r in bad["reasons"]}
    certlib.require("overlap" in kinds and kinds <= {"overlap", "midpoint-overlap"},
                    f"locked checker: the overlap control is not rejected by the overlap tests alone: {sorted(kinds)}")
    certlib.must_fail("locked checker: " + cname, bad["valid"], log)
    rep = _independent(Lc, Xc)
    certlib.require(rep["overlap_violations"] > 0 and chk.accepts(dict(rep, overlap_violations=0), interp=False),
                    "independent checker: the overlap control is not rejected by the overlap test alone")
    certlib.must_fail("independent checker: " + cname, chk.accepts(rep, interp=False), log)
    # Ls longer than Xs: 34 extra lattice-only states, each 0.25% smaller than the last. chk_independent.py does
    # not yet compare the lengths of Ls and Xs (it accepts this file), so only the locked checker is tested here.
    L2 = np.concatenate([d["Ls"], d["Ls"][-1] * 0.9975 ** np.arange(1, 35)[:, None, None]])
    bad = locked.verify_arrays(1, L2, d["Xs"], want_rigidity=False)
    certlib.require(all(r.startswith("format: bad shapes") for r in bad["reasons"]),
                    f"locked checker: the length-mismatch control is not rejected by the format test alone: {bad['reasons']}")
    certlib.must_fail("locked checker: Ls longer than Xs by 34 lattice-only states (format)", bad["valid"], log)
    # The jamming test on a state that is not jammed: the fcc path's next-to-last state.
    I, J, S = locked.bar_graph(1)
    jam = locked.end_state_report(1, d["Ls"][-2], d["Xs"][-2], I, J, S)
    certlib.must_fail("locked checker: first-order jamming test on the next-to-last fcc state",
                      jam["first_order_jammed"] is not False, log)
    cert["must_fail_controls_rejected"] = log

    cert["inputs"] = certlib.inputs(PATHS["hcp"], PATHS["fcc"], FOLD / "verify.py",
                                    ROOT / "results" / "fold-close-packing" / "chk_independent.py",
                                    ROOT / "results" / "fold-end-exact" / "exact_end.py",
                                    HERE / "search_hcp_m1.py", HERE / "search_fcc_m1.py",
                                    HERE / "verify.py", ROOT / "tools" / "certlib.py")
    certlib.write_json(HERE / "certificate.json", cert)
    print("fold-m1: OK (both checkers and the exact check accept the m = 1 hcp and fcc paths; controls rejected: %d)" % len(log))


def _independent(Ls, Xs):
    """The independent checker, without interpolation, on arrays (it reads a path file)."""
    with tempfile.TemporaryDirectory() as td:
        p = pathlib.Path(td) / "bad.npz"
        np.savez(p, m=1, Ls=Ls, Xs=Xs)
        return chk.check(p, interp=False)


def _set(a, idx, val):
    a = a.copy()
    a[idx] = val
    return a


if __name__ == "__main__":
    main()
