#!/usr/bin/env python3
"""Replay for results/fold-close-packing.

Runs BOTH checkers on the hcp fold path and on the fcc extension of it:
  1. autoresearch-fold/verify.py          (locked checker: bars, overlaps, continuity, jamming)
  2. results/fold-close-packing/chk_independent.py   (independent re-implementation + CNA)
Regenerates the fcc extension from the hcp path and compares it with the stored file,
runs must-fail controls, and writes certificate.json. The checkers are imported and run in
this process, so their sha256 pins are checked first, with the standard library only.

Dependencies: numpy, scipy (scipy only through autoresearch-fold/verify.py).
Slow target (about 3 minutes): the independent checker checks 9 interior sub-states
(u = 0.1, ..., 0.9) of every step.   Usage: python verify.py
"""
import hashlib
import itertools
import pathlib
import sys
import tempfile

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FOLD = ROOT / "autoresearch-fold"


def _pins_ok(sumfile, folder, names):
    """True if sumfile ('sha256  name' lines) pins every file in `names` and every file it pins
    in folder matches. Standard library only: it runs before any checker is imported."""
    pins = {}
    for ln in sumfile.read_text().splitlines():
        if ln.strip():
            digest, _, name = ln.partition("  ")
            pins[name] = digest
    return set(names) <= set(pins) and all(
        (folder / name).is_file() and hashlib.sha256((folder / name).read_bytes()).hexdigest() == digest
        for name, digest in pins.items())


for _sums, _folder, _names in ((FOLD / "verify.sha256", FOLD, ["verify.py"]),
                               (HERE / "checkers.sha256", HERE, ["chk_independent.py", "make_fcc_path.py"])):
    if not _pins_ok(_sums, _folder, _names):
        raise SystemExit(f"REFUSING TO RUN: {_sums.relative_to(ROOT)} does not pin {', '.join(_names)}, "
                         "or a file it pins does not match its sha256")

sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(HERE))
import certlib  # noqa: E402
import chk_independent as chk  # noqa: E402
import make_fcc_path  # noqa: E402
locked = make_fcc_path.verify  # autoresearch-fold/verify.py (hash-locked), loaded by path

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
    # Each control above breaks two or three properties at once. Each one below breaks one, and the
    # replay asserts which clauses reject it, not only that it is rejected.
    m = int(d["m"])
    Ls_c, Xs_c = cubic_path(m, 1.25, 625)  # steps of 0.002 rad; past 60 degrees the balls overlap, bars stay 2R
    single_property_control("cubic [111] path continued to 1.25 rad (overlaps only)", m, Ls_c, Xs_c,
                            {"overlap"}, {"overlap", "midpoint-overlap"}, log)
    # the locked checker's overlap search includes bonded pairs, so it reports the shortened bar as an overlap too
    single_property_control("one dot displaced by 2e-6 (bar lengths broken, far below the step limit)",
                            m, d["Ls"], _set(d["Xs"], (400, 0, 0), d["Xs"][400, 0, 0] + 2e-6),
                            {"bars"}, {"bar", "overlap"}, log)
    # Ls longer than Xs: 34 lattice-only states appended, each 0.9975 times the one before.
    # Locked checker only: chk_independent.py does not compare the lengths of Ls and Xs.
    Ls_x = np.concatenate([d["Ls"], d["Ls"][-1] * 0.9975 ** np.arange(1, 35)[:, None, None]])
    cname = "Ls longer than Xs by 34 lattice-only states (format)"
    bad = locked.verify_arrays(m, Ls_x, d["Xs"], want_rigidity=False)
    certlib.require(_reason_kinds(bad) == {"format"}, f"control '{cname}': locked checker reasons {bad['reasons']}")
    certlib.must_fail("locked checker: " + cname, bad["valid"], log)
    cert["must_fail_controls_rejected"] = log

    cert["inputs"] = certlib.inputs(HCP, FCC, FOLD / "verify.py", FOLD / "test_verify.py", HERE / "chk_independent.py",
                                    HERE / "make_fcc_path.py", HERE / "verify.py", ROOT / "tools" / "certlib.py")
    certlib.write_json(HERE / "certificate.json", cert)
    print("fold-close-packing: OK (both checkers accept the hcp and fcc paths; controls rejected: %d)" % len(log))


def _set(a, idx, val):
    a = a.copy()
    a[idx] = val
    return a


def _reason_kinds(rep):
    """The kinds ('bar', 'overlap', 'midpoint-overlap', 'continuity', 'format', ...) of the locked checker's reasons."""
    return {r.split(":")[0] for r in rep["reasons"]}


def _independent_fails(rep):
    """The clauses of chk.accepts(rep, interp=False) that rep fails, one by one."""
    held = {"start state": rep["canonical_bar_err"] < 1e-12 and rep["state0_dev_from_canonical"] < 1e-12
            and rep["lattice0_dev_from_canonical"] < 1e-12,
            "bars": rep["max_bar_rel_err"] < 1e-9, "overlap": rep["overlap_violations"] == 0,
            "dot step": rep["max_dot_step"] <= 0.01, "lattice step": rep["max_lattice_step"] <= 0.01}
    return {k for k, ok in held.items() if not ok}


def single_property_control(cname, m, Ls, Xs, independent_fails, locked_kinds, log):
    """Must-fail control for an input that breaks one property: both checkers must reject it, the
    independent checker failing exactly the clauses `independent_fails`, and the locked checker giving
    reasons of exactly the kinds `locked_kinds`."""
    bad = locked.verify_arrays(m, Ls, Xs, want_rigidity=False)
    certlib.require(_reason_kinds(bad) == locked_kinds, f"control '{cname}': locked checker reasons {bad['reasons']}")
    certlib.must_fail("locked checker: " + cname, bad["valid"], log)
    with tempfile.TemporaryDirectory() as td:
        p = pathlib.Path(td) / "bad.npz"
        np.savez(p, m=m, Ls=Ls, Xs=Xs)
        rep = chk.check(p, interp=False)
    fails, accepted = _independent_fails(rep), chk.accepts(rep, interp=False)
    certlib.require(fails == independent_fails and accepted == (not fails),
                    f"control '{cname}': independent checker fails {sorted(fails)}, accepts() = {accepted}")
    certlib.must_fail("independent checker: " + cname, accepted, log)


def cubic_path(m, th_end, steps):
    """The cubic [111] family of autoresearch-fold/test_verify.py, rebuilt here: the octahedra of even
    and odd cells turn by +th and -th about [111], and the cubic cell side 2a solves
    (2a - 2Kc - q)^2 + 2q^2 = (2R)^2 (c = cos th, q = 2K(1 - c)/3, a = 1 at th = 0), so every bar is
    2R. th runs from 0 to th_end in `steps` equal steps."""
    K, R = 2 - np.sqrt(2), np.sqrt(2) - 1
    n = np.ones(3) / np.sqrt(3)
    Kx = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    Ls, Xs = [], []
    for th in np.linspace(0, th_end, steps + 1):
        c = np.cos(th)
        q = 2 * K * (1 - c) / 3
        a = (2 * K * c + q + np.sqrt(4 * R ** 2 - 2 * q ** 2)) / 2
        X = np.zeros((6 * m ** 3, 3))
        for ix, iy, iz in itertools.product(range(m), repeat=3):
            t = th if (ix + iy + iz) % 2 == 0 else -th
            rot = np.eye(3) + np.sin(t) * Kx + (1 - np.cos(t)) * Kx @ Kx
            cell = (ix * m + iy) * m + iz
            for ax in range(3):
                for sg in range(2):  # dot 6 cell + 2 axis + (0 for +e_axis, 1 for -e_axis), as in both checkers
                    X[6 * cell + 2 * ax + sg] = 2 * a * np.array([ix, iy, iz]) + a + K * (1 - 2 * sg) * rot[:, ax]
        Ls.append(2 * a * m * np.eye(3))
        Xs.append(X)
    return np.array(Ls), np.array(Xs)


if __name__ == "__main__":
    main()
