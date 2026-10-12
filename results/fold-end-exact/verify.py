#!/usr/bin/env python3
"""Replay for results/fold-end-exact.

The end states of the two stored fold paths (hcp: autoresearch-fold/best/path.npz; fcc:
results/fold-close-packing/path_fcc.npz) are exact. For each one this script
  1. snaps the final state to the integer grid in units K/den (K = 2 - sqrt2) in the cube frame
     of octahedron 0 (exact_end.exactify; fails if any coordinate is more than 1e-7 off the grid);
  2. builds the closed form written in PROOF.md and checks that it is the same integer state;
  3. checks the closed form in integer arithmetic (exact_end.check_exact): every bar exactly 2R, no
     pair of balls closer than 2R over all periodic images, 12 contacts per ball, density exactly
     pi/sqrt18, and the common-neighbour signature (hcp or fcc);
  4. runs must-fail controls and writes certificate.json.

Dependencies: numpy (I/O and the snap only). About 30 seconds.   Usage: python verify.py
"""
import itertools
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(HERE))
import certlib  # noqa: E402
import exact_end  # noqa: E402

import numpy as np  # noqa: E402

HCP = ROOT / "autoresearch-fold" / "best" / "path.npz"
FCC = ROOT / "results" / "fold-close-packing" / "path_fcc.npz"

# Closed forms (PROOF.md), m = 2, in units of K. Cell (ix, iy, iz), 0 <= ix, iy, iz <= 1:
#   centre = ix t_x + iy t_y + iz t_z,  dot (axis a, sign s) = centre + s e_{pi(a)},
#   pi = identity on even cells (ix+iy+iz even), (x, y, z) -> (z, x, y) on odd cells.
T_X, T_Y = (1, -1, 2), (1, 2, -1)
CLOSED = {
    "hcp": {"den": 3, "t_z": ("-4/3", "4/3", "4/3"), "L": [(2, -2, 4), (2, 4, -2), ("-8/3", "8/3", "8/3")]},
    "fcc": {"den": 1, "t_z": (-1, 1, 2), "L": [(2, -2, 4), (2, 4, -2), (-2, 3, 3)]},
}


def closed_form(kind, t_z=None, L=None, permute=True, den=None):
    from fractions import Fraction as F
    c = CLOSED[kind]
    den = den or c["den"]
    t_z = t_z or c["t_z"]
    L = L or c["L"]
    unit = lambda v: int(F(v) * den)  # noqa: E731
    Y = []
    for ix, iy, iz in itertools.product(range(2), repeat=3):
        ctr = [unit(ix * F(T_X[k]) + iy * F(T_Y[k]) + iz * F(t_z[k])) for k in range(3)]
        odd = (ix + iy + iz) % 2 == 1 and permute
        for a in range(3):
            b = (2, 0, 1)[a] if odd else a
            for s in (1, -1):
                p = list(ctr)
                p[b] += s * den
                Y.append(p)
    return Y, [[unit(v) for v in r] for r in L], den


def summary(rep):
    keep = ["m", "balls", "bars", "unit", "lattice_rows", "det_L", "bars_exact", "no_overlap_exact",
            "min_pair_distance_over_2R_squared", "contacts_per_ball", "density_is_pi_over_sqrt18_exactly",
            "cna_pairs", "structure", "exact_ok"]
    return {k: rep[k] for k in keep}


def main():
    cert = {"result": "fold-end-exact", "paths": {}}
    for kind, path in (("hcp", HCP), ("fcc", FCC)):
        z = np.load(path)
        Y, Lu, den = exact_end.exactify(z["Xs"][-1], z["Ls"][-1])
        Yc, Lc, dc = closed_form(kind)
        certlib.require((Y, Lu, den) == (Yc, Lc, dc), f"{kind}: stored end state is not the closed form of PROOF.md")
        rep = exact_end.check_exact(2, Yc, Lc, dc)
        certlib.require(rep["exact_ok"], f"{kind}: closed form fails an exact check")
        certlib.require(rep["structure"] == kind, f"{kind}: common-neighbour analysis gives {rep['structure']}")
        cert["paths"][kind] = {"path": path.relative_to(ROOT).as_posix(), "end_state_on_integer_grid_to_1e-7": True,
                               "end_state_equals_closed_form": True, "exact": summary(rep)}

    log = []
    Y, L, den = closed_form("hcp")
    Y[7] = [Y[7][0] + 1, Y[7][1], Y[7][2]]
    certlib.must_fail("hcp closed form with one dot moved by K/3", exact_end.check_exact(2, Y, L, den)["exact_ok"], log)
    certlib.must_fail("hcp closed form without the cyclic axis permutation on odd cells",
                      exact_end.check_exact(2, *closed_form("hcp", permute=False))["exact_ok"], log)
    certlib.must_fail("fcc closed form with third lattice row (-2, 3, 4) instead of (-2, 3, 3)",
                      exact_end.check_exact(2, *closed_form("fcc", L=[(2, -2, 4), (2, 4, -2), (-2, 3, 4)]))["exact_ok"], log)
    # (-2, 2, 4) = 2 t_z is NOT a control: it is a second exact fcc end state (PROOF.md), checked here
    alt = exact_end.check_exact(2, *closed_form("fcc", L=[(2, -2, 4), (2, 4, -2), (-2, 2, 4)]))
    certlib.require(alt["exact_ok"] and alt["structure"] == "fcc", "fcc with third lattice row 2 t_z is not exact fcc")
    cert["fcc_with_third_row_2_t_z"] = summary(alt)
    rep = exact_end.check_exact(2, *closed_form("fcc", t_z=("-4/3", "4/3", "4/3"), L=CLOSED["hcp"]["L"], den=3))
    certlib.must_fail("hcp closed form labelled fcc (common-neighbour test)", rep["structure"] == "fcc", log)
    cert["must_fail_controls_rejected"] = log

    cert["inputs"] = certlib.inputs(HCP, FCC, HERE / "exact_end.py")
    certlib.write_json(HERE / "certificate.json", cert)
    print("fold-end-exact: OK (hcp and fcc end states are exact closed forms; controls rejected: %d)" % len(log))


if __name__ == "__main__":
    main()
