#!/usr/bin/env python3
"""LOCKED trusted checker for the folding / packing search rig. Do NOT edit.

Usage: python verify.py path.npz        (prints a JSON report)

Path file (.npz) keys:
  m      : int, supercell is m x m x m cells (1 <= m <= 4)
  Ls     : (T,3,3) lattice matrices, ROWS are the three supercell vectors
  Xs     : (T,6*m^3,3) Cartesian dot positions (unwrapped, continuous in time)
Any other keys (e.g. a "bars" entry) are ignored: the bar graph is built here.

Model. Cell side 2, dots at c +- k e_i, k = 2 - sqrt2, ball radius R = sqrt2 - 1.
Dot index = 6*cell + 2*axis + (0 for +e_axis, 1 for -e_axis); cell = (ix*m+iy)*m+iz;
cell centre = 2*(ix,iy,iz) + (1,1,1). Canonical lattice = diag(2m,2m,2m).
Bars (length 2R): 12 same-cell pairs of dots on different axes; for each axis the +e dot of
a cell to the -e dot of the neighbour in +axis direction. Each bar carries a FIXED integer
lattice shift (set by the canonical state), bar vector = x_j + s@L - x_i.
"""
import json
import sys
import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import linprog

R = np.sqrt(2.0) - 1.0
K = 2.0 - np.sqrt(2.0)
D = 2.0 * R                      # bar length / ball diameter
TOL_BAR = 1e-9                   # relative bar-length tolerance
TOL_OVL = 1e-9                   # relative overlap tolerance
TOL_MID = 1e-4                   # overlap tolerance at interpolated midpoints
MAX_STEP_X = 0.01                # max dot displacement between consecutive states
MAX_STEP_L = 0.01                # max Frobenius change of L between consecutive states
TOL_INIT = 1e-9
CONTACT_TOL = 1e-6
MAX_STATES = 200000
MAX_REASONS = 40


def canonical_state(m):
    L = 2.0 * m * np.eye(3)
    X = np.zeros((6 * m ** 3, 3))
    for ix in range(m):
        for iy in range(m):
            for iz in range(m):
                c = (ix * m + iy) * m + iz
                ctr = 2.0 * np.array([ix, iy, iz]) + 1.0
                for a in range(3):
                    for sg in range(2):
                        e = np.zeros(3)
                        e[a] = 1.0 if sg == 0 else -1.0
                        X[6 * c + 2 * a + sg] = ctr + K * e
    return L, X


def bar_graph(m):
    """Returns (I, J, S): bar b joins dot I[b] to dot J[b] + S[b]@L (S integer (B,3))."""
    I, J, S = [], [], []
    for ix in range(m):
        for iy in range(m):
            for iz in range(m):
                c = (ix * m + iy) * m + iz
                for a in range(3):
                    for b in range(a + 1, 3):
                        for sa in range(2):
                            for sb in range(2):
                                I.append(6 * c + 2 * a + sa)
                                J.append(6 * c + 2 * b + sb)
                                S.append((0, 0, 0))
                idx = [ix, iy, iz]
                for a in range(3):
                    nb = list(idx)
                    nb[a] = (nb[a] + 1) % m
                    shift = [0, 0, 0]
                    if idx[a] + 1 == m:
                        shift[a] = 1
                    c2 = (nb[0] * m + nb[1]) * m + nb[2]
                    I.append(6 * c + 2 * a + 0)
                    J.append(6 * c2 + 2 * a + 1)
                    S.append(tuple(shift))
    return np.array(I), np.array(J), np.array(S, dtype=float)


def _canon_key(i, j, s):
    s = tuple(int(round(v)) for v in s)
    if i > j or (i == j and s < (0, 0, 0)):
        i, j = j, i
        s = tuple(-v for v in s)
    return (int(i), int(j), s)


def bar_vectors(X, L, I, J, S):
    return X[J] + S @ L - X[I]


def _image_cloud(X, L, rcut):
    """Wrap into the cell and return (tree, owner index, shift) for the image cloud that
    contains EVERY image within rcut of every wrapped base dot (rigorous box bound)."""
    Linv = np.linalg.inv(L)
    F = X @ Linv
    F = F - np.floor(F)
    Xw = F @ L
    colnorm = np.linalg.norm(Linv, axis=0)       # |v_k| <= |r| * ||Linv[:,k]||
    n = np.ceil(rcut * colnorm).astype(int) + 1
    if np.prod(2 * n + 1) > 4000:
        raise ValueError("lattice too degenerate (image box too large)")
    rng = [np.arange(-k, k + 1) for k in n]
    shifts = np.array(np.meshgrid(*rng, indexing="ij")).reshape(3, -1).T.astype(float)
    pts = (Xw[None, :, :] + (shifts @ L)[:, None, :]).reshape(-1, 3)
    owner = np.tile(np.arange(len(X)), len(shifts))
    shp = np.repeat(shifts, len(X), axis=0)
    return Xw, pts, owner, shp


def close_pairs(X, L, rcut):
    """All unordered pairs (i,j,s) (s integer shift w.r.t. WRAPPED coordinates is NOT used;
    returned shift is relative to the given unwrapped X) with distance < rcut, i==j allowed
    for nonzero shift. Returns arrays i, j, s, dist."""
    Linv = np.linalg.inv(L)
    F = X @ Linv
    W = np.floor(F)                      # wrap offset of each dot (integer)
    Xw, pts, owner, shp = _image_cloud(X, L, rcut)
    tree = cKDTree(pts)
    base = cKDTree(Xw)
    pairs = base.query_ball_tree(tree, rcut)
    out_i, out_j, out_s, out_d = [], [], [], []
    for i, lst in enumerate(pairs):
        for p in lst:
            j = owner[p]
            s = shp[p]
            if j == i and not np.any(s):
                continue
            # key to deduplicate: keep each unordered pair once
            if j < i or (j == i and tuple(s) < (0, 0, 0)):
                continue
            # convert shift from wrapped to unwrapped coordinates: x_j^w = x_j - W_j L
            su = s + W[i] - W[j]
            out_i.append(i)
            out_j.append(j)
            out_s.append(su)
            out_d.append(np.linalg.norm(pts[p] - Xw[i]))
    if not out_i:
        return (np.zeros(0, int), np.zeros(0, int), np.zeros((0, 3)), np.zeros(0))
    return (np.array(out_i), np.array(out_j), np.array(out_s), np.array(out_d))


def min_pair_distance(X, L, rcut=D * 1.05):
    i, j, s, d = close_pairs(X, L, rcut)
    return (d.min() if len(d) else np.inf)


def verify_arrays(m, Ls, Xs, want_rigidity=True):
    reasons = []
    res = {"valid": False, "reasons": reasons}

    def fail(msg):
        if len(reasons) < MAX_REASONS:
            reasons.append(msg)

    try:
        m = int(m)
    except Exception:
        fail("format: bad m")
        return res
    if not (1 <= m <= 4):
        fail("format: m must be in 1..4")
        return res
    N = 6 * m ** 3
    try:
        Ls = np.asarray(Ls, dtype=float)
        Xs = np.asarray(Xs, dtype=float)
    except Exception:
        fail("format: arrays not numeric")
        return res
    if Ls.ndim != 3 or Ls.shape[1:] != (3, 3) or Xs.ndim != 3 or Xs.shape[1:] != (N, 3) \
            or Ls.shape[0] != Xs.shape[0] or Ls.shape[0] < 1:
        fail("format: bad shapes (need Ls (T,3,3), Xs (T,%d,3))" % N)
        return res
    T = Ls.shape[0]
    if T > MAX_STATES:
        fail("format: too many states")
        return res
    if not (np.all(np.isfinite(Ls)) and np.all(np.isfinite(Xs))):
        fail("nonfinite: NaN or inf in path")
        return res

    L0, X0 = canonical_state(m)
    if np.abs(Ls[0] - L0).max() > TOL_INIT or np.abs(Xs[0] - X0).max() > TOL_INIT:
        fail("initial: first state is not the canonical initial state")

    I, J, S = bar_graph(m)
    dets = np.linalg.det(Ls)
    for t in np.nonzero(dets <= 0)[0][:5]:
        fail("det: det L <= 0 at state %d" % t)
    nbad_bar = 0
    nbad_ovl = 0
    dens = np.full(T, np.nan)
    vol_balls = N * (4.0 / 3.0) * np.pi * R ** 3
    for t in range(T):
        if dets[t] <= 0:
            continue
        dens[t] = vol_balls / dets[t]
        bv = bar_vectors(Xs[t], Ls[t], I, J, S)
        rel = np.abs(np.linalg.norm(bv, axis=1) / D - 1.0).max()
        if rel > TOL_BAR:
            nbad_bar += 1
            if nbad_bar <= 3:
                fail("bar: state %d bar length rel error %.3e > %.0e" % (t, rel, TOL_BAR))
        try:
            md = min_pair_distance(Xs[t], Ls[t])
        except Exception as e:
            fail("overlap: state %d: %s" % (t, e))
            nbad_ovl += 1
            continue
        if md < D * (1.0 - TOL_OVL):
            nbad_ovl += 1
            if nbad_ovl <= 3:
                fail("overlap: state %d min pair distance %.12g < %.12g" % (t, md, D * (1 - TOL_OVL)))
    if nbad_bar > 3:
        fail("bar: %d states total violate bar lengths" % nbad_bar)
    if nbad_ovl > 3:
        fail("overlap: %d states total overlap" % nbad_ovl)

    nbad_c = 0
    nbad_mid = 0
    for t in range(T - 1):
        dx = np.linalg.norm(Xs[t + 1] - Xs[t], axis=1).max()
        dl = np.linalg.norm(Ls[t + 1] - Ls[t])
        if dx > MAX_STEP_X + 1e-12 or dl > MAX_STEP_L + 1e-12:
            nbad_c += 1
            if nbad_c <= 3:
                fail("continuity: step %d->%d dot move %.5f / lattice change %.5f exceed %.3f"
                     % (t, t + 1, dx, dl, MAX_STEP_X))
            continue
        Lm = 0.5 * (Ls[t] + Ls[t + 1])
        Xm = 0.5 * (Xs[t] + Xs[t + 1])
        if np.linalg.det(Lm) <= 0:
            nbad_mid += 1
            fail("det: midpoint of step %d has det <= 0" % t)
            continue
        try:
            md = min_pair_distance(Xm, Lm)
        except Exception as e:
            nbad_mid += 1
            fail("midpoint: step %d: %s" % (t, e))
            continue
        if md < D * (1.0 - TOL_MID):
            nbad_mid += 1
            if nbad_mid <= 3:
                fail("midpoint-overlap: step %d min distance %.8g" % (t, md))
    if nbad_c > 3:
        fail("continuity: %d steps total too large" % nbad_c)

    res["states"] = int(T)
    res["m"] = m
    ok_dens = np.isfinite(dens)
    if ok_dens.any():
        res["final_density"] = float(dens[-1]) if np.isfinite(dens[-1]) else None
        res["max_density"] = float(np.nanmax(dens))
    else:
        res["final_density"] = None
        res["max_density"] = None
    if np.isfinite(dens[-1]) and np.all(np.isfinite(Xs[-1])):
        try:
            res.update(end_state_report(m, Ls[-1], Xs[-1], I, J, S, want_rigidity))
        except Exception as e:
            # Fail closed: an end-state analysis that throws must not report a valid path.
            fail("end-state analysis failed: %s" % e)
    res["valid"] = len(reasons) == 0     # after the end-state analysis, so its failure counts
    return res


def end_state_report(m, L, X, I, J, S, want_rigidity=True):
    out = {}
    N = len(X)
    ci, cj, cs, cd = close_pairs(X, L, D * (1.0 + CONTACT_TOL))
    deg = np.zeros(N, int)
    np.add.at(deg, ci, 1)
    np.add.at(deg, cj, 1)
    out["contacts_per_ball_min"] = int(deg.min())
    out["contacts_per_ball_mean"] = float(deg.mean())
    out["contacts_per_ball_max"] = int(deg.max())
    out["contact_pairs"] = int(len(ci))
    out["contact_degree_histogram"] = {int(k): int((deg == k).sum()) for k in np.unique(deg)}
    out["end_lattice_lengths"] = [float(v) for v in np.linalg.norm(L, axis=1)]
    out["end_lattice_cos_angles"] = [float(np.dot(L[a], L[b]) / np.linalg.norm(L[a]) / np.linalg.norm(L[b]))
                                     for a, b in ((1, 2), (0, 2), (0, 1))]
    if not want_rigidity:
        return out
    bar_keys = {_canon_key(i, j, s) for i, j, s in zip(I, J, S)}
    extra = [(i, j, s) for i, j, s in zip(ci, cj, cs) if _canon_key(i, j, s) not in bar_keys]
    out["extra_contacts_beyond_bars"] = len(extra)
    n = 3 * N + 9

    def rows(triples):
        A = np.zeros((len(triples), n))
        for k, (i, j, s) in enumerate(triples):
            r = X[j] + np.asarray(s, float) @ L - X[i]
            A[k, 3 * i:3 * i + 3] -= r
            A[k, 3 * j:3 * j + 3] += r
            A[k, 3 * N:] += np.outer(np.asarray(s, float), r).ravel()
        return A

    Ab = rows(list(zip(I, J, S)))
    Ac = rows(extra) if extra else np.zeros((0, n))

    def nullity(A):
        sv = np.linalg.svd(A, compute_uv=False)
        tol = 1e-8 * max(sv.max(), 1.0)
        rank = int((sv > tol).sum())
        return n - rank

    out["flex_bars_only"] = nullity(Ab) - 6
    out["flex_bars_plus_contacts"] = nullity(np.vstack([Ab, Ac])) - 6
    # first-order density-ascent LP: maximise -d ln(det L)/dt  s.t. bars fixed, contacts non-decreasing
    Linv = np.linalg.inv(L)
    c = np.zeros(n)
    c[3 * N:] = Linv.T.ravel()          # d ln det = sum_ab Linv[b,a] Ldot[a,b]; minimise it -> density up
    # minimise d ln det (volume) means density up is  -d ln det  ; we MINIMISE d ln det
    kw = dict(bounds=[(-1, 1)] * n, method="highs")
    if len(Ac):
        lp = linprog(c, A_ub=-Ac, b_ub=np.zeros(len(Ac)), A_eq=Ab, b_eq=np.zeros(len(Ab)), **kw)
    else:
        lp = linprog(c, A_eq=Ab, b_eq=np.zeros(len(Ab)), **kw)
    if lp.status == 0:
        out["first_order_density_ascent_rate"] = float(-lp.fun)   # d ln(density)/dt per unit max speed
        out["first_order_jammed"] = bool(-lp.fun < 1e-7)
    else:
        out["first_order_density_ascent_rate"] = None
        out["first_order_jammed"] = None
    return out


def verify_file(path):
    try:
        z = np.load(path, allow_pickle=False)
        m = int(z["m"])
        Ls = z["Ls"]
        Xs = z["Xs"]
    except Exception as e:
        return {"valid": False, "reasons": ["format: cannot load %s: %s" % (path, e)]}
    return verify_arrays(m, Ls, Xs)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: verify.py path.npz", file=sys.stderr)
        sys.exit(2)
    rep = verify_file(sys.argv[1])
    print(json.dumps(rep, indent=1))
    sys.exit(0 if rep["valid"] else 1)
