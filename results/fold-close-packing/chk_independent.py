#!/usr/bin/env python3
"""Independent checker for fold paths (numpy only; does NOT import autoresearch-fold/verify.py).

This is autoresearch-fold/independent/chk.py reorganised into a function that takes the path
file as an argument and returns its measurements, plus a contact-based common-neighbour
analysis (CNA) of the end state. The checks are the same as in chk.py:

  * the bar graph is rebuilt from the model (cell side 2, dots at centre +- (2-sqrt2) e_i,
    120 bars of length 2R = 2(sqrt2-1) for m = 2) with FIXED integer lattice shifts;
  * every bar length is checked in every state;
  * the minimum distance between non-bonded dots, over all periodic images, is checked in
    every state (overlap = distance < 2R), over all image shifts derived from the lattice;
  * step sizes between consecutive states; optional 10-substate linear interpolation
    between every pair of states (the midpoint-style overlap check);
  * the end state: density, contacts per ball, CNA signatures.

Usage: python chk_independent.py PATH.npz [--no-interp]
"""
import collections
import itertools
import sys

import numpy as np

S2 = np.sqrt(2.0)
R = S2 - 1.0
KAP = 2.0 - S2
B = 2.0 * R


def build_bars(m):
    cells = list(itertools.product(range(m), repeat=3))
    dots, idx = [], {}
    for c in cells:
        for ax in range(3):
            for s in (+1, -1):
                idx[(c, ax, s)] = len(dots)
                p = 2 * np.array(c) + 1.0
                p[ax] += s * KAP
                dots.append(p)
    P0 = np.array(dots)
    L0 = 2 * m * np.eye(3)
    bars = []
    for c in cells:
        for a in range(3):
            for b in range(a + 1, 3):
                for sa in (1, -1):
                    for sb in (1, -1):
                        bars.append((idx[(c, a, sa)], idx[(c, b, sb)], (0, 0, 0)))
        for a in range(3):
            c2 = list(c)
            c2[a] += 1
            n = [0, 0, 0]
            if c2[a] == m:
                c2[a] = 0
                n[a] = 1
            bars.append((idx[(c, a, 1)], idx[(tuple(c2), a, -1)], tuple(n)))
    return P0, L0, bars


def shift_range(X, L, rcut):
    """Integer bound rngn such that EVERY lattice shift n with |X_j - X_i + n L| <= rcut for some
    dots i, j lies in [-rngn, rngn]^3 (rigorous, derived from the lattice and the dots).
    With fractional coordinates F = X L^-1, a vector v of length <= rcut has |v_k| <= rcut |Linv[:,k]|
    along lattice axis k, so |n_k| <= max|F_j,k - F_i,k| + rcut |Linv[:,k]|."""
    Linv = np.linalg.inv(L)
    F = X @ Linv
    spread = F.max(0) - F.min(0)
    return int(np.ceil((spread + rcut * np.linalg.norm(Linv, axis=0)).max()))


def bar_lengths(X, L, bars):
    return np.array([np.linalg.norm(X[j] + np.array(n) @ L - X[i]) for i, j, n in bars])


class Mind:
    """Minimum non-bonded distance and the full distance table. The image-shift block is derived
    from the lattice and dots of each state (shift_range with rcut = 2 * 2R, so every image that
    could be a non-bonded overlap or the next-nearest shell is included), never a fixed size.
    Blocks and bond masks are cached per rngn."""

    def __init__(self, bars, N, rcut=2 * B):
        self.bars, self.N, self.rcut = bars, N, rcut
        self._cache = {}

    def _block(self, rngn):
        if rngn not in self._cache:
            ns = np.array(list(itertools.product(range(-rngn, rngn + 1), repeat=3)))
            mask = np.ones((self.N, self.N, len(ns)), bool)
            zero = np.where((ns == 0).all(1))[0][0]
            for a in range(self.N):
                mask[a, a, zero] = False
            for i, j, n in self.bars:
                s = np.where((ns == np.array(n)).all(1))[0][0]
                mask[i, j, s] = False
                s2 = np.where((ns == -np.array(n)).all(1))[0][0]
                mask[j, i, s2] = False
            self._cache[rngn] = (ns, mask)
        return self._cache[rngn]

    def table(self, X, L):
        return self._table_mask(X, L)[0]

    def _table_mask(self, X, L):
        ns, mask = self._block(max(shift_range(X, L, self.rcut), 1))
        sh = ns @ L
        D = X[None, :, None, :] - X[:, None, None, :] + sh[None, None, :, :]
        return np.linalg.norm(D, axis=-1), mask

    def __call__(self, X, L):
        r, mask = self._table_mask(X, L)
        return r[mask].min(), r


def cna_signatures(X, L, tol=1e-6):
    """Common-neighbour analysis from contacts (distance < 2R(1+tol)); returns a Counter of
    (n_common, n_bonds, 2 if the bonds share an atom else 1) over all directed contacts."""
    N = len(X)
    rngn = max(shift_range(X, L, B * (1 + tol)), 1)
    ns = list(itertools.product(range(-rngn, rngn + 1), repeat=3))
    ct = collections.defaultdict(set)
    for i in range(N):
        for j in range(N):
            for n in ns:
                if i == j and n == (0, 0, 0):
                    continue
                if np.linalg.norm(X[j] - X[i] + np.array(n) @ L) < B * (1 + tol):
                    ct[i].add((j, n))

    def nb(i, off):
        return {(j, tuple(np.add(s, off))) for j, s in ct[i]}

    sig = collections.Counter()
    for i in range(N):
        n0 = nb(i, (0, 0, 0))
        for (j, s) in n0:
            com = n0 & nb(j, s)
            bl = [(p, q) for p in com for q in com if p < q and q in nb(*p)]
            shared = len(bl) == 2 and len(set(bl[0]) & set(bl[1])) > 0
            sig[(len(com), len(bl), 2 if shared else 1)] += 1
    return sig


def check(path, interp=True):
    d = np.load(path)
    Ls, Xs, m = d["Ls"], d["Xs"], int(d["m"])
    T, N = len(Xs), Xs.shape[1]
    P0, L0, bars = build_bars(m)
    out = {"states": T, "dots": N, "bars": len(bars)}
    out["canonical_bar_err"] = float(abs(bar_lengths(P0, L0, bars) - B).max())
    out["state0_dev_from_canonical"] = float(abs(Xs[0] - P0).max())
    out["lattice0_dev_from_canonical"] = float(abs(Ls[0] - L0).max())
    bar_err = np.array([abs(bar_lengths(Xs[t], Ls[t], bars) - B).max() for t in range(T)])
    out["max_bar_rel_err"] = float(bar_err.max() / B)
    md = Mind(bars, N)
    ovl = np.array([md(Xs[t], Ls[t])[0] for t in range(T)])
    out["min_nonbonded_over_2R"] = float(ovl.min() / B)
    out["overlap_violations"] = int((ovl < B * (1 - 1e-9)).sum())
    dX = np.linalg.norm(Xs[1:] - Xs[:-1], axis=2).max(1)
    dL = np.abs(Ls[1:] - Ls[:-1]).max((1, 2))
    out["max_dot_step"] = float(dX.max())
    out["max_lattice_step"] = float(dL.max())
    if interp:
        w = 10
        worst, bw = 1e9, 0.0
        for t in range(T - 1):
            for a in range(1, w):
                u = a / w
                X = (1 - u) * Xs[t] + u * Xs[t + 1]
                L = (1 - u) * Ls[t] + u * Ls[t + 1]
                mn, _ = md(X, L)
                worst = min(worst, mn)
                bw = max(bw, abs(bar_lengths(X, L, bars) - B).max())
        out["interp_min_nonbonded_over_2R"] = float(worst / B)
        out["interp_max_bar_dev_over_2R"] = float(bw / B)
    X, L = Xs[-1], Ls[-1]
    dens = N * 4 / 3 * np.pi * R ** 3 / abs(np.linalg.det(L))
    out["end_density"] = float(dens)
    out["end_density_minus_kepler"] = float(dens - np.pi / np.sqrt(18))
    dens_all = N * 4 / 3 * np.pi * R ** 3 / np.abs(np.linalg.det(Ls))
    out["min_density_along_path"] = float(dens_all.min())
    out["start_density"] = float(dens_all[0])
    r = md.table(X, L)
    cnt = ((r < B * (1 + 1e-6)) & (r > 1e-9)).sum(axis=(1, 2))
    out["end_contacts_per_ball"] = {int(k): int(v) for k, v in zip(*np.unique(cnt, return_counts=True))}
    out["end_next_distance_over_2R"] = float(np.sort(r[r > B * (1 + 1e-6)].ravel())[0] / B)
    out["end_cna"] = {str(k): v for k, v in sorted(cna_signatures(X, L).items())}
    return out


def accepts(rep, interp=True):
    """The acceptance rule used by the certificate and by the must-fail controls."""
    ok = (rep["canonical_bar_err"] < 1e-12 and rep["state0_dev_from_canonical"] < 1e-12
          and rep["lattice0_dev_from_canonical"] < 1e-12
          and rep["max_bar_rel_err"] < 1e-9 and rep["overlap_violations"] == 0
          and rep["max_dot_step"] <= 0.01 and rep["max_lattice_step"] <= 0.01)
    if interp:
        ok = ok and rep["interp_min_nonbonded_over_2R"] > 1 - 1e-4
    return bool(ok)


if __name__ == "__main__":
    import json
    rep = check(sys.argv[1], interp="--no-interp" not in sys.argv)
    print(json.dumps(rep, indent=1, sort_keys=True))
