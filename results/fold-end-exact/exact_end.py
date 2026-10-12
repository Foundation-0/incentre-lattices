#!/usr/bin/env python3
"""Exact check of a fold end state (numpy for I/O only; every verdict is integer arithmetic).

The final state of a path is rotated into the cube frame of octahedron 0 and written in units
u = K/den, K = 2 - sqrt2, for the smallest den <= 12 that makes every coordinate an integer
(to 1e-7). In those units the ball diameter is 2R = sqrt2 K, so (2R)^2 = 2 den^2 exactly, and the
whole check is on integers:

  * every bar of the cab bar graph (built here, fixed integer lattice shifts) has |r|^2 = 2 den^2;
  * no pair of balls, over all periodic images, has |r|^2 < 2 den^2 (no overlap);
  * contacts (|r|^2 = 2 den^2) per ball;
  * density pi/sqrt18 exactly  <=>  |det L| = 2 N K^3  <=>  |det L_u| = 2 N den^3;
  * common-neighbour analysis on the exact contact graph: (4,2,1) for fcc, (4,2,2) for the hcp twin pairs.

Usage: python exact_end.py PATH.npz      (prints a JSON report; exit 1 if any exact check fails)
"""
import collections
import itertools
import json
import sys
from fractions import Fraction

import numpy as np

K = 2.0 - np.sqrt(2.0)


def bar_graph(m):
    """(i, j, s): bar from dot i to dot j + s @ L (s integer). Dot = 6*cell + 2*axis + (0: +e, 1: -e)."""
    bars = []
    for ix, iy, iz in itertools.product(range(m), repeat=3):
        c = (ix * m + iy) * m + iz
        for a, b in itertools.combinations(range(3), 2):
            for sa, sb in itertools.product(range(2), repeat=2):
                bars.append((6 * c + 2 * a + sa, 6 * c + 2 * b + sb, (0, 0, 0)))
        idx = [ix, iy, iz]
        for a in range(3):
            nb = list(idx)
            nb[a] = (nb[a] + 1) % m
            s = [0, 0, 0]
            s[a] = 1 if idx[a] + 1 == m else 0
            bars.append((6 * c + 2 * a, 6 * ((nb[0] * m + nb[1]) * m + nb[2]) + 2 * a + 1, tuple(s)))
    return bars


def exactify(X, L, tol=1e-7, den_max=12):
    """Integer coordinates (dots Y, lattice rows Lu) in units K/den, cube frame of octahedron 0."""
    c0 = X[:6].mean(axis=0)
    Q = np.array([(X[2 * a] - X[2 * a + 1]) / (2 * K) for a in range(3)])
    if np.abs(Q @ Q.T - np.eye(3)).max() > 1e-9:
        raise ValueError("octahedron 0 is not a regular octahedron of half-diagonal K")
    if np.linalg.det(Q) < 0:
        raise ValueError("octahedron 0 frame is reflected")
    Yf = (X - c0) @ Q.T / K
    Lf = L @ Q.T / K
    for den in range(1, den_max + 1):
        Y, Lu = np.rint(Yf * den), np.rint(Lf * den)
        if np.abs(Yf * den - Y).max() < tol and np.abs(Lf * den - Lu).max() < tol:
            return [[int(v) for v in r] for r in Y], [[int(v) for v in r] for r in Lu], den
    raise ValueError("end state is not on the (K/den) integer grid for den <= %d" % den_max)


def det3(A):
    return (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
            - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
            + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]))


def check(path):
    z = np.load(path, allow_pickle=False)
    m = int(z["m"])
    Y, Lu, den = exactify(z["Xs"][-1], z["Ls"][-1])
    out = check_exact(m, Y, Lu, den)
    out["path"] = str(path)
    return out


def check_exact(m, Y, Lu, den):
    """All verdicts on an integer state: dots Y (6 m^3 x 3) and lattice rows Lu, in units K/den."""
    N = 6 * m ** 3
    D2 = 2 * den * den
    sq = lambda v: sum(x * x for x in v)
    add = lambda p, s: tuple(p[k] + sum(s[a] * Lu[a][k] for a in range(3)) for k in range(3))
    sub = lambda p, q: tuple(p[k] - q[k] for k in range(3))
    out = {"m": m, "balls": N, "unit": "K/%d, K = 2 - sqrt2" % den,
           "lattice_rows": Lu, "dots": Y}

    bars = bar_graph(m)
    out["bars"] = len(bars)
    out["bars_exact"] = all(sq(sub(add(Y[j], s), Y[i])) == D2 for i, j, s in bars)

    det = det3(Lu)
    out["det_L"] = det
    out["density_is_pi_over_sqrt18_exactly"] = abs(det) == 2 * N * den ** 3

    # all pairs within the contact distance, over periodic images. Floats are used only to wrap the dots
    # into the cell and to bound the image shifts (|shift_a| <= rcut * |Linv[:, a]| + 1, plus 1 spare);
    # every distance is then an exact integer.
    Linv_f = np.linalg.inv(np.array(Lu, float))
    F = np.array(Y, float) @ Linv_f
    W = np.floor(F).astype(int)
    Yw = [tuple(Y[i][k] - sum(int(W[i][a]) * Lu[a][k] for a in range(3)) for k in range(3)) for i in range(N)]
    rcut = np.sqrt(D2) * 1.01
    nmax = [int(np.ceil(rcut * np.linalg.norm(Linv_f[:, a]))) + 2 for a in range(3)]
    shifts = list(itertools.product(*[range(-n, n + 1) for n in nmax]))
    nbr = collections.defaultdict(set)       # i -> {(j, s)} with s relative to wrapped coordinates
    min_d2 = None
    for i in range(N):
        for j in range(N):
            for s in shifts:
                if i == j and not any(s):
                    continue
                d2 = sq(sub(add(Yw[j], s), Yw[i]))
                if min_d2 is None or d2 < min_d2:
                    min_d2 = d2
                if d2 == D2:
                    nbr[i].add((j, s))
    out["no_overlap_exact"] = min_d2 >= D2
    out["min_pair_distance_over_2R_squared"] = str(Fraction(min_d2, D2))
    deg = [len(nbr[i]) for i in range(N)]
    out["contacts_per_ball"] = dict(collections.Counter(deg))

    def nb_of(j, s):
        return {(k, tuple(a + b for a, b in zip(t, s))) for k, t in nbr[j]}

    def pos(j, s):
        return add(Yw[j], s)

    cna = collections.Counter()
    for i in range(N):
        Ni = {(j, s) for j, s in nbr[i]}
        for j, s in Ni:
            com = list(Ni & nb_of(j, s))
            bonds = [(p, q) for p, q in itertools.combinations(com, 2) if sq(sub(pos(*p), pos(*q))) == D2]
            ends = [v for e in bonds for v in e]
            cna[(len(com), len(bonds), 2 if len(ends) != len(set(ends)) else 1)] += 1
    out["cna_pairs"] = {"(%d,%d,%d)" % k: v for k, v in sorted(cna.items())}
    kinds = set(cna)
    out["structure"] = ("fcc" if kinds == {(4, 2, 1)} else
                        "hcp" if kinds == {(4, 2, 1), (4, 2, 2)} and cna[(4, 2, 1)] == cna[(4, 2, 2)] else
                        "close-packed polytype" if kinds <= {(4, 2, 1), (4, 2, 2)} else "other")
    out["exact_ok"] = bool(out["bars_exact"] and out["no_overlap_exact"] and out["density_is_pi_over_sqrt18_exactly"]
                           and set(deg) == {12})
    return out


if __name__ == "__main__":
    rep = check(sys.argv[1])
    print(json.dumps(rep, indent=1, default=str))
    sys.exit(0 if rep["exact_ok"] else 1)
