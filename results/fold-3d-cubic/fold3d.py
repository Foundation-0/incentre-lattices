"""3D remark: octahedra of 6 dots, alternate rotation +-theta about [111], lattice kept cubic.

a(theta) = K1 + sqrt(4R^2 - 2 K2^2); exact at 60 degrees: a = 8c/3, six contacts per ball,
density 27 pi sqrt2/256; start density pi R^3.  Trigonal (a,b,b) lattices: the feasible set on
a 1 degree x 0.01 grid, connected to the start, contains denser states.
Engine functions (numpy, scipy, sympy); verify.py is the certificate driver.
"""
import itertools
from collections import deque

import numpy as np
import sympy as sp
from scipy.spatial import cKDTree

R_ = sp.sqrt(2) - 1
c_ = 2 - sp.sqrt(2)
R = float(R_)
c = float(c_)
nv = np.ones(3) / 3 ** 0.5


def rot(t):
    K = np.array([[0, -nv[2], nv[1]], [nv[2], 0, -nv[0]], [-nv[1], nv[0], 0]])
    return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K


octp = np.array([s * c * e for e in np.eye(3) for s in (1, -1)])


def a_cubic(t):
    K2 = 2 * c * (1 - np.cos(t)) / 3
    K1 = 2 * c * np.cos(t) + K2
    return K1 + np.sqrt(4 * R ** 2 - 2 * K2 ** 2)


def cell3(t, V, rng=2):
    P, lab = [], []
    Oa = octp @ rot(t).T
    Ob = octp @ rot(-t).T
    for ijk in itertools.product(range(-rng, rng + 1), repeat=3):
        O = Oa if sum(ijk) % 2 == 0 else Ob
        for m in range(6):
            P.append(np.array(ijk) @ V + O[m])
            lab.append((ijk, m))
    return np.array(P), lab


def cubic_path():
    """First non-bonded contact along the cubic-restricted path and the largest bar deviation."""
    P0, lab = cell3(0, 2 * np.eye(3))
    D0 = np.linalg.norm(P0[:, None] - P0[None], axis=2)
    linked = np.abs(D0 - 2 * R) < 1e-9
    ci = [lab.index(((0, 0, 0), m)) for m in range(6)]
    first, dev = None, 0.0
    for t in np.linspace(0, np.pi / 3 + 0.01, 3000):
        P, _ = cell3(t, a_cubic(t) * np.eye(3))
        D = np.linalg.norm(P[:, None] - P[None], axis=2)
        dev = max(dev, np.abs(D - D0)[linked].max())
        nb = np.where(linked | np.eye(len(P), dtype=bool), np.inf, D)[ci].min()
        if first is None and nb <= 2 * R + 1e-9:
            first = float(np.degrees(t))
    return first, float(dev)


def exact_60():
    """sympy: a(60) = 8c/3, all six contacts of every dot exactly 2R, density closed form."""
    t60 = sp.pi / 3
    n_ = sp.Matrix([1, 1, 1]) / sp.sqrt(3)
    K = sp.Matrix([[0, -n_[2], n_[1]], [n_[2], 0, -n_[0]], [-n_[1], n_[0], 0]])
    rotS = lambda t: sp.eye(3) + sp.sin(t) * K + (1 - sp.cos(t)) * K * K  # noqa: E731
    K2s = 2 * c_ * (1 - sp.cos(t60)) / 3
    K1s = 2 * c_ * sp.cos(t60) + K2s
    a60 = sp.simplify(K1s + sp.sqrt(4 * R_ ** 2 - 2 * K2s ** 2))
    P60, lab60 = cell3(np.pi / 3, float(a60) * np.eye(3))
    D60 = np.linalg.norm(P60[:, None] - P60[None], axis=2)
    ci = [lab60.index(((0, 0, 0), m)) for m in range(6)]
    Oa = [rotS(t60) * (s * c_ * sp.eye(3)[:, i]) for i in range(3) for s in (1, -1)]
    Ob = [rotS(-t60) * (s * c_ * sp.eye(3)[:, i]) for i in range(3) for s in (1, -1)]
    exact_ok, contacts, nxt = True, [], []
    for m in range(6):
        i0 = ci[m]
        close = np.where((D60[i0] < 2 * R + 1e-6) & (D60[i0] > 1e-9))[0]
        contacts.append(len(close))
        for j in close:
            ijk, mm = lab60[j]
            O = Oa if sum(ijk) % 2 == 0 else Ob
            v = sp.Matrix(ijk) * a60 + O[mm] - Oa[m]
            exact_ok &= sp.simplify(sp.sqrtdenest(sp.simplify(v.dot(v))) - 4 * R_ ** 2) == 0
        nxt.append(float(np.sort(D60[i0][D60[i0] > 2 * R + 1e-6])[0] / (2 * R)))
    dens = sp.simplify(6 * sp.Rational(4, 3) * sp.pi * R_ ** 3 / a60 ** 3)
    return {"a60_minus_8c_over_3": sp.simplify(a60 - 8 * c_ / 3), "contacts": contacts, "next_over_2R": nxt,
            "all_contacts_exact": bool(exact_ok), "a60": float(a60),
            "density_minus_closed_form": sp.N(dens - 27 * sp.pi * sp.sqrt(2) / 256, 30),
            "density": float(sp.N(dens, 15)), "start_density": float(sp.N(sp.pi * R_ ** 3, 15))}


def coordination_sequence(a60, terms=10):
    Pbig, labbig = cell3(np.pi / 3, a60 * np.eye(3), rng=6)
    tree = cKDTree(Pbig)
    adj = [[] for _ in range(len(Pbig))]
    for i, j in tree.query_pairs(2 * R + 1e-7):
        adj[i].append(j)
        adj[j].append(i)
    start = labbig.index(((0, 0, 0), 0))
    dist = {start: 0}
    dq = deque([start])
    while dq:
        v = dq.popleft()
        if dist[v] >= terms:
            continue
        for w in adj[v]:
            if w not in dist:
                dist[w] = dist[v] + 1
                dq.append(w)
    return [sum(1 for v in dist.values() if v == k) for k in range(1, terms + 1)]


def trigonal_scan():
    """Grid scan of trigonal lattices V = [[a,b,b],[b,a,b],[b,b,a]], b = K2 + u, bars kept exact."""
    shifts = list(itertools.product(range(-2, 3), repeat=3))
    bardev = [0.0]

    def feas(t, u):
        K2 = 2 * c * (1 - np.cos(t)) / 3
        K1 = 2 * c * np.cos(t) + K2
        if 4 * R ** 2 - 2 * u ** 2 < 0:
            return None
        b = K2 + u
        a = K1 + np.sqrt(4 * R ** 2 - 2 * u ** 2)
        V = np.array([[a, b, b], [b, a, b], [b, b, a]])
        Oa = octp @ rot(t).T
        Ob = octp @ rot(-t).T
        P = np.array([np.array(s) @ V + (Oa if sum(s) % 2 == 0 else Ob)[m] for s in shifts for m in range(6)])
        ctr = [shifts.index((0, 0, 0)) * 6 + m for m in range(6)]
        D = np.linalg.norm(P[ctr][:, None] - P[None], axis=2)
        D[np.arange(6), ctr] = np.inf
        bar = max(abs(np.linalg.norm(V[i] + Ob[2 * i + 1] - Oa[2 * i]) - 2 * R) for i in range(3))
        bardev[0] = max(bardev[0], bar)
        return abs(np.linalg.det(V)), D.min()

    ts = np.radians(np.arange(0, 90.01, 1.0))
    us = np.arange(-0.6, 0.6001, 0.01)
    F = {}
    for i, t in enumerate(ts):
        for j, u in enumerate(us):
            r = feas(t, u)
            if r is not None and r[1] >= 2 * R - 1e-7:
                F[(i, j)] = 6 * 4 / 3 * np.pi * R ** 3 / r[0]
    j0 = int(np.argmin(abs(us)))
    seen = {(0, j0)}
    dq = deque([(0, j0)])
    while dq:
        i, j = dq.popleft()
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            k = (i + di, j + dj)
            if k in F and k not in seen:
                seen.add(k)
                dq.append(k)
    best = max(seen, key=lambda k: F[k])
    return {"feasible": len(F), "connected": len(seen), "best_density": float(F[best]),
            "best_theta_deg": float(np.degrees(ts[best[0]])), "best_u": float(us[best[1]]), "bar_dev": bardev[0]}
