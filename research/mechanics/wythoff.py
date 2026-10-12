"""Wythoff frameworks for affine Coxeter groups: rigidity census.
For each group and each ringing (set of mirrors NOT containing the generating point) we build
the periodic bar framework whose bars are the Wythoff edges {g p, g s_i p}, i ringed,
and count zero modes / self-stresses."""
import numpy as np
import itertools
from pframe import *

r2, r3 = np.sqrt(2), np.sqrt(3)


def mir(a, c):
    a = np.asarray(a, float)
    nrm = np.linalg.norm(a)
    return (a / nrm, c / nrm)


def D_lattice(n):
    B = [np.eye(n)[0] + np.eye(n)[1]]
    for i in range(1, n):
        B.append(np.eye(n)[i - 1] * 0 + np.eye(n)[0] - np.eye(n)[i] if i == 1 else np.eye(n)[i - 1] + np.eye(n)[i])
    # basis of D_n: e1+e2, e1-e2, e2+e3, ..., e_{n-1}+e_n
    B = [np.eye(n)[0] + np.eye(n)[1], np.eye(n)[0] - np.eye(n)[1]] + [np.eye(n)[i] + np.eye(n)[i + 1] for i in range(1, n - 1)]
    return np.array(B).T


def group_C(n):
    E = np.eye(n)
    ms = [mir(E[n - 1], 0)] + [mir(E[i] - E[i + 1], 0) for i in range(n - 2, -1, -1)] + [mir(E[0], 1)]
    return dict(name=f"C~{n}", mirrors=ms, L=2 * np.eye(n))


def group_B(n):
    E = np.eye(n)
    ms = [mir(E[n - 1], 0)] + [mir(E[i] - E[i + 1], 0) for i in range(n - 2, -1, -1)] + [mir(E[0] + E[1], 1)]
    return dict(name=f"B~{n}", mirrors=ms, L=D_lattice(n))


def group_D(n):
    E = np.eye(n)
    ms = [mir(E[n - 2] + E[n - 1], 0)] + [mir(E[i] - E[i + 1], 0) for i in range(n - 2, -1, -1)] + [mir(E[0] + E[1], 1)]
    return dict(name=f"D~{n}" if n > 3 else "A~3", mirrors=ms, L=D_lattice(n))


def group_G2():
    ms = [mir((0, 1), 0), mir((-0.5, r3 / 2), 0), mir((1, 0), 1)]
    return dict(name="G~2", mirrors=ms, L=np.array([[2, 0], [1, r3]]).T)


def group_A2():
    ms = [mir((0, 1), 0), mir((-r3 / 2, 0.5), 0), mir((r3 / 2, 0.5), r3 / 2)]
    return dict(name="A~2", mirrors=ms, L=np.array([[0, r3], [1.5, r3 / 2]]).T)


def reflect_map(a, c):
    n = len(a)
    M = np.eye(n) - 2 * np.outer(a, a)
    v = 2 * c * a
    return M, v  # x -> M x + v


class Reducer:
    def __init__(self, L):
        self.L = L; self.Li = np.linalg.inv(L)

    def reduce(self, x):
        f = self.Li @ x
        fl = np.floor(f + 1e-9)
        fr = f - fl
        fr[np.abs(fr - 1) < 1e-9] = 0
        return fr, fl.astype(int)

    def key(self, x):
        fr, _ = self.reduce(x)
        return tuple(np.round(fr, 7))


def cosets(G):
    """BFS over chambers g.F modulo translations; returns list of affine maps (M,v)."""
    ms = G["mirrors"]; red = Reducer(G["L"])
    S = [reflect_map(a, c) for a, c in ms]
    n = len(ms[0][0])
    probe = np.array([0.1234, 0.0567, 0.0089, 0.0013, 0.0002][:n]) * 0 + None if False else None
    # key on image of a generic interior point of the fundamental simplex
    verts = simplex_vertices(G)
    gp = verts.mean(0) + 1e-3 * np.arange(1, n + 1) / n
    seen = {}
    frontier = [(np.eye(n), np.zeros(n))]
    seen[red.key(gp)] = frontier[0]
    while frontier:
        nf = []
        for (M, v) in frontier:
            for (Ms, vs) in S:  # g * s
                M2 = M @ Ms; v2 = M @ vs + v
                k = red.key(M2 @ gp + v2)
                if k not in seen:
                    seen[k] = (M2, v2); nf.append((M2, v2))
        frontier = nf
    return list(seen.values())


def simplex_vertices(G):
    ms = G["mirrors"]; n = len(ms[0][0])
    V = []
    for drop in range(n + 1):
        rows = [ms[i] for i in range(n + 1) if i != drop]
        Aa = np.array([r[0] for r in rows]); cc = np.array([r[1] for r in rows])
        V.append(np.linalg.solve(Aa, cc))
    return np.array(V)


def wythoff_point(G, ringed, weights=None):
    """point on unringed mirrors with distances to ringed mirrors = weights (default equal)."""
    ms = G["mirrors"]; n = len(ms[0][0])
    cen = simplex_vertices(G).mean(0)
    rows, rhs = [], []
    ringed = list(ringed)
    w = np.ones(len(ringed)) if weights is None else np.asarray(weights, float)
    # unknowns p (n) and d (1)
    for i, (a, c) in enumerate(ms):
        sg = np.sign(a @ cen - c)
        if i in ringed:
            rows.append(np.r_[sg * a, -w[ringed.index(i)]]); rhs.append(sg * c)
        else:
            rows.append(np.r_[a, 0.0]); rhs.append(c)
    # n+1 equations, n+1 unknowns
    sol = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)[0]
    return sol[:n]


def build(G, ringed, p=None, weights=None):
    ms = G["mirrors"]; n = len(ms[0][0]); red = Reducer(G["L"]); L = G["L"]
    if p is None:
        p = wythoff_point(G, ringed, weights)
    gs = cosets(G)
    S = [reflect_map(a, c) for a, c in ms]
    sites = {}; X = []

    def site(x):
        k = red.key(x)
        if k not in sites:
            fr, _ = red.reduce(x)
            sites[k] = len(X); X.append(L @ fr)
        return sites[k]

    bonds = set(); blist = []
    for (M, v) in gs:
        x = M @ p + v
        for i in ringed:
            Ms, vs = S[i]
            y = M @ (Ms @ p + vs) + v
            ix, iy = site(x), site(y)
            _, cx = red.reduce(x); _, cy = red.reduce(y)
            m = tuple(int(t) for t in (cy - cx))
            key = (ix, iy, m) if (ix, m) <= (iy, tuple(-t for t in m)) else (iy, ix, tuple(-t for t in m))
            key2 = min((ix, iy, m), (iy, ix, tuple(-t for t in m)))
            if key2 not in bonds:
                bonds.add(key2); blist.append(key2)
    X = np.array(X)
    return L.astype(float), X, blist, len(gs)


def census(G, rng, maxsites=2000, randomise=False):
    n = len(G["mirrors"][0][0])
    rows = []
    for k in range(1, n + 2):
        for ringed in itertools.combinations(range(n + 1), k):
            w = rng.uniform(0.5, 1.5, len(ringed)) if randomise else None
            A, X, B, P = build(G, ringed, weights=w)
            if len(X) > maxsites:
                continue
            lens = np.array([np.linalg.norm(bond_vec(A, X, b)) for b in B])
            if lens.min() < 1e-6:
                rows.append((ringed, "degenerate")); continue
            g = counts(A, X, B, random_q(n, rng))
            z = counts(A, X, B)
            nm, sr, _, ss = flexible_analysis(A, X, B)
            V = len(X); E = len(B)
            rows.append(dict(ring=ringed, V=V, E=E, z=round(2 * E / V, 3), maxwell=n * V - E,
                             generic=g, q0=z, flex=(nm, sr, ss), equal=np.allclose(lens, lens[0])))
    return rows


if __name__ == "__main__":
    import sys
    rng = np.random.default_rng(7)
    groups = [group_C(2), group_G2(), group_A2(), group_C(3), group_B(3), group_D(3),
              group_C(4), group_B(4), group_D(4)]
    for G in groups:
        print("==", G["name"], " |W/T| =", len(cosets(G)))
        for rand in (False, True):
            for r in census(G, rng, randomise=rand):
                if isinstance(r, tuple):
                    print("  ", r); continue
                tag = "generic-pt" if rand else "uniform  "
                print(f"   {tag} ring={r['ring']} V={r['V']} E={r['E']} z={r['z']} nV-E={r['maxwell']} "
                      f"genericq(N0,S)={r['generic']} q0(N0,S)={r['q0']} flex(mech,strainDim,SS)={r['flex']} eqlen={r['equal']}")
        sys.stdout.flush()
