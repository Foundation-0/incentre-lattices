"""Framework (b): the whole drawing = contacts + red lines (dot to its pyramid's vertices)
+ cube edges + pyramid lateral edges (centre to cube vertices). Also (c): pyramids/triangles alone."""
import numpy as np, itertools
from pframe import *
from contact_nd import contact_framework, TSTAR

def drawing(n, with_contacts=True, with_red=True, with_pyr=True):
    A, X, B = contact_framework(n)
    X = list(X); B = [] if not with_contacts else list(B)
    c = len(X); X.append(np.zeros(n))                    # cube centre
    v = len(X); X.append(np.ones(n))                     # cube vertex (1,..,1); others are lattice translates
    # cube vertex at sign vector sg = (1,..,1) + lattice offset (sg-1)/2 cells
    def vert(sg):
        m = tuple(int((s - 1) // 2) for s in sg)        # sg=-1 -> -1 cell
        return v, m
    z = (0,) * n
    # cube edges: vertex (1..1) to vertex (1..1)-2e_i
    for i in range(n):
        m = [0] * n; m[i] = -1
        B.append((v, v, tuple(m)))
    if with_pyr:
        for sg in itertools.product((1, -1), repeat=n):
            vi, m = vert(sg); B.append((c, vi, m))
    if with_red:
        for i in range(n):
            for s in (1, -1):
                d = 2 * i + (0 if s == 1 else 1)
                B.append((d, c, z))
                for sg in itertools.product((1, -1), repeat=n):
                    if sg[i] == s:
                        vi, m = vert(sg); B.append((d, vi, m))
    return A, np.array(X), B

rng = np.random.default_rng(3)
for n in (2, 3, 4):
    for label, kw in [("(b) full drawing", {}), ("(c) pyramid skeleton only (cube+pyramid edges)", dict(with_contacts=False, with_red=False)),
                      ("dots + red lines + cube/pyramid edges, no contacts", dict(with_contacts=False))]:
        A, X, B = drawing(n, **kw)
        if not kw.get('with_red', True):
            X = X[2 * n:]; B = [(i - 2 * n, j - 2 * n, m) for (i, j, m) in B]
        g = counts(A, X, B, random_q(n, rng)); z0 = counts(A, X, B)
        nm, sr, _, ss = flexible_analysis(A, X, B)
        print(f"n={n} {label}: sites={len(X)} bars={len(B)} generic(N0,S)={g} q0(N0,S)={z0} flex mech={nm} strainDim={sr} SS={ss}")
