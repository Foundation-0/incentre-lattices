"""Contact framework of the truncated n-cubic honeycomb (dots at +-t e_i in a cube of side 2).
Bonds: cross-polytope edges inside each cube (all pairs except antipodal) + cross-face bonds.
t = 2 - sqrt(2) is the incentre / Wythoff point (all edges equal 2R, R = sqrt2 - 1)."""
import numpy as np
from pframe import *

TSTAR = 2 - np.sqrt(2)


def contact_framework(n, t=TSTAR, internal=True):
    A = 2.0 * np.eye(n)
    X = []
    for i in range(n):
        for s in (+1, -1):
            v = np.zeros(n); v[i] = s * t; X.append(v)
    X = np.array(X)  # site 2i = +t e_i, 2i+1 = -t e_i
    z = (0,) * n
    bonds = []
    if internal:
        for a in range(2 * n):
            for b in range(a + 1, 2 * n):
                if a // 2 != b // 2:
                    bonds.append((a, b, z))
    for i in range(n):
        m = [0] * n; m[i] = 1
        bonds.append((2 * i, 2 * i + 1, tuple(m)))
    return A, X, bonds


if __name__ == "__main__":
    rng = np.random.default_rng(1)
    print("n | sites bonds | Maxwell 2n^2-bonds | generic-q zero modes, self-stress | q=0 zero, SS | n(n-1)/2  n(n-3)/2")
    for n in range(2, 7):
        A, X, B = contact_framework(n)
        e = [np.linalg.norm(bond_vec(A, X, b)) for b in B]
        assert np.allclose(e, 2 * (np.sqrt(2) - 1))
        g = [counts(A, X, B, random_q(n, rng)) for _ in range(3)]
        assert len(set(g)) == 1
        z0 = counts(A, X, B)
        print(n, "|", len(X), len(B), "|", n * len(X) - len(B), "|", g[0], "|", z0, "|", n * (n - 1) // 2, n * (n - 3) // 2)
    # special q planes: q_1 = 0, others generic
    print("\nzero modes with k of the q-components set to 0 (others generic):")
    for n in range(2, 6):
        A, X, B = contact_framework(n)
        row = []
        for k in range(n + 1):
            q = random_q(n, rng); q[:k] = 0
            row.append(counts(A, X, B, q))
        print(n, row)
    # t dependence (dot moved along the axis, unequal edge lengths)
    print("\nt-dependence of generic-q and q=0 counts (n=2,3):")
    for n in (2, 3):
        for t in (0.2, 0.4, TSTAR, 0.7, 0.9):
            A, X, B = contact_framework(n, t)
            print(n, round(t, 4), counts(A, X, B, random_q(n, rng)), counts(A, X, B))
