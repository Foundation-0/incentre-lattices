#!/usr/bin/env python3
"""Symbolic check of the identities used in HAND-PROOF.md (sympy).

Not part of `make verify` and not hashed into certificate.json: the certified proof of the result
is PROOF.md with its exact Bernstein replay. This script only checks that every algebraic identity
quoted in the written, symmetry-reduced argument is correct, and compares the reduced conditions
with the direct minimum over all 19 motif differences on a parameter grid.
Run from the repository root:  python results/fold-continuous-fcc/hand_proof_check.py
"""
import itertools
import sympy as sp

a, t = sp.symbols("a t", real=True)
d = 2 + a
C = sp.Matrix([[0, -1, -1], [1, 0, -1], [1, -1, 0]])
L = d * sp.eye(3) + t * C
RED = {a**2: 2 - 2 * t**2}          # a^2 + 2 t^2 = 2 on the curve


def red(expr):
    """Reduce modulo a^2 = 2 - 2t^2 (expressions here are at most quadratic in a)."""
    e = sp.expand(expr)
    for _ in range(3):
        e = sp.expand(e.subs(a**3, a * (2 - 2 * t**2)).subs(a**2, 2 - 2 * t**2))
    return e


def same(x, y):
    assert red(x - y) == 0, (x, y, red(x - y))


# Gram matrix G = L L^T and the eigenvalue bound lambda_min >= 7/3
G = L * L.T
E = sp.Matrix([[0, 0, 0], [0, 0, -2], [0, -2, 0]])
J = sp.ones(3, 3)
assert sp.expand(G - ((d**2 + t**2) * sp.eye(3) + t**2 * J + d * t * E)) == sp.zeros(3, 3)
same(d**2, 6 + 4 * a - 2 * t**2)
alpha, gamma, beta2 = d**2 + 2 * t**2, d**2 + 3 * t**2 - 2 * d * t, 2 * t**4
same(alpha, 6 + 4 * a)
same(gamma, 6 + 4 * a + t**2 - 4 * t - 2 * a * t)
assert sp.Rational(35, 3) > 6 + 4 * sp.sqrt(2)

# the nine classes: n, claimed w, claimed |w|^2, and the claimed slack decompositions
classes = [
    ((1, 0, 0), (d, -t, -t), 6 + 4*a,
     [(4*d - 2, 0), (2*(d + t), 2*(1 - t) + 2*a)]),
    ((0, 1, 0), (t, d, -t), 6 + 4*a,
     [(4*d - 2, 0), (2*(d + t), 2*(1 - t) + 2*a)]),
    ((1, 1, 0), (d + t, d - t, -2*t), 12 + 8*a + 2*t**2,
     [(4*(d + t) - 2, 2*(1 - t)**2 + 4 + 4*a), (2*(2*d + 2*t), 2*(1 - t)**2 + 2 + 4*a)]),
    ((1, -1, 0), (d - t, -(d + t), 0), 12 + 8*a - 2*t**2,
     [(4*(d + t) - 2, 2*(1 - t)*(3 + t) + 4*a), (2*(2*d), 4 + 4*a - 2*t**2)]),
    ((0, 1, 1), (2*t, d - t, d - t), 12 + 8*a + 2*t**2 - 8*t - 4*a*t,
     [(4*(2*t) - 2, 2*(1 - t)*(7 - t) + a*(8 - 4*t)),
      (4*(d - t) - 2, 2*(1 - t)**2 + 4 + 4*a*(1 - t)),
      (2*(2*d - 2*t), 2*(1 - t)**2 + 2 + 4*a*(1 - t)),
      (2*(d + t), 2*(1 - t)*(4 - t) + a*(6 - 4*t))]),
    ((0, 1, -1), (0, d + t, -(d + t)), 2*(d + t)**2,
     [(4*(d + t) - 2, 2*(d + t - 1)**2), (2*(2*(d + t)), 2*(d + t)*(d + t - 2))]),
    ((1, 1, 1), (d + 2*t, d - 2*t, d - 2*t), 18 + 12*a + 6*t**2 - 8*t - 4*a*t,
     [(4*(d + 2*t) - 2, 6*(t - sp.Rational(4, 3))**2 + sp.Rational(4, 3) + a*(8 - 4*t)),
      (2*(2*d), 6*(t - sp.Rational(2, 3))**2 + sp.Rational(22, 3) + 4*a*(2 - t))]),
    ((1, -1, -1), (d - 2*t, -d, -d), 18 + 12*a - 2*t**2 - 8*t - 4*a*t,
     [(4*d - 2, 2*(1 - t)*(5 + t) + 2 + a*(8 - 4*t)), (2*(2*d), 2*(1 - t)*(5 + t) + a*(8 - 4*t))]),
    ((1, 1, -1), (d, d, -(d + 2*t)), 18 + 12*a - 2*t**2 + 8*t + 4*a*t,
     [(4*(d + 2*t) - 2, 12 + 8*a - 2*t**2 + 4*a*t), (2*(2*d + 2*t), 10 + 8*a - 2*t**2 + 4*t + 4*a*t)]),
]
for n, w, norm2, slacks in classes:
    wv = sp.Matrix([n]) * L
    assert all(sp.expand(x - y) == 0 for x, y in zip(wv, w)), n
    same((wv * wv.T)[0], norm2)
    for subtract, slack in slacks:
        same(norm2 - subtract, slack)
# |n|^2 = 4: n = 2 e_i
for i in range(3):
    n = sp.Matrix([[2 * int(j == i) for j in range(3)]])
    same(((n * L) * (n * L).T)[0], 4 * (6 + 4 * a))

# numerical cross-check of the reduction (A), (B), (C) against the direct minimum
import math
P = [tuple(s * int(j == i) for j in range(3)) for i in range(3) for s in (1, -1)]
V = sorted({tuple(p[k] - q[k] for k in range(3)) for p in P for q in P})
worst = math.inf
for k in range(201):
    s = k / 200
    tt = 2 * s / (1 + s * s); aa = math.sqrt(2) * (1 - s * s) / (1 + s * s)
    Ln = [[(2 + aa) * (i == j) + tt * int(C[i, j]) for j in range(3)] for i in range(3)]
    for n in itertools.product(range(-3, 4), repeat=3):
        if n == (0, 0, 0):
            continue
        w = [sum(n[i] * Ln[i][j] for i in range(3)) for j in range(3)]
        m = sorted((abs(x) for x in w), reverse=True)
        nn = sum(x * x for x in w)
        reduced = min(nn - 4 * m[0] + 2, nn - 2 * (m[0] + m[1]), nn - 2)
        direct = min(sum((w[j] + v[j]) ** 2 for j in range(3)) for v in V) - 2
        assert abs(reduced - direct) < 1e-9
        if sum(x * x for x in n) >= 5:
            assert nn >= 35 / 3 - 1e-12
        worst = min(worst, direct)
assert worst > -1e-12
print("hand proof identities: all checked; grid minimum of squared clearance", round(worst, 12))
