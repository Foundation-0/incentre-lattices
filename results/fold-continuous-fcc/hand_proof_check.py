#!/usr/bin/env python3
"""Symbolic check of the identities used in HAND-PROOF.md (sympy).

Replayed by `make verify-fold-continuous-hand` (so by `make verify`), but not hashed into certificate.json:
the certified proof of the result is PROOF.md with its exact Bernstein replay. This script only checks
the equalities of Sections 3 and 4 of the written, symmetry-reduced argument that its header lists (not
the block form of G in the rotated basis) and the 48-case centraliser of C used in Section 1, and compares
the reduced conditions with the direct minimum over all 19 motif differences on a parameter grid.
Every check remains active under python -O.
Run from the repository root:  python results/fold-continuous-fcc/hand_proof_check.py
"""
import itertools
import sympy as sp

a, t = sp.symbols("a t", real=True)
d = 2 + a
C = sp.Matrix([[0, -1, -1], [1, 0, -1], [1, -1, 0]])
L = d * sp.eye(3) + t * C
RED = {a**2: 2 - 2 * t**2}          # a^2 + 2 t^2 = 2 on the curve


def require(cond, msg):
    if not cond:
        raise AssertionError(msg)


def red(expr):
    """Reduce modulo a^2 = 2 - 2t^2 (expressions here are at most quadratic in a)."""
    e = sp.expand(expr)
    for _ in range(3):
        e = sp.expand(e.subs(a**3, a * (2 - 2 * t**2)).subs(a**2, 2 - 2 * t**2))
    return e


def same(x, y):
    r = red(x - y)
    require(r == 0, (x, y, r))


# Section 1: of the 48 signed permutations g (e_i g = eps_i e_pi(i)), only +-I and +-m commute with C
gs = [sp.Matrix(3, 3, lambda i, j: sg[i] * int(p[i] == j))
      for p in itertools.permutations(range(3)) for sg in itertools.product((1, -1), repeat=3)]
M = sp.Matrix([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
require(len({tuple(g) for g in gs}) == 48, "48 signed permutations")
cent = {tuple(g) for g in gs if g * C == C * g}
require(cent == {tuple(x) for x in (sp.eye(3), -sp.eye(3), M, -M)}, "centraliser of C")

# Gram matrix G = L L^T and the eigenvalue bound lambda_min >= 7/3
G = L * L.T
E = sp.Matrix([[0, 0, 0], [0, 0, -2], [0, -2, 0]])
J = sp.ones(3, 3)
require(sp.expand(G - ((d**2 + t**2) * sp.eye(3) + t**2 * J + d * t * E)) == sp.zeros(3, 3), "Gram matrix")
same(d**2, 6 + 4 * a - 2 * t**2)
alpha, gamma, beta2 = d**2 + 2 * t**2, d**2 + 3 * t**2 - 2 * d * t, 2 * t**4
same(alpha, 6 + 4 * a)
same(gamma, 6 + 4 * a + t**2 - 4 * t - 2 * a * t)
require(sp.Rational(35, 3) > 6 + 4 * sp.sqrt(2), "35/3 > 6 + 4 sqrt2")

# the nine classes: n, claimed w, claimed |w|^2, and the claimed slack decompositions
# (subtracted term, slack), first for condition (A), then for condition (B)
classes = [
    ((1, 0, 0), (d, -t, -t), 6 + 4*a,
     [(4*d - 2, 0)], [(2*(d + t), 2*(1 - t) + 2*a)]),
    ((0, 1, 0), (t, d, -t), 6 + 4*a,
     [(4*d - 2, 0)], [(2*(d + t), 2*(1 - t) + 2*a)]),
    ((1, 1, 0), (d + t, d - t, -2*t), 12 + 8*a + 2*t**2,
     [(4*(d + t) - 2, 2*(1 - t)**2 + 4 + 4*a)], [(2*(2*d + 2*t), 2*(1 - t)**2 + 2 + 4*a)]),
    ((1, -1, 0), (d - t, -(d + t), 0), 12 + 8*a - 2*t**2,
     [(4*(d + t) - 2, 2*(1 - t)*(3 + t) + 4*a)], [(2*(2*d), 4 + 4*a - 2*t**2)]),
    ((0, 1, 1), (2*t, d - t, d - t), 12 + 8*a + 2*t**2 - 8*t - 4*a*t,
     [(4*(2*t) - 2, 2*(1 - t)*(7 - t) + a*(8 - 4*t)),
      (4*(d - t) - 2, 2*(1 - t)**2 + 4 + 4*a*(1 - t))],
     [(2*(2*d - 2*t), 2*(1 - t)**2 + 2 + 4*a*(1 - t)),
      (2*(d + t), 2*(1 - t)*(4 - t) + a*(6 - 4*t))]),
    ((0, 1, -1), (0, d + t, -(d + t)), 2*(d + t)**2,
     [(4*(d + t) - 2, 2*(d + t - 1)**2)], [(2*(2*(d + t)), 2*(d + t)*(d + t - 2))]),
    ((1, 1, 1), (d + 2*t, d - 2*t, d - 2*t), 18 + 12*a + 6*t**2 - 8*t - 4*a*t,
     [(4*(d + 2*t) - 2, 6*(t - sp.Rational(4, 3))**2 + sp.Rational(4, 3) + a*(8 - 4*t))],
     [(2*(2*d), 6*(t - sp.Rational(2, 3))**2 + sp.Rational(22, 3) + 4*a*(2 - t))]),
    ((1, -1, -1), (d - 2*t, -d, -d), 18 + 12*a - 2*t**2 - 8*t - 4*a*t,
     [(4*d - 2, 2*(1 - t)*(5 + t) + 2 + a*(8 - 4*t))], [(2*(2*d), 2*(1 - t)*(5 + t) + a*(8 - 4*t))]),
    ((1, 1, -1), (d, d, -(d + 2*t)), 18 + 12*a - 2*t**2 + 8*t + 4*a*t,
     [(4*(d + 2*t) - 2, 12 + 8*a - 2*t**2 + 4*a*t)], [(2*(2*d + 2*t), 10 + 8*a - 2*t**2 + 4*t + 4*a*t)]),
]
for n, w, norm2, slacks_a, slacks_b in classes:
    wv = sp.Matrix([n]) * L
    require(all(sp.expand(x - y) == 0 for x, y in zip(wv, w)), n)
    same((wv * wv.T)[0], norm2)
    for subtract, slack in slacks_a + slacks_b:
        same(norm2 - subtract, slack)
# |n|^2 = 4: n = 2 e_i
for i in range(3):
    n = sp.Matrix([[2 * int(j == i) for j in range(3)]])
    same(((n * L) * (n * L).T)[0], 4 * (6 + 4 * a))

# the nine classes, with -n and the mirror image n m, cover the 26 vectors with 1 <= |n|^2 <= 3;
# the listed subtracted terms of (A) and (B), as functions of (a, t), for the grid comparison below
candidates = {}
for n, w, norm2, slacks_a, slacks_b in classes:
    f = sp.lambdify((a, t), [[c for c, _ in slacks_a], [c for c, _ in slacks_b]], "math")
    for k in (n, (n[0], n[2], n[1])):
        candidates[k] = candidates[tuple(-x for x in k)] = f
near = {n for n in itertools.product(range(-1, 2), repeat=3) if 1 <= sum(x * x for x in n) <= 3}
require(len(candidates) == 26 and set(candidates) == near, "the nine classes")

# numerical cross-check of the reduction (A), (B), (C) against the direct minimum
import math
P = [tuple(s * int(j == i) for j in range(3)) for i in range(3) for s in (1, -1)]
V = sorted({tuple(p[k] - q[k] for k in range(3)) for p in P for q in P})
worst = math.inf      # all pairs, bars included (bars are at squared clearance 0)
margin = math.inf     # pairs other than bars, 0 <= s < 1
for k in range(201):
    s = k / 200
    tt = 2 * s / (1 + s * s); aa = math.sqrt(2) * (1 - s * s) / (1 + s * s)
    Ln = [[(2 + aa) * (i == j) + tt * int(C[i, j]) for j in range(3)] for i in range(3)]
    nonbar_s = math.inf
    for n in itertools.product(range(-3, 4), repeat=3):
        if n == (0, 0, 0):
            continue
        w = [sum(n[i] * Ln[i][j] for i in range(3)) for j in range(3)]
        m = sorted((abs(x) for x in w), reverse=True)
        nn = sum(x * x for x in w)
        reduced = min(nn - 4 * m[0] + 2, nn - 2 * (m[0] + m[1]), nn - 2)
        vals = sorted(sum((w[j] + v[j]) ** 2 for j in range(3)) - 2 for v in V)
        direct = vals[0]
        require(abs(reduced - direct) < 1e-9, (s, n))
        if sum(x * x for x in n) >= 5:
            require(nn >= 35 / 3 - 1e-12, (s, n))
        if n in candidates:   # the actual 4 m1 - 2 and 2(m1 + m2) never exceed the largest listed candidate
            cand_a, cand_b = candidates[n](aa, tt)
            require(4 * m[0] - 2 <= max(cand_a) + 1e-9 and 2 * (m[0] + m[1]) <= max(cand_b) + 1e-9, (s, n))
        worst = min(worst, direct)
        nonbar_s = min(nonbar_s, vals[1] if sum(x * x for x in n) == 1 else vals[0])   # n = +-e_i: skip the bridge
    if k < 200:
        margin = min(margin, nonbar_s)
    else:
        at_end = nonbar_s
require(worst > -1e-12, worst)
require(margin > 1e-9, margin)
require(abs(at_end) < 1e-9, at_end)
print("hand proof identities: all checked; centraliser of C is {I, -I, m, -m}; grid of 201 values of s: "
      "smallest squared clearance of non-bar pairs", round(margin, 12), "for s < 1 and", round(at_end, 12), "at s = 1")
