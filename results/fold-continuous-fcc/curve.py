#!/usr/bin/env python3
"""Explicit m=1 motion and exact distance certificates. See PROOF.md.

The proof calculation uses only integers and Fraction, in Q(sqrt(2)).
NumPy is imported only by the optional numerical path sampler.
"""
from fractions import Fraction as F
from itertools import product
from math import comb

C = ((0, -1, -1), (1, 0, -1), (1, -1, 0))
P = tuple(tuple(s * int(j == i) for j in range(3))
          for i in range(3) for s in (1, -1))
ZERO = (F(0), F(0))


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def scale(x, a):
    return x[0] * a, x[1] * a


def mul(x, y):
    return x[0] * y[0] + 2 * x[1] * y[1], x[0] * y[1] + x[1] * y[0]


def sign(x):
    """Sign of a+b*sqrt(2), using rational comparisons only."""
    a, b = x
    if not a:
        return (b > 0) - (b < 0)
    if not b or (a > 0) == (b > 0):
        return (a > 0) - (a < 0)
    delta = a * a - 2 * b * b
    return ((delta > 0) - (delta < 0)) * ((a > 0) - (a < 0))


def dot(xs, ys):
    value = ZERO
    for x, y in zip(xs, ys):
        value = add(value, mul(x, y))
    return value


def differences():
    return sorted({tuple(a - b for a, b in zip(p, q)) for p in P for q in P})


def shifts():
    # The proof excludes all |n|^2 >= 12; the remaining coordinates obey |n_i| <= 3.
    return [n for n in product(range(-3, 4), repeat=3) if sum(x*x for x in n) <= 11]


def coefficients(n, v, matrix=C):
    """Power coefficients of (1+s^2)^2 (|n L(s)+v|^2-2)."""
    a0 = [(F(2*x+y), F(x)) for x, y in zip(n, v)]
    a1 = [(F(2*sum(n[j]*matrix[j][i] for j in range(3))), F(0)) for i in range(3)]
    a2 = [(F(2*x+y), F(-x)) for x, y in zip(n, v)]
    return [add(dot(a0, a0), (-2, 0)), scale(dot(a0, a1), 2),
            add(add(dot(a1, a1), scale(dot(a0, a2), 2)), (-4, 0)),
            scale(dot(a1, a2), 2), add(dot(a2, a2), (-2, 0))]


def bernstein(cs):
    """Degree-four Bernstein coefficients, on the entire interval [0,1]."""
    result = []
    for k in range(5):
        b = ZERO
        for j in range(k+1):
            b = add(b, scale(cs[j], F(comb(k, j), comb(4, j))))
        result.append(b)
    return result


def exact_nonoverlap(matrix=C):
    """Exhaust the finite cases left by the written all-image bound.

    This alone is not a proof for other matrices: the all-image bound in PROOF.md
    is established specifically for C. A changed matrix is used only as a rejection control.
    """
    count = identities = coefficient_count = zero_coefficients = 0
    for n in shifts():
        for v in differences():
            if n == (0, 0, 0) and v == (0, 0, 0):
                continue  # A point compared with itself, the sole exclusion.
            cs = coefficients(n, v, matrix)
            bs = bernstein(cs)
            if any(sign(b) < 0 for b in bs):
                raise ValueError(f"negative Bernstein coefficient: n={n}, v={v}, b={bs}")
            count += 1
            identities += all(c == ZERO for c in cs)
            coefficient_count += len(bs)
            zero_coefficients += sum(b == ZERO for b in bs)
    return {"integer_shifts": len(shifts()), "motif_differences": len(differences()),
            "distance_polynomials": count, "bernstein_coefficients": coefficient_count,
            "zero_coefficients": zero_coefficients, "identically_zero_polynomials": identities,
            "all_coefficients_nonnegative_exactly": True, "subdivisions": 0}


def determinant(m):
    a, b, c = m
    return (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0])
            + a[2]*(b[0]*c[1]-b[1]*c[0]))


def endpoint():
    """Integer endpoint checks; the parity/index argument identifying fcc is in PROOF.md."""
    lattice = tuple(tuple(2*int(i == j) + C[i][j] for j in range(3)) for i in range(3))
    contacts = []
    distinct = True
    for p in P:
        count = 0
        for q in P:
            for n in shifts():
                if p == q and n == (0, 0, 0):
                    continue
                v = [q[i]-p[i]+sum(n[j]*lattice[j][i] for j in range(3)) for i in range(3)]
                dsq = sum(x*x for x in v)
                distinct &= dsq > 0
                count += dsq == 2
        contacts.append(count)
    return {"lattice_rows": lattice, "determinant": determinant(lattice),
            "lattice_rows_have_even_sum": all(sum(row) % 2 == 0 for row in lattice),
            "motif_points_have_odd_sum": all(sum(p) % 2 == 1 for p in P),
            "distinct_cosets": distinct, "contacts_per_ball": contacts}


def sample(destination, states=501):
    """Write a derived numerical NPZ for the unchanged existing fold checkers."""
    import numpy as np
    if states < 2:
        raise ValueError("at least two samples are required")
    s = np.linspace(0., 1., states)
    t = 2*s/(1+s*s)
    a = np.sqrt(2.)*(1-s*s)/(1+s*s)
    k = 2-np.sqrt(2.)
    lattice = k*((2+a[:, None, None])*np.eye(3) + t[:, None, None]*np.array(C))
    motif = 1 + k*np.array(P, dtype=float)
    np.savez_compressed(destination, m=np.array(1), Ls=lattice,
                        Xs=np.broadcast_to(motif, (states, 6, 3)).copy())


if __name__ == "__main__":
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", help="write numerical samples to this NPZ path")
    parser.add_argument("--states", type=int, default=501)
    args = parser.parse_args()
    if args.sample:
        sample(args.sample, args.states)
    else:
        print(json.dumps(exact_nonoverlap(), indent=2, sort_keys=True))
