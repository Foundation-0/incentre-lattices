#!/usr/bin/env python3
"""Independent exact non-overlap check for the explicit one-cell fcc motion.

This module uses only the standard library and imports neither curve.py nor the
certificate verifier. It constructs each distance polynomial by convolution,
checks its Bernstein coefficients in Q(sqrt(2)), and reconstructs the power
polynomial from the Bernstein basis. Every check remains active under python -O.

The finite shell is justified by the analytic image bound in PROOF.md. This is
an internal independent AI check, not outside human review or formalization.
"""

from fractions import Fraction
from itertools import product
import json
from math import comb


# A pair (a, b) represents a + b sqrt(2), with exact rational coefficients.
def field(a=0, b=0):
    return Fraction(a), Fraction(b)


def add(left, right):
    return left[0] + right[0], left[1] + right[1]


def multiply(left, right):
    a, b = left
    c, d = right
    return a * c + 2 * b * d, a * d + b * c


def scale(value, rational):
    return value[0] * rational, value[1] * rational


def sign(value):
    """Return the exact sign; isolate opposing signs before squaring."""
    a, b = value
    if a == 0:
        return (b > 0) - (b < 0)
    if b == 0:
        return (a > 0) - (a < 0)
    if a > 0 and b > 0:
        return 1
    if a < 0 and b < 0:
        return -1

    comparison = (a * a > 2 * b * b) - (a * a < 2 * b * b)
    return comparison if a > 0 else -comparison


def require(condition, message):
    if not condition:
        raise ValueError(message)


def dot(left, right):
    result = field()
    for x, y in zip(left, right):
        result = add(result, multiply(x, y))
    return result


def check_sign_controls():
    # Pell approximants exercise cancellation and both opposing-sign branches.
    controls = (
        (field(), 0),
        (field(1), 1),
        (field(-1), -1),
        (field(0, 1), 1),
        (field(0, -1), -1),
        (field(1, 1), 1),
        (field(-1, -1), -1),
        (field(577, -408), 1),
        (field(-577, 408), -1),
        (field(1393, -985), -1),
        (field(-1393, 985), 1),
    )
    for value, expected in controls:
        require(sign(value) == expected, f"Incorrect exact sign for {value}")


def power_coefficients(n, v, connecting_matrix):
    """Construct ||A0 + A1*s + A2*s^2||^2 - 2(1+s^2)^2."""
    vectors = (
        [field(2 * n[j] + v[j], n[j]) for j in range(3)],
        [
            field(2 * sum(n[k] * connecting_matrix[k][j] for k in range(3)))
            for j in range(3)
        ],
        [field(2 * n[j] + v[j], -n[j]) for j in range(3)],
    )
    coefficients = [field() for _ in range(5)]
    for i, left in enumerate(vectors):
        for j, right in enumerate(vectors):
            coefficients[i + j] = add(coefficients[i + j], dot(left, right))
    for degree, constant in ((0, -2), (2, -4), (4, -2)):
        coefficients[degree] = add(coefficients[degree], field(constant))
    return coefficients


def bernstein_coefficients(coefficients):
    degree = 4
    result = []
    for k in range(degree + 1):
        value = field()
        for j in range(k + 1):
            weight = Fraction(comb(k, j), comb(degree, j))
            value = add(value, scale(coefficients[j], weight))
        result.append(value)
    return result


def reconstruct_power_coefficients(bernstein):
    """Expand binom(4,k)*s^k*(1-s)^(4-k) directly."""
    degree = 4
    result = [field() for _ in range(degree + 1)]
    for k, value in enumerate(bernstein):
        for offset in range(degree - k + 1):
            weight = comb(degree, k) * comb(degree - k, offset) * (-1) ** offset
            result[k + offset] = add(result[k + offset], scale(value, weight))
    return result


def endpoint_coset_count(motif):
    endpoint = ((2, -1, -1), (1, 2, -1), (1, -1, 2))
    inverse = (
        (Fraction(1, 4), Fraction(1, 4), Fraction(1, 4)),
        (Fraction(-1, 4), Fraction(5, 12), Fraction(1, 12)),
        (Fraction(-1, 4), Fraction(1, 12), Fraction(5, 12)),
    )

    # Verify the displayed inverse in both orders before using coset witnesses.
    for left, right in ((endpoint, inverse), (inverse, endpoint)):
        for i in range(3):
            for j in range(3):
                entry = sum(left[i][k] * right[k][j] for k in range(3))
                require(entry == int(i == j), "Endpoint inverse does not invert the lattice")

    residues = {
        tuple(sum(p[k] * inverse[k][j] for k in range(3)) % 1 for j in range(3))
        for p in motif
    }
    require(len(residues) == 6, "Endpoint motif points are not six distinct lattice cosets")
    return len(residues)


def check():
    """Check every finite distance case and return a deterministic integer report."""
    check_sign_controls()
    connecting_matrix = ((0, -1, -1), (1, 0, -1), (1, -1, 0))
    motif = [
        tuple(direction if j == axis else 0 for j in range(3))
        for axis in range(3)
        for direction in (-1, 1)
    ]
    differences = sorted(
        {tuple(p[j] - r[j] for j in range(3)) for p in motif for r in motif}
    )
    translations = sorted(
        n for n in product(range(-3, 4), repeat=3) if sum(x * x for x in n) <= 11
    )
    wider_shell = {
        n for n in product(range(-4, 5), repeat=3) if sum(x * x for x in n) <= 11
    }
    # Outside the wider cube, a coordinate already has square >= 25 > 11.
    require(set(translations) == wider_shell, "Finite translation shell is incomplete")
    require(len(translations) == 171, "Unexpected finite translation count")
    require(len(differences) == 19, "Unexpected motif-difference count")

    polynomials = 0
    positive = 0
    zero = 0
    identically_zero = 0
    for n in translations:
        for v in differences:
            if n == (0, 0, 0) and v == (0, 0, 0):
                continue
            polynomials += 1
            coefficients = power_coefficients(n, v, connecting_matrix)
            bernstein = bernstein_coefficients(coefficients)
            require(
                reconstruct_power_coefficients(bernstein) == coefficients,
                f"Bernstein conversion does not reconstruct the polynomial: n={n}, v={v}",
            )
            signs = [sign(value) for value in bernstein]
            require(
                all(value >= 0 for value in signs),
                f"Negative Bernstein coefficient: n={n}, v={v}, coefficients={bernstein}",
            )
            positive += signs.count(1)
            zero += signs.count(0)
            if all(value == 0 for value in signs):
                identically_zero += 1

    require(polynomials == 3248, "Unexpected distance-polynomial count")
    require(identically_zero == 18, "Unexpected identically-zero polynomial count")
    require(positive + zero == 5 * polynomials, "Coefficient totals are inconsistent")
    return {
        "distance_polynomials": polynomials,
        "positive_coefficients": positive,
        "zero_coefficients": zero,
        "identically_zero_polynomials": identically_zero,
        "distinct_endpoint_cosets": endpoint_coset_count(motif),
    }


if __name__ == "__main__":
    print(json.dumps(check(), sort_keys=True, separators=(",", ":")))
