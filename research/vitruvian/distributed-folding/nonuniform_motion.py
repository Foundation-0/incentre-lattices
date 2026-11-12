"""Exact finite clearance certificate for separated, nonuniform fold phases.

The theorem uses row-vector increments chosen independently along each
material-coordinate line: L_i(k)=2e_i+sqrt(2)cos(theta_i(k))e_i+sin(theta_i(k))C_i.
Every theta_i(k) may lie within 5 degrees of 45 degrees. Positions are sums
of these increments along the three coordinates. No dynamics is claimed.
"""
import itertools
import json
import math
from pathlib import Path

C = ((0, -1, -1), (1, 0, -1), (1, -1, 0))
P = [tuple(s if j == i else 0 for j in range(3)) for i in range(3) for s in [-1, 1]]
V = sorted({tuple(a-b for a, b in zip(p, q)) for p in P for q in P})
assert len(P) == len(set(P)) == 6 and len(V) == 19


def sign_quadratic(a, b):
    """Exact sign of a+b*sqrt(2) for integer coefficients."""
    if not a:
        return (b > 0) - (b < 0)
    if not b or (a > 0) == (b > 0):
        return (a > 0) - (a < 0)
    square_difference = a*a - 2*b*b
    return ((square_difference > 0) - (square_difference < 0)) * ((a > 0) - (a < 0))


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def check():
    checks = []
    for n in itertools.product([-1, 0, 1], repeat=3):
        if n == (0, 0, 0):
            continue  # Same-motif distances are fixed, at least sqrt(2).
        z = tuple(sum(n[i]*C[i][j] for i in range(3)) for j in range(3))
        steps = sum(abs(x) for x in n)
        for v in V:
            if steps == 1 and all(v[i] == -2*n[i] for i in range(3)):
                continue  # Prescribed bridge, exactly sqrt(2) even when phase varies.
            u = tuple(3*n[i] + v[i] for i in range(3))
            # 64*(|nL0+v|^2 - (sqrt(2)+|n|_1/8)^2), L0=3I+C/sqrt(2).
            a = 64*dot(u, u) + 32*dot(z, z) - 128 - steps*steps
            b = 64*dot(u, z) - 16*steps
            assert sign_quadratic(a, b) > 0, (n, v, a, b)
            checks.append({'translation': n, 'motif_difference': v, 'a': a, 'b': b})
    # Symmetric part of L0 has minimum eigenvalue 3-1/sqrt(2).
    # For |n| >= 2, use sigma0 > 9/4 and sqrt(3)/8 < 1/4.
    # Then |X(n+m)-X(n)+v| > (9/4-1/4)|m|-2 >= 2 > sqrt(2).
    # The first two rational bounds follow from sqrt(2)>4/3 and sqrt(3)<2.
    assert 2 > (4/3)**2 and 3 < 2**2
    # Every phase deviation <=5deg has chord displacement
    # 2sqrt(2)sin(5deg/2) <= sqrt(2)*pi/36 < 1/8.
    # Prove the last comparison with pi < 22/7 and integer arithmetic.
    assert 2 * (22 * 8)**2 < (7 * 36)**2
    # Sanity checks for the exact sign function, including near cancellation.
    assert sign_quadratic(1, -1) < 0 and sign_quadratic(-1, 1) > 0
    assert sign_quadratic(3, -2) > 0 and sign_quadratic(-3, 2) < 0
    assert len(checks) == 26*19-6 == 488
    weakest = min(checks, key=lambda q: q['a'] + math.sqrt(2)*q['b'])
    result = {
        'statement': 'Any separated angle sequences theta_i(k) in [40deg,50deg] define a globally non-overlapping, bar-preserving infinite framework through integrated cell increments.',
        'evidence': 'Written analytic construction and bounds plus exact integer-pair certificate for nearby cell translations. Not a formal proof-assistant result.',
        'nearby_distance_conditions': len(checks),
        'all_nearby_comparisons_strictly_positive': True,
        'weakest_comparison': weakest,
        'row_perturbation_upper_bound': '1/8',
        'far_translation_argument': '|m|>=2 => distance > (9/4-1/4)|m|-2 >=2>sqrt(2)',
        'positions': 'X(n)=sum_i S_i(n_i), where S_i(k+1)-S_i(k)=L_i(k)',
        'example': 'theta_1(k,t)=45deg+5deg*sin(2*pi*k/m-omega*t), theta_2=theta_3=45deg, any positive integer m.',
        'scope': 'Kinematic traveling modulation, fixed motif orientations; no actuation, material stability, forces, finite rods, or autonomous wave dynamics established.',
        'conditions': checks,
    }
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'conditions'}, indent=2))


if __name__ == '__main__':
    check()
