"""Finite, floating-point probe of independently rotating cab bridges.

Source: Foundation-0/incentre-lattices at commit
2d374fff2a59cc9331566733043a4bb25c1040bf,
results/fold-continuous-fcc/PROOF.md. This is not a proof or a body model.
All lengths use the proof's normalized units. Sphere diameter is sqrt(2).
"""
import itertools
import json
from pathlib import Path

import numpy as np

C = np.array([[0., -1., -1.], [1., 0., -1.], [1., -1., 0.]])
P = np.concatenate([np.eye(3), -np.eye(3)])
V = np.unique((P[:, None] - P[None, :]).reshape(-1, 3), axis=0)
N = np.array(list(itertools.product(range(-3, 4), repeat=3)))
ns = np.repeat(N, len(V), axis=0)
vs = np.tile(V, (len(N), 1))
self_pair = np.all(ns == 0, axis=1) & np.all(vs == 0, axis=1)
internal = np.all(ns == 0, axis=1) & (np.sum(vs * vs, axis=1) == 2)
bridge = (np.sum(ns * ns, axis=1) == 1) & np.all(vs == -2 * ns, axis=1)
variable = ~(self_pair | internal | bridge)


def measure(degrees):
    angles = np.deg2rad(degrees)
    lattice = np.diag(2 + np.sqrt(2) * np.cos(angles)) + np.sin(angles)[:, None] * C
    bridge_vectors = lattice - 2 * np.eye(3)
    assert np.max(np.abs(np.sum(bridge_vectors ** 2, axis=1) - 2)) < 1e-12
    assert np.linalg.det(lattice) > 0
    singular_min = np.linalg.svd(lattice, compute_uv=False)[-1]
    # Outside [-3,3]^3, |n| >= 4. Every motif difference has |v| <= 2.
    # A strict positive margin rules out omitted overlaps for this sample,
    # up to floating-point accuracy of the singular value calculation.
    outside_lower_bound = 4 * singular_min - 2
    assert outside_lower_bound > np.sqrt(2) + 1e-8
    distances = np.linalg.norm(ns @ lattice + vs, axis=1)
    ids = np.flatnonzero(variable)
    closest = ids[np.argmin(distances[variable])]
    active = ids[np.abs(distances[variable] - np.sqrt(2)) < 1e-8]
    gap = float(distances[closest] - np.sqrt(2))
    return {
        'angles_deg': list(map(float, degrees)),
        'minimum_nonbar_gap': gap,
        'feasible_at_tolerance_1e_8': gap >= -1e-8,
        'directed_active_nonbar_difference_classes': len(active),
        'nearest_class': {'translation': ns[closest].astype(int).tolist(),
                          'motif_difference': vs[closest].astype(int).tolist()},
        'outside_distance_lower_bound': float(outside_lower_bound),
        'determinant': float(np.linalg.det(lattice)),
    }


def main():
    diagonal = [measure([angle] * 3) for angle in [0, 15, 30, 45, 60, 75, 89, 90]]
    cubes = []
    for center, delta in [(30, 5), (45, 5), (60, 5), (89, 1)]:
        rows = [measure([center + delta * i for i in offset])
                for offset in itertools.product([-1, 0, 1], repeat=3)]
        cubes.append({'center_deg': center, 'offset_step_deg': delta,
                      'samples': len(rows),
                      'feasible': sum(r['feasible_at_tolerance_1e_8'] for r in rows),
                      'minimum_gap': min(r['minimum_nonbar_gap'] for r in rows)})
    near_fcc = [measure(a) for a in [[89, 90, 90], [90, 89, 90], [90, 90, 89],
                                   [89, 89, 90], [89, 90, 89], [90, 89, 89],
                                   [85, 90, 90], [90, 85, 90], [90, 90, 85]]]
    grid = [measure(a) for a in itertools.product(range(0, 91, 15), repeat=3)]
    rng = np.random.default_rng(20261112)
    random_samples = [measure(a) for a in rng.uniform(0, 90, size=(300, 3))]
    broad_samples = grid + random_samples
    worst = min(broad_samples, key=lambda r: r['minimum_nonbar_gap'])
    result = {
        'evidence': 'Finite floating-point samples, no continuous feasibility claim.',
        'source_commit': '2d374fff2a59cc9331566733043a4bb25c1040bf',
        'model': 'Rigid regular octahedron, spherical sites, independently rotated fixed-length bridges.',
        'omitted_physics': 'No finite connector thickness, joint limits, stability or human anatomy.',
        'exact_counterexample': {
            'angles_deg': [0, 90, 90],
            'lattice_rows': [['2+sqrt(2)', 0, 0], [1, 2, -1], [1, -1, 2]],
            'translation': [-1, 1, 1], 'motif_difference': [1, -1, 0],
            'separation': ['1-sqrt(2)', 0, 1],
            'squared_distance': '4-2*sqrt(2) < 2 = squared sphere diameter',
            'pair_threshold': 'With angles 2 and 3 at 90 degrees, this pair overlaps iff 0 <= angle 1 < arcsin(1/3).',
            'threshold_degrees': float(np.rad2deg(np.arcsin(1 / 3))),
            'limit': 'This pair clearing above the threshold does not establish feasibility of the whole packing.',
        },
        'diagonal': diagonal, 'independent_angle_samples': cubes, 'near_fcc': near_fcc,
        'full_angle_box_probe': {
            'grid_samples': len(grid), 'random_samples': len(random_samples),
            'random_seed': 20261112,
            'feasible': sum(r['feasible_at_tolerance_1e_8'] for r in broad_samples),
            'worst_sample': worst,
            'minimum_outside_distance_lower_bound': min(r['outside_distance_lower_bound'] for r in broad_samples),
            'violations': [r for r in broad_samples if not r['feasible_at_tolerance_1e_8']],
        },
    }
    destination = Path(__file__).with_name('vitruvian_motion_probe.json')
    assert not measure([0, 90, 90])['feasible_at_tolerance_1e_8']
    assert measure([90, 90, 90])['feasible_at_tolerance_1e_8']
    destination.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
