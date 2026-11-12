"""An exact-form local obstruction to a naive spatial wave of common fold phase.

No dynamics or motif rotations are included. This is a closure calculation,
not a claim that the framework cannot support nonuniform deformation.
"""
import json
from pathlib import Path

import numpy as np

C = np.array([[0., -1., -1.], [1., 0., -1.], [1., -1., 0.]])


def cell_rows(degrees):
    theta = np.deg2rad(degrees)
    return (2 + np.sqrt(2) * np.cos(theta)) * np.eye(3) + np.sin(theta) * C


records = []
for base in [0, 20, 45, 70]:
    delta = 10
    # A common local phase theta(n)=base+delta*n_1 changes in direction 1.
    # On the (1,2) square, paths 1 then 2 vs 2 then 1 differ by this vector.
    mismatch = cell_rows(base)[0] + cell_rows(base + delta)[1] - cell_rows(base)[0] - cell_rows(base)[1]
    expected_norm = 2 * np.sqrt(2) * np.sin(np.deg2rad(delta) / 2)
    assert abs(np.linalg.norm(mismatch) - expected_norm) < 1e-12
    records.append({'base_degrees': base, 'phase_step_degrees': delta,
                    'closure_error_vector': mismatch.tolist(),
                    'closure_error_norm': float(np.linalg.norm(mismatch))})

result = {
    'statement': 'With fixed motif orientation, a common phase increasing by delta along one cell coordinate fails the transverse square-loop closure condition.',
    'exact_norm': '2*sqrt(2)*abs(sin(delta/2))',
    'derivation': 'The closure error is L_2(theta+delta)-L_2(theta), the chord of a radius-sqrt(2) circle separated by angle delta.',
    'examples': records,
    'limits': 'Does not rule out nonuniform motions with independently assigned phases, motif rotations, strain, defects, or changed connections. No collision or dynamics conclusion.',
}
Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
