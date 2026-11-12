# Distributed folding with exact clearance

This note constructs spatially varying folds of the bridged octahedron framework. Each angle sequence depends on one lattice coordinate. Integrating the three sequences gives a compatible infinite framework with fixed bar lengths and no sphere overlap for angles between 40 and 50 degrees.

**PROVEN, with internal AI review:** the written construction, distant-pair bound, and 488 exact nearby-pair comparisons establish this restricted family. Historical priority and outside mathematical review remain open.

## Reproduce

From the repository root, with the dependencies in `requirements.txt` installed:

```sh
python research/vitruvian/distributed-folding/nonuniform_motion.py
python research/vitruvian/distributed-folding/posture_closure.py
python research/vitruvian/distributed-folding/vitruvian_motion_probe.py
```

Use ordinary Python execution with assertions enabled. Each script writes JSON beside its source. The first script is an exact integer-arithmetic certificate. The other two record finite calculations and failed broader prescriptions.

## Evidence

- [Written proof](nonuniform_motion_proof.md): closure, bar lengths, all-pairs clearance, and a traveling kinematic modulation.
- [Exact certificate](nonuniform_motion.json): all 488 strict comparisons, with the integer sign rule in the accompanying Python source.
- [Closure calculation](posture_closure.json): assigning the same spatially varying phase to all three directions generally fails to close a lattice square.
- [Independent-angle probe](vitruvian_motion_probe.json): the full independent-angle box contains overlaps. The restricted interval theorem does not assert otherwise.

The proof permits arbitrary sequences in its stated interval. The global periodic cell may change with time. It does not assert a self-propagating wave, fixed boundaries, actuation, material stability, or hardware clearance.

## Credit and prior art

The starting framework and uniform motion are in [the continuous-fold proof](../../../results/fold-continuous-fcc/PROOF.md).

Borcea and Streinu, *Auxetic deformations and elliptic curves* (2021), Section 7, is related prior work on bridged cross-polytopes: https://journals.sagepub.com/doi/abs/10.1177/1081286520948116 . General rigid-unit flexibility is established. No priority claim is made for this construction.

Ahraz Arifuddin proposed exploring a three-dimensional Vitruvian interpretation of the lattice. GPT-6.1 Sol developed the construction, arithmetic replay, and supporting probes. GPT-6 Astra independently reviewed the proof. No human mathematical verification is recorded. Code is MIT. Text is CC BY 4.0, as in the repository.
