# A referee's reading of the elastic theorem and the fold

Dated 15 November 2026. Head read: `2d374ff`. Internal review. The checks are in [`results/independent-checks`](../results/independent-checks/README.md).
Marks: **[F]** fact checked by the programs there, **[I]** inference.

## What was checked

* [F] The elastic constants of `results/rigidity-nd` for n = 2 to 9, by separate exact code over Q(sqrt 2): C44 = 0, C11/C12 = (n-1)^2,
  nu_100 = 1/(n^2-n-1), and (C11 - C12)/2 = 2(n-2)/(n 2^n). Status: VERIFIED for n = 2 to 9.
* [F] The fold of `results/fold-continuous-fcc`: no overlap for 6,516 pairs of translation and motif difference over the whole interval,
  a separate far-field bound, strict increase of density. Status: VERIFIED.
* [I] The written proof of the general-n step reads as complete, given two assumptions it did not state. Both are now stated in
  `results/rigidity-nd/PROOF.md`: the bars carry no stress at rest, and the constants are the relaxed ones, taken after the joints settle.

## One wording error, now corrected

The repository said every shear modulus is zero. [F] That holds only at n = 2. C44 = 0 for every n, but the second shear modulus,
(C11 - C12)/2, equals 2(n-2)/(n 2^n): 0 at n = 2, 1/12 at n = 3, 1/16 at n = 4, and positive for every n above 2. The lattice resists
the shear that changes one axis against another. It does not resist the shears that tilt the cube axes.

## Why C44 = 0 does not contradict C12 = C44

Central forces between points often force C12 = C44 (the Cauchy relation). [F] Before the joints settle, this lattice obeys that rule
exactly for every n from 2 to 9 (the `affine` values in `out_elastic_n2_9.txt`: A12 = A44). [F] The rule needs every point to be a centre of
symmetry. No joint here is one, so the joints move when the lattice is strained. Each cross-polytope turns as a rigid body. That removes
C44 and lowers C12. At n = 3, C12 falls from 3/2 - sqrt 2 = 0.0858 to 1/18. So the zero is a property of the relaxed constants,
and the rule still holds for the unrelaxed ones.

## Zero modes and self-stresses

* [F] Each joint has 2n-1 bars, one fewer than the 2n that counting would need for rigidity.
* [F] At q = 0 the independent code finds n(n+1)/2 zero modes and n(n-1)/2 self-stresses per cell for n = 2 to 9 (3 and 1 at n = 2, 6 and 3
  at n = 3, 45 and 36 at n = 9).
* [F] The relaxed stiffness tensor has rank n for n >= 3 (rank 1 at n = 2). At n = 3 that is rank 3 in a six-dimensional strain space, so three strain directions
  cost no energy. At n = 2 only the bulk modulus is non-zero.
* The generic-q counts, n(n-1)/2 zero modes and n(n-3)/2 self-stresses, are proved in `PROOF.md` and checked by `rigidity_nd.py`.
  The independent code did not recompute them.

## Closest earlier source

[F] Borcea and Streinu 2021, Section 7 (see `results/rigidity-nd/LITERATURE.md`): the same bar framework in any dimension, at a different
bar length, with a sketched fold. They give no spring constants and no proof that balls stay apart.

[I] No earlier report was found of C11/C12 = (n-1)^2 or of nu_100 = 1/(n^2-n-1). The derivation is short, so a specialist may call it a worked
example. [I] No earlier report was found of an explicit bar-preserving motion of this framework to the densest packing with an exact proof
of no overlap. Paths such as the Bain path move bare balls with no bars. Goodwin 2006 (full text), Milton's book and Teshima and Matsumoto 2012
were not read in full. The wording is "no earlier report found". It is not a claim of priority.

## Next steps

1. Ask a specialist in periodic frameworks whether the relaxed constants are the ones to quote.
2. Offer the fold to the authors of the 2021 paper as an explicit, certified case of their sketch.
