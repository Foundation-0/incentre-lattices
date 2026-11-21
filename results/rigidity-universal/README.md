# Universal rigidity of the cross-polytope and the periodic framework

This folder adds a short note and an exact check. It asks one question. Can a positive semidefinite stress (a sign condition on bar tensions that certifies a framework cannot be bent into a higher dimension) exist for the cross-polytope frameworks of `results/rigidity-nd`?

Label: PROVEN (written proof in [PROOF.md](PROOF.md)) and VERIFIED (finite exact check, listed below). Every piece of mathematics here is known. The incentre instance in part 4 is, as far as a bounded search found, not written down elsewhere. That is not a priority claim.

## What the note shows

1. Maxwell's identity (1870). For every self-stress of a finite bar framework, the sum over bars of tension times squared length is zero. So the trace of the stress matrix is twice the sum of the tensions.
2. If every stressed bar has the same length and the framework has at least d + 2 vertices in R^d, the only positive semidefinite stress is zero. So no super-stability certificate exists.

   A framework is super stable when a positive semidefinite stress of the right rank exists, plus a mild condition on bar directions. Connelly's theorem then gives universal rigidity.

3. More strongly, a finite framework with all bars of one length and at least d + 2 vertices is not universally rigid. Send vertex k to the k-th unit vector of R^V, scaled to the bar length. Every bar keeps its length and the map is not a congruence.

   This covers the finite cross-polytope for every n >= 2.
4. With the lattice held fixed, the periodic incentre framework behaves differently. At zero wavevector there is an explicit stress whose quadratic form is a sum of squares with kernel exactly the constant displacements.

   Then every periodic configuration with the same lattice and the same bar lengths is a translate. This holds for every n >= 2 and every 0 < c < 1. It agrees with Theorem 5.8 of Dewar, Schulze, Tanigawa and Theran (arXiv:2510.19529).

## What is not shown

* The case where the lattice may deform. A remark in PROOF.md says what the stress identities give there. It is not a theorem.
* Whether the stress matrix stays positive semidefinite at every wavevector.
* Competitors that are periodic only on a larger cell, or not periodic at all.
* Anything for the exact value c = 2 - sqrt(2) by computation. The script uses rational c. The proof works for every real c in (0, 1), so it covers 2 - sqrt(2).

Open question: with the lattice fixed, is the Bloch stress matrix of this stress positive semidefinite for every wavevector, with kernel only at zero wavevector on the constants? If yes, the infinite framework resists every bounded displacement, not only periodic ones.

## What the script checks

Exact arithmetic only (`fractions`). No floating point and no tolerance. Seed 20261121 is printed.

| Part | Check | Range |
|---|---|---|
| A | Maxwell's identity on random rational frameworks | d = 1..4, 60 stresses |
| B | Cross-polytope stress space has dimension max(0, n(n-3)/2), the explicit family spans it, Omega = -(W kron J2), every nonzero stress is indefinite | n = 2..7 |
| C | The simplex map keeps every bar length and changes the antipodal distance from 2 to sqrt(2) | n = 2..7 |
| D | Periodic stress is in equilibrium, S = T, Omega(0) is positive semidefinite of rank 2n-1 with kernel the constants, equals the sum of squares, satisfies the fixed-lattice energy identity, and the framework has n(n-1)/2 nontrivial zero modes at q = 0 | n = 2..6, c = 1/2, 1/3, 2/5, 3/4 |

The symbols S, T, W and Omega(0) are defined in PROOF.md.

## Controls that must fail

Each control runs the same check function as the real input on a corrupted input. The run stops if the corrupted input is accepted.

* A wrong bridging tension (two different shifts, n = 3 and n = 4).
* One intra bar with its sign flipped.
* The whole stress negated. It balances, and the positive semidefinite test rejects it.
* A balanced stress with one negative weight.
* A stress on a disconnected weight graph. The kernel is larger than the constants.
* A corrupted Omega, rejected by the sum-of-squares match and by the energy identity.
* Maxwell's identity on a stress that is not in equilibrium.
* A simplex image with one vertex moved, the identity image, and the regular tetrahedron (V = d + 1, where the simplex map is a congruence).
* Unequal bars: the unit square with both diagonals has a nonzero positive semidefinite stress. The equal-bar obstruction does not apply, as it should not.

## Run

```
python3 -I -B results/rigidity-universal/universal_rigidity_check.py
python3 -I -O -B results/rigidity-universal/universal_rigidity_check.py
python3 -I -B -W error results/rigidity-universal/universal_rigidity_check.py
python3 -I -B results/rigidity-universal/universal_rigidity_check.py --corrupt   # must exit 1
python3 -I -B results/rigidity-universal/test_universal_rigidity.py
python3 -I -O -B results/rigidity-universal/test_universal_rigidity.py
python3 -I -B -W error results/rigidity-universal/test_universal_rigidity.py
python3 tools/sha256check.py results/rigidity-universal/SHA256SUMS
```

The check takes about one second. Standard library only. The three modes print identical output.

## Where this sits in the repository

This folder is outside `make verify`. The script is named `universal_rigidity_check.py`, not `verify.py`, because `tools/build_certificates.py` requires a certificate for every `results/*/verify.py`. It has no certificate and is not in `commitment/`.

`SHA256SUMS` binds the files. A maintainer who wants it in CI can add a Makefile target that runs the two commands above. That is a change to an existing file and is left out here.

The open drafts for the exact replay of the rigidity counts and for the folds touch other folders. This note touches neither. The fold drafts study motions where the lattice changes. This note shows that, with the lattice fixed, every periodic configuration that keeps every bar length is a translate.

## AI assistance

Claude Sonnet 5.5 (Anthropic) wrote this note, the script and the tests, and re-derived each formula before writing it. Claude Fable 5.2 directed the work. Separate Claude Opus 5.5 runs checked the results independently and reviewed this folder. They reached the same results. No human has reviewed the proofs.
