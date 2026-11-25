# Corrections to the version 2 preprint

This unissued erratum is for review against the version 2 PDF, SHA-256 `632d7debf3dc06a55e18de70dfaee90da3b664cea9b5be31de8a02599f72a003`. It preserves the published PDF, LaTeX source, DOI and version metadata. Publication as an erratum or incorporation into a later version remains a maintainer decision.

## Reproduction scope

`make verify` replays the registered results, their finite checks and exact certificates. It does not rerun the searches that generated the logged intermediate jam densities. Research scripts outside the registered targets have their own commands and dependencies. The additions since v1 include the continuous fcc motion, one-cell sampled paths, exact end states and the CaB6 attainment argument.

## Rigidity and elasticity

The count formulas and the stress and rotation bases are for n >= 3; the finite count checks cover n = 3 through 8, with n = 2 recorded separately. The elastic calculations also include n = 2. The replay does not compute variable-lattice counts.

The periodic comparison to Borcea and Streinu concerns the same infinite framework at a different bar length. Counts in this paper are per cell of 2Z^n. In the quadratic energy reduction, the matrix is M = nI + J and the linear term is -2(2e) dot y. The stiffness statement C44 = 0 refers to cube-axis shear components, and nu_100 refers to uniaxial stress along a cube axis.

## Folding evidence

The 60-degree cubic-path contact limit was located numerically on a grid. The endpoint spacing and density have exact expressions; that no non-bonded pair touches earlier was checked numerically, not proved by the replay.

In two dimensions, 2 sqrt2 c cos(theta + pi/4) is the gap of the pair that closes. Opposite dots in one square can be closer at earlier angles, while retaining a constant positive gap. The first-contact assertion beyond this pair remains a numerical all-pairs check.

The sampled-path end states have rational closed forms in units K for fcc and K/3 for hcp. A cluster is unrotated in its own chosen frame by definition; it can rotate in the laboratory frame. The fcc endpoint is a translate of the stated parity lattice, without an additional translation being silently applied twice.

The continuous-motion point-group converse needs the bridge argument: edges within octahedra lie in triangles, while bridges do not. A symmetry therefore preserves both classes, has lattice translation and signed-permutation linear part g, and obeys Bg = gB. For B = aI + tC with t > 0, this reduces to Cg = gC. Only the four stated signed permutations commute with C. Independent verification of this reduction remains part of review of the pending hand-proof correction.

## Radius ratio and CaB6 scope

The H1 family list needs its mirror-wall justification: regular polytopes use ridge mirrors; bipyramids use equatorial and vertical mirror planes; the stated Catalan examples inherit the appropriate dual-edge bisectors. The UB12 identification is a numerical research-script comparison, outside `make verify`. The cuboctahedron is a fourth counterexample control, not a facet-transitive example. The CaB6 positional parameter is taken in 0 < x < 1/2; x and 1-x describe the same structure.

## Credit and links

The literature comparison also used GPT-6.1 Sol. The record does not establish which proofs were personally checked by the author, so the acknowledgment should state that the recorded reviews were by AI models under direction and leave human checking scope unasserted. The current source repository is https://github.com/Foundation-0/incentre-lattices.

A later issued revision needs a rebuilt and visually checked PDF, resolved references, updated citation metadata, a new PDF hash and an exact-revision review record. This file supplies corrections for review and does not certify that those release steps occurred.
