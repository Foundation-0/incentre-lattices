# Dimensional rigidity and elasticity: literature comparison

Audit date: 2026-10-12. Result reviewed: `2ed2ab7`, [PROOF.md](PROOF.md).
This is an internal literature audit assisted by GPT-6.1 Sol, not independent human review.

**The bridged cross-polytope family is prior art.** Borcea and Streinu explicitly describe it in arbitrary dimension.
The earlier family-level **NOT FOUND** wording is withdrawn. Priority of this repository's exact formulas remains
unresolved; this audit establishes neither novelty nor an error in those formulas. The proof and finite replay have
different evidential roles, and a successful replay cannot establish priority.

## Primary sources and reading depth

Only primary articles or their institutional records support the comparisons below. Each entry separates the
source's content from the assessment made here. Access limits are part of the record.

### Goodwin (2006): linearly bridged rigid units

A. L. Goodwin, *Rigid unit modes and intrinsic flexibility in linearly bridged framework structures*,
Physical Review B **74**, 134302 (2006), [doi:10.1103/PhysRevB.74.134302](https://doi.org/10.1103/PhysRevB.74.134302).
[Oxford institutional record and abstract](https://ora.ox.ac.uk/objects/uuid%3A0cef4bf2-97c6-4259-bb51-b4b3d807069d).

**Read:** abstract and metadata through the indexed institutional record; full text unavailable in this pass.
The abstract gives a method for counting rigid unit modes in linearly bridged frameworks, an extensive density of
modes across reciprocal space, and examples including Prussian blue, Zn(CN)2 and extended beta-quartz.

**Assessment:** the three-dimensional bridged-unit flexibility phenomenon is known in substance. The abstract
does not establish the exact mode multiplicities, self-stress basis or arbitrary-n elastic formulas in PROOF.md.
Those comparisons require the full paper; they cannot be treated as a completed full-text search.

### Power (2014): the symbol and mode multiplicity

S. C. Power, *Polynomials for crystal frameworks and the rigid unit mode spectrum*,
Philosophical Transactions of the Royal Society A **372**, 20120030 (2014),
[doi:10.1098/rsta.2012.0030](https://doi.org/10.1098/rsta.2012.0030),
[open article at PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3871295/).

**Read:** full-text rendering, focusing on Sections 2, 4 and 5 and the fixed/flexible-lattice distinction in Section 3.
Equation images were not independently transcribed. Section 2 defines a periodic motif and Maxwell counting equilibrium;
Section 4 constructs the matrix symbol. Theorem 5.1 identifies its phase evaluation with the phase-periodic rigidity
map. Definition 5.2 distinguishes the RUM spectrum from mode multiplicity. Section 5(b) counts supercell modes by
summing phase nullities. Section 3 separates fixed-lattice modes from affinely periodic flexes.

**Assessment:** this supplies established machinery for the Bloch calculation and its interpretation. Merely citing
it does not establish that this particular framework's n-dependent multiplicities or elastic constants were already
tabulated. No exhaustive search of its examples, sequels or cited papers was completed here.

### Borcea and Streinu (2021): the family itself

C. S. Borcea and I. Streinu, *Periodic tilings and auxetic deployments*,
Mathematics and Mechanics of Solids **26**(2), 199–216 (2021; online 2020),
[doi:10.1177/1081286520948116](https://doi.org/10.1177/1081286520948116).
[Publisher article](https://journals.sagepub.com/doi/abs/10.1177/1081286520948116);
[Smith institutional record](https://scholarworks.smith.edu/csc_facpubs/261/).

**Read:** Section 7, *A glimpse into higher dimensions*, in publisher-indexed article text, plus abstract and metadata.
Open PDF: <https://scholarworks.smith.edu/csc_facpubs/313> (Wayback copy of `viewcontent.cgi?article=1262&context=csc_facpubs`), read by a referee. Section 7 describes bridged cross-polytopes in arbitrary
dimension, two body orbits, 2d bridges, nine deformation freedoms in 3D, and folding at bridge/edge ratio 2.

**Assessment:** family-level novelty is withdrawn. It is the same infinite framework at a different bar length: regular cross-polytopes, 2d bars per body, any d. Their bars are twice the cross-polytope edge, that is dots at t = sqrt2 - 1; ours are at 2 - sqrt2. They describe it with a periodicity lattice of two body orbits; our per-cell counts are for the finer lattice 2 Z^n, with one body per cell, and hold for any 0 < t < 1. Per-cell counts depend on the periodicity lattice: for an index-2 lattice the generic-phase counts per cell double (n(n-1) zero modes, n(n-3) self-stresses), and the phase-0 count depends on the sublattice (see [below](#fixed-lattice-modes-versus-flexible-lattice-deformation); which sublattice Section 7 uses was not checked here). Section 7 gives only the 3D count (nine degrees of freedom, citing Borcea, Streinu and Tanigawa) and a 3D folding sketch to a folded state, not a packing. It has no n-dependent flex counts, self-stresses, elastic constants or Poisson ratio. The n-dimensional counts are likely routine with existing machinery (priority unresolved). The unit-spring results C11/C12 = (n-1)^2 and nu_100 = 1/(n^2-n-1) are the strongest remaining candidate: not found in Section 7 or in Goodwin's abstract.

### Li, Refson and Dove (2023): physical hexaborides

L. Li, K. Refson and M. T. Dove, *The contribution of phonons to the thermal expansion of some simple cubic
hexaboride structures: SmB6, CaB6, SrB6 and BaB6*, Physical Chemistry Chemical Physics **25**, 10749–10758 (2023),
[doi:10.1039/D3CP01306E](https://doi.org/10.1039/D3CP01306E),
[publisher abstract](https://pubs.rsc.org/en/content/articlelanding/2023/cp/d3cp01306e).

**Read:** abstract and metadata; full article not reviewed. It studies four compounds by density-functional lattice
dynamics, considers rotations of B6 octahedra, and reports very low flexibility of their connected network.

**Assessment:** this is a relevant physical comparison, not evidence for the ideal all-unit-spring tensor or arbitrary-n
mode counts. No flat-band count was checked in its full text. Its physical result neither validates nor refutes the
ideal spring model without matching interactions and geometry.

## Exact claims and elementary consequences

The following are deductions for the repository's stated model, not paraphrases of the articles and not claims of
priority. In dimension n >= 3 the joints are `p(i,s) = s c e_i`, `c = 2 - sqrt(2)`, repeated over `2 Z^n`.
Every cross-polytope edge and each bridge from `p(i,+)` to the next cell's `p(i,-)` is a bar. There is one body per
chosen cell, `2n` joints, `2n(n-1)` internal bars and `n` bridges. Counts also apply for `0 < c < 1` with this connectivity.

| Claim in PROOF.md | Hypotheses and elementary calculation | Priority assessment |
|---|---|---|
| Internal self-stresses: `n(n-3)/2` | The internal cross-polytope is infinitesimally rigid; rank is `2n^2 - n(n+1)/2`. Subtract rank from its `2n(n-1)` edges. | Elementary consequence of internal rigidity and rank-nullity; no novelty claim. |
| Bloch modes: `n(n-1)/2 + k` | Fixed lattice, complex phase `z_i = exp(2 i q_i)`; exactly k coordinates have `z_i = 1`. Derivation below. | Elementary rigid-unit application; exact-formula priority unresolved. |
| Bloch self-stresses: `n(n-3)/2 + k` | Rigidity matrix has `2n^2` columns and `2n^2-n` rows; mode minus stress dimension is n. | Follows from the mode count by rank-nullity. |
| `C11/C12 = (n-1)^2`, `nu_100 = 1/(n^2-n-1)` | Every internal and linking bar is a unit central-force spring; homogeneous strain, relaxed joints, cell volume `2^n`, no angular stiffness or prestress. | Exact spring-model calculation; priority unresolved. |

For the Bloch count, internal rigidity forces the cell velocity to have the form

    u(i,s) = v + s c Omega e_i,    Omega^T = -Omega.

A rotation gives transverse velocity at each bridge endpoint: `e_i . Omega e_i = 0`. The longitudinal bridge equation is
therefore `(z_i - 1) v_i = 0`. The `n(n-1)/2` independent rotations survive at every phase; each `z_i = 1` adds the
translation `e_i`. At phase zero, the count includes n translations; it is not a quotient by all Euclidean motions.
These are first-order modes, not a proof of finite bar-preserving motion or of sphere non-overlap along a path.
The n = 2 internal square is flexible and requires the separate treatment already in PROOF.md.

The elastic calculation requires finite stiffness of the internal edges as well as the bridges. After relaxation,
PROOF.md obtains, with `e_i = eps_ii`,

    E_min = 2(n-2)/n sum_i e_i^2 + 2/n^2 (sum_i e_i)^2,
    C11 = 4(n-1)^2/(n^2 2^n),    C12 = 4/(n^2 2^n),    C44 = 0.

These values follow from minimizing spring energy. Rigid-body kinematics alone does not supply finite `C11` or `C12`:
it eliminates the internal distortions that carry some of that energy. With C44 = 0 the tensor is only
positive semi-definite; `nu_100` describes normal strains under stress along a cube axis, not an isotropic material
Poisson ratio or a measured property of a real compound.

## Fixed-lattice modes versus flexible-lattice deformation

A Bloch nullity with fixed period vectors is a different quantity from a periodic configuration-space dimension with
variable lattice. The repository's [m = 1 fold](../fold-m1/PROOF.md) allows a flexible lattice and gives the three-dimensional
first-order count 6. Borcea and Streinu's stated 9 must not be compared with the generic fixed-lattice Bloch count 3.

Fixed-lattice counts per cell also depend on the periodicity lattice. A description with two body orbits uses an index-2
sublattice of `2 Z^n`, and every such sublattice is `L' = {2m : eps . m even}` for a non-zero `eps` in `{0,1}^n`. A Bloch
phase of `L'` lifts to two phases of `2 Z^n`, the second with `z_i` multiplied by `(-1)^(eps_i)`, and after a unitary
change of basis of the displacements and of the bar stretches the rigidity matrix of one `L'`-cell (two bodies, `4n^2`
coordinates, `4n^2 - 2n` bars) is the direct sum of the two `2 Z^n` matrices; this is the supercell counting by phase
nullities noted above for Power's Section 5(b). For n >= 3, part (b) of PROOF.md then gives per `L'`-cell `n(n-1)` zero
modes and `n(n-3)` self-stresses at a generic phase, and `n^2 + n - |eps|` zero modes and `n^2 - n - |eps|` self-stresses
at phase 0, where `|eps|` is the number of ones in `eps`. For the checkerboard sublattice `eps = (1, ..., 1)` (fcc in 3D)
this is `n^2` and `n(n-2)`. Which sublattice Section 7 uses was not checked against its text, so no phase-0 number is
attributed to it here. This is an internal audit deduction (Claude, Anthropic, via Claude Code), pending review. A
scratch computation from coordinates, outside the replay, agreed for n = 3, 4, 5, three sublattices and t = 2 - sqrt2,
sqrt2 - 1 and 0.3.

Here is conditional Maxwell bookkeeping for a body model with b body orbits and `nb` bridge orbits. Each body has
`n(n+1)/2` infinitesimal rigid-motion variables; a variable lattice adds `n^2`. After removing overall translations and
rotations, and **if all bridge constraints are independent**, the count is

    b n(n+1)/2 + n^2 - n(n+1)/2 - nb = (b+1)n(n-1)/2.

It gives 6 for b = 1 and 9 for b = 2 in dimension three. This accounts for the bookkeeping difference; it is not a rank
proof for either deployment configuration, nor a proof that an infinitesimal flex integrates to a finite path.
Exceptional geometries can have dependent constraints, and regularity must be established before identifying
first-order counts with local configuration-space dimension.

## Independent review packet

**Status: external human review pending; no review request has been sent.** The appropriate reviewer is a periodic
rigidity or lattice-elasticity specialist. The concrete packet is PROOF.md, this audit, the deterministic certificate,
and replay command `make verify-rigidity PY=python3` in the environment pinned by `requirements.txt`. The checker independently builds
matrices for n = 2..8; it does not prove general n or establish literature priority.

Focused questions:

1. Is the internal-rigidity/block proof complete for every n >= 3, including the explicit stress basis and the phase
   intersections where exactly k coordinates satisfy `z_i = 1`?
2. Section 7 treats the same infinite framework with a different t, described with two body orbits per period; our
   per-cell counts are for `2 Z^n`. Do Goodwin's full text, Borcea–Streinu–Tanigawa or other periodic body-and-bar results already state the n-dimensional counts or the unit-spring C11/C12 and nu_100?
3. Is the homogeneous-strain relaxation and normalization correct for finite unit spring stiffness on **all** bars?
   Are C44 = 0 and the axial Poisson formula stated with sufficient hypotheses?
4. Are these dimensional formulas already published, or routine consequences best presented as a worked example?
   Which exact citations should replace any remaining candidate-priority wording?

A separate continuous-path proof is still required for the folding results. Neither this literature correction nor
the rigidity replay discharges that obligation. A focused paper should claim only what survives these comparisons
and review; there is no guarantee that any candidate formula is new.

## Search and access log

2026-10-12, primary-only evidence retained:

- `"Periodic tilings and auxetic deployments" "connecting bars"`: publisher-indexed Section 7 text reached.
- `site:journals.sagepub.com/doi/abs/10.1177/1081286520948116 "A glimpse into higher dimensions"`: section heading and construction confirmed.
- `site:ora.ox.ac.uk "Rigid unit modes and intrinsic flexibility" "0cef4bf2"`: Goodwin abstract and metadata reached.
- `site:pubs.rsc.org/en/content/articlehtml/2023/cp/d3cp01306e "Abstract" "very low flexibility"`: publisher landing abstract reached.
- Power's DOI-matched PMC article opened directly; Sections 2, 4 and 5 inspected. Smith's Borcea–Streinu record opened for metadata.

This is a bounded comparison of four sources, not a systematic literature review. Indexed article text was used where
direct opening failed; full Goodwin and Li papers, source equation images, and Borcea–Streinu's cited higher-dimensional
rigidity literature remain to be reviewed. Borcea–Streinu Section 7 was subsequently read in the open PDF. Search absence is not evidence of novelty.
