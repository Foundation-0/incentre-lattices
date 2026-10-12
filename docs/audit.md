# Audit: what is checked, by what, and how far it goes

Evidence labels used everywhere in this repository:

* **PROVEN**: a written proof exists (in the result's `PROOF.md`). A proof is only as good as its review; so far only internal review has happened.
* **VERIFIED**: a finite computation that replays. It supports the statement on the cases run, with the tolerances stated. It is not a proof.
* **KNOWN**: published; cited in the `PROOF.md`.
* **NEW-COMPUTATION**: our own computation, not independently reproduced beyond the checkers listed.

"Finite check, not formal proof" applies to every VERIFIED entry. Nothing here has been formalised in a proof assistant (see `docs/formalization.md`).

## Per result

### results/rigidity-nd: counts and nu_100 for n = 2..8

| Claim | Label | Checked by |
|---|---|---|
| Generic-q zero modes n(n-1)/2 and self-stresses n(n-3)/2; cross-polytope rigid; q = 0 and plane counts; C11, C12, zero shear, nu_100 = 1/(n^2-n-1) | PROVEN for all n >= 3 (`PROOF.md`) | the written proof; replay only supports it |
| The same for n = 2..8 | VERIFIED | `rigidity_nd.py`: SVD ranks of the Bloch rigidity matrix at fixed-seed random q, tolerance 1e-9; elastic tensor by full relaxation (Schur complement), 1e-12 |
| Steps 1 to 3 of the proof | VERIFIED symbolically for n = 2..8 | `proof_steps.py` (sympy) |
| Controls | must be rejected | n >= 3 counts at n = 2; wrong Poisson formula; wrong C11/C12; deleting a cross bar; wrong stress-basis dimension |

Known limits: numerical ranks depend on an SVD threshold; the proof has had internal review only; n = 2 is outside the theorem and is reported separately;
the model is ideal central-force springs, so shear moduli of exactly zero say nothing about real crystals.

### results/radius-ratio: hole/ball equals the classical radius ratio

| Claim | Label | Checked by |
|---|---|---|
| Under H1 (mirror walls) or H2 (insphere plus equal dihedral angles) all dots equal, adjacent balls touch, hole/ball = classical ratio, and 1/sin(phi/2) - 1 under H2 | PROVEN (`PROOF.md`) | the written proof |
| 20 polytopes (Platonic, Catalan, bipyramids, 4D regular, snub-cube dual, uniform trapezohedron) | VERIFIED | `radius_ratio.py`: inscribed balls by `linprog`, tolerances 1e-6 |
| The corrected hypotheses: facet-transitivity plus an insphere is not enough | VERIFIED (counterexamples) | trapezohedra at 0.7x and 1.3x: contacts open to 1.017 and 1.006 times 2r; prism and cuboctahedron outside the hypotheses |
| CN7, CN9 and dodecahedral CN8 shapes are not inscribed | VERIFIED for J13, J51, J84; classification of deltahedra is cited, not replayed | `radius_ratio.py` |

Known limits: floating-point LP at 1e-6; the families satisfying H1 or H2 are not classified; the capped CN7 shapes under H2 were not examined.

### results/cab6-packing: closed form and the fcc window

| Claim | Label | Checked by |
|---|---|---|
| phi_max(gamma) as a two-branch closed form; maximum (pi/6)(30 sqrt2 - 41) at gamma = sqrt2 - 1; Ca held by Ca-Ca, Ca-B not touching; window 0.41023 < gamma < 0.41591 | PROVEN upper bound (`PROOF.md`); exact values VERIFIED symbolically | sympy |
| Attainment and no other overlap | VERIFIED (finite) | `cab6_packing.py`: all pairs in a 5 x 5 x 5 block, 15 size ratios, agreement with the closed form to 1e-9 |
| Controls | must be rejected | wrong constant, wrong branch, gamma = 0.40 and 0.42 beating fcc, "Ca touches B", wrong Wyckoff reading of x |

Known limits: the all-pairs non-overlap at the optimum is checked, not argued; the model has one free coordinate and the cubic lattice only.

### results/fold-2d: hinged squares as a ball framework

| Claim | Label | Checked by |
|---|---|---|
| The linkage and the fold | KNOWN (Flores; Stuart) | citation only |
| gap(theta) = 2 sqrt2 c cos(theta + pi/4); first contact at pi/4; L = 4R; strain -(sqrt2-1)^2; density 0.5390 to pi/4 | PROVEN exactly | sympy |
| No closer non-bonded pair before pi/4; the incentre is the only c with a single-distance square lattice | VERIFIED on a scan (six values of c, 0.01-degree grid) | `fold2d.py` |

Known limits: the symmetric path only (the framework has other flexes, including shear); the uniqueness of the incentre is a scan, not a proof.

### results/fold-3d-cubic: the cubic-restricted [111] fold

VERIFIED only. Exact algebra at 60 degrees (a = 8c/3, six contacts, density 27 pi sqrt2/256 = 0.4686) by sympy; the first non-bonded contact on a 3000-point grid;
the coordination sequence (ten terms, equal to RCSR `sms`, not an isomorphism); a grid scan of trigonal lattices (density up to 0.6477) that shows the cubic lock is not
a jam. Grid adjacency is not a continuous path.

### results/fold-close-packing: the path to hcp and fcc (m = 2 only)

| Claim | Label | Checked by |
|---|---|---|
| A 793-state sequence from the open cab state to hcp (density pi/sqrt18, 12 contacts, CNA 288 (4,2,1) + 288 (4,2,2)) keeping every bar and with no overlap | VERIFIED (finite) | locked `autoresearch-fold/verify.py` and independent `chk_independent.py` (a parametrised copy of `autoresearch-fold/independent/chk.py`) |
| A 1193-state extension by two layer slides to pure fcc (576 (4,2,1)) | NEW-COMPUTATION, VERIFIED by the same two checkers | `make_fcc_path.py` generates it; both checkers accept it |
| The checker itself | 16 unit tests including must-fail cases; hash lock | `autoresearch-fold/test_verify.py`, `verify.sha256` |
| Controls | must be rejected | displaced dot; lattice scaled in one state; a half-way blend of two states |

Known limits (the important ones):

* **Continuity is certified only by sampling.** Consecutive states differ by at most 0.01 in dot position and lattice; overlap is checked at the states, at midpoints and at 10
  interpolated sub-states. Linear interpolation between states shortens or lengthens bars by up to 1.5e-5 (relative), so the interpolated curve is not exactly bar-preserving.
  This is not an interval-arithmetic or analytic certificate of a continuous motion.
* **Only the m = 2 family** (48 balls) is known to reach close packing. m = 1, 3 and 4 are unexplored or failed.
* Both checkers were written inside the project, by the same group; they share the model definition. Neither is an outside review.
* The density along the path is not monotone (dips to 0.7294 during the fcc slides).
* First-order rigidity at the end state is not a proof of jamming.
* The search is deterministic; "several seeds" gave one path.
* A geometry and rigidity result, not a prediction for real CaB6 or Prussian-blue materials.

## What is not covered at all

The interactive site (`*.html`), the research notes under `research/` (some of which print numbers by scripts not replayed here), the perception question
about triangle centres (an open experiment, not a result), and the UB12 and perovskite-tolerance computations in `research/chemistry`.
