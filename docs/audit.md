# Audit: what is checked, by what, and how far it goes

Evidence labels used everywhere in this repository:

* **PROVEN**: a written proof exists (in the result's `PROOF.md`). A proof is only as good as its review; so far only internal review has happened.
* **VERIFIED**: a finite computation that replays. It supports the statement on the cases run, with the tolerances stated. It is not a proof.
* **KNOWN**: published; cited in the `PROOF.md`.
* **NEW-COMPUTATION**: our own computation, not independently reproduced beyond the checkers listed.

"Finite check, not formal proof" applies to every VERIFIED entry. Nothing here has been formalised in a proof assistant (see `docs/formalization.md`).

## Per result

### results/fold-continuous-fcc: explicit continuous motion

**PROVEN (written argument and exact finite certificate), internal AI review only.** The motif is fixed and the
lattice follows the rational curve in `PROOF.md`. Bar preservation is an exact identity. A bound on the symmetric
part of the lattice excludes all integer shifts of squared norm at least 12. Every remaining squared-distance
inequality is a quartic with nonnegative Bernstein coefficients on the entire interval: 3248 polynomials and
16240 coefficients in Q(sqrt(2)), without floating arithmetic. `curve.py` constructs the exact certificate;
`independent.py` separately forms the polynomials by convolution and reverses the Bernstein expansion.
The fcc endpoint follows from an integer parity/index argument. A written derivative bound proves strictly
increasing density. These analytic reductions still require mathematical review; the code is not a proof assistant.

`make verify-fold-continuous` replays both arithmetic implementations, the endpoint, and must-fail controls.
All proof and replay sources are bound by the result certificate. A separate 501-state sample was accepted by both
existing, unchanged numerical fold checkers: zero overlap violations, 12 endpoint contacts per ball, fcc CNA
72 (4,2,1), and endpoint density pi/sqrt18. The sampled straight-line interpolants are not the exact analytic motion.
The older sampled paths and hcp result retain their previous evidence labels. Priority and outside human review
remain open; no physical transition or material stability is claimed.

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
| H1 (every wall is a mirror) on the 18 cases claimed under H1 | VERIFIED | `verify.py` reflects the vertex set in every wall and requires it to map onto itself (1e-9); the H2-only cases and the counterexamples are asserted not to have mirror walls |
| The corrected hypotheses: facet-transitivity plus an insphere is not enough | VERIFIED (counterexamples) | trapezohedra at 0.7x and 1.3x: contacts open to 1.017 and 1.006 times 2r; prism and cuboctahedron outside the hypotheses; the premises (an inscribed ball in every pyramid) and the quoted numbers (prism inradii, contact ranges, deltahedra orbit radii) are asserted |
| CN7, CN9 and dodecahedral CN8 shapes are not inscribed | VERIFIED for J13, J51, J84; classification of deltahedra is cited, not replayed | `radius_ratio.py` |

Known limits: floating-point LP at 1e-6; the families satisfying H1 or H2 are not classified; the capped CN7 shapes under H2 were not examined.

### results/cab6-packing: closed form and the fcc window

| Claim | Label | Checked by |
|---|---|---|
| phi_max(gamma) as a two-branch closed form; maximum (pi/6)(30 sqrt2 - 41) at gamma = sqrt2 - 1; Ca held by Ca-Ca, Ca-B not touching; window 0.41023 < gamma < 0.41591 | PROVEN upper bound (`PROOF.md`); exact values VERIFIED symbolically | sympy |
| Attainment and no other overlap | PROVEN (`ATTAINMENT.md`, exact all-pairs classification) | corroborated by `cab6_packing.py`: all pairs in a 5 x 5 x 5 block, 15 size ratios, agreement with the closed form to 1e-9 |
| Controls | must be rejected | wrong constant, wrong branch, gamma = 0.40 and 0.42 beating fcc, "Ca touches B", wrong Wyckoff reading of x |

Known limits: the model has one free coordinate and the cubic lattice only.

### results/fold-2d: hinged squares as a ball framework

| Claim | Label | Checked by |
|---|---|---|
| The linkage and the fold | KNOWN (Flores; Stuart) | citation only |
| gap(theta) = 2 sqrt2 c cos(theta + pi/4), which vanishes at pi/4; L = 4R; strain -(sqrt2-1)^2; density 0.5390 to pi/4 | PROVEN exactly | sympy |
| That pi/4 is the first new contact (no closer non-bonded pair before pi/4); among the six c scanned, the incentre is the only one with a single-distance square lattice; for c < 2 - sqrt2 the fold locks elsewhere (66.37, 56.96, 49.51 degrees) | VERIFIED on a scan (six values of c, 0.01-degree grid), not proved | `fold2d.py` |

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
| The second checker and the fcc path generator are part of the lock | sha256 pins, checked by `make verify-fold-cp` (so by `make verify` and CI) | `results/fold-close-packing/checkers.sha256` pins `chk_independent.py` and `make_fcc_path.py`; the second checker derives its image-shift block from the lattice of each state |
| Controls | must be rejected | displaced dot; lattice scaled in one state; a half-way blend of two states |

Known limits (the important ones):

* **Continuity is certified only by sampling.** Consecutive states differ by at most 0.01 in dot position and lattice; overlap is checked at the states (tolerance 1e-9), at step
  midpoints (locked checker, all pairs) and at 9 interior sub-states u = 0.1, ..., 0.9 of each step (independent checker, non-bonded pairs), with tolerance 1e-4 at the midpoints and
  sub-states. Linear interpolation between states shortens or lengthens bars by up to 1.5e-5 (relative), so the interpolated curve is not exactly bar-preserving, and at the sampled
  sub-states it brings two non-bonded balls to 0.999994 * 2R, accepted under the 1e-4 tolerance. Given the state checks and the dot step bound 0.01, every pair with zero lattice
  shift stays at least 0.999927 * 2R apart on the linear interpolant, so the midpoint and sub-state samples can fail only for pairs across a nonzero lattice shift (internal audit
  work, not yet reviewed; the argument is in `results/fold-close-packing/PROOF.md`, Known limits).
  This is not an interval-arithmetic or analytic certificate of a continuous motion.
* **This result is m = 2 only** (48 balls). m = 1, and through it every m in principle (the lift needs subdivided paths; see `results/fold-m1/PROOF.md`), is in `results/fold-m1`.
* Both checkers were written inside the project, by the same group; they share the model definition. Neither is an outside review.
* The density along the path is not monotone (smallest step -7.2e-5 on the hcp path; in each fcc slide it falls by up to 2.3e-4 per step to 0.7294 and rises back).
* First-order rigidity at the end state is not a proof of jamming.
* The stored `autoresearch-fold/best/path.npz` is the certificate; the search that found it is not a certificate and is not guaranteed to reproduce that exact path on another machine. The original run bounded stage B by wall clock, and one slower machine stopped at a jammed state (density 0.618718). `search.py` is now bounded by a step count with BLAS pinned to one thread; a rerun on one machine gave a different valid path (791 states, diverging from the stored one at state 285) reaching the same end density 0.740480.
* A geometry and rigidity result, not a prediction for real CaB6 or Prussian-blue materials.

### results/fold-end-exact: exact closed forms of the close-packed end states (m = 2)

| Claim | Label | Checked by |
|---|---|---|
| The hcp and fcc closed forms in `PROOF.md` have every bar exactly 2R, no pair closer than 2R, 12 contacts per ball, density exactly pi/sqrt18, and CNA hcp / fcc | PROVEN (exact integer arithmetic) | `exact_end.py` (own bar graph; does not import either fold checker) |
| The final states of `autoresearch-fold/best/path.npz` and `results/fold-close-packing/path_fcc.npz` are these closed forms | VERIFIED (1e-7, after a rigid motion) | `exact_end.exactify` |
| Controls | must be rejected | one dot moved by K/3; axis permutation dropped; wrong lattice row; hcp labelled fcc |

Known limits: the end states only, not the path; the image-shift box for the overlap check is computed in floating point with one spare lattice step; one implementation, internal review only.

### results/fold-m1: close packing in the 1 x 1 x 1 cell

| Claim | Label | Checked by |
|---|---|---|
| A 333-state path from the open cab state to hcp at m = 1 (free density descent, no target) | VERIFIED (finite, sampled), NEW-COMPUTATION | locked `autoresearch-fold/verify.py` and `results/fold-close-packing/chk_independent.py` (9 interior sub-states per step) |
| A 268-state path to pure fcc at m = 1 (steered to an explicit end state), det L strictly decreasing at every step | VERIFIED (finite, sampled), NEW-COMPUTATION | the same two checkers; monotonicity at the samples only |
| Both end states exact (12 contacts, density pi/sqrt18, hcp / fcc) | PROVEN (integer arithmetic) | `results/fold-end-exact/exact_end.py` |
| Controls | must be rejected | displaced dot; scaled lattice; 100 deleted states; monotone test on the hcp path |

Known limits: sampled continuity as for m = 2; the lift to larger m is an argument, not a replayed check, and needs subdivided paths (the locked checker rejects the unsubdivided lifts, only on the continuity step limit); searches reproduce byte for byte on one machine only.

## What is not covered at all

The interactive site (`*.html`), the research notes under `research/` (some of which print numbers by scripts not replayed here), the perception question
about triangle centres (an open experiment, not a result), and the UB12 and perovskite-tolerance computations in `research/chemistry`.
