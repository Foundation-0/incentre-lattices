# Formalization: an open track

Nothing in this repository is formally verified. This page lists which statements are good targets for a proof assistant (Lean 4), what is already
in `formal/`, and what would be needed. The reference repository Swapnil-jain/gf2m-circuits ships Lean 4 kernel proofs for its identities; that is the
bar to aim at.

## What exists

`formal/Arithmetic.lean` (core Lean 4, no Mathlib; checked locally with Lean 4.34.1, not run in CI, `make lean`) proves, for n = 3..64 by `decide`, the arithmetic
bookkeeping around the rigidity result: (zero modes) - (self-stresses) = n per cell, (n-1)^2 + (n-2) = n^2 - n - 1, and the bar count 2n(n-1) + n = 2n^2 - n, plus the six values of
the denominator of nu_100. **This is counting arithmetic only. It says nothing about rigidity, ranks or elasticity.** Its value is that the numbers quoted in the
write-ups are the ones the formulas give.

## Good targets, easiest first

1. **Exact algebraic identities in Q(sqrt2).** (4 pi/3)(1/8 + 6 r_S^3) = (pi/6)(30 sqrt2 - 41) with r_S = (sqrt2 - 1)/2; the fcc-window endpoints ((sqrt2-1)/6)^(1/3) and
   (4 + 7 sqrt2)^(-1/3) (the key facts are gamma^3 = (sqrt2-1)/6 and (5 sqrt2 - 7)(10 + 7 sqrt2) = sqrt2); h* = 1 - sqrt2/2 and x* = (sqrt2-1)/2; c sqrt2 = 2R and c + R = 1 for c = 2 - sqrt2;
   R = sqrt2 - 1 = tan(pi/8). With Mathlib these are `ring_nf`, `nlinarith` and `Real.tan_pi_div_eight` style goals. Also the 2D identities L(pi/4) = 4R and gap(theta) = 2 sqrt2 c cos(theta + pi/4)
   (trigonometric addition formula), and a(60 degrees) = 8c/3 with the density 27 pi sqrt2/256.
2. **The cab6 upper bound.** For 0 < x < 1/2: min((1/2 - x)/sqrt2, x) <= (sqrt2 - 1)/2 with equality only at x = (sqrt2-1)/2, and the two-branch bound on phi from rL <= 1/2 and rS <= m*. Real analysis, a few lines.
   Attainment and all-pairs non-overlap are PROVEN (written proof, internal review only) in `results/cab6-packing/ATTAINMENT.md`: the squared distance to Z^3 splits coordinatewise. Its distance lemma
   (d_CaCa = 1, d_CaB = sqrt((1/2 - h)^2 + 1/2), d_BB = min(sqrt2 h, 1 - 2h) for 0 < h < 1/2) and the attainment step are a short formalisation target.
3. **Rank counts for fixed n by `decide`.** For each fixed n (3 to 8), the rank of the rigidity matrix with entries in Z (starting from the unit bar directions of `results/rigidity-nd/rigidity_nd.py`, multiply each
   intra row by sqrt2 so its entries are 0, +-1, and leave the cross rows, which are already +-e_i; row scaling does not change the rank)
   at q = 0, and at a chosen rational point for the phases (a lower bound for the generic rank; the phases e^(i q.L) can be taken as Gaussian rationals on the unit circle, e.g. (3+4i)/5, so the
   entries lie in Q(i)), with an explicit kernel basis (the rotations) for the upper bound. Exact arithmetic in Q or Q(i), no floating point,
   and it removes the SVD tolerance from the story. Generic-q counts need the phases as indeterminates (rank over the function field) or a Schwartz-Zippel style argument; a point-evaluation plus the
   explicit kernel is enough for equality.
4. **The block identities and the minimisation.** Lemma 1 (the stretch formula) is a polynomial identity in the entries, closed by `ring`; the sum identity (*) likewise; the minimisation
   E_min = (2(n-1)^2/n^2) sum e_i^2 + (2/n^2) sum_{i != j} e_i e_j is a positive-definite quadratic with explicit minimiser y_k = (1/n)(2 e_k - (1/n) sum e_j). For symbolic n these need sums over
   `Finset.range n`; fixed n up to 8 is routine. The all-n statement is a good formalisation exercise (kernel of the unsigned incidence matrix of K_n for n >= 3).
5. **The exact fcc motion** (`results/fold-continuous-fcc/PROOF.md`, PROVEN by written argument and exact finite certificate, internal review only). The bar identity a^2 + 2t^2 = 2; the
   symmetric-part bound |w L| >= |w| (Cauchy-Schwarz with d - t >= 1), which leaves the 171 integer shifts with |n|^2 <= 11; the 3248 quartics whose 16240 Bernstein coefficients in Q(sqrt2) are all
   nonnegative (finite exact arithmetic on rational pairs u + v sqrt2, no subdivision); the parity/index endpoint argument (even-sum integer rows of determinant 12, index 6 in D3, six odd-sum motif
   points in distinct cosets); and D'(t) < 0 for 0 < t < 1. The finite part is at present checked by two Python implementations (`curve.py` and `independent.py` there), not by a proof assistant.
6. **The radius-ratio theorem (H1 and H2).** Elementary Euclidean geometry (reflections, the bisector computation sin(theta) = sin(phi/2)); heavier in Mathlib because it needs polytopes and facets.
   A plausible first step is the planar and three-dimensional cases with explicit coordinates.

## Not good targets (for now)

* The sampled fold paths (793 or 1193 floating-point states at m = 2, 333 or 268 at m = 1). Each is a sampled finite check; a formal version would be a certified (interval) path, not a
  formalisation of what exists. For fcc the exact motion in item 5 is the target instead; a certified hcp path is the open problem listed in `docs/status.md`.
* The prior-art statements (they are literature searches).

## How to contribute

Open a pull request with a Lean file under `formal/`, say which statement it proves (exactly, including hypotheses), the Lean and Mathlib versions, and how to check it. A Lean file that proves a
different statement than the one in `PROOF.md` is welcome but must say so. See `.github/pull_request_template.md`.
