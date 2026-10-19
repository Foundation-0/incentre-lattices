# Non-overlap of the continuous fcc motion: a symmetry-reduced written proof

Evidence: **PROVEN (written argument), internal review only.** This is a second proof of the non-overlap step
(Section 3) of [PROOF.md](PROOF.md). It does not replace PROOF.md, whose exact Bernstein certificate
(`make verify-fold-continuous`) remains the certified replay. The label of the result is unchanged.
Identity check: `make verify-fold-continuous-hand`, which runs `results/fold-continuous-fcc/hand_proof_check.py`
(sympy; replayed by `make verify`, but not hashed into `certificate.json`). It checks the 48-case centraliser of
Section 1 exactly, and symbolically the equalities of Sections 3 and 4: the formula for `G`, `d^2`, `alpha`, `gamma`,
`|nL|^2` for `n = 2e_i`, and every `w`, `|w|^2` and slack in the table. It also checks `35/3 > 6 + 4 sqrt2` and that the
nine classes cover the 26 vectors with `1 <= |n|^2 <= 3`. On a grid of 201 values of s it compares the reduced
conditions with the direct minimum over all motif differences, and checks that the actual `4 m1 - 2` and `2(m1 + m2)`
never exceed the largest candidate listed in the table. The block form of `G` in the rotated basis (with
`beta^2 = 2t^4`), the other inequalities of Section 3, the nonnegativity of the table entries and the contact counts
of Section 5 are checked by inspection only.

Notation as in PROOF.md: motif `P = {+-e_i}`, balls of radius `1/sqrt(2)`, row-vector lattice
`L(s) = d I + t C`, `d = 2 + a`, `a^2 + 2t^2 = 2`, `0 <= t <= 1`, `a >= 0`. So `d^2 = 6 + 4a - 2t^2`.

## 1. Symmetry along the motion

For `0 < s <= 1` the framework (balls and bars) has point group `{I, -I, m, -m}` of order 4, where `m` swaps the
second and third coordinates. Indeed `m` commutes with `L(s)` (it fixes `L_1`, swaps `L_2` and `L_3`), maps the motif
to itself, and maps bridge vectors to bridge vectors: `b_1 m = b_1`, `b_2 m = b_3`, where `b_i = L_i - 2e_i = a e_i + t C_i`.
Conversely, let `x -> x g + lambda` be a symmetry. Octahedron edges lie in triangles of bars and bridges in none (a
common neighbour of `e_i` and `L_i - e_i` would need `|L_i|^2 <= 4`, but `|L_i|^2 = 6 + 4a`), so it maps octahedra to
octahedra and bridges to bridges; hence `lambda` is in `Z^3 L(s)` and `g` is a signed permutation,
`e_i g = eps_i e_pi(i)`. The bridge at a ball `eps e_j + nL` has vector `eps b_j`, so `b_i g = eps_i b_pi(i)`, i.e.
`B g = g B` with `B = a I + t C`; for `t > 0` this gives `C g = g C`, and of the 48 signed permutations only `+-I` and
`+-m` commute with `C` (checked exactly in `hand_proof_check.py`). The space group is `C2/m`
(No. 12); balls occupy two orbits, `{+-e_1}` on the mirror and `{+-e_2, +-e_3}` in general position. At `s = 0`
the group is `Pm-3m`. At `s = 1` the sphere packing is fcc, while the framework keeps the order-4 symmetry.

## 2. Reduction to octahedron pairs

Balls in one octahedron are at least `sqrt(2)` apart. For a nonzero integer row `n`, put `w = nL`, and let
`m1 >= m2 >= m3` be the sorted absolute values of the coordinates of `w`. The balls of the octahedra centred at `0`
and `w` are pairs `(p, q + w)`; their squared distances are `|w + v|^2` with `v = q - p` in `P - P`. Minimising the
linear term `2 w.v` over `v = +-2e_i`, `v = +-e_i +- e_j` and `v = 0` shows that no two of these balls overlap iff

    (A) |w|^2 - 4 m1 + 2 >= 0,     (B) |w|^2 - 2(m1 + m2) >= 0,     (C) |w|^2 - 2 >= 0.

Both conditions are unchanged under `n -> -n` (since `P - P = -(P - P)`) and under `n -> n m` (which swaps two
coordinates of `w`). So it suffices to treat `n` up to sign and up to the mirror.

## 3. Distant octahedra

`|nL|^2 = n G n^T` with `G = L L^T = (d^2 + t^2) I + t^2 J + d t E`, where `J` is the all-ones matrix and `E` has
entries `E_23 = E_32 = -2`, zero elsewhere. In the orthonormal basis `e_1`, `(e_2 + e_3)/sqrt2`, `(e_2 - e_3)/sqrt2`,
`G` has the eigenvalue `(d + t)^2 >= 4` on the last vector and the block `[[alpha, beta], [beta, gamma]]` on the
first two, with

    alpha = d^2 + 2t^2 = 6 + 4a,   gamma = d^2 + 3t^2 - 2dt = 6 + 4a + t^2 - 4t - 2at,   beta^2 = 2t^4.

Then `alpha - 7/3 >= 11/3`, and `gamma - 7/3 >= 11/3 - 4t + t^2 >= 2/3` (use `a(4 - 2t) >= 0`; the quadratic decreases on
`[0,1]`). Hence `(alpha - 7/3)(gamma - 7/3) >= 22/9 > 2 >= beta^2`, so the smallest eigenvalue of `G` is at least `7/3` for
every s. If `|n|^2 >= 5` then `|nL|^2 >= 35/3 > 6 + 4 sqrt2 = (2 + sqrt2)^2` (since `17/3 > 4 sqrt2`, i.e. `289 > 288`), and every
ball pair is more than `|nL| - 2 > sqrt2` apart. If `|n|^2 = 4` then `n = +-2e_i` and `|nL|^2 = 4(6 + 4a) >= 24`.

## 4. The nine near classes

The 26 vectors with `1 <= |n|^2 <= 3` fall into nine classes. For each, `w` and `|w|^2` follow from `L = dI + tC`.
(C) is clear in every case (`|w|^2 >= 6`). The table gives the slack in (A) and (B); where the largest coordinate is not
determined, the slack is given for every candidate, so the condition holds whichever coordinate is largest.

| n (with -n, mirror image) | w | `|w|^2` | slack in (A) | slack in (B) |
|---|---|---|---|---|
| e_1 | (d, -t, -t) | 6 + 4a | 0 (bridge) | 2(1-t) + 2a |
| e_2, e_3 | (t, d, -t) | 6 + 4a | 0 (bridge) | 2(1-t) + 2a |
| (1,1,0), (1,0,1) | (d+t, d-t, -2t) | 12 + 8a + 2t^2 | 2(1-t)^2 + 4 + 4a | >= 2(1-t)^2 + 2 + 4a (using m1+m2 <= 2d+2t) |
| (1,-1,0), (1,0,-1) | (d-t, -(d+t), 0) | 12 + 8a - 2t^2 | 2(1-t)(3+t) + 4a | 4 + 4a - 2t^2 |
| (0,1,1) | (2t, d-t, d-t) | 12 + 8a + 2t^2 - 8t - 4at | m1 = 2t: 2(1-t)(7-t) + a(8-4t); m1 = d-t: 2(1-t)^2 + 4 + 4a(1-t) | m1+m2 = d+t: 2(1-t)(4-t) + a(6-4t); = 2d-2t: 2(1-t)^2 + 2 + 4a(1-t) |
| (0,1,-1) | (0, d+t, -(d+t)) | 2(d+t)^2 | 2(d+t-1)^2 | 2(d+t)(d+t-2) |
| (1,1,1) | (d+2t, d-2t, d-2t) | 18 + 12a + 6t^2 - 8t - 4at | 6(t-4/3)^2 + 4/3 + a(8-4t) | 6(t-2/3)^2 + 22/3 + 4a(2-t) |
| (1,-1,-1) | (d-2t, -d, -d) | 18 + 12a - 2t^2 - 8t - 4at | 2(1-t)(5+t) + 2 + a(8-4t) | 2(1-t)(5+t) + a(8-4t) |
| (1,1,-1), (1,-1,1) | (d, d, -(d+2t)) | 18 + 12a - 2t^2 + 8t + 4at | 12 + 8a - 2t^2 + 4at | 10 + 8a - 2t^2 + 4t + 4at |

Every entry is nonnegative for `0 <= t <= 1`, `a >= 0` (and `d + t > 2` because `a` and `t` do not vanish together).
With Sections 2 and 3 this proves non-overlap for every pair of balls and every `s` in `[0,1]`.

## 5. Contacts

The slack vanishes only (i) identically, for the three bridge orbits (class `e_i`, condition (A)), and (ii) at `t = 1`, `a = 0`
(that is `s = 1`), in classes `e_i` (B), `(1,-1,0)` (A), `(0,1,1)` (A with `m1 = 2t`, and B with `m1 + m2 = d + t`) and
`(1,-1,-1)` (B). So for `0 <= s < 1` the only touching pairs are the 15 bars per cell (5 per ball), and at `s = 1`
21 further contacts per cell appear at once, giving 12 per ball. These are exactly the 24 vanishing distance
polynomials, beyond the 18 identically zero ones, that give the 114 zero Bernstein coefficients of PROOF.md
(18 x 5 + 24 = 114).

## Credit

Written by Claude Opus 5.5 (Claude Code) for the v2 preprint, from the curve and proof of PROOF.md (GPT-6 Codex,
checked by GPT-6.1 Sol). Internal work only; no person has checked it yet.
The converse step of Section 1 (only `+-I` and `+-m`), which replaces an unrecorded numerical check at several s,
and the centraliser and candidate checks in `hand_proof_check.py` were added in an internal audit by Claude
(Anthropic), via Claude Code. They are pending review; no person has checked them either.
