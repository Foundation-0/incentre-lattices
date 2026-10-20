# The close-packed end states of the fold are exact (m = 2)

Evidence label: **PROVEN (exact integer arithmetic) for the two closed forms below; VERIFIED (finite check, 1e-7) that the final
states of the stored paths are these closed forms.** This says nothing new about the paths themselves, which stay VERIFIED by
sampling (`results/fold-close-packing/PROOF.md`).
Replay: `make verify-fold-exact` (about 30 seconds; numpy only). Certificate: `certificate.json`.

## Statement

Model as in `results/fold-close-packing/PROOF.md`: dots at cell centre ± K e_i with K = 2 − √2, balls of radius R = √2 − 1,
bars of length 2R (12 octahedron edges per cell, and the +e_i dot of each cell to the −e_i dot of the next), 2 × 2 × 2 supercell
(48 balls, 120 bars). Since 2R = √2 K, a bar has squared length 2K².

**Closed form.** Work in units of K, in the cube frame of octahedron (0,0,0) with the origin at its centre. For a cell
(ix, iy, iz) with 0 ≤ ix, iy, iz ≤ 1:

    centre(ix, iy, iz) = ix t_x + iy t_y + iz t_z,      t_x = (1, −1, 2),  t_y = (1, 2, −1)
    dot (axis a, sign s) = centre + s e_π(a),           π = identity if ix + iy + iz is even,
                                                         π = (x, y, z) → (z, x, y) if it is odd

so every octahedron keeps its orientation up to a cyclic relabelling of its axes, in a checkerboard. The two end states differ
only in t_z and the third lattice row:

| End state | t_z | lattice rows (supercell) | det L |
|---|---|---|---|
| hcp (end of `autoresearch-fold/best/path.npz`) | (−4/3, 4/3, 4/3) | 2t_x, 2t_y, 2t_z = (−8/3, 8/3, 8/3) | 96 K³ |
| fcc (end of `results/fold-close-packing/path_fcc.npz`) | (−1, 1, 2) | 2t_x, 2t_y, (−2, 3, 3) | 96 K³ |

**Claim (PROVEN).** In both states:

1. all 120 bars have length exactly 2R;
2. no two balls, over all periodic images, are closer than 2R;
3. every ball touches exactly 12 others;
4. the density is exactly π/√18 (det L = 2 N K³ with N = 48 balls is equivalent to N (4/3)π R³ / det L = π/√18);
5. the common-neighbour signature of the contact graph is 288 × (4,2,1) + 288 × (4,2,2) for the first state (hcp) and
   576 × (4,2,1) for the second (fcc).

**Proof.** Multiply every coordinate by den = 3 (hcp) or den = 1 (fcc). All dot coordinates and lattice rows become integers in
units K/den, every squared distance becomes an integer times (K/den)², and (2R)² = 2 den² in those units. Items 1 to 5 are then
comparisons of integers, carried out by `exact_end.check_exact`. Items 2, 3 and 5 use every pair whose image shift lies in a box
|n_a| ≤ ⌈2.02R · |L⁻¹ column a|⌉ + 2 around the wrapped cell. That box is computed in floating point, so here is why it is
complete. The wrap subtracts W L from the dots Y, with W the floor of the floating-point product Y L⁻¹, so the wrapped dots Yw are
integer vectors whose fractional coordinates Yw L⁻¹ lie in [−ε, 1 + ε), where ε bounds the rounding error of that product (in
exact arithmetic the interval would be [0, 1)). A pair r = Yw_j + n L − Yw_i has n = r L⁻¹ − (Yw_j − Yw_i) L⁻¹. Each
coordinate of (Yw_j − Yw_i) L⁻¹ is less than 1 + 2ε in absolute value, and |(r L⁻¹)_a| ≤ |r| · |L⁻¹ column a| by
Cauchy–Schwarz, so every pair with |r| ≤ 2R has |n_a| < 2R · |L⁻¹ column a| + 1 + 2ε. The box adds 1% to the radius and one
further full step, so it contains every such n as long as ε ≤ 1 and the computed value of 2.02R · |L⁻¹ column a| is not below
the exact 2R · |L⁻¹ column a|. For the three states in this note ε is below 3e-16 and the computed column norms are within
2e-16 (relative) of the exact ones. ∎ (Internal audit, not yet reviewed: Claude (Anthropic), via Claude Code, wrote out this
bound and measured the two rounding errors against exact rational arithmetic with sympy.)

**Finding (VERIFIED, 1e-7).** After a rigid motion that puts octahedron (0,0,0) in its cube frame, the final state of each
stored path agrees with its closed form to better than 1e-7 in every coordinate (`exact_end.exactify`). So the stored paths end
exactly at these states, up to the round-off of the search.

**A second fcc form.** Replacing the third lattice row of the fcc state by 2t_z = (−2, 2, 4) gives another exact fcc state that
satisfies items 1 to 5 (also checked in the replay). The two differ only in which image the wrap-around z bars attach to. So the
m = 2 fcc end state is not unique, and the simplest fcc closed form is the fully linear one, L = 2(t_x, t_y, t_z) with
t_z = (−1, 1, 2).

## How the octahedra sit

Not part of the replay. A separate one-off analysis of the integer hcp state agrees with the description in
`results/fold-close-packing/PROOF.md`. In the hcp state the balls lie in four close-packed layers per c. All eight octahedra are parallel and fill one third of the
octahedral holes in alternate gaps. Of the 24 linking bars, 16 (all x and y bars) lie inside a layer and the 8 z bars cross the
two empty gaps.

## Controls (must be rejected)

* the hcp closed form with one dot moved by K/3 (bars broken);
* the hcp closed form without the cyclic axis permutation on odd cells (bars broken);
* the fcc closed form with third lattice row (−2, 3, 4) (bars and density broken);
* the hcp state labelled fcc (the common-neighbour test must say hcp).

## Known limits

* The exactness is of the end states only. The path to them is still a sampled path (steps ≤ 0.01 with midpoint and
  sub-state checks); a certified path is still open.
* The space group, Wyckoff sites and an orbit description of the octahedral-hole pattern are not written out here.
* Checked by one implementation (`exact_end.py`, written for this result; it rebuilds the bar graph and does not import either
  existing checker). It has not been reviewed outside the project.

## Credit

The end states were found by the fold search in `autoresearch-fold/` and the layer slides in `results/fold-close-packing/`.
This folder adds only the exact form. AI assistance: Claude Opus 5.5 (Claude Code) wrote `exact_end.py`, `verify.py` and this
note; no person has checked them yet.
