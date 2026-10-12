# Status: the current best per track

As of 2026-10-12. Labels as in `docs/audit.md`. No paper has been published; a preprint is pending.

| Track | Current best | Label | Replay |
|---|---|---|---|
| n-dimensional rigidity of the contact framework | n(n-1)/2 zero modes and n(n-3)/2 self-stresses per cell at generic q; nu_100 = 1/(n^2-n-1) (1/5, 1/11, 1/19, 1/29, 1/41, 1/55 for n = 3..8); all shear moduli 0 | PROVEN n >= 3; VERIFIED n = 2..8 | `make verify-rigidity` |
| Radius ratio of the dot construction | hole/ball = classical radius ratio under H1 (mirror walls) or H2 (insphere, equal dihedral angles); 20 cases; facet-transitive plus insphere is not enough | PROVEN under H1 or H2; VERIFIED on 20 polytopes | `make verify-radius-ratio` |
| CaB6 hard-sphere packing | maximum (pi/6)(30 sqrt2 - 41) = 0.746865 at gamma = sqrt2 - 1; beats fcc for 0.41023 < gamma < 0.41591 | PROVEN (upper bound); all-pairs attainment VERIFIED | `make verify-cab6` |
| 2D fold (hinged squares) | density 0.5390 to pi/4 at theta = pi/4 for the incentre | KNOWN linkage; corollaries PROVEN exactly | `make verify-fold2d` |
| 3D fold, cubic-restricted [111] | density 27 pi sqrt2/256 = 0.4686 at 60 degrees, six contacts; a lock of the cubic path, not a jam | VERIFIED | `make verify-fold3d` |
| **Fold density** (`autoresearch-fold`) | **pi/sqrt18 = 0.74048 (Kepler bound)**: sampled bar-preserving, overlap-free path in the 2x2x2 supercell to hcp (793 states) and, by two layer slides, to pure fcc (1193 states); both checkers accept both paths | VERIFIED (sampled continuity, step <= 0.01); fcc is NEW-COMPUTATION | `make verify-fold-cp`, `make verify-fold-tests` |

## Prior-art framing (updated 2026-10-12)

* Corner-sharing ReO3 reaches hcp by octahedral rotation alone, at 30 degrees about <111>: KNOWN (Johnsson and Lemmens 2006, textbook; Shan et al., under pressure). Only the cab result is "not found": with octahedra joined by bars, the cubic tilt stalls at 0.4686 (0.2637 at 30 degrees, against 0.7405 for hcp), and close packing needs a symmetry-breaking continuation with lattice shear.
* Closest physical study of the cab framework: Li, Refson and Dove, PCCP 25, 10749 (2023), rigid unit modes of the B6 octahedra in CaB6 (abstract read).
* CaB6 packing: the location of the maximum, gamma = sqrt2 - 1, may be known from Travesset, ACS Nano 11, 5375 (2017), which could not be read; the closed form is not confirmed new.
* 2D fold, related work with a different linkage: Tarnai, Fowler, Guest and Kovacs, Symmetry 14, 232 (2022).

## Open problems (fold)

* A certified path: interval arithmetic or an analytic piecewise path in place of sampled continuity.
* ~~A closed-form end state~~: exact integer closed forms for the hcp and fcc end states of the m = 2 paths, PROVEN (`results/fold-end-exact`, `make verify-fold-exact`). Still open: space group, Wyckoff sites and an orbit description of the hole pattern.
* The smallest supercell: m = 1 and odd m (m = 3); only m = 2 is known to reach close packing.
* A monotone-density path, and the minimum density dip.
* The set of close-packed end states and their connectivity; comparison nets (diamond, sodalite, ReO3 with bars).
* An outside review of both checkers, or a third independent implementation.

## Open problems (other)

* A written argument for the all-pairs non-overlap at the CaB6 optimum, and a proof that no non-bonded pair closes before pi/4 in 2D.
* Identify the 60-degree 3D net with RCSR `sms` by isomorphism.
* Classify the polytopes satisfying H1 or H2; the capped CN7 shapes under H2.
* Formalisation of the algebraic statements (`docs/formalization.md`).
* The kaleidoscope census, the boride generator and the triangle-centre perception experiment (see `tracks/README.md`, not replayed here).

## Housekeeping

* Not yet public: the paper sources (`paper/`, ignored by git). Not published: the commitment in `commitment/`.
* The site is separate from these checked results and is not covered by `make verify`.
