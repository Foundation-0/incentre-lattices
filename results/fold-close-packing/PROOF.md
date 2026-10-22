# The cab contact framework folds to close packing (m = 2 supercell)

Evidence label: **VERIFIED (finite check). The path is certified only by sampled continuity (consecutive states differ by at most 0.01
in every dot position and in the lattice) with a midpoint check and a 10-sub-state linear interpolation between every pair of states.
It is not a proof of a continuous overlap-free motion. The fcc extension is NEW WORK, checked by two checkers.**
Replay: `make verify-fold-cp` (slow, about 3 minutes; numpy, scipy) and `make verify-fold-tests` (the locked checker's unit tests).
Certificate: `certificate.json`. Stored paths: `autoresearch-fold/best/path.npz` (hcp, 793 states) and `path_fcc.npz` (this folder, 1193 states).

## Statement

Model. Cube cells of side 2 each hold six dots at the cell centre +- (2 - sqrt2) e_i; balls have radius R = sqrt2 - 1; bars of length
2R join the 12 octahedron edges of each cell and each +e_i dot to the -e_i dot of the neighbouring cell (the CaB6 / RCSR `cab`,
sphere-packing type 5/3/c3 contact framework, density pi R^3 = 0.2233). The lattice may shear. Work in a 2 x 2 x 2 periodic supercell:
48 balls, 120 bars. The generic Maxwell count is 27 (flexes minus self-stresses: 48*3 + 9 - 120 - 3 - 3). At the hcp and fcc end states the locked checker's
replay finds 29 first-order flexes of the bars alone (`flex_bars_only` in `certificate.json`, both end states), hence 2 self-stresses (120 bars minus rank 118).

**Finding (hcp, found by the autoresearch search in `autoresearch-fold/`).** There is a sequence of 793 states from the canonical
open cubic state, keeping every bar length to 1e-9 (relative) and with no two balls overlapping, ending in a hexagonal close packing:
density pi/sqrt18 = 0.74048 (the maximum for any sphere packing, Hales), 12 contacts per ball, common-neighbour signature 288 x (4,2,1)
and 288 x (4,2,2) (hcp), lattice 90/90/120 degrees, first-order jammed.

**Finding (fcc, NEW).** Appending two rigid layer slides (`make_fcc_path.py`, 200 states each, deterministic) to the hcp end state gives a
1193-state path ending in pure fcc: density pi/sqrt18, 12 contacts per ball, 576 x (4,2,1), no (4,2,2). Each slide moves the upper layers by
one close-packed stacking vector, lifted over the bridge site so the chosen interface bar keeps length 2R. The density dips to 0.7294
during the slides. A 4-layer mixed stacking (dhcp-like, 432 x (4,2,1) + 144 x (4,2,2)) is also reachable by other shift choices (see
`autoresearch-fold/PRIOR-ART-FOLD.md`; not part of the certificate).

## Why the end state is plausible (explanation, not proof)

In the hcp end state the balls lie in four close-packed layers per c. Each octahedron spans an adjacent layer pair, so octahedra sit in alternate
gaps only; in each such gap a third of the octahedral holes are chosen so that every ball belongs to exactly one octahedron. Of the 24 linking
bars, 16 lie inside a layer and the other 8 cross the two empty gaps, and all 8 point the same lateral way. That is why the layer pair above
a gap can slide by one stacking vector without losing those contacts. The end state is a combinatorial embedding of the cab bar graph into a
close-packed contact graph, easy to check once written down. The real content is the path: start and end lie in one connected component of the
overlap-free configuration space, which is not automatic. Greedy density ascent from other starting folds jams well below Kepler (0.58, 0.61,
0.62, 0.71), and the basin that works, fold angles about 25 to 35 degrees, is fragile (28 degrees fails). The stored `path.npz` is the certificate. The search is
not guaranteed to reproduce it: it is a numerical continuation with kicks (seeded), its original stage B stopped on a wall-clock deadline, and one slower machine
ended at a jammed state of density 0.618718. The current `search.py` stops on a step count (BLAS pinned to one thread) and reports "incomplete" if a wall-clock
safety limit is hit; one rerun gave a different valid path (791 states, not 793, diverging at state 285) that ends at the same density 0.740480. Rerunning it
adds no independent confirmation of the certificate.

## Two checkers

1. **Locked checker** `autoresearch-fold/verify.py` (sha256 pinned in `autoresearch-fold/verify.sha256`; `make verify-fold-tests` checks the pin and runs its 16
   unit tests, including must-fail cases: overlap, stretched bar, jump, wrong initial state, NaN, negative determinant, altered topology).
   It builds the bar graph itself, checks bar lengths, overlaps (with a midpoint check), step size at most 0.01, the initial state, and
   the end-state contacts and first-order rigidity.
2. **Independent checker** `chk_independent.py`, a re-implementation (numpy only, no import of `verify.py`) of
   `autoresearch-fold/independent/chk.py`, taking the path as an argument. It rebuilds the bar graph with fixed integer lattice shifts, checks all 120
   bar lengths in all states, the minimum non-bonded distance over all periodic images in all states (the image-shift block is derived from each state's lattice and dots, as in the locked checker's `_image_cloud`, not fixed), step sizes, a
   10-sub-state interpolation between every pair of states, and the density, contacts and common-neighbour analysis of the end state. On the hcp path
   it reproduces the numbers of the original `chk.py` (same maximum bar error 8.1e-14 relative, same minimum interpolated distance 0.999994 * 2R).
   Both this file and the path generator `make_fcc_path.py` are sha256-pinned in `checkers.sha256`; `make verify-fold-cp` checks the pins (and the locked checker's) before it runs, so weakening either file and rerunning `make certificates` fails `make verify` and CI. A deliberate change to either file means re-pinning by hand (`shasum -a 256 chk_independent.py make_fcc_path.py > checkers.sha256`) in the same commit.

## Results (all in `certificate.json`)

| Path | States | Locked checker | Independent checker | End state |
|---|---|---|---|---|
| hcp (`autoresearch-fold/best/path.npz`) | 793 | valid | accepted | density pi/sqrt18, 12 contacts, CNA 288 (4,2,1) + 288 (4,2,2) |
| fcc extension (`path_fcc.npz`) | 1193 | valid | accepted | density pi/sqrt18, 12 contacts, CNA 576 (4,2,1) |

**The independent check of the fcc path is new work.** Before this repository the fcc extension had been accepted only by `verify.py` plus a
CNA written by the same author. It now passes the independent checker too: bar errors below 1e-9 relative, zero overlaps in 1193 states, the
interpolated minimum non-bonded distance 0.999994 * 2R (attained in the hcp part, the slides stay above it), pure fcc by CNA. It is still
the same kind of evidence as the hcp path (finite sampling), and both checkers were written within one project; neither is an outside review.

Must-fail controls in the replay (each must be rejected by both checkers where applicable): one dot displaced by 0.05; the lattice scaled by
0.97 in one state; a state replaced by a half-way blend that breaks bars. The stored `path_fcc.npz` is regenerated from the hcp path in every
replay and must agree to 1e-9.

## Known limits

* Continuity is sampled: steps at most 0.01, overlap checked at the states, at midpoints and at 10 interpolation sub-states. Linear
  interpolation between states changes bar lengths by up to 1.5e-5 (relative to 2R), so the interpolated path is not bar-preserving; the
  states are. A rigorous path needs interval arithmetic, regularity of the constraint Jacobian between samples with a margin argument, or an
  analytic piecewise path.
* The density along the path is not monotone (smallest step -7e-5; dips to 0.7294 during the fcc slides).
* **This result is m = 2 only.** For m = 1 (and, by the lift argument, every m) see `results/fold-m1`.
* First-order rigidity at the end state is not a proof of jamming.
* This is a geometry and rigidity result. Real CaB6 and Prussian-blue analogues have bond-angle stiffness and do not behave like this
  model (CaB6 shows no structural transition to about 40 GPa in one study); no prediction about real materials is made.

## Prior art

Full notes in `autoresearch-fold/PRIOR-ART-FOLD.md` (opus novelty judgment, 2026-10-12). Summary:

* **NOT FOUND**: a pin-jointed contact framework of equal balls with the cab topology moving, with every bar length kept and no overlaps, from
  its open cubic state to a close packing. Not the same as proven new. This label covers only the cab result. The corner-sharing version is KNOWN
  and textbook (next bullets).
* **KNOWN, corner-sharing analogue**: in ReO3, rotating the octahedra alone reaches hcp of the anions (RhF3 type), at a rotation of 30 degrees
  about <111>. Johnsson and Lemmens, "Crystallography and Chemistry of Perovskites", Handbook of Magnetism and Advanced Magnetic Materials (Wiley 2006),
  doi:10.1002/9780470022184.hmm411, arXiv:cond-mat/0506606 (full text read): "The ReO3 structure type can be converted to a more dense packing by
  rotating the octahedra until a hexagonal close packing is obtained of the RhF3 type." Observed under pressure: Shan et al., arXiv:2304.09011,
  Natl. Sci. Rev. (2026), doi:10.1093/nsr/nwag425 (arXiv full text read): ReO3 reaches "a hexagonal close-packed arrangement at a rotation angle of 30
  degrees above 30 GPa", and further compression cannot come from rotation alone.
* **What differs in cab (ANALYSIS, this repo's data)**: rotation alone does not reach close packing. On the cubic-restricted [111] path
  (`results/fold-3d-cubic`) the density is 0.2637 at 30 degrees, below the hcp value 0.7405, and the path locks at 0.4686 at 60 degrees. Close packing
  needs a symmetry-breaking continuation with lattice shear (the search folds to 30 degrees, then lets the lattice go free). The 30 degrees
  of the search is a hand-over point, not the ReO3 hcp angle.
* **Closest physical study, KNOWN**: Li, Refson and Dove, Phys. Chem. Chem. Phys. 25, 10749 (2023), doi:10.1039/d3cp01306e (abstract read, full text not):
  rigid unit modes involving rotations of the B6 octahedra in CaB6 and the other cubic hexaborides, with very low flexibility of the connected
  octahedra network. A near-equilibrium lattice-dynamics calculation, not a finite-amplitude hard-contact path, so it does not change the NOT FOUND above.
* **Closest, KNOWN**: corner-sharing MX6 octahedra of ReO3 tilt about [111] in the a-a-a- pattern to the VF3 type with hcp anions (O'Keeffe and Hyde,
  Acta Cryst. B 33, 3802 (1977); Glazer, Acta Cryst. B 28, 3384 (1972)); FeF3 and GaF3 follow it under pressure (Jorgensen and Smith, Acta Cryst.
  B 62 (2006)). Same tilt as our stage A, but in the cab framework that stalls at 0.4686 and the symmetry-breaking continuation is needed.
* Method, KNOWN: sequential linear programming with a deformable cell (Torquato and Jiao, Phys. Rev. E 82, 061302 (2010)); first-order jamming test
  (Donev, Torquato, Stillinger and Connelly, J. Comput. Phys. 197, 139 (2004)); flexibility windows (Sartbaeva et al., Nature Mater. 5, 962 (2006)).
* Not searched, worth a library check before any public claim: Sowa's deformation papers, Koch and Fischer on deformations of 5/3 types,
  O'Keeffe and Hyde's *Crystal Structures I* on cab/ReO3, Conway and Sloane SPLAG chapter 1 on fcc to hcp layer slides.

## Open obligations

* A certified (interval or analytic) path; a closed-form end state (space group, Wyckoff sites, hole pattern).
* The minimal supercell (m = 1, 3), a monotone path, comparison nets (dia, sodalite, ReO3 with bars), connectivity of the set of close-packed end states.
* Review by someone outside the project; a third implementation of the checks.
