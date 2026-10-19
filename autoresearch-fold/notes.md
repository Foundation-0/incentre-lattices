# notes.md

## Iteration log (all numbers from results.tsv / verify.py)
- baseline (e270829): analytic cubic [111] fold to 60 deg, m=2: valid, 0.468585 = 27*pi*sqrt2/256, 6 contacts, flex 18.
- iter1 (d8144dc): fold to 60 deg, then LP-guided free-lattice continuation: valid, 0.710861, 11.5 contacts,
  flex 0, trigonal lattice (|a|=|b|=|c|=2.7476, cos = -3/22). KEPT. First-order jammed local optimum.
- iter2 (876afa7): skip fold, start continuation from canonical state with random kicks (seed 1): valid, 0.572893,
  6.08 contacts, flex 1. DISCARDED (stuck in a poor basin).
- iter3 (6e9ead9): iter1 on m=4: INVALID. Cause: stage-A step (0.002 rad) gives lattice change > 0.01 for m=4 (L scales
  with m); the checker caught it. Also m=4 is slow (~200 engine steps in 240 s). DISCARDED; fixed in iter4.
- iter4 (cae3bdb): stage-A step scaled 1/m; fold only to 30 deg then free continuation (seed 3), m=2: valid,
  0.740480 = pi/sqrt18 (Kepler density), 12 contacts per ball, flex 0, first-order jammed. KEPT (best).

## Characterisation of the iter4 end state
- Lattice (m=2 supercell, 48 balls): |a|=|b|=2.86976, |c|=2.70563, angles 90/90/120 -> hexagonal; det 19.297.
- Every ball has 12 contacts; common-neighbour analysis on the contact graph: each ball has 6 neighbours of type (4,2,1)
  and 6 of type (4,2,2) -> hexagonal close packing (hcp) for all 48 balls (fcc would be all (4,2,1)).
- Balls lie in 4 close-packed layers per c, layer spacing 0.6764 = sqrt(2/3)*2R: consistent with ABAB stacking.
- Density pi/sqrt18 = 0.74048 is the maximum for ANY sphere packing (Hales), so within this model nothing above it exists;
  iter4 reaches it with a valid continuous bar-preserving, overlap-free path.
- Not done: coordination sequence and RCSR lookup (hcp net expected); m >= 3: the only run in this log (iter3, m = 4) was INVALID
  on the step limit (later m = 3 and m = 4 free descents did not finish, see results/fold-m1/PROOF.md); for m = 1 and its lift to
  every m see results/fold-m1; first-order rigidity is not a proof of jamming.
