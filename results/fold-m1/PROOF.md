# The fold reaches close packing in the smallest cell (m = 1), and to fcc with density rising at every step

Evidence label: **VERIFIED (finite check, sampled continuity), the same kind of evidence as `results/fold-close-packing`;
end-state closed forms PROVEN exact (integer arithmetic, `results/fold-end-exact/exact_end.py`); that the stored final states
are these closed forms is VERIFIED (1e-7 after a rigid motion, `exact_end.exactify`). NEW-COMPUTATION.**
Replay: `make verify-fold-m1` (about 10 seconds; numpy, scipy). Certificate: `certificate.json`.

## Statement

Model as in `results/fold-close-packing/PROOF.md`, with the 1 × 1 × 1 supercell: one octahedron (6 balls), 15 bars (12 edges
and 3 linking bars, each joining the +e_i dot to the −e_i dot of the next cell, that is, to its own image). The lattice may shear.
First-order count: 6·3 + 9 − 15 − 3 − 3 = 6 degrees of freedom.

**Finding 1 (hcp, m = 1).** `path_hcp_m1.npz`: 333 states from the canonical open cubic state to hcp at density π/√18, with 12
contacts per ball, common-neighbour signature 36 × (4,2,1) + 36 × (4,2,2), and first-order jammed. It was found by free density
descent, with no fold stage and no target, so the search was not told what the end state should be.

**Finding 2 (fcc, m = 1, density rising at every sample).** `path_fcc_m1.npz`: 268 states from the canonical state to pure fcc
(72 × (4,2,1)) at density π/√18, with 12 contacts per ball and first-order jammed. det L strictly decreases at every one of the 267
steps, so the density rises at every sample. This path was found by steering the same LP engine toward an explicit fcc end
state written down in advance (below). The end state was the target; the path to it was not.

Both paths are accepted by the locked checker `autoresearch-fold/verify.py` and by the independent checker
`results/fold-close-packing/chk_independent.py`, with 9 interior sub-states per step (linear interpolation at spacing 1/10):

| Path | States | Max step (dots / lattice) | Min interpolated distance / 2R | Max interpolated bar error / 2R |
|---|---|---|---|---|
| hcp | 333 | 0.00421 / 0.00482 | 0.999996 | 1.9e-5 |
| fcc | 268 | 0.00421 / 0.00367 | 1.00011 | 1.6e-5 |

**End states (PROVEN, integer arithmetic).** In the cube frame of the octahedron, with the origin at its centre and units
K = 2 − √2, the dots sit at ±e_i, by definition of that frame. In the lab frame the octahedron does rotate (by about 0.44 rad on the hcp path and 0.61 rad on the fcc path, from a best-fit rotation of the 6 dots, start to end). The lattice rows are:

| End state | lattice rows | det L |
|---|---|---|
| fcc | (2, −1, −1), (1, 2, −1), (1, −1, 2) | 12 K³ |
| hcp | (1/3)·(7, −1, 4), (1/3)·(−4, 7, −1), (1/3)·(−3, 3, 6) | 12 K³ |

For fcc, the linking bars are T_a − 2e_a = (0, −1, −1), (1, 0, −1), (1, −1, 0), each of length √2 K = 2R. In both states,
`exact_end.check_exact` confirms in integers that every bar is exactly 2R, no pair is closer, there are 12 contacts per ball, the
density is exactly π/√18 (det L = 2·6·K³), and the structure is hcp or fcc by common-neighbour analysis.

**Stored final states (VERIFIED, 1e-7).** The closed forms above are PROVEN; that the stored paths end at them is a
floating-point match. `exact_end.exactify` moves the final state of each stored path by a rigid motion into the cube frame of
its octahedron and rounds it to the grid of units K/den (den = 3 for hcp, 1 for fcc), accepting it only if every coordinate is
within 1e-7 of that grid. `check_exact` then works on the rounded integer state.

## Consequences

* **Every supercell size, including odd m.** A state of the m = 1 framework is also a state of the m-periodic framework for every m.
  It is the same infinite configuration, and the m-supercell bar graph is the lift of the m = 1 graph. So a continuous m = 1 motion
  is a continuous m-periodic motion. Close packing (hcp and fcc) is therefore reachable for every m, including m = 3, with the
  same sampled-continuity caveat. The lifted paths are exact lifts and carry the same sampled-continuity evidence. Dot and lattice
  steps both grow with m, so the locked checker rejects the unsubdivided lifts at m = 2..4, on the step limit only. To be accepted
  they would need re-sampling, with new states projected back onto the bar constraints (linear subdivision breaks bars by up to
  1.9e-5 relative, against the checker's 1e-9). That is not done here, and m ≥ 5 is outside the locked checker's range (1..4).
* **The [111] checkerboard tilt is not needed.** Neither path uses stage A, and m = 1 cannot carry the checkerboard.
* **A density-monotone path exists** (at the sampled level), answering the question in `tracks/README.md`.

## Reproduce the paths

From the repository root, with `autoresearch-fold/.venv` set up as in its README (this overwrites `autoresearch-fold/search.py` and appends a line to
`autoresearch-fold/results.tsv`; restore both with `git checkout autoresearch-fold/search.py autoresearch-fold/results.tsv`):

    cp results/fold-m1/search_fcc_m1.py autoresearch-fold/search.py && (cd autoresearch-fold && ./run.sh)
    autoresearch-fold/.venv/bin/python -c "import numpy as n;a=n.load('autoresearch-fold/out/path.npz');b=n.load('results/fold-m1/path_fcc_m1.npz');assert int(a['m'])==int(b['m']) and n.array_equal(a['Xs'],b['Xs']) and n.array_equal(a['Ls'],b['Ls']);print('identical')"

The same goes for `search_hcp_m1.py` and `path_hcp_m1.npz`. Both search files pin BLAS to one thread. Each search is bounded by
a step count (`MAX_STEPS`) and a stall count (`STALL_KICKS`), not by time, and finishes in a few seconds. `WALL_LIMIT` (3000 s)
is a safety limit only: if it is hit, the search writes no path and exits with code 2 (`search INCOMPLETE`), so a slow machine
cannot silently produce a truncated path (`run.sh` stops any search at 300 s, also with no path written). On one 2-core Linux
machine (numpy 2.5.4, scipy 1.18.1) both reruns reproduced the stored paths byte for byte. That was with the earlier versions
of these files, which stopped on a 240 s time budget instead; on another machine (numpy 2.3.5 and 2.5.4) the current and earlier
versions wrote byte-identical paths. Other BLAS builds may round differently and give a different path, so the stored paths,
not the searches, are what is certified.

Search engine: the `autoresearch-fold` search.py pipeline (LP step that keeps bars and lets contacts separate, Gauss–Newton
retraction, sliding contacts), with a one-thread BLAS, a cycle breaker and, in the hcp file, three seeded start kicks. The fcc
file has no start kicks, and its run makes no random kick at all: it steers deterministically, with a stage whose objective is
the mismatch to the fcc lattice above, in metric and fractional coordinates, and skips the density descent once that target is
reached.

## What did not work (m = 2, 3, 4, same session)

* Free descent at m = 2 reached π/√18 for 4 of 10 seeds (all hcp). The others jammed at 0.6187 (4 seeds), 0.5977 and 0.5766. None of
  these free runs ended in fcc (session notes; not recorded in `autoresearch-fold/results.tsv` or replayed).
* Matching only the lattice metric of fcc, without fractional coordinates, fell into hcp.
* Free descent at m = 3 and m = 4 had not finished within the 300 s cap (0.549 and 0.297, still rising).

## Known limits

* Continuity is sampled (steps ≤ 0.01, midpoint check and 9 interior sub-states per step), as for the m = 2 paths. Linear
  interpolation changes bar lengths by up to 1.9e-5 relative, so the interpolated curve is not exactly bar-preserving. On the hcp
  path the interpolants also bring two non-bonded balls to 0.999996 · 2R (accepted under the checkers' 1e-4 interpolation
  tolerance), so the interpolated curve is not overlap-free either; the fcc interpolants stay above 2R (minimum 1.00011 · 2R).
  A certified path through these stored samples, and to hcp, is still open; an exact continuous motion to the same fcc end
  state (PROVEN, internal review only) is in `results/fold-continuous-fcc`.
* "Monotone" holds at the 268 samples. Between samples it is not proven.
* Both fold checkers were written inside the project; `exact_end.py` was written for these results. No outside review.
* Prior art: see `autoresearch-fold/PRIOR-ART-FOLD.md`. No separate search was run for the m = 1 case.

## Credit and AI assistance

The search engine and both checkers are the project's (`autoresearch-fold/`, `results/fold-close-packing/`). The m = 1 runs, the fcc
steering stage and this note were made with Claude Opus 5.5 (Claude Code) in a copy of `autoresearch-fold/`. No person has checked them yet.
