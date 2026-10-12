# The fold reaches close packing in the smallest cell (m = 1), and to fcc with density rising at every step

Evidence label: **VERIFIED (finite check, sampled continuity), the same kind of evidence as `results/fold-close-packing`;
end states PROVEN exact (integer arithmetic, `results/fold-end-exact/exact_end.py`). NEW-COMPUTATION.**
Replay: `make verify-fold-m1` (about 10 seconds; numpy, scipy). Certificate: `certificate.json`.

## Statement

Model as in `results/fold-close-packing/PROOF.md`, with the 1 × 1 × 1 supercell: one octahedron (6 balls), 15 bars (12 edges
and 3 linking bars, each joining the +e_i dot to the −e_i dot of the next cell, that is, to its own image). The lattice may shear.
First-order count: 6·3 + 9 − 15 − 3 − 3 = 6 degrees of freedom.

**Finding 1 (hcp, m = 1).** `path_hcp_m1.npz`: 333 states from the canonical open cubic state to hcp at density π/√18, with 12
contacts per ball, common-neighbour signature 36 × (4,2,1) + 36 × (4,2,2), and first-order jammed. It was found by free density
descent, with no fold stage and no target, so the search was not told what the end state should be.

**Finding 2 (fcc, m = 1, density never falls).** `path_fcc_m1.npz`: 268 states from the canonical state to pure fcc
(72 × (4,2,1)) at density π/√18, with 12 contacts per ball and first-order jammed. det L strictly decreases at every one of the 267
steps, so the density rises at every sample. This path was found by steering the same LP engine toward an explicit fcc end
state written down in advance (below). The end state was the target; the path to it was not.

Both paths are accepted by the locked checker `autoresearch-fold/verify.py` and by the independent checker
`results/fold-close-packing/chk_independent.py`, with 10 interpolated sub-states per step:

| Path | States | Max step (dots / lattice) | Min interpolated distance / 2R | Max interpolated bar error / 2R |
|---|---|---|---|---|
| hcp | 333 | 0.00421 / 0.00482 | 0.999996 | 1.9e-5 |
| fcc | 268 | 0.00421 / 0.00367 | 1.00011 | 1.6e-5 |

**End states (PROVEN, integer arithmetic).** In the cube frame of the octahedron, with the origin at its centre and units
K = 2 − √2, the dots sit at ±e_i, so the octahedron ends unrotated in the octahedron's own frame. In the lab frame it does rotate (by about 0.44 rad on the hcp path and 0.61 rad on the fcc path, from a best-fit rotation of the 6 dots, start to end). The lattice rows are:

| End state | lattice rows | det L |
|---|---|---|
| fcc | (2, −1, −1), (1, 2, −1), (1, −1, 2) | 12 K³ |
| hcp | (1/3)·(7, −1, 4), (1/3)·(−4, 7, −1), (1/3)·(−3, 3, 6) | 12 K³ |

For fcc, the linking bars are T_a − 2e_a = (0, −1, −1), (1, 0, −1), (1, −1, 0), each of length √2 K = 2R. In both states,
`exact_end.check_exact` confirms in integers that every bar is exactly 2R, no pair is closer, there are 12 contacts per ball, the
density is exactly π/√18 (det L = 2·6·K³), and the structure is hcp or fcc by common-neighbour analysis.

## Consequences

* **Every supercell size, including odd m.** A state of the m = 1 framework is also a state of the m-periodic framework for every m.
  It is the same infinite configuration, and the m-supercell bar graph is the lift of the m = 1 graph. So a continuous m = 1 motion
  is a continuous m-periodic motion. Close packing (hcp and fcc) is therefore reachable for every m, including m = 3, with the
  same sampled-continuity caveat. The locked checker's absolute step limit (0.01 on L, which scales with m) means the stored
  m = 1 paths would have to be subdivided before that checker accepts them at larger m. That is not done here.
* **The [111] checkerboard tilt is not needed.** Neither path uses stage A, and m = 1 cannot carry the checkerboard.
* **A density-monotone path exists** (at the sampled level), answering the question in `tracks/README.md`.

## Reproduce the paths

From the repository root, with `autoresearch-fold/.venv` set up as in its README (this overwrites `autoresearch-fold/search.py` and appends a line to
`autoresearch-fold/results.tsv`; restore both with `git checkout autoresearch-fold/search.py autoresearch-fold/results.tsv`):

    cp results/fold-m1/search_fcc_m1.py autoresearch-fold/search.py && (cd autoresearch-fold && ./run.sh)
    cmp <(python3 -c "import numpy as n;d=n.load('autoresearch-fold/out/path.npz');print(d['Xs'].tobytes().hex())") \
        <(python3 -c "import numpy as n;d=n.load('results/fold-m1/path_fcc_m1.npz');print(d['Xs'].tobytes().hex())")

The same goes for `search_hcp_m1.py` and `path_hcp_m1.npz`. Both search files pin BLAS to one thread, and each finishes in about
10 seconds, well inside its time budget. On one 2-core Linux machine (numpy 2.5.4, scipy 1.18.1) both reruns reproduced the
stored paths byte for byte. Other BLAS builds may round differently and give a different path, so the stored paths, not the
searches, are what is certified.

Search engine: the `autoresearch-fold` search.py pipeline (LP step that keeps bars and lets contacts separate, Gauss–Newton
retraction, sliding contacts), with a one-thread BLAS, a cycle breaker and three seeded start kicks. The fcc file adds a steering
stage whose objective is the mismatch to the fcc lattice above, in metric and fractional coordinates.

## What did not work (m = 2, 3, 4, same session)

* Free descent at m = 2 reached π/√18 for 4 of 10 seeds (all hcp). The others jammed at 0.6187 (4 seeds), 0.5977 and 0.5766. None of
  14 free runs ended in fcc.
* Matching only the lattice metric of fcc, without fractional coordinates, fell into hcp.
* Free descent at m = 3 and m = 4 had not finished within the 300 s cap (0.549 and 0.297, still rising).

## Known limits

* Continuity is sampled (steps ≤ 0.01, midpoint and 10 sub-state checks), as for the m = 2 paths. Linear interpolation changes
  bar lengths by up to 1.9e-5 relative, so the interpolated curve is not exactly bar-preserving. A certified path is still open.
* "Monotone" holds at the 268 samples. Between samples it is not proven.
* Both fold checkers were written inside the project; `exact_end.py` was written for these results. No outside review.
* Prior art: see `autoresearch-fold/PRIOR-ART-FOLD.md`. No separate search was run for the m = 1 case.

## Credit and AI assistance

The search engine and both checkers are the project's (`autoresearch-fold/`, `results/fold-close-packing/`). The m = 1 runs, the fcc
steering stage and this note were made with Claude Opus 5.5 (Claude Code) in a copy of `autoresearch-fold/`. No person has checked them yet.
