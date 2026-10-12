# Open tracks

| Track | Question | Score | Checker |
|---|---|---|---|
| `fold-density` (folder `autoresearch-fold/`) | The 3D contact framework (rigid boron-like octahedra joined by bars) can fold. Starting from the incentre lattice, what is the highest packing density a continuous, bar-preserving, overlap-free motion can reach, and what structure is the end state? Baseline: 0.4686 on the cubic-restricted path. **Current best: π/√18 ≈ 0.74048 (Kepler bound), ending in hexagonal close packing with 12 contacts per ball, on a 2×2×2-cell supercell; accepted by the locked checker and confirmed by an independent from-scratch checker (`autoresearch-fold/independent/`). No prior description found (searches listed in `autoresearch-fold/PRIOR-ART-FOLD.md`); the closest relative is the ReO3 to VF3 tilt toward hcp (O'Keeffe & Hyde 1977), which stalls at 0.4686 in this framework. Two rigid layer slides then reach fcc (accepted by two independent checkers). Geometry result, not a predicted phase transition in real materials. Open: see the list below.** | max packing density | `autoresearch-fold/verify.py` |

## Open problems on the fold

- ~~**A closed-form end state**~~: done for m = 2 (`results/fold-end-exact/`): both end states are exact, with every coordinate an integer multiple of (2 − √2)/3, and are checked in integer arithmetic. Remaining: space group and Wyckoff sites.
- **A certified path**: replace the numerical path with one whose bar lengths and non-overlap are proven, for example by interval arithmetic or an explicit construction.
- ~~**The smallest supercell**~~: m = 1 works (`results/fold-m1/`): 6 balls, free descent to hcp and a steered path to fcc, accepted by both checkers. An m = 1 motion is also an m-periodic motion for every m, so odd m is covered.
- ~~**Monotone density**~~: the m = 1 fcc path increases density at every one of its 268 samples (`results/fold-m1/`). Open: monotone between samples, and a monotone path to hcp.
- ~~**fcc second check**~~: done (`results/fold-close-packing/chk_independent.py`). Remaining: a certified, non-sampled path.

## Proposed next tracks (checker not written yet; contributions welcome)

- **kaleidoscope-census**: run the equidistant-dot rule over every finite and affine Coxeter group in 2D to 4D. Which contact frameworks fold into a uniform packing, and is the incentre always the only position that does? Score: number of verified uniform end states.
- **boride-generator**: use the rule on Voronoi cells of other lattices and score candidate boron frameworks by their match to known structures in open crystal databases.
- **centre-perception**: an online experiment. Where do people tap when asked for "the middle" of a triangle: the incentre or the centroid? Needs a protocol, ethics notes and a data format before any data is collected.
