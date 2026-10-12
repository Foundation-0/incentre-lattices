# Open tracks

| Track | Question | Score | Checker |
|---|---|---|---|
| `fold-density` (folder `autoresearch-fold/`) | The 3D contact framework (rigid boron-like octahedra joined by bars) can fold. Starting from the incentre lattice, what is the highest packing density a continuous, bar-preserving, overlap-free motion can reach, and what structure is the end state? Baseline: 0.4686 on the cubic-restricted path. **Current best: π/√18 ≈ 0.74048 (Kepler bound), ending in hexagonal close packing with 12 contacts per ball, on a 2×2×2-cell supercell; accepted by the locked checker and confirmed by an independent from-scratch checker (`autoresearch-fold/independent/`). Prior-art search in progress. Open: larger supercells, fcc, a proof by explicit construction.** | max packing density | `autoresearch-fold/verify.py` |

## Proposed next tracks (checker not written yet; contributions welcome)

- **kaleidoscope-census**: run the equidistant-dot rule over every finite and affine Coxeter group in 2D to 4D. Which contact frameworks fold into a uniform packing, and is the incentre always the only position that does? Score: number of verified uniform end states.
- **boride-generator**: use the rule on Voronoi cells of other lattices and score candidate boron frameworks by their match to known structures in open crystal databases.
- **centre-perception**: an online experiment. Where do people tap when asked for "the middle" of a triangle: the incentre or the centroid? Needs a protocol, ethics notes and a data format before any data is collected.
