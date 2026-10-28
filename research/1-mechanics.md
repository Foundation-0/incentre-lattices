# 1 · Mechanics and rigidity of the "dot" frameworks

Code: `research/mechanics/` (plain numpy/scipy/sympy, no venv needed).

| file | what it does |
|---|---|
| `pframe.py` | periodic rigidity matrix R(q) in n dimensions, Calladine counts, flexible-lattice (Guest–Hutchinson / Borcea–Streinu) analysis, supercells |
| `contact_nd.py` | contact framework of the truncated n-cubic honeycomb; mode counts for n = 2…6; t-scan |
| `strain_modes.py` | which macroscopic strains the pin-jointed contact framework allows at first order |
| `elastic.py` | effective elastic tensor of the spring network, including non-affine relaxation |
| `exact_checks.py` | sympy proofs of the elastic constants (n ≤ 7) and of the end states of the fold |
| `mechanism.py`, `final_states.py` | finite counter-rotation ("tilt") mechanisms, ball-collision limits, end-state contact numbers and densities |
| `wythoff.py`, `census.txt`, `census_summary.txt` | rigidity census of every Wythoff framework of the affine Coxeter groups C̃2, G̃2, Ã2, C̃3, B̃3, Ã3, C̃4, B̃4, D̃4 (159 frameworks) |
| `full_drawing.py` | frameworks (b) and (c): the whole drawing with the red lines, cube edges and pyramid edges |

Setup: cube side 2 centred at 0, lattice 2ℤⁿ. The dots sit at ±t eᵢ with t* = 2 − √2. **Contact framework** = bars between dots whose balls touch. These are the edges of the cross-polytope (octahedron in 3D) formed by the 2n dots of one cube, plus one "cross-face" bar per face. Every bar has length 2R = 2√2 − 2. Labels: **PROVEN**, **VERIFIED NUMERICALLY** (SVD rank, tolerance 1e-8, one or more random q: 3 for the n = 2…6 mode counts in `contact_nd.py`, one per framework in `wythoff.py` and `full_drawing.py`), **CONJECTURE**, **KNOWN**.

---

## Ranked findings

### 1. The drawing folds into the square ball packing: 4.8.8 → 4.4.4.4 at exactly 45° (2D), and a [111] fold in 3D that ends at 60° with every ball six-coordinated

Treat the contact graph as a pin-jointed mechanism (each ball is a joint at its centre and each touching pair is a rigid bar of length 2R between the centres; the balls may roll on one another. Pins at the contact points would give a different linkage, a = 2 cos θ; see `results/fold-2d/PROOF.md`). Each square of 4 dots is a cluster. Rotate the clusters alternately by +θ and −θ in a checkerboard. The cross-face bars stay parallel to the axes, and the lattice spacing becomes a(θ) = 2 − 2t(1 − cos θ). Every bar length is preserved exactly: this is a finite motion, not just an infinitesimal one (bar-length error about 1e-15 at the sampled angles 0°, 5°, …, 45°). The macroscopic strain is isotropic, εxx = εyy, so Poisson's ratio is −1 along the whole path.

* **PROVEN** (sympy, `results/fold-2d/fold2d.py`, asserted by its `verify.py`; also `exact_checks.py`): the gap between the balls of neighbouring clusters along the bar axis is 2√2·t·cos(θ + 45°), which vanishes at **θ = 45° exactly**. **VERIFIED** (finite check: 0.01° grid in `results/fold-2d/fold2d.py`; `mechanism.py` agrees on a ≈0.22° scan plus bisection): no other non-bonded pair touches before that, so 45° is the first new contact. At 45°, a = 4√2 − 4 = 4R and the strain is **−(√2−1)² = −R²**. At that moment every ball touches 4 others (it had 3), and the dots form a perfect **square lattice of spacing 2R**. The packing density rises from π(√2−1)² ≈ 0.539 to **π/4**, and the area shrinks by the factor 4R² = 12 − 8√2 ≈ 0.686.
* **PROVEN** (algebra): if a fold ends at 45°, its two edge lengths are √2·t and 2 − 2t, equal only at the incentre t* = 2 − √2.
* **VERIFIED** (finite scan, `results/fold-2d`: six dot positions c = 0.45, 0.5, 0.55, 2 − √2, 0.62, 0.7, 0.01° grid; not proved for every t): for t > 2 − √2 (0.62, 0.7) the fold ends at 45° with two different contact distances. For t < 2 − √2 it does **not** end at 45°: the fold continues and locks at 66.37°, 56.96° and 49.51° (c = 0.45, 0.5, 0.55) with three contacts per ball. So among the scanned values only the incentre ends in the regular square lattice, and in that sense the Wythoff point *is* mechanically singled out: the "uniform" start state folds onto the "uniform" end state.
* **VERIFIED NUMERICALLY** (3D, octahedral clusters, checkerboard tilts R / Rᵀ, `mechanism.py`, `final_states.py`):

| tilt axis | first new contact at | strain at end | volume change | contacts per ball at end | density at end |
|---|---|---|---|---|---|
| [001] (stacked 2D fold) | 45° | (−R², −R², 0) | −0.3137 | 6,6,6,6,5,5 | 0.3253 |
| **[111]** | **60.000°** | isotropic, −0.21895 each axis | −0.5235 | **6 for every ball** | **0.4686** |
| [011] | 54.7356° = arccos(1/√3) | (−0.2476, −0.1427, −0.1427) | −0.4470 | 5,5,6,6,6,6 | 0.4038 |

  For the [111] tilt the end lattice spacing is exactly **8t*/3** (PROVEN). The identity behind it, 4t/3 = √(L² − 2t²/9), is equivalent to √2·t = L, which is exactly the incentre condition. So the 3D fold also closes neatly only at the Wythoff point. The [111] fold is the analogue of the perovskite a⁻a⁻a⁻ tilt, and it is isotropic (ν = −1, cubic). Its end state has z = 6 = 2n, the Maxwell isostatic count.
* Novelty search: "truncated square tiling mechanism", "rotating squares bars/struts/ligaments fully closed", "square lattice of disks to truncated square tiling packing", "a⁻a⁻a⁻ maximal tilt". The 2D linkage is **KNOWN**: A. Flores, "Hinged tilings", North American GeoGebra Journal 6(1), §2.9, Fig. 9 ("hinged squares with rods"), crediting D. Stuart (1963); see `results/fold-2d/PROOF.md`. Its ball-packing numbers (45°, −R², π/4) are elementary corollaries (PARTIALLY KNOWN). **Not found**: the 3D [111] lock state (60°, density 0.4686; `results/fold-3d-cubic`). Closely related **KNOWN** ideas: rotating-squares auxetics (Grima & Evans 2000, ν = −1); octahedral tilting with limits set by anion–anion contact in perovskites (Glazer 1972; O'Keeffe & Hyde 1977, whom I did not retrieve); contraction by rotation as a negative-thermal-expansion mechanism in linearly bridged frameworks (Goodwin 2006). The fold is very much "rotating units with split hinges".

### 2. Exact n-dimensional mode count and exact elastic constants of the contact framework

**Mode count (PROVEN for n ≥ 3; VERIFIED NUMERICALLY for n = 2…6):**
* At a generic wavevector there are exactly **n(n−1)/2 zero modes per cell**, and they are the infinitesimal rotations of the 2n-dot cross-polytopes. They form **flat bands**. There are **n(n−3)/2 states of self-stress** per cell (0 in 3D, 2 in 4D, 5 in 5D, 9 in 6D), all living inside the cross-polytopes.
* On each plane qᵢ = 0 there is one extra mode (rigid sliding of rows along axis i). At q = 0 there are n(n+1)/2 zero modes (the Euclidean motions of the clusters) and n(n−1)/2 self-stresses.
* Proof: for n ≥ 3 the cross-polytope skeleton is infinitesimally rigid (it is a simplicial polytope; Cauchy–Dehn–Alexandrov–Whiteley). So in a zero mode each cluster moves rigidly, as a translation v_c plus a rotation ω_c. The cross bar along eᵢ lies on the line joining the two cluster centres, so ω drops out (eᵢ·ω eᵢ = 0). The bar condition becomes eᵢ·(v_{c+2eᵢ} − v_c) = 0, which in Bloch form is (e^{2iqᵢ} − 1)vᵢ = 0. The self-stress count then follows from Calladine. n = 2 is the exception because the 4-bar square is not rigid: there are 2 zero modes per cell and no self-stress.
* The counts do not depend on where the dot sits on the axis (t ∈ (0, 1)). **The incentre does not change the mode count.**
* **KNOWN in 3D**: this is the Prussian-blue topology (octahedra joined by linear bridges). Goodwin, *PRB* 74, 134302 (2006), found an O(V) density of rigid unit modes, i.e. flat bands. Borcea and Streinu (2021), Section 7, treat the same periodic graph in any dimension (bars at t = √2−1, ours at 2−√2; our counts hold for any 0<t<1) but give only the 3D count. The n-dimensional counts are likely routine with existing machinery (priority unresolved); the unit-spring C11/C12 = (n−1)² and ν100 are the strongest remaining candidate.

**Strains (VERIFIED NUMERICALLY, `strain_modes.py`):** let the lattice deform. In n ≥ 3 the first-order mechanisms produce exactly the n(n−1)/2 pure shears, and every normal strain is blocked. In 2D they produce the 2-dimensional deviatoric space, so only the bulk modulus is non-zero: a fluid-like, "2D pentamode" material.

**Elastic constants (PROVEN for every n ≥ 3 by an explicit minimisation, `results/rigidity-nd`; the same algebra gives ν = 1 formally at n = 2; numerics agree for n = 2…8):** take unit springs on every bar and the symmetric relaxation in which vertices move only along their own axes. Then

  C₁₁ = 4(n−1)²/(n²2ⁿ), C₁₂ = 4/(n²2ⁿ), C₄₄ = 0, so **ν₁₀₀ = 1/(n² − n − 1)**: 1, 1/5, 1/11, 1/19, 1/29, 1/41.

In 3D, **ν = 1/5 exactly**, and it does not depend on t for unit springs. For same-material bars (stiffness proportional to 1/L), ν(t) rises monotonically from 0.035 to 0.32. The value at the incentre equals the unit-spring value 1/5 only because all bars there have equal length; it is not an extremum. The general-n step is a linear solve: s_i = (2εᵢ − Σε/n)/n. **No prior reference found** for these closed forms.

### 3. Wythoff census: the incentre is mechanically ordinary, and simplicial vertex figures give pentamodes (`wythoff.py`)

For 9 of the 11 irreducible affine Coxeter groups in dimensions 2 to 4 (C̃2, G̃2, Ã2, C̃3, B̃3, Ã3, C̃4, B̃4, D̃4; Ã4 and F̃4 were not run) and every nonempty ringing of their diagrams (159 frameworks = 3·7 + 3·15 + 3·31), I compared the uniform (incentre-type) point with a random point on the same face of the fundamental simplex. **The counts were identical in every case**: generic-q zero modes and self-stresses, q = 0 counts, and flexible-lattice strain dimension (VERIFIED NUMERICALLY). How much that shows needs care. In the 36 single-ringed cases the face is a single vertex of the simplex, so the "random" point is the incentre-type point and the agreement is automatic. Ringings related by a symmetry of the diagram are counted separately (for C̃2, the rings (0,) and (2,) have identical rows in `census_summary.txt`). The non-trivial comparison is therefore the 123 multi-ringed cases, which drop to 70 once each symmetry class is counted once (91 classes in all, of which 21 are single-ringed). Moving the dot off the incentre without changing which mirrors it lies on never changes the linear mechanics. The incentre extremises packing quantities, not mode counts.

Patterns:
* Every framework with **z = n+1** (all 2D truncated and omnitruncated tilings; in 3D the omnitruncated honeycombs, C̃3 rings (1,2) and (1,2,3), and the runcicantic cubic honeycomb, B̃3 rings (0,1,2) and (0,1,3); in 4D the z = 5 cases) has 0 self-stress at generic q and **exactly 1 at q = 0**. Its achievable strain space has dimension n(n+1)/2 − 1: only dilation is blocked, i.e. a pentamode (3D) or bulk-only (2D) material.
  * PROVEN for the omnitruncated frameworks: the stress with tension proportional to the area of the facet of the fundamental simplex that each bar crosses is in equilibrium at every vertex, by Minkowski's relation Σ|Fᵢ|nᵢ = 0. Uniqueness is VERIFIED NUMERICALLY.
  * The pentamode principle itself (z = n+1 junctions, one stress) is **KNOWN**: Milton & Cherkaev 1995; Kadic et al. 2012. Its systematic appearance across all Wythoff omnitruncations: no census found.
* Unusual flat-band Maxwell cases: G̃2 ring (0,2) (the 3.4.6.4 tiling, the only 2D case), C̃3 rings (0,2) and (1,3) (cantellated cubic, which recurs as B̃3 (0,2,3)), B̃3 rings (0,2) and (0,3) (runcic cubic, a different honeycomb), C̃4 (1,3), C̃4 (0,2,4), B̃4 (1,3,4) and D̃4 (0,1,3,4) have z = 2n, but N₀ = S > 0 at *every* q. They are "Maxwell lattices with flat zero-mode/self-stress bands". The pyrochlore and kagome cases (Ã3 (0,2), G̃2 (1)) are KNOWN to sit on the Maxwell line, but they share only z = 2n with these: their generic-q counts are N₀ = S = 0.
* No compact formula in n emerged for the whole census beyond the contact family above and the z = n+1 family, where generic zero modes per cell = nV − E = (n−1)|W/T|/2 for omnitruncations.

Novelty search ("rigidity Coxeter/kaleidoscope/uniform honeycomb census", "Archimedean lattices Maxwell 4.8.8 3.12.12"): no census of uniform honeycomb rigidity found. Related KNOWN work: flexibility of periodic tilings by regular polygons (Grasegger & Legerský, arXiv 2305.01570); symmetry-adapted periodic counts (Power, arXiv 1103.1914; Ross; Borcea & Streinu).

### 4. The other frameworks you asked about are rigid and over-braced (`full_drawing.py`; VERIFIED NUMERICALLY)

* (b) The full drawing (contacts + red lines + cube edges + pyramid edges) is rigid with 0 generic zero modes in 2D, 3D and 4D. Generic self-stresses per cell: 12 in 2D, 32 in 3D, 80 in 4D. With a deformable lattice there are no mechanisms at all.
* (c) The pyramid / triangle skeleton is rigid even without face diagonals. Generic self-stresses per cell: 2 in 2D, 5 in 3D, 12 in 4D.
* (d) The dot octahedra as rigid units joined by cross-face bars is exactly the 3D contact framework of finding 2, because a triangulated octahedron is rigid.

### 5. Predictions for a 3D-printed version (CONJECTURE, from the above)

Suppose balls are fused at small necks, so the necks act as flexural hinges.
* **Extreme cubic anisotropy**: C₄₄ comes only from neck bending, while C₁₁ and C₁₂ come from neck stretching, so C₄₄/C₁₁ ~ (neck radius / ball radius)². The Zener ratio A = 2C₄₄/(C₁₁ − C₁₂) should be ≪ 1.
* **Tension/compression asymmetry**: the axes are straight bar chains, stiff in tension. Under compression the fold mode has zero first-order cost, so the lattice should buckle into the counter-rotated pattern and switch to ν ≈ −1, the 3D analogue of Bertoldi's buckled holey sheet.
* In 2D the printed 4.8.8 sheet should start with ν near +1 (bulk-only) and then, after buckling, lock into the square packing at a −17.2% (= −R²) linear strain.

---

## Most promising open question (PhD-sized)

**Map the full configuration space of the ball/bar contact framework, as a real algebraic variety, for every Wythoff orbit.** Questions:
* Which tilt paths, in Glazer-like notation for split-hinge frameworks, are finite?
* Where do they end in new contact states? In 2D the end state is the square packing; in 3D the cubic-lattice [111] fold stalls at a z = 6 packing of density 0.4686, a waypoint and not a jam: if the lattice may shear, the machine search in `autoresearch-fold/` continues to hexagonal close packing, density π/√18 = 0.74048 (fcc follows by two layer slides; `results/fold-close-packing`, VERIFIED, first-order rigid). A trigonal grid scan alone already reaches 0.6477 (`results/fold-3d-cubic`).
* Is the incentre always the unique point whose folds end in a uniform (single edge length) packing?

The 2D result ("4.8.8 folds isometrically onto the square packing at 45°; among the scanned dot positions, only at the incentre") suggests a general principle: Wythoff orbits and the contact graphs they create are linked by finite mechanisms. That would connect the radius-ratio and packing side of this project to rigid-unit-mode crystallography and to transformable metamaterials.

## Honesty notes
* All 3D tilt results use a 2×2×2 checkerboard with lattice vectors kept on the axes. Letting the lattice shear is a separate branch of the mechanism: its continuation is the hcp and fcc path of `results/fold-close-packing`, which is checked on the 2×2×2 supercell only (larger blocks and other stackings are untested).
* The ball-collision test checks the nearest non-bonded neighbours of one cluster in a 5ⁿ block. That block is enough for the end states found, but it is not a proof of global non-overlap along the path. For 2D the gap formula and its zero at 45° are proven analytically; that no other pair touches earlier is checked on a 0.01° grid.
* Prior-art searches were web searches only, not a full literature review. The 2D linkage was later found in print (Flores, after Stuart 1963) and is KNOWN; see `results/fold-2d/PROOF.md`.
