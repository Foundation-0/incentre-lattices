# Exact Pure-$\mathbb{Z}/\mathbb{Q}$ & $\mathbb{Q}[\sqrt{2},\sqrt{3}]$ Maraldi Rhombic-Dodecahedron Certifier (`Foundation-0/incentre-lattices`)

Strictly add-only, zero-float, zero-external-dependency exact algebraic certifier complementing `results/radius-ratio/radius_ratio.py` (lines `73, 134, 153`), `research/chemistry/borides_ideal.py` (lines `59-78`), and the Maxwell-trace prestress obstruction in Draft `PR #109` on `Foundation-0/incentre-lattices` (`origin/main @ 2d374fff2a59cc93`).

> **Scope & Novelty Discipline:** Zero existing repository files are modified (`0` overlap with `#105` / `#109`). Classical geometry builds on Maraldi (1712) and Fejes Tóth (1964) for the rhombic-dodecahedral honeycomb base angle ($\cos\theta = -1/3$), paired with exact Farkas LP primal-dual certification of the incentre-lattice rhombic-dodecahedron pyramid insphere (`hole/ball = 1`, `UB12 x_ideal = 1/6`) replacing float64 `scipy.optimize.linprog` / `ConvexHull` approximations.

## What Is Proven Exactly (Zero Floats, Zero Bare `assert`s)

1. **Part A — Exact 14-Vertex Rhombic Dodecahedron Geometry, `H1` Mirror Walls & Pyramid-Insphere Farkas Duality:**
   - Exact vertex set (`8` cube vertices `(±1,±1,±1)` + `6` octahedron vertices `(±2,0,0)` permutations), `24` edges (`|e|^2 = 3`), `12` faces (`n_F · x <= 2`, `|n_F|^2 = 2`, Euler $\chi = 14 - 24 + 12 = 2$), orthogonal face diagonals (`|d_long|^2 = 8`, `|d_short|^2 = 4`, exact ratio `2`), obtuse Maraldi angle $\cos\theta = -1/3$, adjacent normal dihedral $\cos\phi = 1/2$, and all `24` `H1` ridge mirror reflections $v' = v - (v \cdot w) w$ preserving the vertex set identically.
   - **Primal + Farkas Dual Pyramid-Insphere Certificate (replaces `scipy.optimize.linprog` at `radius_ratio.py:73`):** For every face $F$, dot center $d_F = \frac{2}{3} n_F$ satisfies exact base slack $2 - n_F \cdot d_F = \frac{2}{3}$ and all four ridge-wall slacks $-(n_g - n_F) \cdot d_F = \frac{2}{3}$ with $|n_F|^2 = |n_g - n_F|^2 = 2$, certifying exact squared insphere radius $r^2 = 2/9$ and dot center squared radius $D^2 = 8/9 = 4 r^2$. Dual multipliers $\lambda_{\text{base}} = 1/3$, $\lambda_{\text{wall}} = 1/6$ ($\times 4$) satisfy $\sum \lambda_i = 1$ and $\frac{1}{3} n_F + \frac{1}{6} \sum_{k=1}^4 w_k = (0,0,0)$ identically — proving global optimality with zero numerical tolerance!
   - **Exact `hole/ball` & `UB12` Ideal Parameter (`borides_ideal.py:59-78`):** $D^2 / r^2 = 4 \implies D = 2r \implies (D - r)/r = 1$ exactly, cubic `UB12` boron coordinate $x_{\text{ideal}} = (2/3)/4 = 1/6$ exactly, `24` exact adjacent kissing contacts ($|d_a - d_g|^2 = 8/9 = 4 r^2$), and `42` non-adjacent pairs bounded below by exact $16/9 = 8 r^2$ (`0` overlaps).
2. **Part B — Calculus-Free $\mathbb{Q}[\sqrt{2},\sqrt{3}]$ Maraldi Surface-Area Optimality:**
   - Exact field multiplication and exact sign decision over $\mathbb{Q}[\sqrt{2},\sqrt{3}]$ certify flat base area $S_{\text{flat}} = \frac{3}{2}\sqrt{3}$, Maraldi three-rhombus base area $S_{\text{maraldi}} = \frac{3}{2}\sqrt{2}$, positive area saving $\Delta S = \frac{3}{2}(\sqrt{3} - \sqrt{2}) > 0$, exact squared ratio $S_{\text{maraldi}}^2 / S_{\text{flat}}^2 = 2/3$, and strict positive excess $1 + 8x^2 - 4x\sqrt{2} > 0$ across rational apex perturbations $x \ge 0$.
3. **Part C — Unpinned Maxwell-Trace Obstruction (`PR #109`) vs. Pinned-Boundary 3D Tension-Only RD Cable Web:**
   - Unpinned `14`-vertex / `24`-edge RD has Maxwell trace $\sum_{e} \omega_e |e|^2 = 3 \sum_e \omega_e = 72\,\omega_{\text{rd}} \ne 0$ for any uniform tension $\omega_{\text{rd}} > 0$ — algebraically obstructing free-floating tension-only self-stress.
   - Tethering all `14` interior RD vertices radially to `14` pinned outer anchors at $2v$ (`38` tension-only cables: $\omega_{\text{rd}} = 1$, $\omega_{\text{cube}} = 1$, $\omega_{\text{oct}} = 2$) yields exact zero force residual at all `14` interior nodes and strictly positive-definite interior stress matrix $\Omega_{\text{int}}$ verified via exact-$\mathbb{Q}$ Bunch–Kaufman Sylvester inertia `(14, 0, 0)` (contrasted with compression negative control `(-1, -1, -2)` giving `(0, 14, 0)`).

## Files

- `maraldi_rd_exact_audit.py` — Pure-$\mathbb{Z}/\mathbb{Q}$ and $\mathbb{Q}[\sqrt{2},\sqrt{3}]$ certifier, canonical zero-float certificate emitter, and fail-closed `--check` verifier.
- `test_maraldi_rd_exact_r5.py` — Assertion-free `unittest` R5 suite (positive proofs, negative controls, adversarial JSON/symlink/cap/CLI hardening).
- `maraldi_rd_exact_certificate.json` — Committed canonical zero-float JSON certificate (byte-identical to `--standalone --json`).

## Independent Reproduction (Isolated Mode)

```bash
python3 -I -B studies/maraldi-rhombic-dodecahedron-exact/maraldi_rd_exact_audit.py --standalone --json | cmp -s - studies/maraldi-rhombic-dodecahedron-exact/maraldi_rd_exact_certificate.json
python3 -I -O -B studies/maraldi-rhombic-dodecahedron-exact/maraldi_rd_exact_audit.py --check studies/maraldi-rhombic-dodecahedron-exact/maraldi_rd_exact_certificate.json
python3 -I -O -B studies/maraldi-rhombic-dodecahedron-exact/test_maraldi_rd_exact_r5.py -v
```
