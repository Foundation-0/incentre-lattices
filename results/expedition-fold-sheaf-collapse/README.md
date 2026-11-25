# Universal Rank-1 Elastic Collapse, 4D Pure-Shear Duality, & 8D `F_8 -> E_8` Sheaf Bifurcation

**Status:** Independently **`VERIFIED`** (`AUDIT-019`, oracle: `exact`, novelty: **`N2`**) — 100% zero-float exact rational (`fractions.Fraction`, pure Python 3 stdlib, zero external dependencies).

## Summary of Exact Results

Building directly on the all-dimension `F_n -> D_n` (`C_n`, `C_4^nilp`) and canonical Sylvester–Hadamard `F_8 -> E_8` (`C_8 = (1/2) H_8`) bar-preserving folds, [`fold_sheaf_collapse.py`](./fold_sheaf_collapse.py) computes exact $q=0$ cellular sheaf cohomology and zero-float elastic Schur complements across dimensions $n = 3..8$:

1. **Exact Pure-$\mathbb{Q}$ Curve Parametrizations (Zero $\sqrt{2}$):**
   - For $F_n \to D_n$: $a(k) = \frac{4k}{2k^2+1}$, $t(k) = \frac{2k^2-1}{2k^2+1}$ satisfies $a^2 + 2t^2 \equiv 2$ identically in $\mathbb{Q}(k)$.
   - For all-spinor $c_{ii} = -1/2$ $F_8 \to E_8$: $a(m) = \frac{1-4m}{2m^2-m+1}$, $t(m) = \frac{1-2m^2}{2m^2-m+1}$ satisfies $a^2 - a t + 2 t^2 \equiv 2$ identically in $\mathbb{Q}(m)$.
2. **Universal Instantaneous Rank-1 Elastic Collapse Along $F_n \to D_n$ ($\forall n \ge 3$):**
   - Unfolded $F_n(0)$ has elastic rank $n$. Because $C_n$ ($\tau^* = (-1,1,\dots,1)$) and $C_4^{\mathrm{nilp}}$ ($\tau^* = (1,1,-1,-1)$) are uniquely diagonally symmetrizable, $\operatorname{rk}(M_{\mathrm{rot}}) = n - 1$, $q=0$ sheaf Betti numbers drop discontinuously to $H^0 = \frac{n(n-1)}{2}+1, H^1 = \frac{n(n-3)}{2}+1$, and macroscopic stiffness $S_{\mathrm{cell}}(s)$ collapses from rank $n$ to **exact rank 1** (`inertia = [1, 0, n(n+1)/2 - 1]`) for all $s \in (0, 1]$.
   - **Exact closed-form $s=1$ trace:** `4` at $n=3$ and $\frac{6n-2}{n}$ for $n \ge 4$ (`11/2, 28/5, 17/3, 40/7` for $n=4..7$, and `23/4` at $n=8$), with $S_{00} = 1/n$.
   - **Exact Master Rank-1 Macroscopic Stress Formula (Maxwell Periodic Stress at $L = 2I + C$):**
     $$\sigma^* = C^\top \operatorname{diag}(\tau^*) C + 2 \operatorname{diag}(\tau^*) C.$$
3. **4D Clifford-Nilpotent `C_4^nilp` Pure-Shear / Zero-Axial Exact Duality:**
   - At $s=1$, every axial component of $\sigma^*$ cancels identically (`axial_block_identically_zero = True`, `mode = [0,0,0,0,1,1,0,0,1,1]`, `tr(S_cell) = 4` vs $C_4$'s `11/2`) — 100% pure shear / zero bulk modulus, orthogonal dual to pure-axial $F_4(0)$.
4. **Full Jamming Transition & Exact 4D / 8D Elastic Isotropy:**
   - $D_3$ (`36` bars/cell): `(6, 0, 0)`, $C_{11}=1, C_{12}=C_{44}^{\mathrm{Voigt}}=1/2$ (Cauchy holds), Zener $A=2$.
   - $D_4$ (`96` bars/cell, identical for $C_4$ and $C_4^{\mathrm{nilp}}$): `(10, 0, 0)`, $C_{11}=3/2, C_{12}=C_{44}^{\mathrm{Voigt}}=1/2$, **exact 4D $O(4)$-isotropy $A(D_4) = 1$**.
5. **8D $F_8 \to E_8$ Elastic Bifurcation & $O(8)$-Isotropic Jamming:**
   - **Cyclic $\mathbb{Z}/16\mathbb{Z}$ $\ker(\phi)$ all-spinor bases:** Exhaustively **`0`** diagonally symmetrizable $\implies \operatorname{rk}(M_{\mathrm{rot}}) = 8/8$, $H^0=28, H^1=20$, **complete rank-0 floppy liquid collapse** (`S_cell === 0`) at both pure-$\mathbb{Q}$ interior `(16/11, 7/11)` and `s=1`.
   - **Canonical Sylvester–Hadamard $C_8 = \frac{1}{2} H_8$:** Symmetric ($\tau^* = \mathbf{1}$) $\implies \operatorname{rk}(M_{\mathrm{rot}}) = 7/8$, $H^0=29, H^1=21$, **exact rank-1 auxetic extremal collapse** with unique stress tensor $\sigma^* = 2 I_8 + H_8$ (`inertia = (4, 4, 0)`, eigenvalues $2 \pm 2\sqrt{2}$), persisting along the entire interior $s \in (0, 1]$.
   - **Jammed 240-root $E_8$ (spherical 7-design, Bravais vol = 1):** Full rank `(36, 0, 0)`, $C_{11}=9, C_{12}=C_{44}^{\mathrm{Voigt}}=3$, **exact 36x36 8D $O(8)$-isotropy $A(E_8) = 1$** (vs 112-root $D_8$'s $A(D_8) = 1/3$).
   - **Negative control:** Non-symmetrizable 3-cycle (`C3[0][1]=+1`) yields `rk(M_rot) = 3/3`, `H0 = 3, H1 = 0`, `C_inertia = [0, 0, 6]` (rank 0).

## Reproduction

```bash
python3 -I -O -B -W error results/expedition-fold-sheaf-collapse/fold_sheaf_collapse.py
python3 -I -O -B -W error results/expedition-fold-sheaf-collapse/fold_sheaf_collapse.py --json
```
