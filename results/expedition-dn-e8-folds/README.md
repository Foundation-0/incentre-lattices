# All-Dimension Continuous Bar-Preserving Folds $F_n \to D_n$ ($\forall n \ge 3$), Optimal 4D $D_4$, Canonical Sylvester–Hadamard 8D $F_8 \to E_8$, and Exact 16D Barnes–Wall $\Lambda_{16}$ & 24D Leech $\Lambda_{24}/\sqrt{2}$ Classification

Self-contained companion study (`results/expedition-dn-e8-folds/`) extending `Foundation-0/incentre-lattices` (`results/fold-continuous-fcc/`) from 3D $F_3 \to D_3$ (FCC) to **all dimensions $n \ge 3$ ($F_n \to D_n$)**, **optimal 4D $D_4$ via Clifford-nilpotent $C_4^{\mathrm{nilp}}$**, **optimal 8D $E_8$ (Viazovska 2016) via canonical Sylvester Hadamard $H_8 = H_2^{\otimes 3}$**, and the **complete exact feasibility-and-obstruction classification through 16D ($\Lambda_{16}$ / $E_8\oplus E_8$) and 24D optimal Leech $\Lambda_{24}/\sqrt{2}$**.

All results below are pure Python 3 stdlib, zero floating-point arithmetic on any proof/verification path, zero bare `assert`s, zero CRLF, and independently audited as **`VERIFIED` | `exact` | `N2` CONFIRMED** (`AUDIT-009`, `AUDIT-017`, and `AUDIT-020R`).

## 1. Universal All-Dimension Continuous Folding Theorem $F_n \to D_n$ ($\forall n \ge 3$) & 4D Clifford-Nilpotent $C_4^{\mathrm{nilp}}$ (`AUDIT-009`, `N2`)

- **Proof:** [`ALL_N_DN_FOLD_PROOF.md`](./ALL_N_DN_FOLD_PROOF.md) (complete computer-free analytical proof of Lemmas 1–5 for every $n \ge 3$).
- **Verifier & Certificate:** [`verify_all_n_dn_fold.py`](./verify_all_n_dn_fold.py) $\to$ [`all_n_dn_fold_certificate.json`](./all_n_dn_fold_certificate.json).
- **Key Theorems:**
  - **Twisted $n$-cycle $C_n$:** Exact Schur–Chebyshev determinant $D_n(d, t) = d p_{n-1} + 2 t^2 p_{n-2} + 2 t^n$ ($\det L_n(1) = 4n$) strictly decreasing on $[0, 1]$ ($D_n'(t) \le (1 - 3\sqrt{2}) t p_{n-2} < 0$), exact rational $LDL^T$ min pivot $n/(n-1)$ on $2I_n + S_n$, dual characters $(\chi_1, \chi_2)$ proving $P_n + \mathbb{Z}^n L_n(1) = D_n + e_0$ (kissing number $2n(n-1)$), and exact discrete path Sobolev + SOS proof of degree-4 $\mathbb{Z}[\sqrt{2}]$ Bernstein nonnegativity for **all $n \ge 3$**.
  - **Clifford-Nilpotent Optimal 4D $D_4$ ($C_4^{\mathrm{nilp}}$):** Unique dimension $n = 4$ solving $2^n = 4n$; $C_4^2 = 0$, $I_4 + S_4 = 2\Pi_{\mathrm{orth}}$ (rank-2 projector), $\det L_4(s) = d(s)^4$, kissing number $24$ = optimal 4D $D_4$ density $\pi^2/16$.
  - **Exhaustive zero-float Bernstein certificates:** `universal_C3` (`3,248` pairs / `16,240` coeffs, exact match to `fold-continuous-fcc`), `clifford_nilpotent_C4` (`21,944` pairs / `109,720` coeffs, `0` negative), `universal_C4` (`35,144` pairs / `175,720` coeffs, `0` negative), plus $n = 3..9, 12, 16$ checks and `5/5` must-fail controls.

## 2. Canonical Sylvester–Hadamard Continuous Bar-Preserving Fold $F_8 \to E_8$ + Exact $E_6 / E_7 / D_n^+$ Uniqueness (`AUDIT-017`, `N2`)

- **Proof:** [`E8_SPINOR_FOLD_PROOF.md`](./E8_SPINOR_FOLD_PROOF.md) (Lemmas 1–5 & Theorem 2).
- **Verifier & Certificate:** [`verify_e8_spinor_fold.py`](./verify_e8_spinor_fold.py) $\to$ [`e8_spinor_fold_certificate.json`](./e8_spinor_fold_certificate.json) (bit-identical SHA-256 `e4f14b928a3b08897e562ed56e334c6471ccf12e6bd55f52b6273f2577b0ef31` across all 4 isolated R4 modes).
- **Key Theorems:**
  - **Canonical Sylvester Hadamard $H_8 = H_2^{\otimes 3}$:** Exact zero-trigonometric $\mathbb{Q}(\sqrt{2})(s)$ curve $L_8(s)$ preserving all 8 bars $\|r_i(s) - 2e_i\|^2 \equiv 2$ identically from $F_8$ ($L_8(0) = (2+\sqrt{2})I_8$) to exact optimal 8D $E_8$ ($L_8(1) = 2 I_8 + \frac{1}{2} H_8$, $L_8(1)^{-1} = I_8 - \frac{1}{4} H_8$, $\det = 16$, optimal density $\pi^4/384 = 2\times \Delta(D_8)$, quotient $E_8/\mathbb{Z}^8 L_8(1) \cong (\mathbb{Z}/2\mathbb{Z})^4$, sublattice $\mathbb{Z}^8 L_8(1) \cong \sqrt{2}\,E_8$, exact kissing number **240** on all 16 motif balls).
  - **Exact $2\times 2$ block-diagonalization & strict volume monotonicity:** Pure-integer orthogonal $Q_4$ ($Q_4^T Q_4 = 4I_4$) block-diagonalizes $P^T L_8(s) P$ identically into $3\times M_{+2}(s) \oplus M_{-2}(s)$; $-\frac{d}{ds}\det L_8(s) / s$ (deg 10) has all 11 Bernstein coeffs strictly $> 0$.
  - **Exact Loewner monotonicity & zero-subdivision global $\mathbb{Z}^8$ certificate:** Exact $2\times 2$ Sylvester/Schur check proves $G(s) = L_8(s) L_8(s)^T \succeq G(1) = 6I_8 + 2H_8$ on $[0, 1]$. Pure-rational $LDL^T$ (`pivots = [8, 7/2, 20/7, 5, 4/5, 7/5, 8/7, 1/2]`) + `math.isqrt` Fincke–Pohst over $m^T G(1) m \le \lfloor(2+\sqrt{2})^2\rfloor = 11$ yields `120` ($\text{norm}^2=4$) and `1,080` ($\text{norm}^2=8$) canonical $\pm m$ pairs; exhaustive degree-8 $\mathbb{Z}[\sqrt{2}]$ Bernstein certification over all **`32,640` pairs (`293,760` coeffs)** finds **`0` negative, `0` subdivisions, `976` exact zeros** (`72` persistent bar zeros + `904` $s=1$ contacts, `0` interior zeros on all `32,632` non-bar pairs $\implies$ strict interior separation on $[0, 1)$), plus `5/5` $B_8$ gauge trials and `7/7` must-fail controls.
  - **Exact $E_6 / E_7 / D_n^+$ Uniqueness Theorem:** Determinant-square obstruction rules out $F_6 \not\to E_6$ ($4/3 \notin \mathbb{Z}^2$) and $F_7 \not\to E_7$ ($2 \notin \mathbb{Z}^2$); spinor-norm bar trapping rules out $F_n \not\to D_n^+$ for all even $n \ne 8$.

## 3. Complete Optimal-Lattice Classification ($n = 3, 4, 5, 6, 7, 8, 16, 24$): Direct-Sum Obstruction + 16D Barnes–Wall $\Lambda_{16}$ & 24D Leech $\Lambda_{24}/\sqrt{2}$ (`AUDIT-020R`, `N2`)

- **Proof:** [`BW16_LEECH24_UNJAMMING_PROOF.md`](./BW16_LEECH24_UNJAMMING_PROOF.md) (Theorems 1–3 & Lemmas A–B).
- **Verifier & Certificate:** [`verify_bw16_leech24.py`](./verify_bw16_leech24.py) $\to$ [`bw16_leech24_certificate.json`](./bw16_leech24_certificate.json) (bit-identical SHA-256 `5dbdeb6f470483242705bb632f6dcaa41ec58bb4f8ee7c3d28f2a6b5ab25cadc` across all 4 isolated R4 modes, `9/9` controls rejected).
- **Key Theorems:**
  - **Irreducible Root-Graph Obstruction ($E_8 \oplus E_8$ & ALL nontrivial direct sums $\Gamma_1 \oplus \Gamma_2$):** $D_{16}$'s 480-root graph is connected (`[480]`, rank 16) whereas $E_8 \oplus E_8$'s splits into two orthogonal rank-8 components (`[240, 240]`), proving $F_{16} \not\to E_8 \oplus E_8$ even without bars and $\forall R \in O(16)$.
  - **16D Barnes–Wall $\Lambda_{16}$ (kissing 4,320; $32\times\Delta(D_{16}) = 16\times\Delta(E_8^2)$):** Exact $\mathrm{RM}(1,4)$ $(\mathbb{Z}/2\mathbb{Z})^5$ quotient + explicit positive-orientation bar-preserving endpoint ($\det = +2$) + closed-form exact rational strain $X_{16} = \Pi_{\mathrm{bar}}(I_{16})$ ($X_{16,\mathrm{sym}} \succeq \frac{1}{5} I_{16}$, `0` negative contacts across all `137,280` inter-cell contacts, `32` exact bar zeros, $\delta_{16} > 1/7$) + certified global-$\mathbb{Z}^{16}$ expanding bar-preserving Cayley arc $\tau \in [0, 1/4000]$ (accounting for integral-norm shells $\|w\|^2 \ge 3$ via $\|A(\tau)\|_2 \le 1/8$).
  - **24D Optimal Leech $\Lambda_{24}/\sqrt{2}$ (Cohn–Kumar–Miller–Radchenko–Viazovska 2017: kissing 196,560; $8,192\times\Delta(D_{24})$):** Golay self-duality $(\mathbb{Z}/2\mathbb{Z})^4 \times \mathbb{Z}/3\mathbb{Z}$ quotient (`V_GOLAY_4` $\subset \mathcal{G}_{24}^\perp = \mathcal{G}_{24}$) + explicit positive-orientation bar-preserving endpoint ($\det = +3/256$) + exact rational projected strain $X_{24}$ ($X_{24,\mathrm{sym}} \succeq \frac{1}{5} I_{24}$, `0` negative contacts across ALL `9,432,672` inter-cell contacts, `48` exact bar zeros, $\delta_{24} \ge 2/25$) + certified global-$\mathbb{Z}^{24}$ expanding bar-preserving Cayley arc $\tau \in [0, 10^{-10}]$. *(Honest scope: global extension of the finite expanding arcs all the way to $F_{16}/F_{24}$ remains open.)*

## Reproducibility

```bash
cd results/expedition-dn-e8-folds
python3 -I -B verify_all_n_dn_fold.py > /dev/null
python3 -I -O -B -W error verify_all_n_dn_fold.py > /dev/null
python3 -I -B verify_e8_spinor_fold.py
python3 -I -O -B -W error verify_e8_spinor_fold.py
python3 -I -B verify_bw16_leech24.py
python3 -I -O -B -W error verify_bw16_leech24.py
```
