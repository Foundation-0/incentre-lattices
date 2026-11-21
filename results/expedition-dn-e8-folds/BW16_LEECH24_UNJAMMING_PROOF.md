# Exact Feasibility-and-Obstruction Classification of Optimal Sphere-Packing Lattices from the Bridged Cross-Polytope Framework $F_n$ ($n = 3, 4, 5, 6, 7, 8, 16, 24$)

**Status:** Complete zero-float proof + independent pure Python 3 stdlib verifier [`verify_bw16_leech24.py`](./verify_bw16_leech24.py) and bit-deterministic certificate [`bw16_leech24_certificate.json`](./bw16_leech24_certificate.json).

---

## 1. Master Classification Table

For the $2n$-ball unit cross-polytope motif $P_n = \{\pm e_0, \dots, \pm e_{n-1}\}$ (radius $\sqrt{2}/2$, adjacent contact distance $\sqrt{2}$) and bridged cubic framework $F_n = P_n + (2+\sqrt{2})\mathbb{Z}^n$, we classify every canonical optimal/dense lattice $\Gamma$ normalized to minimum squared norm $2$:

| Dim $n$ | Target Lattice $\Gamma$ | Kissing $\|\Phi(\Gamma)\|$ | Density vs. $D_n$ | Periodic Tiling $P_n + \Lambda = \Gamma + v_0$? | Bar-Preserving Endpoint $L_n(1)=2I_n+C_n$ ($\det>0$)? | Continuous Bar-Preserving Motion Status |
|---|---|---|---|---|---|---|
| $\forall n \ge 3$ | $D_n$ (opt. $n=3,4,5$) | $2n(n-1)$ | $1\times$ | **YES** (`AUDIT-009`) | **YES** ($\det=4n$) | **Global fold $F_n \to D_n$ proved** (`AUDIT-009`) |
| $6$ | $E_6$ (opt. 6D) | $72$ | $2/\sqrt{3}\times$ | **NO** ($4/3 \notin \mathbb{Z}^2$) | **NO** | **Obstructed** (even without bars) |
| $7$ | $E_7$ (opt. 7D) | $126$ | $\sqrt{2}\times$ | **NO** ($2 \notin \mathbb{Z}^2$) | **NO** | **Obstructed** (even without bars) |
| $8$ | **$E_8$ (opt. 8D)** | **$240$** | **$2\times$** | **YES** (`AUDIT-017`) | **YES** ($H_8 = H_2^{\otimes 3}$, $\det=16$) | **Global fold $F_8 \to E_8$ proved** (`AUDIT-017`) |
| even $n\ne 8$ | $D_n^+$ | — | $2\times$ ($n\ge 8$) | Overlap ($n=4,6$) / YES ($n\ge 10$) | **NO** ($n/4 \ne 2$) | **Obstructed** (bar generators trapped in $D_n$) |
| $16$ | $E_8 \oplus E_8$ | $480$ | $2\times$ | **NO** (root graph disconnected!) | **NO** | **Obstructed** (even without bars, $\forall R \in O(16)$) |
| $\forall n\ge 3$ | Any $\Gamma_1 \oplus \Gamma_2$ | $\|\Phi_1\|+\|\Phi_2\|$ | — | **NO** (Thm 1b) | **NO** | **Obstructed** (even without bars, $\forall R \in O(n)$) |
| **$16$** | **Barnes–Wall $\Lambda_{16}$** | **$4{,}320$** | **$32\times$ ($16\times E_8^2$)** | **YES** (Thm 2, $(\mathbb{Z}/2)^5$) | **YES** ($\det = +2$, all $\mathrm{RM}(1,4)$) | **Strictly unjammed + global-$\mathbb{Z}^{16}$ expanding bar-preserving arc $\tau \in [0, \frac{1}{4000}]$ proved**; global connection to $F_{16}$ open |
| **$24$** | **Leech $\Lambda_{24}/\sqrt{2}$ (opt. 24D)** | **$196{,}560$** | **$8{,}192\times$** | **YES** (Thm 3, $(\mathbb{Z}/2)^4\times\mathbb{Z}/3$) | **YES** ($\det = +\frac{3}{256}$, Golay) | **Strictly unjammed + global-$\mathbb{Z}^{24}$ expanding bar-preserving arc $\tau \in [0, 10^{-10}]$ proved**; global connection to $F_{24}$ open |

---

## 2. Theorem 1 (Exact Obstructions: $E_6, E_7$, Direct Sums $E_8 \oplus E_8$, and $D_n^+$)

1. **Determinant-Square Obstruction ($E_6, E_7$):** Any periodic packing $P_n + \Lambda = \Gamma + v_0$ forces $D_n = \langle P_n - P_n \rangle \subseteq \Gamma$ as a finite-index sublattice, hence $\det_{\mathrm{Gram}}(D_n)/\det_{\mathrm{Gram}}(\Gamma) = [\Gamma : D_n]^2 \in \mathbb{Z}^2$. At minimum squared norm $2$, $\det_{\mathrm{Gram}}(D_n) = 4$ while $\det_{\mathrm{Gram}}(E_6) = 3$ and $\det_{\mathrm{Gram}}(E_7) = 2$; since $4/3 \notin \mathbb{Z}^2$ and $4/2 = 2 \notin \mathbb{Z}^2$, neither $E_6$ nor $E_7$ is reachable even without bars.
2. **Irreducible Root-Graph Obstruction ($\Gamma_1 \oplus \Gamma_2$ and $E_8 \oplus E_8$):** Suppose $\Gamma = \Gamma_1 \oplus \Gamma_2$ is any nontrivial orthogonal direct sum ($\dim \Gamma_1, \dim \Gamma_2 \ge 1$) normalized to minimum squared norm $2$. Every minimal vector $w \in \Phi(\Gamma)$ satisfies $\|w\|^2 = \|w_1\|^2 + \|w_2\|^2 = 2$, so either $w_2 = 0$ or $w_1 = 0$; hence $\Phi(\Gamma) = \Phi(\Gamma_1) \sqcup \Phi(\Gamma_2)$ with $\Phi(\Gamma_1) \perp \Phi(\Gamma_2)$. Define the root graph $G(\Phi)$ with edges $(u, v)$ whenever $|u \cdot v| = 1$. Since $D_n$ ($n \ge 3$) is an irreducible root system, $G(\Phi(D_n))$ is **connected** and spans $\mathbb{R}^n$ (for $n=16$, a single component of size $480$); by contrast, $G(\Phi(\Gamma_1 \oplus \Gamma_2))$ is disconnected with every connected component lying inside a proper subspace $\mathbb{R}^{\dim\Gamma_k}$ (for $E_8 \oplus E_8$, two disjoint components of size $240$ and rank $8$). Therefore no orthogonal rotation $R \in O(n)$ can embed $D_n$ into $\Gamma_1 \oplus \Gamma_2$, proving **$F_{16} \not\to E_8 \oplus E_8$** (and $F_n \not\to \Gamma_1 \oplus \Gamma_2$) even without bar constraints. *(Note: the block-diagonal fold $F_8 \oplus F_8 \to E_8 \oplus E_8$ uses the $256$-ball product motif $P_8 \times P_8$, not the $32$-ball cross-polytope $P_{16}$.)*
3. **Spinor-Norm Bar-Trapping Obstruction ($D_n^+$, even $n \ne 8$):** Proved in [`E8_SPINOR_FOLD_PROOF.md`](./E8_SPINOR_FOLD_PROOF.md) §3 and cross-checked in [`verify_bw16_leech24.py`](./verify_bw16_leech24.py). $\square$

---

## 3. Universal Quantitative Unjamming & Global-$\mathbb{Z}^n$ Expanding Cayley Arc Lemma

Let $\Gamma \subset \mathbb{R}^n$ be any **integral-norm lattice** ($\|w\|^2 \in \mathbb{Z}_{\ge 2}$ for all $w \in \Gamma \setminus \{0\}$, so every non-contact shell satisfies $\|w\|^2 \ge 3$) containing $D_n$, and let $L = 2 I_n + C$ ($\|c_i\|^2 = 2$, $\det L > 0$) satisfy $P_n + \mathbb{Z}^n L = \Gamma + e_0$. *(Note: Standard $\Lambda_{16}$ and $\Lambda_{24}$ are even at minimum squared norm $4$, with shells $4, 6, 8, \dots$; rescaling by $1/\sqrt{2}$ to minimum squared norm $2$ gives integral-norm lattices $\Gamma_{16}, \Gamma_{24}$ with shells $2, 3, 4, \dots$, proved explicitly in §§4–5.)*

### Lemma A (Universal Bar-Free Isotropic Unfolding)
For every $t \ge 0$, every $m \in \mathbb{Z}^n \setminus \{0\}$, and every $p, q \in P_n$ ($v = q - p$, $\|v\|^2 \le 4$), letting $x = m L \in \Gamma \setminus \{0\}$ and $w = x + v \in \Gamma \setminus \{0\}$ (so $\|x\|^2 \ge 2, \|w\|^2 \ge 2$):
$$2 w \cdot x = \|w\|^2 + \|x\|^2 - \|v\|^2 \ge 2 + 2 - 4 = 0 \implies \|(1 + t) m L + v\|^2 = \|w\|^2 + 2 t (w \cdot x) + t^2 \|x\|^2 \ge 2 + 2 t^2.$$
Thus without bar constraints, every periodic cross-polytope packing ($D_n, E_8, \Lambda_{16}, \Lambda_{24}/\sqrt{2}$) unfolds globally and overlap-free via pure scaling. $\square$

### Lemma B (Zero-Float Global-$\mathbb{Z}^n$ Bar-Preserving Cayley Arc)
Suppose a rational strain matrix $X \in \mathbb{Q}^{n\times n}$ ($V := L X$) satisfies:
1. **(Exact bar tangency)** $r_i X c_i^T = v_i \cdot c_i = 0$ for all $i = 0, \dots, n-1$;
2. **(Strict symmetric positivity)** $X_{\mathrm{sym}} := \frac{1}{2}(X + X^T) \succeq \lambda_0 I_n$ for a rational $\lambda_0 > 0$ (certified via exact rational $LDL^T$);
3. **(Exact non-bar contact positivity)** For every inter-cell kissing contact $(p, q, m \ne 0)$ with $w = m L + q - p \in \Phi(\Gamma)$ ($\|w\|^2 = 2$), the first-order separation velocity $w \cdot (m V) = (w - (q-p))^T X w$ is $= 0$ on the $2n$ directed bar contacts and $\ge \delta > 0$ on all other inter-cell contacts.

Define the **exact rational bar-preserving Cayley great-circle curve** $L(\tau) = 2 I_n + u(\tau)$ ($\tau \ge 0$) by
$$u_i(\tau) := \frac{\bigl(4 - \frac{1}{2}\|v_i\|^2 \tau^2\bigr) c_i + 4 \tau v_i}{4 + \frac{1}{2}\|v_i\|^2 \tau^2}.$$
Then $\|u_i(\tau)\|^2 \equiv 2$ identically (exact bar preservation for all $\tau$, since $2(4 - \frac{1}{2}\|v_i\|^2\tau^2)^2 + 16\tau^2\|v_i\|^2 = 2(4 + \frac{1}{2}\|v_i\|^2\tau^2)^2$), $u_i(0) = c_i$, $u_i'(0) = v_i$, and $L(\tau) = L + \tau V + R(\tau)$ where $E(\tau) := L^{-1} R(\tau)$ and $A(\tau) := \tau X + E(\tau)$. Every row remainder $R_i(\tau) = -\frac{\frac{1}{2}\|v_i\|^2\tau^2(2 c_i + \tau v_i)}{4 + \frac{1}{2}\|v_i\|^2\tau^2}$ obeys
$$\|R_i(\tau)\|^2 = \frac{\frac{1}{4}\|v_i\|^4 \tau^4}{\bigl(4 + \frac{1}{2}\|v_i\|^2 \tau^2\bigr)^2}(8 + \|v_i\|^2 \tau^2) \le \frac{1}{8}\|v_i\|^4 \tau^4 \qquad \forall \tau \ge 0$$
(since $2(8 + z) \le (4 + z/2)^2 = 16 + 4z + z^2/4$ for all $z \ge 0$). Setting the pure-rational constants $K_X := \|X\|_F^2$, $K_R := \frac{1}{8}\sum_{i=0}^{n-1} \|v_i\|^4$, and $K_E := \|L^{-1}\|_F^2 K_R$, we have $\|E(\tau)\|_2^2 \le \|E(\tau)\|_F^2 \le K_E \tau^4$ for all $\tau \ge 0$.

If a rational $\tau_* > 0$ satisfies the **three pure-rational inequalities**:
$$(1)\; K_E \tau_*^2 \le \frac{\lambda_0^2}{4}, \qquad (2)\; K_X \tau_*^2 \le \frac{1}{144} \text{ and } \frac{\lambda_0 \tau_*}{2} \le \frac{1}{24}, \qquad (3)\; 25 K_E \tau_*^2 < \delta^2,$$
then for **every $\tau \in (0, \tau_*]$** and **every pair $(m, p, q) \ne (0, p, p)$ across all infinitely many vectors $m \in \mathbb{Z}^n$**:
- **Bar pairs ($m = e_i, v = -2 e_i$):** $\|r_i(\tau) - 2 e_i\|^2 = \|u_i(\tau)\|^2 \equiv 2$.
- **All infinitely many non-contact pairs ($\|w\|^2 \ge 3 \implies \|w\| \ge \sqrt{3}$):** By (1), for every unit $y \in \mathbb{R}^n$, $\|(I_n + A(\tau)) y\| \ge y^T(I_n + \tau X_{\mathrm{sym}} + E(\tau)) y \ge 1 + \lambda_0 \tau - \|E(\tau)\|_2 \ge 1 + \frac{\lambda_0}{2}\tau > 1$. By (2), $\|A(\tau)\|_2 \le \sqrt{K_X}\tau + \frac{\lambda_0}{2}\tau \le \frac{1}{12} + \frac{1}{24} = \frac{1}{8}$. Hence
  $$\|m L(\tau) + v\| = \|w(I_n + A(\tau)) - v A(\tau)\| \ge \sigma_{\min}(I_n + A(\tau))\|w\| - \|A(\tau)\|_2\|v\| > 1\cdot \sqrt{3} - \frac{1}{8}\cdot 2 = \sqrt{3} - \frac{1}{4} > \sqrt{2},$$
  since $\sqrt{3} - \sqrt{2} = \frac{1}{\sqrt{3}+\sqrt{2}} > \frac{1}{2+2} = \frac{1}{4}$ (zero-float proof!).
- **All non-bar contact pairs ($\|w\|^2 = 2$, $\|x\| = \|w - v\| \le 2 + \sqrt{2}$):** Since $2\sqrt{2}(2+\sqrt{2}) = 4 + 4\sqrt{2} < 10$, condition (3) gives $(4+4\sqrt{2})\|E(\tau)\|_2 < 10\sqrt{K_E}\tau^2 < 2\delta\tau$, so
  $$\|m L(\tau) + v\|^2 = 2 + 2\tau w\cdot(x X) + 2 w\cdot(x E(\tau)) + \|x A(\tau)\|^2 \ge 2 + \bigl(2\delta\tau - 10\sqrt{K_E}\tau^2\bigr) > 2.$$
- **Strict volume expansion:** Since $I_n + A(\tau)$ has strictly positive symmetric part on $(0, \tau_*]$, $\det L(\tau) = \det L \cdot \det(I_n + A(\tau)) > \det L$ strictly (with $\left.\frac{d}{d\tau}\right|_{\tau=0} \ln \det L(\tau) = \operatorname{tr}(X) \ge n \lambda_0 > 0$). $\square$

---

## 4. Theorem 2 (16D Barnes–Wall Lattice $\Lambda_{16}$)

Normalized to minimum squared norm $2$, the Barnes–Wall lattice is
$$\Gamma_{16} := \frac{1}{2}\left\{ x \in \mathbb{Z}^{16} : x \bmod 2 \in \mathrm{RM}(1, 4),\; \sum_{i=0}^{15} x_i \equiv 0 \pmod 4 \right\},$$
with covolume $2^{-4} = \frac{1}{16}$, $D_{16} \subset \Gamma_{16}$ of index $32$, integral squared norms ($\|x\|^2 \equiv \operatorname{wt}(x\bmod 2) \in \{0, 8, 16\} \equiv 0 \pmod 4 \implies \|x/2\|^2 \in \mathbb{Z}$), minimum squared norm $2$ (zero vectors of norm$^2 < 2$ since every nonzero codeword of $\mathrm{RM}(1,4)$ has weight $\ge 8$), and exact kissing number $480 + 30 \cdot 128 = \mathbf{4{,}320}$ (density $\Delta(\Lambda_{16}) = \frac{\pi^8}{16\cdot 8!} = 32\times\Delta(D_{16}) = 16\times\Delta(E_8^2)$).

- **Exact $(\mathbb{Z}/2\mathbb{Z})^5$ Quotient:** On $D_{16}$, $\psi_{D_{16}}(x) = \bigl(\bigoplus_{i:\, x_i \text{ odd}} i,\; (\tfrac{1}{2}\sum x_i)\bmod 2\bigr) \in \mathbb{F}_2^5$ satisfies $\psi_{D_{16}}(2g) = 0$ for every $g \in \Gamma_{16}$ (every weight-8 hyperplane of $\mathbb{F}_2^4$ sums to $0$ in $\mathbb{F}_2^4$ and has sum $\equiv 0 \pmod 4$). Extending by $u = (0, 16, 16, 16, 16)$ on the 5 generators of $\mathrm{RM}(1, 4)$ yields an exact surjective homomorphism $\Psi_{16} : \Gamma_{16} \to (\mathbb{Z}/2\mathbb{Z})^5$ mapping $P_{16} - e_0$ bijectively onto all $32$ cosets (verified on all $21 \times 21$ generator pairs).
- **Exact Positive-Orientation Bar-Preserving Endpoint:** `C2_16_WITNESS` in [`verify_bw16_leech24.py`](./verify_bw16_leech24.py) has all $16$ rows equal to weight-8 $\mathrm{RM}(1,4)$ minimal vectors with $c_{ii} = -1/2$, every $r_i = 2 e_i + \frac{1}{2}(C_2^{(16)})_i \in \ker\Psi_{16}$, and exact `Fraction` determinant $\det L_{16}(1) = +2 = 32 \cdot \frac{1}{16}$, proving $\mathbb{Z}^{16} L_{16}(1) = \ker\Psi_{16}$ and $P_{16} + \mathbb{Z}^{16} L_{16}(1) = \Gamma_{16} + e_0$.
- **Exact Closed-Form Unjamming & Global-$\mathbb{Z}^{16}$ Expanding Arc:** The exact orthogonal projection $X_{16} = \Pi_{\mathrm{bar}}(I_{16}) \in \mathbb{Q}^{16\times 16}$ (computed directly from $I_{16}$ with zero hand-tuning!) satisfies exact bar tangency $r_i X_{16} c_i^T \equiv 0$, exact $LDL^T$ positivity $X_{16,\mathrm{sym}} \succeq \frac{1}{5} I_{16}$, and across **all $137{,}280$ inter-cell contacts** ($32 \times 4{,}320 - 32 \times 30$): **$0$ negative**, **$32$ exact bar zeros**, and minimum non-bar velocity $\delta_{16} > \frac{1}{7}$ ($\approx 0.1658$). Lemma B certifies global-$\mathbb{Z}^{16}$ non-overlap on $\tau \in [0, \frac{1}{4000}]$. $\square$

---

## 5. Theorem 3 (24D Optimal Leech Lattice $\Lambda_{24}/\sqrt{2}$)

Using the self-dual $[24, 12, 8]$ extended binary Golay code $\mathcal{G}_{24}$ (weight distribution $\{0:1, 8:759, 12:2576, 16:759, 24:1\}$), the Leech lattice rescaled to minimum squared norm $2$ around a $D_{24}$ deep hole is
$$\Gamma_{24} = \frac{\Lambda_{24}}{\sqrt{2}} = \frac{1}{4}\left\{ u \in \mathbb{Z}^{24} : \begin{array}{l}\text{(a) all } u_i \text{ even, } \frac{u}{2}\bmod 2 \in \mathcal{G}_{24},\; \sum u_i \equiv 0 \pmod 8,\text{ or}\\\text{(b) all } u_i \text{ odd, } \frac{u-1}{2}\bmod 2 \in \mathcal{G}_{24},\; \sum u_i \equiv 4 \pmod 8\end{array}\right\},$$
with covolume $2^{-12}$, $D_{24} \subset \Gamma_{24}$ of index $2^{13} = 8{,}192$, integral squared norms (for type (a), $u = 2v$ has $\|v\|^2 \equiv \operatorname{wt}(v\bmod 2) \equiv 0 \pmod 4 \implies \|u\|^2 \equiv 0 \pmod{16}$; for type (b), writing $u_i = 1 + 2c_i + 4y_i$ with $c \in \mathcal{G}_{24}$ gives $u_i^2 \equiv 1 + 8c_i + 8y_i \pmod{16}$, while $\sum u_i = 24 + 2\operatorname{wt}(c) + 4\sum y_i \equiv 4 \pmod 8$ forces $\sum y_i \equiv 1 \pmod 2$, whence $\sum u_i^2 \equiv 24 + 8\operatorname{wt}(c) + 8 \equiv 0 \pmod{16} \implies \|u/4\|^2 \in \mathbb{Z}$), zero vectors of norm$^2 < 2$ (for type (b), all-$(\pm 1)$ gives $\sum u_i = 24 - 2\operatorname{wt}(c) \equiv 0 \not\equiv 4 \pmod 8$, forcing some $|u_i| \ge 3$ and $\|u\|^2 \ge 9 + 23 = 32$), and exact kissing number $1{,}104 + 97{,}152 + 98{,}304 = \mathbf{196{,}560}$ (optimal 24D density $\Delta(\Lambda_{24}) = \pi^{12}/12! = 8{,}192\times\Delta(D_{24})$, Cohn–Kumar–Miller–Radchenko–Viazovska 2017).

- **Golay Self-Duality Quotient Lemma:** Any homomorphism $\psi_2 : \Gamma_{24} \to \mathbb{F}_2^4$ must vanish on $2 x_O = \mathbf{1}_O$ for all $759$ octads $O \in \mathcal{G}_{24}$, forcing every coordinate bit-slice of $m_i = \psi_2(e_i - e_0)$ into $\mathcal{G}_{24}^\perp = \mathcal{G}_{24}$. The explicit 4-tuple `V_GOLAY_4` $\subset \mathcal{G}_{24}$ (whose first three codewords span a constant-weight-12 3D subspace partitioning $\{0,\dots,23\}$ into eight 3-point fibers, split by the fourth codeword `v3` and $\mathbb{Z}/3\mathbb{Z}$ assignment `a_z3`) defines an exact surjective homomorphism $\Psi_{24} : \Gamma_{24} \to (\mathbb{Z}/2\mathbb{Z})^4 \times \mathbb{Z}/3\mathbb{Z}$ (verified on all $37 \times 37$ generator pairs) mapping $P_{24} - e_0$ bijectively onto all $48$ cosets!
- **Exact Positive-Orientation Bar-Preserving Endpoint:** `U24_WITNESS` in [`verify_bw16_leech24.py`](./verify_bw16_leech24.py) has all $24$ rows equal to exact minimal Leech vectors ($\|u_i\|^2 = 32 \iff \|c_i\|^2 = 2$), every $r_i = 2 e_i + \frac{1}{4} u_i \in \ker\Psi_{24}$, and exact `Fraction` determinant $\det L_{24}(1) = +\frac{3}{256} = 48 \cdot 2^{-12}$, proving $\mathbb{Z}^{24} L_{24}(1) = \ker\Psi_{24}$ and $P_{24} + \mathbb{Z}^{24} L_{24}(1) = (\Lambda_{24}/\sqrt{2}) + e_0$.
- **Exact Zero-Float Unjamming Over All $9{,}432{,}672$ Inter-Cell Contacts & Global-$\mathbb{Z}^{24}$ Arc:** Exact rational projection of `X0_24_INT_256 / 256` onto $\{X : r_i X c_i^T = 0\}$ yields $X_{24} \in \mathbb{Q}^{24\times 24}$ with exact zero bar residual, exact $LDL^T$ positivity $X_{24,\mathrm{sym}} \succeq \frac{1}{5} I_{24}$, and across **all $9{,}432{,}672$ inter-cell contacts** ($48 \times 196{,}560 - 48 \times 46$, checked in pure integer arithmetic via Gray-code/incremental updates + exact sufficient bound $u^T M u - 8\|M u\|_\infty \ge \lceil \frac{16 D}{10}\rceil$): **$0$ negative**, **$48$ exact bar zeros**, and exact minimum non-bar velocity $\delta_{24} > \frac{2}{25}$ ($\approx 0.08193$). Lemma B certifies global-$\mathbb{Z}^{24}$ non-overlap on $\tau \in [0, 10^{-10}]$. $\square$

---

## 6. Honest Scope Boundary & Reproducibility

- **What is proved with zero gaps:** Complete obstructions ($E_6, E_7, E_8\oplus E_8$, all $\Gamma_1\oplus\Gamma_2$, $D_n^+$ for $n\ne 8$); exact positive-orientation bar-preserving realizations of $\Lambda_{16}$ and $\Lambda_{24}/\sqrt{2}$; exact zero-float strict bar-preserving unjamming over all $137{,}280$ / $9{,}432{,}672$ contacts; and exact global-$\mathbb{Z}^n$ non-overlap on certified expanding bar-preserving Cayley arcs $\tau \in [0, \tau_*]$ (plus global bar-free isotropic unfolding Lemma A).
- **What remains open:** Extending the local expanding bar-preserving arcs $\tau \in [0, \tau_*]$ globally all the way to $F_{16}$ / $F_{24}$ without intermediate collisions (note that the $R^2=11$ ellipsoid contains $\sim 4\times 10^7$ points in 16D and $\sim 1.1\times 10^{13}$ points in 24D, requiring a dimension-free analytical invariant rather than ellipsoid enumeration).
- **Repro (all 4 isolated R4 modes):**
  ```bash
  cd results/expedition-dn-e8-folds
  python3 -I -B verify_bw16_leech24.py
  python3 -I -O -B verify_bw16_leech24.py
  python3 -I -B -W error verify_bw16_leech24.py
  python3 -I -O -B -W error verify_bw16_leech24.py
  ```
