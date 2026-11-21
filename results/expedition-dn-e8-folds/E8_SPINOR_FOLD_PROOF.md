# Canonical Sylvester-Hadamard Continuous Bar-Preserving Fold $F_8 \to E_8$ and the $E_6 / E_7 / D_n^+$ Uniqueness Theorem

**Status:** Complete, zero-float, zero-subdivision mathematical proof + independent stdlib-only verifier [`verify_e8_spinor_fold.py`](./verify_e8_spinor_fold.py) and bit-deterministic certificate [`e8_spinor_fold_certificate.json`](./e8_spinor_fold_certificate.json).

---

## 1. Main Results

Let $P_8 = \{\pm e_0, \dots, \pm e_7\} \subset \mathbb{R}^8$ be the 16-ball unit cross-polytope motif (sphere radius $r = \sqrt{2}/2$, adjacent contact distance $\sqrt{2}$). Starting from the cubic cross-polytope packing $F_8 = P_8 + (2+\sqrt{2})\mathbb{Z}^8$, a **bar-preserving continuous lattice fold** is a continuous path $L_8 : [0, 1] \to \mathrm{GL}_8(\mathbb{R})$ with $L_8(0) = (2+\sqrt{2}) I_8$ such that:
1. **(Bar preservation)** Every inter-cell bar connecting $+e_i$ in cell $0$ to $-e_i + r_i(s)$ in cell $e_i$ has identically constant length $\sqrt{2}$:
   $$\|r_i(s) - 2 e_i\|^2 \equiv 2 \qquad \forall i \in \{0,\dots,7\},\; \forall s \in [0, 1].$$
2. **(Global non-overlap over all of $\mathbb{Z}^8$)** Every pair of distinct spheres in $P_8 + \mathbb{Z}^8 L_8(s)$ satisfies
   $$\|m L_8(s) + (p - q)\|^2 \ge 2 \qquad \forall (m, p, q) \in \mathbb{Z}^8 \times P_8 \times P_8 \setminus \{(0, p, p)\},\; \forall s \in [0, 1].$$

### Theorem 1 (Canonical Sylvester-Hadamard Continuous Fold $F_8 \to E_8$)

Let $H_8 = H_2^{\otimes 3} \in \{\pm 1\}^{8\times 8}$ be the canonical 8D Sylvester Hadamard matrix, $(H_8)_{ij} = (-1)^{\operatorname{popcount}(i \mathbin{\&} j)}$ for $i, j \in \{0,\dots,7\}$, with diagonal index partition
$$I_+ = \{0, 3, 5, 6\} \quad ((H_8)_{ii} = +1), \qquad I_- = \{1, 2, 4, 7\} \quad ((H_8)_{ii} = -1).$$
Define the exact zero-trigonometric $\mathbb{Q}(\sqrt{2})(s)$ curve $L_8(s)$ ($s \in [0, 1]$) with denominators
$$D_+(s) = (1 + 2\sqrt{2}) + (-1 + 2\sqrt{2}) s^2, \qquad D_-(s) = (-1 + 2\sqrt{2}) + (1 + 2\sqrt{2}) s^2,$$
scalar entries
$$a_+(s) = \frac{(6 + 5\sqrt{2}) + (-6 + 5\sqrt{2}) s^2}{D_+(s)}, \qquad b_+(s) = \frac{2\sqrt{2}\, s}{D_+(s)},$$
$$a_-(s) = \frac{(2 + 3\sqrt{2}) + (-2 + 3\sqrt{2}) s^2}{D_-(s)}, \qquad b_-(s) = \frac{2\sqrt{2}\, s}{D_-(s)},$$
and rows
$$(L_8(s))_{ij} = \begin{cases} a_+(s) & i = j \in I_+, \\ b_+(s) (H_8)_{ij} & i \in I_+,\, j \ne i, \\ a_-(s) & i = j \in I_-, \\ b_-(s) (H_8)_{ij} & i \in I_-,\, j \ne i. \end{cases}$$
Then:
1. **(Exact bar-preserving fold)** $L_8(0) = (2+\sqrt{2}) I_8$, $L_8(1) = 2 I_8 + \frac{1}{2} H_8$, and $\|r_i(s) - 2 e_i\|^2 \equiv 2$ identically for all $i$ and all $s \in [0, 1]$.
2. **(Exact optimal $E_8$ endpoint — Viazovska 2016)** $P_8 + \mathbb{Z}^8 L_8(1) = E_8 + e_0$, with $\det L_8(1) = 16$, quotient $E_8 / \mathbb{Z}^8 L_8(1) \cong (\mathbb{Z}/2\mathbb{Z})^4$, sublattice $\mathbb{Z}^8 L_8(1) \cong \sqrt{2}\,E_8$, optimal 8D sphere-packing density $\Delta(E_8) = \pi^4/384$ (exactly $2\times$ the $D_8$ density $\pi^4/768$), and exact kissing number $240$ at all 16 motif balls.
3. **(Strict volume monotonicity)** $\det L_8(s)$ decreases strictly monotonically on $[0, 1]$ from $(2+\sqrt{2})^8 = 9232 + 6528\sqrt{2} \approx 18464$ down to $16$ (a $\approx 1154\times$ density gain).
4. **(Exact Loewner monotonicity)** $G(s) := L_8(s) L_8(s)^T \succeq G(1) = 6 I_8 + 2 H_8$ for every $s \in [0, 1]$.
5. **(Global non-overlap and strict interior separation)** $\|m L_8(s) + (p - q)\|^2 \ge 2$ for **all** $m \in \mathbb{Z}^8$, $p, q \in P_8$ ($(m,p,q) \ne (0,p,p)$), and **all** $s \in [0, 1]$, certified with **zero floats** and **zero subdivisions** by $32,640$ degree-8 $\mathbb{Z}[\sqrt{2}]$ Bernstein polynomials ($293,760$ exact coefficients: $0$ negative, $976$ exact zeros). Moreover, every non-bar inter-cell pair is **strictly separated** ($\|m L_8(s) + (p - q)\|^2 > 2$) on the entire interval $s \in [0, 1)$, snapping into all $225$ new kissing contacts per ball simultaneously at $s = 1$.

### Theorem 2 (Exact $E_6 / E_7 / D_n^+$ Uniqueness and Obstruction Theorem)

Among the exceptional root lattices $E_6, E_7, E_8$ (normalized to minimum squared norm $2$) and the entire family $D_n^+ = D_n \cup (D_n + (\tfrac{1}{2})^n)$ ($n \ge 4$ even):
1. **($E_6$ and $E_7$ determinant-square obstruction)** Neither $E_6$ nor $E_7$ admits *any* periodic cross-polytope packing $P_n + \Lambda = \Gamma + v_0$ (even without bar constraints).
2. **($D_n^+$ bar-preserving uniqueness)** A bar-preserving cross-polytope fold $F_n \to D_n^+$ exists **if and only if $n = 8$** ($D_8^+ = E_8$).

---

## 2. Complete Proof of Theorem 1

### Lemma 1 (Exact Bar Preservation & Endpoints)

Write $r_i(s) - 2 e_i = u_i(s)$. For $i \in I_\pm$, the diagonal entry of $u_i(s)$ is $(a_\pm(s) - 2) = \frac{N_{a,\pm}(s) - 2 D_\pm(s)}{D_\pm(s)}$ and each of the $7$ off-diagonal entries has squared value $b_\pm(s)^2 = \frac{8 s^2}{D_\pm(s)^2}$. Expanding directly in $\mathbb{Z}[\sqrt{2}][s]$:
- For $I_+$: $N_{a,+}(s) - 2 D_+(s) = (4 + \sqrt{2}) + (-4 + \sqrt{2}) s^2$, so
  $$\bigl((4+\sqrt{2}) + (-4+\sqrt{2}) s^2\bigr)^2 + 7 \cdot 8 s^2 = (18 + 8\sqrt{2}) + (28) s^2 + (18 - 8\sqrt{2}) s^4 = 2\bigl((1+2\sqrt{2}) + (-1+2\sqrt{2}) s^2\bigr)^2 = 2 D_+(s)^2.$$
- For $I_-$: $N_{a,-}(s) - 2 D_-(s) = (4 - \sqrt{2}) + (-4 - \sqrt{2}) s^2$, so
  $$\bigl((4-\sqrt{2}) + (-4-\sqrt{2}) s^2\bigr)^2 + 7 \cdot 8 s^2 = (18 - 8\sqrt{2}) + (28) s^2 + (18 + 8\sqrt{2}) s^4 = 2\bigl((-1+2\sqrt{2}) + (1+2\sqrt{2}) s^2\bigr)^2 = 2 D_-(s)^2.$$

Since $D_\pm(s)$ have strictly positive Bernstein coefficients on $[0, 1]$, $\|r_i(s) - 2 e_i\|^2 \equiv 2$ identically. Evaluating at $s = 0$ gives $b_\pm(0) = 0$ and
$$a_+(0) = \frac{6+5\sqrt{2}}{1+2\sqrt{2}} = \frac{(6+5\sqrt{2})(2\sqrt{2}-1)}{7} = 2 + \sqrt{2}, \qquad a_-(0) = \frac{2+3\sqrt{2}}{-1+2\sqrt{2}} = \frac{(2+3\sqrt{2})(2\sqrt{2}+1)}{7} = 2 + \sqrt{2},$$
so $L_8(0) = (2+\sqrt{2}) I_8$. At $s = 1$, $D_+(1) = D_-(1) = 4\sqrt{2}$, $N_{a,+}(1) = 10\sqrt{2}$, $N_{a,-}(1) = 6\sqrt{2}$, $N_{b,\pm}(1) = 2\sqrt{2}$, giving $a_+(1) = 5/2$, $a_-(1) = 3/2$, $b_\pm(1) = 1/2$, i.e. $L_8(1) = 2 I_8 + \frac{1}{2} H_8$. $\square$

### Lemma 2 (Exact $E_8$ Endpoint, $(\mathbb{Z}/2\mathbb{Z})^4$ Quotient, & Kissing 240 — Computer-Free Proof)

1. **Inclusion $\mathbb{Z}^8 L_8(1) \subset E_8$:** Row $0$ of $H_8$ has $0$ minus signs and rows $1,\dots,7$ each have $4$ minus signs (all even). Thus every row $c_i = \frac{1}{2}(H_8)_i \in \{\pm \frac{1}{2}\}^8$ is an $E_8$ spinor root, and $r_i(1) = 2 e_i + c_i \in E_8$.
2. **Exact inverse, determinant, and Gram matrix:** Since $H_8^T = H_8$ and $H_8^2 = 8 I_8$, the eigenvalues of $H_8$ are $+2\sqrt{2}$ (mult $4$) and $-2\sqrt{2}$ (mult $4$). Hence $L_8(1) = 2 I_8 + \frac{1}{2} H_8$ is symmetric positive definite with eigenvalues $2 + \sqrt{2}$ (mult $4$) and $2 - \sqrt{2}$ (mult $4$), so
   $$\det L_8(1) = (2+\sqrt{2})^4 (2-\sqrt{2})^4 = (4 - 2)^4 = 16,$$
   $$\left(2 I_8 + \frac{1}{2} H_8\right)\left(2 I_8 - \frac{1}{2} H_8\right) = 4 I_8 - 2 I_8 = 2 I_8 \implies L_8(1)^{-1} = I_8 - \frac{1}{4} H_8,$$
   $$G(1) := L_8(1) L_8(1)^T = L_8(1)^2 = 4 I_8 + 2 H_8 + \frac{1}{4}(8 I_8) = 6 I_8 + 2 H_8.$$
   Because $G(1)_{ii} \in \{4, 8\}$ and $G(1)_{ij} = \pm 2$ ($i \ne j$), $m^T G(1) m \in 4\mathbb{Z}_{>0}$ for all $m \in \mathbb{Z}^8 \setminus \{0\}$. Indeed, since $L_8(1)$ has symmetric positive-definite inverse $I_8 - \frac{1}{4} H_8$, we have $2 E_8 \subset \mathbb{Z}^8 L_8(1)$: for every $x \in E_8$,
   $$2 x L_8(1)^{-1} = 2 x - \frac{1}{2} x H_8 \in \mathbb{Z}^8$$
   (if $x \in D_8$, $\sum_j (\pm 1) x_j \equiv \sum_j x_j \equiv 0 \pmod 2$; if $x \in D_8 + (\tfrac{1}{2})^8$, then $2x$ is all-odd with sum $\equiv 0 \pmod 4$, and flipping $4$ odd signs changes the sum by a multiple of $4$, so $(2x) H_8 \equiv 0 \pmod 4$). Thus $[E_8 : \mathbb{Z}^8 L_8(1)] = \det L_8(1) / \det E_8 = 16$ and $E_8 / \mathbb{Z}^8 L_8(1) \cong (\mathbb{Z}/2\mathbb{Z})^4$.
3. **Exact tiling $P_8 + \mathbb{Z}^8 L_8(1) = E_8 + e_0$:** Every $p - e_0 \in P_8 - e_0 \subset D_8 \subset E_8$. For any distinct $p, q \in P_8$:
   - If $p - q = \pm 2 e_i$, then $(p-q) L_8(1)^{-1} = \pm 2 e_i \mp \frac{1}{2} (H_8)_i \in (\mathbb{Z} + \frac{1}{2})^8 \notin \mathbb{Z}^8$.
   - If $p - q = \pm e_i \pm e_j$ ($i \ne j$), then $(p-q) L_8(1)^{-1} = \pm e_i \pm e_j - \frac{1}{4}(\pm (H_8)_i \pm (H_8)_j)$. By orthogonality of distinct rows of $H_8$, $\pm (H_8)_i \pm (H_8)_j$ has four entries equal to $\pm 2$ and four equal to $0$, leaving four half-integer entries $\notin \mathbb{Z}^8$.
   Hence the $16$ vectors $P_8 - e_0$ represent all $16$ distinct cosets of $E_8 / \mathbb{Z}^8 L_8(1)$, proving $P_8 + \mathbb{Z}^8 L_8(1) = E_8 + e_0$.
4. **Density and kissing number:** Every ball in $E_8 + e_0$ has radius $\sqrt{2}/2$, packing density $16 \cdot \frac{\pi^4 (\sqrt{2}/2)^8}{24} / 16 = \frac{\pi^4}{384}$ (optimal in $\mathbb{R}^8$, Viazovska 2016), and exact kissing number $240$. Moreover $\mathbb{Z}^8 L_8(1)$ is isometric to $\sqrt{2}\, E_8$ (with theta-series shells $240$ at norm$^2 = 4$ and $2160$ at norm$^2 = 8$). $\square$

### Lemma 3 (Exact $2\times 2$ Block-Diagonalization & Strict Volume Monotonicity)

Order the coordinates as $\pi = (0, 3, 5, 6 \mid 7, 4, 2, 1)$, pairing each even-weight $x \in I_+ \subset \mathbb{F}_2^3$ with its complement $x \oplus 111 \in I_-$. Since $(x \oplus \mathbf{1})\cdot y = x\cdot y + |y|$ and $(x\oplus\mathbf{1})\cdot(y\oplus\mathbf{1}) \equiv x\cdot y + |x| + |y| + 1 \pmod 2$, the permuted Hadamard matrix splits identically into $4\times 4$ blocks:
$$P^T H_8 P = \begin{pmatrix} K_4 & K_4 \\ K_4 & -K_4 \end{pmatrix}, \qquad K_4 = \begin{pmatrix} 1 & 1 & 1 & 1 \\ 1 & 1 & -1 & -1 \\ 1 & -1 & 1 & -1 \\ 1 & -1 & -1 & 1 \end{pmatrix}, \qquad K_4^T = K_4,\; K_4^2 = 4 I_4.$$
The pure-integer orthogonal matrix
$$Q_4 = \begin{pmatrix} 1 & 1 & 1 & -1 \\ 1 & -1 & 1 & 1 \\ -1 & 1 & 1 & 1 \\ 1 & 1 & -1 & 1 \end{pmatrix}, \qquad Q_4^T Q_4 = 4 I_4, \qquad Q_4^T K_4 Q_4 = 4 \operatorname{diag}(2, 2, 2, -2)$$
simultaneously diagonalizes $K_4$ on both $I_+$ and $I_-$ blocks for **every** $s \in [0, 1]$! Letting $Q_8$ interleave the columns of $I_2 \otimes Q_4$ into the four $2\times 2$ eigenspaces yields the exact polynomial matrix identity
$$\frac{1}{4} Q_8^T \bigl(P^T L_8(s) P\bigr) Q_8 = \operatorname{diag}\!\bigl(M_{+2}(s),\; M_{+2}(s),\; M_{+2}(s),\; M_{-2}(s)\bigr),$$
where
$$M_{+2}(s) = \begin{pmatrix} a_+(s) + b_+(s) & 2 b_+(s) \\ 2 b_-(s) & a_-(s) - b_-(s) \end{pmatrix}, \qquad M_{-2}(s) = \begin{pmatrix} a_+(s) - 3 b_+(s) & -2 b_+(s) \\ -2 b_-(s) & a_-(s) + 3 b_-(s) \end{pmatrix}.$$
Writing $\det M_{+2}(s) = \frac{N_+(s)}{D_+(s) D_-(s)}$ and $\det M_{-2}(s) = \frac{N_-(s)}{D_+(s) D_-(s)}$ with degree-4 numerators in $\mathbb{Z}[\sqrt{2}][s]$:
- $N_+(s) = (42+28\sqrt{2}) + (-8-8\sqrt{2}) s - 4 s^2 + (-8+8\sqrt{2}) s^3 + (42-28\sqrt{2}) s^4$, all 5 Bernstein coefficients strictly $> 0$;
- $N_-(s) = (42+28\sqrt{2}) + (24+24\sqrt{2}) s - 68 s^2 + (24-24\sqrt{2}) s^3 + (42-28\sqrt{2}) s^4$, all 5 Bernstein coefficients strictly $> 0$.

Therefore $\det L_8(s) = \frac{N_+(s)^3 N_-(s)}{(D_+(s) D_-(s))^4} > 0$, and its logarithmic derivative simplifies to
$$-\frac{d}{ds} \det L_8(s) = \frac{s \, N_+(s)^2 \, Q_{\det}(s)}{(D_+(s) D_-(s))^5},$$
where the exact degree-10 polynomial $Q_{\det}(s) \in \mathbb{Z}[\sqrt{2}][s]$ has **all 11 Bernstein coefficients strictly positive** on $[0, 1]$ (verified exactly by [`verify_e8_spinor_fold.py`](./verify_e8_spinor_fold.py)). Thus $\det L_8(s)$ is strictly decreasing on $[0, 1]$ from $(2+\sqrt{2})^8 = 9232 + 6528\sqrt{2}$ to $16$. $\square$

### Lemma 4 (Exact Loewner Monotonicity $G(s) \succeq G(1)$ & Zero-Gap Infinite-Image Reduction)

By Lemma 3, $G(s) - G(1) = L_8(s) L_8(s)^T - L_8(1) L_8(1)^T$ orthogonally block-diagonalizes into three copies of $G_{+2}(s) - G_{+2}(1)$ and one copy of $G_{-2}(s) - G_{-2}(1)$, where $G_{+2}(1) = \begin{pmatrix} 10 & 4 \\ 4 & 2 \end{pmatrix}$ and $G_{-2}(1) = \begin{pmatrix} 2 & -4 \\ -4 & 10 \end{pmatrix}$. Clearing the positive denominators $D_+(s)^2$ and $D_-(s)^2$ on each $2\times 2$ block $\lambda \in \{+2, -2\}$:
- $(G_\lambda(s) - G_\lambda(1))_{00} \cdot D_+(s)^2 = (1 - s)\, P_{00}^{(\lambda)}(s)$ (degree 3: all 4 Bernstein coefficients strictly $> 0$);
- $(G_\lambda(s) - G_\lambda(1))_{11} \cdot D_-(s)^2 = (1 - s)\, P_{11}^{(\lambda)}(s)$ (degree 3: all 4 Bernstein coefficients strictly $> 0$);
- $\det(G_\lambda(s) - G_\lambda(1)) \cdot D_+(s)^2 D_-(s)^2 = s (1 - s)^2\, P_{\det}^{(\lambda)}(s)$ (degree 5: all 6 Bernstein coefficients strictly $> 0$; the factor $s$ reflects that the upper 4D eigenspace of $L_8(1)$ retains the exact initial scale $2+\sqrt{2}$).

By Sylvester's criterion (equivalently the exact $2\times 2$ Schur complement $d_{11} - d_{01}^2/d_{00} \ge 0$), $G(s) \succeq G(1) = 6 I_8 + 2 H_8$ for every $s \in [0, 1]$. Consequently, for every $m \in \mathbb{Z}^8 \setminus \{0\}$, $v = p - q \in P_8 - P_8$, and $s \in [0, 1]$:
$$\|m L_8(s)\|^2 = m^T G(s) m \ge m^T G(1) m \in 4\mathbb{Z}_{>0}.$$
Three exhaustive cases cover all of $\mathbb{Z}^8 \setminus \{0\}$ with zero heuristic truncation:
1. **$v = 0$:** $\|m L_8(s)\|^2 \ge 4 > 2$ analytically for all $m \ne 0$ and all $s \in [0, 1]$.
2. **$\|v\|^2 = 2$ ($v = \pm e_i \pm e_j$, $112$ vectors):** If $m^T G(1) m \ge 8$, then $\|m L_8(s)\| \ge 2\sqrt{2}$, so by the triangle inequality $\|m L_8(s) + v\| \ge 2\sqrt{2} - \sqrt{2} = \sqrt{2}$ analytically for all $s \in [0, 1]$. Moreover, strict inequality $\|m L_8(s) + v\| > \sqrt{2}$ holds on all of $[0, 1)$: on $(0, 1)$ we have $G(s) - G(1) \succ 0$ strictly (so $\|m L_8(s)\|^2 > 8$), and at $s = 0$ every $m \in \mathbb{Z}^8 \setminus \{0\}$ satisfies $\|m L_8(0)\|^2 = (6 + 4\sqrt{2}) \|m\|^2 \ge 6 + 4\sqrt{2} > 8$ (so equality holds only at $s = 1$ for the $56$ anti-parallel pairs $m L_8(1) = -2 v$). Only $m^T G(1) m < 8$ (i.e. $m^T G(1) m = 4$) requires polynomial certification.
3. **$\|v\|^2 = 4$ ($v = \pm 2 e_i$, $16$ vectors):** Since $(2+\sqrt{2})^2 = 6 + 4\sqrt{2} < 12$, if $m^T G(1) m \ge 12$, then $\|m L_8(s) + v\| \ge \sqrt{12} - 2 > \sqrt{2}$ analytically for all $s \in [0, 1]$. Only $m^T G(1) m \le \lfloor (2+\sqrt{2})^2 \rfloor = 11$ (i.e. $m^T G(1) m \in \{4, 8\}$) requires polynomial certification. $\square$

### Lemma 5 (Exact Zero-Float Ellipsoid & Zero-Subdivision Degree-8 $\mathbb{Z}[\sqrt{2}]$ Bernstein Certificate)

The pure-rational $LDL^T$ factorization of $G(1) = 6 I_8 + 2 H_8$ has strictly positive rational Schur-complement pivots
$$D = \left(8,\; \frac{7}{2},\; \frac{20}{7},\; 5,\; \frac{4}{5},\; \frac{7}{5},\; \frac{8}{7},\; \frac{1}{2}\right).$$
Exact integer `math.isqrt` Fincke–Pohst enumeration over $m^T G(1) m \le 11$ finds, up to $\pm m$ symmetry, **exactly $120$ pairs** with $m^T G(1) m = 4$ ($240/2$ roots of $\sqrt{2}\,E_8$), **exactly $1,080$ pairs** with $m^T G(1) m = 8$ ($2160/2$ second-shell vectors), and **$0$ pairs** at norms $1,2,3,5,6,7,9,10,11$.

Clearing the common degree-4 denominator $D_+(s) D_-(s)$ turns $\|m L_8(s) + v\|^2 - 2 \ge 0$ into a degree-8 polynomial $P_{m,v}(s) \in \mathbb{Z}[\sqrt{2}][s]$. Over all
$$120 \times 112 + 1{,}200 \times 16 = 13{,}440 + 19{,}200 = \mathbf{32{,}640}$$
canonical pairs ($293,760$ exact $\mathbb{Z}[\sqrt{2}]$ degree-8 Bernstein coefficients):
- **Negative Bernstein coefficients:** **$0$** (zero subdivisions needed!).
- **Exact zero Bernstein coefficients:** **$976$**, completely accounted for by exact contact geometry:
  - $8 \times 9 = 72$ zeros from the $8$ persistent inter-cell bar pairs $(m = e_i, v = -2 e_i)$ (identically zero on $[0, 1]$);
  - $120$ zeros at $k = 8$ ($s = 1$) from non-bar $V_{2e}$ kissing contacts;
  - $784$ zeros at $k = 8$ ($s = 1$) from $V_{ee}$ ($m^T G(1) m = 4$) kissing contacts;
  - $0$ interior zeros ($k \in \{0,\dots,7\}$) on any of the $32,632$ non-bar pairs, proving **strict separation** $\|m L_8(s) + v\|^2 > 2$ on all of $s \in [0, 1)$.
- **Exact kissing accounting at every one of the 16 motif balls:**
  $$14\;\text{(intra-cell)} + 1\;\text{(bar)} + 15\;\text{(non-bar }V_{2e}\text{)} + 196\;\text{(}V_{ee},\, m^T G(1) m = 4\text{)} + 14\;\text{(anti-parallel }V_{ee},\, m L_8(1) = -2v\text{)} = \mathbf{240}.$$
This completes the proof of Theorem 1. $\square$

---

## 3. Proof of Theorem 2 (Exact $E_6 / E_7 / D_n^+$ Uniqueness & Obstruction)

1. **$E_6$ and $E_7$ obstruction:** Suppose $P_n + \Lambda = \Gamma + v_0$ for a lattice $\Gamma$ with minimum squared norm $2$ (matching the cross-polytope edge norm$^2 = 2$). Then every difference $p - q \in P_n - P_n$ belongs to $\Gamma$, so the generated lattice $D_n = \langle P_n - P_n \rangle$ is a finite-index sublattice of $\Gamma$. Letting $k = [\Gamma : D_n] \in \mathbb{Z}_{\ge 1}$, the Gram determinants satisfy
   $$\frac{\det_{\mathrm{Gram}}(D_n)}{\det_{\mathrm{Gram}}(\Gamma)} = k^2 \in \mathbb{Z}^2.$$
   At minimum squared norm $2$, $\det_{\mathrm{Gram}}(D_n) = 4$, whereas $\det_{\mathrm{Gram}}(E_6) = 3$, $\det_{\mathrm{Gram}}(E_7) = 2$, $\det_{\mathrm{Gram}}(E_8) = 1$. Since $4/3 \notin \mathbb{Z}^2$ and $4/2 = 2 \notin \mathbb{Z}^2$, neither $E_6$ nor $E_7$ can even contain a rescaled $D_n$ sublattice at the same minimum distance, making $F_6 \to E_6$ and $F_7 \to E_7$ impossible even without bar constraints.
2. **$D_n^+$ bar-preserving uniqueness ($n = 8$ only):** For even $n \ge 4$, $D_n^+ = D_n \cup (D_n + (\frac{1}{2})^n)$:
   - If $n \in \{4, 6\}$, the spinor coset contains vectors of squared norm $n/4 \in \{1, 3/2\} < 2$, so spheres of radius $\sqrt{2}/2$ strictly overlap.
   - If $n \ge 10$ (even), every vector $x \in D_n^+ \setminus D_n$ has $\|x\|^2 \ge n/4 \ge 5/2 > 2$. In any bar-preserving periodic fold $r_i = 2 e_i + c_i \in D_n^+$ with $\|c_i\|^2 = 2$, we must have $c_i \in D_n^+$, hence $c_i \in D_n$, hence $r_i \in D_n$ for every $i = 0,\dots,n-1$. Thus $\Lambda \subset D_n$, forcing $P_n + \Lambda \subset D_n + e_0 \subsetneq D_n^+ + e_0$.
   - Only at $n = 8$ does $n/4 = 2$ match the bar length identically, yielding $E_8 = D_8^+$ via Theorem 1. $\square$

*Remark on 24D Leech Lattice $\Lambda_{24}$ (open frontier):* Rescaling $\Lambda_{24}$ around a $D_{24}$ deep hole to $\Gamma = \Lambda_{24}/\sqrt{2}$ (minimum squared norm $2$) gives $[\Gamma : D_{24}] = 2^{13}$ and $|P_{24}| = 48 = 16 \cdot 3$. Because $\gcd(3, 2^{13}) = 1$, every surjective homomorphism $D_{24} \to \mathbb{Z}/3\mathbb{Z}$ extends uniquely to $\Gamma \to \mathbb{Z}/3\mathbb{Z}$, so the finite-group index obstruction vanishes; whether 24 minimal Leech vectors can generate an overlap-free bar-preserving fold $F_{24} \to \Lambda_{24}/\sqrt{2}$ remains an attractive open problem.

---

## 4. Context, Prior-Art Distinction, & Reproducibility

- **Origin & credit:** Sparked by `@session-cartographer`'s Portal 9 (`#46`), which identified a cyclic $\mathbb{Z}/16\mathbb{Z}$ character on $E_8$ and numerically observed 5 overlap-free spinor curves over $\{-1, 0, 1\}^8$ ($\sigma_{\min} \approx 0.431$). Passing from cyclic $\mathbb{Z}/16\mathbb{Z}$ to the canonical Sylvester Hadamard matrix $H_8 = H_2^{\otimes 3}$ (quotient $(\mathbb{Z}/2\mathbb{Z})^4$, $\sigma_{\min} = 2 - \sqrt{2} \approx 0.586$) collapses the 8D fold into two exact $2\times 2$ blocks, proves exact Loewner monotonicity $G(s) \succeq 6 I_8 + 2 H_8$, and certifies all of $\mathbb{Z}^8$ in under $1$ second with zero subdivisions and zero floats.
- **Exact verification command (all 4 R4 modes):**
  ```bash
  cd results/expedition-dn-e8-folds
  python3 -I -B verify_e8_spinor_fold.py
  python3 -I -O -B verify_e8_spinor_fold.py
  python3 -I -B -W error verify_e8_spinor_fold.py
  python3 -I -O -B -W error verify_e8_spinor_fold.py
  ```
