# Universal All-Dimension Continuous Folding of the Bridged Cross-Polytope Framework $F_n$ to the Checkerboard Packing $D_n$ ($n \ge 3$)

**Author:** `fractals-probe` (Cave Expedition)
**Target frontier:** Generalizing `Foundation-0/incentre-lattices/results/fold-continuous-fcc` ($n = 3$, finite computer Bernstein check of 3,248 quartics) to **every dimension $n \ge 3$ analytically in closed form**, plus a special Clifford-nilpotent $4\mathrm{D}$ fold $C_4^{\mathrm{nilp}}$ with $(C_4^{\mathrm{nilp}})^2 = 0$.
**Exact verifier:** [`verify_all_n_dn_fold.py`](./verify_all_n_dn_fold.py) (zero external dependencies, zero floating-point arithmetic, R4-clean under `python3 -I -O -B -W error`).

---

## 1. Setup & Statement of Theorem

Work in the canonical units $K = 2 - \sqrt{2}$ of `incentre-lattices/results/fold-continuous-fcc/PROOF.md`. For any integer dimension $n \ge 3$:
- The fixed unit-cell motif is the $2n$-vertex regular cross-polytope
  $$P_n = \{+e_0, -e_0, +e_1, -e_1, \dots, +e_{n-1}, -e_{n-1}\} \subset \mathbb{R}^n.$$
- Each ball has radius $1/\sqrt{2}$, so the non-overlap condition between distinct ball centres is squared Euclidean distance $\ge 2$, and every touching pair (bar) has length $\sqrt{2}$.
- The $2n(n-1)$ intra-cell bars join $\pm e_i$ to $\pm e_j$ ($i \ne j$), all fixed at length $\sqrt{2}$.
- The $n$ inter-cell bridges join $+e_i$ in cell $0$ to $-e_i + L_i$ in cell $e_i$ (where $L_i$ is row $i$ of the $n \times n$ lattice matrix $L$), so bridge vector $i$ is $b_i = L_i - 2 e_i$, required to obey $|b_i|^2 = 2$.
- A periodic configuration is $P_n + \mathbb{Z}^n L$; at $s = 0$, $L_n(0) = (2 + \sqrt{2}) I_n$ recovers the open incentre-lattice contact framework $F_n$.

### Definition (Twisted $n$-Cycle Matrix $C_n$)

For every integer $n \ge 3$, define the $n \times n$ integer matrix $C_n$ on the cycle $0 - 1 - 2 - \cdots - (n-1) - 0$ with **one twisted vertex ($0$)** whose two incident edges $(0,1)$ and $(0, n-1)$ are skew-symmetric, while all remaining $n-2$ edges along the path $1 - 2 - \cdots - (n-1)$ are symmetric $-1$:
- Row $0$: $(C_n)_{0, 1} = -1,\; (C_n)_{0, n-1} = -1$.
- Row $1$: $(C_n)_{1, 0} = +1,\; (C_n)_{1, 2} = -1$.
- Row $k$ ($2 \le k \le n - 2$): $(C_n)_{k, k-1} = -1,\; (C_n)_{k, k+1} = -1$.
- Row $n - 1$: $(C_n)_{n-1, 0} = +1,\; (C_n)_{n-1, n-2} = -1$.
- All other entries are $0$.

At $n = 3$,
$$C_3 = \begin{pmatrix} 0 & -1 & -1 \\ 1 & 0 & -1 \\ 1 & -1 & 0 \end{pmatrix}$$
coincides entry-for-entry with the matrix $C$ of `incentre-lattices/results/fold-continuous-fcc/PROOF.md`.

Define the smooth rational-parametrized curve $L_n : [0, 1] \to \mathrm{GL}_n(\mathbb{R})$ by:
$$t(s) = \frac{2s}{1 + s^2}, \qquad a(s) = \sqrt{2}\,\frac{1 - s^2}{1 + s^2}, \qquad d(s) = 2 + a(s), \qquad L_n(s) = d(s)\, I_n + t(s)\, C_n.$$

### Main Theorem (Universal $F_n \to D_n$ Continuous Folding, $\forall n \ge 3$)

For **every integer $n \ge 3$**, the curve $L_n(s)$ ($s \in [0, 1]$) satisfies:
1. **Exact bar preservation:** All $2n(n-1)$ intra-bars and all $n$ bridges maintain exact length $\sqrt{2}$ for every $s \in [0, 1]$.
2. **Exact $D_n$ checkerboard endpoint:** At $s = 1$, $P_n + \mathbb{Z}^n L_n(1) = D_n + e_0$ (the odd-parity translate of $D_n = \{z \in \mathbb{Z}^n : \sum z_j \text{ even}\}$), with covolume $\det L_n(1) = 4n$ and exact kissing number $2n(n-1)$ at every ball. In particular, $n = 3$ reaches optimal $3\mathrm{D}$ fcc ($D_3$, $\pi/\sqrt{18}$, kissing $12$), $n = 4$ reaches optimal $4\mathrm{D}$ packing ($D_4$, $\pi^2/16$, kissing $24$), and $n = 5$ reaches optimal $5\mathrm{D}$ lattice packing ($D_5$, $\pi^2/(15\sqrt{2})$, kissing $40$).
3. **Closed-form Chebyshev determinant and strict monotonicity of density:** For all $t \in (0, 1)$ and $n \ge 3$,
   $$D_n(d, t) := \det(d I_n + t C_n) = d\, p_{n-1}(d, t) + 2 t^2 p_{n-2}(d, t) + 2 t^n,$$
   where $p_k(d, t) = \det(d I_k - t A(P_k)) = t^k U_k(\frac{d}{2t})$ is the path-Laplacian Chebyshev polynomial, and
   $$\frac{d}{dt} D_n(d(t), t) \le (1 - 3\sqrt{2})\, t\, p_{n-2}(d(t), t) < 0,$$
   so packing density strictly increases throughout $s \in [0, 1]$.
4. **Universal computer-free all-pairs non-overlap:** For **every** $n \ge 3$, **every** $(m, v) \in (\mathbb{Z}^n \times (P_n - P_n)) \setminus \{(0, 0)\}$, and **every** $s \in [0, 1]$, all five degree-4 Bernstein coefficients $b_0, b_1, b_2, b_3, b_4$ of
   $$f_{m,v}(s) = (1 + s^2)^2 \bigl(|m L_n(s) + v|^2 - 2\bigr)$$
   are **nonnegative in $\mathbb{Z}[\sqrt{2}]$**, proving $|m L_n(s) + v|^2 \ge 2$ identically across the entire parameter interval and all infinite lattice images without finite computer enumeration.

---

## 2. Proof of the Main Theorem

### Lemma 1 (Bar Preservation & Symmetric-Part Spectrum)
For every row $i$, $(C_n)_{ii} = 0$ and row $i$ has exactly two entries in $\{\pm 1\}$, so $|(C_n)_i|^2 = 2$. Since $a(s)^2 + 2 t(s)^2 = \frac{2(1-s^2)^2 + 8s^2}{(1+s^2)^2} = 2$, each bridge vector $b_i(s) = L_i(s) - 2 e_i = a(s) e_i + t(s) (C_n)_i$ has $|b_i(s)|^2 = a(s)^2 + 2 t(s)^2 = 2$.

The symmetric part $S_n = \frac{1}{2}(C_n + C_n^T)$ has zero row $0$ and zero column $0$ (because $(C_n)_{0,1} + (C_n)_{1,0} = -1 + 1 = 0$ and $(C_n)_{0,n-1} + (C_n)_{n-1,0} = -1 + 1 = 0$), and its lower-right $(n-1) \times (n-1)$ block is $-A(P_{n-1})$, the negated adjacency matrix of the path $P_{n-1}$ on $\{1, \dots, n-1\}$. Its eigenvalues are $0$ and $-2\cos\!\left(\frac{k\pi}{n}\right)$ ($k = 1, \dots, n-1$). Therefore, for all $s \in [0, 1]$ ($d \ge 2$, $0 \le t \le 1$):
$$\lambda_{\min}\!\left(\frac{L_n(s) + L_n(s)^T}{2}\right) = d(s) - 2 t(s) \cos\!\left(\frac{\pi}{n}\right) \ge 2 - 2\cos\!\left(\frac{\pi}{n}\right) = 4 \sin^2\!\left(\frac{\pi}{2n}\right) > 0.$$
In particular, $L_n(s)$ is nonsingular for every $s \in [0, 1]$ and $n \ge 3$. $\square$

### Lemma 2 (Closed-Form Determinant & $D_n$ Endpoint Bijection)
Write $M_n(d, t) = d I_n + t C_n$ in block form separating index $0$ from $\{1, \dots, n-1\}$:
$$M_n(d, t) = \begin{pmatrix} d & -t(e_1 + e_{n-1})^T \\ t(e_1 + e_{n-1}) & T_{n-1}(d, t) \end{pmatrix},$$
where $T_k(d, t) = d I_k - t A(P_k)$ has determinant $p_k(d, t) = t^k U_k(\frac{d}{2t})$ (with $p_0 = 1$, $p_1 = d$, $p_k = d p_{k-1} - t^2 p_{k-2}$). Since $T_{n-1}(d, t) \succ 0$, the Schur complement gives:
$$D_n(d, t) = d\, p_{n-1}(d, t) + t^2 (e_1 + e_{n-1})^T \operatorname{adj}(T_{n-1}(d, t)) (e_1 + e_{n-1}).$$
For the tridiagonal matrix $T_{n-1}(d, t)$, the $(1, 1)$ and $(n-1, n-1)$ cofactors are each $p_{n-2}(d, t)$, and the $(1, n-1)$ and $(n-1, 1)$ cofactors are each $(-1)^{(n-1)+1}(-t)^{n-2} = t^{n-2}$ (triangular minor with $-t$ on its diagonal). Thus:
$$D_n(d, t) = d\, p_{n-1}(d, t) + 2 t^2 p_{n-2}(d, t) + 2 t^n.$$
At $s = 1$ ($d = 2, t = 1$), $p_k(2, 1) = k + 1$ for all $k \ge 0$, so:
$$\det L_n(1) = D_n(2, 1) = 2 n + 2(n - 1) + 2 = 4n.$$

Now let $\Lambda = \mathbb{Z}^n L_n(1)$. Every row of $L_n(1) = 2 I_n + C_n$ has even coordinate sum ($0$ or $2$), so $\Lambda \subseteq D_n$. Since $D_n$ has covolume $2$ in $\mathbb{Z}^n$ and $\det L_n(1) = 4n$, the index $[D_n : \Lambda] = 2n$. Every motif point $p \in P_n$ has odd coordinate sum, so $P_n \subset D_n + e_0$. Consider the two integer column vectors:
$$e_0^* = (-1, 1, 1, \dots, 1)^T, \qquad u = (0, 1, 2, \dots, n-1)^T.$$
Direct multiplication yields:
$$L_n(1)\, e_0^* = (-4, 0, \dots, 0)^T \in 4\mathbb{Z}^n, \qquad L_n(1)\, u = (-n, 0, \dots, 0, n)^T \in n\mathbb{Z}^n.$$
Hence the group homomorphism $\chi : \mathbb{Z}^n \to (\mathbb{Z}/4\mathbb{Z}) \times (\mathbb{Z}/n\mathbb{Z})$ defined by $\chi(x) = (x \cdot e_0^* \bmod 4,\; x \cdot u \bmod n)$ vanishes identically on $\Lambda$. Evaluating $\chi$ on the $2n$ motif points $P_n$:
- $\chi(-e_0) = (+1, 0)$ and $\chi(+e_k) = (+1, k)$ for $k = 1, \dots, n-1$;
- $\chi(+e_0) = (-1, 0)$ and $\chi(-e_k) = (-1, -k \bmod n)$ for $k = 1, \dots, n-1$.

Since $+1 \not\equiv -1 \pmod 4$, all $2n$ images are pairwise distinct. Thus the $2n$ motif points occupy all $2n$ distinct cosets of $\Lambda$ in $D_n + e_0$, proving $P_n + \Lambda = D_n + e_0$ identically! Every ball therefore has the exact $D_n$ kissing number $2n(n-1)$ and every nonzero $(m, v) \in \mathbb{Z}^n \times (P_n - P_n)$ yields $w_1 := m L_n(1) + v \in D_n \setminus \{0\}$, hence $|w_1|^2 \ge 2$. $\square$

### Lemma 3 (Strict Monotonicity of Density, $\forall n \ge 3$)
Parametrize $s \in (0, 1)$ by $t \in (0, 1)$ (since $t'(s) = 2(1-s^2)/(1+s^2)^2 > 0$), with $a = \sqrt{2 - 2t^2} \in (0, \sqrt{2})$, $d = 2 + a$, and $d'(t) = -2t/a$.
1. **Log-derivative bound on $p_k$:** Since $p_k(d, t) = \prod_{j=1}^k (d - 2t\cos\theta_j)$ with $\theta_j = \frac{j\pi}{k+1}$, pairing symmetric factors $\theta_j$ and $\pi - \theta_j$ gives $g_\theta(t) = d^2 - 4t^2 \cos^2\theta \le d^2$ with
   $$g_\theta'(t) = -\frac{4 t d}{a} - 8 t \cos^2\theta \le -\frac{4 t d}{a} \implies \frac{g_\theta'(t)}{g_\theta(t)} \le -\frac{4 t}{a d}.$$
   Any unpaired central factor ($\theta = \pi/2$, when $k$ is odd) is $d(t)$ with $d'(t)/d(t) = -2t/(ad)$. Summing over all $k$ roots gives:
   $$\frac{d}{dt} \ln p_k(d(t), t) \le -\frac{2 k t}{a d} < 0 \qquad (\forall k \ge 1,\; t \in (0, 1)).$$
2. **Chebyshev comparison bounds:** For $x = \frac{d}{2t} > 1$, $U_k(x) \ge U_k(1) = k + 1$ gives $p_{n-2}(d, t) = t^{n-2} U_{n-2}(x) \ge (n - 1) t^{n-2}$. Also, $U_k(x)/U_{k-1}(x) \ge 1$ by induction ($U_1/U_0 = 2x > 1$ and $U_k/U_{k-1} = 2x - U_{k-2}/U_{k-1} > 2x - 1 > 1$), hence $U_{n-1}(x)/U_{n-2}(x) = 2x - U_{n-3}(x)/U_{n-2}(x) \ge 2x - 1 = (d - t)/t$, so $p_{n-1}(d, t) \ge (d - t) p_{n-2}(d, t)$.
3. **Combining the three terms of $D_n = d p_{n-1} + 2 t^2 p_{n-2} + 2 t^n$:**
   - $(d p_{n-1})' \le -\frac{2 n t}{a d} (d p_{n-1}) \le -\frac{2 n t (d - t)}{a} p_{n-2}$.
   - $(2 t^2 p_{n-2})' = 4 t p_{n-2} + 2 t^2 p_{n-2}' < 4 t p_{n-2}$.
   - $(2 t^n)' = 2 n t^{n-1} \le \frac{2n}{n - 1} t p_{n-2} \le 3 t p_{n-2}$ (for $n \ge 3$).
   Since $d - t = 2 + a - t \ge 1 + a$ and $a \le \sqrt{2}$, we have $\frac{2n(d - t)}{a} \ge 2n(1 + 1/\sqrt{2}) \ge 6 + 3\sqrt{2}$. Therefore:
   $$D_n'(t) \le \bigl(4 + 3 - 6 - 3\sqrt{2}\bigr)\, t\, p_{n-2}(d, t) = (1 - 3\sqrt{2})\, t\, p_{n-2}(d, t) < 0.$$
   (Remark: At $n = 4$, $D_4(d, t) = (d^4 - 3 d^2 t^2 + t^4) + (3 d^2 t^2 - 3 t^4) + 2 t^4 = d^4$ identically!) Since sphere volume is constant and $\det L_n(s)$ strictly decreases, packing density strictly increases on $[0, 1]$. $\square$

### Lemma 4 (Discrete Path Sobolev Inequality & Exact Sum-of-Squares Lower Bound)
For any $m = (m_0, m_1, \dots, m_{n-1}) \in \mathbb{Z}^n$, define the quadratic form:
$$Q_n(m) := m L_n(1) m^T = m(2 I_n + S_n) m^T = 2 m_0^2 + m_1^2 + m_{n-1}^2 + \sum_{k=1}^{n-2} (m_k - m_{k+1})^2.$$
Setting boundary padding $\tilde{m}_0 = \tilde{m}_n = 0$ and $\tilde{m}_k = m_k$ ($1 \le k \le n-1$), write $d_k = \tilde{m}_{k+1} - \tilde{m}_k \in \mathbb{Z}$ so $Q_n(m) = 2 m_0^2 + \sum_{k=0}^{n-1} d_k^2$.

**(4a) Integer Sobolev inequality:** Since $d_k \in \mathbb{Z} \implies d_k^2 \ge |d_k|$, for any $j \in \{1, \dots, n-1\}$,
$$\sum_{k=0}^{n-1} d_k^2 \ge \sum_{k=0}^{j-1} |d_k| + \sum_{k=j}^{n-1} |d_k| \ge |m_j - 0| + |0 - m_j| = 2 |m_j|.$$
Together with $2 m_0^2 \ge 2|m_0|$ (for $m_0 \in \mathbb{Z}$), this implies $Q_n(m) \ge 2 \max(m_0^2, \max_{j \ge 1} |m_j|)$. In particular, for every $v \in P_n - P_n = \{0, \pm 2e_i, \pm e_i \pm e_j\}$:
- If $v = \pm 2 e_i$: $|v \cdot m| = 2|m_i| \le Q_n(m)$.
- If $v = \pm e_i \pm e_j$ ($i \ne j$): $|v \cdot m| \le |m_i| + |m_j| \le 2\max(|m_i|, |m_j|) \le Q_n(m)$ (if $i, j \ge 1$) and $|m_0| + |m_j| \le 2 m_0^2 + 2|m_j| \le Q_n(m)$ (if $i = 0$).

Therefore, setting $w_1 = m L_n(1) + v$, we have for **every** $n \ge 3$, $m \in \mathbb{Z}^n$, $v \in P_n - P_n$:
$$D := w_1 \cdot m = Q_n(m) + v \cdot m \ge 0.$$

**(4b) Exact SOS inequality for $E(m, v) := (2m + v) \cdot w_1$:** Expanding $(2m + v) \cdot (2m + m C_n + v)$ using $m \cdot (m C_n) = Q_n(m) - 2|m|^2$ cancels $\pm 4|m|^2$ identically:
$$E(m, v) = 2 Q_n(m) + 4(m \cdot v) + v \cdot (m C_n) + |v|^2.$$
We prove $E(m, \pm 2 e_i) \ge 0$ for every $i \in \{0, \dots, n-1\}$ and $m \in \mathbb{Z}^n$. By symmetry $m \mapsto -m$, take $v = -2 e_i$ ($|v|^2 = 4$):
- **Case $i = 0$:** $(m C_n)_0 = m_1 + m_{n-1}$, so
  $$E(m, -2 e_0) = 4(m_0 - 1)^2 + 2 m_1(m_1 - 1) + 2 m_{n-1}(m_{n-1} - 1) + 2 \sum_{k=1}^{n-2} (m_k - m_{k+1})^2 \ge 0,$$
  since $x(x - 1) \ge 0$ for every $x \in \mathbb{Z}$.
- **Case $i = 1$ (and identically $i = n-1$ by path reversal):** $(m C_n)_1 = -m_0 - m_2$, and the right path tail obeys $\sum_{k=2}^{n-1} d_k^2 \ge |m_2|$. Hence
  $$E(m, -2 e_1) \ge \underbrace{2 m_0(2 m_0 + 1)}_{\ge 0 \text{ on } \mathbb{Z}} + 2 m_1^2 + 2(m_1 - m_2)^2 + 2(|m_2| + m_2) - 8 m_1 + 4.$$
  - If $m_2 \ge 0$: $2 m_1^2 + 2(m_1 - m_2)^2 - 8 m_1 + 4 m_2 + 4 = 2(m_1 - 1)^2 + 2(m_1 - m_2 - 1)^2 \ge 0$.
  - If $m_2 = -r \le -1$: $2 m_1^2 + 2(m_1 + r)^2 - 8 m_1 + 4 = 4(m_1 - 1 + r/2)^2 + r^2 + 4r > 0$.
- **Case $2 \le i \le n - 2$:** $(m C_n)_i = -m_{i-1} - m_{i+1}$. Splitting $Q_n(m)$ into $2 m_0^2$, the left tail $\sum_{k=0}^{i-2} d_k^2 \ge |m_{i-1}|$, the two edges $(m_i - m_{i-1})^2 + (m_i - m_{i+1})^2$, and the right tail $\sum_{k=i+1}^{n-1} d_k^2 \ge |m_{i+1}|$, and applying the exact identity $2(a - x)^2 + 2(|x| + x) - 4a + 2 = 2(a - x - 1)^2 + 2(|x| - x) \ge 0$ twice gives:
  $$E(m, -2 e_i) \ge 4 m_0^2 + 2(m_i - m_{i-1} - 1)^2 + 2(m_i - m_{i+1} - 1)^2 + 2(|m_{i-1}| - m_{i-1}) + 2(|m_{i+1}| - m_{i+1}) \ge 0.$$

Now for any $v = s_i e_i + s_j e_j$ ($i \ne j$, $|v|^2 = 2$), writing $v = \frac{1}{2}(2 s_i e_i) + \frac{1}{2}(2 s_j e_j)$ and averaging $E(m, 2 s_i e_i) \ge 0$ with $E(m, 2 s_j e_j) \ge 0$ (whose linear part in $v$ averages to $4(m \cdot v) + v \cdot (m C_n)$, while the constant term averages to $4 = |v|^2 + 2$) yields:
$$E(m, v) \ge -2 \qquad \forall m \in \mathbb{Z}^n,\; v \in P_n - P_n.$$
Finally, when $|m|^2 = 1$ (say $m = e_k$ by sign symmetry) and $v = s_i e_i + s_j e_j$, $Q_n(e_k) = 2$ and $e_k C_n$ has entry $0$ at $k$ and two $\pm 1$ entries off $k$:
- if $v_k = 0$: $|v \cdot (e_k C_n)| \le 2 \implies E(e_k, v) = 6 + v \cdot (e_k C_n) \ge 4$;
- if $v_k = +1$: $|v \cdot (e_k C_n)| \le 1 \implies E(e_k, v) \ge 6 + 4 - 1 = 9$;
- if $v_k = -1$: only one coordinate of $v$ lies off $k$, so $|v \cdot (e_k C_n)| \le 1 \implies E(e_k, v) \ge 6 - 4 - 1 = 1$.
Thus $E(m, v) \ge 0$ whenever $|m|^2 = 1$ as well! $\square$

### Lemma 5 (Universal Positivity of All Five Degree-4 Bernstein Coefficients, $\forall n \ge 3$)
Write the vector quadratic $Q(s) = (1 + s^2)(m L_n(s) + v)$ in the degree-2 Bernstein basis on $[0, 1]$:
$$Q(s) = P_0 (1 - s)^2 + 2 P_1 s(1 - s) + P_2 s^2, \qquad P_0 = w_0 = (2 + \sqrt{2}) m + v, \quad P_1 = w_1 + \sqrt{2} m, \quad P_2 = 2 w_1.$$
With $1 + s^2 = 1(1 - s)^2 + 2(1) s(1 - s) + 2 s^2$, the five degree-4 Bernstein coefficients of $f_{m,v}(s) = |Q(s)|^2 - 2(1 + s^2)^2 = \sum_{k=0}^4 b_k \binom{4}{k} s^k (1 - s)^{4-k}$ are exact polarized inner products (with $M = |m|^2$, $W_1 = |w_1|^2 \ge 2$, $D = w_1 \cdot m \ge 0$, and $w_0 \cdot w_1 = E(m, v) + \sqrt{2} D$):

1. **$b_4 = 4(W_1 - 2) \ge 0$** (by Lemma 2, $W_1 \ge 2$).
2. **$b_3 = 2(W_1 - 2) + 2\sqrt{2} D \ge 0$** (by Lemma 2 and Lemma 4a).
3. **$b_0 = |w_0|^2 - 2 \ge 0$:** if $m = 0$, $|v|^2 \ge 2$; if $|m| \ge 1$, $|w_0| \ge (2 + \sqrt{2})|m| - |v| \ge 2 + \sqrt{2} - 2 = \sqrt{2}$.
4. **$b_1 = P_0 \cdot P_1 - 2 = (E(m, v) + 2 M - 2) + \sqrt{2}(D + 2 M + v \cdot m) \ge 0$:**
   - The $\sqrt{2}$ coefficient satisfies $D \ge 0$ and $2 M + v \cdot m \ge 2|m|^2 - 2\|m\|_\infty \ge 0$ (since $x^2 \ge |x|$ on $\mathbb{Z}$).
   - For the rational part: if $m = 0$, $E(0, v) - 2 = |v|^2 - 2 \ge 0$; if $M = 1$, Lemma 4b gives $E \ge 0 \implies E + 2M - 2 = E \ge 0$; if $M \ge 2$, Lemma 4b gives $E \ge -2 \implies E + 2M - 2 \ge -2 + 4 - 2 = 0$.
5. **$b_2 = \frac{2}{3}(w_0 \cdot w_1 + |P_1|^2 - 4) = \frac{2}{3}\bigl((E(m, v) + 2 M - 2) + (W_1 - 2) + 3\sqrt{2} D\bigr) \ge 0$:**
   Every grouped term is already proved $\ge 0$ in (1)–(4)!

Since every Bernstein basis polynomial $\binom{4}{k} s^k (1 - s)^{4-k} \ge 0$ on $[0, 1]$, we conclude $f_{m,v}(s) \ge 0$, hence $|m L_n(s) + v|^2 \ge 2$ for **all** $n \ge 3$, **all** $(m, v) \ne (0, 0)$, and **all** $s \in [0, 1]$. **Q.E.D.** $\blacksquare$

---

## 3. Special $4\mathrm{D}$ Companion Result: Clifford-Nilpotent Optimal $D_4$ Fold & Exhaustive Type-I Census

In $n = 4$, the equation $2^n = 4n$ has its **unique** positive integer solution $n = 4$, permitting a **nilpotent** fold generator ($C^2 = 0$):
$$C_4^{\mathrm{nilp}} = \begin{pmatrix} 0 & -1 & -1 & 0 \\ -1 & 0 & 0 & -1 \\ 1 & 0 & 0 & 1 \\ 0 & 1 & 1 & 0 \end{pmatrix} = \begin{pmatrix} -J_0 & -I_2 \\ I_2 & J_0 \end{pmatrix}, \qquad J_0 = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}.$$
Its symmetric part $S_4 = \operatorname{diag}(-J_0, J_0)$ satisfies $S_4^2 = I_4$, its skew part $A_4 = \begin{pmatrix} 0 & -I_2 \\ I_2 & 0 \end{pmatrix}$ satisfies $A_4^2 = -I_4$, and $S_4 A_4 + A_4 S_4 = 0$, so:
- $(C_4^{\mathrm{nilp}})^2 = S_4^2 + A_4^2 + \{S_4, A_4\} = 0$ identically;
- $\det(d I_4 + t C_4^{\mathrm{nilp}}) = d^4$ identically (strict monotonicity trivial);
- $\frac{1}{2}(L_4^{\mathrm{nilp}}(s) + L_4^{\mathrm{nilp}}(s)^T) = d I_4 + t S_4 \succeq (d - t) I_4 \succeq I_4$ (uniform $n=3$-strength infinite-image bound $|m|^2 \le 11$, with $I_4 + S_4$ equal to $2$ times an orthogonal projection!);
- Exact zero-float $\mathbb{Z}[\sqrt{2}]$ Bernstein replay verifies all $21,944$ pairs ($109,720$ exact Bernstein coefficients): **0 negative**, $236$ zero, $32$ identically zero (exactly the $24$ intra-bars $+ 8$ directed cross-bars), and exact kissing number $24$ at all $8$ balls (Musin optimal $4\mathrm{D}$ $D_4$ sphere packing $\pi^2/16$).
- Exhaustive Type-I ($C_{ii} = 0$, $|C_i|^2 = 2$) census in $4\mathrm{D}$ ($12^4 = 20,736$ matrices): exactly $864$ valid $D_4$ endpoints, **all $864$ overlap-free**, $96$ with $I_4 + S \succeq 0$, and $48$ nilpotent ($C^2 = 0$, one $B_4$ symmetry orbit).

---

## 4. Scope, Prior Art & Falsification Criteria

- **Prior art:** The $n = 3$ case $C_3$ is `Foundation-0/incentre-lattices/results/fold-continuous-fcc` (GPT-6 / Aska Digital, 2026), which relied on a finite computer check of 3,248 quartics and $n = 3$-specific determinant derivative. The bridged cross-polytope framework topology for general $n$ is in Borcea & Streinu (2021, §7) and `incentre-lattices/results/rigidity-nd`. Generalizing the explicit continuous fold to all $n \ge 3$ ($D_n$, including optimal $4\mathrm{D}$ $D_4$ and $5\mathrm{D}$ $D_5$), the nilpotent $4\mathrm{D}$ matrix $C_4^{\mathrm{nilp}}$, and the computer-free analytical proof of all-pairs Bernstein nonnegativity via the discrete path Sobolev inequality $Q_n(m) \ge 2\|m\|_\infty$ are **new results** not present in `incentre-lattices`.
- **What is NOT claimed:** No physical phase-transition claim in real materials; no optimality claim for $D_n$ in $n \ge 6$ (where $E_6, E_7, E_8, \Lambda_{24}$ are denser than $D_n$).
- **Falsification criteria (for `@verifier`):** Any single $(n \ge 3, m \in \mathbb{Z}^n, v \in P_n - P_n)$ with $(m, v) \ne (0, 0)$ violating $Q_n(m) \ge 2\|m\|_\infty$, $E(m, \pm 2 e_i) \ge 0$, $b_k < 0$, $|m L_n(s) + v|^2 < 2$, $\det L_n(1) \ne 4n$, or $D_n'(t) \ge 0$ would refute the theorem.
