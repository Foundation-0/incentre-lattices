# Universal rigidity: the cross-polytope and the periodic framework

Evidence label: **PROVEN** (written proofs below) and **VERIFIED** (finite exact check, `universal_rigidity_check.py`, n = 2..7 for the finite parts and n = 2..6 with four rational c for the periodic part). Known results are marked **KNOWN** with a citation.
Replay: see [README.md](README.md). Notation follows [results/rigidity-nd/PROOF.md](../rigidity-nd/PROOF.md).

## Terms

A bar framework (G, p) has V vertices p_1, ..., p_V in R^d and a set of bars. A **self-stress** is a number omega_ij on each bar with sum over j of omega_ij (p_i - p_j) = 0 at every vertex. The **stress matrix** Omega has Omega_ij = -omega_ij off the diagonal and zero row sums.

A framework is **universally rigid** if every configuration in any R^D with the same bar lengths is congruent to it. It is **super stable** if some positive semidefinite (PSD) stress has rank V - d - 1 and the bar directions do not lie on a conic at infinity.

Connelly proved that a super stable framework is universally rigid (KNOWN: Connelly 1982, "Rigidity and energy", Inventiones Mathematicae 66). Connelly, Gortler and Theran (2017, "Affine rigidity and conics at infinity", IMRN) give the conic condition.

## 1. Maxwell's identity

**Proposition 1 (KNOWN: Maxwell 1870, "On reciprocal figures, frames, and diagrams of forces").** For every self-stress, sum over bars of omega_ij |p_i - p_j|^2 = 0. So trace Omega = 2 sum of omega_ij.

*Proof.* Take the dot product of the equilibrium equation at vertex i with p_i and sum over i. A bar ij contributes omega_ij p_i . (p_i - p_j) at i and omega_ij p_j . (p_j - p_i) at j. These add to omega_ij |p_i - p_j|^2. The total is zero. The diagonal of Omega is the sum of omega_ij over the bars at i, so trace Omega counts each bar twice. QED

## 2. No certificate for equal bars

**Proposition 2.** Let every bar with nonzero stress have the same length, and let V >= d + 2. Then trace Omega = 0 and the only PSD stress is zero. So the framework has no super-stability certificate.

*Proof.* By Proposition 1, l^2 times the sum of omega_ij is zero, so the sum is zero and trace Omega = 0. A PSD matrix has non-negative eigenvalues. If their sum is zero, all are zero, so Omega = 0. With no loops and no repeated bars, the off-diagonal entries give omega = 0. A certificate needs a PSD stress of rank V - d - 1 >= 1. That is not zero. QED

The hypothesis V >= d + 2 matters. A regular simplex has V = d + 1. The zero stress has the required rank 0, and the simplex is universally rigid. The script's regular tetrahedron control shows this.

## 3. Equal bars are never universally rigid

**Proposition 3.** A finite framework with all bars of length l > 0 and V >= d + 2 is not universally rigid.

*Proof.* Send vertex k to (l / sqrt 2) e_k in R^V. Every pair of images is at distance l, so every bar keeps its length.

Suppose the image were congruent to p. Then all pairs of vertices of p would be at distance l. Translate one vertex to the origin. The other V - 1 vectors v_k have |v_k|^2 = l^2 and v_j . v_k = l^2 / 2. Their Gram matrix is (l^2 / 2)(I + J), which is positive definite of size V - 1. A Gram matrix of vectors in R^d has rank at most d, so V - 1 <= d. This contradicts V >= d + 2. QED

This needs no stress. It rules out every stress-based route, including higher-order ones. The idea of a regular-simplex realisation is folklore (KNOWN, in the same spirit as Pak's colour-class simplex embeddings).

**Corollary (finite cross-polytope).** The cross-polytope with vertices +-e_i in R^n has V = 2n >= n + 2 exactly when n >= 2. All its 2n(n-1) bars have length sqrt 2. So it is not universally rigid for any n >= 2. The map sends vertex k to the unit vector e_k of R^(2n). The antipodal squared distance changes from 4 to 2.

## 4. Stresses of the finite cross-polytope

**Proposition 4.** The self-stresses of the cross-polytope are exactly the following. Choose weights w_ij for i < j with sum over j of w_ij = 0 for every i. Put tension w_ij on all four bars between axes i and j. The space has dimension max(0, n(n-3)/2).

The stress matrix is Omega = -(W kron J_2), where W is the symmetric weight matrix with zero diagonal and J_2 is the 2 by 2 all-ones matrix. Its trace is 0 and trace Omega^2 = 8 times the sum over i < j of w_ij^2. Every nonzero stress is indefinite.

*Proof.* For n >= 3 the description and the dimension are part (a) of `results/rigidity-nd/PROOF.md`. For n = 2 the framework is a four-cycle and has no stress. The script confirms the dimension by exact rank for n = 2..7 and checks that the family is independent.

The vertices (i, +) and (i, -) of one axis share no bar, so the diagonal block of Omega for axis i is 2 r_i I with r_i = sum of w_ij over j, which is 0. Between axes i and j every entry is -w_ij. This is -(W kron J_2). Its trace is 0. The trace of Omega^2 is tr(W^2) tr(J_2^2) = 2 (sum of w_ij^2 over i < j) times 4.

A nonzero stress has trace Omega^2 > 0 and trace Omega = 0. A PSD matrix with trace 0 is zero, so Omega is not PSD. For an explicit witness, take a pair with w_ij nonzero. The vector with 1 on the four vertices of axes i and j gives the value -8 w_ij. The vector with 1 on axis i and -1 on axis j gives +8 w_ij. QED

The last claim of Proposition 4 also follows from Proposition 2. All bars have one length and V = 2n >= n + 2, so neither omega nor -omega is PSD unless omega = 0. The explicit form is what the script checks.

## 5. The periodic framework with the lattice fixed

Take the periodic framework of `results/rigidity-nd/PROOF.md`. The dots are p(i, s) = s c e_i for i = 1..n and s = +-1, with the lattice 2 Z^n. Here 0 < c < 1 is any value. The incentre uses c = 2 - sqrt 2.

Bars are the 2n(n-1) intra bars between dots of different axes in one cell, and n cross bars joining p(i, +) to p(i, -) + 2 e_i of the next cell. A bar e from dot a to dot b with lattice shift L has vector l_e = p_b + L - p_a.

A **stress at zero wavevector** gives each bar a tension omega_e so that the pulls omega_e l_e on the two ends of every bar balance at every dot. The matrix Omega(0) has x^T Omega(0) x = sum over bars of omega_e (x_b - x_a)^2.

Here omega_e is force per unit length. The rigidity-nd proof uses the force itself, which is omega_e |l_e|. At c = 2 - sqrt 2 every bar has length 2R, so the two scales differ by the factor 2R. The cross value c (n-1) / (1-c) below is then sqrt 2 (n-1), as in part (b) there with w = 1.

**Proposition 5.** Let w_ij = 1 for all i < j and r_i = n - 1. Put tension 1 on every intra bar and tension tau = c (n-1) / (1-c) on every cross bar. This is a stress at zero wavevector. For x in R^(2n) put m_i = (x(i,+) + x(i,-)) / 2 and d_i = (x(i,+) - x(i,-)) / 2. Then

    x^T Omega(0) x = 4 [ sum over i<j of (m_i - m_j)^2 + (n-1) sum over i of d_i^2 / (1-c) ].

So Omega(0) is PSD of rank 2n - 1, and its kernel is exactly the constant vectors.

*Proof.* Balance at p(i, +). The intra bars to the dots (j, t) pull with the sum over j and t of w_ij c (t e_j - e_i). The t e_j terms cancel in pairs. The total is -2 c r_i e_i. The cross bar pulls with tau (2 - 2c) e_i. Balance gives tau = c r_i / (1 - c). At p(i, -) the signs reverse and the same value balances.

Write x(i, +-) = m_i +- d_i. For one pair of axes, the sum over s and t of (x(i,s) - x(j,t))^2 equals 4 (m_i - m_j)^2 + 4 d_i^2 + 4 d_j^2. Summing with weight w_ij over i < j gives 4 times the sum of w_ij (m_i - m_j)^2 plus 4 times the sum of r_i d_i^2. A cross bar contributes tau_i (2 d_i)^2 = 4 c r_i d_i^2 / (1-c). Since r_i (1 + c/(1-c)) = r_i / (1-c), the formula follows.

With w = 1 the first sum is the quadratic form of the Laplacian of the complete graph K_n. It is zero only when all m_i are equal. The second sum is zero only when all d_i are zero, because 1 - c > 0 and n - 1 > 0. So the form is non-negative and vanishes exactly when x(i, +) = x(i, -) = m for every i, which means x is constant. QED

**Proposition 6.** For every n >= 2 and 0 < c < 1, every configuration of the dots in any R^D, periodic with the same lattice 2 Z^n (embedded in R^D) and with the same length for every bar, is a translate of p.

*Proof.* Let the new dot positions be p(a) + u(a) with u(a) in R^D. Pad p and the lattice with zeros. For the stress of Proposition 5 compute the stress energy E = sum over bars of omega_e |l_e|^2. After the move each bar vector is l_e + u(b) - u(a). So

    E(p + u) = E(p) + 2 sum over bars of omega_e l_e . (u(b) - u(a)) + sum over bars of omega_e |u(b) - u(a)|^2.

The middle term is 2 times the sum over dots a of u(a) . (net pull on a), and this is zero by balance. The last term is the sum over the D coordinates k of u_k^T Omega(0) u_k, where u_k is the k-th coordinate of u. If every bar keeps its length, then E(p + u) = E(p). So the last term is zero. By Proposition 5, each u_k is constant. So u is a constant vector. QED

**Stated hypotheses of Theorem 5.8 in Dewar, Schulze, Tanigawa and Theran** (arXiv:2510.19529, "Uniquely realizable crystalline structures"). This is KNOWN. An independent check read the statement as follows: the lattice is fixed, there is an equilibrium stress of the periodic framework, its stress matrix at zero wavevector is positive semidefinite with a one-dimensional kernel, and no conic condition is needed.

The conclusion is universal rigidity among periodic configurations with that lattice.

What was checked, in exact arithmetic, for this framework: equilibrium, positive semidefiniteness, and kernel exactly the constants. The script does this for n = 2..6 and c = 1/2, 1/3, 2/5, 3/4. Proposition 5 does it for all n >= 2 and all 0 < c < 1 by hand.

What was taken from the paper by that check and not re-read for this note: the exact wording of the theorem, any standing assumptions on the framework (for example on the bar set or on regularity), and the paper's precise definition of universal rigidity. Proposition 6 does not depend on any of these. It proves the conclusion for periodic configurations with the same lattice directly.

**Zero modes.** The framework still has n(n-1)/2 nontrivial zero modes at zero wavevector, besides the n translations (`results/rigidity-nd/PROOF.md`, part (b)). The script confirms the count n(n+1)/2 by exact rank for n = 2..6. These modes are first-order flexes. Proposition 6 shows none of them extends to a motion that keeps the lattice and the bar lengths.

## 6. A remark on a free lattice (not a theorem)

Let the lattice be free to deform. The balance condition then has an extra term. As an independent check reads it from the literature (Dewar et al., equation 2), the matrix S = sum over bars of omega_e l_e kron l_e must vanish. This condition is taken from that reading and not re-derived here.

For the stresses of Proposition 5 with general weights w_ij, a direct sum gives S = diag(4 c r_i). The script checks this with exact arithmetic. So S = 0 forces every r_i = 0, and then every cross-bar tension is zero. The remaining stressed bars are intra bars of one length c sqrt 2, and Proposition 2 applies to them.

This says the stresses above cannot give a PSD certificate with a free lattice. It does not show that the framework fails to be universally rigid with a free lattice, and it does not show that it is rigid. That case is not shown.

## 7. What is not shown

* The free-lattice case, as stated in section 6.
* Whether the stress matrix of Proposition 5 is positive semidefinite at every wavevector. If it is, with kernel only at zero wavevector on the constants, the infinite framework resists every bounded displacement. This is the open question.
* Competitors periodic only on a larger cell, and non-periodic competitors.
* A computation at the irrational value c = 2 - sqrt 2. The proofs hold for it.

## Relation to the fold drafts

The open fold drafts (pull requests 106 and 107) study bar-preserving motions in which the lattice changes. This note concerns the fixed lattice and finds rigidity there. The two statements do not conflict.

## Novelty

Default N0, meaning known. Maxwell 1870, Connelly 1982, Connelly, Gortler and Theran 2017, and Dewar, Schulze, Tanigawa and Theran (2025) cover the mathematics.

Of these, only the last was read for this note, in part: its Theorem 5.8 and its equation 2. The others are cited as standard references.

The bridged cross-polytope family is prior art (see `results/rigidity-nd/LITERATURE.md`). A bounded search by mechanism did not find the incentre instance of Propositions 5 and 6 written down. That is not a priority claim.
