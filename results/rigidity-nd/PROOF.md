# Rigidity and elasticity of the n-dimensional contact framework

Evidence label: **PROVEN for every n >= 3 (written proof below); VERIFIED (finite check) for n = 2, ..., 8.**
Replay: `make verify-rigidity` (a few seconds; numpy and sympy). Certificate: `certificate.json`.

## Statement

Take the cube Q = [-1, 1]^n and cut it into the 2n pyramids over its facets. The pyramid over the facet x_i = 1 is
{x in Q : x_i >= |x_j| for all j}. Its facets are the base x_i = 1 and the walls x_i = +-x_j. A point t e_i is at
distance 1 - t from the base and t/sqrt2 from every wall, so the incentre is c e_i with

    c = 2 - sqrt2,   R = 1 - c = sqrt2 - 1,   c sqrt2 = 2R.

The dots are p(i, s) = s c e_i for i = 1..n and s = +-1, repeated over the lattice 2 Z^n. Two dots with different axes
in one cell are at distance c sqrt2 = 2R. The dot p(i,+) and the dot p(i,-) + 2 e_i of the next cell are at distance
2 - 2c = 2R. All other distances are larger.

The **contact framework** F_n has these dots as joints and a bar of length 2R for every touching pair. Per cell there
are 2n joints (2n^2 degrees of freedom) and 2n(n-1) + n = 2n^2 - n bars: the 2n(n-1) *intra* bars (the edges of the
cross-polytope) and n *cross* bars X_i joining neighbouring cells. Each ball touches 2n - 1 others.

**Theorem.** Let n >= 3.

* (a) The cross-polytope framework alone (intra bars of one cell) is infinitesimally rigid in R^n and has exactly
  n(n-3)/2 independent self-stresses. A basis: for any w = (w_ij), i < j, with sum over j != i of w_ij = 0 for every i,
  put tension w_ij on all four bars between the axes i and j.
* (b) At a wavevector q with no q_i in pi Z, F_n has exactly n(n-1)/2 zero modes per cell (the infinitesimal rotations of
  each cross-polytope, with the common Bloch phase) and exactly n(n-3)/2 self-stresses (those of (a), zero on cross bars).
  If exactly k of the q_i lie in pi Z there are n(n-1)/2 + k zero modes and n(n-3)/2 + k self-stresses. At q = 0 there are
  n(n+1)/2 zero modes and n(n-1)/2 self-stresses: for every w, tension w_ij on the four bars of the pair (i,j) and
  tension tau_i = sqrt2 (sum over j of w_ij) on X_i.
* (c) Give every bar unit spring constant and relax all joints under a homogeneous strain eps. The energy per unit volume is
  W(eps) = 2^-n [ (2(n-1)^2/n^2) sum_i eps_ii^2 + (2/n^2) sum_{i != j} eps_ii eps_jj ], independent of the off-diagonal
  strains. So C11 = 4(n-1)^2/(n^2 2^n), C12 = 4/(n^2 2^n), every shear modulus is 0, C11/C12 = (n-1)^2, and the Poisson
  ratio for uniaxial stress along a cube axis is **nu_100 = 1/(n^2 - n - 1)**: 1/5, 1/11, 1/19, 1/29, 1/41, 1/55 for n = 3..8.

Scope. Bars are ideal central-force springs, the framework is the infinite periodic one with one-cell periodic strain, and
"zero mode" means a first-order (infinitesimal) flex. Nothing here is a statement about finite motions or real materials.
The counting statements use only bar directions, so they hold for dots at +-t e_i for any 0 < t < 1. The elastic
formula with "identical bars" needs t = c (otherwise cross bars have length 2 - 2t and c + R = 1 is replaced by
t + (1 - t) = 1; the formula stays the same).

By rank-nullity, (zero modes) - (self-stresses) = 2n^2 - (2n^2 - n) = n at every q. This is the periodic Maxwell-Calladine
count, so each number below determines the other.

## Proof

### 1. Splitting the problem into blocks

Fix a cell and write the joint displacements as u(i,s). Set a_i = (u(i,+) + u(i,-))/2 and b_i = (u(i,+) - u(i,-))/2, and
write v^k for the k-th component of a vector v. For i < j define four scalars:

    alpha_ij = a_i^i - a_j^i,   beta_ij = a_i^j - a_j^j,   kappa_ij = b_i^i + b_j^j,   sigma_ij = b_i^j + b_j^i.

**Lemma 1.** Under a homogeneous strain eps plus displacements u, the first-order stretch of the intra bar from p(j,t) to
p(i,s) is

    delta = (1/sqrt2) [ c (eps_ii + eps_jj) + kappa_ij + s alpha_ij - t beta_ij - s t (2c eps_ij + sigma_ij) ].

For a one-cell-periodic displacement the stretch of the cross bar X_i is 2R eps_ii - 2 b_i^i. For a Bloch displacement
with eps = 0 it is e_i . (z_i u(i,-) - u(i,+)), where z_i = exp(2 i q_i).

*Proof.* The bar vector is l = c (s e_i - t e_j), with unit vector l/(c sqrt2), and delta = (l/|l|) . (eps l + u(i,s) - u(j,t)).
Substitute u(i,s) = a_i + s b_i and use s^2 = t^2 = 1. For X_i the bar vector is 2R e_i and the displacement difference is
-2 b_i in the periodic case. QED

The four functions 1, s, t, st on {+-1}^2 are orthogonal. So the four bars between the axes i and j impose exactly the
four independent scalar conditions kappa = alpha = beta = sigma = 0 when eps = 0, and

    (1/2) sum over (s,t) of delta^2 = (c(eps_ii + eps_jj) + kappa_ij)^2 + alpha_ij^2 + beta_ij^2 + (2c eps_ij + sigma_ij)^2.   (*)

The linear system has split into three independent pieces: one in the n numbers b_i^i (through kappa), one in the n^2
entries of the a_i (through alpha, beta), and one in the n(n-1) entries b_i^j with i != j (through sigma).

### 2. Part (a): rigidity of the cross-polytope and its stresses

* *The a-block.* alpha_ij = beta_ij = 0 for all i < j says that for each k all the numbers a_j^k (j = 1..n) equal a_k^k.
  So a_j = v for every j: an n-dimensional space of translations. The map a -> (alpha, beta) goes from R^(n^2) to
  R^(n(n-1)) with an n-dimensional kernel, so it is onto and its n(n-1) conditions are independent.
* *The sigma-block.* Each sigma_ij involves the two variables b_i^j and b_j^i, which occur nowhere else. The C(n,2)
  conditions are independent, and the solutions are the matrices B = (b_i^j) whose symmetric off-diagonal part vanishes.
* *The kappa-block.* The conditions b_i^i + b_j^j = 0 for all i < j are x_i + x_j = 0 along the edges of the complete graph
  K_n. For n >= 3 a triangle i, j, k gives x_i = -x_j = x_k = -x_i, so x = 0. The coefficient matrix is the unsigned
  incidence matrix of K_n, which has rank n; its left kernel (edge weights w with sum over j of w_ij = 0 at every vertex)
  has dimension C(n,2) - n = n(n-3)/2.

So the zero modes are a_j = v, b_i^i = 0 and B antisymmetric. That is u(i,s) = v + s B' e_i with B' antisymmetric, and with
Omega = B'/c these are the translations plus the rotations u(i,s) = v + Omega p(i,s). The kernel has dimension
n + n(n-1)/2, which is exactly the trivial motions of a rigid body, so the cross-polytope is infinitesimally rigid. The
rank is 2n^2 - n(n+1)/2 = n(3n-1)/2, so the self-stress space has dimension 2n(n-1) - n(3n-1)/2 = n(n-3)/2.

For the explicit basis, expand a stress on the four bars of the pair (i,j) in the characters 1, s, t, st. By Lemma 1 its
pairing with a displacement is a combination of kappa, alpha, beta, sigma with those coefficients. The a- and
sigma-blocks have no relations, so the s, t and st coefficients vanish, and the character-1 coefficients w_ij satisfy
sum over i<j of w_ij kappa_ij = sum_i b_i^i (sum over j of w_ij) = 0 for all b, that is sum over j of w_ij = 0.
Physically, equal tension on the four bars of the square spanned by p(i,+-) and p(j,+-) pulls p(i,s) inward along e_i
with force sqrt2 w_ij, and these pulls cancel at every joint. QED (a)

### 3. Part (b): counts at every wavevector

Intra bars join joints of one cell, so their rows carry no Bloch phase. A zero mode of F_n therefore restricts, on each
cell, to an infinitesimal motion of the cross-polytope, which by (a) is u(i,s) = v + Omega p(i,s) with complex v and Omega.
Because e_i . Omega e_i = 0, both e_i . u(i,+) and e_i . u(i,-) equal v^i, and Lemma 1 shows that the cross bar X_i forces
(z_i - 1) v^i = 0.

If z_i != 1 for all i then v = 0 and Omega is free: n(n-1)/2 zero modes. If exactly k of the z_i equal 1, the k matching v^i
are free: n(n-1)/2 + k zero modes. The self-stress counts follow from rank-nullity (self-stresses = zero modes - n).
The stresses of (a) remain stresses of F_n with zero cross-bar tension (intra rows carry no phase), and at generic q they
are all of them. At q = 0, take any w and the tensions in the statement. At p(i,s) the intra tensions give force
-sqrt2 s e_i (sum over j of w_ij) and the cross bar gives s tau_i e_i, so tau_i = sqrt2 (sum over j of w_ij) balances every
joint. These stresses form a space of dimension C(n,2) = n(n-1)/2, the full count. QED (b)

### 4. Part (c): the elastic constants

A homogeneous strain stretches every cell's bars identically. In a Bloch decomposition of the relaxation over any periodic
supercell, only the q = 0 component couples linearly to eps; the other components add non-negative energy and are optimally
zero. So it is enough to relax over one-cell-periodic displacements. By Lemma 1 and (*), the energy per cell is

    E = sum over i<j of [ (c(eps_ii + eps_jj) + kappa_ij)^2 + alpha_ij^2 + beta_ij^2 + (2c eps_ij + sigma_ij)^2 ] + 2 sum_i (R eps_ii - b_i^i)^2.

The a-variables appear only in alpha, beta; choose a = 0. The off-diagonal b_i^j appear only in sigma_ij, and each sigma_ij
can be chosen freely: take sigma_ij = -2c eps_ij. This removes every off-diagonal strain from the energy, so all shear
moduli vanish. (The relaxed motion keeps each cross-polytope undeformed; shear is taken up by first-order rotations.)

Put e_i = eps_ii, x_i = b_i^i and y_i = x_i + c e_i. Using c + R = 1 the remaining energy is

    E = sum over i<j of (y_i + y_j)^2 + 2 sum_i (e_i - y_i)^2.

Its Hessian in y is 2 (n I + J), with J the all-ones matrix, which is positive definite. Stationarity gives
n y_k + sum_j y_j = 2 e_k, hence sum_j y_j = (1/n) sum_j e_j and y_k = (1/n)(2 e_k - (1/n) sum_j e_j). Substituting,

    E_min = (2(n-2)/n) sum_i e_i^2 + (2/n^2) (sum_i e_i)^2 = (2(n-1)^2/n^2) sum_i e_i^2 + (2/n^2) sum_{i != j} e_i e_j.

Divide by the cell volume 2^n to get W. For uniaxial stress along e_1 the conditions dW/de_j = 0 (j >= 2) read
C12 e_1 + (C11 + (n-2) C12) e_perp = 0, so nu_100 = C12/(C11 + (n-2) C12). With C11 + (n-2) C12 = (4/(n^2 2^n)) ((n-1)^2 + n - 2)
this is 1/(n^2 - n - 1). QED (c)

**The case n = 2.** The proof of (c) does not use n >= 3, so W holds for n = 2: C11 = C12 = 1/4, all shear moduli vanish,
only the bulk modulus is non-zero, and nu_100 = 1 formally. The counts in (a) and (b) change: the kappa-block of K_2 is the
single condition b_1^1 + b_2^2 = 0, so the four-bar square has one internal flex, which couples to the translations through
the cross bars at generic q. Result: 2 zero modes and 0 self-stresses at generic q, 3 and 1 at q = 0. The replay confirms
exactly these numbers.

## What the replay adds (VERIFIED, n = 2..8)

`rigidity_nd.py` builds the Bloch rigidity matrix and the spring-network elastic tensor directly from coordinates, without
using the block decomposition. For each n = 2..8 (fixed seed) it checks: the counts at three random generic q, at q = 0 and
on the planes q_1 = ... = q_k = 0; that the cross-polytope is rigid with n(n-3)/2 stresses; that the explicit stress basis
spans the generic stresses; the explicit q = 0 stresses; that the generic zero modes span the cluster rotations; and the
elastic constants (C11, C12 to 1e-12, shear and normal-shear coupling below 1e-12, nu_100 to 1e-10) for c = 2 - sqrt2 and
c = 0.3, 0.75. `proof_steps.py` verifies Lemma 1, equation (*) and the minimisation symbolically (sympy) for n = 2..8.
Must-fail controls: the n >= 3 counts at n = 2, a wrong Poisson formula, a wrong C11/C12 ratio, deleting a cross bar
(which must add a zero mode), and a wrong stress-basis dimension are all rejected.

These are finite checks. They corroborate the proof for n <= 8; they do not replace it for larger n.

## Prior art

The source comparison, reading depths, elementary deductions and independent-review questions are recorded in
[LITERATURE.md](LITERATURE.md). The earlier family-level **NOT FOUND** claim is withdrawn.

* The 3D bridged-unit flexibility phenomenon is **KNOWN in substance**: A. L. Goodwin, "Rigid unit modes and intrinsic
  flexibility in linearly bridged framework structures", Phys. Rev. B 74, 134302 (2006), doi:10.1103/PhysRevB.74.134302.
  Abstract read; exact multiplicities not checked against the full paper. S. C. Power, "Polynomials for crystal frameworks
  and the rigid unit mode spectrum", Phil. Trans. R. Soc. A 372, 20120030 (2014), doi:10.1098/rsta.2012.0030, supplies
  the general symbol/phase-periodic machinery (Sections 2, 4, 5 read), not proof of priority for these exact formulas.
* **The bridged cross-polytope family is prior art**: C. S. Borcea and I. Streinu, "Periodic tilings and auxetic deployments",
  Math. Mech. Solids 26(2), 199–216 (2021), doi:10.1177/1081286520948116, Section 7. Its periodicity and spring-model
  comparison with this realization is detailed in LITERATURE.md: the same periodic graph at a different bar length (t = sqrt2 - 1 versus 2 - sqrt2); the counts here hold for any 0 < t < 1 and so cover it.
* **Closest physical study of the cab framework, KNOWN**: L. Li, K. Refson and M. T. Dove, "The contribution of phonons to the thermal expansion of some
  simple cubic hexaboride structures: SmB6, CaB6, SrB6 and BaB6", Phys. Chem. Chem. Phys. 25, 10749 (2023), doi:10.1039/d3cp01306e (abstract read, full
  text not). DFT lattice dynamics focused on rigid unit modes involving rotations of the B6 octahedra, as in Zn(CN)2, Prussian blue and Si(NCN)2; it reports
  very low flexibility of the connected octahedra. This studies four real compounds; its abstract does not establish the mode counts or elastic
  formulas for the ideal bar framework in n dimensions. I have not checked whether its full text counts the flat bands the way this result does.
* The periodic Maxwell-Calladine count: J. C. Maxwell, Phil. Mag. 27, 294 (1864); C. R. Calladine, Int. J. Solids Struct.
  14, 161 (1978).
* **Priority unresolved** for the exact closed forms in n (n(n-1)/2 and n(n-3)/2 per cell, C11/C12 = (n-1)^2,
  nu_100 = 1/(n^2-n-1)). The counts are elementary consequences of internal rigidity, transverse rotational bridge
  velocities and rank-nullity. The elastic formulas additionally require the stated finite spring stiffness on all bars.
  This bounded audit is not an exhaustive full-text search, and no novelty is claimed.
* The ingredients are classical; an independent specialist may judge the result routine. External human review remains
  pending, and no review request has been sent.

## Open obligations

* Independent human review of the written proof (internal pre-publication review only so far; see `docs/review-ledger.md`).
* A formal (Lean 4) check of the block identities for general n is an open track (`docs/formalization.md`).
* Physical meaning: the elastic tensor is only positive semi-definite (zero shear moduli); real materials add angular
  stiffness.
