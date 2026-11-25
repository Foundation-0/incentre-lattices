# Fixed-angle force law and small-fold response

The specified angle springs bound the ideal internal-motion singularity. The leading maximum grows as the inverse square root of their stiffness, while its location moves toward the open state. The derivation below applies to the distributed dead-load test in [PROTOCOL.md](PROTOCOL.md).

Status: written derivation with symbolic and numerical checks, pending independent review. The formulas are not a historical-priority claim.

## Forces and geometric stiffness

For one angle let `r=|b|`, `l=|d|`, `n=b/r`, `m=d/l`, and `c=n.m`. Its fixed preferred cosine is `c0=-1/sqrt(2)`. The derivatives with respect to the two bond vectors are

```text
g_b = (m-c*n)/r
g_d = (n-c*m)/l
H_bb = (-m*n.T-n*m.T-c*I+3*c*n*n.T)/r^2
H_dd = (-n*m.T-m*n.T-c*I+3*c*m*m.T)/l^2
H_bd = (I-n*n.T)(I-m*m.T)/(r*l)
H_db = H_bd.T.
```

If `g=(g_b,g_d)` and Hc is the block matrix above, the angle contributes

```text
grad U_angle = kappa*(c-c0)*g
Hess U_angle = kappa*[g*g.T + (c-c0)*Hc].
```

For a central spring, with `v/|v|=n` and `l0=sqrt(2)`,

```text
grad U_bar = (|v|-l0)*n
Hess U_bar = n*n.T + (1-l0/|v|)*(I-n*n.T).
```

Each bond vector depends linearly on joint displacements and on the six symmetric strain coordinates. Pulling these derivatives back through that linear map supplies forces and the full tangent matrix. The reference bars are unstretched. The angles are stressed at every nonzero fold.

Summing the 24 angles on the exact path gives

```text
U(s) = 12*kappa*s^2*(1+2*s^2)/(1+s^2)^2
U'(s) = 24*kappa*s*(1+3*s^2)/(1+s^2)^3.
```

The energy increases for s > 0. Thus no nonzero point on this path is an unloaded equilibrium of the chosen material. A dead load equal to each energy derivative establishes equilibrium. Its work is linear in the declared coordinates, so it adds no tangent stiffness. These loads must be changed when moving to a different reference fold, then kept fixed during its local test.

## Why the response has a finite peak

Let R0 be the fifteen-bar compatibility matrix at s=0, B0 its strain-extension matrix, and P the orthogonal projector onto the three-dimensional nullspace of R0.T. Let Q map a rotation vector to the six joint displacements `omega cross p_i`. Then `Q.T Q=4I` and `R0 Q=0`.

The twelve stiff internal deformations can be eliminated at leading order by applying P. Write R1 for the derivative of the compatibility matrix at zero fold. The reduced central energy uses `A=P R1 Q`. The angle Hessian on rotations at s=0 is G. Exact arithmetic in [asymptotic.py](asymptotic.py) gives

```text
A.T A = [[16/3, 8/3,   8/3],
         [8/3, 80/9, -56/9],
         [8/3,-56/9,  80/9]]
A.T P B0[:,xx] = (20/9)*(1+sqrt(2))*(0,-1,1)
G = 2*(2+sqrt(2))^2*I.
```

The leading rotation response to xx strain is therefore

```text
omega/e_xx = -s*(s^2*A.T*A+kappa*G)^(-1)*A.T*P*B0[:,xx].
```

The forcing vector is an eigenvector of `A.T A`, with eigenvalue 136/9. Defining `L0=2+sqrt(2)`, the divergent part of the internal displacement norm has the form

```text
||du/de_xx|| ~ (40*L0/9)*s / [(136/9)*s^2 + 2*L0^2*kappa].
```

The precise asymptotic statement keeps `x=s/sqrt(kappa)>0` fixed as kappa approaches zero. The product `sqrt(kappa)*||du/de_xx||` tends to the corresponding finite function of x. Stiff internal modes contribute bounded corrections. Angular geometric stiffness changes higher-order terms in this limit, and remains in every full calculation.

Differentiating the limiting function gives

```text
s_peak/sqrt(kappa) -> 3*L0/(2*sqrt(17)) = 1.2421026305365
sqrt(kappa)*gain_peak -> 5/(3*sqrt(17)) = 0.4042260417272.
```

For `sqrt(kappa) << s << 1`, the ideal `5*L0/(17*s)` response is recovered to leading order. At fixed positive kappa, taking s to zero removes the divergent rotational part. Other internal deformations remain finite.

This scaling comes with a shrinking range of linear strain. Rotation is of order `e/sqrt(kappa)` near the peak. Its quadratic geometric contribution becomes comparable to the imposed strain when e is of order kappa. Thus `|e| << kappa` is the asymptotic condition for trusting the linear gain. The finite-strain checks test this cost directly.

## Stability beyond uniform strain

At wavevector q, every neighbour displacement has phase `exp(i*q.n)`. Assemble the Hermitian Bloch matrix by using that phase in the vector maps and retaining the complete angular Hessian. The external dead loads add no matrix term. The phases are measured per lattice step, not in fixed laboratory inverse-length units.

Uniform-strain positivity does not decide stability of these modes. Near q=0, evaluate the phase derivatives directly:

```text
H(rho*d) = H0 + rho*H1 + rho^2*H2 + O(rho^3).
```

Let T contain the three normalized translations and Z span their orthogonal complement. The acoustic coefficient is

```text
A_acoustic(d) = T.H * [H2 - H1*Z*(Z.H*H0*Z)^(-1)*Z.H*H1] * T.
```

This expression follows from eliminating the internal modes in the full Bloch matrix. It retains rotations and prestress terms that a calculation from symmetric strain alone could miss. The code checks a finite set of directions. Positive values on that set do not prove positivity for every direction or every nonzero q.

At the fully folded endpoint the saved alternating-cell mode has negative curvature. Direct energy changes in a 2-by-2-by-2 periodic supercell reproduce it. Omitting angular geometric stiffness makes this same wavevector appear positive, so the prestress term changes the stability verdict.

The model has point joints and central springs. It does not add the nonbar contact constraints that finite balls would acquire at the packing endpoint. Its endpoint instability therefore concerns this specified mechanical design, not a hard-sphere packing.
