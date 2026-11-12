# Fixed-angle springs along the octahedral fold

This calculation tests whether the large internal motion of an ideal folded framework survives a fixed bending law. It describes an ideal mechanical design with specified joints and holding loads. No material constants or fabrication tolerance are fitted.

Status: protocol fixed, 12 November 2026. Calculations and independent review are pending. AI assistance: gpt-forge using the requested gpt-6-astra model.

The source geometry is Foundation-0/incentre-lattices at `2d374fff2a59cc9331566733043a4bb25c1040bf`, verified as its current main commit on this date.

Its [continuous fold proof](https://github.com/Foundation-0/incentre-lattices/blob/2d374fff2a59cc9331566733043a4bb25c1040bf/results/fold-continuous-fcc/PROOF.md) supplies the path. The earlier [ideal-elasticity investigation](https://github.com/Foundation-0/gravity/pull/11) supplies the question being tested.

## Fixed design and reference

Use six joints p = (+x,-x,+y,-y,+z,-z), twelve edges joining different axes, and three bridges joining +axis to -axis in the adjacent cell. All fifteen central springs have unit stiffness and rest length sqrt(2). The motif radius is the length unit. With 0 <= s <= 1,

```text
C = [[0,-1,-1], [1,0,-1], [1,-1,0]]
a = sqrt(2)*(1-s*s)/(1+s*s)
t = 2*s/(1+s*s)
L(s) = (2+a)*I + t*C
```

Every joint has one outward bridge and four intra-octahedron edges. Add an angular interaction between its outward bridge vector b and each of those four edge vectors d. There are 24 interactions per cell, each counted once. The fixed energy is

```text
c(b,d) = b.d / (|b| |d|)
c0 = -1/sqrt(2)
U = (1/2) sum_bars (|v|-sqrt(2))^2
    + (kappa/2) sum_24_angles (c(b,d)-c0)^2.
```

The same graph, spring constants and preferred angle apply at every s. The angular law is invariant under rigid rotation. Near its preferred angle it is an ordinary torsional spring with angular stiffness kappa/2. Away from that angle it is the specified cosine law, not a claim about any manufactured hinge.

All original bars remain at rest along the exact geometric path. The angular energy on that path is `3*kappa*((sqrt(2)-a)^2+t^2)`. The open state is unloaded. A folded state generally is not an unloaded equilibrium.

## Holding loads and tangent response

At each selected s and kappa, evaluate the exact joint forces and the six cell forces conjugate to symmetric homogeneous strain. Apply opposite internal forces as external dead loads. Their directions and magnitudes remain fixed during a local perturbation. The total potential in joint displacement u and symmetric cell strain e is

```text
Pi(u,e) = U(u,e) - grad_u U(0,0).u - grad_e U(0,0).e.
```

This defines force balance at the chosen state without changing the material. Different points on the fold require different applied loads. It is a model for a controlled test with distributed joint loading, not a free-standing specimen held only at its boundary. Cell rotation is fixed by the test frame.

Dead loads have zero Hessian in these coordinates. The internal angular Hessian is `kappa*(grad(c) grad(c)^T + (c-c0)*Hess(c))`. The second term is the geometric stiffness from the stored angular stress. Dropping it changes the model. A follower force whose direction turns with a joint would require its own load derivative and is outside this protocol.

For a stable cell, eliminate the fifteen internal coordinates after removing translations. Compute the relaxed six-component tangent stiffness by the Schur complement. Report internal displacement per unit xx strain, its maximum over the scan, the holding load norms, and the smallest internal and relaxed eigenvalues. This stiffness belongs to the specified loaded test, rather than an unloaded material tensor.

## Checks and stop rules

1. Derive angular gradients and Hessians explicitly. Compare them with finite changes in energy and force at both open and folded states. Check rigid-rotation invariance and net internal force.
2. Scan kappa from 1e-8 to 1e-2 and s from zero through one, with extra points around s proportional to sqrt(kappa). Retain the geometric stiffness in every result.
3. At promising points, calculate the full 18-by-18 Bloch Hessian on a reciprocal grid, plus small wavevectors and random directions. Remove only the three exact translations at q=0. A sampled negative eigenvalue rejects stability at that point. Positive sampled values are a finite check, not proof over every wavevector.
4. At any stable enhanced-response point, minimize the full nonlinear potential under small positive and negative imposed xx strain. Compare the resulting joint motion and energy with the tangent prediction, with smaller strain and tighter tolerances as controls.
5. Distinguish a broad interval of bounded enhancement from a near-instability divergence. If no useful stable regime remains, report that negative result. A second interaction model requires a physical reason.

No historical priority is claimed. General rank-one elasticity and stiffness regularization by restoring torques occur in [Nassar, Chen and Huang (2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6283905/), including equations 2.11 and 2.12. The comparison will target this fold and this fixed angle law, rather than the general existence of rank-one media.
