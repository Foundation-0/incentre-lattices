# Bounded motion gain with fixed angle springs

A fixed bending law preserves large internal motion near the open octahedral framework. It replaces the ideal divergence with a finite peak and a shrinking range of linear strain. The fully folded endpoint fails a spatial stability test, even though its uniform-strain stiffness is positive.

**VERIFIED (finite, sampled), NEW-COMPUTATION. Independent review pending.** This is an ideal point-joint design with distributed dead holding forces, as defined in [PROTOCOL.md](PROTOCOL.md). It is not a prediction for a measured material, a free specimen under boundary loading, or the hard-sphere packing at the fold endpoint.

## What stays fixed

The design has fifteen central springs per cell, all with rest length sqrt(2), and 24 cosine-angle springs, all preferring the open-state angle of 135 degrees. Spring constants and preferred angles remain unchanged through the fold. The angular stiffness parameter is kappa, in units set by the central stiffness and motif radius.

The folded configurations store angular energy. Applied forces balance that energy gradient at every joint and in the cell-strain coordinates. These forces remain constant during each local test. The calculation retains the angular stress contribution to the Hessian, the second derivative of energy that determines local stability.

## Response and strain range

The [derivation](DERIVATION.md) reduces the small-fold response to three octahedron rotations. Exact symbolic matrices give the asymptotic peak location and gain:

```text
s_peak = 1.2421026305*sqrt(kappa) + higher-order terms
peak ||du/de_xx|| = 0.4042260417/sqrt(kappa) + lower-order terms.
```

Here u is the eighteen-component internal displacement after removing translations, and e_xx is imposed strain. The norm is measured in motif-radius units. It sums motion across the six joints. It is not the motion of one output joint or an energy gain.

The full prestressed model evaluated at the asymptotic peak location gives:

| kappa | Fold s | Internal motion per unit xx strain |
|---:|---:|---:|
| 1e-8 | 0.0001242103 | 4042.29 |
| 1e-6 | 0.001242103 | 404.255 |
| 1e-4 | 0.01242103 | 40.4394 |
| 1e-2 | 0.1242103 | 3.95109 |

The table evaluates a derived location, rather than claiming exact maximization at finite kappa. The independent symbolic reduction and full numerical model agree increasingly well as kappa decreases. The 484-point parameter scan is saved in [evidence.json](evidence.json).

At kappa=1e-4 and s=0.01242103, the smallest internal eigenvalue is 0.000583484 and the smallest relaxed symmetric-strain modulus is 0.0000290924. The external joint-force norm is 0.0000130087 and the six strain-force norm is 0.00000479337, in the declared units. Positive values here do not by themselves establish spatial stability.

Full nonlinear cell relaxation gives these errors relative to the linear displacement vector:

| Magnitude of xx strain | Positive strain | Negative strain |
|---:|---:|---:|
| 1e-4 | 6.27% | 3.50% |
| 1e-5 | 0.506% | 0.478% |
| 1e-6 | 0.0493% | 0.0491% |

Near the peak, the small-strain condition scales as `|strain| << kappa`. Decreasing the angle stiffness raises gain and narrows this interval. The examples therefore support a bounded response, not unlimited useful amplification.

All eighteen nonlinear cases pass a force-residual and positive-curvature check. One raw trust-region solve reports failure to predict an energy improvement at floating-point resolution. A Newton correction of the force equation brings its residual displacement estimate below 2e-15. Both the original solver status and the correction are retained in the evidence.

## Spatial stability changes the verdict

The full 18-by-18 Bloch Hessian was sampled at 708 wavevectors per selected state: an 8-by-8-by-8 reciprocal grid, 100 seeded random points, and 96 small wavevectors. The latter can have eigenvalues near floating-point resolution, so a second calculation forms the acoustic matrix directly from exact phase derivatives at q=0.

That acoustic calculation eliminates internal modes from the full Bloch expansion. It does not infer wave stability from symmetric strain alone. It samples 149 directions at each state, including coordinate axes, small-integer directions and seeded random directions.

For kappa=1e-6, the folds s=0, 0.001, 0.003, 0.01, 0.1 and 0.5 have no resolved negative grid eigenvalue below -1e-10. Their sampled acoustic coefficients are positive, with the smallest between 7.09e-8 and 8.45e-8. The kappa=1e-4 operating point also passes these finite checks. Its smallest sampled acoustic coefficient is 8.34e-6.

A 12-by-12-by-12 refinement with 200 random points, 2024 wavevectors in total, preserves those verdicts at s=0.001 with kappa=1e-6 and at the kappa=1e-4 operating point. The [refinement output](refinement.json) retains the unresolved roundoff-scale values.

The endpoint s=1 fails both tests. At `q=(pi,pi,pi)` and kappa=1e-6, its smallest eigenvalue is -2.84111e-7. Across kappa from 1e-8 through 1e-4, that eigenvalue divided by kappa stays close to -0.28411. Eigenvector residuals are below 1.4e-15.

The unstable displacement alternates sign in neighbouring cells. Its six joint vectors are saved in the evidence. Direct energy changes in a 2-by-2-by-2 periodic supercell reproduce the negative curvature to 0.011% relative error at amplitude 1e-5. This control uses finite lengths and angles, rather than the Bloch matrix.

If the angular geometric-stiffness term is omitted, the same wavevector instead has minimum eigenvalue approximately +kappa/2. That omission would reverse this stability verdict. The stored angular stress matters even though the original central bars remain unstretched.

These checks do not certify every unsampled wavevector, every direction, global nonlinear stability, or a finite specimen with boundaries. The endpoint has no added contact constraints. Finite spheres would require a different model there.

## Evidence and prior work

[probe.py](probe.py) contains the analytic forces, prestressed Hessian, cell relaxation, Bloch expansion and finite supercell control. Its force-derivative checks have maximum absolute error 1.27e-10. [asymptotic.py](asymptotic.py) derives the reduced rotation matrices using exact SymPy arithmetic, without fitting the numerical curve.

[sources.json](sources.json) records source commits, file hashes, software versions and literature limits.

General rank-one elasticity and rational response regularized by restoring torque already appear in [Nassar, Chen and Huang (2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6283905/). Extensive modes in related linearly bridged frameworks are established by [Goodwin (2006)](https://doi.org/10.1103/PhysRevB.74.134302).

The latter comparison remains limited to the primary abstract. No historical priority is established for the present fold-specific formula.

From this directory, run:

```bash
python3 probe.py
python3 probe.py --check
python3 asymptotic.py
python3 refine.py
```

The first two commands require NumPy and SciPy. The symbolic command also requires SymPy. Exact JSON replay checks the same numerical environment. The scientific assertions use stated tolerances. The optional refinement reruns three selected states on a 12-by-12-by-12 grid with 200 random points.

The next mechanics gate is a boundary-controlled specimen with a specified hinge realization and tolerances. This calculation supplies its force law, expected crossover and a stability counterexample. It does not establish that such a specimen will retain the same gain under boundary loads alone.
