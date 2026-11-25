# Independent physical review of the fixed-angle fold calculation

**Disposition: no blocking scientific defect found in the frozen calculation.** The bounded internal-motion asymptote and the prescribed dead-load interpretation are consistent. The endpoint instability admits an exact negative-curvature proof, supplied separately in [ENDPOINT-PROOF.md](ENDPOINT-PROOF.md).

Reviewed 12 November 2026 by gpt-forge using gpt-6-astra, independently of the author. The reviewed files are those listed in the author's `checkpoint.json`. Every SHA-256 was verified before the checks. Exact hashes and executable results are retained in [independent-results.json](independent-results.json).

## Checks performed independently

1. Rebuilt the exact bond vectors from the specified fold geometry and differentiated their unit directions symbolically at \(s=0\). Eliminating the central compatibility range reproduces the reported rotation matrix, forcing vector and angular restoring matrix.
2. Recomputed the Euclidean displacement normalization \(Q^TQ=4I\), the forced eigenvalue \(136/9\), and the extremum of the response. This gives \(s_{\rm peak}/\sqrt\kappa=3(2+\sqrt2)/(2\sqrt{17})\) and \(\sqrt\kappa\,g_{\rm peak}=5/(3\sqrt{17})\). No numerical curve was fitted.
3. Constructed an alternating rigid-octahedron rotation at \(q=(\pi,\pi,\pi)\). Direct angle-energy expansion, including stored angular stress, gives an exact negative Rayleigh quotient for every \(\kappa>0\).
4. Checked that quotient with direct finite energies at 70-digit precision, without calling the author's force, Hessian or Bloch routines. The smallest-step relative error is \(1.19\times10^{-11}\).
5. Inspected the full Bloch assembly and its long-wavelength Schur reduction. The phase factors consistently act on neighbour displacements in both vectors of each angle. The acoustic reduction retains the full prestressed Hessian and eliminates the internal modes of the same matrix.

The executable independent checks are in [independent_check.py](independent_check.py). They do not rerun or import the author's model code. Numerical claims outside these checks were inspected against the stored evidence rather than independently resimulated across the full grid.

## Load and stability interpretation

The material law is unchanged along the fold. Since its angular energy increases away from the open state, folded states need applied forces. The declared dead loads provide those forces and contribute no Hessian in the stated Cartesian displacement and symmetric-strain coordinates.

This is consistent with the calculation, but it is a specific controlled test. Distributed joint forces are essential to that interpretation. The result does not establish an unloaded material modulus or the response of a specimen held only at its boundary. Cell rotation is constrained by the test frame, as stated in the protocol.

The omission of angular geometric stiffness would change the physical model. The independent endpoint energy calculation confirms that this term can reverse the stability verdict. Positive internal and uniform-strain matrices alone would not establish Bloch stability.

The author correctly limits positive Bloch and acoustic results to the sampled points and directions. These finite checks leave global wavevector positivity and global nonlinear stability unproved.

The new endpoint proof establishes one instability for all positive \(\kappa\), under the original model assumptions. It does not extend to contact-stabilized spheres.

## Prior work and next action

Restoring-torque regularization and rational rotational response are established mechanisms. Full primary text for [Nassar, Chen and Huang (2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6283905/) was accessible during this review and supports that overlap. Its grounded rotational design is different from this three-dimensional fixed-angle fold.

The primary abstract of [Goodwin (2006)](https://doi.org/10.1103/PhysRevB.74.134302) supports extensive rigid-unit flexibility in linearly bridged frameworks. Full-text priority comparison remains incomplete. Neither generic soft-mode regularization nor historical priority for the exact constants is established by the present computation.

The current ideal-model question is answered sufficiently for publication: a specified fixed angular law bounds the response, narrows its linear range, and destabilizes the endpoint.

More sampling of the same ideal model has lower value than a boundary-controlled realization with stated hinge behavior and tolerances. That next calculation needs a concrete physical loading and hinge design. No additional arbitrary spring variant is justified by this review.
