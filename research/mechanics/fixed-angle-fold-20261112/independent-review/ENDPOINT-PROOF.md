# An exact negative-curvature direction at the fold endpoint

The fixed-angle point-joint model is locally unstable at the fully folded endpoint for every positive angular stiffness \(\kappa\), under the prescribed distributed dead loads. A single alternating rotation gives a negative second variation of its total potential. This strengthens the author's finite sampled endpoint result.

This independent addition uses the model in the reviewed protocol. It changes no spring, reference angle, contact or external-load assumption. It is not a statement about finite spheres with contact constraints or a specimen held only at its boundaries.

## A permitted alternating displacement

At \(s=1\), the six motif points are \(p_i\in\{\pm e_x,\pm e_y,\pm e_z\}\), and the lattice-row matrix is

\[
L=\begin{pmatrix}2&-1&-1\\1&2&-1\\1&-1&2\end{pmatrix}.
\]

In cell \(n\), choose the straight displacement path

\[
u_i(n;h)=(-1)^{n_x+n_y+n_z}h(\omega\times p_i).
\]

This is the Bloch wavevector \((\pi,\pi,\pi)\). It is periodic in an even supercell. The first-order changes of every central-spring length vanish: intra-octahedron edges undergo infinitesimal rigid rotation, while each bridge vector is unchanged to first order because it joins opposite motif points in opposite-sign cells.

The central springs therefore contribute zero quadratic energy along this path. They begin to contribute at fourth order in \(h\). The squared displacement norm per cell is \(4h^2|\omega|^2\).

## Direct angular-energy expansion

For each angle, let \(b\) be its endpoint bridge vector and \(d\) its edge vector. Both have length \(\sqrt2\). The bridge is unchanged along the chosen straight path. Write \(d(h)=d+hD\omega\), where \(D\omega=\omega\times d\). Then \(d\cdot D\omega=0\), so

\[
c(h)=\frac{b\cdot(d+hD\omega)}{2\sqrt{1+h^2|D\omega|^2/2}}
=c+\frac h2 b\cdot D\omega-\frac{h^2}4c|D\omega|^2+O(h^3).
\]

Substituting into the fixed energy \(\kappa(c(h)+1/\sqrt2)^2/2\), and summing the 24 angles, gives the rotation-space Hessian

\[
H_{\rm rot}=\kappa
\begin{pmatrix}
2&2&2\\2&2&-2+\sqrt2\\2&-2+\sqrt2&2
\end{pmatrix}.
\]

The external dead-load potential is linear in Cartesian displacements and has zero second derivative along this straight path. Its work also sums to zero over the alternating supercell. Thus this is the total-potential second variation for an admissible perturbation.

Choose \(\omega=(-\sqrt2,1,1)\). Dividing by the displacement norm gives the exact Rayleigh quotient

\[
\frac{\omega^TH_{\rm rot}\omega}{4|\omega|^2}
=\kappa\left(\frac14-\frac{3\sqrt2}{8}\right)
=-0.2803300859\,\kappa<0.
\]

By the variational characterization of the least eigenvalue, the complete Bloch Hessian has a negative eigenvalue for every \(\kappa>0\). Coupling to other displacements can only lower the minimum below this trial value. No assertion about the other wavevectors is needed to establish instability.

Minimizing within the three rotation coordinates gives

\[
\lambda_{\rm rot,min}/\kappa
=\frac{2+\sqrt2-\sqrt{38-4\sqrt2}}8
=-0.2841112370.
\]

This projected value matches the small-\(\kappa\) numerical coefficient. It is not asserted to be the exact least eigenvalue of the full eighteen-coordinate matrix at finite \(\kappa\).

## Independent checks

[independent_check.py](independent_check.py) constructs the geometry and energy without importing the author's model. It checks the central nullspace, derives the angular matrix symbolically and evaluates the finite energy at 70-digit precision. [independent-results.json](independent-results.json) records the source hashes and results.

At \(\kappa=10^{-6}\), finite second differences along the trial path have relative errors \(1.18\times10^{-3}\), \(1.18\times10^{-7}\), and \(1.18\times10^{-11}\) at amplitudes \(10^{-5}\), \(10^{-7}\), and \(10^{-9}\). Each curvature is negative. The convergence includes the central springs' positive fourth-order contribution.
