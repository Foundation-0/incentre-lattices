# Hole/ball ratio of the pyramid-incentre ("dot") construction equals the classical radius ratio

Evidence label: **PROVEN under hypothesis H1 or H2 (written proof below); VERIFIED (finite check) on 20 polytopes, with 4
counterexample controls.** Replay: `make verify-radius-ratio` (a few seconds; numpy, scipy). Certificate: `certificate.json`.

## Statement

Let P be a convex polytope in R^n with the origin O in its interior. For a facet F let Pi_F = conv(O, F), the pyramid over F.
For a ridge E = F n G let W_E be the hyperplane through O and E (the *wall*). Suppose each Pi_F has an *inscribed ball*, a ball
inside Pi_F tangent to every facet of Pi_F. Its centre d_F is the *dot* and its radius is r_F.

**Theorem (corrected hypotheses).** Assume one of

* **H1.** Every wall W_E is a mirror of P (a reflection symmetry of P), and Pi_F has an inscribed ball for one, hence every, facet.
* **H2.** P has an insphere centred at O, and the angle phi between the normals of two adjacent facets is the same for all ridges.

Then all dots have the same radius r and the same distance D from O; the balls of two facets sharing a ridge touch, at a point
of the common wall; no two balls overlap; and the ball of radius D - r centred at O touches every ball. Consequently

    hole/ball = (D - r)/r = 2D/e - 1,   e = 2r,

which is the classical radius-ratio limit of the polytope Q formed by the dots: anions of radius e/2 on a sphere of radius D,
touching their neighbours and a central cation. Under H1, P is facet-transitive and Q vertex-transitive. Under H2, Q is a scaled
polar dual of P, inscribed with all edges equal to 2r, and

    (D - r)/r = 1/sin(phi/2) - 1.

Conversely, every convex inscribed polytope with equal edges and its circumcentre inside arises this way.

**Why the hypotheses are what they are.** Facet-transitivity plus an insphere is NOT enough. For tetragonal trapezohedra polar to
square antiprisms of height 0.7 and 1.3 times the uniform height, every pyramid has an inscribed ball, yet the distances between
adjacent dots range over [1, 1.017] * 2r and [1, 1.006] * 2r, so some contacts open. (Replay controls 1 and 2.) The theorem
holds for the trapezohedron only at the single shape dual to the uniform antiprism, and only under H2. A triangular prism with
square sides is not facet-transitive: its pyramids have different inradii (0.1667 and 0.1340) and the theorem does not apply
(control 3). The cuboctahedron is not facet-transitive either (control 4).

## Proof

**Lemma (uniqueness).** If a convex polytope has a ball tangent to all its facet hyperplanes from inside, it is the unique
largest inscribed ball.

*Proof.* Write the facets as n_k . x <= h_k with unit n_k and let B(x, r) be the tangent ball, so n_k . x + r = h_k. Suppose
B(y, r') lies inside with r' >= r. Then n_k . (y - x) <= r - r' <= 0 for all k. The polytope is bounded, so there are
lambda_k > 0 with sum lambda_k n_k = 0; applying it forces n_k . (y - x) = 0 for all k and r' = r. The n_k span R^n, so y = x. QED

**H1.** Let E = F n G and let rho_E be the reflection in W_E. It is a symmetry of P fixing O and E pointwise, so it permutes the
two facets through E. It cannot fix F: the plane of F contains E, which lies in W_E, and does not contain O, so if it were
rho_E-invariant it would be perpendicular to W_E and F would lie on both sides of E inside its own plane. So rho_E(F) = G and
rho_E(Pi_F) = Pi_G. The image of an inscribed ball of Pi_F is an inscribed ball of Pi_G, and by the Lemma it is THE one:
rho_E(B_F) = B_G. The ball B_F is tangent to W_E at a point that rho_E fixes, so B_F and B_G touch there, and r_F = r_G and
|O d_F| = |O d_G|. The ridge graph of a polytope is connected, so r and D are constant, and the reflections generate a group acting
transitively on facets. The balls lie in pyramids with disjoint interiors, so they do not overlap. The remaining claims
follow, since the central ball of radius D - r touches each ball along the ray O d_F.

**H2.** Let varrho be the inradius, n_F the unit outer normal of F and T_F = varrho n_F the tangency point on F. Because O is
equidistant from the planes of F and G and lies on the inner side of both, the wall W_E is the internal bisector of the
dihedral angle at E, with normal proportional to n_F - n_G. The angle theta between n_F and W_E satisfies
sin(theta) = n_F . (n_F - n_G) / |n_F - n_G| = sin(phi/2). The point d = (varrho - r) n_F on the segment O T_F has distance r
from the base and (varrho - r) sin(phi/2) from each wall of Pi_F. With r = varrho sin(phi/2) / (1 + sin(phi/2)) these agree, so
d = d_F and D = varrho - r = r / sin(phi/2). The reflection in W_E exchanges n_F and n_G, hence d_F and d_G, whose distance is
2 D sin(phi/2) = 2r. The dots D n_F are D varrho times the vertices n_F/varrho of the polar dual, whose edges correspond to
ridges of P. Conversely, if Q is inscribed in a sphere about O with equal edges, its polar has all facets at one distance from
O, and the angle between the normals at the ridge dual to an edge qq' is the angle qOq', the same for every edge. QED

**Unreachable shapes.** The CN7 and CN9 polyhedra of the radius-ratio tables and the dodecahedral CN8 are not
vertex-transitive, so no polytope satisfying H1 produces them. Under H2, Q must be inscribed with all edges equal. For the
triangulated shapes (pentagonal bipyramid CN7, snub disphenoid CN8, tricapped trigonal prism CN9) all faces would be equilateral,
and by the classical enumeration of the eight convex deltahedra the only ones with these combinatorics are the Johnson solids
J13, J84, J51. None is inscribed: for unit edge their vertex orbits lie at distances 0.851 and 0.526, 0.930 and 0.677, 0.764 and
0.996 from the centre. So these limits are unreachable by either route. The capped CN7 shapes under H2 were not examined. (The
replay recomputes these orbit radii; the classification of deltahedra is cited, not replayed.)

## What the replay adds (VERIFIED)

`radius_ratio.py` finds each inscribed ball by linear programming (largest ball in the pyramid, scipy `linprog`) and reports
whether it touches every facet, the spread of r and D, the range of adjacent-dot distances over 2r, the minimum pair distance
and the hole/ball ratio. `verify.py` requires, for all 20 cases that satisfy H1 or H2 (Platonic solids, rhombic dodecahedron,
two Catalan solids, bipyramids of six heights over k = 3, 4, 6, the pentagonal icositetrahedron, the uniform trapezohedron and
four regular 4-polytopes), that every pyramid has an inscribed ball, r and D are constant to 1e-6, all adjacent dots are at
2r, no pair is closer, and hole/ball equals the classical value (and 1/sin(phi/2) - 1 for the H2 cases). The 4 counterexamples
above must violate that conclusion. These are finite checks at floating-point precision; they do not prove the theorem.

## Prior art

**NOT FOUND as a statement.** The ingredients are KNOWN: the Wythoff construction and Catalan duals, the classical
radius-ratio rule for coordination polyhedra, and polar duality; the proof is elementary. The result explains one structural
limitation of the radius-ratio rule (it applies naturally only to mirror-generated, equal-edged polyhedra); it does not address
the softness and size effects behind most of the rule's failures (see George et al., 2020, cited in the preprint). Absence of
evidence only.

## Open obligations

* Independent human review of the H1 and H2 proofs (internal pre-publication review only so far).
* Which polytopes satisfy H1 or H2 is not classified here beyond the listed families.
* The capped CN7 shapes under H2.
