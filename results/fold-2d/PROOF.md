# The hinged-squares fold as a ball framework (2D)

Evidence label: **KNOWN linkage (Flores; Stuart). Corollary numbers PROVEN exactly (sympy). Uniqueness of the incentre among dot
positions VERIFIED on a scan only (finite check).** Replay: `make verify-fold2d` (seconds; numpy, sympy). Certificate:
`certificate.json`. This result is deliberately modest: the numbers are elementary corollaries of a known mechanism.

## Statement

Take squares of four dots at +-c e_i (i = 1, 2) in cells of side 2, joined by cross bars between neighbouring cells, with bars
between touching balls of radius R = sqrt2 - 1 and c = 2 - sqrt2 (the incentre of the 4.8.8 construction), so every bar has length
2R. Balls are joined centre to centre (they may roll on one another); a model with pins at the contact points behaves differently
and gives L = 2 cos(theta), which shortens the bars.

Rotate alternate squares by +-theta about their centres. The cross bars stay on the axes and the period in both directions is

    L(theta) = 2R + 2c cos(theta).

* The closest pair of non-bonded balls has gap 2 sqrt2 c cos(theta + pi/4), so the first new contact is at **theta = pi/4**.
* There L = 4R, the linear strain is L/2 - 1 = -R^2 = -(3 - 2 sqrt2), the dots form a square lattice of spacing 2R, and the packing
  fraction has risen from pi (sqrt2 - 1)^2 = 0.5390 to pi/4.
* The motion is isotropic along the path (equal period in both directions), the symmetric path of the mechanism. The periodic
  framework has other infinitesimal flexes, including a shear.
* Among dots on the axis mirror lines at +-t e_i, only t = 2 - sqrt2 ends in a square lattice with a single contact distance
  (at theta = pi/4 the dot spacings along each axis alternate between t sqrt2 and 2 - 2t, equal only at the incentre). For t < 2 - sqrt2
  the motion continues past pi/4 and stops at a state with three contacts per ball (numerical scan).

Scope: the plane, rolling contacts, the symmetric path. "Dots on the diagonal mirror lines" is a different family with no fold
at all (two bar constraints force theta = 0); do not conflate the two.

## Derivation

Gap along y between the dots of adjacent counter-rotated squares: L(theta) - 2 c sin(theta) - 2R = 2R + 2c cos(theta) - 2c sin(theta) -
2R = 2c (cos theta - sin theta) = 2 sqrt2 c cos(theta + pi/4). It vanishes first at theta = pi/4, and sympy confirms
L(pi/4) - 4R = 0 and the strain identity. That no other non-bonded pair comes closer for theta < pi/4 is checked numerically
(minimum distance 1.0003 * 2r at theta = 44.99 degrees in the internal review), not proved in general; the replay finds the lock
angle on a 0.01-degree grid. The statement about t is clear at theta = pi/4 for the spacing argument above; the case t < 2 - sqrt2
is numerical.

## What the replay adds

Exact simplification of the gap and the values at pi/4; a scan over c in {0.45, 0.5, 0.55, 2 - sqrt2, 0.62, 0.7} giving the lock
angle (66.37, 56.96, 49.51, 45.00, 45.00, 45.00 degrees) and the nearest-neighbour distances at the lock; must-fail controls (wrong
gap formula, 45 degrees claimed for c = 0.5, a non-incentre c claimed to end in a single-distance square lattice, a wrong end
density).

## Prior art

**KNOWN (the mechanism).** A. Flores, "Hinged tilings", North American GeoGebra Journal 6(1), section 2.9, Fig. 9, "hinged squares with
rods" (squares linked by rigid rods; the 4.8.8 tiling and the square tiling are its end states), crediting D. Stuart, "Polyhedral and
mosaic transformations", Student Publications of the School of Design, NC State 12(1), 2-38 (1963); the term "hinged tessellations" is
from D. Wells, Penguin Dictionary of Curious and Interesting Geometry (1991). Rotating squares as an auxetic mechanism: J. N. Grima
and K. E. Evans, J. Mater. Sci. Lett. 19, 1563 (2000). The equi-auxetic criterion: Mitschke et al., EPL 102, 66005 (2013) (this predicts
nu = -1 along the symmetric path). The flexibility-window notion: Sartbaeva, Wells, Treacy and Thorpe, Nature Mater. 5, 962 (2006); the
kagome analogue (triangles twisted into the triangular packing) is of the same type (Donev et al. 2004; Sun et al. 2012).
**PARTIALLY KNOWN (ball-packing numbers): trivial corollaries. NOT FOUND: the uniqueness of the incentre as an elementary statement.**

## Open obligations

* A proof that no non-bonded pair is closer than 2r before pi/4 (currently numerical).
* The case t < 2 - sqrt2 and the continuation after pi/4 are numerical.
