# The cubic-restricted [111] fold of the 3D contact framework

Evidence label: **Exact values at 60 degrees PROVEN (sympy identity); first contact at 60 degrees, coordination sequence and
trigonal scan VERIFIED (finite check). Not a proof of any path.** Replay: `make verify-fold3d` (slow, about 1 minute; numpy,
scipy, sympy). Certificate: `certificate.json`. This folder is an addition because the numbers belong to the same fold remarks
as `fold-2d`.

## Statement

Take the 3D contact framework (octahedra of six dots at +-c e_i, c = 2 - sqrt2, in cells of side 2, joined by bars of length
2R, R = sqrt2 - 1). Rotate the octahedra alternately by +-theta about [111] (the analogue of the a-a-a- tilt) and **keep the lattice
cubic**. The bars then fix the lattice parameter

    a(theta) = K1 + sqrt(4R^2 - 2 K2^2),   K2 = (2/3) c (1 - cos theta),   K1 = 2c cos theta + K2.

* At theta = 60 degrees, a = 8c/3 exactly, every ball has exactly six contacts at 2R (the next distance is 1.3744 * 2R), and the
  packing fraction is 27 pi sqrt2/256 = 0.4686, up from pi R^3 = 0.2233 at theta = 0. Numerically no earlier non-bonded contact occurs
  (first contact at 60.0 degrees on a 3000-point grid; bar-length deviation along the path below 1e-12).
* This is **not a jammed state.** If the lattice may become trigonal, with rows (a,b,b), (b,a,b), (b,b,a), every bar length is
  preserved, and on a grid of 1 degree by 0.01 the set of non-overlapping states connected (by grid adjacency) to theta = 0 contains
  states of packing fraction up to 0.648 (near theta = 31 degrees, with large shear). So 60 degrees / 0.4686 is the lock of the
  *cubic-restricted* path only.
* The six-contact net at 60 degrees has coordination sequence 6, 17, 42, 68, 110, 153, 212, 272, 348, 425, which matches RCSR `sms`
  (sphere-packing type 6/3/h2) in its first ten terms. We checked the coordination sequence only, not an isomorphism.

Scope. Cubic-restricted path (an extra constraint); grid connectivity is not a continuous path; the 60-degree lock angle is not
unique to the incentre but the six-contacts-per-ball closed form is. For c > 2 - sqrt2 the same angle appears with only two
contacts, for c < 2 - sqrt2 the lock comes past 60 degrees with five contacts (internal review, scanned 0.45 to 0.8).

## Derivation of the exact values

At theta = pi/3 the lattice parameter simplifies to a = 8c/3 (sympy: a(60) - 8c/3 = 0 exactly). Each of the six dots of a cell has
six contacts whose squared distance simplifies exactly to 4R^2 (sympy, sqrtdenest). The density is 6 (4 pi R^3/3)/a^3, which
simplifies to 27 pi sqrt2/256 (difference below 1e-25 at 30 digits). Exactly: at a = 8c/3 the density is 27 pi R^3/(64 c^3), and
R/c = 1/sqrt2, which gives 27 pi sqrt2/256 (this one-line reduction was added in an internal audit, pending review; the replay
checks the density numerically, as above).

## What the replay adds

The cubic path scan, the exact 60-degree algebra, the coordination sequence against the RCSR `sms` sequence (hardcoded
in `verify.py`; checked against the RCSR data file `3dall.txt`, which is not vendored in this repo), and the trigonal grid scan with its connected component. Must-fail controls: the 60-degree state claimed
to be a jam; a shifted lattice constant claimed to keep six contacts; the coordination sequence claimed to be fcc's.

## Prior art

**PARTIALLY KNOWN.** The framework is the Prussian-blue topology (rigid octahedra joined by linear bridges; Goodwin, Phys. Rev. B 74,
134302 (2006)); octahedral tilts in Prussian-blue analogues are an established topic and the a-a-a- tilt is Glazer's notation
(Acta Cryst. B 28, 3384 (1972)). The end net matches RCSR `sms` (R-3m, 6/3/h2), a catalogued packing type; its appearance as the lock
state of this mechanism, the 60-degree angle and the density 0.4686 were NOT FOUND. Before claiming the density, the density of 6/3/h2 in
Fischer's trigonal sphere-packing tables should be checked.

## Open obligations

* Identify the 60-degree net with `sms` by an isomorphism, not by ten terms.
* A continuous (not grid-adjacent) path for the trigonal states; superseded in practice by `results/fold-close-packing`.
