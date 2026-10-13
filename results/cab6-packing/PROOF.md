# Maximal hard-sphere packing fraction of the CaB6 structure type

Evidence label: **PROVEN (elementary upper-bound argument below); attainment and non-overlap VERIFIED (finite check: all pairs in a
5 x 5 x 5 block of cells, 15 size ratios).** Replay: `make verify-cab6` (a few seconds; numpy, scipy, sympy). Certificate:
`certificate.json`.

## Statement

CaB6 has Ca on the 1a site (0,0,0) and B on the 6f sites (x, 1/2, 1/2) and permutations (space group Pm-3m). Put the lattice
parameter a = 1. Hard spheres of radius rL sit on Ca and of radius rS on B; x and a are free; gamma = rS/rL. The packing fraction is
phi = (4 pi/3)(rL^3 + 6 rS^3).

**Theorem.** The maximal packing fraction of the CaB6 type at size ratio gamma is

    phi_max(gamma) = (pi/6)(1 + 6 gamma^3)                  for 0 < gamma <= sqrt2 - 1,
    phi_max(gamma) = (pi/6)(sqrt2 - 1)^3 (6 + gamma^-3)     for gamma >= sqrt2 - 1.

The global maximum is **phi* = (pi/6)(30 sqrt2 - 41) = 0.746865**, attained only at gamma = sqrt2 - 1, x = (sqrt2 - 1)/2 (equivalently
the B offset h = 1 - sqrt2/2 = 0.2929 from the octahedron centre), a = 2 rL. There the Ca spheres touch each other and each B sphere
touches five B spheres, but Ca and B do **not** touch: d(Ca,B) = 0.73681 a > rL + rS = 0.70711 a (gap 0.0297 a). Moreover
phi_max(gamma) beats fcc, pi/(3 sqrt2) = 0.740480, exactly when

    ((sqrt2 - 1)/6)^(1/3) = 0.41023 < gamma < (4 + 7 sqrt2)^(-1/3) = 0.41591,

a narrow, asymmetric window (-0.0040 and +0.0017 about sqrt2 - 1). At gamma = sqrt2 - 1 the NaCl type reaches pi(5 - 3 sqrt2)/3 = 0.7931,
so CaB6 is not the densest binary packing at this size ratio.

Scope. Model: two sphere sizes on the CaB6 positions, one free coordinate x and the lattice parameter; no other structure
type, no softness, no relaxation of the cubic symmetry.

## Proof

Set a = 1. Three of the non-overlap constraints are: Ca-Ca, 2 rL <= 1; B-B inside an octahedron, 2 rS <= sqrt2 (1/2 - x); B-B between
octahedra, 2 rS <= 2x. The last two give

    rS <= m(x) = min( (1/2 - x)/sqrt2, x ) <= m* = (sqrt2 - 1)/2,

with equality only at x* = (sqrt2 - 1)/2 (the two terms are equal exactly there, and one increases while the other decreases).

* From rL <= 1/2: phi <= (4 pi/3)(1/8)(1 + 6 gamma^3) = (pi/6)(1 + 6 gamma^3).
* From rS <= m*, writing phi = (4 pi/3) rS^3 (gamma^-3 + 6): phi <= (4 pi/3) m*^3 (6 + gamma^-3) = (pi/6)(sqrt2 - 1)^3 (6 + gamma^-3).

The first bound is smaller for gamma <= sqrt2 - 1, the second for gamma >= sqrt2 - 1 (they coincide at gamma = sqrt2 - 1).

*Attainment.* Both bounds are attained at x = x*: for gamma <= sqrt2 - 1 take rL = 1/2, rS = gamma/2 <= m*; for gamma >= sqrt2 - 1
take rS = m*, rL = m*/gamma <= 1/2. What remains is that no other pair of spheres overlaps. Here the Ca-B distance is
sqrt(x*^2 + 1/2) = 0.73681 and rL + rS <= 0.70711, so Ca-B is not even touching at the optimum, and Ca is held by Ca-Ca contact.
**This remaining step (every other pair in the infinite structure is at least as far apart as the listed constraints demand) is
checked here by exhaustive enumeration of all pairs in a 5 x 5 x 5 block of cells, not by a written argument.** It is a finite
check of an elementary geometric fact, and a short written argument would close it (open obligation below).

A proposed exact replacement for this finite-only step is given in [ATTAINMENT.md](ATTAINMENT.md), awaiting review.

*Value and window.* (4 pi/3)(1/8 + 6 rS^3) with rS = (sqrt2 - 1)/2 equals (pi/6)(30 sqrt2 - 41) exactly (sympy). The left window
end solves (pi/6)(1 + 6 gamma^3) = pi/(3 sqrt2), i.e. gamma^3 = (sqrt2 - 1)/6. At the right end, with gamma^-3 = 4 + 7 sqrt2,
(sqrt2 - 1)^3 = 5 sqrt2 - 7 and (5 sqrt2 - 7)(10 + 7 sqrt2) = sqrt2, so (pi/6) sqrt2 = pi/(3 sqrt2). Both verified symbolically. QED

**Corrections recorded during internal review** (they are the reason the statement above reads as it does): in the Wyckoff
convention B at (x, 1/2, 1/2), x = (sqrt2 - 1)/2 = 0.2071 is the distance from the Ca corner. As an offset h from the octahedron centre
the equal-bond value is h = 1 - sqrt2/2 = 0.2929 (at h = 0.2071 the intra- and inter-octahedron B-B distances are 0.293 and
0.586). And "Ca touches B" is false (above). Replay control 6 and control 5 reject these two misstatements.

## What the replay adds (VERIFIED)

* Exact (sympy): the closed form, the optimum, the gap, both window ends, the NaCl value.
* Brute force (numpy, scipy): for each gamma in a grid of 15 ratios (0.1 to 1.5, including sqrt2 - 1 and both sides of the window), maximise
  phi over h (bounded scalar minimisation on ten sub-intervals plus the exact optimum) with rL set by ALL pair distances in the
  5 x 5 x 5 block. The maximum difference from the closed form is below 1e-9.
* Must-fail controls: a wrong constant, the wrong branch, gamma = 0.40 and 0.42 claimed to beat fcc, "Ca touches B", and the wrong
  Wyckoff reading of x.

## Prior art

**Travesset 2017, not retrieved.** A. Travesset, "Soft Skyrmions, Spontaneous Valence and Selection Rules in Nanoparticle Superlattices", ACS Nano 11, 5375
(2017), doi:10.1021/acsnano.7b02219, is paywalled and I could not read it; the abstract (PubMed) has no CaB6 numbers. It is cited by S. Ren et al.,
arXiv:2006.01860, for superlattices sitting at the peak of the packing fraction and for the MgZn2 (sqrt(2/3)) and CsCl (sqrt3 - 1) peaks, but their
text quotes no CaB6 peak, and their Table 5 gives CaB6 only as ranges (experiment 0.43 to 0.47, binary LJ 0.3 to 0.5). So the **location of the maximum,
gamma = sqrt2 - 1, may already be known from Travesset 2017 (and, as below, can be read off Horst and Travesset 2016, Fig. 3); I could not confirm it.** The
closed form (pi/6)(30 sqrt2 - 41) is not confirmed as new either: it was not found in the sources that I could read. A library copy of Travesset 2017 is
the one check that would settle both.

**PARTIALLY KNOWN.** The curve phi_max(gamma) was plotted for this and 23 other lattices by N. Horst and A. Travesset, "Prediction of binary
nanoparticle superlattices from soft potentials", J. Chem. Phys. 144, 014502 (2016) (Fig. 3); the peak position can be read from
their figure, and the closed forms were not found in text. CaB6 is absent from the densest binary packings of A. B. Hopkins, F. H.
Stillinger and S. Torquato, Phys. Rev. E 85, 021130 (2012), and from the Murray-Sanders list of structures beating fcc (Phil. Mag. A 42,
721 (1980)), consistent with the narrow window. The equal-bond boron position is the Pauling-Weinbaum model (Z. Kristallogr. 87, 181
(1934)); the net is RCSR `cab`. CaB6-type nanoparticle superlattices are observed (Shevchenko et al., Nature 439, 55 (2006); author list
and title not verified).

## Open obligations

* Review the proposed exact all-pairs argument in [ATTAINMENT.md](ATTAINMENT.md). The existing evidence label remains unchanged pending review.
* Independent human review (internal pre-publication review only so far).
