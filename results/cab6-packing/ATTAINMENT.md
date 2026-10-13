# Exact all-pairs distances and attainment for the CaB6 packing bound

Status: **proposed written proof, submitted for review**. This note closes the all-pairs obligation identified in
[PROOF.md](PROOF.md) by minimizing over the entire integer lattice. The existing result labels and numerical certificates
remain unchanged pending review. This is an elementary completion of the project's argument, with no claim of literature novelty.

## Model and exact distance lemma

Normalize the cubic period to 1. Write b = (1/2, 1/2, 1/2), and let e_1, e_2, e_3 be the coordinate unit vectors. The two sets of centres are

    Ca = Z^3,
    B  = { z + b + s h e_i : z in Z^3, i in {1,2,3}, s in {-1,+1} },

where 0 < h < 1/2 and the Wyckoff coordinate is x = 1/2 - h. The representatives 0 < x < 1/2 cover the nondegenerate family:
x and 1 - x describe the same six-site set, while the boundary cases have coincident sites and cannot accommodate positive radii.

**Lemma.** The minimum distances between distinct centres in the infinite periodic structure are

    d_CaCa = 1,
    d_CaB  = sqrt((1/2 - h)^2 + 1/2),
    d_BB   = min(sqrt(2) h, 1 - 2h).

These are minima over all lattice translations, not just a finite block of cells.

**Proof.** For a vector v, the squared distance to Z^3 is the sum of the squared distances of its coordinates to Z:
each integer coordinate can be chosen independently. A nonzero integer vector has length at least 1, proving the Ca-Ca case.

For Ca-B, two coordinates are congruent to 1/2 modulo Z, while the remaining coordinate is congruent to 1/2 + s h.
Their respective minimum absolute values are 1/2, 1/2 and 1/2 - h. This proves the stated Ca-B minimum, which is attained.

Every B-B displacement is z + s h e_i - t h e_j, for an arbitrary z in Z^3. There are three exhaustive cases:

1. **Same axis and same sign:** i = j, s = t. Distinctness excludes z = 0, so the minimum distance is 1.
2. **Same axis and opposite signs:** i = j, s = -t. The affected coordinate is congruent to 2h or -2h.
   Since 0 < 2h < 1, its distance to Z is min(2h, 1 - 2h); the other two coordinates can be zero.
3. **Different axes:** i != j. The two affected coordinates are congruent to h or -h.
   Since h < 1/2, each has minimum absolute value h, so the minimum distance is sqrt(2) h.

Each minimum is attained by an allowed pair. Taking the minimum of the three cases gives

    min(1, 2h, 1 - 2h, sqrt(2) h) = min(sqrt(2) h, 1 - 2h),

because sqrt(2) h < 1 and sqrt(2) h < 2h. This exhausts every pair and every integer translation. QED.

In particular, for any fixed h and gamma = rS/rL > 0, the exact largest admissible radius is

    rL_max(h, gamma)
      = min(1/2,
            sqrt((1/2 - h)^2 + 1/2) / (1 + gamma),
            min(sqrt(2) h, 1 - 2h) / (2 gamma)).

The three conditions are necessary and sufficient, because they cover Ca-Ca, Ca-B and B-B pairs.

## Attainment for every positive size ratio

Put

    x* = (sqrt(2) - 1)/2,
    h* = 1/2 - x* = 1 - sqrt(2)/2,
    gamma* = sqrt(2) - 1 = 2x*.

At h = h* the lemma gives exactly

    d_CaCa = 1,
    d_BB   = sqrt(2) h* = 1 - 2h* = 2x*,
    d_CaB  = sqrt(x*^2 + 1/2).

Thus radii rL = 1/2 and rS = x* are feasible: the same-species inequalities are equalities, and the mixed-species inequality is strict,

    (rL + rS)^2 = (1/2 + x*)^2 = 1/2 < x*^2 + 1/2 = d_CaB^2.

This is an exact comparison with no decimal rounding. Reducing either radius preserves non-overlap. Therefore both radius assignments
used in the upper-bound proof are feasible for every gamma > 0:

- If 0 < gamma <= gamma*, choose rL = 1/2 and rS = gamma/2 <= x*.
- If gamma >= gamma*, choose rS = x* and rL = x*/gamma <= 1/2.

The upper bounds in PROOF.md are consequently attained throughout their stated domains:

    phi_max(gamma) = (pi/6)(1 + 6 gamma^3)                       for 0 < gamma <= gamma*,
    phi_max(gamma) = (pi/6)(sqrt(2) - 1)^3 (6 + gamma^(-3))     for gamma >= gamma*.

The first branch is strictly increasing and the second strictly decreasing. Their unique maximum as a function of gamma is at gamma*.
Equality there requires rL = 1/2 and rS = x*. The inequalities

    rS <= x,     rS <= (1/2 - x)/sqrt(2)

then force x = x*. Thus the global optimum is unique in the chosen Wyckoff representative, up to the usual lattice symmetries and overall scale.
This does not assert uniqueness of the optimizing coordinate at other values of gamma.

## The five B-B contacts at the optimum

Since 1/4 < h* < 1/2, a B site has one closest opposite-sign site on its own axis, across the adjacent cell boundary:
1 - 2h* < 2h*. It also has four closest sites on the other two axes, one for each sign on each axis; for these the minimizing
integer translation is zero. All five distances equal 2x*. Same-axis, same-sign translates have distance at least 1 > 2x*.
The coordinatewise minima are unique in these cases, so there are exactly five B contacts, not merely five exhibited contacts.
Ca has its six simple-cubic Ca contacts, and there are no Ca-B contacts.

## Evidence, scope and credit

The proof is the coordinate classification above. Finite replay of the existing packing verifier corroborates the result but is not
the justification for the infinite-lattice claim. The model remains two sphere radii on the fixed cubic CaB6 site family, with one internal
coordinate; no assertion concerns a real-material transition or optimality among other binary structures.

This note builds on the construction and upper-bound argument in this repository, directed by Ahraz Arifuddin (Aska Digital).
The existing prior-art qualifications in PROOF.md, including Horst and Travesset's plotted packing curve, remain in force.
Prepared with OpenAI Codex (GPT-6); no independent human review is claimed.
