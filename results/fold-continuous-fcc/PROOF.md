# An explicit continuous motion from the cab contact framework to fcc

Evidence: **PROVEN (written construction and exact finite certificate), internal review only.**
This is a proposed research result awaiting independent specialist review, not a claim of historical novelty.
Replay: `make verify-fold-continuous`. The exact replay uses the Python standard library only.

## Statement and model

The one-cell cab contact framework admits a smooth, bar-preserving motion from its open cubic configuration to
face-centred cubic sphere packing. Every pair of distinct balls remains non-overlapping throughout the motion.
Density increases strictly from the open-state value to `pi/sqrt(18)`. The lattice is allowed to shear.

Use units `K = 2 - sqrt(2)`. In these units the six motif centres are
`P = {+e_1,-e_1,+e_2,-e_2,+e_3,-e_3}`, each ball has radius `1/sqrt(2)`, and every bar has length `sqrt(2)`.
The 12 internal edges join different coordinate axes. The three bridges join `+e_i` to `-e_i` in the cell translated
by lattice row `L_i`. A configuration consists of all `p + n L`, with `p in P`, `n in Z^3`.
Multiplying by K and adding `(1,1,1)` gives exactly the repository's canonical initial motif and ball radius.

This is a new explicit witness. It does not certify the interpolants of the previously stored numerical paths,
which retain their sampled-continuity labels. It makes no assertion about a physical transition in CaB6.

## 1. The curve

Let

    C = [ 0 -1 -1 ]
        [ 1  0 -1 ]
        [ 1 -1  0 ].

For `0 <= s <= 1`, set

    t = 2s/(1+s^2),
    a = sqrt(2)(1-s^2)/(1+s^2),
    d = 2+a,
    L(s) = d I + t C.

Then `a >= 0`, `0 <= t <= 1`, and `a^2 + 2t^2 = 2`. At the start,
`L(0) = (2+sqrt(2)) I`, hence `K L(0) = 2 I`. At the endpoint,

    L(1) = [ 2 -1 -1 ]
           [ 1  2 -1 ]
           [ 1 -1  2 ].

The motif stays fixed in this frame. Internal edges therefore retain their lengths. Each bridge vector is
`L_i - 2e_i = a e_i + t C_i`. Since `C_ii = 0` and `|C_i|^2 = 2`, its squared length is
`a^2 + 2t^2 = 2` for every real s in the interval. All 15 bars are preserved exactly.

## 2. A uniform bound for all distant periodic images

The symmetric part of L is

    (L+L^T)/2 = [ d  0  0 ]
                 [ 0  d -t ]
                 [ 0 -t  d ],

with eigenvalues `d`, `d+t`, `d-t`. In particular `d-t >= 1`. For a real row vector w, Cauchy--Schwarz gives

    |w L| |w| >= w L w^T >= |w|^2,

so `|w L| >= |w|`. Thus the lattice is nonsingular throughout the motion. Any motif difference `v = p-q` has
`|v| <= 2`. If the integer translation satisfies `|n|^2 >= 12`, then

    |n L + v| >= |n L| - |v| >= sqrt(12)-2 > sqrt(2).

The last strict inequality follows from `12 > (2+sqrt(2))^2 = 6+4sqrt(2)`.
It remains to check exactly the 171 integer vectors with `|n|^2 <= 11`; they all lie in `[-3,3]^3`.
The set `P-P` has 19 distinct vectors. Excluding only `(n,v)=(0,0)` leaves `171*19-1 = 3248` distance conditions.

## 3. Exact non-overlap on the entire parameter interval

For each of those pairs define vectors with coordinates in `Q(sqrt(2))`:

    A0 = (2+sqrt(2)) n + v,
    A1 = 2 n C,
    A2 = (2-sqrt(2)) n + v.

The squared-distance condition, after multiplying by the positive denominator, is the quartic inequality

    f(s) = (1+s^2)^2 (|n L(s)+v|^2 - 2)
         = |A0 + A1 s + A2 s^2|^2 - 2(1+s^2)^2 >= 0.

Its power coefficients are

    c0 = |A0|^2 - 2,
    c1 = 2 A0.A1,
    c2 = |A1|^2 + 2 A0.A2 - 4,
    c3 = 2 A1.A2,
    c4 = |A2|^2 - 2.

Write the same polynomial in the degree-four Bernstein basis:

    f(s) = sum(k=0..4) b_k binom(4,k) s^k (1-s)^(4-k),
    b_k = sum(j=0..k) c_j binom(k,j)/binom(4,j).

Every basis function is nonnegative on `[0,1]`. The exact enumeration in `curve.py` finds all **16,240 coefficients
nonnegative**, including 114 zero coefficients; 18 distance polynomials vanish identically. No interval subdivision
or numerical tolerance is used. Coefficients are rational pairs representing `u+v sqrt(2)`. Their signs are determined
by rational comparisons and, when the terms have opposite signs, by comparing `u^2` with `2v^2`.

This is an exhaustive finite certificate after the infinite-image bound, not a parameter scan. Together the bound
and coefficient signs prove the distance inequality for every pair and every real `s in [0,1]`.
The code is part of the proof's trusted arithmetic implementation; it has not been formalised in a proof assistant.

## 4. The endpoint is fcc

The endpoint lattice has integer rows of even coordinate sum and determinant 12. It is therefore a sublattice of
`D3 = {z in Z^3 : sum(z_i) is even}`, which has covolume 2, with index 6.
The six motif points have odd coordinate sum and occupy distinct lattice cosets: equality of two translated centres
would contradict Section 3's strictly positive minimum separation. They consequently fill all six cosets of the
odd-parity translate of D3. This is the fcc lattice, whose shortest distance is `sqrt(2)`.

The exact endpoint replay also finds 12 contacts at each of the six motif points. The packing fraction is

    6 * (4pi/3) * (1/sqrt(2))^3 / 12 = pi/sqrt(18).

The endpoint matrix is the fcc matrix already recorded in [fold-m1](../fold-m1/PROOF.md).
The contribution here is the explicit motion and its all-parameter certificate.

## 5. Density is strictly increasing

In the t parametrization, `a = sqrt(2-2t^2)`. Direct expansion gives

    D(t) = det L = d^3 + d t^2 + 2t^3.

For `0 < t < 1`, differentiation and `a^2+2t^2=2` give

    D'(t) = (2t/a) [-16-10a+3t^2+3at].

Since `t <= 1` and `a >= 0`, the bracket is at most `-13-7a < 0`.
The determinant decreases strictly on the open interval and, by continuity, between any two distinct parameter
values in the closed interval. Also `t(s)` is strictly increasing on `[0,1]`. Sphere volume is constant, so density
is strictly increasing throughout the motion, including between any stored samples.

## Reproduction, review, and limits

`verify.py` replays the 3248 exact quartics through two separately written implementations, checks the integer endpoint
and bridge identities, tests Bernstein basis conversion, and rejects an overlapping sign-changed motion, a polynomial
positive at its endpoints but negative inside, and an incorrectly accepted negative coefficient near cancellation.
`independent.py` uses polynomial convolution and reverses every Bernstein expansion; it also independently checks
the endpoint cosets and image enumeration. The certificate hashes this proof, both implementations, the verifier,
and shared certificate helper. The image bound and monotonicity argument are written proofs, not automatically
formalised by that replay.

To compare a derived sample with the existing, unchanged numerical checkers:

    python results/fold-continuous-fcc/curve.py --sample /tmp/continuous-fcc.npz
    python autoresearch-fold/verify.py /tmp/continuous-fcc.npz
    python results/fold-close-packing/chk_independent.py /tmp/continuous-fcc.npz

The default has 501 samples. Their straight-line interpolation is not the analytic curve and is not asserted to
preserve bars exactly. Those numerical checks are supplemental; the continuous claim rests on Sections 1--5.
An m=1 motion lifts to any integer supercell by repeating it. This proves existence for all supercell sizes, without
claiming a fixed sampling meets every larger-cell checker's step-size limit.

Prior art: bridged rigid-unit flexibility and cross-polytope deployment have precedents in
[Goodwin (2006)](https://doi.org/10.1103/PhysRevB.74.134302) and
[Borcea and Streinu (2021), Section 7](https://doi.org/10.1177/1081286520948116).
The existing project's [fold prior-art survey](../../autoresearch-fold/PRIOR-ART-FOLD.md) and
[dimensional literature comparison in PR #27](https://github.com/aska-digital/incentre-lattice/pull/27)
explain related models. Whether this exact equal-bar, non-overlapping curve was previously published remains open.
This work does not establish material stability, energy barriers, or a physical mechanism in real hexaborides.

AI assistance: GPT-6 derived the curve and proof and implemented the exact replay; GPT-6.1 Sol independently checked
the polynomial enumeration and analytic steps. These are internal reviews. No independent human review is claimed.
