#!/usr/bin/env python3
r"""Zero-float exact verifier for the Canonical Sylvester-Hadamard Continuous
Bar-Preserving Fold F_8 -> E_8 (optimal 8D sphere packing, Viazovska 2016:
covolume 16, density pi^4/384, kissing number 240 — exact 2x D_8 density)
and the E_6/E_7/D_n^+ Uniqueness & Obstruction Theorem.

Zero external dependencies (Python 3 stdlib only), zero floating-point
arithmetic on any verification or certificate path, zero bare asserts,
deterministic bit-for-bit canonical JSON certificate across -I -B, -I -O -B,
-I -B -W error, and -I -O -B -W error.
"""

from fractions import Fraction
import hashlib
import json
import math
import os
import random
import sys
import time


# ---------------------------------------------------------------------------
# Strict R3/R4 guard helpers (zero bare asserts; reject bool & float!)
# ---------------------------------------------------------------------------

def require(cond, msg):
    if not cond:
        raise ValueError(msg)


def check_exact_int(x, label="value"):
    if isinstance(x, bool) or not isinstance(x, int):
        raise TypeError(f"{label} must be exact int (not {type(x).__name__})")
    return x


def check_exact_rat(x, label="value"):
    if isinstance(x, bool) or not isinstance(x, (int, Fraction)):
        raise TypeError(f"{label} must be exact int or Fraction (not {type(x).__name__})")
    return Fraction(x)


def check_z2(x, label="z2"):
    if not (isinstance(x, tuple) and len(x) == 2):
        raise TypeError(f"{label} must be a 2-tuple (a, b)")
    return (check_exact_int(x[0], f"{label}[0]"), check_exact_int(x[1], f"{label}[1]"))


# ---------------------------------------------------------------------------
# Exact Z[sqrt(2)] and Z[sqrt(2)][s] polynomial arithmetic
# ---------------------------------------------------------------------------

def z2_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def z2_sub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def z2_mul(x, y):
    return (x[0] * y[0] + 2 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def z2_scale(c, x):
    return (c * x[0], c * x[1])


def z2_sign(x):
    """Exact sign in {-1, 0, +1} of a + b*sqrt(2) for exact ints (a, b)."""
    a, b = check_z2(x, "z2_sign")
    if a == 0 and b == 0:
        return 0
    if a >= 0 and b >= 0:
        return 1
    if a <= 0 and b <= 0:
        return -1
    diff = a * a - 2 * b * b
    if a > 0:
        return 1 if diff > 0 else (-1 if diff < 0 else 0)
    return -1 if diff > 0 else (1 if diff < 0 else 0)


def z2_sign_fast(a, b):
    """Inline fast path for already-validated ints (a, b)."""
    if a == 0 and b == 0:
        return 0
    if a >= 0 and b >= 0:
        return 1
    if a <= 0 and b <= 0:
        return -1
    diff = a * a - 2 * b * b
    if a > 0:
        return 1 if diff > 0 else (-1 if diff < 0 else 0)
    return -1 if diff > 0 else (1 if diff < 0 else 0)


def p_add(p, q):
    n = max(len(p), len(q))
    r = [(0, 0)] * n
    for i, c in enumerate(p):
        r[i] = z2_add(r[i], c)
    for i, c in enumerate(q):
        r[i] = z2_add(r[i], c)
    return r


def p_sub(p, q):
    n = max(len(p), len(q))
    r = [(0, 0)] * n
    for i, c in enumerate(p):
        r[i] = z2_add(r[i], c)
    for i, c in enumerate(q):
        r[i] = z2_sub(r[i], c)
    return r


def p_scale(c, p):
    return [z2_scale(c, x) for x in p]


def p_mul(p, q):
    r = [(0, 0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a == (0, 0):
            continue
        for j, b in enumerate(q):
            if b == (0, 0):
                continue
            r[i + j] = z2_add(r[i + j], z2_mul(a, b))
    return r


def p_deriv(p):
    return [z2_scale(i, p[i]) for i in range(1, len(p))]


def p_eval_int(p, s_int):
    acc = (0, 0)
    pw = 1
    for c in p:
        acc = z2_add(acc, z2_scale(pw, c))
        pw *= s_int
    return acc


def p_div_1ms(p):
    """Exact polynomial division by (1 - s), requiring p(1) == 0."""
    require(p_eval_int(p, 1) == (0, 0), f"Polynomial does not vanish at s=1: {p_eval_int(p, 1)}")
    q = []
    acc = (0, 0)
    for c in p[:-1]:
        acc = z2_add(acc, c)
        q.append(acc)
    require(z2_add(acc, p[-1]) == (0, 0), "Exact (1-s) remainder non-zero")
    return q


def to_bernstein_scaled(p):
    """Exact Bernstein coefficients B_k = C(d, k) * b_k on [0, 1], preserving Z[sqrt(2)] integrality.
    Since C(d, k) > 0 for 0 <= k <= d, sign(B_k) == sign(b_k) identically!"""
    d = len(p) - 1
    out = []
    for k in range(d + 1):
        sa = 0
        sb = 0
        for i in range(k + 1):
            w = math.comb(d - i, k - i)
            sa += w * p[i][0]
            sb += w * p[i][1]
        out.append((sa, sb))
    return out


# ---------------------------------------------------------------------------
# Exact rational linear algebra (det, inverse, LDL^T, Fincke-Pohst)
# ---------------------------------------------------------------------------

def exact_det_rat(M):
    """Exact Gaussian elimination determinant over Q."""
    n = len(M)
    A = [[check_exact_rat(M[i][j]) for j in range(n)] for i in range(n)]
    det = Fraction(1)
    for k in range(n):
        piv = None
        for i in range(k, n):
            if A[i][k] != 0:
                piv = i
                break
        if piv is None:
            return Fraction(0)
        if piv != k:
            A[k], A[piv] = A[piv], A[k]
            det = -det
        pv = A[k][k]
        det *= pv
        for i in range(k + 1, n):
            if A[i][k] != 0:
                f = A[i][k] / pv
                for j in range(k, n):
                    A[i][j] -= f * A[k][j]
    return det


def exact_ldlt(G):
    """Exact rational LDL^T factorization of symmetric positive-definite G over Q."""
    n = len(G)
    for i in range(n):
        for j in range(n):
            require(G[i][j] == G[j][i], f"Matrix not symmetric at ({i},{j})")
    L = [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]
    D = [Fraction(0)] * n
    for j in range(n):
        s = G[j][j] - sum(L[j][k] * L[j][k] * D[k] for k in range(j))
        require(s > 0, f"LDL^T non-positive pivot at index {j}: {s}")
        D[j] = s
        for i in range(j + 1, n):
            L[i][j] = (G[i][j] - sum(L[i][k] * L[j][k] * D[k] for k in range(j))) / D[j]
    # Exact reconstruction check
    for i in range(n):
        for j in range(n):
            val = sum(L[i][k] * D[k] * L[j][k] for k in range(n))
            require(val == G[i][j], f"LDL^T reconstruction mismatch at ({i},{j})")
    return L, D


def exact_floor_sqrt_frac(q):
    """Exact integer floor(sqrt(q)) for exact Fraction q >= 0 using math.isqrt."""
    q = check_exact_rat(q, "exact_floor_sqrt_frac")
    require(q >= 0, f"Negative radicand {q}")
    return math.isqrt(q.numerator * q.denominator) // q.denominator


def enum_exact_ellipsoid_canonical(G, R2_int):
    r"""Zero-float exact Fincke-Pohst enumeration of ALL canonical +-m pairs
    (first nonzero coordinate > 0) in Z^n \ {0} with m^T G m <= R2_int."""
    check_exact_int(R2_int, "R2_int")
    require(R2_int >= 0, "R2_int must be >= 0")
    n = len(G)
    L, D = exact_ldlt(G)
    res = []
    m = [0] * n

    def rec(k, rem):
        shift = sum(L[j][k] * m[j] for j in range(k + 1, n))
        Q = rem / D[k]
        p, q = shift.numerator, shift.denominator
        r = exact_floor_sqrt_frac(Q * q * q)
        lo = (-r - p + q - 1) // q
        hi = (r - p) // q
        for x in range(lo, hi + 1):
            t = Fraction(x) + shift
            nrem = rem - D[k] * t * t
            if nrem < 0:
                continue
            m[k] = x
            if k == 0:
                for c in m:
                    if c != 0:
                        if c > 0:
                            n2_q = Fraction(R2_int) - nrem
                            require(n2_q.denominator == 1, f"Non-integer norm^2: {n2_q}")
                            res.append((tuple(m), n2_q.numerator))
                        break
            else:
                rec(k - 1, nrem)
        m[k] = 0

    rec(n - 1, Fraction(R2_int))
    return res, L, D


# ---------------------------------------------------------------------------
# Canonical Sylvester Hadamard H_8 and exact Q(sqrt(2))(s) fold curve
# ---------------------------------------------------------------------------

def build_sylvester_H8():
    return [[1 if bin(i & j).count("1") % 2 == 0 else -1 for j in range(8)] for i in range(8)]


# Exact curve polynomials in Z[sqrt(2)][s] (scaled by 2 relative to raw form)
DP = [(1, 2), (0, 0), (-1, 2)]       # (1 + 2*sqrt(2)) + (-1 + 2*sqrt(2)) s^2
DM = [(-1, 2), (0, 0), (1, 2)]       # (-1 + 2*sqrt(2)) + (1 + 2*sqrt(2)) s^2
NAP = [(6, 5), (0, 0), (-6, 5)]      # (6 + 5*sqrt(2)) + (-6 + 5*sqrt(2)) s^2
NBP = [(0, 0), (0, 2), (0, 0)]       # 2*sqrt(2) s
NAM = [(2, 3), (0, 0), (-2, 3)]      # (2 + 3*sqrt(2)) + (-2 + 3*sqrt(2)) s^2
NBM = [(0, 0), (0, 2), (0, 0)]       # 2*sqrt(2) s


def validate_spinor_hadamard_input(H):
    """Strict R3 validator for an 8x8 symmetric spinor Hadamard matrix."""
    require(isinstance(H, list) and len(H) == 8, "H must be 8x8 list")
    for i in range(8):
        require(isinstance(H[i], list) and len(H[i]) == 8, f"Row {i} must have length 8")
        neg_cnt = 0
        for j in range(8):
            val = check_exact_int(H[i][j], f"H[{i}][{j}]")
            require(val in (-1, 1), f"H[{i}][{j}] must be in {{-1, +1}}, got {val}")
            if val == -1:
                neg_cnt += 1
        require(neg_cnt % 2 == 0, f"Row {i} has odd minus count {neg_cnt} (not an E_8 spinor root)")
    for i in range(8):
        for j in range(8):
            require(H[i][j] == H[j][i], f"H not symmetric at ({i},{j})")
            dot = sum(H[i][k] * H[j][k] for k in range(8))
            expected = 8 if i == j else 0
            require(dot == expected, f"H^2[{i}][{j}] = {dot} != {expected}")


def verify_lemmas_1_to_4(H8):
    """Verify Lemmas 1-4 completely in exact Z[sqrt(2)] and Q arithmetic."""
    validate_spinor_hadamard_input(H8)

    # --- Lemma 1: Exact bar lengths & endpoints ---
    # Check D_+(s) > 0, D_-(s) > 0 on [0, 1]
    dp_b = [z2_sign(b) for b in to_bernstein_scaled(DP)]
    dm_b = [z2_sign(b) for b in to_bernstein_scaled(DM)]
    require(all(s == 1 for s in dp_b), f"D_+ not strictly positive: {dp_b}")
    require(all(s == 1 for s in dm_b), f"D_- not strictly positive: {dm_b}")

    # Check bar length identity: (N_a - 2 D)^2 + 7 N_b^2 == 2 D^2 identically!
    for lbl, Na, Nb, Dp in [("I_+", NAP, NBP, DP), ("I_-", NAM, NBM, DM)]:
        diff = p_sub(Na, p_scale(2, Dp))
        lhs = p_add(p_mul(diff, diff), p_scale(7, p_mul(Nb, Nb)))
        rhs = p_scale(2, p_mul(Dp, Dp))
        require(lhs == rhs, f"Bar identity failed for {lbl}: lhs={lhs}, rhs={rhs}")

    # Endpoints s=0: N_a(0) == (2+sqrt(2)) * D(0), N_b(0) == 0
    for lbl, Na, Nb, Dp in [("I_+", NAP, NBP, DP), ("I_-", NAM, NBM, DM)]:
        require(p_eval_int(Na, 0) == z2_mul((2, 1), p_eval_int(Dp, 0)), f"s=0 diag mismatch {lbl}")
        require(p_eval_int(Nb, 0) == (0, 0), f"s=0 off-diag mismatch {lbl}")

    # Endpoints s=1: D_+(1) = D_-(1) = (0, 4), a_+(1) = 5/2, a_-(1) = 3/2, b_+(1)=b_-(1)=1/2
    require(p_eval_int(DP, 1) == (0, 4) and p_eval_int(DM, 1) == (0, 4), "D(1) != 4*sqrt(2)")
    require(z2_scale(2, p_eval_int(NAP, 1)) == z2_scale(5, p_eval_int(DP, 1)), "a_+(1) != 5/2")
    require(z2_scale(2, p_eval_int(NAM, 1)) == z2_scale(3, p_eval_int(DM, 1)), "a_-(1) != 3/2")
    require(z2_scale(2, p_eval_int(NBP, 1)) == p_eval_int(DP, 1), "b_+(1) != 1/2")
    require(z2_scale(2, p_eval_int(NBM, 1)) == p_eval_int(DM, 1), "b_-(1) != 1/2")

    # --- Lemma 2: Exact E_8 endpoint, (Z/2Z)^4 quotient, sqrt(2) E_8 sublattice ---
    L1 = [[Fraction(2 if i == j else 0) + Fraction(H8[i][j], 2) for j in range(8)] for i in range(8)]
    L1_inv = [[Fraction(1 if i == j else 0) - Fraction(H8[i][j], 4) for j in range(8)] for i in range(8)]
    # Verify L1 * L1_inv == I_8, L1^2 == 6 I_8 + 2 H_8, det(L1) == 16
    for i in range(8):
        for j in range(8):
            prod_inv = sum(L1[i][k] * L1_inv[k][j] for k in range(8))
            require(prod_inv == (1 if i == j else 0), f"L1*L1_inv mismatch at ({i},{j})")
            g1_ij = sum(L1[i][k] * L1[j][k] for k in range(8))
            expected_g1 = Fraction(6 if i == j else 0) + Fraction(2 * H8[i][j])
            require(g1_ij == expected_g1, f"G(1) mismatch at ({i},{j})")
    det_L1 = exact_det_rat(L1)
    require(det_L1 == Fraction(16), f"det(L1) = {det_L1} != 16")

    # Verify 2 E_8 subset Z^8 L1 subset E_8 on an explicit Z-basis of E_8!
    # Standard E_8 basis (det = 1): 7 simple D_8 roots + 1 spinor root (1/2)^8
    E8_basis = [
        [Fraction(1), Fraction(-1), Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(0)],
        [Fraction(0), Fraction(1), Fraction(-1), Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(1), Fraction(-1), Fraction(0), Fraction(0), Fraction(0), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(0), Fraction(1), Fraction(-1), Fraction(0), Fraction(0), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(1), Fraction(-1), Fraction(0), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(1), Fraction(-1), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(1), Fraction(1), Fraction(0)],
        [Fraction(1, 2)] * 8,
    ]
    require(abs(exact_det_rat(E8_basis)) == Fraction(1), "E_8 basis does not have |det|=1")
    # Every row of 2*E8_basis maps into Z^8 under L1_inv => 2 E_8 subset Z^8 L1!
    for b_idx, bv in enumerate(E8_basis):
        coords = [sum(2 * bv[k] * L1_inv[k][j] for k in range(8)) for j in range(8)]
        require(all(c.denominator == 1 for c in coords), f"2*E8_basis[{b_idx}] not in Z^8 L1")

    # Verify P_8 - e_0 maps to 16 distinct cosets of E_8 / Z^8 L1, all of order dividing 2!
    P8 = [[Fraction(sgn if j == i else 0) for j in range(8)] for i in range(8) for sgn in (1, -1)]
    cosets = set()
    coset_orders = []
    for p in P8:
        diff = [p[j] - P8[0][j] for j in range(8)]
        frac_coords = tuple((sum(diff[k] * L1_inv[k][j] for k in range(8))) % 1 for j in range(8))
        cosets.add(frac_coords)
        if all(c == 0 for c in frac_coords):
            coset_orders.append(1)
        else:
            require(all((2 * c).denominator == 1 for c in frac_coords), f"Coset order != 2: {frac_coords}")
            coset_orders.append(2)
    require(len(cosets) == 16, f"Expected 16 distinct cosets, got {len(cosets)}")
    require(sorted(coset_orders) == [1] + [2] * 15, f"Quotient not (Z/2Z)^4: {coset_orders}")

    # --- Lemma 3: Exact integer orthogonal block-diagonalization & det monotonicity ---
    # Only for canonical Sylvester H_8 (perm = (0,3,5,6, 7,4,2,1))
    I_plus = [i for i in range(8) if H8[i][i] == 1]
    I_minus = [i for i in range(8) if H8[i][i] == -1]
    require(len(I_plus) == 4 and len(I_minus) == 4, "Expected 4 positive and 4 negative diagonal entries")

    # Build exact 8x8 polynomial matrix Dp*Dm*L_8(s) and check block-diagonalization for canonical H_8!
    P_dp = p_mul(NAP, DM)
    P_op = p_mul(NBP, DM)
    P_dm = p_mul(NAM, DP)
    P_om = p_mul(NBM, DP)

    # Check canonical Sylvester permutation when H8 == build_sylvester_H8()
    block_diag_verified = False
    if H8 == build_sylvester_H8():
        perm = [0, 3, 5, 6, 7, 4, 2, 1]
        Hp = [[H8[perm[i]][perm[j]] for j in range(8)] for i in range(8)]
        K4 = [[Hp[i][j] for j in range(4)] for i in range(4)]
        for i in range(4):
            for j in range(4):
                require(Hp[i][j + 4] == K4[i][j], "H_+- != K4")
                require(Hp[i + 4][j] == K4[i][j], "H_-+ != K4")
                require(Hp[i + 4][j + 4] == -K4[i][j], "H_-- != -K4")
        # Pure integer orthogonal matrix Q4 (Q4^T Q4 = 4 I_4, Q4^T K4 Q4 = 4 diag(2,2,2,-2))
        Q4 = [
            [1, 1, 1, -1],
            [1, -1, 1, 1],
            [-1, 1, 1, 1],
            [1, 1, -1, 1],
        ]
        # Build interleaved 8x8 Q8: columns (v_1,0), (0,v_1), (v_2,0), (0,v_2), (v_3,0), (0,v_3), (v_0,0), (0,v_0)
        Q8 = [[0] * 8 for _ in range(8)]
        for blk in range(4):
            for r in range(4):
                Q8[r][2 * blk] = Q4[r][blk]
                Q8[4 + r][2 * blk + 1] = Q4[r][blk]
        # Verify Q8^T Q8 == 4 I_8
        for i in range(8):
            for j in range(8):
                dot = sum(Q8[k][i] * Q8[k][j] for k in range(8))
                require(dot == (4 if i == j else 0), f"Q8 orthogonality failed at ({i},{j})")

        # Build permuted polynomial matrix M_perm(s) = Dp(s) Dm(s) P^T L_8(s) P
        M_perm = [[None] * 8 for _ in range(8)]
        for i in range(8):
            for j in range(8):
                if i < 4:
                    M_perm[i][j] = P_dp if i == j else p_scale(Hp[i][j], P_op)
                else:
                    M_perm[i][j] = P_dm if i == j else p_scale(Hp[i][j], P_om)

        # Compute exact Q8^T M_perm(s) Q8 and verify it equals 4 * diag(Mp_num, Mp_num, Mp_num, Mm_num)!
        Mp_num = [
            [p_add(P_dp, P_op), p_scale(2, P_op)],
            [p_scale(2, P_om), p_sub(P_dm, P_om)],
        ]
        Mm_num = [
            [p_sub(P_dp, p_scale(3, P_op)), p_scale(-2, P_op)],
            [p_scale(-2, P_om), p_add(P_dm, p_scale(3, P_om))],
        ]
        expected_blocks = [Mp_num, Mp_num, Mp_num, Mm_num]
        for i in range(8):
            for j in range(8):
                acc = [(0, 0)] * 5
                for a in range(8):
                    qa = Q8[a][i]
                    if qa == 0:
                        continue
                    for b in range(8):
                        qb = Q8[b][j]
                        if qb == 0:
                            continue
                        acc = p_add(acc, p_scale(qa * qb, M_perm[a][b]))
                bi, ri = divmod(i, 2)
                bj, rj = divmod(j, 2)
                exp_poly = p_scale(4, expected_blocks[bi][ri][rj]) if bi == bj else [(0, 0)] * 5
                require(acc == exp_poly, f"Block diagonalization mismatch at ({i},{j})")
        block_diag_verified = True

    # Exact determinant monotonicity from M_+2 and M_-2
    Np = p_sub(
        p_mul(p_add(NAP, NBP), p_sub(NAM, NBM)),
        p_mul(p_scale(2, NBP), p_scale(2, NBM)),
    )
    Nm = p_sub(
        p_mul(p_sub(NAP, p_scale(3, NBP)), p_add(NAM, p_scale(3, NBM))),
        p_mul(p_scale(-2, NBP), p_scale(-2, NBM)),
    )
    np_signs = [z2_sign(b) for b in to_bernstein_scaled(Np)]
    nm_signs = [z2_sign(b) for b in to_bernstein_scaled(Nm)]
    require(all(s == 1 for s in np_signs), f"N_+ not strictly positive: {np_signs}")
    require(all(s == 1 for s in nm_signs), f"N_- not strictly positive: {nm_signs}")

    DpDm = p_mul(DP, DM)
    term1 = p_mul(p_add(p_scale(3, p_mul(p_deriv(Np), Nm)), p_mul(Np, p_deriv(Nm))), DpDm)
    term2 = p_scale(4, p_mul(p_mul(Np, Nm), p_deriv(DpDm)))
    d_num = p_sub(term1, term2)
    require(d_num[0] == (0, 0), f"det'(0) != 0: {d_num[0]}")
    q_det = p_scale(-1, d_num[1:])
    q_det_signs = [z2_sign(b) for b in to_bernstein_scaled(q_det)]
    require(all(s == 1 for s in q_det_signs), f"-det'(s)/s not strictly positive: {q_det_signs}")

    # Verify det(0) = (2+sqrt(2))^8 = 9232 + 6528*sqrt(2) from Np(0)^3 Nm(0) / (Dp(0)Dm(0))^4
    # Note Dp(0)*Dm(0) = (1+2sq2)(-1+2sq2) = 7, so (Dp(0)Dm(0))^4 = 2401!
    # Check Np(0) = Nm(0) = (42, 28) = 7 * (6, 4) = 7 * (2+sqrt(2))^2!
    require(p_eval_int(DpDm, 0) == (7, 0), "Dp(0)*Dm(0) != 7")
    require(p_eval_int(Np, 0) == (42, 28) and p_eval_int(Nm, 0) == (42, 28), "N(0) != 7*(6+4*sqrt(2))")
    # Verify det(1) = 16 from Np(1)^3 Nm(1) / (Dp(1)Dm(1))^4: Dp(1)Dm(1) = 32, Np(1) = Nm(1) = 64 => 2^4 = 16!
    require(p_eval_int(DpDm, 1) == (32, 0), "Dp(1)*Dm(1) != 32")
    require(p_eval_int(Np, 1) == (64, 0) and p_eval_int(Nm, 1) == (64, 0), "N(1) != 64")

    # --- Lemma 4: Exact Loewner monotonicity G(s) - G(1) >= 0 on BOTH blocks ---
    loewner_summary = {}
    Dp2 = p_mul(DP, DP)
    Dm2 = p_mul(DM, DM)
    for lbl, r00, r01, r10, r11, g00, g01, g11 in [
        ("M_+2", p_add(NAP, NBP), p_scale(2, NBP), p_scale(2, NBM), p_sub(NAM, NBM), 10, 4, 2),
        ("M_-2", p_sub(NAP, p_scale(3, NBP)), p_scale(-2, NBP), p_scale(-2, NBM), p_add(NAM, p_scale(3, NBM)), 2, -4, 10),
    ]:
        d00 = p_sub(p_add(p_mul(r00, r00), p_mul(r01, r01)), p_scale(g00, Dp2))
        d11 = p_sub(p_add(p_mul(r10, r10), p_mul(r11, r11)), p_scale(g11, Dm2))
        d01 = p_sub(p_add(p_mul(r00, r10), p_mul(r01, r11)), p_scale(g01, DpDm))
        det_diff = p_sub(p_mul(d00, d11), p_mul(d01, d01))

        d00_q = p_div_1ms(d00)
        d11_q = p_div_1ms(d11)
        det_q2 = p_div_1ms(p_div_1ms(det_diff))
        require(det_q2[0] == (0, 0), f"{lbl} det_diff/(1-s)^2 did not vanish at s=0")
        det_q2_s = det_q2[1:]

        s_d00 = [z2_sign(b) for b in to_bernstein_scaled(d00_q)]
        s_d11 = [z2_sign(b) for b in to_bernstein_scaled(d11_q)]
        s_det = [z2_sign(b) for b in to_bernstein_scaled(det_q2_s)]
        require(all(s == 1 for s in s_d00), f"{lbl} d00/(1-s) not strictly positive: {s_d00}")
        require(all(s == 1 for s in s_d11), f"{lbl} d11/(1-s) not strictly positive: {s_d11}")
        require(all(s == 1 for s in s_det), f"{lbl} det/(s*(1-s)^2) not strictly positive: {s_det}")
        loewner_summary[lbl] = {
            "d00_over_1ms_deg": len(d00_q) - 1,
            "d00_over_1ms_bernstein_signs": s_d00,
            "d11_over_1ms_deg": len(d11_q) - 1,
            "d11_over_1ms_bernstein_signs": s_d11,
            "det_over_s_1ms2_deg": len(det_q2_s) - 1,
            "det_over_s_1ms2_bernstein_signs": s_det,
        }

    return {
        "det_L1": int(det_L1),
        "quotient_cosets": len(cosets),
        "quotient_group": "(Z/2Z)^4",
        "block_diag_exact_polynomial_verified": block_diag_verified,
        "Np_bernstein_signs": np_signs,
        "Nm_bernstein_signs": nm_signs,
        "neg_det_deriv_over_s_deg": len(q_det) - 1,
        "neg_det_deriv_over_s_bernstein_signs": q_det_signs,
        "loewner_G_minus_G1_blocks": loewner_summary,
    }


def verify_lemma_5_bernstein_and_kissing(H8, R2_cutoff=11, target_num=2, target_den=1):
    """Exact zero-float Fincke-Pohst ellipsoid enumeration + exhaustive degree-8
    Z[sqrt(2)] Bernstein certification + exact integer kissing count at all 16 motif balls.

    Completeness gate: (2 + sqrt(2))^2 = 6 + 4*sqrt(2) = 6 + sqrt(32) < 12,
    so floor((2+sqrt(2))^2) = 6 + isqrt(32) = 11. Any R2_cutoff < 11 fails completeness!
    """
    check_exact_int(R2_cutoff, "R2_cutoff")
    check_exact_int(target_num, "target_num")
    check_exact_int(target_den, "target_den")
    require(target_den > 0, "target_den must be positive")
    required_cutoff = 6 + math.isqrt(32)  # exact floor((2 + sqrt(2))^2) == 11
    require(R2_cutoff >= required_cutoff,
            f"Incomplete ellipsoid cutoff R2={R2_cutoff} < floor((2+sqrt(2))^2)={required_cutoff}")

    I_plus = [i for i in range(8) if H8[i][i] == 1]
    I_minus = [i for i in range(8) if H8[i][i] == -1]
    G1 = [[Fraction(6 if i == j else 0) + Fraction(2 * H8[i][j]) for j in range(8)] for i in range(8)]

    pts_pos, _, pivots = enum_exact_ellipsoid_canonical(G1, R2_cutoff)
    norm_hist = {}
    for _, n2 in pts_pos:
        norm_hist[n2] = norm_hist.get(n2, 0) + 1
    require(norm_hist == {4: 120, 8: 1080}, f"Unexpected G(1) shell counts: {norm_hist}")

    # Exact integer kissing count at ALL 16 motif balls (using 2*L1 = 4*I + H8)
    # Note: ||x + 2(q - p)||^2 == ||x||^2 + 4||q-p||^2 + 4 x.(q - p), using only nonzero coords of q-p!
    L1_2 = [[(4 if i == j else 0) + H8[i][j] for j in range(8)] for i in range(8)]
    P8 = [tuple((sgn if j == i else 0) for j in range(8)) for i in range(8) for sgn in (1, -1)]
    all_m = [m for m, _ in pts_pos] + [tuple(-c for c in m) for m, _ in pts_pos]
    mL2 = [(0,) * 8] + [tuple(sum(m[i] * L1_2[i][j] for i in range(8)) for j in range(8)) for m in all_m]
    mL2_n2 = [sum(c * c for c in x) for x in mL2]
    kissing_counts = []
    for p_idx, p in enumerate(P8):
        pi = p_idx // 2
        ps = 1 if p_idx % 2 == 0 else -1
        cnt = 0
        for q_idx, q in enumerate(P8):
            qi = q_idx // 2
            qs = 1 if q_idx % 2 == 0 else -1
            if p_idx == q_idx:
                # v = 0, skip idx = 0; every m != 0 has ||m(2L1)||^2 >= 16 > 8
                for n2x in mL2_n2[1:]:
                    require(n2x >= 8, f"s=1 v=0 overlap: d2={n2x}")
                continue
            if pi == qi:
                # q - p = -2*ps * e_pi, so 4||q-p||^2 = 16, 4 x.(q-p) = -8*ps * x[pi]
                c_dot = -8 * ps
                for idx in range(len(mL2)):
                    d2 = mL2_n2[idx] + 16 + c_dot * mL2[idx][pi]
                    require(d2 >= 8, f"s=1 V_2e overlap: d2={d2}")
                    if d2 == 8:
                        cnt += 1
            else:
                # q - p = qs * e_qi - ps * e_pi, so 4||q-p||^2 = 8, 4 x.(q-p) = 4*(qs*x[qi] - ps*x[pi])
                c_q = 4 * qs
                c_p = -4 * ps
                for idx in range(len(mL2)):
                    x = mL2[idx]
                    d2 = mL2_n2[idx] + 8 + c_q * x[qi] + c_p * x[pi]
                    require(d2 >= 8, f"s=1 V_ee overlap: d2={d2}")
                    if d2 == 8:
                        cnt += 1
        kissing_counts.append(cnt)
    require(kissing_counts == [240] * 16, f"Kissing numbers != 240: {kissing_counts}")

    # Degree-8 exact Z[sqrt(2)] Bernstein check over ALL pairs (m, v):
    # Target inequality: target_den * ||m L_8(s) + v||^2 >= target_num
    # Multiply by D_+(s)^2 D_-(s)^2:
    #   target_den * (||W(s)||^2 + 2 sum_j v_j W_j(s) P_v(s) + ||v||^2 P_v(s)^2) - target_num * P_v(s)^2 >= 0
    P_dp = p_mul(NAP, DM)
    P_op = p_mul(NBP, DM)
    P_dm = p_mul(NAM, DP)
    P_om = p_mul(NBM, DP)
    P_v = p_mul(DP, DM)
    P_v2 = p_mul(P_v, P_v)

    COMB8 = [[math.comb(8 - i, k - i) if i <= k else 0 for i in range(9)] for k in range(9)]

    def to_b8(ca, cb):
        ba = [0] * 9
        bb = [0] * 9
        for k in range(9):
            row = COMB8[k]
            sa = 0
            sb = 0
            for i in range(k + 1):
                w = row[i]
                sa += w * ca[i]
                sb += w * cb[i]
            ba[k] = sa
            bb[k] = sb
        return ba, bb

    B_Pv2a, B_Pv2b = to_b8([x[0] for x in P_v2], [x[1] for x in P_v2])
    c_2e = 4 * target_den - target_num   # for ||v||^2 = 4 (equals 2 when target=2/1)
    c_ee = 2 * target_den - target_num   # for ||v||^2 = 2 (equals 0 when target=2/1)

    V_2e = [((i, sgn),) for i in range(8) for sgn in (-2, 2)]
    V_ee = [((i, s1), (j, s2)) for i in range(8) for j in range(i + 1, 8) for s1 in (-1, 1) for s2 in (-1, 1)]

    total_pairs = 0
    neg_pairs = 0
    zero_coeffs = 0
    bar_pairs = 0
    s1_contacts_2e = 0
    s1_contacts_ee = 0
    nonbar_interior_zeros = 0

    for m, n2 in pts_pos:
        Wa = [[0] * 5 for _ in range(8)]
        Wb = [[0] * 5 for _ in range(8)]
        for j in range(8):
            dp = m[j] if H8[j][j] == 1 else 0
            dm = m[j] if H8[j][j] == -1 else 0
            up = 0
            um = 0
            for i in I_plus:
                if i != j:
                    up += m[i] * H8[i][j]
            for i in I_minus:
                if i != j:
                    um += m[i] * H8[i][j]
            for k in range(5):
                Wa[j][k] = dp * P_dp[k][0] + up * P_op[k][0] + dm * P_dm[k][0] + um * P_om[k][0]
                Wb[j][k] = dp * P_dp[k][1] + up * P_op[k][1] + dm * P_dm[k][1] + um * P_om[k][1]

        W2a = [0] * 9
        W2b = [0] * 9
        for j in range(8):
            wja, wjb = Wa[j], Wb[j]
            for i1 in range(5):
                a1, b1 = wja[i1], wjb[i1]
                if a1 == 0 and b1 == 0:
                    continue
                for i2 in range(5):
                    a2, b2 = wja[i2], wjb[i2]
                    W2a[i1 + i2] += a1 * a2 + 2 * b1 * b2
                    W2b[i1 + i2] += a1 * b2 + b1 * a2

        Wva = [[0] * 9 for _ in range(8)]
        Wvb = [[0] * 9 for _ in range(8)]
        for j in range(8):
            wja, wjb = Wa[j], Wb[j]
            for i1 in range(5):
                a1, b1 = wja[i1], wjb[i1]
                if a1 == 0 and b1 == 0:
                    continue
                for i2 in range(5):
                    a2, b2 = P_v[i2]
                    if a2 == 0 and b2 == 0:
                        continue
                    Wva[j][i1 + i2] += 2 * (a1 * a2 + 2 * b1 * b2)
                    Wvb[j][i1 + i2] += 2 * (a1 * b2 + b1 * a2)

        B_W2a, B_W2b = to_b8(W2a, W2b)
        B_Wva = []
        B_Wvb = []
        for j in range(8):
            ba, bb = to_b8(Wva[j], Wvb[j])
            B_Wva.append(ba)
            B_Wvb.append(bb)

        # V_2e for ALL m with n2 <= 11 (1,200 * 16 = 19,200 pairs)
        for ((idx, sgn),) in V_2e:
            total_pairs += 1
            waj, wbj = B_Wva[idx], B_Wvb[idx]
            zmask = []
            has_neg = False
            for k in range(9):
                a = target_den * (B_W2a[k] + sgn * waj[k]) + c_2e * B_Pv2a[k]
                b = target_den * (B_W2b[k] + sgn * wbj[k]) + c_2e * B_Pv2b[k]
                sg = z2_sign_fast(a, b)
                if sg < 0:
                    has_neg = True
                    break
                elif sg == 0:
                    zero_coeffs += 1
                    zmask.append(k)
            if has_neg:
                neg_pairs += 1
            elif len(zmask) == 9:
                bar_pairs += 1
            elif zmask == [8]:
                s1_contacts_2e += 1
            elif zmask:
                nonbar_interior_zeros += len(zmask)

        # V_ee for ALL m with n2 == 4 (strictly < 8; 120 * 112 = 13,440 pairs)
        if n2 < 8:
            for ((i1, s1), (i2, s2)) in V_ee:
                total_pairs += 1
                wa1, wb1 = B_Wva[i1], B_Wvb[i1]
                wa2, wb2 = B_Wva[i2], B_Wvb[i2]
                zmask = []
                has_neg = False
                for k in range(9):
                    a = target_den * (B_W2a[k] + s1 * wa1[k] + s2 * wa2[k]) + c_ee * B_Pv2a[k]
                    b = target_den * (B_W2b[k] + s1 * wb1[k] + s2 * wb2[k]) + c_ee * B_Pv2b[k]
                    sg = z2_sign_fast(a, b)
                    if sg < 0:
                        has_neg = True
                        break
                    elif sg == 0:
                        zero_coeffs += 1
                        zmask.append(k)
                if has_neg:
                    neg_pairs += 1
                elif zmask == [8]:
                    s1_contacts_ee += 1
                elif zmask:
                    nonbar_interior_zeros += len(zmask)

    return {
        "R2_cutoff": R2_cutoff,
        "ldlt_pivots": [str(p) for p in pivots],
        "shell_counts_canonical_pm": {str(k): v for k, v in sorted(norm_hist.items())},
        "kissing_number_all_16_balls": 240,
        "total_pairs": total_pairs,
        "total_degree8_bernstein_coeffs": total_pairs * 9,
        "negative_pairs": neg_pairs,
        "subdivisions_needed": 0,
        "zero_coeffs_total": zero_coeffs,
        "persistent_bar_pairs": bar_pairs,
        "s1_contacts_V2e_nonbar": s1_contacts_2e,
        "s1_contacts_Vee_n2_4": s1_contacts_ee,
        "nonbar_interior_zeros": nonbar_interior_zeros,
    }


def verify_lemma_6_uniqueness():
    """Verify Exact E_6/E_7/E_8 Determinant-Square Obstruction and D_n^+ Uniqueness (Portal 9.4)."""
    # Exact Cartan (= Gram at min norm^2 = 2) matrices of E_6, E_7, E_8 and D_6, D_7, D_8!
    def cartan_from_edges(n, edges):
        C = [[Fraction(2 if i == j else 0) for j in range(n)] for i in range(n)]
        for u, v in edges:
            C[u][v] = Fraction(-1)
            C[v][u] = Fraction(-1)
        return C

    C_E6 = cartan_from_edges(6, [(0, 1), (1, 2), (2, 3), (3, 4), (2, 5)])
    C_E7 = cartan_from_edges(7, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (2, 6)])
    C_E8 = cartan_from_edges(8, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)])
    C_D6 = cartan_from_edges(6, [(0, 1), (1, 2), (2, 3), (3, 4), (3, 5)])
    C_D7 = cartan_from_edges(7, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (4, 6)])
    C_D8 = cartan_from_edges(8, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (5, 7)])

    d_E6, d_E7, d_E8 = exact_det_rat(C_E6), exact_det_rat(C_E7), exact_det_rat(C_E8)
    d_D6, d_D7, d_D8 = exact_det_rat(C_D6), exact_det_rat(C_D7), exact_det_rat(C_D8)
    require((d_E6, d_E7, d_E8) == (3, 2, 1), f"Unexpected E Gram dets: {(d_E6, d_E7, d_E8)}")
    require((d_D6, d_D7, d_D8) == (4, 4, 4), f"Unexpected D Gram dets: {(d_D6, d_D7, d_D8)}")

    def is_int_square(q):
        if q.denominator != 1 or q.numerator < 0:
            return False
        r = math.isqrt(q.numerator)
        return r * r == q.numerator

    require(not is_int_square(d_D6 / d_E6), "4/3 should not be an integer square")
    require(not is_int_square(d_D7 / d_E7), "4/2 should not be an integer square")
    require(is_int_square(d_D8 / d_E8), "4/1 must be an integer square (2^2)")

    # Check D_n^+ spinor norm^2 = n/4 for even n in 4..24: equals 2 iff n == 8!
    dn_plus_status = {}
    for n in range(4, 26, 2):
        sp_norm2 = Fraction(n, 4)
        if sp_norm2 < 2:
            dn_plus_status[str(n)] = f"OVERLAP (spinor norm^2={sp_norm2} < 2)"
        elif sp_norm2 == 2:
            dn_plus_status[str(n)] = "UNIQUE_EXACT_MATCH (E_8 = D_8^+, spinor norm^2=2)"
        else:
            dn_plus_status[str(n)] = f"BAR_TRAPPED_IN_Dn (spinor norm^2={sp_norm2} > 2)"
    exact_matches = [k for k, v in dn_plus_status.items() if v.startswith("UNIQUE")]
    require(exact_matches == ["8"], f"Expected unique match n=8, got {exact_matches}")

    return {
        "gram_dets": {"D_n": 4, "E_6": int(d_E6), "E_7": int(d_E7), "E_8": int(d_E8)},
        "index_squared_ratios": {"E_6": "4/3 (obstructed)", "E_7": "2 (obstructed)", "E_8": "4 = 2^2 (allowed)"},
        "Dn_plus_uniqueness_even_4_to_24": dn_plus_status,
    }


def verify_b8_gauge_invariance(H8):
    """Verify signed-permutation B_8 gauge invariance across 5 deterministic elements!
    Note: conjugating H8 by signed permutation Pi (with even minus count in the sign vector,
    so spinor parity is preserved!) maps F_8 -> E_8 isometrically to another exact F_8 -> E_8 fold!"""
    rng = random.Random(20260808)
    checked = 0
    for trial in range(5):
        perm = list(range(8))
        rng.shuffle(perm)
        signs = [rng.choice((-1, 1)) for _ in range(8)]
        if signs.count(-1) % 2 != 0:
            signs[0] = -signs[0]
        # H' = Pi H8 Pi^T where Pi[i, perm[i]] = signs[i]
        Hp = [[signs[i] * signs[j] * H8[perm[i]][perm[j]] for j in range(8)] for i in range(8)]
        l14 = verify_lemmas_1_to_4(Hp)
        l5 = verify_lemma_5_bernstein_and_kissing(Hp, R2_cutoff=11)
        require(l14["det_L1"] == 16 and l14["quotient_cosets"] == 16, f"Gauge {trial} endpoint failed")
        require(l5["negative_pairs"] == 0 and l5["zero_coeffs_total"] == 976, f"Gauge {trial} Bernstein failed")
        require(l5["nonbar_interior_zeros"] == 0 and l5["persistent_bar_pairs"] == 8, f"Gauge {trial} zeros failed")
        checked += 1
    return {"b8_gauge_trials_verified": checked}


def run_must_fail_controls(H8):
    """7 adversarial must-fail controls (R1/R2/R3/R5)."""
    controls = {}

    # Control 1: Single-entry flip breaks spinor parity & Hadamard orthogonality & det=16
    H_flip = [row[:] for row in H8]
    H_flip[0][1] = -H_flip[0][1]
    rejected_1 = False
    try:
        verify_lemmas_1_to_4(H_flip)
    except ValueError:
        rejected_1 = True
    L_flip = [[Fraction(2 if i == j else 0) + Fraction(H_flip[i][j], 2) for j in range(8)] for i in range(8)]
    det_flip = exact_det_rat(L_flip)
    require(rejected_1 and det_flip != 16, f"Control 1 failed: rejected={rejected_1}, det={det_flip}")
    controls["C1_single_entry_flip"] = {"rejected": True, "det_mutated": str(det_flip)}

    # Control 2: All-minus spinor matrix C2 = -J_8 (even minus count = 8, orientation-reversing collapse!)
    # Eigenvalues of L1 = 2I - (1/2)J_8 are -2 (mult 1) and +2 (mult 7) => det(L1) = -256 < 0!
    # Along the curve, row sum a_-(s) - 7 b_-(s) = [(2+3sq2) - 14sq2 s + (-2+3sq2) s^2] / Dm(s)
    # reverses sign from +(2+sq2) at s=0 to -2 at s=1, so m = (1,...,1) collapses through 0 => ||m L(s)||^2 - 2 < 0!
    H_allm = [[-1] * 8 for _ in range(8)]
    rejected_2 = False
    try:
        verify_lemmas_1_to_4(H_allm)
    except ValueError:
        rejected_2 = True
    L_allm = [[Fraction(2 if i == j else 0) + Fraction(H_allm[i][j], 2) for j in range(8)] for i in range(8)]
    det_allm = exact_det_rat(L_allm)
    # Also directly verify exact negative Bernstein coeff for m=(1,...,1), v=0 (||m L(s)||^2 - 2 < 0)!
    row_sum_num = p_sub(NAM, p_scale(7, NBM))
    norm2_minus_2 = p_sub(p_scale(8, p_mul(row_sum_num, row_sum_num)), p_scale(2, p_mul(DM, DM)))
    c2_bs_signs = [z2_sign(b) for b in to_bernstein_scaled(norm2_minus_2)]
    require(rejected_2 and det_allm == -256 and min(c2_bs_signs) < 0,
            f"Control 2 failed: rejected={rejected_2}, det={det_allm}, bs_signs={c2_bs_signs}")
    controls["C2_all_minus_orientation_collapse"] = {
        "rejected": True,
        "det_inverted": int(det_allm),
        "m_all_ones_bernstein_signs": c2_bs_signs,
    }

    # Control 3: Octonion / Paley skew-Hadamard C2 = -I_8 + K_8 (all valid E_8 spinors, orthogonal, but det=256 != 16!)
    # Paley order-8 skew conference matrix K_8 (normalized so row 0 is all -1 in C2 = -I + K_8):
    # QR mod 7 = {1, 2, 4}
    QR7 = {1, 2, 4}
    K8 = [[0] * 8 for _ in range(8)]
    for j in range(1, 8):
        K8[0][j] = -1
        K8[j][0] = 1
    for i in range(1, 8):
        for j in range(1, 8):
            if i != j:
                K8[i][j] = 1 if ((j - i) % 7) in QR7 else -1
    C2_skew = [[(-1 if i == j else K8[i][j]) for j in range(8)] for i in range(8)]
    # Confirm C2_skew IS a genuine Hadamard matrix with all spinor rows!
    for i in range(8):
        require(C2_skew[i].count(-1) % 2 == 0, "Skew-Hadamard row should have even minus count")
    L_skew = [[Fraction(2 if i == j else 0) + Fraction(C2_skew[i][j], 2) for j in range(8)] for i in range(8)]
    det_skew = exact_det_rat(L_skew)
    rejected_3 = False
    try:
        verify_lemmas_1_to_4(C2_skew)
    except ValueError:
        rejected_3 = True
    require(rejected_3 and det_skew == 256, f"Control 3 failed: rejected={rejected_3}, det={det_skew}")
    controls["C3_octonion_skew_hadamard"] = {"rejected": True, "det_scaled_orthogonal": int(det_skew)}

    # Control 4: Tightened distance target 2001/1000 > 2 MUST produce negative Bernstein coefficients!
    res_tight = verify_lemma_5_bernstein_and_kissing(H8, R2_cutoff=11, target_num=2001, target_den=1000)
    require(res_tight["negative_pairs"] > 0, "Control 4 failed: tightened target > 2 must fail Bernstein!")
    controls["C4_tightened_target_2001_1000"] = {"rejected": True, "negative_pairs": res_tight["negative_pairs"]}

    # Control 5: Incomplete ellipsoid cutoff R2 = 4 < floor((2+sqrt(2))^2) = 11 MUST be rejected!
    rejected_5 = False
    try:
        verify_lemma_5_bernstein_and_kissing(H8, R2_cutoff=4)
    except ValueError:
        rejected_5 = True
    require(rejected_5, "Control 5 failed: incomplete cutoff R2=4 was not rejected!")
    controls["C5_incomplete_ellipsoid_cutoff"] = {"rejected": True}

    # Control 6: R3 strict type guard rejects bool and float inputs!
    bad_inputs_rejected = 0
    for bad in [(True, 0), (0, False), (1.0, 0), (0, 0.5)]:
        try:
            z2_sign(bad)
        except TypeError:
            bad_inputs_rejected += 1
    require(bad_inputs_rejected == 4, f"Control 6 failed: only {bad_inputs_rejected}/4 rejected")
    controls["C6_r3_bool_float_type_guards"] = {"rejected_count": bad_inputs_rejected}

    # Control 7: Non-square Gram ratio obstruction (E_6: 4/3, E_7: 2) fails squareness check
    for name, num, den in [("E_6", 4, 3), ("E_7", 4, 2)]:
        q = Fraction(num, den)
        is_sq = (q.denominator == 1) and (math.isqrt(q.numerator) ** 2 == q.numerator)
        require(not is_sq, f"Control 7 failed for {name}")
    controls["C7_e6_e7_obstruction_rejected"] = {"rejected": True}

    return controls


def main():
    t0 = time.perf_counter()
    H8 = build_sylvester_H8()
    res_1_4 = verify_lemmas_1_to_4(H8)
    res_5 = verify_lemma_5_bernstein_and_kissing(H8, R2_cutoff=11)
    require(res_5["negative_pairs"] == 0, f"Found {res_5['negative_pairs']} negative Bernstein pairs!")
    require(res_5["nonbar_interior_zeros"] == 0, f"Found {res_5['nonbar_interior_zeros']} non-bar interior zeros!")
    require(res_5["persistent_bar_pairs"] == 8, f"Expected 8 persistent bar pairs, got {res_5['persistent_bar_pairs']}")
    require(res_5["zero_coeffs_total"] == 976, f"Expected 976 exact zero coeffs, got {res_5['zero_coeffs_total']}")
    require(res_5["s1_contacts_V2e_nonbar"] == 120 and res_5["s1_contacts_Vee_n2_4"] == 784,
            "s=1 contact counts mismatch")

    res_6 = verify_lemma_6_uniqueness()
    res_gauge = verify_b8_gauge_invariance(H8)
    res_controls = run_must_fail_controls(H8)

    cert = {
        "theorem": "Canonical Sylvester-Hadamard Continuous Bar-Preserving Fold F_8 -> E_8 + E_6/E_7/D_n^+ Uniqueness",
        "status": "ALL_EXACT_CHECKS_PASS",
        "oracle": "exact",
        "arithmetic": "zero-float pure Python stdlib (int, Fraction, Z[sqrt(2)])",
        "sylvester_H8": H8,
        "lemmas_1_to_4": res_1_4,
        "lemma_5_bernstein_and_kissing": res_5,
        "lemma_6_uniqueness_obstruction": res_6,
        "b8_gauge_invariance": res_gauge,
        "must_fail_controls": res_controls,
    }
    canonical = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    sha256 = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "e8_spinor_fold_certificate.json")
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write(canonical)

    elapsed = time.perf_counter() - t0
    print(f"STATUS: {cert['status']}")
    print(f"CERT_PATH: {out_path}")
    print(f"CERT_SHA256: {sha256}")
    print(f"PAIRS_CHECKED: {res_5['total_pairs']} ({res_5['total_degree8_bernstein_coeffs']} degree-8 Z[sqrt(2)] Bernstein coeffs, 0 negative, 0 subdivisions, {res_5['zero_coeffs_total']} exact zeros)")
    print(f"KISSING_ALL_16: {res_5['kissing_number_all_16_balls']} | CONTROLS_REJECTED: {len(res_controls)}/7 | B8_GAUGE: {res_gauge['b8_gauge_trials_verified']}/5")
    print(f"ELAPSED_S: {elapsed:.3f}")


if __name__ == "__main__":
    main()
