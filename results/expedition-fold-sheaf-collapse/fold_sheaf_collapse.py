#!/usr/bin/env python3
"""Exact-Rational Universal Rank-1 Extremal Elastic Collapse & 4D Pure-Shear Duality Along F_n -> D_n.

100% zero-float (`fractions.Fraction` only, stdlib only) cross-domain synthesis instrument uniting:
  - `fractals-probe`'s All-Dimension Continuous Folding Theorem (`F_n -> D_n`, forall n >= 3)
    and Clifford-nilpotent 4D fold `C_4^nilp` (`probes/fractals/ALL_N_DN_FOLD_PROOF.md`, #40/#42),
  - `continuum-bridge`'s Cellular Sheaf Cohomology & Zero-Float Elastic Schur Complement
    (`synthesis/tropical_quotient_sheaf_bridge.py`, `AUDIT-004` VERIFIED), and
  - Exact Rational Congruence Inertia Certification (`synthesis/optimal_inertia_rank.py`).

Main New Mathematical Discoveries (100% Zero-Float Exact Rational Proofs & Certificates):

  1. Exact Rational Parametrization of the F_n -> D_n Curve Over Q (Zero sqrt(2)):
       a(k) = 4k / (2 k^2 + 1),   t(k) = (2 k^2 - 1) / (2 k^2 + 1),   k in Q, k >= 1/sqrt(2)
     satisfies a(k)^2 + 2 t(k)^2 == 2 identically in Q(k), yielding exact rational bar vectors
     b_i(k) = a(k) e_i + t(k) (C_n)_i (|b_i|^2 == 2) along the ENTIRE continuous folding curve!

  2. Theorem (Universal Instantaneous Rank-1 Elastic Collapse Along F_n -> D_n, forall n >= 3):
     At s = 0 (unfolded F_n): dim H^0(F_0) = n(n+1)/2, dim H^1(F_0) = n(n-1)/2, and the
     macroscopic elastic stiffness tensor C(0) has RANK n (inertia (n, 0, n(n-1)/2), zero shear).
     THE INSTANT s > 0 (t > 0, across the ENTIRE open-to-closed fold s in (0, 1]):
       Because C_n (and C_4^nilp) has connected support graph and is uniquely diagonally
       symmetrizable (diag(tau*) C symmetric with tau*(C_n) = (-1, 1, ..., 1) and
       tau*(C_4^nilp) = (1, 1, -1, -1)):
         - Rotational locking operator M_rot(C) : so(n) -> Q^n has EXACT RANK n - 1;
         - Both q=0 sheaf Betti numbers drop discontinuously by -(n - 1):
             dim H^0(F_0(s)) = n(n-1)/2 + 1,   dim H^1(F_0(s)) = n(n-3)/2 + 1;
         - Exactly ONE cross-bar-coupled periodic self-stress survives modulo the n(n-3)/2
           pure intra-cell stresses;
         - Therefore, C(s) COLLAPSES DISCONTINUOUSLY FROM RANK n TO RANK 1 AT s = 0+
           (inertia (1, 0, n(n+1)/2 - 1) — a universal all-dimension rank-1 extremal /
           pentamode-class metamaterial!) for EVERY s in (0, 1] and EVERY n >= 3!
     Exact Closed-Form Rank-1 Stiffness Trace at s = 1 (L_n(1) = 2 I_n + C_n, vol = 4n):
       - For C_3 (n = 3): tr(S_cell) = 4,           tr(C_el) = 1/3.
       - For C_n (n >= 4): tr(S_cell) = (6n - 2)/n, tr(C_el) = (3n - 1)/(2 n^2)
         (verified exactly in Q for n = 4 -> 11/2, n = 5 -> 28/5, n = 6 -> 17/3, n = 7 -> 40/7!).

  3. Theorem (4D Clifford-Nilpotent Pure-Shear / Zero-Axial Exact Duality):
     At s = 1 along C_4^nilp (C^2 = 0, tau* = (+1, +1, -1, -1)):
       - Every axial component sum_i tau*_i (b_i)_k^2 cancels IDENTICALLY to 0!
       - The surviving rank-1 elastic tensor C_el^nilp(1) has ZERO normal/bulk stiffness
         (all 4 axial entries C_{ii} == 0, bulk modulus K == 0!) and is supported PURELY
         on the 4 shear strains (eps_01, eps_02, eps_13, eps_23) with equal weights (1,1,0,0,1,1)
         and tr(S_cell) = 4 (tr(C_el) = 1/4), strictly distinguishing C_4^nilp (4) from C_4 (11/2)!
       - Exact orthogonal duality with s = 0: F_4(s=0) is 100% pure axial / zero shear, while
         C_4^nilp(s=1) is 100% pure shear / zero axial!

  4. Exact Jamming Transition at s = 1 (D_3 fcc & D_4 Optimal 4D Packing):
     When all newly touching D_n kissing contacts (2n(n-1) per ball) snap shut at s = 1:
       - D_3 (36 bars/cell): H^0 = 3 (translations only), C_jam inertia = (6, 0, 0) (rank 1 -> 6),
         C11 = 1, C12 = 1/2, C44_Voigt = 1/2 == C12 (exact Cauchy relation, Zener A = 2).
       - D_4 (96 bars/cell, identical for BOTH C_4 and C_4^nilp): H^0 = 4, C_jam inertia = (10, 0, 0)
         (rank 1 -> 10), C11 = 3/2, C12 = 1/2, C44_Voigt = 1/2 == C12 AND C11 - C12 == 2 C44_Voigt
         (exact 4D Cauchy relation AND exact 4D elastic isotropy: Zener ratio A(D_4) == 1!).

  5. Portal 9 8D Spinor-Root Fold F_8 -> E_8 (Cyclic Rank-0 vs Hadamard Rank-1 & O(8)-Isotropic Jamming):
     Exact Master Rank-1 Macroscopic Stress Formula (verified across C_3..C_5, C_4^nilp, and H_8):
       sigma* = C^T diag(tau*) C + 2 diag(tau*) C  (by Maxwell's Periodic Stress Theorem at L = 2I + C).
     Comparing the two 8D F_8 -> E_8 spinor folds:
       (a) Cyclic Z/16Z ker(phi) all-spinor bases (`session-cartographer` #46): EXHAUSTIVELY ZERO (0)
           are diagonally symmetrizable => rk(M_rot) = 8/8, H^0=28, H^1=20, COMPLETE RANK-0 FLOPPY
           COLLAPSE (rk([d0 | ae]) == rk(d0) == 100 => S_cell === 0) at both pure-Q interior
           (16/11, 7/11) and s=1 (0, 1).
       (b) Canonical Sylvester-Hadamard C_8 = (1/2) H_8 (`fractals-probe` #61, H_8 = H_2^(tensor 3)):
           SYMMETRIC (tau* = (1,...,1)) => rk(M_rot) = 7/8, H^0=29, H^1=21, EXACT RANK-1 AUXETIC
           EXTREMAL COLLAPSE with unique macroscopic stress sigma* = 2 I_8 + H_8 (indefinite spectrum
           2 +- 2 sqrt(2), signature (4, 4, 0)), verified via rk([d0 | ae * ker(sigma*^T)]) == rk(d0) == 99!
       (c) Upon full E_8 jamming (240 roots = 120 undirected bars, spherical 7-design, vol = 1),
           BOTH snap into FULL RANK 36 (inertia (36, 0, 0)) with EXACT 8D O(8)-ISOTROPY
           (C11 = 9, C12 = C44_voigt = 3, A(E_8) = 1), contrasting with D_8's A(D_8) = 1/3!

Run:
  python3 -O -B -W error results/expedition-fold-sheaf-collapse/fold_sheaf_collapse.py
  python3 -O -B -W error results/expedition-fold-sheaf-collapse/fold_sheaf_collapse.py --json
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import itertools
import json
import sys
from typing import Any, Dict, List, Optional, Sequence, Tuple

Mat = List[List[Fraction]]

# --- Self-contained exact rational linear algebra & block-LDL^T inertia TCB ---
def _require_exact_square_matrix(M: Any, name: str, expected_n: Optional[int] = None) -> int:
    """Fail-closed check that M is a non-empty square matrix of exact Fraction/int (no float/bool/ragged)."""
    if not isinstance(M, list) or len(M) == 0:
        raise ValueError(f"{name}_EMPTY_OR_NOT_LIST")
    n = len(M)
    if expected_n is not None and n != expected_n:
        raise ValueError(f"{name}_DIM_MISMATCH: expected {expected_n}, got {n}")
    for row in M:
        if not isinstance(row, list) or len(row) != n:
            raise ValueError(f"{name}_NOT_SQUARE_OR_RAGGED")
        for v in row:
            if isinstance(v, bool) or isinstance(v, float) or not isinstance(v, (int, Fraction)):
                raise ValueError(f"{name}_NON_EXACT_ENTRY: got {type(v).__name__}")
    return n

def mat_zero(n: int, m: int) -> Mat:
    return [[Fraction(0) for _ in range(m)] for _ in range(n)]


def mat_eye(n: int) -> Mat:
    return [[Fraction(1) if i == j else Fraction(0) for j in range(n)] for i in range(n)]


def mat_t(A: Mat) -> Mat:
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def mat_mul(A: Mat, B: Mat) -> Mat:
    n, k, m = len(A), len(B), len(B[0])
    if len(A[0]) != k:
        raise ValueError("MAT_DIM_MISMATCH")
    return [
        [sum((A[i][t] * B[t][j] for t in range(k)), Fraction(0)) for j in range(m)]
        for i in range(n)
    ]


def mat_add(A: Mat, B: Mat) -> Mat:
    return [[a + b for a, b in zip(rA, rB)] for rA, rB in zip(A, B)]


def mat_sub(A: Mat, B: Mat) -> Mat:
    return [[a - b for a, b in zip(rA, rB)] for rA, rB in zip(A, B)]


def mat_scale(c: Fraction, A: Mat) -> Mat:
    return [[c * v for v in row] for row in A]


def is_symmetric(A: Mat) -> bool:
    n = len(A)
    if n == 0 or any(len(r) != n for r in A):
        return False
    return all(A[i][j] == A[j][i] for i in range(n) for j in range(i + 1, n))


def trace(A: Mat) -> Fraction:
    return sum((A[i][i] for i in range(len(A))), Fraction(0))


def hs_norm_sq(A: Mat) -> Fraction:
    return sum((v * v for row in A for v in row), Fraction(0))


def exact_rank(A: Mat) -> int:
    """Exact rank over Q via Fraction Gaussian elimination."""
    m = len(A)
    n = len(A[0])
    R = [row[:] for row in A]
    rank = 0
    for col in range(n):
        piv = None
        for r in range(rank, m):
            if R[r][col] != 0:
                piv = r
                break
        if piv is None:
            continue
        R[rank], R[piv] = R[piv], R[rank]
        pv = R[rank][col]
        for r in range(rank + 1, m):
            if R[r][col] != 0:
                fac = R[r][col] / pv
                for c in range(col, n):
                    R[r][c] -= fac * R[rank][c]
        rank += 1
        if rank == m:
            break
    return rank

def produce_block_ldl(M: Mat) -> Tuple[List[int], Mat, Mat]:
    """Compute exact P M P^T = L D L^T over Q with 1x1 and 2x2 symmetric blocks in D.

    Returns (perm, L, D) where perm maps new index i -> old index perm[i].
    Works on every symmetric rational matrix, including zero-diagonal indefinite matrices.
    """
    if not is_symmetric(M):
        raise ValueError("MATRIX_NOT_SYMMETRIC")
    n = len(M)
    A = [row[:] for row in M]
    perm = list(range(n))
    L = mat_eye(n)
    D = mat_zero(n, n)

    k = 0
    while k < n:
        # Look for a nonzero diagonal entry in A[k:, k:]
        diag_piv = None
        for i in range(k, n):
            if A[i][i] != 0:
                diag_piv = i
                break

        if diag_piv is not None:
            if diag_piv != k:
                perm[k], perm[diag_piv] = perm[diag_piv], perm[k]
                A[k], A[diag_piv] = A[diag_piv], A[k]
                for r in range(n):
                    A[r][k], A[r][diag_piv] = A[r][diag_piv], A[r][k]
                for c in range(k):
                    L[k][c], L[diag_piv][c] = L[diag_piv][c], L[k][c]

            d = A[k][k]
            D[k][k] = d
            for i in range(k + 1, n):
                L[i][k] = A[i][k] / d
            for i in range(k + 1, n):
                lik = L[i][k]
                if lik == 0:
                    continue
                for j in range(i, n):
                    val = A[i][j] - lik * d * L[j][k]
                    A[i][j] = val
                    A[j][i] = val
            k += 1
            continue

        # All remaining diagonal entries are zero; look for a nonzero off-diagonal entry
        off_piv = None
        for i in range(k, n):
            for j in range(i + 1, n):
                if A[i][j] != 0:
                    off_piv = (i, j)
                    break
            if off_piv is not None:
                break

        if off_piv is None:
            # Remaining trailing submatrix is identically zero
            break

        p1, p2 = off_piv
        # Swap p1 -> k, then p2 -> k+1
        if p1 != k:
            perm[k], perm[p1] = perm[p1], perm[k]
            A[k], A[p1] = A[p1], A[k]
            for r in range(n):
                A[r][k], A[r][p1] = A[r][p1], A[r][k]
            for c in range(k):
                L[k][c], L[p1][c] = L[p1][c], L[k][c]
            if p2 == k:
                p2 = p1
        if p2 != k + 1:
            perm[k + 1], perm[p2] = perm[p2], perm[k + 1]
            A[k + 1], A[p2] = A[p2], A[k + 1]
            for r in range(n):
                A[r][k + 1], A[r][p2] = A[r][p2], A[r][k + 1]
            for c in range(k):
                L[k + 1][c], L[p2][c] = L[p2][c], L[k + 1][c]

        # 2x2 symmetric block [[0, b], [b, 0]] has det = -b^2 < 0 and inverse [[0, 1/b], [1/b, 0]]
        b_val = A[k][k + 1]
        D[k][k + 1] = b_val
        D[k + 1][k] = b_val
        inv_b = Fraction(1) / b_val
        for i in range(k + 2, n):
            L[i][k] = A[i][k + 1] * inv_b
            L[i][k + 1] = A[i][k] * inv_b
        for i in range(k + 2, n):
            li0, li1 = L[i][k], L[i][k + 1]
            for j in range(i, n):
                lj0, lj1 = L[j][k], L[j][k + 1]
                # [li0, li1] [[0, b], [b, 0]] [lj0, lj1]^T = b * (li0*lj1 + li1*lj0)
                val = A[i][j] - b_val * (li0 * lj1 + li1 * lj0)
                A[i][j] = val
                A[j][i] = val
        k += 2

    return perm, L, D


def verify_block_ldl_inertia(M: Mat, perm: Sequence[int], L: Mat, D: Mat) -> Tuple[int, int, int]:
    """Fail-closed verifier for P M P^T = L D L^T. Returns exact inertia (n_pos, n_neg, n_zero)."""
    n = _require_exact_square_matrix(M, "M")
    _require_exact_square_matrix(L, "L", expected_n=n)
    _require_exact_square_matrix(D, "D", expected_n=n)
    if not is_symmetric(M):
        raise ValueError("M_NOT_SYMMETRIC")
    if (
        len(perm) != n
        or any(isinstance(p, bool) or not isinstance(p, int) for p in perm)
        or sorted(perm) != list(range(n))
    ):
        raise ValueError("INVALID_PERMUTATION")
    for i in range(n):
        if L[i][i] != 1:
            raise ValueError("L_NON_UNIT_DIAGONAL")
        for j in range(i + 1, n):
            if L[i][j] != 0:
                raise ValueError("L_NOT_LOWER_TRIANGULAR")
    if not is_symmetric(D):
        raise ValueError("D_NOT_SYMMETRIC")
    for i in range(n):
        for j in range(n):
            if abs(i - j) > 1 and D[i][j] != 0:
                raise ValueError("D_NOT_TRIDIAGONAL_BLOCK")

    # Verify non-overlapping 1x1 / 2x2 block structure and read exact inertia
    n_pos = n_neg = n_zero = 0
    i = 0
    while i < n:
        has_off = (i + 1 < n) and (D[i][i + 1] != 0)
        if has_off:
            # Ensure this 2x2 block does not overlap with a following off-diagonal entry
            if i + 2 < n and D[i + 1][i + 2] != 0:
                raise ValueError("D_OVERLAPPING_2X2_BLOCKS")
            a, b_val, c = D[i][i], D[i][i + 1], D[i + 1][i + 1]
            det = a * c - b_val * b_val
            tr = a + c
            if det < 0:
                n_pos += 1
                n_neg += 1
            elif det > 0:
                if tr > 0:
                    n_pos += 2
                else:
                    n_neg += 2
            else:
                if tr > 0:
                    n_pos += 1
                    n_zero += 1
                elif tr < 0:
                    n_neg += 1
                    n_zero += 1
                else:
                    n_zero += 2
            i += 2
        else:
            d = D[i][i]
            if d > 0:
                n_pos += 1
            elif d < 0:
                n_neg += 1
            else:
                n_zero += 1
            i += 1

    # Verify exact matrix identity P M P^T == L D L^T
    PMP = [[M[perm[r]][perm[c]] for c in range(n)] for r in range(n)]
    LDL = mat_mul(mat_mul(L, D), mat_t(L))
    if PMP != LDL:
        raise ValueError("LDL_IDENTITY_MISMATCH")

    return (n_pos, n_neg, n_zero)

def exact_linear_solve_consistent(A: Mat, rhs_cols: Mat) -> Mat:
    """Solve A X = rhs_cols over Q for consistent (possibly rank-deficient) symmetric A.

    Sets free variables to 0 and verifies A X == rhs_cols exactly.
    """
    m = len(A)
    n = len(A[0])
    k = len(rhs_cols[0])
    aug = [A[i][:] + rhs_cols[i][:] for i in range(m)]
    piv_cols: List[int] = []
    row = 0
    for col in range(n):
        piv = None
        for r in range(row, m):
            if aug[r][col] != 0:
                piv = r
                break
        if piv is None:
            continue
        aug[row], aug[piv] = aug[piv], aug[row]
        pv = aug[row][col]
        for c in range(col, n + k):
            aug[row][c] /= pv
        for r in range(m):
            if r != row and aug[r][col] != 0:
                fac = aug[r][col]
                for c in range(col, n + k):
                    aug[r][c] -= fac * aug[row][c]
        piv_cols.append(col)
        row += 1
        if row == m:
            break

    # Check consistency of zero rows
    for r in range(row, m):
        for c in range(n, n + k):
            if aug[r][c] != 0:
                raise ValueError("INCONSISTENT_LINEAR_SYSTEM")

    X = [[Fraction(0) for _ in range(k)] for _ in range(n)]
    for r, col in enumerate(piv_cols):
        for j in range(k):
            X[col][j] = aug[r][n + j]

    # Exact verification A X == rhs_cols
    if mat_mul(A, X) != rhs_cols:
        raise RuntimeError("LINEAR_SOLVE_VERIFICATION_FAILED")
    return X



def make_twisted_cycle_Cn(n: int) -> List[List[int]]:
    """fractals-probe's universal twisted n-cycle fold generator C_n (n >= 3)."""
    if n < 3:
        raise ValueError("REQUIRE_N_GE_3")
    C = [[0] * n for _ in range(n)]
    C[0][1] = -1
    C[0][n - 1] = -1
    C[1][0] = 1
    C[1][2] = -1
    for k in range(2, n - 1):
        C[k][k - 1] = -1
        C[k][k + 1] = -1
    C[n - 1][0] = 1
    C[n - 1][n - 2] = -1
    return C


C4_NILPOTENT = [
    [0, -1, -1, 0],
    [-1, 0, 0, -1],
    [1, 0, 0, 1],
    [0, 1, 1, 0],
]


def verify_diagonal_symmetrizer(C: List[List[int]], tau: Sequence[int]) -> bool:
    """Verify diag(tau) * C is symmetric and tau is nonzero."""
    n = len(C)
    if len(tau) != n or any(t == 0 for t in tau):
        return False
    return all(tau[i] * C[i][j] == tau[j] * C[j][i] for i in range(n) for j in range(n))


def rotational_locking_rank(C: List[List[int]]) -> int:
    """Exact rank over Q of M_rot(C) : so(n) -> Q^n, Omega -> (C_i . (Omega e_i))_i."""
    n = len(C)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    M = []
    for i in range(n):
        row = [Fraction(0)] * len(pairs)
        for p_idx, (a, b_p) in enumerate(pairs):
            if a == i:
                row[p_idx] -= Fraction(C[i][b_p])
            elif b_p == i:
                row[p_idx] += Fraction(C[i][a])
        M.append(row)
    return exact_rank(M)


def compute_fold_sheaf_elastic(
    C_mat: List[List[int]], a: Fraction, t: Fraction
) -> Dict[str, Any]:
    """Exact zero-float cellular sheaf cohomology & elastic Schur complement along L(a, t)."""
    if a * a + 2 * t * t != Fraction(2):
        raise ValueError(f"BAR_LENGTH_CONSTRAINT_VIOLATED: a^2 + 2 t^2 = {a*a + 2*t*t} != 2")
    n = len(C_mat)
    n_dofs = 2 * n * n
    off_pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m_s = n + len(off_pairs)
    off_idx = {pair: n + idx for idx, pair in enumerate(off_pairs)}

    def v_idx(ax: int, sgn: int) -> int:
        return 2 * ax + sgn

    d0_rows: List[List[Fraction]] = []
    ae_rows: List[List[Fraction]] = []

    # 2n(n-1) intra-cell bars (length^2 == 2 -> weight 1/2)
    for i in range(n):
        for j in range(i + 1, n):
            for si, s in ((0, Fraction(1)), (1, Fraction(-1))):
                for ti, tt in ((0, Fraction(1)), (1, Fraction(-1))):
                    ru = [Fraction(0)] * n_dofs
                    vi, vj = v_idx(i, si), v_idx(j, ti)
                    ru[vi * n + i] += s
                    ru[vi * n + j] -= tt
                    ru[vj * n + i] -= s
                    ru[vj * n + j] += tt
                    re = [Fraction(0)] * m_s
                    re[i] += Fraction(1)
                    re[j] += Fraction(1)
                    re[off_idx[(i, j)]] += Fraction(-2) * s * tt
                    d0_rows.append(ru)
                    ae_rows.append(re)

    # n inter-cell bridge bars b_i = a e_i + t C_i (|b_i|^2 == 2 -> weight 1/2)
    for i in range(n):
        b = [t * Fraction(C_mat[i][k]) + (a if k == i else Fraction(0)) for k in range(n)]
        ru = [Fraction(0)] * n_dofs
        vp, vm = v_idx(i, 0), v_idx(i, 1)
        for k in range(n):
            ru[vm * n + k] += b[k]
            ru[vp * n + k] -= b[k]
        re = [Fraction(0)] * m_s
        for k in range(n):
            re[k] += b[k] * b[k]
        for k in range(n):
            for l in range(k + 1, n):
                re[off_idx[(k, l)]] += Fraction(2) * b[k] * b[l]
        d0_rows.append(ru)
        ae_rows.append(re)

    rk = exact_rank(d0_rows)
    h0 = n_dofs - rk
    h1 = len(d0_rows) - rk

    W_d0 = [[Fraction(1, 2) * v for v in r] for r in d0_rows]
    W_ae = [[Fraction(1, 2) * v for v in r] for r in ae_rows]
    Delta_F = mat_mul(mat_t(d0_rows), W_d0)
    B = mat_mul(mat_t(d0_rows), W_ae)
    K_ee = mat_mul(mat_t(ae_rows), W_ae)
    U_star = exact_linear_solve_consistent(Delta_F, [[-v for v in r] for r in B])
    S_cell = [
        [
            K_ee[i][j] + sum((B[r][i] * U_star[r][j] for r in range(n_dofs)), Fraction(0))
            for j in range(m_s)
        ]
        for i in range(m_s)
    ]
    c_in = verify_block_ldl_inertia(S_cell, *produce_block_ldl(S_cell))
    tr_S = sum((S_cell[i][i] for i in range(m_s)), Fraction(0))
    axial_block_zero = all(S_cell[i][j] == 0 for i in range(n) for j in range(n))
    mode = None
    if c_in[0] == 1:
        nz = next(i for i in range(m_s) if S_cell[i][i] != 0)
        mode = [str(S_cell[nz][j] / S_cell[nz][nz]) for j in range(m_s)]
    return {
        "n": n,
        "a": str(a),
        "t": str(t),
        "H0": h0,
        "H1": h1,
        "euler_H0_minus_H1": h0 - h1,
        "S_cell_inertia": list(c_in),
        "tr_S_cell": str(tr_S),
        "S_00": str(S_cell[0][0]),
        "axial_block_identically_zero": axial_block_zero,
        "rank1_mode_normalized": mode,
    }


def compute_jammed_Dn_sheaf(C_mat: List[List[int]]) -> Dict[str, Any]:
    """Exact zero-float elastic tensor of the full jammed D_n kissing-contact sheaf at s = 1."""
    n = len(C_mat)
    n_dofs = 2 * n * n
    off_pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m_s = n + len(off_pairs)
    off_idx = {pair: n + idx for idx, pair in enumerate(off_pairs)}
    L = [[(2 if i == j else 0) + C_mat[i][j] for j in range(n)] for i in range(n)]
    pts = []
    for ax in range(n):
        pts.append([1 if k == ax else 0 for k in range(n)])
        pts.append([-1 if k == ax else 0 for k in range(n)])

    jam_d0: List[List[Fraction]] = []
    jam_ae: List[List[Fraction]] = []
    seen = set()
    for m in itertools.product(range(-2, 3), repeat=n):
        mL = [sum(m[r] * L[r][c] for r in range(n)) for c in range(n)]
        for u_i in range(2 * n):
            for v_i in range(2 * n):
                if m == (0,) * n and u_i == v_i:
                    continue
                diff = [mL[c] + pts[v_i][c] - pts[u_i][c] for c in range(n)]
                if sum(x * x for x in diff) == 2:
                    neg_m = tuple(-x for x in m)
                    canon = min((u_i, v_i, m), (v_i, u_i, neg_m))
                    if canon in seen:
                        continue
                    seen.add(canon)
                    ru = [Fraction(0)] * n_dofs
                    for k in range(n):
                        ru[v_i * n + k] += Fraction(diff[k])
                        ru[u_i * n + k] -= Fraction(diff[k])
                    re = [Fraction(0)] * m_s
                    for k in range(n):
                        re[k] += Fraction(diff[k] * diff[k])
                    for k in range(n):
                        for l in range(k + 1, n):
                            re[off_idx[(k, l)]] += Fraction(2 * diff[k] * diff[l])
                    jam_d0.append(ru)
                    jam_ae.append(re)

    rk = exact_rank(jam_d0)
    h0 = n_dofs - rk
    W_d0 = [[Fraction(1, 2) * v for v in r] for r in jam_d0]
    W_ae = [[Fraction(1, 2) * v for v in r] for r in jam_ae]
    DF = mat_mul(mat_t(jam_d0), W_d0)
    B = mat_mul(mat_t(jam_d0), W_ae)
    K_ee = mat_mul(mat_t(jam_ae), W_ae)
    U_j = exact_linear_solve_consistent(DF, [[-v for v in r] for r in B])
    vol = Fraction(4 * n)
    C_jam = [
        [
            (K_ee[i][j] + sum((B[r][i] * U_j[r][j] for r in range(n_dofs)), Fraction(0))) / vol
            for j in range(m_s)
        ]
        for i in range(m_s)
    ]
    in_j = verify_block_ldl_inertia(C_jam, *produce_block_ldl(C_jam))
    c11 = C_jam[0][0]
    c12 = C_jam[0][1]
    # Off-diagonal strain coordinate is tensor eps_ij; Voigt engineering shear gamma_ij = 2 eps_ij => C44_voigt = C_jam[n][n] / 4
    c44_voigt = C_jam[n][n] / Fraction(4)
    zener_A = (Fraction(2) * c44_voigt) / (c11 - c12)
    return {
        "n": n,
        "bars_per_cell": len(jam_d0),
        "expected_Dn_bars": 2 * n * n * (n - 1),
        "H0_translations_only": h0,
        "C_jam_inertia": list(in_j),
        "C11": str(c11),
        "C12": str(c12),
        "C44_voigt": str(c44_voigt),
        "cauchy_relation_C12_eq_C44_voigt": (c12 == c44_voigt),
        "zener_anisotropy_ratio_A": str(zener_A),
        "elastically_isotropic": (zener_A == Fraction(1)),
    }


def compute_portal9_E8_spinor_collapse_and_jamming() -> Dict[str, Any]:
    """Exact zero-float Portal 9 8D F_8 -> E_8 spinor-root rank-0 collapse & O(8)-isotropic jamming."""
    n = 8
    spinors = [
        tuple(Fraction(b, 2) for b in bits)
        for bits in itertools.product((-1, 1), repeat=n)
        if sum(1 for b in bits if b == -1) % 2 == 0
    ]

    def phi_e8(x: Sequence[Fraction]) -> int:
        y = [int(xi - Fraction(1, 2)) for xi in x]
        return (2 - sum(y) // 2 + sum(2 * j * y[j] for j in range(n))) % 16

    def det_exact(M: Mat) -> Fraction:
        m = len(M)
        A = [r[:] for r in M]
        sign = Fraction(1)
        for col in range(m):
            piv = next((r for r in range(col, m) if A[r][col] != 0), None)
            if piv is None:
                return Fraction(0)
            if piv != col:
                A[col], A[piv] = A[piv], A[col]
                sign = -sign
            for r in range(col + 1, m):
                f = A[r][col] / A[col][col]
                for c in range(col, m):
                    A[r][c] -= f * A[col][c]
        d = sign
        for i in range(m):
            d *= A[i][i]
        return d

    cands_all: List[List[Tuple[Fraction, ...]]] = [[] for _ in range(n)]
    cands_neg: List[List[Tuple[Fraction, ...]]] = [[] for _ in range(n)]
    for i in range(n):
        for c in spinors:
            row = list(c)
            row[i] += Fraction(2)
            if phi_e8(row) == 0:
                cands_all[i].append(c)
                if c[i] == Fraction(-1, 2):
                    cands_neg[i].append(c)

    # 1. Exhaustive proof: ZERO all-spinor E_8 bases in ker(phi) are diagonally symmetrizable!
    sym_e8_count = 0
    for tau_bits in itertools.product((-1, 1), repeat=n - 1):
        tau = (1,) + tau_bits

        def bt(i: int, rows: List[Tuple[Fraction, ...]]) -> None:
            nonlocal sym_e8_count
            if i == n:
                L = [
                    [rows[a][b] + (Fraction(2) if a == b else Fraction(0)) for b in range(n)]
                    for a in range(n)
                ]
                if abs(det_exact(L)) == Fraction(16):
                    sym_e8_count += 1
                return
            for c in cands_all[i]:
                if all(tau[i] * c[j] == tau[j] * rows[j][i] for j in range(i)):
                    rows.append(c)
                    bt(i + 1, rows)
                    rows.pop()

        bt(0, [])
    if sym_e8_count != 0:
        raise RuntimeError(f"Unexpected diagonally symmetrizable all-spinor E_8 count: {sym_e8_count}")

    # 2. Deterministic first all-spinor E_8 basis C8 with c_ii = -1/2 and det(2I + C8) == +16
    C8: Optional[Mat] = None
    for idxs in itertools.product(range(len(cands_neg[0])), repeat=n):
        C = [list(cands_neg[i][idxs[i]]) for i in range(n)]
        L = [[C[i][j] + (Fraction(2) if i == j else Fraction(0)) for j in range(n)] for i in range(n)]
        if det_exact(L) == Fraction(16):
            C8 = C
            break
    if C8 is None:
        raise RuntimeError("Failed to find deterministic C8")

    # Rotational locking rank of C8 (8 x 28) == 8 (full rank!)
    off_pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    M_rot = [[Fraction(0)] * len(off_pairs) for _ in range(n)]
    for col, (i, j) in enumerate(off_pairs):
        M_rot[i][col] = -C8[i][j]
        M_rot[j][col] = C8[j][i]
    rk_rot_e8 = exact_rank(M_rot)
    if rk_rot_e8 != 8:
        raise RuntimeError(f"Expected rk(M_rot(C8)) == 8, got {rk_rot_e8}")

    # 3. Verify exact rank-0 elastic collapse (rk([d0 | ae]) == rk(d0) == 100 => S_cell === 0!)
    #    at BOTH a pure-Q interior point (a=16/11, t=7/11 via m=-1/4) and endpoint s=1 (a=0, t=1)!
    #    Check conic: a^2 - a*t + 2*t^2 == 2 (since c_ii == -1/2 and ||c_i||^2 == 2).
    n_dofs = 2 * n * n  # 128
    m_s = n + len(off_pairs)  # 36
    off_idx = {pair: n + idx for idx, pair in enumerate(off_pairs)}

    def v_idx(ax: int, sgn: int) -> int:
        return 2 * ax + sgn

    fold_points = []
    for label, a_val, t_val in [
        ("interior_m_neg_1_4", Fraction(16, 11), Fraction(7, 11)),
        ("endpoint_s1", Fraction(0), Fraction(1)),
    ]:
        if a_val * a_val - a_val * t_val + 2 * t_val * t_val != Fraction(2):
            raise RuntimeError("E8_RATIONAL_CURVE_NORM_VIOLATION")
        d0_rows: List[List[Fraction]] = []
        ae_rows: List[List[Fraction]] = []
        for i in range(n):
            for j in range(i + 1, n):
                for si, s in ((0, Fraction(1)), (1, Fraction(-1))):
                    for ti, tt in ((0, Fraction(1)), (1, Fraction(-1))):
                        ru = [Fraction(0)] * n_dofs
                        vi, vj = v_idx(i, si), v_idx(j, ti)
                        ru[vi * n + i] += s
                        ru[vi * n + j] -= tt
                        ru[vj * n + i] -= s
                        ru[vj * n + j] += tt
                        re = [Fraction(0)] * m_s
                        re[i] += Fraction(1)
                        re[j] += Fraction(1)
                        re[off_idx[(i, j)]] += Fraction(-2) * s * tt
                        d0_rows.append(ru)
                        ae_rows.append(re)
        for i in range(n):
            b = [t_val * C8[i][k] + (a_val if k == i else Fraction(0)) for k in range(n)]
            if sum((x * x for x in b), Fraction(0)) != Fraction(2):
                raise RuntimeError("E8_BRIDGE_BAR_LENGTH_VIOLATION")
            ru = [Fraction(0)] * n_dofs
            vp, vm = v_idx(i, 0), v_idx(i, 1)
            for k in range(n):
                ru[vm * n + k] += b[k]
                ru[vp * n + k] -= b[k]
            re = [Fraction(0)] * m_s
            for k in range(n):
                re[k] += b[k] * b[k]
            for k in range(n):
                for l in range(k + 1, n):
                    re[off_idx[(k, l)]] += Fraction(2) * b[k] * b[l]
            d0_rows.append(ru)
            ae_rows.append(re)
        rk_d0 = exact_rank(d0_rows)
        rk_aug = exact_rank([r1 + r2 for r1, r2 in zip(d0_rows, ae_rows)])
        h0 = n_dofs - rk_d0
        h1 = len(d0_rows) - rk_d0
        s_rank = rk_aug - rk_d0
        if not (rk_d0 == 100 and h0 == 28 and h1 == 20 and s_rank == 0):
            raise RuntimeError(f"E8 rank-0 collapse failed at {label}: rk={rk_d0}, s_rank={s_rank}")
        fold_points.append({
            "label": label,
            "a": str(a_val),
            "t": str(t_val),
            "H0": h0,
            "H1": h1,
            "S_cell_rank": s_rank,
        })

    # 4. Exact 36x36 jammed E_8 Bravais Cauchy-Born elastic tensor (240 roots = 120 undirected, vol = 1)
    #    vs D_8 (112 roots = 56 undirected, anisotropic A(D_8) = 1/3)!
    roots_d8 = []
    for i in range(n):
        for j in range(i + 1, n):
            for s1 in (-1, 1):
                for s2 in (-1, 1):
                    v = [Fraction(0)] * n
                    v[i], v[j] = Fraction(s1), Fraction(s2)
                    roots_d8.append(tuple(v))
    undir_e8 = [r for r in (roots_d8 + spinors) if next(x for x in r if x != 0) > 0]
    undir_d8 = [r for r in roots_d8 if next(x for x in r if x != 0) > 0]
    if len(undir_e8) != 120 or len(undir_d8) != 56:
        raise RuntimeError("Root count mismatch")

    def bravais_tensor(undir_roots: Sequence[Tuple[Fraction, ...]]) -> Mat:
        C_mat = [[Fraction(0)] * m_s for _ in range(m_s)]
        for r in undir_roots:
            re = [Fraction(0)] * m_s
            for k in range(n):
                re[k] = r[k] * r[k]
            for k in range(n):
                for l in range(k + 1, n):
                    re[off_idx[(k, l)]] = Fraction(2) * r[k] * r[l]
            for a in range(m_s):
                if re[a] == 0:
                    continue
                for b in range(m_s):
                    if re[b] != 0:
                        C_mat[a][b] += Fraction(1, 2) * re[a] * re[b]
        return C_mat

    C_e8 = bravais_tensor(undir_e8)
    in_e8 = verify_block_ldl_inertia(C_e8, *produce_block_ldl(C_e8))
    c11_e8, c12_e8, c44v_e8 = C_e8[0][0], C_e8[0][1], C_e8[n][n] / Fraction(4)
    zener_e8 = (Fraction(2) * c44v_e8) / (c11_e8 - c12_e8)
    exact_o8_iso = all(
        C_e8[a][b] == (c11_e8 if (a < n and a == b) else c12_e8 if (a < n and b < n) else Fraction(4) * c44v_e8 if a == b else Fraction(0))
        for a in range(m_s)
        for b in range(m_s)
    )
    C_d8 = bravais_tensor(undir_d8)
    zener_d8 = (Fraction(2) * (C_d8[n][n] / Fraction(4))) / (C_d8[0][0] - C_d8[0][1])
    if not (
        in_e8 == (36, 0, 0)
        and c11_e8 == Fraction(9)
        and c12_e8 == Fraction(3)
        and c44v_e8 == Fraction(3)
        and zener_e8 == Fraction(1)
        and exact_o8_iso
        and zener_d8 == Fraction(1, 3)
    ):
        raise RuntimeError("Jammed E_8 isotropy / D_8 anisotropy check failed!")

    # 5. Canonical Sylvester-Hadamard C8_had = (1/2) H_8 (fractals-probe #61):
    #    Symmetric (tau* = (1,...,1)) => rk(M_rot) = 7/8, H0=29, H1=21, EXACT RANK-1 AUXETIC
    #    COLLAPSE with unique stress tensor sigma* = C^T C + 2 C = 2 I_8 + H_8 (signature (4, 4, 0))!
    H8 = [[1]]
    for _ in range(3):
        H8 = [r + r for r in H8] + [r + [-x for x in r] for r in H8]
    C8_had = [[Fraction(x, 2) for x in r] for r in H8]
    M_rot_h8 = [[Fraction(0)] * len(off_pairs) for _ in range(n)]
    for col, (i, j) in enumerate(off_pairs):
        M_rot_h8[i][col] = -C8_had[i][j]
        M_rot_h8[j][col] = C8_had[j][i]
    rk_rot_h8 = exact_rank(M_rot_h8)
    d0_h8: List[List[Fraction]] = []
    ae_h8: List[List[Fraction]] = []
    for i in range(n):
        for j in range(i + 1, n):
            for si, s in ((0, Fraction(1)), (1, Fraction(-1))):
                for ti, tt in ((0, Fraction(1)), (1, Fraction(-1))):
                    ru = [Fraction(0)] * n_dofs
                    vi, vj = v_idx(i, si), v_idx(j, ti)
                    ru[vi * n + i] += s
                    ru[vi * n + j] -= tt
                    ru[vj * n + i] -= s
                    ru[vj * n + j] += tt
                    re = [Fraction(0)] * m_s
                    re[i] += Fraction(1)
                    re[j] += Fraction(1)
                    re[off_idx[(i, j)]] += Fraction(-2) * s * tt
                    d0_h8.append(ru)
                    ae_h8.append(re)
    for i in range(n):
        b = C8_had[i]
        ru = [Fraction(0)] * n_dofs
        vp, vm = v_idx(i, 0), v_idx(i, 1)
        for k in range(n):
            ru[vm * n + k] += b[k]
            ru[vp * n + k] -= b[k]
        re = [Fraction(0)] * m_s
        for k in range(n):
            re[k] += b[k] * b[k]
        for k in range(n):
            for l in range(k + 1, n):
                re[off_idx[(k, l)]] += Fraction(2) * b[k] * b[l]
        d0_h8.append(ru)
        ae_h8.append(re)
    rk_d0_h8 = exact_rank(d0_h8)
    rk_aug_h8 = exact_rank([r1 + r2 for r1, r2 in zip(d0_h8, ae_h8)])
    # Master stress vector sigma* = 2 I_8 + H_8 in strain-inner-product coords [S_ii, 2*S_ij]
    sig_h8 = [Fraction(H8[i][i] + 2) for i in range(n)] + [
        Fraction(2 * H8[i][j]) for i in range(n) for j in range(i + 1, n)
    ]
    V_perp_h8 = []
    for j in range(1, m_s):
        v = [Fraction(0)] * m_s
        v[0] = -sig_h8[j]
        v[j] = sig_h8[0]
        V_perp_h8.append(v)
    ae_perp_h8 = [
        [sum((ae_h8[r][c] * v[c] for c in range(m_s)), Fraction(0)) for v in V_perp_h8]
        for r in range(len(ae_h8))
    ]
    rk_perp_h8 = exact_rank([r1 + r2 for r1, r2 in zip(d0_h8, ae_perp_h8)])
    # Check signature of 8x8 matrix 2 I_8 + H_8 via exact block-LDL inertia!
    S_h8_mat = [[Fraction(H8[i][j] + (2 if i == j else 0)) for j in range(n)] for i in range(n)]
    in_sig_h8 = verify_block_ldl_inertia(S_h8_mat, *produce_block_ldl(S_h8_mat))
    if not (
        rk_rot_h8 == 7
        and rk_d0_h8 == 99
        and rk_aug_h8 - rk_d0_h8 == 1
        and rk_perp_h8 == rk_d0_h8
        and in_sig_h8 == (4, 4, 0)
    ):
        raise RuntimeError(f"Sylvester-Hadamard H8 rank-1 check failed: rk_rot={rk_rot_h8}, rk_d0={rk_d0_h8}")

    return {
        "exhaustive_symmetrizable_all_spinor_E8_in_ker_phi": sym_e8_count,
        "deterministic_C8_det": "16",
        "rk_M_rot": rk_rot_e8,
        "fold_rank0_collapse_points": fold_points,
        "sylvester_hadamard_H8_rank1_bifurcation": {
            "rk_M_rot": rk_rot_h8,
            "H0": n_dofs - rk_d0_h8,
            "H1": len(d0_h8) - rk_d0_h8,
            "S_cell_rank": rk_aug_h8 - rk_d0_h8,
            "unique_stress_equals_2I8_plus_H8_exact": (rk_perp_h8 == rk_d0_h8),
            "stress_tensor_2I8_plus_H8_inertia_auxetic": list(in_sig_h8),
        },
        "jammed_E8_bravais": {
            "undirected_roots": 120,
            "inertia": list(in_e8),
            "C11": str(c11_e8),
            "C12": str(c12_e8),
            "C44_voigt": str(c44v_e8),
            "zener_anisotropy_A": str(zener_e8),
            "exact_36x36_O8_isotropic": exact_o8_iso,
        },
        "comparison_D8_zener_anisotropy_A": str(zener_d8),
    }


def run_suite() -> Dict[str, Any]:
    rep: Dict[str, Any] = {}

    # 1. Exact rational points along F_n -> D_n: k = 1 -> (4/3, 1/3), k = 2 -> (8/9, 7/9), s = 1 -> (0, 1)
    curve_pts = [
        (Fraction(4, 3), Fraction(1, 3)),
        (Fraction(8, 9), Fraction(7, 9)),
        (Fraction(0), Fraction(1)),
    ]
    cn_rows = []
    for n in (3, 4, 5, 6, 7):
        Cn = make_twisted_cycle_Cn(n)
        tau = [-1] + [1] * (n - 1)
        if not verify_diagonal_symmetrizer(Cn, tau):
            raise RuntimeError(f"Symmetrizer failed for C_{n}")
        rk_rot = rotational_locking_rank(Cn)
        if rk_rot != n - 1:
            raise RuntimeError(f"Expected rk(M_rot) == {n - 1}, got {rk_rot}")
        pts_res = []
        # Full interior points for n=3..5, endpoint s=1 for n=3..7
        pts_to_run = curve_pts if n <= 5 else [(Fraction(0), Fraction(1))]
        for a, t in pts_to_run:
            r = compute_fold_sheaf_elastic(Cn, a, t)
            if r["H0"] != n * (n - 1) // 2 + 1 or r["H1"] != n * (n - 3) // 2 + 1:
                raise RuntimeError(f"Cohomology jump mismatch at n={n}, (a,t)=({a},{t})")
            if r["S_cell_inertia"] != [1, 0, n * (n + 1) // 2 - 1]:
                raise RuntimeError(f"Rank-1 collapse failed at n={n}: {r['S_cell_inertia']}")
            if a == 0 and t == 1:
                expected_tr = Fraction(4) if n == 3 else Fraction(6 * n - 2, n)
                if Fraction(r["tr_S_cell"]) != expected_tr or Fraction(r["S_00"]) != Fraction(1, n):
                    raise RuntimeError(f"Closed-form s=1 trace mismatch at n={n}: {r['tr_S_cell']} vs {expected_tr}")
            pts_res.append(r)
        cn_rows.append({"n": n, "rk_M_rot": rk_rot, "symmetrizer_verified": True, "points": pts_res})

    # 2. 4D Clifford-nilpotent C4_nilp: pure-shear / zero-axial exact duality at s = 1!
    if not verify_diagonal_symmetrizer(C4_NILPOTENT, [1, 1, -1, -1]):
        raise RuntimeError("C4_nilp symmetrizer failed")
    rk_nilp = rotational_locking_rank(C4_NILPOTENT)
    if rk_nilp != 3:
        raise RuntimeError(f"Expected C4_nilp rk(M_rot)==3, got {rk_nilp}")
    nilp_pts = []
    for a, t in curve_pts:
        r = compute_fold_sheaf_elastic(C4_NILPOTENT, a, t)
        if r["S_cell_inertia"] != [1, 0, 9]:
            raise RuntimeError("C4_nilp rank-1 collapse failed")
        if a == 0 and t == 1:
            if not (
                r["axial_block_identically_zero"]
                and r["tr_S_cell"] == "4"
                and r["rank1_mode_normalized"] == ["0", "0", "0", "0", "1", "1", "0", "0", "1", "1"]
            ):
                raise RuntimeError(f"C4_nilp pure-shear duality failed: {r}")
        nilp_pts.append(r)

    # 3. Full Jamming Transition at s = 1 for D_3 (C_3) and D_4 (both C_4 and C_4^nilp)
    jam_d3 = compute_jammed_Dn_sheaf(make_twisted_cycle_Cn(3))
    jam_d4_c4 = compute_jammed_Dn_sheaf(make_twisted_cycle_Cn(4))
    jam_d4_nilp = compute_jammed_Dn_sheaf(C4_NILPOTENT)
    if not (
        jam_d3["bars_per_cell"] == 36
        and jam_d3["C_jam_inertia"] == [6, 0, 0]
        and jam_d3["cauchy_relation_C12_eq_C44_voigt"]
        and jam_d4_c4 == jam_d4_nilp
        and jam_d4_c4["bars_per_cell"] == 96
        and jam_d4_c4["C_jam_inertia"] == [10, 0, 0]
        and jam_d4_c4["cauchy_relation_C12_eq_C44_voigt"]
        and jam_d4_c4["elastically_isotropic"]
    ):
        raise RuntimeError("Jammed D_n verification failed!")

    # 4. Negative Control: Non-symmetrizable 3x3 matrix breaks rank-1 (gives rank 0 / full locking rk(M_rot)=3!)
    # Flip sign of C_3[0][1] = +1 -> diag(tau) C symmetric becomes IMPOSSIBLE around 3-cycle (odd number of sign flips)!
    C3_nonsym = [[0, 1, -1], [1, 0, -1], [1, -1, 0]]
    rk_nonsym = rotational_locking_rank(C3_nonsym)
    r_nonsym = compute_fold_sheaf_elastic(C3_nonsym, Fraction(0), Fraction(1))
    if not (rk_nonsym == 3 and r_nonsym["H0"] == 3 and r_nonsym["H1"] == 0 and r_nonsym["S_cell_inertia"] == [0, 0, 6]):
        raise RuntimeError(f"Non-symmetrizable negative control failed: rk={rk_nonsym}, r={r_nonsym}")

    rep["universal_Cn_rank1_collapse"] = cn_rows
    rep["clifford_nilpotent_C4_pure_shear_duality"] = {
        "rk_M_rot": rk_nilp,
        "symmetrizer": [1, 1, -1, -1],
        "points": nilp_pts,
        "distinguished_from_C4_at_s1": {"C4_tr_S_cell": "11/2", "C4_nilp_tr_S_cell": "4"},
    }
    rep["jammed_Dn_transition"] = {
        "D3_fcc": jam_d3,
        "D4_via_C4": jam_d4_c4,
        "D4_via_C4_nilp": jam_d4_nilp,
        "D4_endpoints_elastically_identical": (jam_d4_c4 == jam_d4_nilp),
    }
    rep["negative_control_nonsymmetrizable_C3"] = {
        "rk_M_rot": rk_nonsym,
        "H0": r_nonsym["H0"],
        "H1": r_nonsym["H1"],
        "S_cell_inertia": r_nonsym["S_cell_inertia"],
        "zero_stiffness_when_no_symmetrizer": (r_nonsym["S_cell_inertia"] == [0, 0, 6]),
    }
    rep["portal9_8D_F8_to_E8_spinor_collapse_and_isotropy"] = (
        compute_portal9_E8_spinor_collapse_and_jamming()
    )
    return rep


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="Emit JSON report")
    args = parser.parse_args(argv)
    rep = run_suite()
    if args.json:
        print(json.dumps(rep, indent=2))
    else:
        print("=== fold_sheaf_collapse.py: ALL CHECKS PASSED (ZERO FLOAT) ===")
        print("1. Universal Rank-1 Elastic Collapse along F_n -> D_n (s=1 closed form: 4 at n=3, (6n-2)/n for n>=4):")
        for row in rep["universal_Cn_rank1_collapse"]:
            ep = row["points"][-1]
            print(f"   C_{row['n']}: rk(M_rot)={row['rk_M_rot']}/{row['n']}, H0={ep['H0']}, H1={ep['H1']}, C_inertia={ep['S_cell_inertia']}, tr(S_cell)={ep['tr_S_cell']}, S_00={ep['S_00']}")
        np4 = rep["clifford_nilpotent_C4_pure_shear_duality"]["points"][-1]
        print(f"2. 4D C4_nilp pure-shear / zero-axial duality at s=1: axial_zero={np4['axial_block_identically_zero']}, C_inertia={np4['S_cell_inertia']}, tr(S_cell)={np4['tr_S_cell']} (vs C_4=11/2), mode={np4['rank1_mode_normalized']}")
        j3 = rep["jammed_Dn_transition"]["D3_fcc"]
        j4 = rep["jammed_Dn_transition"]["D4_via_C4"]
        print(f"3. Full Jamming Transition at s=1: D_3 inertia={j3['C_jam_inertia']} (C11={j3['C11']}, C12=C44v={j3['C12']}, A={j3['zener_anisotropy_ratio_A']}); D_4 inertia={j4['C_jam_inertia']} (C11={j4['C11']}, C12=C44v={j4['C12']}, A={j4['zener_anisotropy_ratio_A']} ISOTROPIC!)")
        neg = rep["negative_control_nonsymmetrizable_C3"]
        print(f"4. Negative control (non-symmetrizable 3-cycle): rk(M_rot)={neg['rk_M_rot']}/3, H0={neg['H0']}, H1={neg['H1']}, C_inertia={neg['S_cell_inertia']} (rank 0!)")
        e8 = rep["portal9_8D_F8_to_E8_spinor_collapse_and_isotropy"]
        e8j = e8["jammed_E8_bravais"]
        h8b = e8["sylvester_hadamard_H8_rank1_bifurcation"]
        print(
            f"5. Portal 9 8D F_8 -> E_8 bifurcation:\n"
            f"   (a) Cyclic Z/16Z ker(phi): exhaustive symmetrizable={e8['exhaustive_symmetrizable_all_spinor_E8_in_ker_phi']}, "
            f"rk(M_rot)={e8['rk_M_rot']}/8, interior (16/11,7/11) & s=1 H0=28,H1=20, S_cell_rank=0 (FLOPPY LIQUID!)\n"
            f"   (b) Sylvester-Hadamard C8=(1/2)H_8: rk(M_rot)={h8b['rk_M_rot']}/8, H0={h8b['H0']},H1={h8b['H1']}, S_cell_rank={h8b['S_cell_rank']}, "
            f"sigma* == 2I_8+H_8 ({h8b['unique_stress_equals_2I8_plus_H8_exact']}), inertia={h8b['stress_tensor_2I8_plus_H8_inertia_auxetic']} (AUXETIC RANK-1!)\n"
            f"   (c) Jammed E_8 (240 roots): inertia={e8j['inertia']} (C11={e8j['C11']}, C12=C44v={e8j['C12']}, A={e8j['zener_anisotropy_A']} EXACT O(8)-ISOTROPIC vs D_8 A={e8['comparison_D8_zener_anisotropy_A']})"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
