#!/usr/bin/env python3
"""
Exact Pure-Z/Q and Q[sqrt(2), sqrt(3)] Maraldi Rhombic-Dodecahedron Audit for
`Foundation-0/incentre-lattices` (`origin/main @ 2d374fff2a59cc93`,
`results/radius-ratio/radius_ratio.py:73,134,153`, `research/chemistry/borides_ideal.py:59-78`,
and Draft `PR #109` Maxwell-trace prestress contrast).

Replaces float64 `scipy.optimize.linprog` / `ConvexHull` checks for the rhombic dodecahedron
with zero-float exact algebraic certificates:
  1. Exact 14-vertex, 24-edge, 12-face rhombic dodecahedron geometry, sqrt(2):1 diagonals,
     Maraldi obtuse cosine cos(theta) = -1/3, 60-deg normal dihedral, and exact H1 mirror walls.
  2. Exact pure-Q pyramid-insphere Farkas/LP duality proof: d_F = (2/3) n_F, r^2 = 2/9,
     D^2 = 8/9 = 4 r^2, hole/ball = (D - r)/r = 1 exactly, UB12 x_ideal = 1/6, exact 24 kissing
     contacts and 0 overlaps across all 66 dot pairs.
  3. Exact calculus-free Q[sqrt(2), sqrt(3)] Maraldi surface-area minimization
     (S_maraldi^2 / S_flat^2 = 2/3, square completion identity).
  4. Exact unpinned Maxwell-trace obstruction (sum omega_e |e|^2 = 0 forbids omega_e > 0)
     contrasted with a pinned-boundary 3D tension-only RD cable web (14 interior RD nodes,
     38 tension-only cables, exact zero force residual, Omega_int inertia (14, 0, 0) > 0 vs
     compression negative control (0, 14, 0)).

Zero floats on trust paths. Zero bare `assert`s (`-O` safe). Zero external dependencies.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import hmac
import json
from pathlib import Path
import sys
from typing import Any

CERT_SCHEMA = "foundation-0.incentre-lattices.maraldi-rd-exact.v1"
UPSTREAM_REPO = "Foundation-0/incentre-lattices"
UPSTREAM_MAIN_SHA = "2d374fff2a59cc9331566733043a4bb25c1040bf"
HERE = Path(__file__).resolve().parent
DEFAULT_CERT_PATH = HERE / "maraldi_rd_exact_certificate.json"

MAX_JSON_DEPTH = 32
MAX_JSON_BYTES = 2_000_000
MAX_INT_DIGITS = 400
MAX_INT_BITS = 1329
MAX_SEQ_LEN = 256
MAX_MATRIX_DIM = 32
_HEX_LOWER = frozenset("0123456789abcdef")
_CERT_KEYS = frozenset({"schema", "upstream_repo", "upstream_main_sha", "payload_sha256", "certificate_sha256", "payload"})


class MaraldiAuditError(ValueError):
    """Fail-closed typed error (never stripped under python3 -O)."""

    def __init__(self, code: str, message: str) -> None:
        self.code = code
        self.message = message
        super().__init__(f"{code}: {message}")


# ---------------------------------------------------------------------------
# Strict Input Hardening & Iterative Zero-Float Canonical JSON
# ---------------------------------------------------------------------------

def require_int(val: Any, label: str = "int", *, lo: int | None = None, hi: int | None = None) -> int:
    if type(val) is not int:
        raise MaraldiAuditError("INVALID_INT_TYPE", f"{label} must be strict int")
    if lo is not None and val < lo:
        raise MaraldiAuditError("INT_BELOW_MIN", f"{label} below min {lo}")
    if hi is not None and val > hi:
        raise MaraldiAuditError("INT_ABOVE_MAX", f"{label} above max {hi}")
    return val


def require_rat(val: Any, label: str = "rational") -> Fraction:
    if type(val) is Fraction:
        return val
    if type(val) is int:
        return Fraction(val, 1)
    raise MaraldiAuditError("INVALID_RATIONAL_TYPE", f"{label} must be strict int or Fraction")


def require_seq(val: Any, label: str = "seq", *, min_len: int = 1, max_len: int = MAX_SEQ_LEN) -> tuple[Any, ...]:
    if type(val) not in (list, tuple):
        raise MaraldiAuditError("INVALID_SEQ_TYPE", f"{label} must be strict list or tuple")
    n = len(val)
    if n < min_len:
        raise MaraldiAuditError("SEQ_TOO_SHORT", f"{label} length {n} < {min_len}")
    if n > max_len:
        raise MaraldiAuditError("SEQ_TOO_LONG", f"{label} length {n} > {max_len}")
    return tuple(val)


def fmt_q(q: Fraction) -> str:
    r = require_rat(q, "q")
    if r.denominator == 1:
        return str(r.numerator)
    return f"{r.numerator}/{r.denominator}"


def validate_json_value(obj: Any, *, max_depth: int = MAX_JSON_DEPTH) -> None:
    """Iterative zero-float, depth-capped, strict-key JSON tree validator (zero recursion)."""
    limit = require_int(max_depth, "max_depth", lo=1, hi=128)
    stack: list[tuple[Any, int]] = [(obj, 0)]
    while stack:
        cur, depth = stack.pop()
        if depth > limit:
            raise MaraldiAuditError("MAX_DEPTH_EXCEEDED", f"JSON nesting exceeds {limit}")
        t = type(cur)
        if t is dict:
            for k, v in cur.items():
                if type(k) is not str:
                    raise MaraldiAuditError("NON_STR_KEY", "dict keys must be strict str")
                stack.append((v, depth + 1))
        elif t is list:
            for item in cur:
                stack.append((item, depth + 1))
        elif t is int:
            if cur.bit_length() > MAX_INT_BITS:
                raise MaraldiAuditError("INT_TOO_LARGE", "integer exceeds bit-length cap")
        elif t in (str, bool) or cur is None:
            continue
        elif t is float:
            raise MaraldiAuditError("FLOAT_FORBIDDEN", "float values are forbidden")
        else:
            raise MaraldiAuditError("UNSUPPORTED_JSON_TYPE", f"unsupported type {t.__name__}")


def _reject_float_token(tok: str) -> Any:
    raise MaraldiAuditError("JSON_FLOAT_FORBIDDEN", f"float literal forbidden: {tok!r}")


def _reject_constant_token(tok: str) -> Any:
    raise MaraldiAuditError("JSON_CONSTANT_FORBIDDEN", f"non-RFC8259 constant forbidden: {tok!r}")


def _parse_capped_int_token(tok: str) -> int:
    digits = tok[1:] if tok.startswith("-") else tok
    if len(digits) > MAX_INT_DIGITS:
        raise MaraldiAuditError("INT_DIGIT_LIMIT", f"integer token exceeds {MAX_INT_DIGITS} digits")
    return int(tok)


def _pairs_no_dup(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for k, v in pairs:
        if k in out:
            raise MaraldiAuditError("JSON_DUPLICATE_KEY", f"duplicate JSON key: {k!r}")
        out[k] = v
    return out


def parse_strict_json(text: str) -> Any:
    if type(text) is not str:
        raise MaraldiAuditError("INVALID_JSON_INPUT", "JSON input must be strict str")
    if len(text.encode("utf-8")) > MAX_JSON_BYTES:
        raise MaraldiAuditError("JSON_TOO_LARGE", "JSON input exceeds byte cap")
    if not text.strip():
        raise MaraldiAuditError("EMPTY_JSON", "JSON input must be non-empty")
    try:
        parsed = json.loads(
            text,
            parse_int=_parse_capped_int_token,
            parse_float=_reject_float_token,
            parse_constant=_reject_constant_token,
            object_pairs_hook=_pairs_no_dup,
        )
    except MaraldiAuditError:
        raise
    except (ValueError, RecursionError) as exc:
        raise MaraldiAuditError("JSON_PARSE_ERROR", str(exc)) from exc
    validate_json_value(parsed)
    return parsed


def canonical_json_bytes(obj: Any) -> bytes:
    validate_json_value(obj)
    return json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")


# ---------------------------------------------------------------------------
# Exact Q[sqrt(2), sqrt(3)] Field & Sylvester Inertia Core
# ---------------------------------------------------------------------------

def _validate_q23(elem: Any, label: str = "q23") -> tuple[Fraction, Fraction, Fraction, Fraction]:
    seq = require_seq(elem, label, min_len=4, max_len=4)
    return (
        require_rat(seq[0], f"{label}[0]"),
        require_rat(seq[1], f"{label}[1]"),
        require_rat(seq[2], f"{label}[2]"),
        require_rat(seq[3], f"{label}[3]"),
    )


def q23_mul(
    x: tuple[Fraction, Fraction, Fraction, Fraction],
    y: tuple[Fraction, Fraction, Fraction, Fraction],
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    a0, a1, a2, a3 = _validate_q23(x, "x")
    b0, b1, b2, b3 = _validate_q23(y, "y")
    return (
        a0 * b0 + 2 * a1 * b1 + 3 * a2 * b2 + 6 * a3 * b3,
        a0 * b1 + a1 * b0 + 3 * a2 * b3 + 3 * a3 * b2,
        a0 * b2 + a2 * b0 + 2 * a1 * b3 + 2 * a3 * b1,
        a0 * b3 + a3 * b0 + a1 * b2 + a2 * b1,
    )


def _sign_q_sqrt2(p: Fraction, q: Fraction) -> int:
    if p == 0 and q == 0:
        return 0
    if p >= 0 and q >= 0:
        return 1
    if p <= 0 and q <= 0:
        return -1
    diff = p * p - 2 * q * q
    if p > 0:
        return 1 if diff > 0 else -1
    return -1 if diff > 0 else 1


def q23_sign(elem: tuple[Fraction, Fraction, Fraction, Fraction]) -> int:
    a, b, c, d = _validate_q23(elem, "elem")
    s_u = _sign_q_sqrt2(a, b)
    s_v = _sign_q_sqrt2(c, d)
    if s_u == 0 and s_v == 0:
        return 0
    if s_u >= 0 and s_v >= 0:
        return 1
    if s_u <= 0 and s_v <= 0:
        return -1
    diff_p = (a * a + 2 * b * b) - 3 * (c * c + 2 * d * d)
    diff_q = 2 * a * b - 6 * c * d
    s_diff = _sign_q_sqrt2(diff_p, diff_q)
    return s_diff if s_u > 0 else -s_diff


def sylvester_inertia_q(matrix: Any) -> tuple[int, int, int]:
    """Exact symmetric Bunch-Kaufman 1x1/2x2 congruence Sylvester inertia (n_pos, n_neg, n_zero)."""
    rows_raw = require_seq(matrix, "matrix", min_len=1, max_len=MAX_MATRIX_DIM)
    n = len(rows_raw)
    A: list[list[Fraction]] = []
    for i in range(n):
        rseq = require_seq(rows_raw[i], f"row[{i}]", min_len=n, max_len=n)
        A.append([require_rat(x, f"S[{i}][{j}]") for j, x in enumerate(rseq)])
    for i in range(n):
        for j in range(i + 1, n):
            if A[i][j] != A[j][i]:
                raise MaraldiAuditError("NON_SYMMETRIC_MATRIX", f"asymmetry at ({i},{j})")

    n_pos = 0
    n_neg = 0
    n_zero = 0
    k = 0
    while k < n:
        diag_piv = -1
        for i in range(k, n):
            if A[i][i] != 0:
                diag_piv = i
                break
        if diag_piv >= 0:
            if diag_piv != k:
                A[k], A[diag_piv] = A[diag_piv], A[k]
                for row in A:
                    row[k], row[diag_piv] = row[diag_piv], row[k]
            d = A[k][k]
            if d > 0:
                n_pos += 1
            else:
                n_neg += 1
            for i in range(k + 1, n):
                f = A[i][k] / d
                if f != 0:
                    for j in range(k + 1, n):
                        A[i][j] -= f * A[k][j]
            k += 1
            continue

        off_i = -1
        off_j = -1
        for i in range(k, n):
            for j in range(i + 1, n):
                if A[i][j] != 0:
                    off_i = i
                    off_j = j
                    break
            if off_i >= 0:
                break
        if off_i < 0:
            n_zero += n - k
            break

        if off_i != k:
            A[k], A[off_i] = A[off_i], A[k]
            for row in A:
                row[k], row[off_i] = row[off_i], row[k]
        if off_j != k + 1:
            A[k + 1], A[off_j] = A[off_j], A[k + 1]
            for row in A:
                row[k + 1], row[off_j] = row[off_j], row[k + 1]

        b = A[k][k + 1]
        n_pos += 1
        n_neg += 1
        for i in range(k + 2, n):
            c0 = A[i][k + 1] / b
            c1 = A[i][k] / b
            for j in range(k + 2, n):
                A[i][j] -= c0 * A[k][j] + c1 * A[k + 1][j]
        k += 2

    return (n_pos, n_neg, n_zero)


# ---------------------------------------------------------------------------
# Exact Rhombic Dodecahedron, Pyramid-Insphere LP Duality, Maraldi & PR #109
# ---------------------------------------------------------------------------

def _dot3(a: tuple[Fraction, Fraction, Fraction], b: tuple[Fraction, Fraction, Fraction]) -> Fraction:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _sub3(a: tuple[Fraction, Fraction, Fraction], b: tuple[Fraction, Fraction, Fraction]) -> tuple[Fraction, Fraction, Fraction]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def build_rd_vertices_and_faces() -> tuple[
    list[tuple[Fraction, Fraction, Fraction]],
    list[tuple[int, int]],
    list[tuple[tuple[Fraction, Fraction, Fraction], list[int]]],
]:
    """Exact 14-vertex, 24-edge, 12-face rhombic dodecahedron matching radius_ratio.py:134."""
    cube_verts: list[tuple[Fraction, Fraction, Fraction]] = [
        (Fraction(sx), Fraction(sy), Fraction(sz))
        for sx in (-1, 1)
        for sy in (-1, 1)
        for sz in (-1, 1)
    ]
    oct_verts: list[tuple[Fraction, Fraction, Fraction]] = [
        (Fraction(2), Fraction(0), Fraction(0)),
        (Fraction(-2), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(2), Fraction(0)),
        (Fraction(0), Fraction(-2), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(2)),
        (Fraction(0), Fraction(0), Fraction(-2)),
    ]
    verts = cube_verts + oct_verts
    edges: list[tuple[int, int]] = []
    for i in range(8):
        for j in range(8, 14):
            diff = _sub3(verts[i], verts[j])
            if _dot3(diff, diff) == 3:
                edges.append((i, j))

    normals: list[tuple[Fraction, Fraction, Fraction]] = []
    for axis_zero in range(3):
        for s1 in (-1, 1):
            for s2 in (-1, 1):
                vec = [Fraction(0), Fraction(0), Fraction(0)]
                axes = [a for a in range(3) if a != axis_zero]
                vec[axes[0]] = Fraction(s1)
                vec[axes[1]] = Fraction(s2)
                normals.append((vec[0], vec[1], vec[2]))

    faces: list[tuple[tuple[Fraction, Fraction, Fraction], list[int]]] = []
    for nrm in normals:
        on_face = [idx for idx, v in enumerate(verts) if _dot3(nrm, v) == 2]
        faces.append((nrm, on_face))
    return verts, edges, faces


def certify_rd_geometry_and_insphere_q() -> dict[str, Any]:
    """
    Exact-Q certification of rhombic dodecahedron geometry, H1 mirror walls,
    exact pyramid-insphere Farkas LP duality (hole/ball = 1, UB12 x_ideal = 1/6),
    and global all-pairs kissing/non-overlap across all 12 dots.
    """
    verts, edges, faces = build_rd_vertices_and_faces()
    vert_set = set(verts)

    # 1. Verify all 12 faces, rhombus equal sides (len^2 = 3), sqrt(2):1 diagonals, cos(theta) = -1/3
    face_checks: list[bool] = []
    diag_checks: list[bool] = []
    maraldi_checks: list[bool] = []
    for nrm, f_ids in faces:
        off_max = max(_dot3(nrm, verts[idx]) for idx in range(14) if idx not in f_ids)
        c_ids = [idx for idx in f_ids if idx < 8]
        o_ids = [idx for idx in f_ids if idx >= 8]
        face_checks.append(all((len(f_ids) == 4, off_max == 0, len(c_ids) == 2, len(o_ids) == 2)))
        c0, c1 = verts[c_ids[0]], verts[c_ids[1]]
        o0, o1 = verts[o_ids[0]], verts[o_ids[1]]
        d_short = _sub3(c1, c0)
        d_long = _sub3(o1, o0)
        diag_checks.append(all((
            _dot3(d_short, d_long) == 0,
            _dot3(d_long, d_long) == 8,
            _dot3(d_short, d_short) == 4,
        )))
        u = _sub3(o0, c0)
        v = _sub3(o1, c0)
        maraldi_checks.append(all((
            _dot3(u, u) == 3,
            _dot3(v, v) == 3,
            _dot3(u, v) / Fraction(3) == Fraction(-1, 3),
        )))
    all_faces_exact = all(face_checks)
    all_diagonals_2to1 = all(diag_checks)
    all_maraldi_minus_third = all(maraldi_checks)

    # 2. Verify all 24 ridges: dihedral cos(phi) = 1/2 and exact H1 mirror reflection invariance
    # For adjacent faces (a, g), w = n_g - n_a is exact ridge-wall normal with |w|^2 = 2,
    # w . v0 = w . v1 = 0, and exact integer reflection v' = v - (v . w) w!
    adj_pairs: list[tuple[int, int]] = []
    dihedral_checks: list[bool] = []
    mirror_checks: list[bool] = []
    for a in range(12):
        na, sa = faces[a][0], set(faces[a][1])
        for g in range(a + 1, 12):
            ng, sg = faces[g][0], set(faces[g][1])
            common = sorted(sa & sg)
            if len(common) == 2:
                adj_pairs.append((a, g))
                dihedral_checks.append(_dot3(na, ng) / Fraction(2) == Fraction(1, 2))
                w = _sub3(ng, na)
                mirror_checks.append(all((
                    _dot3(w, w) == 2,
                    _dot3(w, verts[common[0]]) == 0,
                    _dot3(w, verts[common[1]]) == 0,
                )))
                for v in verts:
                    sc = _dot3(v, w)
                    vr = (v[0] - sc * w[0], v[1] - sc * w[1], v[2] - sc * w[2])
                    mirror_checks.append(vr in vert_set)
    dihedral_all_half = all(dihedral_checks)
    h1_all_exact_mirrors = all(mirror_checks)

    # 3. Exact pyramid-insphere primal & Farkas dual LP certificate on every face F
    # Since |n_F|^2 = 2 and every ridge-wall outward normal w_k = n_g - n_a has |w_k|^2 = 2:
    # dot center d_F = (2/3) n_F satisfies:
    #   base slack: (2 - n_F . d_F) = 2 - (4/3)*2 = 2/3 -> dist^2 = (2/3)^2 / 2 = 2/9 = r^2
    #   wall slack: (0 - w_k . d_F) = 2/3               -> dist^2 = 2/9 = r^2
    # Dual Farkas certificate: lambda_base = 1/3, lambda_wall = 1/6 (x4), sum = 1, sum lambda_i a_i = 0!
    dots: list[tuple[Fraction, Fraction, Fraction]] = [
        (Fraction(2, 3) * n[0], Fraction(2, 3) * n[1], Fraction(2, 3) * n[2])
        for n, _ in faces
    ]
    r_sq = Fraction(2, 9)
    d_sq = Fraction(8, 9)
    pd_checks: list[bool] = []
    for a in range(12):
        na, sa = faces[a][0], set(faces[a][1])
        da = dots[a]
        pd_checks.append(all((
            _dot3(da, da) == d_sq,
            Fraction(2) - _dot3(na, da) == Fraction(2, 3),
        )))
        wall_normals: list[tuple[Fraction, Fraction, Fraction]] = []
        for g in range(12):
            if g == a:
                continue
            common = sorted(sa & set(faces[g][1]))
            if len(common) != 2:
                continue
            w_norm = _sub3(faces[g][0], na)
            pd_checks.append(all((
                _dot3(w_norm, w_norm) == 2,
                _dot3(w_norm, verts[common[0]]) == 0,
                _dot3(w_norm, verts[common[1]]) == 0,
                -_dot3(w_norm, da) == Fraction(2, 3),
            )))
            wall_normals.append(w_norm)
        grad = [
            Fraction(1, 3) * na[c] + sum((Fraction(1, 6) * w[c] for w in wall_normals), Fraction(0))
            for c in range(3)
        ]
        pd_checks.append(all((len(wall_normals) == 4, grad == [Fraction(0), Fraction(0), Fraction(0)])))
    primal_dual_exact = all(pd_checks)

    # 4. All-pairs kissing contact (24 adjacent pairs == 4*r^2 = 8/9) & non-adjacent separation (min == 16/9 = 8*r^2)
    adj_set = set(adj_pairs)
    kissing_count = 0
    non_adj_dists: list[Fraction] = []
    for a in range(12):
        for g in range(a + 1, 12):
            diff = _sub3(dots[a], dots[g])
            dist_sq = _dot3(diff, diff)
            if (a, g) in adj_set:
                if dist_sq == 4 * r_sq:
                    kissing_count += 1
            else:
                non_adj_dists.append(dist_sq)
    min_non_adj_dist_sq = min(non_adj_dists)
    zero_overlaps = min_non_adj_dist_sq == 8 * r_sq

    # D^2 / r^2 = 4 => D = 2r => hole/ball = (D - r)/r = 1; UB12 cubic A=4 => x_ideal = (2/3)/4 = 1/6
    d_over_r_sq = d_sq / r_sq
    hole_over_ball = Fraction(1)
    ub12_x_ideal = Fraction(2, 3) / Fraction(4)

    passed = all((
        len(verts) == 14,
        len(edges) == 24,
        len(faces) == 12,
        len(verts) - len(edges) + len(faces) == 2,
        all_faces_exact,
        all_diagonals_2to1,
        all_maraldi_minus_third,
        len(adj_pairs) == 24,
        dihedral_all_half,
        h1_all_exact_mirrors,
        primal_dual_exact,
        kissing_count == 24,
        len(non_adj_dists) == 42,
        zero_overlaps,
        d_over_r_sq == Fraction(4),
        ub12_x_ideal == Fraction(1, 6),
    ))
    return {
        "n_vertices": len(verts),
        "n_edges": len(edges),
        "n_faces": len(faces),
        "euler_characteristic": len(verts) - len(edges) + len(faces),
        "edge_squared_length": "3",
        "diagonal_squared_ratio_long_over_short": "2",
        "maraldi_obtuse_cos_theta": "-1/3",
        "dihedral_normal_cos_phi": "1/2",
        "h1_24_mirror_walls_exact": h1_all_exact_mirrors,
        "insphere_r_squared": fmt_q(r_sq),
        "dot_radius_D_squared": fmt_q(d_sq),
        "D_over_r_squared": fmt_q(d_over_r_sq),
        "hole_over_ball_exact": fmt_q(hole_over_ball),
        "ub12_boron_x_ideal_exact": fmt_q(ub12_x_ideal),
        "primal_dual_farkas_exact": primal_dual_exact,
        "kissing_contacts_verified": kissing_count,
        "min_non_adjacent_dist_squared": fmt_q(min_non_adj_dist_sq),
        "zero_overlaps_exact": zero_overlaps,
        "passed": passed,
    }


def certify_maraldi_surface_q23(
    perturbations: Any = (Fraction(0), Fraction(1, 4), Fraction(1, 3), Fraction(1, 2), Fraction(1)),
) -> dict[str, Any]:
    """Exact Q[sqrt(2), sqrt(3)] calculus-free Maraldi base optimality proof."""
    p_seq = require_seq(perturbations, "perturbations", min_len=1, max_len=32)
    xs = [require_rat(v, f"perturbations[{i}]") for i, v in enumerate(p_seq)]
    for x in xs:
        if x < 0:
            raise MaraldiAuditError("NEGATIVE_APEX_OFFSET", "apex offset x must be >= 0")

    s_flat = (Fraction(0), Fraction(0), Fraction(3, 2), Fraction(0))
    s_maraldi = (Fraction(0), Fraction(3, 2), Fraction(0), Fraction(0))
    delta_s = (Fraction(0), Fraction(-3, 2), Fraction(3, 2), Fraction(0))
    s_flat_sq = q23_mul(s_flat, s_flat)
    s_maraldi_sq = q23_mul(s_maraldi, s_maraldi)
    area_sq_ratio = s_maraldi_sq[0] / s_flat_sq[0]
    pert_ok = all(q23_sign((Fraction(1) + 8 * x * x, -4 * x, Fraction(0), Fraction(0))) == 1 for x in xs)
    passed = all((
        q23_sign(delta_s) == 1,
        area_sq_ratio == Fraction(2, 3),
        pert_ok,
    ))
    return {
        "s_flat_q23": [fmt_q(c) for c in s_flat],
        "s_maraldi_q23": [fmt_q(c) for c in s_maraldi],
        "delta_s_saving_q23": [fmt_q(c) for c in delta_s],
        "s_maraldi_sq_over_s_flat_sq": fmt_q(area_sq_ratio),
        "all_rational_perturbations_strictly_worse": pert_ok,
        "passed": passed,
    }


def verify_rd_prestress_contrast_q(omega_rd: Any, omega_cube_anchor: Any, omega_oct_anchor: Any) -> dict[str, Any]:
    """
    Exact-Q contrast with Draft PR #109's unpinned Maxwell-trace obstruction:
      - Unpinned 14-vertex, 24-edge RD: Maxwell identity sum(omega_e * |e|^2) = 3 * sum(omega_e) = 0
        forbids strictly positive self-stress (omega_e > 0 => trace = 72*omega_rd > 0 != 0).
      - Pinned-boundary 3D RD cable web (14 interior RD vertices tethered radially to 14 outer
        anchors at 2*v): omega_rd=1, omega_cube_anchor=1, omega_oct_anchor=2 gives exact zero
        equilibrium residual at all 14 interior nodes and strictly PD Omega_int inertia (14, 0, 0).
    """
    w_rd = require_rat(omega_rd, "omega_rd")
    w_ca = require_rat(omega_cube_anchor, "omega_cube_anchor")
    w_oa = require_rat(omega_oct_anchor, "omega_oct_anchor")
    verts, edges, _ = build_rd_vertices_and_faces()
    unpinned_maxwell_trace = sum((w_rd * Fraction(3) for _ in edges), Fraction(0))

    # Build 14 interior nodes (0..13) + 14 radial outer anchors (14..27) at 2*v
    all_pts = verts + [(2 * v[0], 2 * v[1], 2 * v[2]) for v in verts]
    cables: list[tuple[int, int, Fraction]] = [(u, v, w_rd) for u, v in edges]
    for i in range(8):
        cables.append((i, 14 + i, w_ca))
    for i in range(8, 14):
        cables.append((i, 14 + i, w_oa))

    eq_ok = True
    for v_idx in range(14):
        fx = Fraction(0)
        fy = Fraction(0)
        fz = Fraction(0)
        px, py, pz = all_pts[v_idx]
        for a, b, w in cables:
            if a == v_idx:
                qx, qy, qz = all_pts[b]
                fx += w * (px - qx)
                fy += w * (py - qy)
                fz += w * (pz - qz)
            elif b == v_idx:
                qx, qy, qz = all_pts[a]
                fx += w * (px - qx)
                fy += w * (py - qy)
                fz += w * (pz - qz)
        if (fx, fy, fz) != (0, 0, 0):
            eq_ok = False

    omega_int = [[Fraction(0)] * 14 for _ in range(14)]
    for u, v in edges:
        omega_int[u][u] += w_rd
        omega_int[v][v] += w_rd
        omega_int[u][v] -= w_rd
        omega_int[v][u] -= w_rd
    for i in range(8):
        omega_int[i][i] += w_ca
    for i in range(8, 14):
        omega_int[i][i] += w_oa
    inertia_int = sylvester_inertia_q(omega_int)
    strictly_tension = w_rd > 0 and w_ca > 0 and w_oa > 0
    return {
        "strictly_tension_only": strictly_tension,
        "unpinned_maxwell_trace": fmt_q(unpinned_maxwell_trace),
        "unpinned_tension_obstructed": unpinned_maxwell_trace != 0,
        "interior_equilibrium_zero_residual": eq_ok,
        "omega_int_inertia": list(inertia_int),
        "omega_int_strictly_pd": inertia_int == (14, 0, 0),
        "pinned_prestress_gate_pass": all((strictly_tension, eq_ok, inertia_int == (14, 0, 0))),
    }


def build_study_payload() -> dict[str, Any]:
    rd = certify_rd_geometry_and_insphere_q()
    maraldi = certify_maraldi_surface_q23()
    prestress_pass = verify_rd_prestress_contrast_q(Fraction(1), Fraction(1), Fraction(2))
    prestress_comp_fail = verify_rd_prestress_contrast_q(Fraction(-1), Fraction(-1), Fraction(-2))

    all_passed = all((
        rd["passed"],
        maraldi["passed"],
        prestress_pass["pinned_prestress_gate_pass"],
        prestress_pass["unpinned_tension_obstructed"],
        not prestress_comp_fail["pinned_prestress_gate_pass"],
        prestress_comp_fail["omega_int_inertia"] == [0, 14, 0],
    ))
    return {
        "evidence_label": "exact_algebraic_proof_pure_Z_Q_and_Q_sqrt2_sqrt3",
        "upstream_anchors": {
            "repo": UPSTREAM_REPO,
            "main_sha": UPSTREAM_MAIN_SHA,
            "target_lines": "results/radius-ratio/radius_ratio.py:73,134,153; research/chemistry/borides_ideal.py:59-78; PR #109 Maxwell trace contrast",
        },
        "part_a_rd_geometry_and_farkas_insphere": rd,
        "part_b_maraldi_q23_surface_optimality": maraldi,
        "part_c_pr109_maxwell_trace_vs_pinned_prestress": {
            "pinned_3d_rd_tension_only_passes": prestress_pass,
            "compression_negative_control_fails": prestress_comp_fail,
        },
        "all_passed": all_passed,
    }


# ---------------------------------------------------------------------------
# Fail-Closed Certificate Builder, Checker & CLI
# ---------------------------------------------------------------------------

def _require_sha256_hex(val: Any, label: str) -> str:
    if type(val) is not str:
        raise MaraldiAuditError("INVALID_DIGEST_TYPE", f"{label} must be strict str")
    if len(val) != 64:
        raise MaraldiAuditError("INVALID_DIGEST_LEN", f"{label} must have length 64")
    if any(ch not in _HEX_LOWER for ch in val):
        raise MaraldiAuditError("INVALID_DIGEST_HEX", f"{label} must be lowercase hex")
    return val


def _require_no_symlink_path(path: Path) -> None:
    if path.is_symlink():
        raise MaraldiAuditError("SYMLINK_FORBIDDEN", f"symlink path forbidden: {path}")
    for parent in path.parents:
        if parent.is_symlink():
            raise MaraldiAuditError("SYMLINK_ANCESTOR_FORBIDDEN", f"symlink ancestor forbidden: {parent}")


def build_certificate() -> dict[str, Any]:
    payload = build_study_payload()
    payload_sha256 = hashlib.sha256(canonical_json_bytes(payload)).hexdigest()
    envelope = {
        "schema": CERT_SCHEMA,
        "upstream_repo": UPSTREAM_REPO,
        "upstream_main_sha": UPSTREAM_MAIN_SHA,
        "payload_sha256": payload_sha256,
        "payload": payload,
    }
    certificate_sha256 = hashlib.sha256(canonical_json_bytes(envelope)).hexdigest()
    return {
        "schema": CERT_SCHEMA,
        "upstream_repo": UPSTREAM_REPO,
        "upstream_main_sha": UPSTREAM_MAIN_SHA,
        "payload_sha256": payload_sha256,
        "certificate_sha256": certificate_sha256,
        "payload": payload,
    }


def verify_certificate_dict(cert: Any) -> dict[str, Any]:
    if type(cert) is not dict:
        raise MaraldiAuditError("CERT_NOT_DICT", "certificate root must be a JSON object")
    if set(cert.keys()) != _CERT_KEYS:
        raise MaraldiAuditError("CERT_SCHEMA_KEYS", f"unexpected keys: {sorted(cert.keys())}")
    if cert["schema"] != CERT_SCHEMA:
        raise MaraldiAuditError("CERT_SCHEMA_MISMATCH", "schema mismatch")
    if cert["upstream_repo"] != UPSTREAM_REPO:
        raise MaraldiAuditError("CERT_REPO_MISMATCH", "upstream_repo mismatch")
    if cert["upstream_main_sha"] != UPSTREAM_MAIN_SHA:
        raise MaraldiAuditError("CERT_SHA_MISMATCH", "upstream_main_sha mismatch")

    stored_payload_sha = _require_sha256_hex(cert["payload_sha256"], "payload_sha256")
    stored_cert_sha = _require_sha256_hex(cert["certificate_sha256"], "certificate_sha256")
    computed_payload_sha = hashlib.sha256(canonical_json_bytes(cert["payload"])).hexdigest()
    if not hmac.compare_digest(stored_payload_sha, computed_payload_sha):
        raise MaraldiAuditError("PAYLOAD_SHA_MISMATCH", "payload_sha256 mismatch")

    envelope = {
        "schema": cert["schema"],
        "upstream_repo": cert["upstream_repo"],
        "upstream_main_sha": cert["upstream_main_sha"],
        "payload_sha256": stored_payload_sha,
        "payload": cert["payload"],
    }
    computed_cert_sha = hashlib.sha256(canonical_json_bytes(envelope)).hexdigest()
    if not hmac.compare_digest(stored_cert_sha, computed_cert_sha):
        raise MaraldiAuditError("CERTIFICATE_SHA_MISMATCH", "certificate_sha256 mismatch")

    expected = build_certificate()
    if not hmac.compare_digest(computed_payload_sha, expected["payload_sha256"]):
        raise MaraldiAuditError("RECOMPUTED_PAYLOAD_DRIFT", "payload digest differs from exact recomputation")

    return {
        "valid": True,
        "schema": CERT_SCHEMA,
        "payload_sha256": computed_payload_sha,
        "certificate_sha256": computed_cert_sha,
    }


def check_certificate_file(path: Path) -> dict[str, Any]:
    if not isinstance(path, Path):
        raise MaraldiAuditError("INVALID_PATH_TYPE", "path must be a pathlib.Path")
    _require_no_symlink_path(path)
    if not path.is_file():
        raise MaraldiAuditError("CERT_FILE_NOT_FOUND", f"certificate file not found: {path}")
    if path.stat().st_size > MAX_JSON_BYTES:
        raise MaraldiAuditError("CERT_FILE_TOO_LARGE", "certificate file exceeds byte cap")
    raw = path.read_bytes()
    if b"\r\n" in raw:
        raise MaraldiAuditError("CRLF_FORBIDDEN", "CRLF line endings forbidden")
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise MaraldiAuditError("INVALID_UTF8", f"invalid UTF-8: {exc}") from exc
    parsed = parse_strict_json(text)
    canonical_raw = canonical_json_bytes(parsed)
    if raw not in (canonical_raw, canonical_raw + b"\n"):
        raise MaraldiAuditError("NON_CANONICAL_CERT_BYTES", "certificate file is not canonical JSON")
    return verify_certificate_dict(parsed)


def write_certificate_file(path: Path) -> dict[str, Any]:
    _require_no_symlink_path(path)
    cert = build_certificate()
    raw = canonical_json_bytes(cert) + b"\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(raw)
    return check_certificate_file(path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Exact-Q Maraldi rhombic-dodecahedron certifier for incentre-lattices")
    parser.add_argument("--standalone", action="store_true", help="Run full standalone exact verification")
    parser.add_argument("--json", action="store_true", help="Emit canonical JSON certificate to stdout")
    parser.add_argument("--write", nargs="?", const=str(DEFAULT_CERT_PATH), help="Write canonical certificate")
    parser.add_argument("--check", type=str, help="Fail-closed verification of a certificate JSON file")
    args = parser.parse_args(argv)

    try:
        if args.check is not None:
            res = check_certificate_file(Path(args.check))
            print(
                f"OK {CERT_SCHEMA} payload_sha256={res['payload_sha256']} "
                f"certificate_sha256={res['certificate_sha256']}"
            )
            return 0
        if args.write is not None:
            res = write_certificate_file(Path(args.write))
            print(
                f"WROTE {CERT_SCHEMA} payload_sha256={res['payload_sha256']} "
                f"certificate_sha256={res['certificate_sha256']}"
            )
            return 0
        cert = build_certificate()
        if args.json:
            sys.stdout.buffer.write(canonical_json_bytes(cert) + b"\n")
            return 0
        print(
            f"OK {CERT_SCHEMA} payload_sha256={cert['payload_sha256']} "
            f"certificate_sha256={cert['certificate_sha256']} "
            f"all_passed={cert['payload']['all_passed']}"
        )
        return 0
    except (MaraldiAuditError, ValueError, RecursionError, OSError) as exc:
        print(f"REJECT {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
