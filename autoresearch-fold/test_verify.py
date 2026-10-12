#!/usr/bin/env python3
"""Tests for verify.py (the trusted checker). Run: .venv/bin/python test_verify.py"""
import os
import tempfile
import unittest
import numpy as np
import verify
from verify import canonical_state, verify_arrays, verify_file, R, K


def rot(n, th):
    n = np.asarray(n, float) / np.linalg.norm(n)
    Kx = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return np.eye(3) + np.sin(th) * Kx + (1 - np.cos(th)) * Kx @ Kx


def cubic_state(m, th):
    """Independent reference: checkerboard rotation by +-th about [111], cubic lattice."""
    c = np.cos(th)
    q = 2 * K * (1 - c) / 3
    # (2a - 2Kc - q)^2 + 2 q^2 = (2R)^2, take the root continuous with a=1 at th=0
    two_a = 2 * K * c + q + np.sqrt(4 * R ** 2 - 2 * q ** 2)
    a = two_a / 2
    L = 2 * a * m * np.eye(3)
    X = np.zeros((6 * m ** 3, 3))
    for ix in range(m):
        for iy in range(m):
            for iz in range(m):
                cc = (ix * m + iy) * m + iz
                par = (ix + iy + iz) % 2
                Rm = rot([1, 1, 1], th if par == 0 else -th)
                ctr = 2 * a * np.array([ix, iy, iz]) + a
                for ax in range(3):
                    for sg in range(2):
                        e = np.zeros(3)
                        e[ax] = 1 if sg == 0 else -1
                        X[6 * cc + 2 * ax + sg] = ctr + K * Rm @ e
    return L, X


def cubic_path(m, th_end, n=None, th0=0.0):
    ths = []
    th = th0
    prev = cubic_state(m, th)
    ths.append(th)
    Ls, Xs = [prev[0]], [prev[1]]
    d = 0.002
    while th < th_end - 1e-15:
        th = min(th + d, th_end)
        L, X = cubic_state(m, th)
        Ls.append(L)
        Xs.append(X)
    return np.array(Ls), np.array(Xs)


def run(m, Ls, Xs):
    return verify_arrays(m, Ls, Xs)


def has(rep, key):
    return any(r.startswith(key) for r in rep["reasons"])


class TestMustPass(unittest.TestCase):
    def test_trivial_path(self):
        L, X = canonical_state(2)
        rep = run(2, L[None], X[None])
        self.assertTrue(rep["valid"], rep)
        self.assertAlmostEqual(rep["final_density"], np.pi * R ** 3, places=9)
        rep = run(2, np.array([L, L]), np.array([X, X]))
        self.assertTrue(rep["valid"], rep)

    def test_trivial_all_m(self):
        for m in (1, 3, 4):
            L, X = canonical_state(m)
            self.assertTrue(run(m, L[None], X[None])["valid"])

    def test_cubic_111_to_60(self):
        Ls, Xs = cubic_path(2, np.pi / 3)
        rep = run(2, Ls, Xs)
        self.assertTrue(rep["valid"], rep["reasons"])
        self.assertAlmostEqual(rep["final_density"], 27 * np.pi * np.sqrt(2) / 256, delta=1e-4)
        self.assertAlmostEqual(rep["final_density"], 0.4686, delta=1e-4)
        self.assertEqual(rep["max_density"], rep["final_density"])

    def test_file_roundtrip(self):
        Ls, Xs = cubic_path(2, 0.3)
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "p.npz")
            np.savez(p, m=2, Ls=Ls, Xs=Xs)
            self.assertTrue(verify_file(p)["valid"])
            self.assertFalse(verify_file(os.path.join(d, "missing.npz"))["valid"])


class TestMustFail(unittest.TestCase):
    def setUp(self):
        self.Ls, self.Xs = cubic_path(2, 0.6)

    def test_jump(self):
        keep = list(range(0, 100)) + list(range(150, len(self.Ls)))
        rep = run(2, self.Ls[keep], self.Xs[keep])
        self.assertFalse(rep["valid"])
        self.assertTrue(has(rep, "continuity"), rep["reasons"])

    def test_stretched_bar(self):
        Xs = self.Xs.copy()
        Xs[50, 7] += np.array([2e-6, 0, 0])
        rep = run(2, self.Ls, Xs)
        self.assertFalse(rep["valid"])
        self.assertTrue(has(rep, "bar"), rep["reasons"])

    def test_overlap(self):
        # the cubic [111] family is a legal motion up to 60 deg; beyond it balls interpenetrate
        Ls, Xs = cubic_path(2, 1.25)
        rep = run(2, Ls, Xs)
        self.assertFalse(rep["valid"])
        self.assertTrue(has(rep, "overlap") or has(rep, "midpoint-overlap"), rep["reasons"])

    def test_overlap_with_good_bars_but_wrong_branch(self):
        # lattice squeezed past the 60-degree state by rescaling everything: bars break AND balls overlap
        L, X = canonical_state(2)
        f = 0.99
        rep = run(2, np.array([L, L * f]), np.array([X, X * f]))
        self.assertFalse(rep["valid"])
        self.assertTrue(has(rep, "overlap") or has(rep, "midpoint-overlap"), rep["reasons"])

    def test_wrong_initial(self):
        Ls, Xs = cubic_path(2, 0.5, th0=0.1)
        rep = run(2, Ls, Xs)
        self.assertFalse(rep["valid"])
        self.assertTrue(has(rep, "initial"), rep["reasons"])

    def test_translated_initial(self):
        L, X = canonical_state(2)
        rep = run(2, L[None], (X + 1e-6)[None])
        self.assertFalse(rep["valid"])
        self.assertTrue(has(rep, "initial"))

    def test_altered_topology(self):
        # rotate ONE octahedron (cell 0) about [111]: its own 12 bars stay intact, the
        # inter-cell bars (absent from a fake, smaller bar list we also ship in the npz) break.
        L, X = canonical_state(2)
        Ls, Xs = [L], [X]
        ctr = X[:6].mean(axis=0)
        for th in np.arange(0.002, 0.2, 0.002):
            Xn = X.copy()
            Xn[:6] = ctr + (X[:6] - ctr) @ rot([1, 1, 1], th).T
            Ls.append(L)
            Xs.append(Xn)
        Ls, Xs = np.array(Ls), np.array(Xs)
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "p.npz")
            np.savez(p, m=2, Ls=Ls, Xs=Xs, bars=np.zeros((0, 2), int))   # fake topology ignored
            rep = verify_file(p)
        self.assertFalse(rep["valid"])
        self.assertTrue(has(rep, "bar"), rep["reasons"])

    def test_negative_det(self):
        L, X = canonical_state(2)
        Lf = L[[1, 0, 2]]
        rep = run(2, np.array([L, Lf]), np.array([X, X]))
        self.assertFalse(rep["valid"])
        self.assertTrue(has(rep, "det"), rep["reasons"])

    def test_lattice_shrinks_dots_dont_follow(self):
        Ls = self.Ls.copy()
        Ls[60:] *= 0.995
        Ls[60:] = Ls[60:]
        rep = run(2, Ls, self.Xs)
        self.assertFalse(rep["valid"])
        self.assertTrue(has(rep, "bar") or has(rep, "continuity"), rep["reasons"])

    def test_lattice_shrinks_tiny(self):
        Ls = self.Ls.copy()
        Ls[:, 0, 0] *= (1 - 1e-6)
        rep = run(2, Ls, self.Xs)
        self.assertFalse(rep["valid"])

    def test_nan(self):
        Xs = self.Xs.copy()
        Xs[10, 3, 1] = np.nan
        rep = run(2, self.Ls, Xs)
        self.assertFalse(rep["valid"])
        self.assertTrue(has(rep, "nonfinite"))
        Ls = self.Ls.copy()
        Ls[5, 0, 0] = np.inf
        self.assertFalse(run(2, Ls, self.Xs)["valid"])

    def test_bad_format(self):
        L, X = canonical_state(2)
        self.assertFalse(run(5, L[None], X[None])["valid"])
        self.assertFalse(run(2, L[None], X[None, :10])["valid"])
        self.assertFalse(run(2, L[None][:0], X[None][:0])["valid"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
