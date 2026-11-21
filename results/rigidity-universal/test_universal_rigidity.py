#!/usr/bin/env python3
"""Tests for universal_rigidity_check.py. Standard library only.

Run: python3 -I -B results/rigidity-universal/test_universal_rigidity.py
The known-answer numbers below are written as literals, not computed by the module under test.
"""
import itertools
import pathlib
import random
import sys
import unittest
from fractions import Fraction as F

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import universal_rigidity_check as u  # noqa: E402


class KnownAnswers(unittest.TestCase):
    def test_cross_polytope_stress_dimensions(self):
        for n, want in {2: 0, 3: 0, 4: 2, 5: 5, 6: 9, 7: 14}.items():
            pts, bars = u.cross_polytope(n)
            self.assertEqual(len(u.stress_space(pts, bars)), want, n)

    def test_periodic_zero_modes_and_kernel(self):
        for n in (2, 3, 4):
            c = F(2, 5)
            self.assertEqual(u.check_periodic(random.Random(1), n, c), {2: 1, 3: 3, 4: 6}[n])
            pts, bars = u.incentre(n, c)
            Om = u.omega0(pts, bars, u.periodic_stress(n, c))
            self.assertEqual(u.psd_rank(Om), (True, 2 * n - 1))

    def test_bridging_tension_literal(self):
        # n = 3, c = 1/2: c (n-1)/(1-c) = 2
        om = u.periodic_stress(3, F(1, 2))
        self.assertEqual(om[-3:], [2, 2, 2])
        self.assertEqual(om[:-3], [1] * 12)

    def test_square_with_diagonals(self):
        sq = [(F(0), F(0)), (F(1), F(0)), (F(1), F(1)), (F(0), F(1))]
        bars = [(0, 1), (1, 2), (2, 3), (0, 3), (0, 2), (1, 3)]
        om = [F(1)] * 4 + [F(-1)] * 2
        self.assertEqual(u.maxwell_defects(sq, bars, om), [])
        self.assertEqual(u.psd_rank(u.stress_matrix(4, bars, om)), (True, 1))

    def test_simplex_map(self):
        pts, bars = u.cross_polytope(3)
        self.assertEqual(u.simplex_map_defects(pts, bars, u.simplex_image(6), 2), [])


class MustFail(unittest.TestCase):
    def test_wrong_bridging_tension(self):
        for shift in (F(1), F(-1, 2)):
            om = u.periodic_stress(3, F(1, 2), tau_shift=shift)
            self.assertTrue(u.periodic_defects(random.Random(1), 3, F(1, 2), om))

    def test_sign_flip(self):
        om = u.periodic_stress(3, F(1, 3))
        om[5] = -om[5]
        self.assertTrue(u.periodic_defects(random.Random(1), 3, F(1, 3), om))

    def test_negated_stress_is_equilibrium_but_not_psd(self):
        om = [-x for x in u.periodic_stress(3, F(1, 3))]
        pts, bars = u.incentre(3, F(1, 3))
        self.assertEqual(u.periodic_equilibrium_defects(pts, bars, om), [])
        self.assertEqual(u.kernel_defects(u.omega0(pts, bars, om), 1), ["not positive semidefinite"])

    def test_psd_test_rejects_indefinite(self):
        Om = [[F(0), F(1)], [F(1), F(0)]]
        self.assertEqual(u.psd_rank(Om)[0], False)
        pts, bars = u.cross_polytope(4)
        w = u.zero_row_sum_weights(4)[0]
        Om = u.stress_matrix(8, bars, u.axis_pair_weights(4, w))
        self.assertFalse(u.psd_rank(Om)[0])

    def test_wrong_family_dimension_is_noticed(self):
        # n = 4: a framework with one bar removed has fewer stresses than the cross-polytope
        pts, bars = u.cross_polytope(4)
        self.assertEqual(len(u.stress_space(pts, bars[1:])), 1)

    def test_corrupt_simplex_image(self):
        pts, bars = u.cross_polytope(3)
        img = u.simplex_image(6)
        img[2] = tuple(x * 2 for x in img[2])
        self.assertTrue(u.simplex_map_defects(pts, bars, img, 2))

    def test_equal_bars_at_d_plus_one_is_a_congruence(self):
        tet = [tuple(map(F, p)) for p in ((0, 0, 0), (1, 1, 0), (1, 0, 1), (0, 1, 1))]
        bars = list(itertools.combinations(range(4), 2))
        self.assertEqual(u.simplex_map_defects(tet, bars, u.simplex_image(4), 2), ["image is a congruence"])

    def test_corrupt_flag_exits_one(self):
        with self.assertRaises(u.CheckFailed):
            u.check_periodic(random.Random(1), 3, F(1, 2), tau_shift=F(1))


class FullRun(unittest.TestCase):
    def test_controls_and_main(self):
        rng = random.Random(u.SEED)
        u.check_maxwell(rng)
        u.run_controls(rng)


if __name__ == "__main__":
    unittest.main(verbosity=1)
