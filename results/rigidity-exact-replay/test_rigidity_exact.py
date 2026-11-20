#!/usr/bin/env python3
"""Tests for rigidity_exact.py and modp_check.py: known-answer values and controls that must fail.

Standard library only (unittest). Uses unittest methods, not bare assert statements, so it also runs under python -O.
Run:  python3 -I -B -W error results/rigidity-exact-replay/test_rigidity_exact.py
      python3 -I -O -B -W error results/rigidity-exact-replay/test_rigidity_exact.py

The known answers below are written out as literal numbers (taken from the formulas in results/rigidity-nd/PROOF.md,
evaluated by hand), not computed from the formulas in the scripts, so a wrong formula in a script cannot also be wrong here.
"""
import importlib.util
import json
import pathlib
import subprocess
import sys
import unittest

HERE = pathlib.Path(__file__).resolve().parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / (name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


RE = load("rigidity_exact")
MP = load("modp_check")

# n:      zero modes, self-stresses of the lone cross-polytope (n(n+1)/2 and n(n-3)/2; n = 2: one internal flex)
CP = {2: (4, 0), 3: (6, 0), 4: (10, 2), 5: (15, 5), 6: (21, 9), 7: (28, 14), 8: (36, 20)}
# generic phase: n(n-1)/2 and n(n-3)/2 (n = 2: 2 and 0)
GEN = {2: (2, 0), 3: (3, 0), 4: (6, 2), 5: (10, 5), 6: (15, 9), 7: (21, 14), 8: (28, 20)}
# q = 0: n(n+1)/2 and n(n-1)/2 (n = 2: 3 and 1)
Q0 = {2: (3, 1), 3: (6, 3), 4: (10, 6), 5: (15, 10), 6: (21, 15), 7: (28, 21), 8: (36, 28)}


def run(script, *args):
    return subprocess.run([sys.executable, "-I", "-B", str(HERE / script), *args], capture_output=True, text=True, timeout=600)


class ExactKnownAnswers(unittest.TestCase):
    def test_counts_match_literal_tables(self):
        for n in range(2, 6):
            pos, ib = RE.positions(n), RE.intra_bars(n)
            cp = RE.counts(n, pos, ib)
            self.assertEqual((cp["zero_modes"], cp["self_stresses"]), CP[n], f"cross-polytope n={n}")
            q0 = RE.counts(n, pos, ib, [RE.ONE] * n, True)
            self.assertEqual((q0["zero_modes"], q0["self_stresses"]), Q0[n], f"q=0 n={n}")
            gen = RE.counts(n, pos, ib, RE.generic_z(n), True)
            self.assertEqual((gen["zero_modes"], gen["self_stresses"]), GEN[n], f"generic n={n}")

    def test_stored_file_matches_literal_tables(self):
        stored = json.loads((HERE / "exact_counts.json").read_text())
        rows = {r["n"]: r for r in stored["counts_by_n"]}
        self.assertEqual(sorted(rows), list(range(2, 9)))
        for n, r in rows.items():
            self.assertEqual((r["crosspolytope"]["zero_modes"], r["crosspolytope"]["self_stresses"]), CP[n])
            self.assertEqual((r["q0"]["zero_modes"], r["q0"]["self_stresses"]), Q0[n])
            self.assertEqual(tuple(r["generic_z_minus1_and_pythagorean_zero_modes_selfstresses"][1]), GEN[n])
            self.assertEqual(r["crosspolytope_rigid"], n >= 3)
            self.assertEqual(r["bars_per_cell"], 2 * n * n - n)
            for k, zm, ss in r["planes_k_zeromodes_selfstresses"]:
                if n >= 3:
                    self.assertEqual((zm, ss), (GEN[n][0] + k, GEN[n][1] + k))

    def test_recompute_equals_stored_rows(self):
        stored = json.loads((HERE / "exact_counts.json").read_text())
        fresh = json.loads(RE.dump(RE.build_certificate(nmax=5)))
        self.assertEqual(fresh["counts_by_n"], stored["counts_by_n"][:4])
        self.assertEqual(fresh["control_details"]["delete_one_bar_n3"], stored["control_details"]["delete_one_bar_n3"])

    def test_stress_basis_dimension(self):
        for n, dim in ((3, 0), (4, 2), (5, 5), (6, 9)):
            got, ok = RE.exact_stress_basis_check(n)
            self.assertEqual((got, ok), (dim, True), f"n={n}")

    def test_exact_rank_known_matrices(self):
        Q = RE.QI
        self.assertEqual(RE.rank([[Q(1), Q(0), Q(0)], [Q(0), Q(1), Q(0)], [Q(0), Q(0), Q(1)]]), 3)
        self.assertEqual(RE.rank([[Q(1), Q(2)], [Q(2), Q(4)]]), 1)
        # rows (1, i) and (i, -1): the second is i times the first, so rank 1 (needs exact Q(i) arithmetic)
        self.assertEqual(RE.rank([[Q(1), Q(0, 1)], [Q(0, 1), Q(-1)]]), 1)
        self.assertEqual(RE.rank([[Q(1), Q(0, 1)], [Q(0, 1), Q(1)]]), 2)


class ExactMustFail(unittest.TestCase):
    def test_require_raises_under_any_flags(self):
        with self.assertRaises(SystemExit):
            RE.require(False, "x")
        with self.assertRaises(SystemExit):
            MP.require(False, "x")
        RE.require(True, "x")

    def test_deleting_a_bar_changes_the_counts(self):
        for n in (3, 4, 5):
            pos, ib = RE.positions(n), RE.intra_bars(n)
            base, less = RE.counts(n, pos, ib), RE.counts(n, pos, ib[1:])
            self.assertNotEqual((base["rank"], base["self_stresses"]), (less["rank"], less["self_stresses"]))

    def test_square_is_not_rigid(self):
        cp = RE.counts(2, RE.positions(2), RE.intra_bars(2))
        self.assertNotEqual(cp["zero_modes"], 3)

    def test_sign_twisted_stress_is_rejected(self):
        for n in (4, 5):
            self.assertEqual(RE.exact_stress_basis_check(n, twist=True)[1], False)

    def test_wrong_phase_family_gives_wrong_counts(self):
        # all phases 1 is the q = 0 point: it must NOT give the generic counts
        n = 4
        c = RE.counts(n, RE.positions(n), RE.intra_bars(n), [RE.ONE] * n, True)
        self.assertNotEqual((c["zero_modes"], c["self_stresses"]), GEN[n])

    def test_must_fail_helper_rejects_an_accepting_verdict(self):
        with self.assertRaises(SystemExit):
            RE.must_fail("broken input was accepted", True, [])
        log = []
        RE.must_fail("broken input was rejected", False, log)
        self.assertEqual(log, ["broken input was rejected"])

    def test_corrupt_runs_exit_1(self):
        for script in ("rigidity_exact.py", "modp_check.py"):
            r = run(script, "--corrupt", *(("--nmax", "4") if script == "rigidity_exact.py" else ()))
            self.assertEqual(r.returncode, 1, script)
            self.assertIn("MISMATCH", r.stderr + r.stdout, script)


class ModPKnownAnswers(unittest.TestCase):
    def test_square_roots(self):
        for p in MP.PRIMES:
            self.assertEqual(p % 8, 1)
            s2, ii = MP.sqrt_mod(2, p), MP.sqrt_mod(p - 1, p)
            self.assertEqual(s2 * s2 % p, 2)
            self.assertEqual(ii * ii % p, p - 1)

    def test_rank_mod_known_matrices(self):
        p = 998244353
        self.assertEqual(MP.rank_mod([[1, 0], [0, 1]], p), 2)
        self.assertEqual(MP.rank_mod([[1, 2], [2, 4]], p), 1)
        self.assertEqual(MP.rank_mod([[p, 0], [0, p]], p), 0)

    def test_counts_match_literal_tables(self):
        for p in MP.PRIMES:
            c = (2 - MP.sqrt_mod(2, p)) % p
            for n in range(2, 6):
                rows, sq = MP.system(n, p, c)
                self.assertEqual(MP.counts(rows, n, p), CP[n], f"p={p} cross-polytope n={n}")
                rows, sq = MP.system(n, p, c, [1] * n)
                self.assertEqual(MP.counts(rows, n, p), Q0[n], f"p={p} q=0 n={n}")
                self.assertEqual(len(sq), 1, "all bars one length")
                rows, _ = MP.system(n, p, c, [p - 1] * n)
                self.assertEqual(MP.counts(rows, n, p), GEN[n], f"p={p} z=-1 n={n}")

    def test_agrees_with_the_exact_method(self):
        p = MP.PRIMES[0]
        c = (2 - MP.sqrt_mod(2, p)) % p
        for n in (3, 4, 5):
            ex = RE.counts(n, RE.positions(n), RE.intra_bars(n), [RE.QI(-1)] * n, True)
            rows, _ = MP.system(n, p, c, [p - 1] * n)
            self.assertEqual(MP.counts(rows, n, p), (ex["zero_modes"], ex["self_stresses"]))


class ModPMustFail(unittest.TestCase):
    def test_must_fail_helper_rejects_an_accepting_verdict(self):
        with self.assertRaises(SystemExit):
            MP.must_fail("broken input was accepted", True, [])
        log = []
        MP.must_fail("broken input was rejected", False, log)
        self.assertEqual(log, ["broken input was rejected"])

    def test_wrong_c_breaks_equal_length(self):
        p = MP.PRIMES[0]
        _, sq = MP.system(4, p, 1, [1] * 4)
        self.assertGreater(len(sq), 1)

    def test_moved_joint_breaks_equal_length(self):
        p = MP.PRIMES[0]
        c = (2 - MP.sqrt_mod(2, p)) % p
        _, sq = MP.system(4, p, c, off_axis=True)
        self.assertGreater(len(sq), 1)

    def test_deleting_a_bar_changes_the_counts(self):
        p = MP.PRIMES[1]
        c = (2 - MP.sqrt_mod(2, p)) % p
        for n in (3, 4):
            base = MP.counts(MP.system(n, p, c)[0], n, p)
            less = MP.counts(MP.system(n, p, c, drop=1)[0], n, p)
            self.assertNotEqual(base, less)


if __name__ == "__main__":
    unittest.main(verbosity=2)
