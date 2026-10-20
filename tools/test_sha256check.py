#!/usr/bin/env python3
"""Exercise tools/sha256check.py on manifests and files in a temporary directory only."""
import contextlib
import hashlib
import io
import os
import pathlib
import subprocess
import sys
import tempfile
import types
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tools/sha256check.py"
checker = types.ModuleType("sha256check")
checker.__file__ = str(SCRIPT)
exec(compile(SCRIPT.read_text(), checker.__file__, "exec"), checker.__dict__)

MALFORMED = "expected '<64 lowercase hex digits>  <path>'"


def line(rel, data):
    return f"{hashlib.sha256(data).hexdigest()}  {rel}"


class Sha256CheckTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.tmp = pathlib.Path(tmp.name)
        self.base = self.tmp / "base"
        (self.base / "sub").mkdir(parents=True)
        self.data = {"a.txt": b"first\n", "sub/b.txt": b"second\n"}
        for rel, data in self.data.items():
            (self.base / rel).write_bytes(data)
        self.outside = self.tmp / "outside.txt"
        self.outside.write_bytes(b"outside\n")
        self.sums = self.tmp / "SUMS"

    def check(self, lines, expect=None):
        self.sums.write_text("".join(ln + "\n" for ln in lines))
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            checker.main(str(self.sums), str(self.base), expect)
        return out.getvalue()

    def good_lines(self):
        return [line(rel, data) for rel, data in self.data.items()]

    def test_valid_manifest_passes(self):
        self.assertTrue(self.check(self.good_lines()).endswith("2 checksums OK\n"))
        self.assertTrue(self.check(["", *self.good_lines(), "  "], expect=["sub/b.txt", "a.txt"]).endswith("2 checksums OK\n"))

    def test_empty_manifest_rejected(self):
        for lines in ([], ["", "   "]):
            with self.subTest(lines=lines), self.assertRaisesRegex(SystemExit, "no checksums listed"):
                self.check(lines)

    def test_malformed_line_rejected(self):
        digest = hashlib.sha256(self.data["a.txt"]).hexdigest()
        for bad in ("deadbeef", f"{digest} a.txt", f"{digest}\ta.txt", f"{digest.upper()}  a.txt",
                    f"{digest[:-1]}  a.txt", f"{digest}0  a.txt", f"{digest}  ", f" {digest}  a.txt"):
            with self.subTest(line=bad), self.assertRaisesRegex(SystemExit, MALFORMED):
                self.check([self.good_lines()[0], bad])

    def test_absolute_path_rejected(self):
        for path, data in ((self.outside, b"outside\n"), (self.base / "a.txt", self.data["a.txt"])):
            with self.subTest(path=path), self.assertRaisesRegex(SystemExit, "absolute path"):
                self.check([line(str(path), data)])

    def test_escaping_path_rejected(self):
        for rel in ("../outside.txt", "sub/../../outside.txt"):
            with self.subTest(rel=rel), self.assertRaisesRegex(SystemExit, "resolves outside"):
                self.check([line(rel, b"outside\n")])
        try:
            os.symlink(self.outside, self.base / "link.txt")
        except OSError:
            self.skipTest("symbolic links not available")
        with self.assertRaisesRegex(SystemExit, "resolves outside"):
            self.check([line("link.txt", b"outside\n")])

    def test_duplicate_path_rejected(self):
        for second in ("a.txt", "./a.txt", "sub/../a.txt"):
            with self.subTest(second=second), self.assertRaisesRegex(SystemExit, "listed twice"):
                self.check([line("a.txt", self.data["a.txt"]), line(second, self.data["a.txt"])])

    def test_expected_set_mismatch_rejected(self):
        # A deleted line no longer drops a file from the lock silently.
        with self.assertRaisesRegex(SystemExit, r"missing \['sub/b.txt'\], unexpected \[\]"):
            self.check(self.good_lines()[:1], expect=["a.txt", "sub/b.txt"])
        with self.assertRaisesRegex(SystemExit, r"missing \[\], unexpected \['sub/b.txt'\]"):
            self.check(self.good_lines(), expect=["a.txt"])

    def test_changed_missing_or_directory_entry_fails(self):
        good = self.good_lines()
        for case, first in (("changed", line("a.txt", b"changed\n")), ("missing", line("gone.txt", b"")),
                            ("directory", line("sub", b""))):
            with self.subTest(case=case), self.assertRaisesRegex(SystemExit, "^1 of 2 checksums failed$"):
                self.check([first, good[1]])

    def test_command_line(self):
        self.sums.write_text("".join(ln + "\n" for ln in self.good_lines()))
        cmd = [sys.executable, "-B", str(SCRIPT), str(self.sums), str(self.base), "--expect"]
        ok = subprocess.run(cmd + ["a.txt", "sub/b.txt"], capture_output=True, text=True)
        self.assertEqual(ok.returncode, 0, ok.stderr)
        self.assertTrue(ok.stdout.endswith("2 checksums OK\n"))
        bad = subprocess.run(cmd + ["a.txt"], capture_output=True, text=True)
        self.assertNotEqual(bad.returncode, 0)
        self.assertIn("listed paths differ from --expect", bad.stderr)


if __name__ == "__main__":
    unittest.main()
