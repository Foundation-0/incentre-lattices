#!/usr/bin/env python3
"""Exercise the real certificate generator using an in-memory filesystem only."""
import contextlib
import copy
import hashlib
import io
import json
import pathlib
import types
import unittest
from unittest.mock import patch

ROOT = pathlib.Path(__file__).resolve().parents[1]
builder = types.ModuleType("build_certificates")
builder.__file__ = str(ROOT / "tools/build_certificates.py")
exec(compile(pathlib.Path(builder.__file__).read_text(), builder.__file__, "exec"), builder.__dict__)


class MemoryFiles:
    def __init__(self, files):
        self.files = dict(files)

    def read_text(self, path, *args, **kwargs):
        return self.files[pathlib.Path(path)].decode()

    def write_text(self, path, text, *args, **kwargs):
        self.files[pathlib.Path(path)] = text.encode()
        return len(text)

    def exists(self, path):
        path = pathlib.Path(path)
        return path in self.files or any(path in p.parents for p in self.files)

    def glob(self, path, pattern):
        return [p for p in self.files if path in p.parents and p.relative_to(path).match(pattern)] if path == ROOT / "results" else []

    def open(self, path, mode="r", *args, **kwargs):
        if mode != "rb":
            raise AssertionError(f"unexpected open mode: {mode}")
        return io.BytesIO(self.files[pathlib.Path(path)])

    @contextlib.contextmanager
    def mounted(self):
        with contextlib.ExitStack() as stack:
            stack.enter_context(patch.object(pathlib.Path, "read_text", lambda p, *a, **k: self.read_text(p, *a, **k)))
            stack.enter_context(patch.object(pathlib.Path, "write_text", lambda p, t, *a, **k: self.write_text(p, t, *a, **k)))
            stack.enter_context(patch.object(pathlib.Path, "exists", lambda p: self.exists(p)))
            stack.enter_context(patch.object(pathlib.Path, "glob", lambda p, g: self.glob(p, g)))
            stack.enter_context(patch.object(pathlib.Path, "mkdir", lambda *a, **k: None))
            stack.enter_context(patch("builtins.open", self.open))
            stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
            yield


class CertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Discover fixtures independently of the generator's registries.
        cls.certificates = sorted((ROOT / "results").glob("*/certificate.json"))
        cls.seed = {p: p.read_bytes() for p in cls.certificates}
        cls.paths = sorted({ROOT / rel for p in cls.certificates
                            for rel in json.loads(cls.seed[p])["inputs"]
                            if rel.endswith(".npz")})
        cls.seed.update({p: p.read_bytes() for p in cls.paths})

    def setUp(self):
        self.fs = MemoryFiles(self.seed)
        self.stack = contextlib.ExitStack()
        self.addCleanup(self.stack.close)
        self.stack.enter_context(patch.object(builder, "RESULTS", copy.deepcopy(builder.RESULTS)))
        self.stack.enter_context(patch.object(builder, "SUM_FILES", list(builder.SUM_FILES)))
        self.stack.enter_context(self.fs.mounted())

    def generate(self):
        builder.main("--write")
        builder.main("--check")

    def test_round_trip_covers_all_certificates_and_paths(self):
        self.generate()
        selected = json.loads(self.fs.read_text(ROOT / "certificates/selected-results.json"))
        self.assertEqual({r["certificate"] for r in selected["results"]},
                         {p.relative_to(ROOT).as_posix() for p in self.certificates})
        required = {p.relative_to(ROOT).as_posix() for p in self.certificates + self.paths}
        for name in ("certificates/SHA256SUMS", "commitment/SHA256SUMS"):
            lines = self.fs.read_text(ROOT / name).splitlines()
            listed = [line.split("  ", 1)[1] for line in lines]
            self.assertTrue(required <= set(listed))
            self.assertEqual(len(listed), len(set(listed)))
            for line in lines:
                digest, rel = line.split("  ", 1)
                self.assertEqual(digest, hashlib.sha256(self.fs.files[ROOT / rel]).hexdigest())
        self.assertEqual(self.fs.read_text(ROOT / "commitment/ROOT").strip(),
                         hashlib.sha256(self.fs.files[ROOT / "commitment/SHA256SUMS"]).hexdigest())

    def test_unregistered_new_certificate_rejected_before_write(self):
        path = ROOT / "results/future-result/certificate.json"
        self.fs.files[path] = b'{"result":"future-result","inputs":{}}'
        for mode in ("--write", "--check"):
            with self.subTest(mode=mode), self.assertRaisesRegex(SystemExit, "unregistered certificates=.*future-result"):
                builder.main(mode)
        self.assertNotIn(ROOT / "certificates/selected-results.json", self.fs.files)

    def test_removed_registry_entry_rejected(self):
        omitted = builder.RESULTS.pop()["id"]
        with self.assertRaisesRegex(SystemExit, f"unregistered certificates=.*{omitted}"):
            builder.main("--write")

    def test_extra_registry_entry_rejected(self):
        builder.RESULTS.append({"id": "missing-result"})
        with self.assertRaisesRegex(SystemExit, "missing certificates=.*missing-result"):
            builder.main("--write")

    def test_duplicate_registry_id_rejected(self):
        builder.RESULTS.append(dict(builder.RESULTS[0]))
        with self.assertRaisesRegex(SystemExit, "duplicate IDs=\\[.+\\]"):
            builder.main("--write")

    def test_missing_certificate_rejected(self):
        del self.fs.files[self.certificates[0]]
        with self.assertRaisesRegex(SystemExit, "missing certificates=\\[.+\\]"):
            builder.main("--write")

    def test_removed_path_registration_rejected(self):
        for path in self.paths:
            with self.subTest(path=path):
                rel = path.relative_to(ROOT).as_posix()
                with patch.object(builder, "SUM_FILES", [p for p in builder.SUM_FILES if p != rel]):
                    with self.assertRaisesRegex(SystemExit, "unregistered path inputs"):
                        builder.main("--write")

    def test_future_path_input_requires_direct_manifest_coverage(self):
        path = self.certificates[0]
        cert = json.loads(self.fs.files[path])
        cert["inputs"]["results/future-path.npz"] = "0" * 64
        self.fs.files[path] = json.dumps(cert).encode()
        with self.assertRaisesRegex(SystemExit, "unregistered path inputs.*future-path.npz"):
            builder.main("--write")

    def test_each_certificate_tamper_rejected(self):
        self.generate()
        for path in self.certificates:
            with self.subTest(path=path):
                original = self.fs.files[path]
                self.fs.files[path] = original + b" "
                with self.assertRaisesRegex(SystemExit, "selected-results.json is out of date"):
                    builder.main("--check")
                self.fs.files[path] = original

    def test_each_path_tamper_rejected(self):
        self.generate()
        for path in self.paths:
            with self.subTest(path=path):
                original = self.fs.files[path]
                self.fs.files[path] = original + b"tampered"
                with self.assertRaisesRegex(SystemExit, "certificates/SHA256SUMS does not match"):
                    builder.main("--check")
                self.fs.files[path] = original

    def test_each_generated_output_tamper_rejected(self):
        self.generate()
        for rel in ("certificates/selected-results.json", "certificates/SHA256SUMS",
                    "commitment/SHA256SUMS", "commitment/ROOT"):
            with self.subTest(path=rel):
                path = ROOT / rel
                original = self.fs.files[path]
                self.fs.files[path] = original + b"tampered"
                with self.assertRaises(SystemExit):
                    builder.main("--check")
                self.fs.files[path] = original


if __name__ == "__main__":
    unittest.main()
