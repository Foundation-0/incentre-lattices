#!/usr/bin/env python3
"""Exercise the real certificate generator using an in-memory filesystem only,
and check that `make verify` runs the verifier of every registered result."""
import contextlib
import copy
import hashlib
import io
import json
import pathlib
import re
import types
import unittest
from unittest.mock import patch

ROOT = pathlib.Path(__file__).resolve().parents[1]
builder = types.ModuleType("build_certificates")
builder.__file__ = str(ROOT / "tools/build_certificates.py")
exec(compile(pathlib.Path(builder.__file__).read_text(), builder.__file__, "exec"), builder.__dict__)
OUTPUTS = ("certificates/selected-results.json", "certificates/SHA256SUMS", "commitment/SHA256SUMS", "commitment/ROOT")


class MemoryFiles:
    def __init__(self, files):
        self.files = dict(files)
        self.mkdirs = []

    def read_bytes(self, path):
        try:
            return self.files[pathlib.Path(path)]
        except KeyError:
            raise FileNotFoundError(path) from None

    def read_text(self, path, *args, **kwargs):
        return self.read_bytes(path).decode()

    def write_text(self, path, text, *args, **kwargs):
        self.files[pathlib.Path(path)] = text.encode()
        return len(text)

    def exists(self, path):
        path = pathlib.Path(path)
        return path in self.files or any(path in p.parents for p in self.files)

    def is_file(self, path):
        return pathlib.Path(path) in self.files

    def glob(self, path, pattern):
        return [p for p in self.files if path in p.parents and p.relative_to(path).match(pattern)] if path == ROOT / "results" else []

    def open(self, path, mode="r", *args, **kwargs):
        if mode != "rb":
            raise AssertionError(f"unexpected open mode: {mode}")
        return io.BytesIO(self.read_bytes(path))

    @contextlib.contextmanager
    def mounted(self):
        with contextlib.ExitStack() as stack:
            stack.enter_context(patch.object(pathlib.Path, "read_text", lambda p, *a, **k: self.read_text(p, *a, **k)))
            stack.enter_context(patch.object(pathlib.Path, "write_text", lambda p, t, *a, **k: self.write_text(p, t, *a, **k)))
            stack.enter_context(patch.object(pathlib.Path, "exists", lambda p: self.exists(p)))
            stack.enter_context(patch.object(pathlib.Path, "is_file", lambda p: self.is_file(p)))
            stack.enter_context(patch.object(pathlib.Path, "glob", lambda p, g: self.glob(p, g)))
            stack.enter_context(patch.object(pathlib.Path, "mkdir", lambda p, *a, **k: self.mkdirs.append(p)))
            stack.enter_context(patch("builtins.open", self.open))
            stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
            yield


class CertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Discover fixtures independently of the generator's registries.
        cls.certificates = sorted((ROOT / "results").glob("*/certificate.json"))
        certs = {p: json.loads(p.read_bytes()) for p in cls.certificates}
        # Every recorded input (sources and stored data); `paths` are the stored data files. Inputs get
        # placeholder contents and the fixture certificates record their hashes, so these tests do not depend
        # on the working tree's inputs (make verify may run them before the verifiers rewrite the certificates).
        cls.inputs = sorted({ROOT / rel for cert in certs.values() for rel in cert["inputs"]})
        cls.paths = [p for p in cls.inputs if not p.name.endswith((".py", ".md"))]
        cls.seed = {p: f"contents of {p.relative_to(ROOT).as_posix()}\n".encode() for p in cls.inputs}
        for p, cert in certs.items():
            cert["inputs"] = {rel: hashlib.sha256(cls.seed[ROOT / rel]).hexdigest() for rel in cert["inputs"]}
            cls.seed[p] = (json.dumps(cert, indent=2, sort_keys=True) + "\n").encode()

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
        for suffix in (".npz", ".npy", ".json", ".csv"):
            with self.subTest(suffix=suffix):
                rel = f"results/future-path{suffix}"
                self.fs.files[ROOT / rel] = b"stored data"
                cert = json.loads(self.seed[path])
                cert["inputs"][rel] = hashlib.sha256(b"stored data").hexdigest()
                self.fs.files[path] = json.dumps(cert).encode()
                with self.assertRaisesRegex(SystemExit, f"unregistered path inputs.*{re.escape(rel)}"):
                    builder.main("--write")

    def test_future_source_input_bound_through_certificate(self):
        path = self.certificates[0]
        cert = json.loads(self.seed[path])
        sources = ("results/future-source.py", "results/future-notes.md")
        for rel in sources:
            self.fs.files[ROOT / rel] = b"# source\n"
            cert["inputs"][rel] = hashlib.sha256(b"# source\n").hexdigest()
        self.fs.files[path] = json.dumps(cert).encode()
        self.generate()
        sums = self.fs.read_text(ROOT / "certificates/SHA256SUMS")
        self.assertEqual([rel for rel in sources if rel in sums], [])

    def test_verifier_without_certificate_rejected(self):
        self.fs.files[ROOT / "results/future-result/verify.py"] = b'print("x")\n'
        self.fs.files[ROOT / "results/future-result/PROOF.md"] = b"# PROOF\n"
        for mode in ("--write", "--check"):
            with self.subTest(mode=mode), self.assertRaisesRegex(
                    SystemExit, r"verifiers without a certificate=\['future-result'\]"):
                builder.main(mode)
        self.assertNotIn(ROOT / "certificates/selected-results.json", self.fs.files)

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
                with self.assertRaisesRegex(SystemExit, "stale certificate|certificates/SHA256SUMS does not match"):
                    builder.main("--check")
                self.fs.files[path] = original

    def test_each_input_tamper_rejected(self):
        self.generate()
        for path in self.inputs:
            rel = path.relative_to(ROOT).as_posix()
            with self.subTest(path=rel):
                original = self.fs.files[path]
                self.fs.files[path] = original + b"\n"
                with self.assertRaisesRegex(SystemExit, f"stale certificate .*: input {re.escape(rel)} changed"):
                    builder.main("--check")
                self.fs.files[path] = original

    def test_source_drift_rejected_without_rewrite(self):
        # One comment line appended to an engine after its certificate was written: neither mode accepts it,
        # and --write does not republish the stale input hash.
        self.generate()
        before = {rel: self.fs.files[ROOT / rel] for rel in OUTPUTS}
        self.fs.files[ROOT / "results/fold-2d/fold2d.py"] += b"# drift\n"
        for mode in ("--check", "--write"):
            with self.subTest(mode=mode), self.assertRaisesRegex(
                    SystemExit, "stale certificate results/fold-2d/certificate.json: input results/fold-2d/fold2d.py"):
                builder.main(mode)
        self.assertEqual({rel: self.fs.files[ROOT / rel] for rel in OUTPUTS}, before)

    def test_missing_input_rejected(self):
        self.generate()
        for path in self.inputs:
            rel = path.relative_to(ROOT).as_posix()
            with self.subTest(path=rel):
                original = self.fs.files.pop(path)
                with self.assertRaisesRegex(SystemExit, f"missing {re.escape(rel)}, an input of"):
                    builder.main("--check")
                self.fs.files[path] = original

    def test_missing_registered_path_file_rejected(self):
        extra = "results/extra-path.npz"
        builder.SUM_FILES.append(extra)
        self.fs.files[ROOT / extra] = b"stored data"
        self.generate()
        del self.fs.files[ROOT / extra]
        for mode in ("--check", "--write"):
            with self.subTest(mode=mode), self.assertRaisesRegex(SystemExit, f"missing files: .*{re.escape(extra)}"):
                builder.main(mode)

    def test_each_generated_output_tamper_rejected(self):
        self.generate()
        for rel in OUTPUTS:
            with self.subTest(path=rel):
                path = ROOT / rel
                original = self.fs.files[path]
                self.fs.files[path] = original + b"tampered"
                with self.assertRaises(SystemExit):
                    builder.main("--check")
                self.fs.files[path] = original

    def test_each_generated_output_missing_rejected(self):
        self.generate()
        for rel in OUTPUTS:
            with self.subTest(path=rel):
                original = self.fs.files.pop(ROOT / rel)
                with self.assertRaisesRegex(SystemExit, f"missing {re.escape(rel)} \\(run make certificates\\)"):
                    builder.main("--check")
                self.fs.files[ROOT / rel] = original

    def test_check_has_no_side_effects(self):
        with self.assertRaisesRegex(SystemExit, "missing certificates/selected-results.json"):
            builder.main("--check")
        self.assertEqual((self.fs.files, self.fs.mkdirs), (self.seed, []))
        self.generate()
        before = dict(self.fs.files)
        self.fs.mkdirs.clear()
        builder.main("--check")
        self.assertEqual((self.fs.files, self.fs.mkdirs), (before, []))


def make_rules(text):
    """Variables and rules of a plain Makefile: {name: words}, {target: (prerequisites, recipe lines)}."""
    variables, rules, target = {}, {}, None
    for line in text.splitlines():
        if line.startswith("\t") and target is not None:
            rules[target][1].append(line.strip())
            continue
        line, target = line.split("#", 1)[0].rstrip(), None
        if m := re.fullmatch(r"(\w+)\s*[:?]?=\s*(.*)", line):
            variables[m[1]] = m[2].split()
        elif m := re.fullmatch(r"([^\s:=]+)\s*:(?!=)\s*(.*)", line):
            target = m[1]
            rules[target] = (m[2].split(), [])
    return variables, rules


def unverified(makefile, results):
    """Registered results whose verifier `make verify` does not run, or whose verify_command does not name
    a prerequisite of `verify` that runs it (empty when every result is covered)."""
    variables, rules = make_rules(makefile)

    def expand(words):
        out = []
        for w in words:
            m = re.fullmatch(r"\$\((\w+)\)", w)
            out += expand(variables[m[1]]) if m else [w]
        return out

    targets, todo = [], expand(rules["verify"][0])
    while todo:
        t = todo.pop(0)
        if t not in targets and t in rules:
            targets.append(t)
            todo += expand(rules[t][0])
    problems = []
    for r in results:
        runs = re.compile(rf"\$\(PY\)( -\S+)* results/{re.escape(r['id'])}/verify\.py(\s|$)")
        if not any(runs.search(line) for t in targets for line in rules[t][1]):
            problems.append(f"{r['id']}: make verify does not run results/{r['id']}/verify.py")
        named = re.match(r"make (\S+)", r["verify_command"])
        if not (named and named[1] in targets and any(runs.search(line) for line in rules[named[1]][1])):
            problems.append(f"{r['id']}: {r['verify_command']!r} is not a prerequisite of verify that runs it")
    return problems


class MakefileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.makefile = (ROOT / "Makefile").read_text()

    def test_make_verify_runs_every_registered_verifier(self):
        self.assertEqual(unverified(self.makefile, builder.RESULTS), [])

    def test_dropped_verify_target_detected(self):
        # Must-fail control: verify-fold-continuous removed from the FAST or SLOW list.
        text = re.sub(r"(?m)^((?:FAST|SLOW) =.*?) verify-fold-continuous(?=\s|$)", r"\1", self.makefile)
        self.assertTrue(text != self.makefile, "control expects verify-fold-continuous in FAST or SLOW")
        self.assertEqual([p for p in unverified(text, builder.RESULTS) if p.startswith("fold-continuous-fcc:")], [
            "fold-continuous-fcc: make verify does not run results/fold-continuous-fcc/verify.py",
            "fold-continuous-fcc: 'make verify-fold-continuous' is not a prerequisite of verify that runs it"])

    def test_wrong_verify_command_detected(self):
        results = copy.deepcopy(builder.RESULTS)
        next(r for r in results if r["id"] == "fold-2d")["verify_command"] = "make verify-cab6"
        self.assertEqual([p for p in unverified(self.makefile, results) if p.startswith("fold-2d:")],
                         ["fold-2d: 'make verify-cab6' is not a prerequisite of verify that runs it"])


if __name__ == "__main__":
    unittest.main()
