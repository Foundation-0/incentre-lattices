#!/usr/bin/env python3
"""Check a sha256sum-format file (two-space separator), paths relative to a base directory.

Usage: python tools/sha256check.py SUMS_FILE [BASE_DIR] [--expect NAME ...]
  BASE_DIR defaults to the repository root. With --expect, the listed paths must be exactly NAME ...
It fails before hashing anything if the file has no entries, if a non-blank line is not 64 lowercase
hex digits, two spaces and a path, or if a path is absolute, resolves outside BASE_DIR or names a file
already listed. Blank lines are ignored.
Standard library only, so it works wherever Python does (no need for sha256sum or shasum).
"""
import argparse
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import certlib  # noqa: E402

LINE = re.compile(r"([0-9a-f]{64})  (.+)")


def main(sums, base, expect=None):
    base_dir = pathlib.Path(base).resolve()
    entries, seen = [], set()
    for n, ln in enumerate(pathlib.Path(sums).read_text().splitlines(), 1):
        if not ln.strip():
            continue
        m = LINE.fullmatch(ln)
        if not m:
            raise SystemExit(f"{sums}:{n}: expected '<64 lowercase hex digits>  <path>', got {ln!r}")
        digest, rel = m.groups()
        if pathlib.PurePath(rel).is_absolute():
            raise SystemExit(f"{sums}:{n}: absolute path {rel!r}")
        p = (base_dir / rel).resolve()
        if not p.is_relative_to(base_dir):
            raise SystemExit(f"{sums}:{n}: path {rel!r} resolves outside {base_dir}")
        if p in seen:
            raise SystemExit(f"{sums}:{n}: {rel!r} is listed twice")
        seen.add(p)
        entries.append((digest, rel, p))
    if not entries:
        raise SystemExit(f"{sums}: no checksums listed")
    listed = {rel for _, rel, _ in entries}
    if expect is not None and listed != set(expect):
        raise SystemExit(f"{sums}: listed paths differ from --expect: "
                         f"missing {sorted(set(expect) - listed)}, unexpected {sorted(listed - set(expect))}")
    bad = 0
    for digest, rel, p in entries:
        ok = p.is_file() and certlib.sha256_file(p) == digest
        bad += not ok
        print(("OK      " if ok else "FAILED  ") + rel)
    if bad:
        raise SystemExit(f"{bad} of {len(entries)} checksums failed")
    print(f"{len(entries)} checksums OK")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sums", metavar="SUMS_FILE")
    ap.add_argument("base", metavar="BASE_DIR", nargs="?", default=certlib.ROOT)
    ap.add_argument("--expect", metavar="NAME", nargs="+", help="the exact set of paths the file must list")
    args = ap.parse_args()
    main(args.sums, args.base, args.expect)
