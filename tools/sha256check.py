#!/usr/bin/env python3
"""Check a sha256sum-format file (two-space separator), paths relative to a base directory.

Usage: python tools/sha256check.py SUMS_FILE [BASE_DIR]     (BASE_DIR defaults to the repository root)
Standard library only, so it works wherever Python does (no need for sha256sum or shasum).
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import certlib  # noqa: E402


def main(sums, base):
    bad = 0
    lines = [ln for ln in pathlib.Path(sums).read_text().splitlines() if ln.strip()]
    for ln in lines:
        digest, _, rel = ln.partition("  ")
        p = pathlib.Path(base) / rel
        ok = p.exists() and certlib.sha256_file(p) == digest
        bad += not ok
        print(("OK      " if ok else "FAILED  ") + rel)
    if bad:
        raise SystemExit(f"{bad} of {len(lines)} checksums failed")
    print(f"{len(lines)} checksums OK")


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2] if len(sys.argv) == 3 else certlib.ROOT)
