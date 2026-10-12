"""Small helpers shared by every verifier: hashing, rounding, canonical JSON.

Standard library only. Certificates must be byte-identical across reruns, so floats are
recorded to a fixed number of significant digits and tolerance results are recorded as
booleans, never as raw round-off residuals.
"""
import hashlib
import json
import math
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def inputs(*paths):
    """{repo-relative path: sha256} for the files a verifier reads."""
    out = {}
    for p in paths:
        p = pathlib.Path(p).resolve()
        out[p.relative_to(ROOT).as_posix()] = sha256_file(p)
    return dict(sorted(out.items()))


def sig(x, digits=8):
    """Round to `digits` significant digits (0.0 and non-finite pass through)."""
    x = float(x)
    if x == 0.0 or not math.isfinite(x):
        return x
    return float(f"{x:.{digits - 1}e}")


def write_json(path, obj):
    path = pathlib.Path(path)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")


def fail(msg):
    raise SystemExit("VERIFY FAILED: " + msg)


def require(cond, msg):
    if not cond:
        fail(msg)


def must_fail(name, predicate_holds, log):
    """A must-fail control: `predicate_holds` is the verdict of the checker on a deliberately
    broken input. The control passes only if the checker REJECTED the input."""
    require(not predicate_holds, f"control '{name}' was wrongly accepted")
    log.append(name)
