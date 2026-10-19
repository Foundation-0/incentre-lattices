# Reproducibility

Every checked result in `results/` has a deterministic replay. This page says how to run them, what "reproduced" means here, and what is
pinned. For what each replay does and does not establish, read `docs/audit.md`.

## Run everything

```
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt     # numpy, scipy, sympy, mpmath, all pinned
make verify                         # about 4 to 6 minutes on a laptop
git diff --exit-code -- results certificates commitment   # must print nothing
```

`make verify` runs the verifiers below, rewriting each `results/<name>/certificate.json`, then checks that `certificates/selected-results.json`,
`certificates/SHA256SUMS` and `commitment/` still match those files. If anything differs, the command fails. `make verify-fast` runs the quick
targets only (about 30 seconds). `make certificates` is the same as `make verify` but then regenerates the certificate index, the checksum list and
the commitment; run it only after a deliberate change, and review the diff.

| Target | Result | Tools | Time | Slow |
|---|---|---|---|---|
| `verify-rigidity` | `results/rigidity-nd` | numpy, sympy | a few seconds | |
| `verify-radius-ratio` | `results/radius-ratio` | numpy, scipy | a few seconds | |
| `verify-cab6` | `results/cab6-packing` | numpy, scipy, sympy | a few seconds | |
| `verify-fold2d` | `results/fold-2d` | numpy, sympy | seconds | |
| `verify-fold-tests` | the locked fold checker's hash lock and 16 unit tests | numpy, scipy | 12 seconds | |
| `verify-fold3d` | `results/fold-3d-cubic` | numpy, scipy, sympy | about 1 minute | yes |
| `verify-fold-cp` | `results/fold-close-packing` (both checkers, hcp and fcc paths) | numpy, scipy | about 3 minutes | yes |

Times are from a recent laptop. The continuous-integration job has a 15-minute limit and runs `make verify` on Python 3.11, 3.12 and 3.13
(`.github/workflows/verify.yml`).

## What "reproduced" means

Each verifier recomputes its values from scratch (it does not read earlier outputs), asserts every claim it can check, runs must-fail controls
(deliberately broken inputs that the checker must reject), and writes a `certificate.json`. A certificate holds:

* the statement's values, as exact strings where they are exact (`"1/55"`, `"(pi/6)*(30*sqrt(2) - 41)"`) and as numbers rounded to a stated number of
  significant digits where they are numerical;
* tolerance results as booleans (`"bar_error_below_1e-9_relative": true`), never raw round-off residuals, so the file does not depend on the
  last bits of a particular BLAS;
* the names of the must-fail controls that were rejected;
* the sha256 of every input file the verifier read (`inputs`), including the checkers themselves and the stored path files.

`certificates/selected-results.json` has one record per result: statement, value, evidence label, verifier command, input hashes, and the
hash of the certificate. `certificates/SHA256SUMS` lists every certificate and path file. Check it with
`python tools/sha256check.py certificates/SHA256SUMS` (standard library only, no `sha256sum` needed).

A green `make verify` therefore means: the finite computations replay, the certificates are byte-for-byte what the code produces, and the
inputs are the files whose hashes are recorded. It does not mean that any statement beyond those computations is proved.

## Pinned inputs

* `requirements.txt` pins numpy 2.3.5, scipy 1.16.3, sympy 1.14.0 and mpmath 1.3.0. These ship wheels for Python 3.11, 3.12 and 3.13. (An earlier
  pin, scipy 1.15.3, failed to import on a recent macOS because of a linker-section problem in that wheel; scipy 1.16.3 imports fine. Linux was
  not affected, but the pin must work on the machines of people who replay it.)
* `autoresearch-fold/verify.py` is the locked checker. Its sha256 is in `autoresearch-fold/verify.sha256` and is checked by `verify-fold-tests`.
  Do not edit it; changes are reviewed separately from results.
* The stored paths are `autoresearch-fold/best/path.npz` (hcp, 793 states, 0.97 MB) and `results/fold-close-packing/path_fcc.npz` (hcp plus two layer
  slides, 1193 states, 0.77 MB). Together they are 1.7 MB. `path_fcc.npz` is regenerated from the hcp path (`make_fcc_path.py`, no random numbers) in
  every replay and must agree with the stored file to 1e-9; the stored file is what the checkers read. Nothing was decimated.
* Random choices: `rigidity_nd.py` samples random wavevectors with a fixed seed; nothing else uses random numbers. The search that produced the hcp path
  (`autoresearch-fold/search.py`) is not part of `make verify` (it takes minutes and is a search, not a check). It is step-bounded with BLAS pinned to one
  thread, but it is not guaranteed to regenerate the stored path on another machine; `best/path.npz` is the certificate (see `docs/audit.md`).

## Commitment

`commitment/SHA256SUMS` lists the hashes of every certificate, `certificates/*` and all four stored path files; `commitment/ROOT` is the sha256
of that list. `commitment/README.md` explains how to check it. An earlier commitment (root `213275970171...`, commit `e84a128`) was posted
publicly on 2026-10-12 as a gist. The current root covers more files, some files bound by the earlier root have changed since, and it is not
separately timestamped (see `commitment/README.md`).

## Optional: Lean

`make lean` (not part of `make verify`, not run in CI) checks `formal/Arithmetic.lean` with a Lean 4 toolchain, if you have one. It proves only
arithmetic bookkeeping. See `docs/formalization.md`.

## If a replay fails

First check the Python and package versions against `requirements.txt`. If the certificate differs in a digit, that is a finding: open a
"research review" issue (`.github/ISSUE_TEMPLATE/research-review.md`) with the platform, versions and the diff.
