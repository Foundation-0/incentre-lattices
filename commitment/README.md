# Commitment

`SHA256SUMS` lists the sha256 of every file that carries a result: each `results/*/certificate.json`, everything in `certificates/` (the index of selected results and its own checksum list),
and all four stored path files: `autoresearch-fold/best/path.npz` and `results/fold-close-packing/path_fcc.npz` (m = 2 hcp / fcc),
and `results/fold-m1/path_hcp_m1.npz` and `results/fold-m1/path_fcc_m1.npz` (m = 1 hcp / fcc). `ROOT` is the sha256 of the `SHA256SUMS` file itself,
so one 64-character string commits to all of it.

Purpose: a public record of the exact artefacts bound by a root at a given time. A public timestamp binds only the root and files in the posted copy; regenerating this folder does not timestamp the new root.

## Verify

From the repository root, with Python 3 (no packages needed):

```
python3 tools/sha256check.py commitment/SHA256SUMS          # every listed file matches its hash
python3 - <<'PY'
import hashlib, pathlib
root = hashlib.sha256(pathlib.Path("commitment/SHA256SUMS").read_bytes()).hexdigest()
print(root, root == pathlib.Path("commitment/ROOT").read_text().strip())
PY
```

or with standard tools: `shasum -a 256 -c commitment/SHA256SUMS` (macOS) or `sha256sum -c commitment/SHA256SUMS` (Linux), and `shasum -a 256 commitment/SHA256SUMS` compared with `commitment/ROOT`.

To compare against a posted copy, check that its root equals `commitment/ROOT` at the commit being audited, then run the commands above at that commit.

## What it does not show

A matching hash shows that the files are the ones committed to. It does not show that the claims in them are correct (see `docs/audit.md`), and it does not bind anything outside the listed files
(a verifier source is bound only if its hash is among the `inputs` recorded inside a certificate; each certificate's `inputs` lists what it binds). The paper manuscript is not part of this commitment.

## Regenerate (after a deliberate change)

`make certificates` reruns the verifiers and rewrites `certificates/` and this folder. Review the diff: every change to a certificate changes `ROOT`.

## Public timestamp

An earlier commitment was posted 2026-10-12 as a public gist: https://gist.github.com/ahrazzle/46adf76f15bf766e24b8df626f9c5706

The root posted in that gist, `213275970171a80d301fb7e54e262d82cd18611b3ffa21754e7b370b327444bf`, is the `commitment/ROOT` at commit `e84a128`. Later roots, including the expanded root introduced by registering every result certificate and stored path, are not separately timestamped; they do not retroactively timestamp the added certificates and m = 1 paths. To determine what the earlier timestamp binds, check out `e84a128` and compare its `commitment/ROOT` and `SHA256SUMS` with the posted root.

The gist also records the sha256 of the then-unpublished preprint PDF: `9d6dbf351d62c91f17a4b0da10e5e8ef3fe0185374f31586fcf074565068e48a`.

## Published preprint

Preprint v1 was published on 2026-10-12 as `papers/incentre-lattice.pdf` (commit `4d8488e`, DOI 10.5281/zenodo.23322080). Since v2 that file is kept unchanged as `papers/incentre-lattice-v1.pdf`, sha256 `64a59f3c6fdd15419f0819ee897091dcb02ebff1bfc1b7b2a42c1b156404ee26`. It differs from the PDF whose hash was timestamped earlier the same day only by bibliography corrections and the placement of the origin note. Both hashes are recorded in the gist.

`papers/incentre-lattice.pdf` now holds v2 (commit `38bbf45`, sha256 `632d7debf3dc06a55e18de70dfaee90da3b664cea9b5be31de8a02599f72a003`, DOI 10.5281/zenodo.23340326); v2 is not covered by the public timestamp.
