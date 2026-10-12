# Commitment

`SHA256SUMS` lists the sha256 of every file that carries a result: each `results/*/certificate.json`, everything in `certificates/` (the index of selected results and its own checksum list),
and the two path files (`autoresearch-fold/best/path.npz`, the hcp path, and `results/fold-close-packing/path_fcc.npz`, the fcc extension). `ROOT` is the sha256 of the `SHA256SUMS` file itself,
so one 64-character string commits to all of it.

Purpose: a public record that these exact artefacts existed in this form at a given time, to be posted before any paper that relies on them. **It has not been posted anywhere.** Where and when
to post it (a gist, a signed tag, a timestamp service) is for the maintainers to decide; this repository does not publish it for you.

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
(the verifier sources are bound only through the input hashes recorded inside each certificate). The paper manuscript is not part of this commitment.

## Regenerate (after a deliberate change)

`make certificates` reruns the verifiers and rewrites `certificates/` and this folder. Review the diff: every change to a certificate changes `ROOT`.

## Public timestamp

Posted 2026-10-12 as a public gist: https://gist.github.com/ahrazzle/46adf76f15bf766e24b8df626f9c5706

The gist also records the sha256 of the then-unpublished preprint PDF: `9d6dbf351d62c91f17a4b0da10e5e8ef3fe0185374f31586fcf074565068e48a`.

## Published preprint

The preprint was published on 2026-10-12 as `papers/incentre-lattice.pdf`, sha256 `64a59f3c6fdd15419f0819ee897091dcb02ebff1bfc1b7b2a42c1b156404ee26`. It differs from the hash timestamped earlier the same day only by bibliography corrections and the placement of the origin note. Both hashes are recorded in the gist.
