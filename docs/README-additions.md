# Text for README.md (to be merged by the person editing README.md)

README.md and tracks/README.md are being edited elsewhere, so this text is kept here. Paste the sections you want. Paths are relative to the repository root.

---

## Results

Every result below has a folder in `results/` with a `PROOF.md` (statement, scope, evidence label, proof, prior art, open obligations), a deterministic `verify.py`, and a
`certificate.json`. Labels: **PROVEN** = written proof; **VERIFIED** = finite check, not a proof; **KNOWN** = published and cited. Prior-art entries mean "not found in the
searches listed", not "new".

| Result | Value | Label | Verify | Prior art |
|---|---|---|---|---|
| [n-dimensional rigidity](results/rigidity-nd/PROOF.md) | generic zero modes n(n-1)/2, self-stresses n(n-3)/2 per cell; nu_100 = 1/(n^2-n-1); shear moduli 0 | PROVEN n >= 3; VERIFIED n = 2..8 | `make verify-rigidity` | 3D case KNOWN (Goodwin 2006); n-D closed forms not found |
| [Radius ratio](results/radius-ratio/PROOF.md) | hole/ball of the dot construction = classical radius ratio, under mirror walls or insphere plus equal dihedral angles | PROVEN under H1 or H2; VERIFIED on 20 polytopes | `make verify-radius-ratio` | not found as a statement; ingredients known |
| [CaB6 packing](results/cab6-packing/PROOF.md) | maximum (pi/6)(30 sqrt2 - 41) = 0.746865 at gamma = sqrt2 - 1; beats fcc for 0.4102 < gamma < 0.4159 | PROVEN upper bound; attainment VERIFIED | `make verify-cab6` | curve plotted by Horst and Travesset 2016; closed form not found |
| [2D fold](results/fold-2d/PROOF.md) | density 0.5390 to pi/4 at theta = pi/4 | KNOWN linkage (Flores; Stuart); corollaries PROVEN | `make verify-fold2d` | known mechanism |
| [3D fold, cubic-restricted](results/fold-3d-cubic/PROOF.md) | 27 pi sqrt2/256 = 0.4686 at 60 degrees; a lock of the cubic path, not a jam | VERIFIED | `make verify-fold3d` | Prussian-blue tilts known; end net matches RCSR sms (ten terms) |
| [Fold to close packing](results/fold-close-packing/PROOF.md) | pi/sqrt18 = 0.74048: sampled bar-preserving, overlap-free path in the 2x2x2 supercell to hcp (793 states) and fcc (1193 states) | VERIFIED (sampled continuity, step <= 0.01; two checkers); fcc extension is new work | `make verify-fold-cp` | not found; closest is the ReO3 to VF3 tilt (O'Keeffe and Hyde 1977) |

Only the m = 2 supercell is known to reach close packing. The fold is a geometry and rigidity result, not a predicted phase transition of real CaB6 or Prussian-blue materials.

## Reproduce

```
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
make verify          # about 4 to 6 minutes; make verify-fast for the quick targets
```

`make verify` reruns every verifier, including must-fail controls, and the CI matrix (Python 3.11, 3.12, 3.13) fails if a regenerated certificate differs from the committed one.
See [docs/reproducibility.md](docs/reproducibility.md), [docs/audit.md](docs/audit.md) (what each check does and does not show), [docs/status.md](docs/status.md) (current best
per track), [docs/review-ledger.md](docs/review-ledger.md) (reviews so far; internal only) and [docs/formalization.md](docs/formalization.md) (Lean targets, open).

## Certificates and commitment

`certificates/selected-results.json` has one record per result (statement, value, label, verifier command, input hashes); `certificates/SHA256SUMS` lists every certificate and path
file; `commitment/` holds a SHA-256 manifest and its root hash, with instructions to check it ([commitment/README.md](commitment/README.md)). The commitment has not been published.

## Credit, notices, citation

[CONTRIBUTORS.md](CONTRIBUTORS.md) (including AI assistance and which models did which roles), [NOTICE](NOTICE) (third-party data and licences; note the mixed licences of the
COD extracts), [CITATION.cff](CITATION.cff) (a paper preprint is pending). Pull requests use `.github/pull_request_template.md`; reviews use the "Research review" issue template.

## Suggested edits to existing text

* In the fold track description (`tracks/README.md`), the line "second check pending" for the fcc slides can now read: "accepted by the locked checker and by an independent checker
  (`results/fold-close-packing/chk_independent.py`); the open item is a certified, non-sampled path".
* The open-problem "fcc second check" in `tracks/README.md` is done (see `results/fold-close-packing`); the remaining items are unchanged.
* CONTRIBUTING.md labels claims PROVEN, VERIFIED NUMERICALLY, CONJECTURE or KNOWN; this build uses PROVEN, VERIFIED (finite check), KNOWN and NEW-COMPUTATION. Pick one vocabulary.
