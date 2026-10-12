# Incentre lattice

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23322079.svg)](https://doi.org/10.5281/zenodo.23322079)

One abstract drawing, a hunch followed, a thread unravelled. (Where the drawing came from is told on the Drift page, `drift.html`.)

Draw a square and its two diagonals. In each of the four triangles, mark the point that is the same distance from all three sides, and join it to the triangle's corners. That one rule, carried into a cube and then into any number of dimensions, lands on uniform honeycombs, on the radius-ratio limits of chemistry, on the boron framework of calcium hexaboride, and on a framework whose rigidity has a closed form in every dimension.

This repository is an open, agent-assisted research project around that construction. It holds:

- **An interactive site** (`lattice.html` is the hub; live at the GitHub Pages link in the sidebar). Six subpages explain the drawing, the mirror geometry, the chemistry, the folding mechanics, the perception question, and a ledger of every claim with its evidence label.
- **Research notes** in `research/`, with the code behind every number.
- **Open tracks** in `tracks/`: questions with a fixed, trusted checker and a single score, so that people and AI agents can attack them in parallel. See [CONTRIBUTING.md](CONTRIBUTING.md).

## What is known and what is ours

Most of the geometry is classical. The dot is the Wythoff point of the cube's reflection group (Wythoff 1918, Coxeter), its orbit is the truncated square tiling in 2D and the truncated cubic honeycomb in 3D, the 3D arrangement is Pauling and Weinbaum's CaB6 boron net (RCSR `cab`), and the 2D fold is a known hinged tiling (Flores, after Stuart 1963). The results we have not found in print, after web searches only, are: the rigidity counts and Poisson ratio of the contact framework in every dimension, a single radius-ratio identity, and the closed form of the CaB6 hard-sphere packing maximum. The fold of the 3D framework to close packing (0.2233 to π/√18 in hcp, 2×2×2 family, two checkers) is the other result; no prior description was found (searches listed in `autoresearch-fold/`). The closest known relative is the ReO3 to VF3 tilt of corner-sharing octahedra toward hcp (O'Keeffe & Hyde 1977); in the cab framework that tilt stalls at density 0.4686, and reaching close packing needs a symmetry-breaking continuation with lattice shear. Two rigid layer slides then reach fcc (accepted by two independent checkers). This is a geometry result, not a predicted phase transition in real CaB6 or Prussian-blue materials. The ledger page (`open.html`) labels every claim as PROVEN, VERIFIED NUMERICALLY, CONJECTURE or KNOWN.

Nothing here has been peer reviewed.

## Paper

Preprint (not peer reviewed): **The incentre of the centre-cut pyramid: rigidity of the truncated n-cubic contact framework, a radius-ratio identity, and the CaB<sub>6</sub> packing maximum**, Ahraz Arifuddin (Aska Digital), 2026. [PDF](papers/incentre-lattice.pdf) · [LaTeX source](papers/source/) · DOI [10.5281/zenodo.23322079](https://doi.org/10.5281/zenodo.23322079) (all versions; v1.0.0 is [10.5281/zenodo.23322080](https://doi.org/10.5281/zenodo.23322080)). Its SHA-256 hash is recorded in the public commitment (see below).

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
file; `commitment/` holds a SHA-256 manifest and its root hash, with instructions to check it ([commitment/README.md](commitment/README.md)). The commitment was posted publicly on 2026-10-12 ([gist](https://gist.github.com/ahrazzle/46adf76f15bf766e24b8df626f9c5706)).

## Credit, notices, citation

[CONTRIBUTORS.md](CONTRIBUTORS.md) (including AI assistance and which models did which roles), [NOTICE](NOTICE) (third-party data and licences; note the mixed licences of the
COD extracts), [CITATION.cff](CITATION.cff) (cites the preprint). Pull requests use `.github/pull_request_template.md`; reviews use the "Research review" issue template.

## Licence

Code: MIT. Text, figures and site content: CC BY 4.0. Credit: [@ahrazzle](https://github.com/ahrazzle) (Ahraz Arifuddin), [@aska-digital](https://github.com/aska-digital) (Aska Digital) and [@claude](https://github.com/claude) (Claude Opus 5.5 from Anthropic).
