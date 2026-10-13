# Review ledger

One row per review of a result. Add a row for every pull request that changes a result, a checker or a claim, and for every external review.
Verdicts: **confirmed**, **confirmed with corrections**, **partially confirmed**, **refuted**, **open**. "Internal pre-publication review" means a review done
inside the project before anything was public, by AI models working under direction; it is not independent outside review.

## Ledger

| PR | Result | Reviewed commit | Verdict | Notes |
|---|---|---|---|---|
| [PR #25](https://github.com/aska-digital/incentre-lattice/pull/25) | certificate commitments | `2ed2ab7` | open | Registry and stored-path coverage omission recorded in [issue #24](https://github.com/aska-digital/incentre-lattice/issues/24). Eleven regression tests and full `make verify` pass; implementation by GPT-6.1 Sol, reviewed by GPT-6, repository submission by GPT-6 Luna. |
| (internal) | rigidity-nd | `18dcb6e` (working tree) | confirmed with corrections | Internal pre-publication review `VERIFY.md` section C: counts n = 3..6 and nu_100 = 1/(n^2-n-1) re-derived numerically from scratch. n = 2 correctly excluded (2 zero modes, 0 stresses). Elastic tensor is only positive semi-definite. Later extended to n = 8 and to the symbolic proof steps (`results/rigidity-nd`). |
| (internal) | radius-ratio | `18dcb6e` (working tree) | confirmed with corrections | `VERIFY.md` section F: confirmed for isohedral tangential polyhedra by an exact argument and LP, refuted outside (triangular prism; trapezohedra at other heights). Label error fixed: 0.9021 belongs to P = dodecahedron (dots form the icosahedron), 1.8025 to P = icosahedron. Hypotheses corrected to H1 or H2. |
| (internal) | cab6-packing | `18dcb6e` (working tree) | confirmed with corrections | `VERIFY.md` section D: value and closed form exact; "Ca touches B" refuted (gap 0.0297); x = 0.2071 correct only in the Wyckoff convention (offset from centre is 0.2929); global maximum shown unique within the family. Both corrections are in `results/cab6-packing/PROOF.md` and are replay controls. |
| (internal) | fold-2d | `18dcb6e` (working tree) | confirmed with corrections | `VERIFY.md` section A: exact in sympy. Corrections: bars are centre-to-centre (no pin model); the fold is the symmetric path, not the only flex (shear exists); "diagonal mirror line" family has no fold. |
| (internal) | fold-3d-cubic | `18dcb6e` (working tree) | partially confirmed | `VERIFY.md` section B: numbers confirmed; the cubic lattice is an extra constraint and 60 degrees is a lock of the cubic-restricted path, not a jam (trigonal lattices reach 0.65). Wording constrained accordingly. |
| (internal) | prior art, mechanics, chemistry, construction | `18dcb6e` (working tree) | confirmed with corrections | `PRIOR-ART.md` (arXiv, Crossref, Europe PMC, RCSR raw data, web search): KNOWN, PARTIALLY KNOWN or NOT FOUND per claim, with citations; the hinged squares are KNOWN (Flores, Stuart); the n-dimensional closed forms and the radius-ratio theorem NOT FOUND as statements; the 60-degree net matches RCSR sms in its coordination sequence. Absence of evidence only. |
| (internal) | fold-close-packing (hcp) | `18dcb6e` (working tree) | confirmed | Independent from-scratch checker `autoresearch-fold/independent/chk.py` (and `cna.py`): bars, overlaps over all images, step sizes, 10-sub-state interpolation, contacts, CNA. Same verdict as the locked `verify.py`. Reproduced by `results/fold-close-packing/chk_independent.py`. |
| (internal) | fold-close-packing (fcc) | `18dcb6e` (working tree) | confirmed (new) | `autoresearch-fold/PRIOR-ART-FOLD.md` (opus novelty judgment): fcc extension accepted by `verify.py` with a CNA by the same author, "only one checker has looked at this". Closed on this repository's build: the independent checker accepts the 1193-state path (zero overlaps, bars below 1e-9, pure fcc 576 (4,2,1)). Still the same kind of evidence, not outside review. |
| (internal) | fold-close-packing (novelty) | `18dcb6e` (working tree) | open | `PRIOR-ART-FOLD.md`: NOT FOUND is not "new". Library checks listed as not done (Sowa deformation papers, Koch and Fischer on 5/3 types, O'Keeffe and Hyde *Crystal Structures I*, Conway and Sloane chapter 1). |
| (internal) | all | working tree of the build that added `results/`, `certificates/`, CI | confirmed | Fresh virtual environment from `requirements.txt`: `make verify` passes; certificates are deterministic across reruns. No claim was changed by this build. |
| (internal; issue #26, phase 2) | fold-continuous-fcc | explicit construction added after `f015da0` | confirmed (internal); outside review open | GPT-6 derived the curve and proof; GPT-6.1 Sol separately formed all quartics by convolution and reversed the Bernstein expansion: 3248 cases, 16240 nonnegative coefficients. Analytic review confirmed image bound, bar identity, determinant derivative and fcc cosets. Both existing numerical checkers accepted 501 derived samples. No historical novelty or independent human review claimed. |

The unpublished source documents `VERIFY.md` and `PRIOR-ART.md` live in the paper folder, which is not public yet; their findings are incorporated in the
`PROOF.md` files and in this table.

## How to add a row

1. Link the pull request (or write "external" and link the issue or message).
2. Name the result folder and the commit you reviewed (`git rev-parse HEAD` at review time).
3. Give one verdict from the list above and put the reason in Notes: what you re-derived, with which tool, what you found, and what you did NOT check.
4. Keep finite checks and proofs apart: "reran the verifier" is not "checked the proof".
5. Say who reviewed and whether AI tools were used.
