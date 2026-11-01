# Corrections for the published v1 and v2 archive descriptions

These corrections are prepared for review. They do not edit a published release or Zenodo record, replace an archived PDF, or issue a new DOI.

## v1.0.0 known replay issue

On some machines `make verify` in the v1 archive stops at `verify-rigidity` with `VERIFY FAILED: n=7, c=0.3000: elastic constants`. The numerical relaxation used the default cutoff of `numpy.linalg.pinv`, which could retain a round-off singular value. The exact formulas are unaffected. PR #8 (landed `81aa067`) fixes the numerical projection; the maintained replay is in v2.0.0 and current main. This is a numerical conditioning issue, rather than evidence that the exact elastic formula fails.

## v2.0.0 result description

The unit-spring elastic result is for every dimension n >= 3. C44 = 0 means that every stiffness component involving a cube-axis shear strain vanishes. The axial ratio is C11/C12 = (n-1)^2, and the Poisson ratio for uniaxial stress along a cube axis is nu_100 = 1/(n^2-n-1). This wording does not claim that every rotated shear experiment has zero stiffness.

Borcea and Streinu (2021), Section 7, describe the same infinite framework at a different bar length and with a different periodicity lattice. The counts here are per cell of 2Z^n. The specific formulas and explicit packing motion retain the priority qualifications in the repository; a failed literature search does not establish novelty.

## AI assistance and review

Computation, proofs, checking and writing used Claude models (Anthropic), GPT-6 via Codex, GPT-6.1 Sol and GPT-6 Luna (OpenAI), under the author's direction. The roles are recorded in CONTRIBUTORS.md. Unconfirmed provider attribution for the historical name Muse and unrecorded human proof-checking scope are not inferred from account names or co-author trailers. Independent human mathematical review remains pending.

## Licence wording for both archive records

Code is MIT licensed. Text, figures and site content are CC BY 4.0. Third-party material retains its own terms in NOTICE, including 25 COD CIFs provided by IUCr Journals for scientific use with attribution. A single archive-level licence field does not relicense these files.

The published files retain their original hashes. Whether future source archives retain the 25 IUCr-provided CIFs requires a separate distribution decision. This correction does not remove files or add export rules.
