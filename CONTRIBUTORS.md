# Contributors and credit

## Originator and direction

**Ahraz Arifuddin, Aska Digital** originated the construction (the incentre-and-spokes reading of the 4.8.8 tiling and its generalisation
to pyramids over the facets of a polytope), chose the questions, and directed the project. He decides what is claimed and what is published.

## AI assistance

Most of the computation, proofs, checking and writing in this repository was done by AI models under his direction: **Claude Opus 5.5
and other Claude models by Anthropic, PBC.** Their roles, in general terms:

| Model family | Role |
|---|---|
| Claude Opus | proofs and derivations, research and novelty judgment, prior-art review, decisions about what a result can honestly be called |
| Claude Sonnet | builds and checks: verifiers, replay scripts, tests, certificates, continuous integration, documentation |
| Claude Haiku | operations: running known scripts, searching and retrieval, file moves and repository housekeeping |

The October 2026 site update was prepared with **GPT-6** (research synthesis, interactive fcc illustration and content),
**GPT-6.1 Sol** (shared site controls, claim ledger, accessibility fixes and code review), and **GPT-6 Luna** (link checks and repository publication).
The research update separately identifies AI contributions to the open result and paper drafts. Publishing their summaries on the site does not constitute acceptance or outside mathematical review.

OpenAI Codex (GPT-6) prepared the exact CaB6 attainment proof in `results/cab6-packing/ATTAINMENT.md`; GPT-6 Luna handled its repository submission. No independent human review is claimed for this contribution.

No AI model is an author in the legal or scholarly sense, and none can vouch for a result. The labels PROVEN, VERIFIED and KNOWN in
`docs/audit.md` say what was actually done; a person has not independently re-derived every proof. Internal reviews so far were performed
by AI models under direction (see `docs/review-ledger.md`); that is "internal pre-publication review" and not outside review.

For the explicit continuous fcc motion (`results/fold-continuous-fcc`), GPT-6 (Codex) derived the curve, written proof
and first exact arithmetic implementation. GPT-6.1 Sol independently checked the analytic reductions and supplied
a separate polynomial-convolution replay. GPT-6 Luna handled repository publishing. These are internal AI
contributions under the user's direction; independent human review and historical priority remain open.

The `Co-authored-by` trailers on the squash-merge commits on `main` are copied from the commits of each pull request, so they name
only models whose commits carried a trailer. Commits made with Claude Code add one; the commits of work that the PR bodies credit to
GPT-6, GPT-6.1 Sol, GPT-6 Luna or Muse have none (for PR #33, the body of PR #28, which it replaced). So the squash commits of PRs #11,
#23, #25, #27, #31 and #33 name only Claude, and that of PR #30 names no model. The attribution record is this file, the pull-request
bodies, the credit sections of the result notes and `docs/review-ledger.md`, not the trailers.

## Prior work this repository builds on

Credit for the ideas below belongs to the cited people; where we reuse a construction we say so in the `PROOF.md` of the result.

* The hinged-squares linkage: Antonio Flores ("Hinged tilings"), crediting Donald Stuart (1963); Daniel Wells for the term "hinged tessellations".
* Rigid unit modes of linearly bridged frameworks and the rigid-unit-mode spectrum: Andrew L. Goodwin; Stephen C. Power.
* Periodic tilings and auxetic deployments, including bridged cross-polytopes in arbitrary dimensions: Ciprian S. Borcea and Ileana Streinu.
* Maxwell-Calladine counting: James Clerk Maxwell; C. R. Calladine.
* Sphere packing nets and catalogues: the Reticular Chemistry Structure Resource (RCSR); Fischer, Koch and Sowa on homogeneous sphere packings.
* Packing of the CaB6 type and binary packings: Horst and Travesset; Hopkins, Stillinger and Torquato; Murray and Sanders; Pauling and Weinbaum.
* The ReO3 to VF3 octahedral tilt towards hcp anions: O'Keeffe and Hyde; Glazer's tilt notation.
* Packing and jamming methods: Torquato and Jiao; Donev, Torquato, Stillinger and Connelly.
* Perovskite tolerance-factor data (used in `research/chemistry`): Christopher J. Bartel and co-authors.
* Crystal structure data: the Crystallography Open Database and its contributors.
* The shape of this repository (hash-pinned certificates, a CI matrix that reruns every verifier, scoped claims, review ledger) follows
  CrocSwap/integer-mult-bounds (Douglas Colkitt and community contributors) and Swapnil-jain/gf2m-circuits (per-result folders,
  standard-library verifier, must-fail controls, public commitment of artefacts). We thank them for setting the standard.

Full citations are in the `PROOF.md` of each result. See `NOTICE` for third-party material and licences.

## How to be added

Open a pull request using `.github/pull_request_template.md`. Say what you contributed, disclose AI assistance, and your credit will be recorded here
and in `docs/review-ledger.md` for reviews.
