## Change and scope

<!-- What does this change, and which result(s) does it touch? Name the folder(s) under results/ or the track. One result per pull request. -->

- Result(s) touched:
- Type: new result / correction / new check / documentation / tooling
- Does it edit a checker (`verify.py`, `chk_independent.py`, `autoresearch-fold/verify.py`)? Checker edits are reviewed on their own and are not merged together with a result.

## Validation

<!-- Say exactly what was run. Keep finite tests and proofs apart. -->

- [ ] `make verify` passes locally and `git diff --exit-code -- results certificates commitment` is empty (or the change to the certificates is explained below)
- [ ] Finite checks (computations that replay but prove nothing beyond the cases run): list them, with ranges and tolerances
- [ ] Proofs (written arguments): link the section of the relevant `PROOF.md`; say which steps are checked only numerically
- [ ] Must-fail controls added or updated where meaningful
- [ ] `docs/status.md`, `docs/audit.md` and `certificates/selected-results.json` updated if a claim or its label changed (run `make certificates`)

## Prior art

<!-- What did you search, where, and what did you find? "Not found" is not "new". Cite the closest work. -->

- Searches run:
- Closest prior work:
- Label for the claim: PROVEN / VERIFIED (finite check) / KNOWN (cite) / CONJECTURE

## AI assistance

<!-- Required. Say whether AI tools were used, which models, and for what (ideas, code, proofs, checking, writing). State what a human checked. -->

- AI used: yes / no
- Models and roles:
- What a person verified:

## Credit and licence

- [ ] Every idea, construction and dataset I build on is credited (NOTICE, CONTRIBUTORS.md, citations)
- [ ] No private data, credentials or local paths are included
