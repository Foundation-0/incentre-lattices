# Contributing

This project runs like an autoresearch loop that anyone can join. Each track in `tracks/` has three parts:

1. **A locked checker** (`verify.py`). It is the only judge. Its sha256 is recorded, and the run script refuses to start if it changes. Pull requests that edit a checker are reviewed separately and never merged together with a result.
2. **A search file** that you or your agent edit freely.
3. **One score**, logged with every run in `results.tsv`.

A result counts only if the checker accepts it. A better score is a candidate, not a theorem: maintainers re-derive it independently and search for prior work before it is recorded as a finding.

## Workflow

1. Pick a track and read its `program.md` and `notes.md`, including the failed ideas.
2. Fork, branch, and run the loop. Keep your log in `notes.md`.
3. Open one pull request per result. Include: the path or certificate file, the checker output, the exact command to reproduce, and what you tried that did not work.
4. Do not merge your own result, and do not claim it is new. Say what you searched for.

## Rules

- Credit every idea, construction and person you build on. Keep licences.
- Disclose AI assistance and which model.
- Label claims: PROVEN, VERIFIED NUMERICALLY, CONJECTURE, or KNOWN (with citation).
- No hype. Exact statements, honest limits.

## Research prompt

Paste this into Claude Code, Codex or a similar agent from inside your fork:

```
Read README.md, CONTRIBUTING.md and tracks/README.md. Choose one open track. Read its program.md,
notes.md and results.tsv, and the open pull requests, so you do not repeat earlier attempts. Reproduce
the current best result with ./run.sh. Then run the loop in program.md: form a hypothesis, edit only the
search file, run, keep the change if the checker accepts it and the score improves, otherwise revert.
Never edit verify.py. Log every attempt, including failures, in notes.md. When you have a checker-accepted
improvement, characterise it (what structure, why it works, how it compares with known results), search
for prior work, and prepare one pull request with the certificate, the checker output, the reproduction
command and an honest novelty statement. Do not merge, and do not claim acceptance.
```
