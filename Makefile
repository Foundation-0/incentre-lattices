# Incentre lattice: replay every checked result.
#   make verify          run every verifier, rewrite each results/*/certificate.json, then check that
#                        certificates/ and commitment/ still match (about 6 minutes)
#   make verify-fast     the quick targets only (about 30 seconds)
#   make certificates    verify, then regenerate certificates/selected-results.json,
#                        certificates/SHA256SUMS and commitment/ (run this after a deliberate change)
#   make verify-<result> one result: rigidity, radius-ratio, cab6, fold2d, fold3d, fold-cp, fold-exact, fold-tests
# Setup: python -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt
# After `make verify`, `git diff --exit-code results certificates commitment` must be empty.
# All finite checks are labelled as such; see docs/audit.md for what each target does and does not show.

PY ?= python3

FAST = verify-rigidity verify-radius-ratio verify-cab6 verify-fold2d verify-fold-tests
SLOW = verify-fold3d verify-fold-cp verify-fold-exact   # fold3d about 1 minute; fold-cp about 3 minutes (independent checker interpolates 10 sub-states per step); fold-exact about 30 seconds

.PHONY: verify verify-fast certificates check-hashes clean help $(FAST) $(SLOW)

help:
	@sed -n '2,10p' Makefile

verify: $(FAST) $(SLOW)
	$(PY) tools/build_certificates.py --check
	$(PY) tools/sha256check.py certificates/SHA256SUMS
	$(PY) tools/sha256check.py commitment/SHA256SUMS

verify-fast: $(FAST)

certificates: $(FAST) $(SLOW)
	$(PY) tools/build_certificates.py --write
	$(PY) tools/sha256check.py certificates/SHA256SUMS

check-hashes:
	$(PY) tools/build_certificates.py --check
	$(PY) tools/sha256check.py certificates/SHA256SUMS
	$(PY) tools/sha256check.py commitment/SHA256SUMS

verify-rigidity:
	$(PY) results/rigidity-nd/verify.py

verify-radius-ratio:
	$(PY) results/radius-ratio/verify.py

verify-cab6:
	$(PY) results/cab6-packing/verify.py

verify-fold2d:
	$(PY) results/fold-2d/verify.py

verify-fold3d:   # slow
	$(PY) results/fold-3d-cubic/verify.py

verify-fold-cp:  # slow: runs both checkers on the hcp path and on the fcc extension
	$(PY) results/fold-close-packing/verify.py

verify-fold-exact:  # exact (integer) closed forms of the hcp and fcc end states
	$(PY) results/fold-end-exact/verify.py

verify-fold-tests:  # the locked checker's own unit tests, and its hash lock
	$(PY) tools/sha256check.py autoresearch-fold/verify.sha256 autoresearch-fold
	cd autoresearch-fold && $(PY) test_verify.py

clean:
	find . -name __pycache__ -type d -prune -exec rm -rf {} +

.PHONY: lean
lean:  # optional: needs a Lean 4 toolchain on PATH (checked with v4.34.1); not part of `make verify`
	lean formal/Arithmetic.lean
