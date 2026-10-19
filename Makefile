# Incentre lattice: replay every checked result.
#   make verify          run every verifier, rewrite each results/*/certificate.json, then check that
#                        certificates/ and commitment/ still match (about 6 minutes)
#   make verify-fast     the quick targets only (about 30 seconds)
#   make certificates    verify, then regenerate certificates/selected-results.json,
#                        certificates/SHA256SUMS and commitment/ (run this after a deliberate change)
#   make verify-<result> one result: rigidity, radius-ratio, cab6, fold2d, fold3d, fold-cp, fold-exact, fold-m1, fold-tests
# Setup: python -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt
# After `make verify`, `git diff --exit-code results certificates commitment` must be empty.
# All finite checks are labelled as such; see docs/audit.md for what each target does and does not show.

PY ?= python3

# One BLAS thread: multi-threaded reductions round differently from run to run and machine to machine.
export OPENBLAS_NUM_THREADS = 1
export OMP_NUM_THREADS = 1

FAST = verify-certificates-tests verify-rigidity verify-radius-ratio verify-cab6 verify-fold2d verify-fold-tests verify-fold-continuous verify-fold-continuous-hand
SLOW = verify-fold3d verify-fold-cp verify-fold-exact verify-fold-m1   # fold3d about 1 minute; fold-cp about 3 minutes (independent checker interpolates 10 sub-states per step); fold-exact about 30 seconds; fold-m1 about 10 seconds

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

check-hashes: verify-certificates-tests
	$(PY) tools/build_certificates.py --check
	$(PY) tools/sha256check.py certificates/SHA256SUMS
	$(PY) tools/sha256check.py commitment/SHA256SUMS

verify-certificates-tests:  # standard-library regression tests; in-memory filesystem only
	$(PY) -B tools/test_build_certificates.py

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

verify-fold-cp:  # slow: runs both checkers on the hcp path and on the fcc extension; the second checker and the path generator are hash-locked
	$(PY) tools/sha256check.py results/fold-close-packing/checkers.sha256 results/fold-close-packing
	$(PY) tools/sha256check.py autoresearch-fold/verify.sha256 autoresearch-fold
	$(PY) results/fold-close-packing/verify.py

verify-fold-exact:  # exact (integer) closed forms of the hcp and fcc end states
	$(PY) results/fold-end-exact/verify.py

verify-fold-m1:  # m = 1 paths to hcp and fcc (both checkers, exact end states, monotone fcc path)
	$(PY) results/fold-m1/verify.py

verify-fold-continuous:  # exact all-parameter certificate; standard library only
	$(PY) results/fold-continuous-fcc/verify.py

verify-fold-continuous-hand:  # sympy identity check of the written proof HAND-PROOF.md; replayed, not hashed into the certificate
	$(PY) results/fold-continuous-fcc/hand_proof_check.py

verify-fold-tests:  # the locked checker's own unit tests, and its hash lock
	$(PY) tools/sha256check.py autoresearch-fold/verify.sha256 autoresearch-fold
	cd autoresearch-fold && $(PY) test_verify.py

clean:
	find . -name __pycache__ -type d -prune -exec rm -rf {} +

.PHONY: lean
lean:  # optional: needs a Lean 4 toolchain on PATH (checked with v4.34.1); not part of `make verify`
	lean formal/Arithmetic.lean
