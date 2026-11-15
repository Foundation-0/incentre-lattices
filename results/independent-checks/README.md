# Independent exact checks (November 2026 review)

Two programs written from the statements alone. They import nothing from this repository and share no code with its verifiers.
Both use exact arithmetic over Q(sqrt 2), built on `python-flint` (`fmpq`). This package is not in `requirements.txt`, and these
checks are not part of `make verify`. They have no certificate and are not listed in `certificates/` or `commitment/`.
Their files are bound by `SHA256SUMS` in this folder.

Evidence label: VERIFIED (a separately written program accepts the statement). Review: internal, by the review team that wrote the
programs. No outside human review. Priority is not claimed: the closest earlier source found is Borcea and Streinu 2021,
Section 7 (same bar framework, no spring constants, a sketched fold with no packing proof). The wording is "no earlier report found".

## 1. Elastic constants of results/rigidity-nd, n = 2 to 9

`elastic_indep.py` builds the lattice 2Z^n with dots at +-c e_i, c = 2 - sqrt 2, finds the bars by a distance search, puts a unit
spring with no pre-stress on every bar, relaxes the joints exactly, and compares with the closed forms. For every n from 2 to 9
it gives C44 = 0, C11 = 4(n-1)^2/(n^2 2^n), C12 = 4/(n^2 2^n), C11/C12 = (n-1)^2, nu_100 = 1/(n^2-n-1), and
(C11 - C12)/2 = 2(n-2)/(n 2^n). The last value is 0 at n = 2, 1/12 at n = 3 and positive for every n >= 3. So "every shear modulus is 0"
holds only at n = 2. The two assumptions are no stress in the bars at rest, and relaxed constants taken after the joints settle.
Before the joints settle (the printed `affine` values) the tensor has C12 = C44 exactly for every n.

It also prints the zero modes and self-stresses at q = 0 (n(n+1)/2 and n(n-1)/2), the rank of the relaxed tensor (n for n >= 3), and the
coordination (2n-1 bars at every joint, checked).

```
python3 -I results/independent-checks/elastic_indep.py --nmax 9      # about 4 minutes; output: out_elastic_n2_9.txt
python3 -I results/independent-checks/elastic_indep.py --nmax 4 --corrupt   # must exit 1; output: out_elastic_corrupt.txt
```

`--corrupt` multiplies the expected C12 by 1.01. The run must report MISMATCH and exit 1, also under `python3 -O`.

## 2. The cab to fcc fold of results/fold-continuous-fcc

`fold_indep.py` uses only the formula for L(s) and the motif from `PROOF.md`, Section 1. It proves in exact arithmetic that

* all 15 bars have squared length 2 for every s (polynomial identity);
* no two balls overlap for 6,516 pairs of lattice translation and motif difference (box [-3,3]^3), over the whole interval [0,1],
  by Bernstein coefficients with exact subdivision;
* translations outside the box are bounded by a different argument from the repository's: L(s) L(s)^T - I is positive
  semidefinite on [0,1] (all seven principal minors certified nonnegative), so the distance is at least |n| - 2;
* the end state has determinant 12, 12 contacts per ball and density pi/sqrt 18;
* det L(s) is strictly decreasing, so the density rises strictly.

```
python3 -I results/independent-checks/fold_indep.py                # 2 seconds; output: out_fold.txt
python3 -I results/independent-checks/fold_indep.py --corrupt 2>&1    # must exit 1; output: out_fold_corrupt.txt
```

`--corrupt` flips the sign of one bridge. Every bar keeps its length and balls overlap, so the run must fail, also under `python3 -O`.
Proof steps raise an error through `check()`; there is no bare `assert`.

## Check the files

```
python3 tools/sha256check.py results/independent-checks/SHA256SUMS
```

## What this does not show

The programs check the statements, not the written proofs of the general-n step. The elastic check stops at n = 9. The 6,516-pair
box and the PSD bound are a proof of non-overlap for the motion given in `PROOF.md`, not a check of the sampled paths of
`results/fold-close-packing`. The model is ideal central-force springs, so these constants say nothing about real crystals.
