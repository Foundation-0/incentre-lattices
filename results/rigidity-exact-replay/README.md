# Exact replay of the rigidity counts

This folder replays the counts in `results/rigidity-nd` with exact arithmetic, for n = 2 to 8. It uses no floating point and no tolerance. `rigidity-nd/rigidity_nd.py` does the same counts with numerical ranks (SVD, tolerance 1e-9). Here every rank is exact.

Label: VERIFIED (finite check). It covers the listed n and the listed wavevectors. It is not a proof for all n or all wavevectors. The proof stays in `results/rigidity-nd/PROOF.md`.

## What is replayed

A rank of the rigidity matrix gives the zero modes (motions that cost no bar stretch) and the self-stresses (bar tensions that balance at every joint). Both come from one rank.

| Statement | Where it is stated | Replayed for |
|---|---|---|
| The lone cross-polytope is infinitesimally rigid for n >= 3, with n(n+1)/2 zero modes and n(n-3)/2 self-stresses. For n = 2 the square has one internal flex. | `PROOF.md` (a), and "The case n = 2" | n = 2..8 |
| At generic wavevector: n(n-1)/2 zero modes and n(n-3)/2 self-stresses per cell | `PROOF.md` (b), and "The case n = 2" | n = 3..8, two phase families. For n = 2 the stated values are 2 and 0, and those are checked. |
| With k phases equal to 1: both counts rise by k | `PROOF.md` (b) | n = 3..8, k = 0..n. `PROOF.md` states no plane counts for n = 2. Those values are written to `exact_counts.json` but not compared. |
| At q = 0: n(n+1)/2 zero modes and n(n-1)/2 self-stresses | `PROOF.md` (b) | n = 2..8 |
| The explicit stress basis of `PROOF.md` (a) is a self-stress and has dimension n(n-3)/2 | `PROOF.md` (a) | n = 3..8 |

All of these agree with the stated numbers. Nothing differs.

## Two methods

`rigidity_exact.py` works over Q(i) with Python `fractions`. Every bar vector is a nonzero multiple of an integer direction, so the irrational factor c = 2 - sqrt2 divides out of each row. Scaling a row by a nonzero number does not change the rank or the number of self-stresses. The Bloch phase exp(2iq) is taken at exact points: 1, -1, and Pythagorean points such as (3+4i)/5.

`modp_check.py` keeps the true coordinates, with c = 2 - sqrt2, and works in the field GF(p) for the primes 1000000009, 998244353 and 2013265921. Each is 1 mod 8, so sqrt2 and i exist in GF(p). It also checks that the intra and cross bars have equal length. The two scripts share no code.

Reducing a matrix modulo p can lower its rank but never raise it. So a zero-mode count mod p is at least the true count. Three primes agreeing with the formulas is evidence, not proof.

## How it differs from the open draft

Draft pull request 102 adds `results/independent-checks/elastic_indep.py`. It computes the elastic constants over Q(sqrt2) and also prints the q = 0 zero modes and self-stresses by exact rank, for n = 2..9. That overlaps with the q = 0 row here and with nothing else.

This folder adds the cross-polytope rigidity, the generic-wavevector counts, the planes, the self-stress count n(n-3)/2, the explicit stress basis, and the check modulo a prime. It does not compute elastic constants and does not use python-flint. It touches no file that the draft touches.

## Run

```
python3 -I -O -B results/rigidity-exact-replay/rigidity_exact.py --check
python3 -I -O -B results/rigidity-exact-replay/modp_check.py
python3 -I -B -W error results/rigidity-exact-replay/test_rigidity_exact.py
python3 -I -O -B -W error results/rigidity-exact-replay/test_rigidity_exact.py
python3 -I -O -B results/rigidity-exact-replay/rigidity_exact.py --corrupt   # must exit 1
python3 -I -O -B results/rigidity-exact-replay/modp_check.py --corrupt       # must exit 1
python3 tools/sha256check.py results/rigidity-exact-replay/SHA256SUMS
```

`rigidity_exact.py` without `--check` rewrites `exact_counts.json`. The run takes about 25 seconds. Standard library only. This folder is outside `make verify` and has no certificate in `commitment/`. `SHA256SUMS` binds the files.

## Controls that must fail

Each control feeds a broken input to the same code that accepts the real framework, and the run stops if the broken input is accepted.

They include: delete one bar (the counts must change), cut one joint to degree n-1 (it must stop being rigid), move a joint off its axis (the bar lengths must stop being equal), sign-twist the stress basis (it must stop being a stress), delete a cross bar at generic phase, and claim the rigid count for the n = 2 square.

The test file adds known-answer values written out as literal numbers.

## What is not covered

* Elastic constants.
* All wavevectors. Phases are tried at chosen exact points only.
* Equal bar length from c = 2 - sqrt2 in exact arithmetic. `modp_check.py` checks it modulo three primes only.
* Identification of the zero modes as rotations. Only counts are replayed.
* n above 8. The list of Pythagorean points ends there.

## Known limits

A moved joint along its own axis keeps the counts (the theorem holds for any 0 < t < 1) and changes the bar lengths. The counts cannot reject that move. Only the length check can. `exact_counts.json` records it under `control_details`.
