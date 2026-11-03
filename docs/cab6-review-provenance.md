# CaB6 attainment review provenance

The written all-pairs argument is in `results/cab6-packing/ATTAINMENT.md`. Its original submission reports an internal referee and a check over Q(sqrt2), but neither a named reviewer nor the check implementation is in the repository. These reports are historical provenance, not independently replayable evidence.

`make verify-cab6` runs symbolic identities and a floating-point finite pair search. It does not replay the missing exact check or establish that a human reviewed the proof. The written-proof label remains conditional on inspection of the argument, with independent human review pending.

A future exact implementation needs its source, input definition, translation bound, must-fail controls and reproduction command committed together. A future reviewer record needs an identity, exact reviewed revision and explicit verdict. No such evidence is inferred from a GitHub account name.
