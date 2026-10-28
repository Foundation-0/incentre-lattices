# Locked fold checker history

| Revision | SHA-256 | Effect |
|---|---|---|
| Initial search, `18dcb6e` | `5062e649aafe62ad2d08c8faf5255049865cf474c6d14840c24e3166d788f00b` | Original judge; `best/verify.json` is its stored output. |
| PR #11, `aea9dd6` | `05d2385dcd3ec5a4ec7fec5d23971c9df685c2ccedcea2edf43629ac813835cd` | End-state exceptions moved into reasons. Validity was assigned earlier, so this did not fail closed. |
| PR #50, `3d52aa8` | `e449913320b1902175382376691399bd91115aa48b23a49ae87c169a19de9e42` | Validity assigned after end-state analysis; exceptions now invalidate a path. Twenty-four checker tests pass. |

These pins identify different judge revisions. A later pin does not retroactively identify the judge used in the original search. The stored paths must be replayed against the selected revision before their current acceptance is asserted.
