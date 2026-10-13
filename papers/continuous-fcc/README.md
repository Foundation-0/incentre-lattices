# Focused research draft: an explicit cab-to-fcc motion

[Read the draft PDF](continuous-fcc.pdf) or [edit the TeX source](source.tex).
This note concentrates on the explicit continuous motion and its exact non-overlap certificate.
The broader paper remains in `papers/source/`; this is a separate, narrower research draft.

Status: written proof with an exact computational component, internally reviewed with two separately written
arithmetic implementations. Independent human review and historical priority remain unresolved. The linked
proof is in [the result folder](../../results/fold-continuous-fcc/PROOF.md); its complete replay is

```sh
make verify-fold-continuous PY=python3
```

Build from this directory with [Tectonic](https://tectonic-typesetting.github.io/):

```sh
tectonic source.tex
mv source.pdf continuous-fcc.pdf
```

No bibliography tool or figure-generation step is required. The PDF is a review artifact, not a new public release,
DOI deposit, or journal submission. The draft names the existing framework and endpoint precedents and preserves
the sampled-continuity labels of older numerical paths. See research assessment issue #26 and literature PR #27.

Suggested review order:

1. Confirm the model, row-vector convention, and exact bridge identity.
2. Check the symmetric-part bound reducing infinitely many images to the finite shell.
3. Inspect both exact arithmetic implementations and their exhaustive Bernstein checks.
4. Check the fcc coset argument and the determinant derivative.
5. Compare the exact motion with the closest prior literature before making a priority claim.
