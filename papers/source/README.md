# Preprint source

`main.tex`, `refs.bib` and `figures/` are the LaTeX source of the preprint.

| File | Version | Status |
|---|---|---|
| `papers/incentre-lattice-v3.pdf` | 3 | built from this source; not yet released or archived |
| `papers/incentre-lattice.pdf` | 2 | archived, doi:10.5281/zenodo.23340326; kept unchanged, built from the v2.0.0 source |
| `papers/incentre-lattice-v1.pdf` | 1 | archived, doi:10.5281/zenodo.23322080; kept unchanged |

## Build

The v3 PDF was built with Tectonic 0.17.0 (XeTeX engine, xdvipdfmx output):

```
cp -r papers/source /tmp/paper && cd /tmp/paper
tectonic -X compile main.tex        # writes main.pdf; runs BibTeX and reruns as needed
cp main.pdf <repo>/papers/incentre-lattice-v3.pdf
```

The v2 PDF was produced by an xdvipdfmx toolchain; building the v2.0.0 source with Tectonic 0.17.0 gives the same text.
The PDF bytes depend on the toolchain and the build date, so compare text (`pdftotext`), not hashes.

## Figures

From `papers/source/figures`, with the repository virtual environment (needs numpy and matplotlib):

```
../../../.venv/bin/python make_figures.py        # fig1_construction, fig2_fold (unused), fig3_modes, fig4_cab6
cd ../../.. && .venv/bin/python papers/source/figures/make_fold_figure.py   # fig_motion (Figure 3)
```

In v3 only `fig1_construction.pdf` was regenerated (the label "cross bar" became "bridge").
