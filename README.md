# RIVERSEA project

Analysis pipeline and manuscript for the RIVERSEA project, quantifying how
dams and nutrient loading affect fish food web structure along the
river-estuary-sea continuum in France. The project website (built with
[Quarto](https://quarto.org/)) is published at `docs/` from the sources
below.

## Repository structure

- `analysis/` — lab notebooks (Quarto `.qmd`). Exploratory, unpolished
  traces of the analysis as it happens, meant for internal/private use.
  Listed in `analysis/index.qmd`.
- `output/` — the manuscript, and in time synthetic reports meant to
  cleanly document the finalized pipeline for public sharing.
- `presentations/` — slides (Quarto) presented at meetings and conferences.
- `R/` — package functions used throughout the analysis (data
  download/preprocessing, food web construction, SEM fitting, plotting).
  Documented in `man/` (roxygen-generated).
- `targets/` and `_targets.R` — the [`targets`](https://books.ropensci.org/targets/)
  pipeline: `data_targets.R` (data download/processing), `plot_targets.R`
  (figures), `report_targets.R` (rendering notebooks/manuscript to
  `docs/`).
- `data/` — raw and intermediate data (gitignored; not versioned).
- `assets/` and `figures/` — images used in notebooks, the manuscript, and
  presentations.
- `index.qmd` / `_quarto.yml` — the project website: landing page and
  Quarto site configuration.

## Setup

Package dependencies are managed with [`renv`](https://rstudio.github.io/renv/).
Run `renv::restore()` to install them.
