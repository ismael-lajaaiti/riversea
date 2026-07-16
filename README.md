<p align="center">
  <img src="assets/logo.svg" alt="RIVERSEA" width="300">
</p>

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

## Reproducing the pipeline

The pipeline is containerized so it runs the same everywhere: R version,
system libraries, Quarto, and R package versions (via
[`renv`](https://rstudio.github.io/renv/)) are all pinned in the
`Dockerfile`.

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/).
2. Build the image once from the repository root (~15-20 min):
   ```sh
   docker build -t riversea .
   ```
3. Run a container, mounting the repo so data and rendered outputs are
   saved back to your machine:
   ```sh
   docker run --rm -it -v "$(pwd)":/project riversea
   ```
4. In the R console that opens, run the pipeline:
   ```r
   targets::tar_make()
   ```

Have R installed locally and prefer to skip Docker? Run `renv::restore()`
to install the pinned package versions, then step 4 directly (you'll also
need the system libraries and Quarto version listed in the `Dockerfile`).

## Collaborating

No Docker knowledge beyond the commands above is needed:

- Edit code (`R/`, `targets/`, `.qmd` files).
- Re-run `targets::tar_make()` to regenerate affected data, figures, and
  rendered documents.
- Commit code changes and any outputs meant for publishing (`docs/`,
  `output/`); `data/` is local-only and gitignored.
- Open a PR against `main`.

## License

Code (`R/`, `targets/`) is licensed under GPL-3 (see `LICENSE.md`). The
manuscript (`output/manuscript.qmd`) is licensed under
[CC BY-NC](https://creativecommons.org/licenses/by-nc/4.0/).
