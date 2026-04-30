report_targets <- list(
  tar_quarto(
    index_main,
    "index.qmd",
    quarto_args = c("--embed-resources")
  ),
  tar_quarto(
    index_presentation,
    "presentations/index.qmd",
    quarto_args = c("--embed-resources")
  ),
  tar_quarto(report_A1,
    "analysis/A1-model-eel-abundance.qmd",
    quarto_args = c("--embed-resources")
  ),
  tar_quarto(
    report_A2,
    "analysis/A2-define-model.qmd",
    quarto_args = c("--embed-resources")
  ),
  tar_quarto(
    report_B2,
    "analysis/B2-get-diet-fishbase.qmd",
    quarto_args = c("--embed-resources")
  ),
  tar_quarto(
    report_B3,
    "analysis/B3-build-food-webs.qmd",
    quarto_args = c("--embed-resources")
  ),
  tar_quarto(
    report_C1,
    "analysis/C1-explore-coast-data.qmd",
    quarto_args = c("--embed-resources")
  ),
  tar_quarto(
    report_C2,
    "analysis/C2-get-fish-sizes.qmd",
    quarto_args = c("--embed-resources")
  ),
  tar_quarto(
    report_D1,
    "analysis/D1-spatial-distribution.qmd",
    quarto_args = c("--embed-resources")
  ),
  tar_quarto(
    report_D2,
    "analysis/D2-estuary-sea-continuum.qmd",
    quarto_args = c("--embed-resources")
  ),
  # Presentations.
  tar_quarto(
    presentation_wine,
    "presentations/2026-04-14_wine.qmd",
    quarto_args = c("--embed-resources")
  ),
  # Manuscript.
  tar_quarto(
    manuscript,
    "manuscript/manuscript.qmd",
    quarto_args = c("--embed-resources")
  )
)
