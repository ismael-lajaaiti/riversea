library(targets)
library(tarchetypes)
tar_option_set(
  packages = c(
    "zen4R",
    "here",
    "quarto",
    "dagitty",
    "ggdag",
    "ggplot2",
    "tidyr",
    "dplyr",
    "purrr",
    "stringr",
    "rfishbase",
    "truncdist"
  )
)
tar_source()

workshop_dir <- "data/river_workshop"
sea_data_raw <- "data/sea/raw"

list(
  # Workshop river data.
  tar_target(
    workshop_zip,
    download_workshop_data(workshop_dir),
    format = "file"
  ),
  tar_target(
    workshop_unzipped,
    unzip_workshop(workshop_zip, workshop_dir),
    format = "file"
  ),
  tar_target(
    plot_dag,
    create_plot_dag()
  ),
  # Sea survey data.
  tar_target(
    sea_data_tidy,
    preprocess_sea_data(sea_data_raw)
  ),
  # Infer missing size.
  tar_target(
    sea_data_imputed,
    infer_missing_size(sea_data_tidy)
  ),
  # Reports.
  tar_quarto(
    report_A1,
    "analysis/A1-model-eel-abundance.qmd",
    quarto_args = c("--embed-resources")
  ),
  tar_quarto(
    report_A2,
    "analysis/A2-define-model.qmd",
    quarto_args = c("--embed-resources")
  ),
  tar_quarto(
    report_C1,
    "analysis/C1-explore-coast-data.qmd",
    quarto_args = c("--embed-resources"),
    quiet = FALSE
  ),
  tar_quarto(
    report_C2,
    "analysis/C2-get-fish-sizes.qmd",
    quarto_args = c("--embed-resources"),
    quiet = FALSE
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
