report_targets <- list(
  tar_quarto(
    manuscript,
    "output/manuscript.qmd",
    quarto_args = c("--embed-resources")
  )
)
