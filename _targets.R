library(targets)
library(tarchetypes)
tar_option_set(
    packages = c(
        "zen4R",
        "here",
        "quarto",
        "dagitty",
        "ggdag",
        "ggplot2"
    )
)
tar_source()

workshop_dir <- "data/river_workshop"

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
        clean_sea_data(sea_data_raw),
        format = "file"
    )
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
        quarto_args = c("--embed-resources"),
        quiet = FALSE
    )
)
