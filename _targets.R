library(targets)
library(tarchetypes)
library(ggplot2)
set_theme(theme_minimal())

tar_option_set(
  packages = c(
    "zen4R",
    "here",
    "quarto",
    "dagitty",
    "ggdag",
    "ggplot2",
    "ggVennDiagram",
    "GGally",
    "tidyr",
    "tibble",
    "dplyr",
    "purrr",
    "stringr",
    "rfishbase",
    "truncdist",
    "readr",
    "rnaturalearth",
    "rnaturalearthdata",
    "sf",
    "foodwebbuilder",
    "igraph",
    "metanetwork",
    "cheddar",
    "vegan",
    "ggvegan",
    "lubridate",
    "mapview"
  )
)

tar_source()

lapply(list.files("targets", full.names = TRUE), source)

workshop_dir <- "data/river_workshop"
sea_data_raw <- "data/sea/raw"
station_file <- "data/river/station_analysis.csv"

list(
  data_targets,
  plot_targets,
  report_targets
)
