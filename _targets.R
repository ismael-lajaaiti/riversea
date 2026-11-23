library(targets)
library(tarchetypes)
library(ggplot2)
set_theme(theme_minimal())

tar_option_set(
  controller = crew::crew_controller_local(workers = 3, seconds_idle = 60), # 3 parallel workers (INLA uses 4 threads each).
  storage = "worker",
  retrieval = "worker",
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
    "fmesher",
    "metanetwork",
    "cheddar",
    "vegan",
    "ggvegan",
    "lubridate",
    "mapview",
    "patchwork"
  )
)

tar_source()

lapply(list.files("targets", full.names = TRUE), source)

station_file <- "data/river/station_analysis.csv"

list(
  data_targets,
  dam_web_targets,
  presence_targets,
  foodweb_dam_targets,
  nutrient_targets,
  nutrient_driver_targets,
  foodweb_nutrient_targets,
  community_quality_targets,
  functional_diversity_targets,
  sea_community_targets,
  plot_targets,
  report_targets
)
