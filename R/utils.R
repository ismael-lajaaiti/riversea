#' Download the data from the spatiotemporal workshop from Zenodo.
#'
#' @param dir Directory where to store the zip file.
#' @return Path of the downloaded zip file.
#' @export
download_workshop_data <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  zen4R::download_zenodo("10.5281/zenodo.17962542", path = dir)
  file.path(dir, "miste_data.zip")
}

#' Unzip data from the spatiotemporal workshop.
#' Data is stored within zenodo folder.
#'
#' @param zip_file path of the zip file.
#' @param dir directory where to store the unzipped files.
#' @return Path of the unzipped directory.
#' @export
unzip_workshop <- function(zip_file, dir) {
  unzip(zip_file, exdir = dir)
  file.path(dir, "zenodo")
}

#' Create the plot of the DAG using ggplot.
#'
#' This represents the DAG we assume for our study.
#' That is, how environmental variable can shape food web structure.
#'
#' @return ggplot
#' @export
create_plot_dag <- function() {
  dag <- dagify(
    TL ~ S + Comp + Env,
    S ~ Env,
    Comp ~ Env,
    latent = "Comp",
    exposure = "Env",
    outcome = "TL",
    labels = c(
      TL = "Food web structure",
      S = "Richness",
      Comp = "Composition",
      Env = "Environment"
    )
  )
  ggdag(dag, use_labels = "label") + theme_dag()
}

#' Plot venn diagram of common species between surveys
#'
#' @param df data frame
#'
#' @return plot
#' @export
plot_venn_diagram <- function(df) {
  df |>
    distinct(species_valid, survey) |>
    group_split(survey) |>
    setNames(unique(df$survey)) |>
    lapply(\(x) x$species_valid) |>
    ggVennDiagram()
}

#' Get relative path.
#'
#' @param x the absolute path
#'
#' @return relative path
#' @export
to_rel <- function(x) {
  rel <- gsub(here::here(), "", x, fixed = TRUE)
  paste0("..", rel)
}

#' Plot sampling station sea and river survey.
#'
#' @param sea_data tidy data frame of sea data.
#' @param station_file file containing information about river station.
#'
#' @return plot.
#' @export
plot_sampling <- function(sea_data, station_file) {
  d_trait <- sea_data$trait
  d_sf <- d_trait |>
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
  france <- rnaturalearth::ne_countries(
    geounit = "france",
    type = "map_units",
    scale = "medium",
    returnclass = "sf"
  )
  rivers <- ne_download(
    scale = 10, type = "rivers_lake_centerlines", category = "physical",
    returnclass = "sf"
  )
  rivers_fr <- st_intersection(rivers, france)
  d_station <- read.csv(station_file, sep = ";") |>
    st_as_sf(coords = c("lat", "long"), crs = 4326)
  d_sampling <- rbind(
    d_sf |> select(survey, geometry),
    d_station |> mutate(survey = "onema") |> select(survey, geometry)
  )
  ggplot() +
    geom_sf(data = france) +
    geom_sf(data = rivers_fr, color = "grey") +
    geom_sf(data = d_sampling, aes(color = survey))
}

plot_network_groups <- function(diet_resource) {
  adj <- t(as.matrix(diet_resource |> select(-c(light, species, reference))))
  colnames(adj) <- diet_resource$species
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "directed")
  build_metanet(metaweb = g) |>
    compute_TL() |>
    ggmetanet()
}
