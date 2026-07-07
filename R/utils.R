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

read_environment_data <- function(dir, year_start = 1993, year_end = 2023) {
  year_vals <- seq(year_start, year_end)
  purrr::map(year_vals, function(year) {
    read.csv(here(dir, paste0("Environment_QM_", year, ".csv"))) |>
      mutate(year = year)
  }) |>
    purrr::list_rbind()
}

download_hydrographic_files <- function(dir) {
  url_sea_polygones <- c(
    estuary = "https://services.sandre.eaufrance.fr/telechargement/geo/MDO/MasseDEau/vedl_2019/PolygMasseDEauTransition/PolygMasseDEauTransition_FXX-shp.zip",
    coast = "https://services.sandre.eaufrance.fr/telechargement/geo/MDO/MasseDEau/vrap_2022/MasseDEauCotiere/MasseDEauCotiere_FRA-shp.zip"
  )
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  purrr::walk2(
    url_sea_polygones, names(url_sea_polygones),
    \(url, name) {
      zip_file <- here::here(dir, paste0(name, ".zip"))
      download.file(url, destfile = zip_file)
      unzip(zipfile = zip_file, exdir = here::here(dir))
      file.remove(zip_file)
    }
  )
  dir
}

format_hydrographic_area <- function(dir) {
  # District classification from
  # https://www.donnees.statistiques.developpement-durable.gouv.fr/lesessentiels/essentiels/dce-district.html
  file_list <- list.files(path = dir, pattern = ".shp")
  purrr::map(
    here::here(dir, file_list),
    \(x) sf::read_sf(x) |> dplyr::select(CdEuMasseD)
  ) |>
    purrr::map(\(x) sf::st_transform(x, crs = 4326)) |>
    bind_rows() |>
    rename(zone_id = CdEuMasseD) |>
    sf::st_make_valid() |>
    mutate(district = case_when(
      str_detect(zone_id, "^FRA") ~ "Escaut-Somme",
      str_detect(zone_id, "^FRB") ~ "Meuse",
      str_detect(zone_id, "^FRC") ~ "Rhin",
      str_detect(zone_id, "^FRD") ~ "Rhone",
      str_detect(zone_id, "^FRE") ~ "Corse",
      str_detect(zone_id, "^FRF") ~ "Ardour-Garonne",
      str_detect(zone_id, "^FRG") ~ "Loire",
      str_detect(zone_id, "^FRH") ~ "Seine",
      TRUE ~ NA_character_
    ))
}

match_foodweb_district <- function(foodweb, hydro_zone, dist_max = 5) {
  fw_sf <- foodweb |>
    select(-c(
      foodweb,
      log_trophic_richness,
      log_species_richness,
      connectance,
      trophic_length
    )) |>
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
    st_join(hydro_zone, join = sf::st_within)

  # Food web within zone.
  fw_within <- fw_sf |>
    filter(!is.na(zone_id)) |>
    mutate(dist_to_district = 0)
  fw_outside <- fw_sf |>
    filter(is.na(zone_id))

  # Match outside point to nearest district.
  fw_nearest <- fw_outside |>
    mutate(
      nearest_idx = sf::st_nearest_feature(geometry, hydro_zone),
      district = hydro_zone$district[nearest_idx],
      dist_to_district = sf::st_distance(geometry, hydro_zone[nearest_idx, ],
        by_element = TRUE
      ) |> as.numeric() / 1000
    ) |>
    select(-nearest_idx) |>
    filter(dist_to_district < dist_max)

  fw_within |>
    rbind(fw_nearest) |>
    mutate(
      longitude = sf::st_coordinates(geometry)[, 1],
      latitude  = sf::st_coordinates(geometry)[, 2]
    ) |>
    sf::st_drop_geometry()
}

#' Convert a dataframe to sf object
#'
#' Assumes coordinates are given by the longitude and latitude columns.
#'
#' @param x data.frame
#'
#' @return sf object
#' @export
to_sf <- function(x) {
  x |> sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
}

#' Plot data points on the France map
#'
#' @param data data frame
#'
#' @return plot by ggplot
#' @export
plot_data_on_map <- function(data) {
  coastline <- rnaturalearth::ne_coastline(scale = "large", returnclass = "sf")
  ggplot(data, aes(x = longitude, y = latitude)) +
    geom_point(size = 1) +
    geom_sf(data = coastline, inherit.aes = FALSE, linewidth = 0.3) +
    coord_sf(
      xlim = range(data$longitude),
      ylim = range(data$latitude)
    ) +
    labs(x = NULL, y = NULL)
}

#' Theme for plots of the RIVERSEA project
#'
#' A clean, minimalist theme suitable for scientific publication: sober
#' typography, thin axis lines instead of a panel border, faint horizontal
#' gridlines only, and unobtrusive facet strips and legends.
#'
#' @param base_size base font size, in points.
#'
#' @return nothing, sets the theme globally via [ggplot2::theme_set()].
#' @export
nice_theme <- function(base_size = 11) {
  ink <- "grey20"
  muted <- "grey45"
  faint <- "grey85"
  ggplot2::theme_set(
    ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(
        # Text
        plot.title = ggplot2::element_text(
          face = "bold", size = ggplot2::rel(1.15), color = ink,
          margin = ggplot2::margin(b = 4)
        ),
        plot.subtitle = ggplot2::element_text(
          color = muted, size = ggplot2::rel(0.95),
          margin = ggplot2::margin(b = 10)
        ),
        plot.caption = ggplot2::element_text(
          color = muted, size = ggplot2::rel(0.75), hjust = 0
        ),
        axis.title = ggplot2::element_text(color = ink, size = ggplot2::rel(0.95)),
        axis.text = ggplot2::element_text(color = muted, size = ggplot2::rel(0.85)),
        # Axes and grid: thin axis lines instead of a panel border, only
        # faint horizontal gridlines to guide the eye across the plot.
        axis.line = ggplot2::element_line(color = ink, linewidth = 0.3),
        axis.ticks = ggplot2::element_line(color = ink, linewidth = 0.3),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(color = faint, linewidth = 0.3),
        # Facets
        strip.text = ggplot2::element_text(
          face = "bold", color = ink, size = ggplot2::rel(0.9)
        ),
        strip.background = ggplot2::element_blank(),
        # Legend
        legend.position = "bottom",
        legend.title = ggplot2::element_text(color = ink, size = ggplot2::rel(0.9)),
        legend.text = ggplot2::element_text(color = muted, size = ggplot2::rel(0.85)),
        legend.key = ggplot2::element_blank(),
        # Backgrounds and margins
        plot.margin = ggplot2::margin(t = 10, r = 15, b = 10, l = 15)
      )
  )
}
