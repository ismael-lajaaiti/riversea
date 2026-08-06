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
  utils::unzip(zip_file, exdir = dir)
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
  dag <- ggdag::dagify(
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
  ggdag::ggdag(dag, use_labels = "label") + ggdag::theme_dag()
}

#' Plot venn diagram of common species between surveys
#'
#' Styled to match `nice_theme()`: thin ink-grey outlines and labels, a
#' restrained blue fill, no legend (redundant with the in-plot counts).
#'
#' @param df data frame
#'
#' @return plot
#' @export
plot_venn_diagram <- function(df) {
  ink <- "grey20"
  # split() names each group after its own factor level - avoids mismatches
  # that group_split() + setNames(unique(df$survey)) can silently produce.
  distinct_df <- df |> distinct(species_valid, survey)
  sets <- split(distinct_df$species_valid, distinct_df$survey)

  ggVennDiagram::ggVennDiagram(
    sets,
    label = "count",
    label_geom = "text",
    label_color = ink,
    label_alpha = 1,
    label_size = 3.2,
    set_color = ink,
    set_size = 4,
    edge_size = 0.4
  ) +
    ggplot2::scale_fill_gradient(low = "white", high = "#a9c9e0") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.15)) +
    ggplot2::theme(legend.position = "none")
}

#' Get relative path.
#'
#' @param x path relative to the project root
#'
#' @return path relative to `output/`
#' @export
to_rel <- function(x) {
  file.path("..", x)
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
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
  france <- rnaturalearth::ne_countries(
    geounit = "france",
    type = "map_units",
    scale = "medium",
    returnclass = "sf"
  )
  rivers <- rnaturalearth::ne_download(
    scale = 10, type = "rivers_lake_centerlines", category = "physical",
    returnclass = "sf"
  )
  rivers_fr <- sf::st_intersection(rivers, france)
  d_station <- utils::read.csv(station_file, sep = ";") |>
    sf::st_as_sf(coords = c("lat", "long"), crs = 4326)
  d_sampling <- rbind(
    d_sf |> select(survey, geometry),
    d_station |> mutate(survey = "onema") |> select(survey, geometry)
  )
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = france) +
    ggplot2::geom_sf(data = rivers_fr, color = "grey") +
    ggplot2::geom_sf(data = d_sampling, ggplot2::aes(color = survey))
}

plot_network_groups <- function(diet_resource) {
  adj <- t(as.matrix(diet_resource |> select(-c(light, species, reference))))
  colnames(adj) <- diet_resource$species
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "directed")
  metanetwork::build_metanet(metaweb = g) |>
    metanetwork::compute_TL() |>
    metanetwork::ggmetanet()
}

read_environment_data <- function(dir, year_start = 1993, year_end = 2023) {
  year_vals <- seq(year_start, year_end)
  purrr::map(year_vals, function(year) {
    utils::read.csv(here::here(dir, paste0("Environment_QM_", year, ".csv"))) |>
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
    dir.create(dir, recursive = TRUE)
  }
  purrr::walk2(
    url_sea_polygones, names(url_sea_polygones),
    \(url, name) {
      zip_file <- here::here(dir, paste0(name, ".zip"))
      utils::download.file(url, destfile = zip_file)
      utils::unzip(zipfile = zip_file, exdir = here::here(dir))
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

download_hydrographic_basin <- function(dir) {
  url <- "https://services.sandre.eaufrance.fr/telechargement/geo/ETH/BDTopage/2025/BassinHydrographique/BassinHydrographique_FXX-shp.zip"
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  zip_file <- here::here(dir, "bassin.zip")
  utils::download.file(url, destfile = zip_file)
  utils::unzip(zipfile = zip_file, exdir = dir)
  file.remove(zip_file)
  dir
}

format_hydrographic_basin <- function(dir) {
  # LbBH is the basin's official name (BD Topage, SANDRE) - mapped to this
  # project's district scheme. Rhin-Meuse is a single combined polygon
  # upstream (unlike the FRB/FRC split used elsewhere) - not an issue since
  # neither Rhin nor Meuse is a kept district.
  file_list <- list.files(path = dir, pattern = ".shp")
  sf::read_sf(here::here(dir, file_list)) |>
    sf::st_transform(crs = 4326) |>
    sf::st_make_valid() |>
    mutate(district = case_when(
      LbBH == "Artois-Picardie" ~ "Escaut-Somme",
      LbBH == "Seine-Normandie" ~ "Seine",
      LbBH == "Loire-Bretagne" ~ "Loire",
      LbBH == "Adour-Garonne" ~ "Ardour-Garonne",
      LbBH == "Rhône-Méditerranée" ~ "Rhone",
      LbBH == "Corse" ~ "Corse",
      LbBH == "Rhin-Meuse" ~ "Rhin-Meuse",
      TRUE ~ NA_character_
    )) |>
    select(district, geometry)
}

#' Flag whether each operation sits inland or at sea, and restrict `river`
#' operations to the catchments of interest.
#'
#' Point-in-polygon against `hydrographic_basin` (land catchments), no
#' nearest-match fallback - a point just offshore stays "sea" rather than
#' being pulled into the nearest catchment. `district` is kept (not just
#' `inland`) so it can later be joined against `river_distance_to_mouth`.
#' `river` operations outside `district_kept` are dropped - this also
#' drops the handful of border-precision misses (no matched `district` at
#' all), since `NA` can't be in `district_kept` either.
#'
#' @param operation tibble with `longitude`, `latitude`, `survey`.
#' @param hydrographic_basin sf tibble with `district`, `geometry`, e.g.
#'   `format_hydrographic_basin()`'s output.
#' @param district_kept character vector of districts to keep `river`
#'   operations in.
#'
#' @return `operation` with `district` and `inland` added (TRUE if inside a
#' basin polygon), `river` rows restricted to `district_kept`.
#' @export
classify_operation_location <- function(operation, hydrographic_basin, district_kept) {
  basin_valid <- hydrographic_basin |> filter(!is.na(district))
  operation |>
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) |>
    sf::st_join(basin_valid["district"], join = sf::st_within) |>
    mutate(inland = !is.na(district)) |>
    sf::st_drop_geometry() |>
    filter(survey != "river" | district %in% district_kept)
}

#' District of the nearest polygon, within a distance threshold.
#'
#' @param points sf points.
#' @param hydro_zone sf polygons with `district`.
#' @param dist_max maximum distance to accept a match, in km. 0 means the
#'   point must be inside the polygon (`st_distance` to a polygon a point
#'   is inside is 0).
#'
#' @return character vector of `district`, one per row of `points`, `NA`
#' beyond `dist_max`.
nearest_district <- function(points, hydro_zone, dist_max) {
  nearest_idx <- sf::st_nearest_feature(points, hydro_zone)
  dist_km <- as.numeric(sf::st_distance(points, hydro_zone[nearest_idx, ], by_element = TRUE)) / 1000
  ifelse(dist_km <= dist_max, hydro_zone$district[nearest_idx], NA_character_)
}

match_sea_district <- function(foodweb, hydro_zone, dist_max = 5) {
  fw_sf <- foodweb |>
    select(-c(
      foodweb,
      log_trophic_richness,
      log_species_richness,
      connectance,
      trophic_length
    )) |>
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
    sf::st_join(hydro_zone, join = sf::st_within)

  # Food web within zone.
  fw_within <- fw_sf |>
    filter(!is.na(zone_id))
  fw_outside <- fw_sf |>
    filter(is.na(zone_id))

  # Match outside point to nearest district.
  fw_nearest <- fw_outside |>
    mutate(district = nearest_district(geometry, hydro_zone, dist_max)) |>
    filter(!is.na(district))

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
  ggplot2::ggplot(data, ggplot2::aes(x = longitude, y = latitude)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::geom_sf(data = coastline, inherit.aes = FALSE, linewidth = 0.3) +
    ggplot2::coord_sf(
      xlim = range(data$longitude),
      ylim = range(data$latitude)
    ) +
    ggplot2::labs(x = NULL, y = NULL)
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

download_sea_raw_data <- function(dest_dir = "data/sea/raw") {
  parent_dir <- dirname(dest_dir)
  dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)

  zip_path <- file.path(parent_dir, "sea.zip")

  zen4R::download_zenodo(
    doi = "10.5281/zenodo.21395668",
    path = parent_dir,
    files = "sea.zip"
  )

  utils::unzip(zip_path, exdir = parent_dir)
  file.remove(zip_path)

  dest_dir
}

#' Get coastline from rnaturalearth package.
#'
#'
#' @return sf object
#' @export
coastline <- function() {
  rnaturalearth::ne_coastline(scale = "large", returnclass = "sf")
}

#' Filter data such that year is in the specified time window
#'
#' @param data data tibble with `year` column.
#' @param year_min minimum year to be kept.
#' @param year_max maximal year to be kept.
#'
#' @return data tibble filtered.
#' @export
filter_year <- function(data, year_min, year_max) {
  data |>
    filter(year >= year_min, year <= year_max)
}

#' Convert date to year and month
#'
#' @param data tibble with `date` column.
#' @param keep logical, whether or not to keep the original `date` column.
#'
#' @return tibble
#' @export
from_date_to_year_month <- function(data, keep = FALSE) {
  data <- data |>
    mutate(
      year = lubridate::year(date),
      month = lubridate::month(date)
    )
  if (keep) data else data |> select(-date)
}
