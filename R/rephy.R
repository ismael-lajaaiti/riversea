#' Download and unzip the REPHY dataset from SEANOE
#'
#' @param url character. URL of the REPHY zip archive on SEANOE.
#' @param dest_dir character. Destination folder where the archive will be
#'   downloaded and extracted (created if it doesn't exist).
#' @param overwrite logical. If FALSE (default) and the destination folder
#'   already contains files, the download/unzip step is skipped and the
#'   existing folder path is returned as-is.
#'
#' @return character. Path to dest_dir (invisibly usable as a targets
#'   dependency for downstream steps that read files from it).
download_rephy <- function(
  url = "https://www.seanoe.org/data/00361/47248/data/106989.zip",
  dest_dir = "data/rephy",
  overwrite = FALSE
) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  not_empty <- length(list.files(dest_dir)) > 0

  if (not_empty && !overwrite) {
    message("REPHY data already present in ", dest_dir, " - skipping download.")
    return(dest_dir)
  }

  zip_path <- file.path(tempdir(), "rephy_data.zip")

  message("Downloading REPHY archive from SEANOE...")
  utils::download.file(
    url      = url,
    destfile = zip_path,
    mode     = "wb",
    quiet    = FALSE
  )

  message("Unzipping into ", dest_dir, " ...")
  utils::unzip(zip_path, exdir = dest_dir)

  unlink(zip_path) # Remove zip.

  invisible(dest_dir) # Return dest_dir.
}

#' REPHY columns using "." as decimal mark
#'
#' The coordinate and depth columns use "." as decimal mark, unlike the rest
#' of the file (e.g. the measured value column) which follows the French
#' convention ("," as decimal mark) - see [read_rephy()].
#'
#' @return character vector of original REPHY column names
rephy_dot_decimal_cols <- function() {
  c(
    "Coordonnées passage : Coordonnées minx",
    "Coordonnées passage : Coordonnées maxx",
    "Coordonnées passage : Coordonnées miny",
    "Coordonnées passage : Coordonnées maxy",
    "Prélèvement : Immersion"
  )
}

#' Read a REPHY export file
#'
#' REPHY exports from SEANOE are `;`-delimited and ISO-8859-1 (Latin-1)
#' encoded, not comma-delimited UTF-8, hence the dedicated reader instead of
#' `readr::read_csv()`. The two files distributed on SEANOE (Med and
#' Manche-Atlantique) each also contain one spurious data row that repeats
#' the header, left over from concatenating two export batches; that row is
#' dropped here. Coordinate and depth columns are read as character and
#' parsed separately because they use "." as decimal mark while the rest of
#' the file uses "," (the French default) - reading everything with a
#' single locale silently corrupts them (e.g. "3.14" becomes 314).
#'
#' @param path character. Path to a REPHY csv file.
#'
#' @return tibble
#' @export
read_rephy <- function(path) {
  dot_decimal_col_types <- do.call(
    readr::cols,
    stats::setNames(
      rep(list(readr::col_character()), length(rephy_dot_decimal_cols())),
      rephy_dot_decimal_cols()
    )
  )
  data <- readr::read_csv2(
    path,
    locale = readr::locale(encoding = "ISO-8859-1"),
    col_types = dot_decimal_col_types,
    show_col_types = FALSE
  )
  data |>
    dplyr::filter(.data[[names(data)[1]]] != names(data)[1]) |>
    dplyr::mutate(
      dplyr::across(dplyr::all_of(rephy_dot_decimal_cols()), as.numeric)
    )
}

#' REPHY parameter codes for nitrogen, phosphorus and oxygen measurements
#'
#' Note that REPHY reports nitrate and nitrite as a single combined
#' measurement ("NO3+NO2"), unlike AMOBIO which reports them separately
#' (see [rename_amobio()]).
#'
#' @return character vector of REPHY parameter codes
#' @export
rephy_nutrient_codes <- function() {
  c("NH4", "NO3+NO2", "PO4", "OXYGENE")
}

#' Columns kept and their new names when extracting REPHY nutrient data
#'
#' Keeps only what is needed to identify the sampling operation (site,
#' sample, date, depth), its geographical coordinates, and the measured
#' parameter/value/unit.
#'
#' @return named character vector (new name = original REPHY column name)
#' @export
rename_rephy <- function() {
  c(
    site_id = "Lieu de surveillance : Identifiant",
    site_label = "Lieu de surveillance : Libellé",
    sample_id = "Prélèvement : Identifiant interne",
    date = "Passage : Date",
    level = "Prélèvement : Niveau",
    depth = "Prélèvement : Immersion",
    depth_unit = "Prélèvement : Symbole de l'unité d'immersion",
    longitude = "Coordonnées passage : Coordonnées minx",
    latitude = "Coordonnées passage : Coordonnées miny",
    parameter_code = "Résultat : Code paramètre",
    parameter_label = "Résultat : Libellé paramètre",
    value = "Résultat : Valeur de la mesure",
    unit = "Résultat : Symbole unité de mesure associé au quintuplet"
  )
}

#' Extract nutrient (nitrogen, phosphorus, oxygen) measurements from REPHY
#'
#' Filters rows to the nutrient parameters of interest and drops all columns
#' but those identifying the sampling operation, its coordinates, and the
#' measured value (see [rename_rephy()]).
#'
#' Dissolved oxygen is reported in two units ("ml.l-1" and "mg.l-1");
#' "ml.l-1" values are converted to "mg.l-1" (factor 1.42903, the standard
#' oceanographic conversion) so all oxygen values are directly comparable.
#'
#' @param data tibble read with [read_rephy()].
#'
#' @return tibble with one row per nutrient measurement.
#' @export
extract_nutrients_rephy <- function(data) {
  data |>
    dplyr::filter(
      .data[["Résultat : Code paramètre"]] %in% rephy_nutrient_codes()
    ) |>
    dplyr::filter_out(
      .data[["Résultat : Niveau de qualité"]] == "Douteux"
    ) |>
    dplyr::select(dplyr::all_of(rename_rephy())) |>
    dplyr::mutate(
      date = lubridate::dmy(date),
      value = dplyr::if_else(unit == "ml.l-1", value * 1.42903, value),
      unit = dplyr::if_else(unit == "ml.l-1", "mg.l-1", unit)
    )
}

#' Summarise coverage of REPHY nutrient data
#'
#' For each nutrient parameter, report the number of observations, the
#' number and share of missing values, the number of distinct sampling
#' sites, and the covered date range. Useful to compare spatial and temporal
#' coverage across parameters before further analysis.
#'
#' @param data tibble as returned by [extract_nutrients_rephy()].
#'
#' @return tibble with one row per parameter.
#' @export
summarise_rephy_coverage <- function(data) {
  data |>
    dplyr::summarise(
      n = dplyr::n(),
      n_sites = dplyr::n_distinct(site_id),
      year_min = lubridate::year(min(date, na.rm = TRUE)),
      year_max = lubridate::year(max(date, na.rm = TRUE)),
      .by = parameter_code
    ) |>
    dplyr::arrange(parameter_code)
}

#' Summarise yearly observation counts of REPHY nutrient data
#'
#' @param data tibble as returned by [extract_nutrients_rephy()].
#'
#' @return tibble with one row per (year, parameter).
#' @export
summarise_rephy_temporal_coverage <- function(data) {
  data |>
    dplyr::mutate(year = lubridate::year(date)) |>
    dplyr::reframe(n_obs = dplyr::n(), .by = c(year, parameter_code)) |>
    dplyr::arrange(year, parameter_code)
}

#' Plot yearly observation counts of REPHY nutrient data
#'
#' Heatmap of the number of observations for each nutrient parameter, across
#' years. Useful to spot parameters or periods with poor sampling coverage.
#'
#' @param data tibble as returned by [extract_nutrients_rephy()].
#'
#' @return ggplot
#' @export
plot_rephy_temporal_coverage <- function(data) {
  coverage <- summarise_rephy_temporal_coverage(data)
  ggplot2::ggplot(
    coverage,
    ggplot2::aes(x = year, y = parameter_code, fill = n_obs)
  ) +
    ggplot2::geom_tile(color = "#FAFAF8", linewidth = 0.6) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty()) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

#' Map REPHY sampling sites
#'
#' Plots distinct sampling sites (deduplicated on coordinates) over the
#' French coastline, to get a sense of the spatial resolution of the REPHY
#' network. When `facet` is `TRUE` (default), sites are faceted by nutrient
#' parameter to compare spatial coverage across parameters - some sites may
#' only report a subset of the nutrients of interest.
#'
#' @param data tibble as returned by [extract_nutrients_rephy()].
#' @param facet logical. Facet the map by `parameter_code`.
#'
#' @return ggplot
#' @export
plot_rephy_map <- function(data, facet = TRUE) {
  sites <- if (facet) {
    data |>
      dplyr::group_by(site_id, longitude, latitude, parameter_code) |>
      dplyr::summarise(value = mean(value)) |>
      ungroup() |>
      group_by(parameter_code) |>
      group_modify(~ .x |> mutate(value = log10(1 + range_zero_one(value)))) |>
      ungroup()
  } else {
    dplyr::distinct(data, site_id, longitude, latitude)
  }
  coastline <- rnaturalearth::ne_coastline(scale = "large", returnclass = "sf")
  p <- ggplot2::ggplot(sites, ggplot2::aes(x = longitude, y = latitude))
  if (facet) {
    p <- p +
      ggplot2::geom_point(mapping = ggplot2::aes(color = value), alpha = 0.8) +
      ggplot2::scale_color_viridis_c()
  } else {
    p <- p +
      ggplot2::geom_point(alpha = 0.8)
  }
  p <- p +
    ggplot2::geom_sf(data = coastline, inherit.aes = FALSE, linewidth = 0.3) +
    ggplot2::coord_sf(
      xlim = range(sites$longitude),
      ylim = range(sites$latitude)
    ) +
    ggplot2::labs(x = NULL, y = NULL)
  if (facet) {
    p <- p + ggplot2::facet_wrap(~parameter_code) +
      ggplot2::labs(color = "Scaled value [log10 scale]")
  }
  p
}

range_zero_one <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#' Match REPHY sampling sites to their DCE hydrographic district
#'
#' Joins each site to the coastal/transitional water body it falls within
#' (see [format_hydrographic_area()]); sites falling outside every polygon
#' are matched to the nearest district, provided it is within `dist_max` km.
#' A given `site_id` sometimes carries slightly different coordinates across
#' visits (rounding differences); only the first coordinate pair is used so
#' each site is matched exactly once.
#'
#' @param data tibble as returned by [extract_nutrients_rephy()].
#' @param hydro_zone sf tibble as returned by [format_hydrographic_area()].
#' @param dist_max maximum distance (km) for nearest-district fallback
#'   matching.
#'
#' @return the input data with an added `district` column; sites that could
#'   not be matched within `dist_max` are dropped.
#' @export
match_rephy_district <- function(data, hydro_zone, dist_max = 5) {
  sites <- data |>
    dplyr::distinct(site_id, longitude, latitude) |>
    dplyr::slice_head(n = 1, by = site_id) |>
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
    sf::st_join(hydro_zone["district"], join = sf::st_within)

  sites_within <- sites |> dplyr::filter(!is.na(district))
  sites_outside <- sites |> dplyr::filter(is.na(district))

  sites_nearest <- sites_outside |>
    dplyr::mutate(
      nearest_idx = sf::st_nearest_feature(geometry, hydro_zone),
      district = hydro_zone$district[nearest_idx],
      dist_km = sf::st_distance(
        geometry, hydro_zone[nearest_idx, ],
        by_element = TRUE
      ) |> as.numeric() / 1000
    ) |>
    dplyr::filter(dist_km < dist_max) |>
    dplyr::select(-nearest_idx, -dist_km)

  district_by_site <- dplyr::bind_rows(sites_within, sites_nearest) |>
    sf::st_drop_geometry() |>
    dplyr::select(site_id, district)

  dplyr::inner_join(data, district_by_site, by = "site_id")
}

#' Summarise yearly nutrient concentration by DCE district
#'
#' @param data tibble as returned by [match_rephy_district()].
#'
#' @return tibble with one row per (year, district, parameter), giving the
#'   mean concentration and the number of observations it is based on.
#' @export
summarise_rephy_yearly <- function(data) {
  data |>
    dplyr::mutate(year = lubridate::year(date)) |>
    dplyr::summarise(
      value = mean(value, na.rm = TRUE),
      n = dplyr::n(),
      .by = c(year, district, parameter_code)
    ) |>
    dplyr::arrange(parameter_code, district, year)
}

#' Plot yearly nutrient concentration trends by DCE district
#'
#' @param data tibble as returned by [summarise_rephy_yearly()].
#'
#' @return ggplot
#' @export
plot_rephy_temporal_trend <- function(data) {
  ggplot2::ggplot(
    data,
    ggplot2::aes(x = year, y = value, color = district)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::facet_wrap(~parameter_code, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Yearly mean concentration", color = NULL)
}

#' Summarise REPHY measurement method share by year and parameter
#'
#' For each nutrient parameter and year, the share of observations recorded
#' with each measurement method. Expressing counts as a share (rather than
#' raw counts) makes method changes visible independently of the sampling
#' effort trend (see [plot_rephy_temporal_coverage()]), which otherwise
#' dominates raw counts and masks method switches in low-effort years.
#'
#' @param raw_data tibble as returned by [read_rephy()] (must still have the
#'   method column, which [extract_nutrients_rephy()] drops).
#'
#' @return tibble with one row per (year, parameter_code, method), giving the
#'   observation count and its share of that (year, parameter_code) total.
#' @export
summarise_rephy_method_share <- function(raw_data) {
  raw_data |>
    dplyr::transmute(
      parameter_code = .data[["Résultat : Code paramètre"]],
      method = .data[["Résultat : Libellé méthode"]],
      year = .data[["Passage : Année"]]
    ) |>
    dplyr::filter(parameter_code %in% rephy_nutrient_codes()) |>
    dplyr::summarise(n = dplyr::n(), .by = c(parameter_code, method, year)) |>
    dplyr::mutate(share = n / sum(n), .by = c(parameter_code, year)) |>
    dplyr::arrange(parameter_code, year, dplyr::desc(share))
}

summarise_rephy_quality_share <- function(raw_data) {
  raw_data |>
    dplyr::transmute(
      parameter_code = .data[["Résultat : Code paramètre"]],
      quality = .data[["Résultat : Niveau de qualité"]],
      year = .data[["Passage : Année"]]
    ) |>
    dplyr::filter(parameter_code %in% rephy_nutrient_codes()) |>
    dplyr::summarise(n = dplyr::n(), .by = c(parameter_code, quality, year)) |>
    dplyr::mutate(share = n / sum(n), .by = c(parameter_code, year)) |>
    dplyr::arrange(parameter_code, year, dplyr::desc(share))
}

#' Plot REPHY measurement method share by year and parameter
#'
#' Stacked bar chart of the share of observations per method, for each
#' nutrient parameter and year. Useful to spot method changes over time
#' that could introduce artificial breaks in the concentration time series.
#'
#' @param data tibble as returned by [summarise_rephy_method_share()].
#'
#' @return ggplot
#' @export
plot_rephy_method_share <- function(data) {
  ggplot2::ggplot(
    data,
    ggplot2::aes(x = year, y = share, fill = method)
  ) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~parameter_code) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty()) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(x = NULL, y = "Share of observations", fill = "Method")
}

plot_rephy_quality_share <- function(data) {
  ggplot2::ggplot(
    data,
    ggplot2::aes(x = year, y = share, fill = quality)
  ) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~parameter_code) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty()) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(x = NULL, y = "Share of observations", fill = "Quality")
}

#' Filter rephy data keeping samples above depth_max
#'
#' @param rephy_data data.frame
#' @param depth_max double
#'
#' @return data.frame filtered
#' @export
filter_surface_rephy <- function(rephy_data, depth_max = 1) {
  rephy_data |>
    dplyr::filter(depth <= depth_max) |>
    dplyr::select(-level)
}

#' Filter rephy data keeping stations with non-small sampling effort
#'
#' Station with less samples than sampling_min are removed.
#'
#' @param rephy_data data tibble
#' @param sampling_min minimum sampling to be kept
#'
#' @return data.frame filtered
#' @export
filter_sampling_rephy <- function(rephy_data, sampling_min = 20) {
  station_n <- table(rephy_data$site_id)
  rephy_data |>
    dplyr::filter(station_n[as.character(site_id)] >= sampling_min)
}

#' Drop years with little from REPHY data
#'
#' @param data tibble.
#' @param year_min minimum year to be kept
#'
#' @return data tibble filtered.
#' @export
filter_year_rephy <- function(data, year_min = 2000) {
  data |>
    dplyr::mutate(
      year = lubridate::year(date),
      month = lubridate::month(date)
    ) |>
    dplyr::filter(year >= year_min)
}
