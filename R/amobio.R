#' Download the AMOBIO data base
#'
#' From data.gouv at the link
#' https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/EQYVLP
#'
#' @param dir where the data is downloaded.
#'
#' @return nothing
#' @export
download_amobio_data <- function(dir, verbose = FALSE) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  Sys.setenv("DATAVERSE_SERVER" = "entrepot.recherche.data.gouv.fr")
  dataset <- dataverse::get_dataset("doi:10.57745/EQYVLP")
  files <- dataset$files
  paths <- purrr::map2_chr(files$id, files$label, \(id, filename) {
    if (verbose) {
      print(paste("Downloading", filename))
    }
    path <- file.path(dir, filename)
    writeBin(dataverse::get_file(id, progress = verbose), path)
    path
  })
  paths
}

extract_metrics_amobio <- function(paths, average_time) {
  metric_file <- paths[grepl("METRICS_AMOBIO.Rdata", paths)]
  df_metric <- get(base::load(metric_file))
  pc_x_cols <- names(df_metric) |> stringr::str_subset("^PC_X")
  df_metric <- df_metric |>
    dplyr::select(-dplyr::any_of(pc_x_cols)) |>
    dplyr::rename_with(tolower) |>
    dplyr::select(
      node_id, id, date_operation,
      dplyr::any_of(amobio_selected_vars(average_time))
    ) |>
    dplyr::rename(dplyr::any_of(rename_amobio(average_time))) |>
    dplyr::mutate(
      node_id = as.integer(as.character(node_id)),
      sandre_code = stringr::str_pad(as.character(sandre_code), 8, pad = "0")
    )
}

extract_nodes_amobio <- function(paths) {
  network_file <- paths[grepl("NETWORK_AMOBIO.Rdata", paths)]
  network <- get(load(network_file))
  network |>
    sfnetworks::activate("nodes") |>
    tibble::as_tibble() |>
    tibble::as_tibble() |>
    dplyr::filter(IS_FISH == "FISH") |>
    dplyr::select(node_id, geometry, ELEVATION_IGN25M) |>
    dplyr::rename_with(tolower) |>
    dplyr::rename(dplyr::any_of(rename_amobio()))
}

#' Deduplicate AMOBIO metrics on (sandre_code, date).
#'
#' Some nodes have several candidate physicochemical stations
#' (`METADATA_PC_CdStation`) for the same date, which is dropped upstream.
#' This leaves rows that share `(sandre_code, date)` but disagree on the
#' `pc_n_*` nutrient columns, so those are averaged across candidates while
#' every other column (identical across candidates) is kept as is.
#'
#' @param metrics tibble
#'
#' @return tibble with one row per (sandre_code, date)
dedup_amobio_metrics <- function(metrics) {
  pc_cols <- names(metrics) |> stringr::str_subset("^pc_n_")
  other_cols <- setdiff(names(metrics), c(pc_cols, "sandre_code", "date"))

  metrics |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(pc_cols),
        \(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
      ),
      dplyr::across(dplyr::all_of(other_cols), dplyr::first),
      .by = c(sandre_code, date)
    )
}

#' Combine the AMOBIO tibble in a single one.
#'
#' Also select variables of interest.
#'
#' @param metrics tibble
#' @param nodes tibble
#'
#' @return sf object CRS Lambert-93
#' @export
combine_amobio_data <- function(metrics, nodes) {
  inner_join(nodes, metrics, by = join_by(node_id)) |>
    sf::st_as_sf()
}

#' Utility function to visualise AMOBIO data.
#'
#' @param amobio_data dataframe
#' @param metric string specifying the metric name to plot
#' @param title string giving the title of the plot
#' @param subtitle string giving the subtitle of the plot
#' @param legend string giving the legend (color caption) of the plot
#' @param log_transform whether or not to log transform the metric
#' @param year string to select point of a given year
#'
#' @return plot
#' @export
plot_amobio_metric <- function(
  amobio_data,
  metric,
  title = NULL,
  subtitle = NULL,
  legend = NULL,
  f_transform = \(x) x,
  year = NA
) {
  data <- amobio_data |>
    select(all_of(metric), date) |>
    filter(!is.na(.data[[metric]])) |>
    distinct() |>
    mutate(to_plot = f_transform(.data[[metric]]))

  if (!is.na(year)) {
    data <- data |>
      filter(lubridate::year(date) == year)
  }

  coast <- rnaturalearth::ne_coastline(scale = "large", returnclass = "sf") |>
    sf::st_transform(sf::st_crs(data))
  boundaries <- sf::st_bbox(sf::st_as_sf(data))

  ggplot2::ggplot() +
    ggplot2::geom_sf(data = coast) +
    ggplot2::geom_sf(
      data = data,
      ggplot2::aes(color = to_plot),
      size = 0.5
    ) +
    ggplot2::coord_sf(
      xlim = boundaries[c("xmin", "xmax")],
      ylim = boundaries[c("ymin", "ymax")]
    ) +
    hrbrthemes::scale_color_flexoki_continuous() +
    ggplot2::labs(
      x = "Latitude",
      y = "Longitude",
      color = legend,
      title = title,
      subtitle = subtitle
    )
}

#' Environmental variables from the AMOBIO databsed selected for the analysis
#'
#' @param average_time string indicating the duration of the average for
#' chemical variabels. Possible values are "1y" (one year), "3m" (three months),
#' "5y" (five years) and "all" (all available points).
#'
#' @return vector of string
#' @export
amobio_selected_vars <- function(average_time) {
  c(
    "elevation",
    "poe_ddam_l3_dam",
    "poe_ddam_l2_dam",
    "poe_bh5_l6",
    "poe_bh5_l5",
    paste("pc_n_ptotal_s1", average_time, sep = "_"),
    paste("pc_n_orthophosp_s1", average_time, sep = "_"),
    paste("pc_n_no3_s1", average_time, sep = "_"),
    paste("pc_n_no2_s1", average_time, sep = "_"),
    paste("pc_n_nh4_s1", average_time, sep = "_"),
    paste("pc_n_nkj_s1", average_time, sep = "_")
  )
}

#' Rename selected environmental variables of the AMOBIO database.
#'
#' @return vector of strings.
#' @export
rename_amobio <- function(average_time = "1y") {
  c(
    "date" = "date_operation",
    "sandre_code" = "id",
    "elevation" = "elevation_ign25m",
    "dam_distance_upstream" = "poe_ddam_l3_dam",
    "dam_distance_downstream" = "poe_ddam_l2_dam",
    "barrier_downstream" = "poe_bh5_l5",
    "barrier_upstream" = "poe_bh5_l6",
    "total_phosphorus" = paste("pc_n_ptotal_s1", average_time, sep = "_"),
    "phosphate" = paste("pc_n_orthophosp_s1", average_time, sep = "_"),
    "ammonium" = paste("pc_n_nh4_s1", average_time, sep = "_"),
    "nitrite" = paste("pc_n_no2_s1", average_time, sep = "_"),
    "nitrate" = paste("pc_n_no3_s1", average_time, sep = "_"),
    "kjeldahl_nitrogen" = paste("pc_n_nkj_s1", average_time, sep = "_")
  )
}
