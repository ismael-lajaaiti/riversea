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

extract_metrics_amobio <- function(paths) {
  metric_file <- paths[grepl("METRICS_AMOBIO.Rdata", paths)]
  df_metric <- get(load(metric_file))
  df_metric |>
    dplyr::select(
      -matches(
        "B_FISH_|B_DIA_|B_INV_|NC_|H_|HM_|T_|P_|IS_|COMPARTMENT|META"
      ),
      -matches("^PC_(?!N_)", perl = TRUE)
    ) |>
    dplyr::rename_with(tolower) |>
    dplyr::rename(dplyr::any_of(rename_amobio())) |>
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

#' Combine the AMOBIO tibble in a single one.
#'
#' @param metrics tibble
#' @param nodes tibble
#'
#' @return sf object CRS Lambert-93
#' @export
combine_amobio_data <- function(metrics, nodes) {
  left_join(nodes, metrics, by = join_by(node_id)) |>
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

#' Dictionnary of renamed columns in the AMOBIO database.
#'
#' To keep track of the renaming and ensure consistency.
#'
#' @return named vector
#' @export
rename_amobio <- function() {
  c(
    "date" = "date_operation",
    "sandre_code" = "id",
    "elevation" = "elevation_ign25m"
  )
}
