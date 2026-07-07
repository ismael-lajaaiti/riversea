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

get_aspe_data <- function(aspe_file_foodweb, aspe_file_code) {

  aspe_file_foodweb <- here("data", "river", "output_size2webs.rda")
  aspe_file_code <- here("data", "river", "output_individual_fish.rda")

  foodweb <-
    get(base::load(aspe_file_foodweb))$tab_local_foodwebs_summary_metrics |>
    mutate(operation_id = as.integer(operation_id))

  code <- get(base::load(aspe_file_code))
  clean_code <- code$fishing_operation |>
    select(operation_id, site_id, date) |>
      left_join(
        code$station |> select(site_id, sandre_code),
        by = join_by(site_id)
      )

  foodweb |>
    left_join(clean_code, by = join_by(operation_id))
}

#' Join AMOBIO and ASPE datasets
#'
#' @param amobio dataframe
#' @param aspe dataframe
#'
#' @return dataframe
#' @export
join_amobio_aspe <- function(amobio, aspe) {
  inner_join(amobio, aspe, by = join_by(sandre_code, date))
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

#' Summarise completeness of environmental variables
#'
#' For each environmental variable, report the number of observations, the
#' number and share of missing values, and the range of observed (non-NA)
#' values. Useful to spot variables with poor coverage before modelling.
#'
#' @param data tibble containing the environmental variables as columns.
#' @param vars character vector of variable names to summarise. Defaults to
#' the AMOBIO environmental variables (see [rename_amobio()]).
#'
#' @return tibble with one row per variable.
#' @export
summarise_env_na <- function(
  data,
  vars = setdiff(names(rename_amobio()), c("date", "sandre_code"))
) {
  n_total <- nrow(data)
  purrr::map(vars, \(v) {
    x <- data[[v]]
    n_na <- sum(is.na(x))
    tibble::tibble(
      variable = v,
      n = n_total,
      n_na = n_na,
      pct_na = round(100 * n_na / n_total, 1),
      min = suppressWarnings(min(x, na.rm = TRUE)),
      median = suppressWarnings(stats::median(x, na.rm = TRUE)),
      max = suppressWarnings(max(x, na.rm = TRUE))
    )
  }) |>
    purrr::list_rbind()
}

#' Summarise temporal coverage of environmental variables
#'
#' For each year and environmental variable, count the number of
#' non-missing observations, out of the number of points sampled that year.
#' Useful to spot years or variables with poor sampling coverage.
#'
#' @param data tibble containing a `date` column and the environmental
#' variables as columns.
#' @param vars character vector of variable names to summarise. Defaults to
#' the AMOBIO environmental variables (see [rename_amobio()]).
#'
#' @return tibble with one row per (year, variable).
#' @export
summarise_env_temporal_coverage <- function(
  data,
  vars = setdiff(names(rename_amobio()), c("date", "sandre_code"))
) {
  data |>
    sf::st_drop_geometry() |>
    dplyr::mutate(year = lubridate::year(date)) |>
    dplyr::reframe(
      dplyr::across(dplyr::all_of(vars), \(x) sum(!is.na(x))),
      n = dplyr::n(),
      .by = year
    ) |>
    tidyr::pivot_longer(
      dplyr::all_of(vars),
      names_to = "variable",
      values_to = "n_obs"
    ) |>
    dplyr::mutate(
      pct_obs = 100 * n_obs / n,
      variable = factor(variable, levels = rev(vars))
    ) |>
    dplyr::arrange(year, variable)
}

#' Plot temporal coverage of environmental variables
#'
#' Heatmap of the share of non-missing observations for each environmental
#' variable, across years.
#'
#' @param data tibble containing a `date` column and the environmental
#' variables as columns.
#' @param vars character vector of variable names to summarise. Defaults to
#' the AMOBIO environmental variables (see [rename_amobio()]).
#'
#' @return ggplot
#' @export
plot_env_temporal_coverage <- function(
  data,
  vars = setdiff(names(rename_amobio()), c("date", "sandre_code"))
) {
  coverage <- summarise_env_temporal_coverage(data, vars)

  ggplot2::ggplot(
    coverage,
    ggplot2::aes(x = year, y = variable, fill = n_obs)
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

#' Summarise extreme values in environmental variables
#'
#' Flags extreme values using Tukey's fences (values below
#' `Q1 - 1.5 * IQR` or above `Q3 + 1.5 * IQR` are flagged), computed on the
#' log10 scale by default since nutrient concentrations are typically
#' right-skewed. Non-positive values are excluded from the log-scale fences
#' (reported separately) since they cannot be log-transformed.
#'
#' @param data tibble containing the environmental variables as columns.
#' @param vars character vector of variable names to check. Defaults to the
#' nutrient variables (phosphorus and nitrogen derivatives).
#' @param log_transform whether to compute the fences on the log10 scale.
#'
#' @return tibble with one row per variable.
#' @export
summarise_env_outliers <- function(
  data,
  vars = c(
    "total_phosphorus", "phosphate", "ammonium",
    "nitrite", "nitrate", "kjeldahl_nitrogen"
  ),
  log_transform = TRUE
) {
  purrr::map(vars, \(v) {
    x <- data[[v]]
    x <- x[!is.na(x)]
    x_pos <- if (log_transform) x[x > 0] else x
    x_t <- if (log_transform) log10(x_pos) else x_pos

    q <- stats::quantile(x_t, c(0.25, 0.75))
    iqr <- q[[2]] - q[[1]]
    lower <- q[[1]] - 1.5 * iqr
    upper <- q[[2]] + 1.5 * iqr
    is_outlier <- x_t < lower | x_t > upper

    tibble::tibble(
      variable = v,
      n = length(x),
      n_outlier = sum(is_outlier),
      pct_outlier = round(100 * sum(is_outlier) / length(x_t), 2),
      fence_low = if (log_transform) 10^lower else lower,
      fence_high = if (log_transform) 10^upper else upper,
      min = min(x),
      max = max(x)
    )
  }) |>
    purrr::list_rbind()
}

#' Boxplot of environmental variable distributions
#'
#' Boxplots on a log10 scale to visually spot extreme values: points beyond
#' the whiskers are flagged as outliers using Tukey's fences (1.5 * IQR),
#' the same convention as [summarise_env_outliers()] and
#' [ggplot2::geom_boxplot()]. Non-positive values are dropped since they
#' cannot be shown on a log scale. Optionally overlays a per-variable
#' reference threshold, e.g. a regulatory quality-class boundary.
#'
#' @param data tibble containing the environmental variables as columns.
#' @param vars character vector of variable names to plot. Defaults to the
#' nutrient variables (phosphorus and nitrogen derivatives).
#' @param thresholds optional tibble with a `variable` column (matching
#' `vars`) and a `threshold_col` column giving the reference value to
#' overlay for each variable, e.g. the SEQ-Eau grid (see
#' `data/river/seq_eau_nutrient_thresholds.csv`).
#' @param threshold_col name of the column in `thresholds` holding the
#' reference value.
#' @param threshold_label legend label for the threshold reference line.
#'
#' @return ggplot
#' @export
plot_env_boxplot <- function(
  data,
  vars = c(
    "total_phosphorus", "phosphate", "ammonium",
    "nitrite", "nitrate", "kjeldahl_nitrogen"
  ),
  thresholds = NULL,
  threshold_col = "mauvais_min",
  threshold_label = "SEQ-Eau \"mauvais\" threshold"
) {
  p <- data |>
    sf::st_drop_geometry() |>
    dplyr::select(dplyr::all_of(vars)) |>
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to = "variable", values_to = "value"
    ) |>
    dplyr::filter(!is.na(value), value > 0) |>
    ggplot2::ggplot(ggplot2::aes(x = variable, y = value)) +
    ggplot2::geom_boxplot(
      fill = "#cde2fb", color = "#104281",
      outlier.color = "#db8474", outlier.alpha = 0.5
    )

  if (!is.null(thresholds)) {
    thresholds_plot <- thresholds |>
      dplyr::filter(variable %in% vars) |>
      dplyr::rename(
        threshold = dplyr::all_of(threshold_col)
      )

    p <- p +
      ggplot2::geom_errorbar(
        data = thresholds_plot,
        ggplot2::aes(
          x = variable, ymin = threshold, ymax = threshold,
          color = threshold_label
        ),
        inherit.aes = FALSE,
        linewidth = 0.8, width = 0.7
      ) +
      ggplot2::scale_color_manual(
        name = NULL,
        values = stats::setNames("#9f66b3", threshold_label)
      )
  }

  p +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = NULL, y = "Concentration [log10 scale]") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}
