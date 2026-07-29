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

extract_river_mouth_amobio <- function(paths) {
  network_file <- paths[grepl("NETWORK_AMOBIO.Rdata", paths)]
  network <- get(base::load(network_file))
  network |>
    sfnetworks::activate("nodes") |>
    tibble::as_tibble() |>
    tibble::as_tibble() |>
    dplyr::filter(IS_RIVERMOUTH_LEVEL1 == "RIVERMOUTH_LEVEL1") |>
    dplyr::select(node_id, geometry)
}

match_mouth_district <- function(river_mouth, district, dist_max) {
  district_kept <- c("Loire", "Ardour-Garonne", "Seine", "Escaut-Somme")
  hydro_kept <- district |>
    filter(district %in% district_kept)

  river_mouth_sf <- river_mouth |>
    sf::st_as_sf() |>
    sf::st_transform(4326)

  joined <- river_mouth_sf |>
    sf::st_join(hydro_kept["district"], join = sf::st_within)
  matched <- joined |>
    dplyr::filter(!is.na(district)) |>
    dplyr::mutate(dist_km = 0)
  unmatched <- joined |> dplyr::filter(is.na(district))

  nearest_idx <- sf::st_nearest_feature(unmatched, hydro_kept)
  unmatched <- unmatched |>
    dplyr::mutate(
      district = hydro_kept$district[nearest_idx],
      dist_km = as.numeric(units::set_units(
        sf::st_distance(
          unmatched,
          hydro_kept[nearest_idx, ],
          by_element = TRUE
        ),
        "km"
      ))
    )

  river_mouth_district <- dplyr::bind_rows(matched, unmatched) |>
    dplyr::mutate(matched = dist_km <= dist_max)
}

#' Load the AMOBIO river network topology, stripped down to the minimum
#'
#' The raw network is ~1GB on disk / ~2.9GB in memory (1M nodes, 1.1M edges,
#' ~80 columns each - mostly links to unrelated monitoring networks, plus
#' full edge geometries). Keeps only what a shortest-path distance
#' computation needs: node geometry and the river-mouth flag, and edge
#' topology with the precomputed length (dropping edge geometry, which is
#' most of the memory).
#'
#' @param paths character vector of AMOBIO file paths.
#'
#' @return list(nodes, edges): `nodes` has `node_id`, `geometry`,
#'   `is_river_mouth`; `edges` has `from`, `to`, `length` (meters).
#' @export
extract_amobio_network <- function(paths) {
  network_file <- paths[grepl("NETWORK_AMOBIO.Rdata", paths)]
  network <- get(base::load(network_file))

  nodes <- network |>
    sfnetworks::activate("nodes") |>
    tibble::as_tibble() |>
    tibble::as_tibble() |>
    dplyr::transmute(
      node_id,
      geometry,
      is_river_mouth = IS_RIVERMOUTH_LEVEL1 %in% "RIVERMOUTH_LEVEL1"
    )

  edges <- network |>
    sfnetworks::activate("edges") |>
    tibble::as_tibble() |>
    tibble::as_tibble() |>
    dplyr::select(from, to, length = length_meters_num)

  list(nodes = nodes, edges = edges)
}

#' Restrict the AMOBIO network to within `max_dist` of the kept river mouths
#'
#' Grows the network outward (channel length, not straight-line) from every
#' matched river mouth at once, using a single Dijkstra run from a virtual
#' source connected to all of them with zero-weight edges - much cheaper than
#' one shortest-path search per mouth.
#'
#' @param network list(nodes, edges), see [extract_amobio_network()].
#' @param river_mouth tibble with `node_id`, `matched`, e.g. from
#'   [match_mouth_district()].
#' @param max_dist maximum channel-length distance to keep, in meters.
#'
#' @return list(nodes, edges) restricted to nodes within `max_dist`.
#' @export
restrict_amobio_network <- function(network, river_mouth, max_dist) {
  g <- igraph::graph_from_edgelist(
    as.matrix(network$edges[c("from", "to")]),
    directed = FALSE
  )
  igraph::E(g)$weight <- network$edges$length

  seeds <- river_mouth$node_id[river_mouth$matched]
  source <- igraph::vcount(g) + 1L
  g <- igraph::add_vertices(g, 1)
  g <- igraph::add_edges(g, as.vector(rbind(source, as.integer(seeds))), weight = 0)

  dist <- igraph::distances(g, v = source, weights = igraph::E(g)$weight)[1, ]
  kept_id <- setdiff(which(dist <= max_dist), source)

  list(
    nodes = dplyr::filter(network$nodes, node_id %in% kept_id),
    edges = dplyr::filter(network$edges, from %in% kept_id & to %in% kept_id)
  )
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

#' Join AMOBIO and river food web data.
#'
#' @param amobio dataframe.
#' @param river_foodweb river food web structure, as returned by
#' `measure_foodweb_structure()` (filtered to the river survey).
#' @param river_operation river operation metadata, as returned by
#' `get_river_operation()`, providing the `sandre_code`/`date` join keys
#' that `river_foodweb` itself doesn't carry.
#'
#' @return dataframe.
#' @export
join_amobio_aspe <- function(amobio, river_foodweb, river_operation) {
  aspe <- river_foodweb |>
    left_join(
      river_operation |> select(operation_id, sandre_code, date),
      by = join_by(operation_id)
    )
  # dplyr::inner_join() on an sf object currently errors (sf 1.1.1 / dplyr
  # 1.2.1 incompatibility), so geometry is carried through as plain x/y
  # columns across the join and reattached afterwards instead.
  amobio_df <- amobio |>
    mutate(
      geom_x = sf::st_coordinates(geometry)[, 1],
      geom_y = sf::st_coordinates(geometry)[, 2]
    ) |>
    sf::st_drop_geometry()
  inner_join(amobio_df, aspe, by = join_by(sandre_code, date)) |>
    sf::st_as_sf(coords = c("geom_x", "geom_y"), crs = sf::st_crs(amobio))
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
