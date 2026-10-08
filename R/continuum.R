#' Variables describing the position of river operations along the
#' continuum.
#'
#' Site variables of Danet et al. (2021) from the RHT (altitude, slope,
#' distance from source, width, Strahler order) plus the hydrographic
#' distance to the river mouth. Depth is left out: in the RHT it carries the
#' same information as width. Variables with zeros are `log(x + 1)`
#' transformed. Only operations correctly matched to the RHT
#' (`rht_match`) and with all variables known are kept.
#'
#' @param rht_operation tibble, e.g. `check_rht_match()`'s output.
#' @param distance_to_mouth tibble with `operation_id`, `dist_mouth_m`.
#'
#' @return tibble, one row per operation, with `operation_id`, `longitude`,
#' `latitude`, `district` and the transformed variables (see
#' [continuum_variables()]).
#' @export
prepare_continuum_variables <- function(rht_operation, distance_to_mouth) {
  rht_operation |>
    dplyr::filter(rht_match) |>
    dplyr::inner_join(
      distance_to_mouth |> dplyr::select(operation_id, dist_mouth_m),
      by = "operation_id"
    ) |>
    dplyr::transmute(
      operation_id, longitude, latitude, district,
      log_altitude = log1p(altitude),
      log_slope = log1p(pente),
      log_dist_source = log1p(d_source),
      log_width = log(l),
      strahler,
      log_dist_mouth = log(dist_mouth_m / 1000)
    ) |>
    dplyr::filter(dplyr::if_all(
      dplyr::all_of(continuum_variables()), \(x) !is.na(x)
    ))
}

#' Names of the continuum variables used in the PCA.
#'
#' @return character vector.
#' @export
continuum_variables <- function() {
  c(
    "log_altitude", "log_slope", "log_dist_source", "log_width", "strahler",
    "log_dist_mouth"
  )
}

#' PCA of the continuum variables.
#'
#' Fitted on standardised variables, one row per sampling location
#' (repeated visits share the same site variables). PC1 is oriented so that
#' it increases with width (downstream), PC2 so that it increases with
#' distance to mouth.
#'
#' @param continuum tibble, e.g. `prepare_continuum_variables()`'s output.
#'
#' @return `prcomp` object.
#' @export
fit_continuum_pca <- function(continuum) {
  location <- continuum |>
    dplyr::distinct(longitude, latitude, .keep_all = TRUE) |>
    dplyr::select(dplyr::all_of(continuum_variables()))
  pca <- stats::prcomp(location, scale. = TRUE)

  pc_sign <- c(
    sign(pca$rotation["log_width", "PC1"]),
    sign(pca$rotation["log_dist_mouth", "PC2"])
  )
  pca$rotation[, 1:2] <- sweep(pca$rotation[, 1:2], 2, pc_sign, `*`)
  pca$x[, 1:2] <- sweep(pca$x[, 1:2], 2, pc_sign, `*`)
  pca
}

#' Scores of operations on the first two continuum PCA axes.
#'
#' @param continuum tibble, e.g. `prepare_continuum_variables()`'s output.
#' @param pca `prcomp` object, e.g. `fit_continuum_pca()`'s output.
#'
#' @return `continuum` with `PC1` and `PC2` added.
#' @export
score_continuum_pca <- function(continuum, pca) {
  scores <- stats::predict(pca, continuum)[, c("PC1", "PC2")]
  dplyr::bind_cols(continuum, tibble::as_tibble(scores))
}
