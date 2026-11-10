#' Station-year means of the NAIADES parameters.
#'
#' Mean of each parameter per station and year, kept when the parameter has at
#' least 4 measurements over at least 3 seasons (winter December-February,
#' spring, summer, autumn). Nutrients are log10-transformed; temperature and
#' oxygen are not. Only station-years with all parameters are returned.
#'
#' @param clean tibble, e.g. `clean_naiades_pc()`'s output.
#' @param code named character, Sandre codes of the parameters (names used as
#'   columns).
#' @return tibble with `code_station`, `year` and one column per parameter.
#' @export
compute_station_year_mean <- function(clean, code) {
  clean |>
    dplyr::mutate(
      parameter = factor(names(code)[match(code_parametre, code)], names(code)),
      year = as.integer(format(date, "%Y")),
      season = (as.integer(format(date, "%m")) %/% 3) %% 4
    ) |>
    dplyr::summarise(
      value = mean(value), n = dplyr::n(), seasons = dplyr::n_distinct(season),
      .by = c(parameter, code_parametre, code_station, year)
    ) |>
    dplyr::filter(n >= 4, seasons >= 3) |>
    dplyr::mutate(
      value = dplyr::if_else(code_parametre %in% c("1301", "1311", "1312"), value, log10(value))
    ) |>
    dplyr::select(code_station, year, parameter, value) |>
    dplyr::arrange(parameter) |>
    tidyr::pivot_wider(names_from = parameter, values_from = value) |>
    tidyr::drop_na()
}

#' PCA of the station-year means, PC1 increasing with water quality.
#'
#' @param station_year tibble, e.g. `compute_station_year_mean()`'s output.
#' @return `prcomp` object (standardised variables), PC1 signed so that its
#'   oxygen loading is positive.
#' @export
fit_water_quality_pca <- function(station_year) {
  pca <- stats::prcomp(dplyr::select(station_year, -code_station, -year), scale. = TRUE)
  if (pca$rotation["oxygen", "PC1"] < 0) {
    pca$rotation[, 1] <- -pca$rotation[, 1]
    pca$x[, 1] <- -pca$x[, 1]
  }
  pca
}

#' Water quality index (PC1) per station-year, with coordinates.
#'
#' @param station_year tibble, e.g. `compute_station_year_mean()`'s output.
#' @param pca `prcomp`, e.g. `fit_water_quality_pca()`'s output, fitted on
#'   `station_year`.
#' @param station sf points with `code_station`.
#' @return tibble with `code_station`, `year`, `PC1`, `X`, `Y` (Lambert-93,
#'   km).
#' @export
compute_water_quality <- function(station_year, pca, station) {
  coordinate <- station |>
    sf::st_transform(2154) |>
    sf::st_coordinates()
  station <- tibble::tibble(
    code_station = station$code_station,
    X = coordinate[, "X"] / 1000, Y = coordinate[, "Y"] / 1000
  )
  station_year |>
    dplyr::transmute(code_station, year, PC1 = pca$x[, "PC1"]) |>
    dplyr::inner_join(station, by = "code_station")
}

#' Mesh of the water quality model.
#'
#' Edges up to 40 km over the stations (nodes at least 15 km apart), up to
#' 150 km in the outer buffer.
#'
#' @param water_quality tibble with `X`, `Y` (km).
#' @return INLA mesh, coordinates in km (Lambert-93).
#' @export
build_water_quality_mesh <- function(water_quality) {
  INLA::inla.mesh.2d(
    loc = unique(as.matrix(water_quality[c("X", "Y")])), max.edge = c(40, 150), cutoff = 15
  )
}

#' Spatial-block folds of the cross-validation.
#'
#' Stations are grouped in square blocks of 50 km, and blocks are randomly
#' assigned to `k` folds, so that held-out stations have no close neighbour
#' in the training data.
#'
#' @param water_quality tibble with `X`, `Y` (km).
#' @param k integer, number of folds.
#' @return integer vector, fold of each row of `water_quality`.
#' @export
assign_block_fold <- function(water_quality, k) {
  block <- paste(water_quality$X %/% 50, water_quality$Y %/% 50)
  block_fold <- sample(rep_len(seq_len(k), dplyr::n_distinct(block)))
  block_fold[match(block, unique(block))]
}

#' Fit the space-time model of the water quality index.
#'
#' Gaussian model of PC1 with an intercept, a Matern spatial field (SPDE on
#' `mesh`) constant over time, an AR1 year effect shared by all stations and
#' an iid station effect. Rows
#' with `held_out` are not used in the fit and are predicted.
#'
#' @param water_quality tibble, e.g. `compute_water_quality()`'s output.
#' @param mesh INLA mesh, e.g. `build_water_quality_mesh()`'s output.
#' @param held_out logical, rows of `water_quality` left out of the fit.
#' @return list with `summary.fixed`, `summary.hyperpar` and `prediction`:
#'   `water_quality` with the posterior `mean` and `sd` of the linear
#'   predictor and the predictive `sd_observation` (linear predictor and
#'   observation noise), for all rows.
#' @export
fit_water_quality_model <- function(water_quality, mesh, held_out) {
  spde <- INLA::inla.spde2.pcmatern(mesh, prior.range = c(50, 0.5), prior.sigma = c(2, 0.01))
  stack <- INLA::inla.stack(
    data = list(PC1 = dplyr::if_else(held_out, NA, water_quality$PC1)),
    A = list(INLA::inla.spde.make.A(mesh, loc = as.matrix(water_quality[c("X", "Y")])), 1),
    effects = list(
      list(field = seq_len(spde$n.spde), Intercept = 1),
      list(station = water_quality$code_station, year = water_quality$year)
    )
  )
  model <- INLA::inla(
    PC1 ~ -1 + Intercept + f(field, model = spde) + f(year, model = "ar1") +
      f(station, model = "iid"),
    family = "gaussian", data = INLA::inla.stack.data(stack), num.threads = 4,
    control.predictor = list(A = INLA::inla.stack.A(stack)),
    control.inla = list(int.strategy = "eb")
  )
  predictor <- model$summary.linear.predictor[seq_len(nrow(water_quality)), ]
  noise <- 1 / model$summary.hyperpar["Precision for the Gaussian observations", "mean"]
  out <- list(
    summary.fixed = model$summary.fixed,
    summary.hyperpar = model$summary.hyperpar,
    prediction = dplyr::mutate(
      water_quality,
      held_out = held_out,
      mean = predictor$mean,
      sd = predictor$sd,
      sd_observation = sqrt(predictor$sd^2 + noise)
    )
  )
  rm(model)
  gc()
  out
}
