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

#' Build a spatial mesh over REPHY station locations
#'
#' Uses a non-convex hull boundary so the mesh follows the coastline instead of
#' spanning the land between the Atlantic and Channel station clusters.
#'
#' @param data tibble with `longitude`/`latitude` columns.
#'
#' @return an `inla.mesh` object.
#' @export
build_rephy_mesh <- function(data) {
  loc <- data |>
    dplyr::distinct(longitude, latitude) |>
    as.matrix()
  boundary <- INLA::inla.nonconvex.hull(loc, convex = 0.3, resolution = 100)
  INLA::inla.mesh.2d(
    loc = loc,
    boundary = boundary,
    max.edge = c(0.3, 1),
    cutoff = 0.05
  )
}

#' Predict a REPHY parameter with an additive spatio-temporal SPDE model
#'
#' Fits `log(value) ~ spatial Matern field (SPDE) + year (RW1) + cyclic
#' month (RW1)`, additive only (no space-time interaction), and predicts at
#' `newdata`. Exact zeros (below detection limit) are set to half the
#' smallest positive value before fitting, instead of being dropped.
#'
#' @param train tibble with `value`, `longitude`, `latitude`, `year`, `month`.
#' @param newdata tibble with `longitude`, `latitude`, `year`, `month` to
#'   predict at - years should stay within `train`'s range.
#' @param mesh `inla.mesh`, e.g. from [build_rephy_mesh()].
#'
#' @return numeric vector, one prediction per row of `newdata`.
#' @export
predict_rephy_spde <- function(train, newdata, mesh) {
  min_positive <- min(train$value[train$value > 0])
  train <- dplyr::mutate(
    train, 
    value = dplyr::if_else(value == 0, min_positive / 2, value)
  )
  spde <- INLA::inla.spde2.pcmatern(
    mesh,
    prior.range = c(1, 0.5),
    prior.sigma = c(1, 0.5)
  )
  field_idx <- INLA::inla.spde.make.index("field", n.spde = spde$n.spde)

  make_stack <- function(d, tag, y = NA) {
    a <- INLA::inla.spde.make.A(mesh, loc = as.matrix(d[c("longitude", "latitude")]))
    INLA::inla.stack(
      tag = tag,
      data = list(y = y),
      A = list(a, 1),
      effects = list(
        field_idx,
        list(
          intercept = 1,
          year = d$year - min(train$year) + 1,
          month = d$month
        )
      )
    )
  }

  stack <- INLA::inla.stack(
    make_stack(train, "est", y = log(train$value)),
    make_stack(newdata, "pred")
  )

  fit <- INLA::inla(
    y ~ 0 + intercept +
      f(field, model = spde) +
      f(year, model = "rw1") +
      f(month, model = "rw1", cyclic = TRUE),
    data = INLA::inla.stack.data(stack),
    control.predictor = list(A = INLA::inla.stack.A(stack), compute = TRUE),
    family = "gaussian"
  )

  idx_pred <- INLA::inla.stack.index(stack, "pred")$data
  exp(fit$summary.fitted.values$mean[idx_pred])
}

#' Leave-station-out cross-validation for a REPHY parameter
#'
#' Splits stations (`site_id`) into `k` groups; for each, fits
#' [predict_rephy_spde()] on the other stations and predicts at the held-out
#' stations' observations. Tests spatial generalization to unseen locations,
#' rather than leave-one-observation-out, which would let a station's other
#' readings leak into its own training set.
#'
#' @param data tibble for one parameter, e.g. from [filter_year_rephy()].
#' @param mesh `inla.mesh`, e.g. from [build_rephy_mesh()].
#' @param k number of folds.
#'
#' @return `data` with added `pred` and `pred_baseline` columns.
#'   `pred_baseline` is the training-fold's geometric mean, i.e. a naive
#'   prediction using no spatial/temporal information at all - a reference
#'   point for how much `pred` actually improves on.
#' @export
cv_rephy_spde <- function(data, mesh, k = 10) {
  stations <- unique(data$site_id)
  station_fold <- stats::setNames(
    sample(rep_len(seq_len(k), length(stations))),
    as.character(stations)
  )
  data$fold <- station_fold[as.character(data$site_id)]

  purrr::map_dfr(seq_len(k), \(i) {
    train <- dplyr::filter(data, fold != i)
    test <- dplyr::filter(data, fold == i)
    test$pred <- predict_rephy_spde(train, test, mesh)
    test$pred_baseline <- exp(mean(log(train$value[train$value > 0])))
    test
  })
}

#' Cross-validation error metrics, on the log scale
#'
#' Log scale is the relevant one here: concentrations are right-skewed and
#' natural-scale error is dominated by rare extreme events (e.g. a
#' multi-year point-source spike at one station) that no spatially-smoothed
#' model can predict from its neighbors, and that would otherwise swamp the
#' metric regardless of model quality. Zero values are dropped (log
#' undefined).
#'
#' @param cv tibble with `value` and `pred_col` columns, e.g. from
#'   [cv_rephy_spde()].
#' @param pred_col name of the column to score against `value` - `"pred"`
#'   for the model, `"pred_baseline"` for the reference.
#'
#' @return named list: `rmse`, `mae`, `cor` (all on the log scale).
#' @export
rephy_cv_metrics <- function(cv, pred_col = "pred") {
  cv <- dplyr::filter(cv, value > 0)
  log_pred <- log(cv[[pred_col]])
  log_value <- log(cv$value)
  list(
    rmse = sqrt(mean((log_pred - log_value)^2)),
    mae = mean(abs(log_pred - log_value)),
    cor = stats::cor(log_pred, log_value)
  )
}
