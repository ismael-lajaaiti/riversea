
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
#'
#' @param data tibble for one parameter, e.g. from [filter_year()].
#' @param mesh `inla.mesh`, e.g. from [build_rephy_mesh()].
#' @param k number of folds.
#'
#' @return `data` with added `pred` and `pred_baseline` columns.
#'   `pred_baseline` is the training-fold's geometric mean.
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
#'
#' @param cv tibble with `value` and `pred_col` columns
#' @param pred_col name of the column to score against `value` - `"pred"`
#'
#' @return named list: `rmse`, `mae`, `cor` (on the log scale).
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
