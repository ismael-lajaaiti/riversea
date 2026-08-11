
#' Build a spatial mesh over REPHY stations and fishing operations
#'
#' Uses a non-convex hull boundary so the mesh follows the coastline instead of
#' spanning the land between the Atlantic and Channel station clusters. Both
#' point sets are included so the mesh also covers fishing operations, not just
#' REPHY stations.
#'
#' @param rephy tibble with `longitude`/`latitude` columns, e.g. REPHY data.
#' @param foodwebs tibble with `longitude`/`latitude` columns, e.g. sea-survey
#'   fishing operations.
#'
#' @return an `inla.mesh` object.
#' @export
build_rephy_mesh <- function(rephy, foodwebs) {
  loc <- dplyr::bind_rows(
    rephy |> dplyr::distinct(longitude, latitude),
    foodwebs |> dplyr::distinct(longitude, latitude)
  ) |>
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
    control.compute = list(waic = TRUE),
    family = "gaussian"
  )

  idx_pred <- INLA::inla.stack.index(stack, "pred")$data
  result <- exp(fit$summary.fitted.values$mean[idx_pred])
  attr(result, "waic") <- fit$waic$waic
  result
}

#' Predict a REPHY parameter with an additive spatio-temporal SPDE model,
#' with salinity as a covariate
#'
#' Same as [predict_rephy_spde()], with an added fixed effect for salinity.
#' `train` and `newdata` must already be complete in `salinity`.
#'
#' @param train tibble with `value`, `longitude`, `latitude`, `year`,
#'   `month`, `salinity`.
#' @param newdata tibble with `longitude`, `latitude`, `year`, `month`,
#'   `salinity` to predict at - years should stay within `train`'s range.
#' @param mesh `inla.mesh`, e.g. from [build_rephy_mesh()].
#'
#' @return numeric vector, one prediction per row of `newdata`.
#' @export
predict_rephy_spde_salinity <- function(train, newdata, mesh) {
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
          month = d$month,
          salinity = d$salinity
        )
      )
    )
  }

  stack <- INLA::inla.stack(
    make_stack(train, "est", y = log(train$value)),
    make_stack(newdata, "pred")
  )

  fit <- INLA::inla(
    y ~ 0 + intercept + salinity +
      f(field, model = spde) +
      f(year, model = "rw1") +
      f(month, model = "rw1", cyclic = TRUE),
    data = INLA::inla.stack.data(stack),
    control.predictor = list(A = INLA::inla.stack.A(stack), compute = TRUE),
    control.compute = list(waic = TRUE),
    family = "gaussian"
  )

  idx_pred <- INLA::inla.stack.index(stack, "pred")$data
  result <- exp(fit$summary.fitted.values$mean[idx_pred])
  attr(result, "waic") <- fit$waic$waic
  result
}

#' Predict a REPHY parameter with an additive spatio-temporal SPDE model,
#' with salinity as a quadratic covariate
#'
#' Same as [predict_rephy_spde_salinity()], with an added `salinity^2` fixed
#' effect. `train` and `newdata` must already be complete in `salinity`.
#'
#' @param train tibble with `value`, `longitude`, `latitude`, `year`,
#'   `month`, `salinity`.
#' @param newdata tibble with `longitude`, `latitude`, `year`, `month`,
#'   `salinity` to predict at - years should stay within `train`'s range.
#' @param mesh `inla.mesh`, e.g. from [build_rephy_mesh()].
#'
#' @return numeric vector, one prediction per row of `newdata`.
#' @export
predict_rephy_spde_salinity_quadratic <- function(train, newdata, mesh) {
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
          month = d$month,
          salinity = d$salinity,
          salinity_sq = d$salinity^2
        )
      )
    )
  }

  stack <- INLA::inla.stack(
    make_stack(train, "est", y = log(train$value)),
    make_stack(newdata, "pred")
  )

  fit <- INLA::inla(
    y ~ 0 + intercept + salinity + salinity_sq +
      f(field, model = spde) +
      f(year, model = "rw1") +
      f(month, model = "rw1", cyclic = TRUE),
    data = INLA::inla.stack.data(stack),
    control.predictor = list(A = INLA::inla.stack.A(stack), compute = TRUE),
    control.compute = list(waic = TRUE),
    family = "gaussian"
  )

  idx_pred <- INLA::inla.stack.index(stack, "pred")$data
  result <- exp(fit$summary.fitted.values$mean[idx_pred])
  attr(result, "waic") <- fit$waic$waic
  result
}

#' Predict a REPHY parameter with an additive spatio-temporal SPDE model,
#' with salinity as a smooth (non-linear) covariate
#'
#' Same as [predict_rephy_spde_salinity()], but salinity enters as a random
#' walk over binned values (`f(salinity_group, model = "rw1")`) instead of a
#' linear fixed effect. Bins are computed jointly across `train` and
#' `newdata` (via `INLA::inla.group()`) so both share the same index.
#'
#' @param train tibble with `value`, `longitude`, `latitude`, `year`,
#'   `month`, `salinity`.
#' @param newdata tibble with `longitude`, `latitude`, `year`, `month`,
#'   `salinity` to predict at - years should stay within `train`'s range.
#' @param mesh `inla.mesh`, e.g. from [build_rephy_mesh()].
#'
#' @return numeric vector, one prediction per row of `newdata`, with a `sd`
#'   attribute (same length) giving the posterior SD of the linear
#'   predictor - on the log scale, since the model is fit on log(value).
#' @export
predict_rephy_spde_salinity_smooth <- function(train, newdata, mesh) {
  min_positive <- min(train$value[train$value > 0])
  train <- dplyr::mutate(
    train,
    value = dplyr::if_else(value == 0, min_positive / 2, value)
  )

  salinity_group <- INLA::inla.group(
    c(train$salinity, newdata$salinity),
    method = "quantile"
  )
  train$salinity_group <- salinity_group[seq_len(nrow(train))]
  newdata$salinity_group <-
    salinity_group[nrow(train) + seq_len(nrow(newdata))]

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
          month = d$month,
          salinity_group = d$salinity_group
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
      f(month, model = "rw1", cyclic = TRUE) +
      f(salinity_group, model = "rw1"),
    data = INLA::inla.stack.data(stack),
    control.predictor = list(A = INLA::inla.stack.A(stack), compute = TRUE),
    control.compute = list(waic = TRUE),
    family = "gaussian"
  )

  idx_pred <- INLA::inla.stack.index(stack, "pred")$data
  result <- exp(fit$summary.fitted.values$mean[idx_pred])
  attr(result, "waic") <- fit$waic$waic
  attr(result, "sd") <- fit$summary.fitted.values$sd[idx_pred]
  result
}

#' Fit the smooth-salinity REPHY model once, for later reuse
#'
#' Same model as [predict_rephy_spde_salinity_smooth()], fit on `train`
#' alone (`control.compute = list(config = TRUE)`, 100 salinity bins instead
#' of the default 25) so the fit can be reused by [sample_rephy_spde_smooth()]
#' at any later prediction locations without refitting.
#'
#' @param train tibble with `value`, `longitude`, `latitude`, `year`,
#'   `month`, `salinity`.
#' @param mesh `inla.mesh`, e.g. from [build_rephy_mesh()].
#'
#' @return list with `fit` (the `inla` object), `mesh`, `salinity_centers`
#'   (sorted bin centers) and `year_min` - everything
#'   [sample_rephy_spde_smooth()] needs.
#' @export
fit_rephy_spde_smooth <- function(train, mesh) {
  min_positive <- min(train$value[train$value > 0])
  train <- dplyr::mutate(
    train,
    value = dplyr::if_else(value == 0, min_positive / 2, value),
    salinity_group = INLA::inla.group(salinity, method = "quantile", n = 100)
  )

  spde <- INLA::inla.spde2.pcmatern(
    mesh,
    prior.range = c(1, 0.5),
    prior.sigma = c(1, 0.5)
  )
  field_idx <- INLA::inla.spde.make.index("field", n.spde = spde$n.spde)
  a <- INLA::inla.spde.make.A(mesh, loc = as.matrix(train[c("longitude", "latitude")]))
  stack <- INLA::inla.stack(
    tag = "est",
    data = list(y = log(train$value)),
    A = list(a, 1),
    effects = list(
      field_idx,
      list(
        intercept = 1,
        year = train$year - min(train$year) + 1,
        month = train$month,
        salinity_group = train$salinity_group
      )
    )
  )

  fit <- INLA::inla(
    y ~ 0 + intercept +
      f(field, model = spde) +
      f(year, model = "rw1") +
      f(month, model = "rw1", cyclic = TRUE) +
      f(salinity_group, model = "rw1"),
    data = INLA::inla.stack.data(stack),
    control.predictor = list(A = INLA::inla.stack.A(stack), compute = TRUE),
    control.compute = list(config = TRUE),
    family = "gaussian"
  )

  list(
    fit = fit,
    mesh = mesh,
    salinity_centers = sort(unique(train$salinity_group)),
    year_min = min(train$year)
  )
}

#' Linear-interpolation projection matrix onto fixed salinity bin centers
#'
#' Row `i` gets weight 1 at the nearest `centers` entry if `x[i]` falls
#' outside `centers`' range, or linearly-interpolated weights at the two
#' bracketing centers otherwise - the standard way to evaluate a
#' `model = "rw1"` term (a smooth curve known only at discrete points) at an
#' arbitrary `x`. `NA` in `x` gives an all-zero row: no salinity
#' contribution, same as [predict_rephy_spde_salinity_smooth()]'s handling
#' of missing salinity.
#'
#' @param x numeric vector of new salinity values.
#' @param centers sorted numeric vector of fitted bin centers.
#'
#' @return sparse `Matrix`, `length(x)` rows by `length(centers)` columns.
#' @noRd
.salinity_interp_A <- function(x, centers) {
  n <- length(x)
  k <- length(centers)
  has_val <- !is.na(x)
  x_clipped <- pmin(pmax(x, centers[1]), centers[k])
  j <- findInterval(x_clipped, centers, all.inside = TRUE)
  frac <- (x_clipped - centers[j]) / (centers[j + 1] - centers[j])

  Matrix::sparseMatrix(
    i = c(which(has_val), which(has_val)),
    j = c(j[has_val], j[has_val] + 1L),
    x = c((1 - frac)[has_val], frac[has_val]),
    dims = c(n, k)
  )
}

#' Draw posterior samples from an already-fitted smooth-salinity model
#'
#' Projects a [fit_rephy_spde_smooth()] fit onto new locations without
#' refitting: builds fresh A-matrices for the spatial field and for
#' salinity ([.salinity_interp_A()]) and applies them to posterior draws
#' via `INLA::inla.posterior.sample.eval()`.
#'
#' @param fit_obj output of [fit_rephy_spde_smooth()].
#' @param newdata tibble with `longitude`, `latitude`, `year`, `month`,
#'   `salinity` to predict at - years should stay within the original
#'   training range.
#' @param n_draws number of posterior draws.
#'
#' @return matrix, `nrow(newdata)` rows by `n_draws` columns.
#' @export
sample_rephy_spde_smooth <- function(fit_obj, newdata, n_draws) {
  a_spatial <- INLA::inla.spde.make.A(
    fit_obj$mesh, loc = as.matrix(newdata[c("longitude", "latitude")])
  )
  a_salinity <- .salinity_interp_A(newdata$salinity, fit_obj$salinity_centers)
  year_idx <- newdata$year - fit_obj$year_min + 1
  month_idx <- newdata$month

  samples <- INLA::inla.posterior.sample(n_draws, fit_obj$fit)
  # inla.posterior.sample.eval() replaces this function's environment with
  # one exposing each sample's components (field, year, month,
  # salinity_group, intercept) as plain names, parented to .GlobalEnv - it
  # can't see a_spatial/a_salinity/year_idx/month_idx via closure, so they
  # have to come in as arguments via `...` instead.
  proj <- function(a_spatial, a_salinity, year_idx, month_idx) {
    spatial <- as.vector(a_spatial %*% field)
    salinity_contrib <- as.vector(a_salinity %*% salinity_group)
    lp <- intercept[1] + spatial + year[year_idx] + month[month_idx] + salinity_contrib
    exp(lp)
  }
  INLA::inla.posterior.sample.eval(
    proj, samples,
    a_spatial = a_spatial, a_salinity = a_salinity,
    year_idx = year_idx, month_idx = month_idx
  )
}

#' Leave-station-out cross-validation for a REPHY parameter
#'
#' Splits stations (`site_id`) into `k` groups; for each, fits `predict_fn`
#' on the other stations and predicts at the held-out group.
#'
#' @param data tibble for one parameter, e.g. from [filter_year()].
#' @param mesh `inla.mesh`, e.g. from [build_rephy_mesh()].
#' @param predict_fn prediction function, e.g. [predict_rephy_spde()]
#'   (default) or [predict_rephy_spde_salinity()], called as
#'   `predict_fn(train, test, mesh)`.
#' @param k number of folds.
#'
#' @return `data` with added `pred` and `pred_baseline` columns.
#'   `pred_baseline` is the training-fold's geometric mean.
#' @export
cv_rephy_spde <- function(data, mesh, predict_fn = predict_rephy_spde, k = 10) {
  stations <- unique(data$site_id)
  station_fold <- stats::setNames(
    sample(rep_len(seq_len(k), length(stations))),
    as.character(stations)
  )
  data$fold <- station_fold[as.character(data$site_id)]

  purrr::map_dfr(seq_len(k), \(i) {
    train <- dplyr::filter(data, fold != i)
    test <- dplyr::filter(data, fold == i)
    test$pred <- predict_fn(train, test, mesh)
    test$pred_baseline <- exp(mean(log(train$value[train$value > 0])))
    test
  })
}

#' Model-selection criterion from a predict_rephy_spde*() result
#'
#' Reads the `waic` attribute attached by [predict_rephy_spde()] and its
#' salinity variants (INLA's Bayesian analogue of AIC) - meant for a single
#' full-data fit (`newdata` = `train`), not a cross-validation fold.
#'
#' @param pred numeric vector returned by [predict_rephy_spde()] or a
#'   salinity variant.
#'
#' @return one-row data.frame with `waic`.
#' @export
rephy_model_criteria <- function(pred) {
  data.frame(waic = attr(pred, "waic"))
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
