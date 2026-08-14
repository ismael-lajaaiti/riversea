
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
#' Same as [predict_rephy_spde_salinity()], but salinity enters as a
#' P-spline: a cubic B-spline basis (18 basis functions, knots at `train`'s
#' salinity quantiles) with an `f(model = "rw2")` smoothness penalty on the
#' basis coefficients. The basis is fixed from `train` alone and evaluated
#' at `newdata` via [.spline_basis_A()], so no joint binning of `train` and
#' `newdata` is needed.
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

  bs_basis <- splines::bs(train$salinity, df = 18, degree = 3, intercept = TRUE)
  n_basis <- ncol(bs_basis)

  spde <- INLA::inla.spde2.pcmatern(
    mesh,
    prior.range = c(1, 0.5),
    prior.sigma = c(1, 0.5)
  )
  field_idx <- INLA::inla.spde.make.index("field", n.spde = spde$n.spde)

  make_stack <- function(d, tag, y = NA) {
    a_spatial <- INLA::inla.spde.make.A(mesh, loc = as.matrix(d[c("longitude", "latitude")]))
    a_spline <- .spline_basis_A(d$salinity, bs_basis)
    INLA::inla.stack(
      tag = tag,
      data = list(y = y),
      A = list(a_spatial, a_spline, 1),
      effects = list(
        field_idx,
        list(spline = seq_len(n_basis)),
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
      f(month, model = "rw1", cyclic = TRUE) +
      f(spline, model = "rw2"),
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
#' alone (`control.compute = list(config = TRUE)`) so the fit can be reused
#' by [sample_rephy_spde_smooth()] at any later prediction locations without
#' refitting.
#'
#' @param train tibble with `value`, `longitude`, `latitude`, `year`,
#'   `month`, `salinity`.
#' @param mesh `inla.mesh`, e.g. from [build_rephy_mesh()].
#'
#' @return list with `fit` (the `inla` object), `mesh`, `bs_basis` (the
#'   fitted `splines::bs()` object, fixing the basis for later prediction)
#'   and `year_min` - everything [sample_rephy_spde_smooth()] needs.
#' @export
fit_rephy_spde_smooth <- function(train, mesh) {
  min_positive <- min(train$value[train$value > 0])
  train <- dplyr::mutate(
    train,
    value = dplyr::if_else(value == 0, min_positive / 2, value)
  )

  bs_basis <- splines::bs(train$salinity, df = 18, degree = 3, intercept = TRUE)
  n_basis <- ncol(bs_basis)

  spde <- INLA::inla.spde2.pcmatern(
    mesh,
    prior.range = c(1, 0.5),
    prior.sigma = c(1, 0.5)
  )
  field_idx <- INLA::inla.spde.make.index("field", n.spde = spde$n.spde)
  a_spatial <- INLA::inla.spde.make.A(mesh, loc = as.matrix(train[c("longitude", "latitude")]))
  a_spline <- .spline_basis_A(train$salinity, bs_basis)

  stack <- INLA::inla.stack(
    tag = "est",
    data = list(y = log(train$value)),
    A = list(a_spatial, a_spline, 1),
    effects = list(
      field_idx,
      list(spline = seq_len(n_basis)),
      list(
        intercept = 1,
        year = train$year - min(train$year) + 1,
        month = train$month
      )
    )
  )

  fit <- INLA::inla(
    y ~ 0 + intercept +
      f(field, model = spde) +
      f(year, model = "rw1") +
      f(month, model = "rw1", cyclic = TRUE) +
      f(spline, model = "rw2"),
    data = INLA::inla.stack.data(stack),
    control.predictor = list(A = INLA::inla.stack.A(stack), compute = TRUE),
    control.compute = list(config = TRUE),
    family = "gaussian"
  )

  list(
    fit = fit,
    mesh = mesh,
    bs_basis = bs_basis,
    year_min = min(train$year)
  )
}

#' B-spline basis design matrix, with missing salinity zeroed out
#'
#' Evaluates `bs_basis` (a fitted `splines::bs()` object, fixing the knots,
#' degree and boundary) at `x` via `predict()`. `NA` entries in `x` give an
#' all-zero row: no salinity contribution, same convention as elsewhere for
#' missing salinity.
#'
#' @param x numeric vector of salinity values.
#' @param bs_basis a `splines::bs()` object, e.g. from
#'   [fit_rephy_spde_smooth()].
#'
#' @return matrix, `length(x)` rows by `ncol(bs_basis)` columns.
#' @noRd
.spline_basis_A <- function(x, bs_basis) {
  has_val <- !is.na(x)
  result <- matrix(0, nrow = length(x), ncol = ncol(bs_basis))
  result[has_val, ] <- predict(bs_basis, x[has_val])
  result
}

#' Draw posterior samples from an already-fitted smooth-salinity model
#'
#' Projects a [fit_rephy_spde_smooth()] fit onto new locations without
#' refitting: builds fresh A-matrices for the spatial field and for
#' salinity ([.spline_basis_A()]) and applies them to posterior draws via
#' `INLA::inla.posterior.sample.eval()`.
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
  a_spline <- .spline_basis_A(newdata$salinity, fit_obj$bs_basis)
  year_idx <- newdata$year - fit_obj$year_min + 1
  month_idx <- newdata$month

  samples <- INLA::inla.posterior.sample(n_draws, fit_obj$fit)
  proj <- function(a_spatial, a_spline, year_idx, month_idx) {
    spatial <- as.vector(a_spatial %*% field)
    salinity_contrib <- as.vector(a_spline %*% spline)
    lp <- intercept[1] + spatial + year[year_idx] + month[month_idx] + salinity_contrib
    exp(lp)
  }
  INLA::inla.posterior.sample.eval(
    proj, samples,
    a_spatial = a_spatial, a_spline = a_spline,
    year_idx = year_idx, month_idx = month_idx
  )
}

#' Predict REPHY nutrient concentrations, with uncertainty, at fishing
#' operations.
#'
#' Wraps [predict_rephy_spde_salinity_smooth()] (fit-and-predict-in-one-call,
#' same function/model as the rest of the notebook's point predictions) over
#' each parameter. Uses INLA's own posterior summary for `sd` (the linear
#' predictor's sd on the log scale, computed exactly), not an empirical sd
#' across posterior draws - the latter is noisy for a modest draw count and
#' can be inflated by occasional draws with extreme salinity-spline
#' extrapolation (see the `bs()` "beyond boundary knots" warning).
#'
#' @param rephy_data tibble with `parameter_code`, `value`, `longitude`,
#'   `latitude`, `year`, `month`, `salinity`, e.g. `rephy_salinity_complete`.
#' @param operation tibble with `operation_id`, `survey`, `longitude`,
#'   `latitude`, `year`, `month`, `salinity`, e.g. `foodweb_structure`
#'   restricted to `nurse`/`pomet`.
#' @param mesh `inla.mesh`, e.g. `rephy_mesh`.
#' @param parameters REPHY parameter codes to predict.
#'
#' @return tibble, one row per (operation, parameter): `operation_id`,
#' `survey`, `longitude`, `latitude`, `parameter_code`, `pred` (natural
#' scale), `sd` (log-scale sd of the linear predictor).
#' @export
predict_rephy_at_operations <- function(
  rephy_data, operation, mesh, parameters = c("NH4", "NO3+NO2", "PO4")
) {
  purrr::map_dfr(parameters, function(param) {
    train <- rephy_data |> dplyr::filter(parameter_code == param)
    pred <- predict_rephy_spde_salinity_smooth(train, operation, mesh)
    sd_val <- attr(pred, "sd")
    operation |>
      dplyr::select(operation_id, survey, longitude, latitude) |>
      dplyr::mutate(parameter_code = param, pred = as.numeric(pred), sd = sd_val)
  })
}

#' Map of interpolated REPHY nutrient concentrations at fishing operations.
#'
#' Three panels (NH4, NO3+NO2, PO4) of predicted concentration at fishing
#' operations, restricted to `district_kept`, each with its own legend
#' directly beneath it (concentration ranges differ too much between
#' parameters to share one).
#'
#' @param predictions tibble with `operation_id`, `survey`, `longitude`,
#'   `latitude`, `parameter_code`, `pred`, e.g.
#'   `predict_rephy_at_operations()`'s output.
#' @param operation_location tibble with `operation_id`, `district`, e.g.
#'   `classify_operation_location()`'s output.
#' @param district_kept character vector of districts to keep.
#' @param basin sf tibble with `district`, `geometry`, e.g.
#'   `format_basin()`'s output.
#'
#' @return patchwork object, 1 row x 3 columns.
#' @export
plot_rephy_predictions_map <- function(
  predictions, operation_location, district_kept, basin
) {
  nice_theme()

  pred <- predictions |>
    dplyr::left_join(
      operation_location |> dplyr::select(operation_id, district),
      by = "operation_id"
    ) |>
    dplyr::filter(district %in% district_kept)

  france <- rnaturalearth::ne_countries(
    geounit = "france", type = "map_units", scale = "medium",
    returnclass = "sf"
  )
  basin_kept <- basin |> dplyr::filter(district %in% district_kept)
  xlim <- c(-5.3, 8.3)
  ylim <- c(41.2, 51.2)

  params <- c("NH4", "NO3+NO2", "PO4")

  panels <- lapply(params, function(p) {
    pred_p <- pred |> dplyr::filter(parameter_code == p)
    ggplot2::ggplot() +
      ggplot2::geom_sf(data = france, fill = NA, color = "grey80", linewidth = 0.25) +
      ggplot2::geom_sf(
        data = basin_kept, fill = "grey97", color = "grey40", linewidth = 0.4
      ) +
      ggplot2::geom_point(
        data = pred_p, ggplot2::aes(longitude, latitude, color = pred),
        size = 0.6, alpha = 0.75
      ) +
      ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      ggplot2::scale_color_viridis_c(trans = "log10", name = "Concentration (µmol/L)") +
      ggplot2::labs(x = "Longitude", y = "Latitude", title = p) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 13, face = "bold"),
        axis.text = ggplot2::element_text(size = 9),
        axis.title = ggplot2::element_text(size = 10),
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 8, angle = 30, hjust = 1),
        legend.title = ggplot2::element_text(size = 10)
      )
  })

  patchwork::wrap_plots(panels, ncol = 3)
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
