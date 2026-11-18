#' Dam visits with their water quality index.
#'
#' @param visit tibble, e.g. `prepare_dam_visits()`'s output.
#' @param quality tibble with `operation_id` and the predicted PC1 in `mean`,
#'   e.g. `predict_water_quality_visit()`'s output.
#' @return `visit` restricted to the visits with a predicted PC1, with `PC1`
#'   centred and scaled over these visits.
#' @export
prepare_quality_visit <- function(visit, quality) {
  visit |>
    dplyr::inner_join(dplyr::select(quality, operation_id, PC1 = mean), by = "operation_id") |>
    dplyr::mutate(PC1 = as.numeric(scale(PC1)))
}

#' Fit the presence model of a species with dams and water quality.
#'
#' Same model as [fit_presence_model()], with the water quality index `PC1`
#' available to `term` (e.g. "accumulation * PC1 + log_width").
#'
#' @param visit tibble, e.g. `prepare_quality_visit()`'s output.
#' @param presence tibble with `operation_id`, `species_valid`.
#' @param mesh INLA mesh, e.g. `build_presence_mesh()`'s output, or `NULL`
#'   for a model without spatial field.
#' @param species character, one species.
#' @param term character, right-hand side terms of the stressor and position
#'   effects.
#' @return list as [fit_presence_model()]'s output.
#' @export
fit_stressor_presence_model <- function(visit, presence, mesh, species, term) {
  data <- visit |>
    dplyr::transmute(
      present = as.numeric(
        operation_id %in% presence$operation_id[presence$species_valid == species]
      ),
      accumulation, PC1, log_width, log_alt, log_dist, log_slope,
      district = factor(district), protocol = factor(protocol), log_surface,
      site, year
    )
  sampling <- stats::model.matrix(~ district + protocol, data)[, -1, drop = FALSE]
  colnames(sampling) <- make.names(colnames(sampling))
  data <- dplyr::bind_cols(
    dplyr::select(data, -district, -protocol), tibble::as_tibble(sampling)
  )
  covariate <- dplyr::select(data, -present)
  if (is.null(mesh)) {
    A <- list(1)
    effects <- list(dplyr::mutate(covariate, Intercept = 1))
    field_term <- ""
  } else {
    spde <- INLA::inla.spde2.pcmatern(
      mesh, prior.range = c(50, 0.5), prior.sigma = c(1, 0.01)
    )
    A <- list(INLA::inla.spde.make.A(mesh, loc = visit_coordinates(visit)), 1)
    effects <- list(list(spatial = seq_len(spde$n.spde), Intercept = 1), covariate)
    field_term <- "+ f(spatial, model = spde)"
  }
  stack <- INLA::inla.stack(data = list(present = data$present), A = A, effects = effects)
  model <- INLA::inla(
    stats::as.formula(paste(
      "present ~ -1 + Intercept +", term, "+",
      paste(colnames(sampling), collapse = " + "), "+ log_surface +",
      "f(site, model = 'iid', hyper = list(prec = list(initial = 0))) +",
      "f(year, model = 'iid')", field_term
    )),
    family = "binomial", data = INLA::inla.stack.data(stack), num.threads = 4,
    quantiles = c(0.045, 0.5, 0.955),
    control.predictor = list(A = INLA::inla.stack.A(stack)),
    control.compute = list(config = FALSE, waic = TRUE)
  )
  out <- list(
    species = species,
    term = term,
    summary.fixed = model$summary.fixed,
    summary.hyperpar = model$summary.hyperpar,
    summary.random = model$summary.random,
    summary.linear.predictor = model$summary.linear.predictor[seq_len(nrow(visit)), ],
    waic = list(waic = model$waic$waic),
    mlik = model$mlik
  )
  rm(model)
  gc()
  out
}
