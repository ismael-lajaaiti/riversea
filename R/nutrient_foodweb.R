#' Food webs of the dam visits with their water quality index.
#'
#' @param web tibble, e.g. `measure_dam_webs()`'s output.
#' @param quality tibble with `operation_id` and the predicted PC1 of the
#'   water quality in `mean`, e.g. `predict_water_quality_visit()`'s output.
#' @param pca `prcomp` of the web metrics, e.g. `fit_structure_pca()`'s
#'   output.
#' @return `web` restricted to the visits with a predicted water quality,
#'   with `quality` centred and scaled over these visits, and the scores of
#'   the web on the first two axes of `pca` in `PC1` and `PC2`, scaled to
#'   unit variance.
#' @export
prepare_quality_foodweb <- function(web, quality, pca) {
  score <- stats::predict(pca, web)[, c("PC1", "PC2")]
  web |>
    dplyr::mutate(
      PC1 = score[, "PC1"] / stats::sd(score[, "PC1"], na.rm = TRUE),
      PC2 = score[, "PC2"] / stats::sd(score[, "PC2"], na.rm = TRUE)
    ) |>
    dplyr::inner_join(dplyr::select(quality, operation_id, quality = mean), by = "operation_id") |>
    dplyr::mutate(quality = as.numeric(scale(quality)))
}

#' Fit a model of a web structure metric on dams and water quality.
#'
#' Same model as [fit_structure_model()], with the water quality index
#' `quality` available to `term` (e.g. "accumulation + quality + log_width").
#'
#' @param web tibble, e.g. `prepare_quality_foodweb()`'s output.
#' @inheritParams fit_structure_model
#' @return list as [fit_structure_model()]'s output.
#' @export
fit_stressor_structure_model <- function(web, mesh, metric, family, term) {
  data <- web |>
    dplyr::transmute(
      response = .data[[metric]],
      accumulation, quality, log_width, log_alt, log_dist, log_slope,
      log_richness = log(fish_richness) - mean(log(fish_richness)),
      district = factor(district), protocol = factor(protocol), log_surface,
      site, year
    )
  sampling <- stats::model.matrix(~ district + protocol, data)[, -1, drop = FALSE]
  colnames(sampling) <- make.names(colnames(sampling))
  covariate <- dplyr::bind_cols(
    dplyr::select(data, -response, -district, -protocol), tibble::as_tibble(sampling)
  )
  spde <- INLA::inla.spde2.pcmatern(
    mesh, prior.range = c(50, 0.5), prior.sigma = c(1, 0.01)
  )
  stack <- INLA::inla.stack(
    data = list(response = data$response),
    A = list(INLA::inla.spde.make.A(mesh, loc = visit_coordinates(web)), 1),
    effects = list(list(spatial = seq_len(spde$n.spde), Intercept = 1), covariate)
  )
  model <- INLA::inla(
    stats::as.formula(paste(
      "response ~ -1 + Intercept +", term, "+",
      paste(colnames(sampling), collapse = " + "), "+ log_surface +",
      "f(site, model = 'iid', hyper = list(prec = list(initial = 0))) +",
      "f(year, model = 'iid') + f(spatial, model = spde)"
    )),
    family = family, data = INLA::inla.stack.data(stack), num.threads = 4,
    quantiles = c(0.045, 0.5, 0.955),
    control.predictor = list(A = INLA::inla.stack.A(stack)),
    control.compute = list(config = FALSE, waic = TRUE)
  )
  out <- list(
    metric = metric,
    term = term,
    summary.fixed = model$summary.fixed,
    summary.hyperpar = model$summary.hyperpar,
    waic = list(waic = model$waic$waic)
  )
  rm(model)
  gc()
  out
}

#' Change of a PCA axis when one species is removed.
#'
#' Same as [compute_removal_pc1()] for any axis, without the richness
#' correction.
#'
#' @inheritParams compute_removal_pc1
#' @param axis character, PCA axis, e.g. "PC2".
#' @return tibble with `metric` (the axis), `species`, `operation_id`,
#'   `sandre_code`, and the score of the web in `full` and without the
#'   species in `removed`, as input of [summarise_species_removal()].
#' @export
compute_removal_axis <- function(removal, pca, axis) {
  pca_metric <- rownames(pca$rotation)
  removal |>
    dplyr::filter(metric %in% pca_metric) |>
    tidyr::pivot_longer(c(full, removed), names_to = "web", values_to = "value") |>
    tidyr::pivot_wider(names_from = metric, values_from = value) |>
    dplyr::mutate(
      score = stats::predict(pca, dplyr::pick(dplyr::all_of(pca_metric)))[, axis] /
        pca$sdev[as.integer(sub("PC", "", axis))]
    ) |>
    dplyr::select(species, operation_id, sandre_code, web, score) |>
    tidyr::pivot_wider(names_from = web, values_from = score) |>
    dplyr::mutate(metric = axis)
}
