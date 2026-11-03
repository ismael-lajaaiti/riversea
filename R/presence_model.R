#' Coordinates of the visits for the spatial field (Lambert-93, km).
#'
#' @param visit tibble with `longitude`, `latitude` (WGS84).
#' @return matrix with columns `X`, `Y`, one row per visit.
#' @export
visit_coordinates <- function(visit) {
  visit |>
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
    sf::st_transform(2154) |>
    sf::st_coordinates() / 1000
}

#' Mesh of the spatial field of the presence models.
#'
#' Small triangles (edges up to 40 km) over the sites, larger ones (up to
#' 150 km) in the outer buffer.
#'
#' @param visit tibble, e.g. `prepare_dam_visits()`'s output.
#' @return INLA mesh, coordinates in km (Lambert-93).
#' @export
build_presence_mesh <- function(visit) {
  INLA::inla.mesh.2d(
    loc = unique(visit_coordinates(visit)), max.edge = c(40, 150), cutoff = 5
  )
}

#' Fit the presence model of a species with a spatial field.
#'
#' Bernoulli model with logit link of the presence of `species` on `term`,
#' effort (district, protocol, log fished surface), iid site and year
#' effects and a Matern spatial field (SPDE on `mesh`, none if `mesh` is
#' `NULL`). The intercept is coded explicitly in the field block of the
#' stack, so that it enters each visit once. The site log precision starts at 0: from INLA's default start
#' the fit can stay in a much worse mode where a short-range field replaces
#' the site effect.
#'
#' @param visit tibble, e.g. `prepare_dam_visits()`'s output.
#' @param presence tibble with `operation_id`, `species_valid`.
#' @param mesh INLA mesh, e.g. `build_presence_mesh()`'s output, or `NULL`
#'   for a model without spatial field.
#' @param species character, one species.
#' @param term character, right-hand side terms of the dam and position
#'   effects, e.g. "accumulation + log_alt".
#' @return list with `species`, `term` and the parts of the INLA fit used
#'   downstream, under their INLA names: `summary.fixed`,
#'   `summary.hyperpar`, `summary.random` (site, year and spatial if any),
#'   `summary.linear.predictor` (one row per visit, in the order of
#'   `visit`), `waic` (list with `waic`) and `mlik`.
#' @export
fit_presence_model <- function(visit, presence, mesh, species, term) {
  data <- visit |>
    dplyr::transmute(
      present = as.numeric(
        operation_id %in% presence$operation_id[presence$species_valid == species]
      ),
      accumulation, log_width, log_alt, log_dist, log_slope,
      district = factor(district), protocol = factor(protocol), log_surface,
      site, year
    )
  # District and protocol as 0/1 columns: without an implicit intercept, a
  # factor would get all its levels.
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
    # Each row of the projection matrix sums to 1, so the intercept in the
    # field block enters every visit once.
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
