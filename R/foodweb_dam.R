#' Structure metrics of a local size-class food web.
#'
#' Trophic levels are prey-averaged, over the fish size classes: their
#' maximum is the trophic length, their mean tells whether the bulk of the
#' fish community feeds high or low in the web. Fish size classes sharing
#' the same prey and the same predators play the same trophic role.
#'
#' @param web adjacency matrix, rows are prey and columns are predators,
#'   fish nodes named `"species_class"`.
#' @param resource character, resource node names.
#' @param biomass named numeric, biomass of each fish node of `web`.
#' @return one-row tibble: `fish_richness`, `trophic_length`,
#'   `mean_trophic_level`, `frac_piscivorous` (share of fish size classes
#'   eating fish), `fish_diet_overlap`, `trophic_redundancy` (1 - number of
#'   trophic roles / number of fish size classes), `biomass_trophic_level`
#'   (mean trophic level weighted by biomass), `log_biomass` (log total fish
#'   biomass), `piscivore_biomass` (share of
#'   the fish biomass in size classes eating fish) and `biomass_pyramid`
#'   (correlation of log10 biomass with trophic level over the fish size
#'   classes, negative when biomass decreases up the web; `NA` with fewer
#'   than 3 distinct trophic levels).
#' @export
measure_web <- function(web, resource, biomass) {
  fish <- !colnames(web) %in% resource
  level <- compute_prey_averaged_trophic_level(
    igraph::graph_from_adjacency_matrix(web), resource
  )[fish]
  piscivorous <- colSums(web[fish, fish, drop = FALSE]) > 0
  role <- paste(
    apply(web[, fish, drop = FALSE], 2, paste, collapse = ""),
    apply(web[fish, , drop = FALSE], 1, paste, collapse = "")
  )
  mass <- biomass[colnames(web)[fish]]
  tibble::tibble(
    fish_richness = length(unique(sub("_.*", "", colnames(web)[fish]))),
    trophic_length = max(level),
    mean_trophic_level = mean(level),
    frac_piscivorous = mean(piscivorous),
    fish_diet_overlap = get_fish_diet_overlap(web, resource),
    trophic_redundancy = 1 - length(unique(role)) / sum(fish),
    biomass_trophic_level = stats::weighted.mean(level, mass),
    log_biomass = log(sum(mass)),
    piscivore_biomass = sum(mass[piscivorous]) / sum(mass),
    biomass_pyramid = if (length(unique(level)) > 2) {
      stats::cor(level, log10(mass))
    } else {
      NA_real_
    }
  )
}

#' Biomass of the fish size classes of each visit.
#'
#' Individuals are assigned to the size classes of their species (lengths
#' above the last bound go to the last class) and their weights summed.
#'
#' @param size tibble with `operation_id`, `species_valid`, `length` (cm)
#'   and `weight` (g), e.g. `get_river_size()`'s output.
#' @param size_class tibble with `species_code`, `lower_bound`,
#'   `upper_bound_1`, ..., e.g. `web_list$size_class`.
#' @return tibble with `operation_id`, `node` (`"species_class"`) and
#'   `biomass` (g).
#' @export
compute_node_biomass <- function(size, size_class) {
  bound <- size_class |>
    tidyr::pivot_longer(dplyr::starts_with("upper_bound_"), values_to = "upper") |>
    dplyr::summarise(upper = list(sort(upper)), .by = species_code)
  size |>
    dplyr::filter(!is.na(length), !is.na(weight)) |>
    dplyr::mutate(operation_id = as.character(operation_id)) |>
    dplyr::inner_join(bound, by = c(species_valid = "species_code")) |>
    dplyr::mutate(
      class = purrr::map2_int(length, upper, \(x, u) min(findInterval(x, u) + 1L, length(u))),
      node = paste(species_valid, class, sep = "_")
    ) |>
    dplyr::summarise(biomass = sum(weight), .by = c(operation_id, node))
}

#' Observed food webs of the dam visits and their structure.
#'
#' @param foodweb tibble with `operation_id` and the local webs in
#'   `foodweb`, e.g. the `river_foodweb_structure` target.
#' @param visit tibble with `operation_id`, e.g. `prepare_dam_visits()`'s
#'   output.
#' @param node_biomass tibble with `operation_id`, `node`, `biomass`, e.g.
#'   `compute_node_biomass()`'s output.
#' @param resource character, resource node names.
#' @return `visit` restricted to the visits with a web in which every fish
#'   size class has prey (a fish without prey would get the trophic level
#'   of a basal resource), with the web in `foodweb`, the node biomass in `biomass` (named numeric) and the
#'   metrics of [measure_web()].
#' @export
measure_dam_webs <- function(foodweb, visit, node_biomass, resource) {
  biomass <- node_biomass |>
    dplyr::summarise(biomass = list(stats::setNames(biomass, node)), .by = operation_id)
  web <- foodweb |>
    dplyr::filter(operation_id %in% visit$operation_id) |>
    dplyr::select(operation_id, foodweb) |>
    dplyr::inner_join(biomass, by = "operation_id") |>
    dplyr::filter(!purrr::map_lgl(foodweb, \(foodweb) {
      fish <- !colnames(foodweb) %in% resource
      any(colSums(foodweb[, fish, drop = FALSE]) == 0)
    }))
  metric <- parallel::mclapply(
    seq_len(nrow(web)),
    \(i) measure_web(web$foodweb[[i]], resource, web$biomass[[i]]),
    mc.cores = 4
  )
  .stop_on_mclapply_errors(metric, "measure_dam_webs")
  visit |>
    dplyr::inner_join(dplyr::bind_cols(web, dplyr::bind_rows(metric)), by = "operation_id")
}

#' Change of the web structure when one species is removed.
#'
#' Each fish species of each web is removed in turn (all its size classes
#' and their biomass), and the metrics of [measure_web()] are measured again
#' on the remaining web. Webs with fewer than 3 fish species are skipped.
#'
#' @param web tibble with `operation_id`, `sandre_code`, `foodweb` and
#'   `biomass`, e.g. `measure_dam_webs()`'s output.
#' @param resource character, resource node names.
#' @return tibble, one row per web, removed species and metric, with the
#'   metric of the full web (`full`) and of the web without the species
#'   (`removed`).
#' @export
remove_web_species <- function(web, resource) {
  out <- parallel::mclapply(seq_len(nrow(web)), function(i) {
    w <- web$foodweb[[i]]
    node_species <- sub("_.*", "", colnames(w))
    fish <- setdiff(unique(node_species), resource)
    if (length(fish) < 3) {
      return(NULL)
    }
    full <- measure_web(w, resource, web$biomass[[i]])
    purrr::map_dfr(fish, function(species) {
      keep <- node_species != species
      removed <- measure_web(w[keep, keep, drop = FALSE], resource, web$biomass[[i]])
      tibble::tibble(
        metric = names(full), full = unlist(full), removed = unlist(removed),
        species = species
      )
    }) |>
      dplyr::mutate(operation_id = web$operation_id[i], sandre_code = web$sandre_code[i])
  }, mc.cores = 4)
  .stop_on_mclapply_errors(out, "remove_web_species")
  dplyr::bind_rows(out)
}

#' Change of the first PCA axis when one species is removed.
#'
#' The score is scaled to unit variance over the PCA visits. Metric "PC1" is
#' the raw change; "PC1_richness" removes from it the change expected from
#' the loss of richness alone, to compare with a model at equal richness.
#'
#' @param removal tibble, `remove_web_species()`'s output.
#' @param pca `prcomp` object, `fit_structure_pca()`'s output.
#' @param richness_effect numeric, effect of log richness on the scaled PC1
#'   score, from the model of PC1 with richness.
#' @return tibble with `metric`, `species`, `operation_id`, `sandre_code`
#'   and the score of the web in `full` and without the species in
#'   `removed`, as input of [summarise_species_removal()].
#' @export
compute_removal_pc1 <- function(removal, pca, richness_effect) {
  pca_metric <- rownames(pca$rotation)
  score <- removal |>
    dplyr::filter(metric %in% pca_metric) |>
    tidyr::pivot_longer(c(full, removed), names_to = "web", values_to = "value") |>
    tidyr::pivot_wider(names_from = metric, values_from = value) |>
    dplyr::mutate(
      score = stats::predict(pca, dplyr::pick(dplyr::all_of(pca_metric)))[, "PC1"] / pca$sdev[1],
      log_richness = log(fish_richness)
    ) |>
    dplyr::select(species, operation_id, sandre_code, web, score, log_richness) |>
    tidyr::pivot_wider(names_from = web, values_from = c(score, log_richness))
  dplyr::bind_rows(
    dplyr::transmute(
      score, species, operation_id, sandre_code,
      metric = "PC1", full = score_full, removed = score_removed
    ),
    dplyr::transmute(
      score, species, operation_id, sandre_code,
      metric = "PC1_richness", full = score_full,
      removed = score_removed - richness_effect * (log_richness_removed - log_richness_full)
    )
  )
}

#' Mean change of the web structure when a species is removed.
#'
#' Changes are standardised by the standard deviation of the metric over
#' the webs, averaged over the visits of each site, then over the sites,
#' with a 91% bootstrap interval over the sites (500 resamples).
#'
#' @param removal tibble, e.g. `remove_web_species()`'s output.
#' @param min_sites integer, species kept: removed from webs of at least
#'   this number of sites.
#' @return tibble, one row per species and metric: `n_sites`, `mean`,
#'   `lower`, `upper` (standardised change) and `same_sign` (share of the
#'   webs changed in the direction of the mean).
#' @export
summarise_species_removal <- function(removal, min_sites) {
  removal |>
    dplyr::mutate(
      change = (removed - full) / stats::sd(full[!duplicated(operation_id)], na.rm = TRUE),
      .by = metric
    ) |>
    dplyr::filter(!is.na(change)) |>
    dplyr::summarise(
      pos = mean(change > 0), neg = mean(change < 0), change = mean(change), .by = c(species, metric, sandre_code)
    ) |>
    dplyr::filter(dplyr::n() >= min_sites, .by = c(species, metric)) |>
    dplyr::summarise(
      n_sites = dplyr::n(),
      mean = mean(change),
      boot = list(replicate(500, mean(sample(change, replace = TRUE)))),
      same_sign = if (mean > 0) mean(pos) else mean(neg),
      .by = c(species, metric)
    ) |>
    dplyr::mutate(
      lower = purrr::map_dbl(boot, \(x) stats::quantile(x, 0.045)),
      upper = purrr::map_dbl(boot, \(x) stats::quantile(x, 0.955)),
      .after = mean
    ) |>
    dplyr::select(-boot)
}

#' Fit a model of a web structure metric on the dams.
#'
#' Same structure as [fit_presence_model()]: `term`, effort (district,
#' protocol, log fished surface), iid site and year effects and a Matern
#' spatial field on `mesh` (none if `mesh` is `NULL`), with an identity
#' link for Gaussian metrics and a log link for Poisson ones.
#'
#' @param web tibble with the metric and the covariates of the visits,
#'   e.g. `measure_dam_webs()`'s output.
#' @param mesh INLA mesh, e.g. `build_presence_mesh()`'s output, or `NULL`
#'   for a model without spatial field.
#' @param metric character, column of `web` modelled.
#' @param family character, INLA likelihood, "gaussian" or "poisson".
#' @param term character, right-hand side terms of the dam and position
#'   effects, e.g. "accumulation + log_alt".
#' @return list with `metric`, `term` and the parts of the INLA fit used
#'   downstream, as `fit_presence_model()`, with 91% credible intervals
#'   (quantiles 0.045 and 0.955).
#' @export
fit_structure_model <- function(web, mesh, metric, family, term) {
  data <- web |>
    dplyr::transmute(
      response = .data[[metric]],
      accumulation, log_width, log_alt, log_dist, log_slope,
      log_richness = log(fish_richness) - mean(log(fish_richness)),
      district = factor(district), protocol = factor(protocol), log_surface,
      site, year
    )
  sampling <- stats::model.matrix(~ district + protocol, data)[, -1, drop = FALSE]
  colnames(sampling) <- make.names(colnames(sampling))
  covariate <- dplyr::bind_cols(
    dplyr::select(data, -response, -district, -protocol), tibble::as_tibble(sampling)
  )
  if (is.null(mesh)) {
    A <- list(1)
    effects <- list(dplyr::mutate(covariate, Intercept = 1))
    field_term <- ""
  } else {
    spde <- INLA::inla.spde2.pcmatern(
      mesh, prior.range = c(50, 0.5), prior.sigma = c(1, 0.01)
    )
    A <- list(INLA::inla.spde.make.A(mesh, loc = visit_coordinates(web)), 1)
    effects <- list(list(spatial = seq_len(spde$n.spde), Intercept = 1), covariate)
    field_term <- "+ f(spatial, model = spde)"
  }
  stack <- INLA::inla.stack(
    data = list(response = data$response), A = A, effects = effects
  )
  model <- INLA::inla(
    stats::as.formula(paste(
      "response ~ -1 + Intercept +", term, "+",
      paste(colnames(sampling), collapse = " + "), "+ log_surface +",
      "f(site, model = 'iid', hyper = list(prec = list(initial = 0))) +",
      "f(year, model = 'iid')", field_term
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
    summary.random = model$summary.random,
    summary.linear.predictor = model$summary.linear.predictor[seq_len(nrow(web)), ],
    waic = list(waic = model$waic$waic),
    mlik = model$mlik
  )
  rm(model)
  gc()
  out
}

#' Principal component analysis of the web metrics.
#'
#' Run on the scaled metrics of the visits where all are defined. The sign
#' of the first two axes is fixed so that a higher score means, on PC1, a
#' higher diet overlap and, on PC2, more fish species.
#'
#' @param web tibble with the metrics, e.g. `measure_dam_webs()`'s output.
#' @param metric character, columns of `web` used.
#' @return `prcomp` object.
#' @export
fit_structure_pca <- function(web, metric) {
  pca <- web |>
    dplyr::select(dplyr::all_of(metric)) |>
    stats::na.omit() |>
    stats::prcomp(scale. = TRUE)
  flip <- c(
    sign(pca$rotation["fish_diet_overlap", "PC1"]),
    sign(pca$rotation["fish_richness", "PC2"])
  )
  pca$rotation[, 1:2] <- sweep(pca$rotation[, 1:2], 2, flip, "*")
  pca$x[, 1:2] <- sweep(pca$x[, 1:2], 2, flip, "*")
  pca
}

#' Presence of the fish species in the webs.
#'
#' @param web tibble with `operation_id` and `foodweb`, e.g.
#'   `measure_dam_webs()`'s output.
#' @param resource character, resource node names.
#' @return visit x species 0/1 matrix.
#' @export
compute_species_presence <- function(web, resource) {
  species <- lapply(web$foodweb, \(foodweb) {
    unique(sub("_.*", "", setdiff(colnames(foodweb), resource)))
  })
  all_species <- sort(unique(unlist(species)))
  presence <- t(vapply(species, \(x) as.numeric(all_species %in% x), numeric(length(all_species))))
  dimnames(presence) <- list(web$operation_id, all_species)
  presence
}

#' Web structure change expected from the species lost to dams.
#'
#' Linear approximation: one more unit of `log(1 + dams)` changes the
#' presence probability of species s by beta_s p_s (1 - p_s) (logit link),
#' and its presence changes the metric by minus its removal effect. The
#' expected change of the metric is the sum over species of
#' beta_s mean(p_s (1 - p_s)) (-removal_s), with p_s the fitted presence
#' probability at the visits. Uncertainty is propagated by drawing the dam
#' effects (normal, posterior mean and sd) and the removal effects (normal,
#' sd from their 91% bootstrap interval).
#'
#' @param presence_model list of presence model fits, e.g. the
#'   `presence_model` target.
#' @param removal_effect tibble with `species`, `metric`, `mean`, `lower`,
#'   `upper`, e.g. `summarise_species_removal()`'s output.
#' @param n_draws integer, number of draws.
#' @return list with `species`, the contribution of each species to each
#'   metric (mean over draws), and `draw`, the expected change of each
#'   metric at each draw.
#' @export
decompose_dam_structure <- function(presence_model, removal_effect, n_draws) {
  presence <- purrr::map_dfr(presence_model, \(model) tibble::tibble(
    species = model$species,
    beta = model$summary.fixed["accumulation", "mean"],
    beta_sd = model$summary.fixed["accumulation", "sd"],
    slope = mean(stats::dlogis(model$summary.linear.predictor$mean))
  ))
  effect <- removal_effect |>
    dplyr::inner_join(presence, by = "species") |>
    dplyr::mutate(effect_sd = (upper - lower) / (2 * stats::qnorm(0.955)))
  draw <- effect |>
    tidyr::expand_grid(draw = seq_len(n_draws)) |>
    dplyr::mutate(
      beta_draw = stats::rnorm(dplyr::n(), beta, beta_sd),
      contribution = beta_draw * slope * -stats::rnorm(dplyr::n(), mean, effect_sd)
    )
  list(
    species = draw |>
      dplyr::mutate(presence_change = beta_draw * slope) |>
      dplyr::summarise(
        presence = mean(presence_change),
        presence_lower = stats::quantile(presence_change, 0.045),
        presence_upper = stats::quantile(presence_change, 0.955),
        lower = stats::quantile(contribution, 0.045),
        upper = stats::quantile(contribution, 0.955),
        contribution = mean(contribution),
        .by = c(species, metric)
      ),
    draw = draw |>
      dplyr::summarise(change = sum(contribution), .by = c(metric, draw))
  )
}

#' Expected structure change with the dam effects shuffled among species.
#'
#' Null model of [decompose_dam_structure()]: the dam effects on presence
#' are permuted among the species, keeping their values but not which
#' species carries them, so that the loss is as large but falls on random
#' species.
#'
#' @inheritParams decompose_dam_structure
#' @param n_permutations integer, number of permutations.
#' @return tibble with `metric`, `permutation` and the expected `change`
#'   (posterior means of the dam and removal effects).
#' @export
permute_dam_structure <- function(presence_model, removal_effect, n_permutations) {
  presence <- purrr::map_dfr(presence_model, \(model) tibble::tibble(
    species = model$species,
    beta = model$summary.fixed["accumulation", "mean"],
    slope = mean(stats::dlogis(model$summary.linear.predictor$mean))
  ))
  purrr::map_dfr(seq_len(n_permutations), function(k) {
    shuffled <- dplyr::mutate(presence, beta = sample(beta))
    removal_effect |>
      dplyr::inner_join(shuffled, by = "species") |>
      dplyr::summarise(change = sum(beta * slope * -mean), .by = metric) |>
      dplyr::mutate(permutation = k)
  })
}
