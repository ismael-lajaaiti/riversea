#' Build the river visits used to relate dams to species and food webs.
#'
#' One row per river (ASPE) operation with fish records, a continuum score,
#' a fished surface and a catchment, as in D12. Position controls are
#' centred: log river width, log altitude and log distance to the mouth.
#' Dam accumulation is the number of large dams and weirs between the site
#' and the sea.
#'
#' @param continuum_operation tibble with `operation_id`, `district`,
#'   `log_altitude`, e.g. `score_continuum_pca()`'s output.
#' @param continuum_variables tibble with `operation_id`, `log_width`, e.g.
#'   `prepare_continuum_variables()`'s output.
#' @param river_operation tibble with `operation_id`, `sandre_code`, `date`,
#'   `year`.
#' @param river_effort tibble with `operation_id`, `protocol`, `surface`.
#' @param operation_mouth tibble with `operation_id`, `mouth_node`,
#'   `dist_mouth_m`, e.g. `assign_river_mouth()`'s output.
#' @param operation_dam tibble with `operation_id`, `roe_code`, `position`,
#'   e.g. `locate_operations_dams()`'s output.
#' @param large_dam sf tibble with `roe_code`, `kind`, e.g.
#'   `snap_dams_network()`'s output.
#' @param presence tibble with `operation_id`, `species_valid`.
#' @return tibble, one row per visit, with `dam_count`, `weir_count`, `dams`,
#'   the centred position controls, `log_surface` and the INLA indices
#'   `site` and `catchment`.
#' @export
prepare_dam_visits <- function(
  continuum_operation, continuum_variables, river_operation, river_effort,
  operation_mouth, operation_dam, large_dam, presence
) {
  dam_accumulation <- operation_dam |>
    dplyr::filter(position == "upstream") |>
    dplyr::left_join(
      sf::st_drop_geometry(large_dam) |> dplyr::select(roe_code, kind),
      by = "roe_code"
    ) |>
    dplyr::summarise(
      dam_count = sum(kind == "dam"), weir_count = sum(kind == "weir"),
      .by = operation_id
    )
  continuum_operation |>
    dplyr::select(operation_id, longitude, latitude, district, log_altitude) |>
    dplyr::filter(!is.na(log_altitude)) |>
    dplyr::semi_join(presence, by = "operation_id") |>
    dplyr::inner_join(
      dplyr::distinct(river_operation, operation_id, sandre_code, date, year),
      by = "operation_id"
    ) |>
    dplyr::left_join(
      dplyr::select(river_effort, operation_id, protocol, surface),
      by = "operation_id"
    ) |>
    dplyr::filter(!is.na(surface)) |>
    dplyr::left_join(dam_accumulation, by = "operation_id") |>
    dplyr::inner_join(
      dplyr::select(operation_mouth, operation_id, mouth_node, dist_mouth_m),
      by = "operation_id"
    ) |>
    dplyr::left_join(
      dplyr::select(continuum_variables, operation_id, log_width),
      by = "operation_id"
    ) |>
    dplyr::mutate(
      dam_count = dplyr::coalesce(dam_count, 0L),
      weir_count = dplyr::coalesce(weir_count, 0L),
      dams = dam_count + weir_count,
      accumulation = log1p(dams),
      log_dist = log(dist_mouth_m / 1000) - mean(log(dist_mouth_m / 1000)),
      log_alt = log_altitude - mean(log_altitude),
      log_width = log_width - mean(log_width),
      log_surface = log(surface),
      # Indices for INLA.
      catchment = as.integer(factor(mouth_node)),
      site = as.integer(factor(sandre_code))
    )
}

#' Formula terms shared by the dam models of species and food webs.
#'
#' Position controls (log width, log altitude with a quadratic term, log
#' distance to the mouth), effort (district, protocol, log fished surface)
#' and iid random effects of the site, the year and the catchment.
#'
#' @return character, right-hand side terms to add after the dam term.
#' @export
dam_model_terms <- function() {
  paste(
    "log_width + log_alt + I(log_alt^2) + log_dist +",
    "district + protocol + log_surface +",
    "f(site, model = 'iid') + f(year, model = 'iid') + f(catchment, model = 'iid')"
  )
}

#' Species found at enough sites to fit a presence model.
#'
#' @param presence tibble with `operation_id`, `species_valid`.
#' @param visit tibble with `operation_id`, `sandre_code`, e.g.
#'   `prepare_dam_visits()`'s output.
#' @param min_sites minimum number of distinct sites.
#' @return character vector of species, sorted.
#' @export
select_presence_species <- function(presence, visit, min_sites) {
  presence |>
    dplyr::inner_join(
      dplyr::distinct(visit, operation_id, sandre_code),
      by = "operation_id"
    ) |>
    dplyr::distinct(species_valid, sandre_code) |>
    dplyr::count(species_valid) |>
    dplyr::filter(n >= min_sites) |>
    dplyr::pull(species_valid) |>
    sort()
}

#' Fit the presence model of a species and draw from its posterior.
#'
#' Bernoulli model with logit link of the presence of `species` on
#' log(1 + dams) and `dam_model_terms()`. Only the posterior draws needed
#' downstream are kept: the linear predictor of each visit and the dam
#' effect, from the same joint draws.
#'
#' @param visit tibble, e.g. `prepare_dam_visits()`'s output.
#' @param presence tibble with `operation_id`, `species_valid`.
#' @param species character, one species.
#' @param n_draws number of joint posterior draws.
#' @return list with `species`, `operation_id` (order of the rows of `lp`),
#'   `lp` (visits x draws linear predictor), `beta` (draws of the dam
#'   effect), `summary` (posterior summary of the dam effect) and `waic`.
#' @export
sample_presence_dam <- function(visit, presence, species, n_draws) {
  visit$present <- as.numeric(
    visit$operation_id %in% presence$operation_id[presence$species_valid == species]
  )
  model <- INLA::inla(
    stats::as.formula(paste("present ~ accumulation +", dam_model_terms())),
    family = "binomial", data = visit, num.threads = 4,
    control.compute = list(config = TRUE, waic = TRUE)
  )
  set.seed(1)
  draw <- INLA::inla.posterior.sample(n_draws, model)
  name <- rownames(draw[[1]]$latent)
  predictor <- grep("^Predictor:", name)
  dam <- which(name == "accumulation:1")
  out <- list(
    species = species,
    operation_id = visit$operation_id,
    lp = round(sapply(draw, \(x) x$latent[predictor, 1]), 3),
    beta = sapply(draw, \(x) x$latent[dam, 1]),
    summary = model$summary.fixed["accumulation", ],
    waic = model$waic$waic
  )
  rm(model, draw)
  gc()
  out
}

#' Maximum prey-averaged trophic level of a food web.
#'
#' Same value as `get_trophic_length()`, by solving the linear system
#' directly (much faster than cheddar, for the many simulated webs): nodes
#' without prey have a trophic level of 1, others 1 + the mean trophic
#' level of their prey.
#'
#' @param web adjacency matrix, rows are prey and columns are predators.
#' @return numeric, `NA` if the system has no solution: a consumer whose
#'   only prey is itself, e.g. a pike alone (`get_trophic_length()` returns
#'   `-Inf` there).
#' @export
get_trophic_length_fast <- function(web) {
  n_prey <- colSums(web)
  P <- t(web)
  P[n_prey > 0, ] <- P[n_prey > 0, ] / n_prey[n_prey > 0]
  tryCatch(
    max(solve(diag(nrow(web)) - P, rep(1, nrow(web)))),
    error = \(e) NA_real_
  )
}

#' Structure of the species-level food web of a set of species.
#'
#' The web is cut out of the species-level metaweb: the fish species and
#' all resources.
#'
#' @param species character vector of fish species.
#' @param metaweb species-level adjacency matrix (rows prey, columns
#'   predators), e.g. `aggregate_species_web()`'s output.
#' @param resource character vector of resource nodes.
#' @return named numeric: `richness` (fish species), `trophic_length`,
#'   `piscivore` (1 if a fish eats a fish) and `fish_diet_overlap` (`NA`
#'   with fewer than 2 fish).
#' @export
measure_species_web <- function(species, metaweb, resource) {
  fish <- intersect(species, colnames(metaweb))
  web <- metaweb[c(fish, resource), c(fish, resource), drop = FALSE]
  c(
    richness = length(fish),
    trophic_length = get_trophic_length_fast(web),
    piscivore = as.numeric(length(fish) > 0 && sum(web[fish, fish]) > 0),
    fish_diet_overlap = if (length(fish) >= 2) get_fish_diet_overlap(web, resource) else NA_real_
  )
}

#' Simulate species-level food webs from the presence models.
#'
#' For each visit and posterior draw, the presence of each species is drawn
#' from its presence probability, and the structure of the resulting web is
#' measured, with the observed dams and with no dam (dam term set to 0). The
#' same random numbers are used in both scenarios, so that their difference
#' is the dam effect carried by the species composition.
#'
#' @param sample list of `sample_presence_dam()` outputs, one per species.
#' @param visit tibble with `operation_id`, `accumulation`.
#' @param metaweb species-level adjacency matrix.
#' @param resource character vector of resource nodes.
#' @return list with `operation_id` and `web`, an array visits x draws x
#'   metrics x scenarios ("dams", "no_dam").
#' @export
simulate_dam_webs <- function(sample, visit, metaweb, resource) {
  species <- vapply(sample, `[[`, character(1), "species")
  operation_id <- sample[[1]]$operation_id
  stopifnot(all(vapply(sample, \(x) identical(x$operation_id, operation_id), logical(1))))
  accumulation <- visit$accumulation[match(operation_id, visit$operation_id)]
  beta <- sapply(sample, `[[`, "beta")
  n_draws <- nrow(beta)
  metric <- c("richness", "trophic_length", "piscivore", "fish_diet_overlap")
  one_visit <- function(i) {
    lp <- sapply(sample, \(x) x$lp[i, ])
    set.seed(i)
    u <- matrix(stats::runif(length(lp)), nrow(lp))
    webs <- function(lp) {
      t(vapply(
        seq_len(n_draws),
        \(d) measure_species_web(species[u[d, ] < stats::plogis(lp[d, ])], metaweb, resource),
        numeric(length(metric))
      ))
    }
    array(c(webs(lp), webs(lp - beta * accumulation[i])), c(n_draws, length(metric), 2))
  }
  web <- parallel::mclapply(
    seq_along(operation_id), one_visit,
    mc.cores = max(1, parallel::detectCores() - 2)
  )
  web <- aperm(simplify2array(web), c(4, 1, 2, 3))
  dimnames(web) <- list(NULL, NULL, metric, c("dams", "no_dam"))
  list(operation_id = operation_id, web = web)
}

#' Structure of the observed species-level food webs.
#'
#' @param presence tibble with `operation_id`, `species_valid`.
#' @param operation_id visits to measure.
#' @param species modelled species: the observed webs are restricted to
#'   them, as the simulated ones.
#' @param metaweb species-level adjacency matrix.
#' @param resource character vector of resource nodes.
#' @return tibble with `operation_id` and the metrics of
#'   `measure_species_web()`.
#' @export
measure_observed_species_webs <- function(presence, operation_id, species, metaweb, resource) {
  present <- split(presence$species_valid, presence$operation_id)
  web <- t(vapply(
    as.character(operation_id),
    \(id) measure_species_web(intersect(present[[id]], species), metaweb, resource),
    numeric(4)
  ))
  dplyr::bind_cols(tibble::tibble(operation_id = operation_id), tibble::as_tibble(web))
}

#' Posterior mean presence probability of each species at each visit.
#'
#' @param sample list of `sample_presence_dam()` outputs, one per species.
#' @return matrix visits x species, with the visits as row names.
#' @export
mean_presence_probability <- function(sample) {
  probability <- sapply(sample, \(x) rowMeans(stats::plogis(x$lp)))
  dimnames(probability) <- list(
    sample[[1]]$operation_id, vapply(sample, `[[`, character(1), "species")
  )
  probability
}

#' Dam effect on food web structure, from one set of webs.
#'
#' Fits the D12 models of web structure on log(1 + dams) and
#' `dam_model_terms()`: fish richness (Poisson), trophic length (Gaussian,
#' webs with a piscivore) and diet overlap (Gaussian), on the visits with at
#' least 3 fish species. Returns draws of the marginal posterior of the dam
#' effect, so that the fits to several sets of generated webs can be pooled.
#'
#' @param visit tibble, e.g. `prepare_dam_visits()`'s output.
#' @param structure tibble with `operation_id` and the metrics of
#'   `measure_species_web()`.
#' @param n_draws number of draws of the dam effect per response.
#' @return tibble with `response` and `draw`.
#' @export
fit_structure_dam <- function(visit, structure, n_draws) {
  data <- visit |>
    dplyr::inner_join(structure, by = "operation_id") |>
    dplyr::filter(richness >= 3)
  family <- c(richness = "poisson", trophic_length = "gaussian", fish_diet_overlap = "gaussian")
  purrr::map_dfr(names(family), function(response) {
    fit_data <- if (response == "trophic_length") dplyr::filter(data, piscivore == 1) else data
    model <- INLA::inla(
      stats::as.formula(paste(response, "~ accumulation +", dam_model_terms())),
      family = family[[response]], data = fit_data, num.threads = 4,
      control.compute = list(config = FALSE)
    )
    draw <- INLA::inla.rmarginal(n_draws, model$marginals.fixed$accumulation)
    rm(model)
    gc()
    tibble::tibble(response = response, draw = draw)
  })
}

#' Structure of one realisation of the generated webs.
#'
#' @param simulated `simulate_dam_webs()`'s output.
#' @param realisation index of the posterior draw.
#' @param scenario "dams" or "no_dam".
#' @return tibble with `operation_id` and the metrics.
#' @export
extract_generated_webs <- function(simulated, realisation, scenario = "dams") {
  dplyr::bind_cols(
    tibble::tibble(operation_id = simulated$operation_id),
    tibble::as_tibble(simulated$web[, realisation, , scenario])
  )
}

#' Dam effect carried by the species composition, per posterior draw.
#'
#' For each draw, the mean difference in web structure between the webs
#' generated with the observed dams and with no dam, over the visits with at
#' least one dam, overall and by number of dams.
#'
#' @param simulated `simulate_dam_webs()`'s output.
#' @param visit tibble with `operation_id`, `dams`.
#' @return tibble with `dam_class`, `metric`, `draw` and `difference`.
#' @export
summarise_dam_counterfactual <- function(simulated, visit) {
  dams <- visit$dams[match(simulated$operation_id, visit$operation_id)]
  dam_class <- cut(dams, c(0, 1, 3, Inf), c("1", "2-3", "4+"))
  difference <- simulated$web[, , , "dams", drop = FALSE] - simulated$web[, , , "no_dam", drop = FALSE]
  metric <- dimnames(simulated$web)[[3]]
  purrr::map_dfr(c("all", levels(dam_class)), function(class) {
    keep <- if (class == "all") dams > 0 else dam_class %in% class
    purrr::map_dfr(metric, function(m) {
      tibble::tibble(
        dam_class = class, metric = m, draw = seq_len(dim(difference)[2]),
        difference = colMeans(difference[keep, , m, 1], na.rm = TRUE)
      )
    })
  })
}

#' Rank-weighted number of dams.
#'
#' The k-th dam between the site and the sea weighs `weight^(k - 1)`: 1
#' gives the number of dams, 0 only the first dam (presence of a dam),
#' values between them a saturating count and values above 1 an
#' accelerating one. The order of the dams does not matter, only their
#' number.
#'
#' @param dams number of dams.
#' @param weight rank weight, >= 0.
#' @return numeric, the weighted number of dams.
#' @export
rank_dam_impact <- function(dams, weight) {
  if (weight == 1) dams else (1 - weight^dams) / (1 - weight)
}

#' Fit a dam model with rank-weighted dams, for model comparison.
#'
#' @param data tibble with `dams`, the response and the columns of
#'   `dam_model_terms()`.
#' @param response name of the response column.
#' @param family INLA family.
#' @param weight rank weight of the dams, see `rank_dam_impact()`.
#' @return one-row tibble with `response`, `weight`, the posterior summary
#'   of the effect of log(1 + weighted dams), the WAIC and the log marginal
#'   likelihood.
#' @export
fit_dam_weight <- function(data, response, family, weight) {
  data$accumulation <- log1p(rank_dam_impact(data$dams, weight))
  model <- INLA::inla(
    stats::as.formula(paste(response, "~ accumulation +", dam_model_terms())),
    family = family, data = data, num.threads = 4,
    control.compute = list(config = FALSE, waic = TRUE)
  )
  effect <- model$summary.fixed["accumulation", ]
  out <- tibble::tibble(
    response = response, weight = weight, mean = effect$mean,
    lower = effect$`0.025quant`, upper = effect$`0.975quant`,
    waic = model$waic$waic, mlik = model$mlik[1, 1]
  )
  rm(model)
  gc()
  out
}

#' Compare dam weightings for the structure of a set of webs.
#'
#' @param visit tibble, e.g. `prepare_dam_visits()`'s output.
#' @param structure tibble with `operation_id` and the metrics of
#'   `measure_species_web()`.
#' @param weight rank weight of the dams.
#' @return tibble, one row per response (see `fit_dam_weight()`).
#' @export
weight_structure_dams <- function(visit, structure, weight) {
  data <- visit |>
    dplyr::inner_join(structure, by = "operation_id") |>
    dplyr::filter(richness >= 3)
  dplyr::bind_rows(
    fit_dam_weight(data, "richness", "poisson", weight),
    fit_dam_weight(dplyr::filter(data, piscivore == 1), "trophic_length", "gaussian", weight),
    fit_dam_weight(data, "fish_diet_overlap", "gaussian", weight)
  )
}

#' Compare dam weightings for the presence of a species.
#'
#' @param visit tibble, e.g. `prepare_dam_visits()`'s output.
#' @param presence tibble with `operation_id`, `species_valid`.
#' @param species character, one species.
#' @param weight rank weight of the dams.
#' @return one-row tibble (see `fit_dam_weight()`) with `species`.
#' @export
weight_presence_dams <- function(visit, presence, species, weight) {
  visit$present <- as.numeric(
    visit$operation_id %in% presence$operation_id[presence$species_valid == species]
  )
  fit_dam_weight(visit, "present", "binomial", weight) |>
    dplyr::mutate(species = species, .before = 1)
}

#' Simulate species-level food webs with one dam fewer downstream.
#'
#' As `simulate_dam_webs()` (same posterior draws and random numbers, so the
#' difference with its "dams" scenario is paired), with the number of dams
#' of each visit reduced by one. Removing any one of the dams downstream of
#' a site changes its dam count by one, so this single scenario gives the
#' effect of removing each dam at the visits upstream of it.
#'
#' @param sample list of `sample_presence_dam()` outputs, one per species.
#' @param visit tibble with `operation_id`, `dams`.
#' @param metaweb species-level adjacency matrix.
#' @param resource character vector of resource nodes.
#' @return list with `operation_id`, `web` (array visits x draws x metrics;
#'   `NA` at the visits without a dam) and `presence_change` (visits x
#'   species, posterior mean change in presence probability).
#' @export
simulate_one_dam_less <- function(sample, visit, metaweb, resource) {
  species <- vapply(sample, `[[`, character(1), "species")
  operation_id <- sample[[1]]$operation_id
  dams <- visit$dams[match(operation_id, visit$operation_id)]
  shift <- log1p(pmax(dams - 1, 0)) - log1p(dams)
  beta <- sapply(sample, `[[`, "beta")
  n_draws <- nrow(beta)
  metric <- c("richness", "trophic_length", "piscivore", "fish_diet_overlap")
  one_visit <- function(i) {
    if (dams[i] == 0) {
      return(list(web = matrix(NA_real_, n_draws, length(metric)), change = rep(0, length(species))))
    }
    lp <- sapply(sample, \(x) x$lp[i, ])
    lp_less <- lp + beta * shift[i]
    set.seed(i)
    u <- matrix(stats::runif(length(lp)), nrow(lp))
    web <- t(vapply(
      seq_len(n_draws),
      \(d) measure_species_web(species[u[d, ] < stats::plogis(lp_less[d, ])], metaweb, resource),
      numeric(length(metric))
    ))
    list(web = web, change = colMeans(stats::plogis(lp_less) - stats::plogis(lp)))
  }
  out <- parallel::mclapply(
    seq_along(operation_id), one_visit,
    mc.cores = max(1, parallel::detectCores() - 2)
  )
  web <- aperm(simplify2array(lapply(out, `[[`, "web")), c(3, 1, 2))
  dimnames(web) <- list(NULL, NULL, metric)
  change <- t(sapply(out, `[[`, "change"))
  dimnames(change) <- list(operation_id, species)
  list(operation_id = operation_id, web = web, presence_change = change)
}

#' Importance of each dam for the species and the food webs upstream.
#'
#' For each dam with sampled sites upstream, the change expected if it were
#' removed: the paired difference between the webs with one dam fewer and
#' with the observed dams, averaged over the visits of each site, then
#' summed (and averaged) over the sites upstream of the dam, for each
#' posterior draw.
#'
#' @param simulated `simulate_dam_webs()`'s output.
#' @param one_less `simulate_one_dam_less()`'s output.
#' @param operation_dam tibble with `operation_id`, `roe_code`, `position`.
#' @param visit tibble with `operation_id`, `sandre_code`.
#' @return list with `web` (tibble: `roe_code`, `sites`, `visits`, `metric`,
#'   `draw`, `total` and `per_site` change) and `species` (tibble:
#'   `roe_code`, `species`, `sites_gained`, the expected number of sites
#'   upstream where the species would be found in addition).
#' @export
summarise_dam_removal <- function(simulated, one_less, operation_dam, visit) {
  stopifnot(identical(simulated$operation_id, one_less$operation_id))
  index <- match(visit$operation_id, simulated$operation_id)
  upstream <- operation_dam |>
    dplyr::filter(position == "upstream", operation_id %in% visit$operation_id) |>
    dplyr::distinct(roe_code, operation_id) |>
    dplyr::left_join(dplyr::select(visit, operation_id, sandre_code), by = "operation_id")
  difference <- one_less$web - simulated$web[, , , "dams"]
  metric <- dimnames(difference)[[3]]
  site_of <- visit$sandre_code[match(simulated$operation_id, visit$operation_id)]
  # Site means of the differences, sites x draws x metrics.
  site <- sort(unique(site_of[!is.na(site_of)]))
  site_mean <- vapply(metric, function(m) {
    rowsum(difference[, , m], site_of, na.rm = TRUE)[site, , drop = FALSE] /
      as.vector(rowsum(as.numeric(!is.na(difference[, 1, m])), site_of)[site, ])
  }, matrix(0, length(site), dim(difference)[2]))
  change_site <- rowsum(one_less$presence_change, site_of)[site, , drop = FALSE] /
    as.vector(table(site_of)[site])
  dam_site <- split(match(upstream$sandre_code, site), upstream$roe_code) |> lapply(unique)
  web <- purrr::imap_dfr(dam_site, function(s, code) {
    purrr::map_dfr(metric, function(m) {
      x <- site_mean[s, , m, drop = FALSE][, , 1, drop = TRUE]
      x <- if (is.null(dim(x))) matrix(x, nrow = 1) else x
      tibble::tibble(
        roe_code = code, sites = length(s),
        visits = sum(upstream$roe_code == code), metric = m,
        draw = seq_len(ncol(x)), total = colSums(x, na.rm = TRUE),
        per_site = colMeans(x, na.rm = TRUE)
      )
    })
  })
  species <- purrr::imap_dfr(dam_site, function(s, code) {
    tibble::tibble(
      roe_code = code, species = colnames(change_site),
      sites_gained = colSums(change_site[s, , drop = FALSE])
    )
  })
  list(web = web, species = species)
}

#' Species holding the maximum trophic level of each food web.
#'
#' For the size-class webs of the visits and for their species-level webs
#' (modelled species only), the trophic level of each node is computed as
#' in `get_trophic_length_fast()`, and the species of the highest fish node
#' is kept.
#'
#' @param local named list of size-class adjacency matrices, e.g.
#'   `build_foodweb()`'s `local` output.
#' @param presence tibble with `operation_id`, `species_valid`.
#' @param operation_id visits.
#' @param species modelled species.
#' @param metaweb species-level adjacency matrix.
#' @param resource character vector of resource nodes.
#' @return tibble with `operation_id`, `top_size_class`, `top_species`
#'   (species of the top node in each web) and `trophic_length_size_class`,
#'   `trophic_length_species`.
#' @export
find_top_trophic_species <- function(local, presence, operation_id, species, metaweb, resource) {
  trophic_level <- function(web) {
    n_prey <- colSums(web)
    P <- t(web)
    P[n_prey > 0, ] <- P[n_prey > 0, ] / n_prey[n_prey > 0]
    level <- tryCatch(solve(diag(nrow(web)) - P, rep(1, nrow(web))), error = \(e) rep(NA_real_, nrow(web)))
    stats::setNames(level, rownames(web))
  }
  present <- split(presence$species_valid, presence$operation_id)
  purrr::map_dfr(intersect(as.character(operation_id), names(local)), function(id) {
    size_class <- trophic_level(local[[id]])
    size_class <- size_class[!names(size_class) %in% resource]
    fish <- intersect(intersect(present[[id]], species), colnames(metaweb))
    species_level <- trophic_level(metaweb[c(fish, resource), c(fish, resource), drop = FALSE])[fish]
    tibble::tibble(
      operation_id = id,
      top_size_class = sub("_[0-9]+$", "", names(which.max(size_class))),
      top_species = names(which.max(species_level)),
      trophic_length_size_class = max(size_class, na.rm = TRUE),
      trophic_length_species = max(species_level, na.rm = TRUE)
    )
  })
}

#' River reaches of the network, with the covariates of the dam models.
#'
#' RHT reaches of the kept districts, placed on the AMOBIO network by their
#' midpoint (as the operations), with their length, the position controls
#' of the visits (log width, log altitude, log distance to the mouth,
#' centred with the means of the visits) and their catchment. Only the
#' reaches within the range of the visits (1% to 99% quantiles of each
#' control) and within 200 m of the network are kept, to avoid
#' extrapolating.
#'
#' @param rht sf tibble of RHT reaches, e.g. `read_rht()`'s output.
#' @param basin sf tibble with `district`.
#' @param network_geo AMOBIO network with edge geometry, e.g.
#'   `add_edge_geometry()`'s output.
#' @param network AMOBIO network, e.g. `restrict_river_network()`'s output.
#' @param river_mouth tibble with `node_id`, `district`, `matched`.
#' @param visit tibble, e.g. `prepare_dam_visits()`'s output.
#' @param continuum_variables tibble with `operation_id`, `log_width`.
#' @return tibble, one row per reach: `reach_id`, `node_id`, `district`,
#'   `length_km`, `log_width`, `log_alt`, `log_dist` (centred),
#'   `catchment` (index of the visits' catchments, `NA` for others).
#' @export
prepare_network_reaches <- function(rht, basin, network_geo, network, river_mouth, visit, continuum_variables) {
  width_visit <- continuum_variables$log_width[match(visit$operation_id, continuum_variables$operation_id)]
  mean_width <- mean(width_visit)
  mean_alt <- mean(visit$log_altitude)
  mean_dist <- mean(log(visit$dist_mouth_m / 1000))
  district_kept <- unique(visit$district)
  midpoint <- sf::st_line_sample(sf::st_geometry(rht), sample = 0.5) |>
    sf::st_cast("POINT")
  reach <- sf::st_sf(
    reach_id = rht$id_drain, length_km = rht$long_ / 1000,
    log_altitude = log1p(rht$altitude), log_width = log(rht$l),
    geometry = midpoint
  ) |>
    # Simplified district outlines: the detailed coastlines make the join very slow.
    sf::st_join(
      sf::st_transform(basin, sf::st_crs(rht)) |> sf::st_simplify(dTolerance = 200),
      join = sf::st_within
    ) |>
    dplyr::filter(district %in% district_kept, is.finite(log_width), !is.na(log_altitude))
  coordinate <- sf::st_coordinates(sf::st_transform(reach, 4326))
  reach <- sf::st_drop_geometry(reach) |>
    dplyr::mutate(
      operation_id = as.character(reach_id), longitude = coordinate[, 1],
      latitude = coordinate[, 2], survey = "reach"
    )
  snapped <- snap_operations(reach, network_geo) |>
    dplyr::filter(snap_dist_m <= 200)
  mouth <- compute_distance_to_mouth(network, snapped, river_mouth, district_kept) |>
    dplyr::filter(is.finite(dist_mouth_m), dist_mouth_m > 0) |>
    assign_river_mouth(network, operation = _, river_mouth = river_mouth)
  visit_mouth <- levels(factor(visit$mouth_node))
  reach |>
    dplyr::inner_join(dplyr::select(mouth, operation_id, node_id, dist_mouth_m, mouth_node), by = "operation_id") |>
    dplyr::mutate(
      log_width = log_width - mean_width,
      log_alt = log_altitude - mean_alt,
      log_dist = log(dist_mouth_m / 1000) - mean_dist,
      catchment = match(mouth_node, visit_mouth)
    ) |>
    dplyr::filter(
      dplyr::between(log_width, stats::quantile(width_visit - mean_width, 0.01), stats::quantile(width_visit - mean_width, 0.99)),
      dplyr::between(log_alt, stats::quantile(visit$log_alt, 0.01), stats::quantile(visit$log_alt, 0.99)),
      dplyr::between(log_dist, stats::quantile(visit$log_dist, 0.01), stats::quantile(visit$log_dist, 0.99))
    ) |>
    dplyr::select(reach_id, node_id, district, length_km, log_width, log_alt, log_dist, catchment)
}

#' Fit the presence model of a species with rank-weighted dams, and draw
#' from the posterior of its fixed and catchment effects.
#'
#' As `sample_presence_dam()`, with the dams weighted by rank (see
#' `rank_dam_impact()`), keeping what is needed to predict at new places:
#' the fixed effects and the catchment effects, from joint draws. The site
#' and year effects are left out: predictions are for a typical site and
#' year.
#'
#' @param visit tibble, e.g. `prepare_dam_visits()`'s output.
#' @param presence tibble with `operation_id`, `species_valid`.
#' @param species character, one species.
#' @param weight rank weight of the dams.
#' @param n_draws number of joint posterior draws.
#' @return list with `species`, `weight`, `mlik`, `fixed` (draws x fixed
#'   effects) and `catchment` (draws x catchments).
#' @export
sample_presence_weight <- function(visit, presence, species, weight, n_draws) {
  visit$present <- as.numeric(
    visit$operation_id %in% presence$operation_id[presence$species_valid == species]
  )
  visit$accumulation <- log1p(rank_dam_impact(visit$dams, weight))
  model <- INLA::inla(
    stats::as.formula(paste("present ~ accumulation +", dam_model_terms())),
    family = "binomial", data = visit, num.threads = 4,
    control.compute = list(config = TRUE)
  )
  set.seed(1)
  draw <- INLA::inla.posterior.sample(n_draws, model)
  name <- rownames(draw[[1]]$latent)
  fixed_name <- paste0(rownames(model$summary.fixed), ":1")
  catchment_name <- grep("^catchment:", name, value = TRUE)
  latent <- sapply(draw, \(x) x$latent[c(fixed_name, catchment_name), 1])
  out <- list(
    species = species, weight = weight, mlik = model$mlik[1, 1],
    fixed = t(latent[fixed_name, , drop = FALSE]) |> `colnames<-`(rownames(model$summary.fixed)),
    catchment = t(latent[catchment_name, , drop = FALSE])
  )
  rm(model, draw)
  gc()
  out
}

#' Importance of each dam over the whole river network upstream.
#'
#' For each species, the presence probability is predicted on every reach
#' with its observed number of dams downstream and with one dam fewer, for a
#' typical site and year (reference protocol, median fished surface). The
#' weighting of the dams is averaged over the rank weights by their log
#' marginal likelihood (Bayesian model averaging): each posterior draw comes
#' from one weight, drawn with these probabilities. The gain of a dam is the
#' change in presence summed over the reaches upstream of it, weighted by
#' their length (species x km).
#'
#' @param sample list of `sample_presence_weight()` outputs (species x
#'   weights).
#' @param reach tibble, e.g. `prepare_network_reaches()`'s output.
#' @param reach_dam tibble with `operation_id` (reach), `roe_code`,
#'   `position`, e.g. `locate_operations_dams()`'s output on the reaches.
#' @param visit tibble of the visits the models were fitted on.
#' @return list with `richness` (tibble: `roe_code`, `draw`, `gain` in
#'   species x km), `species` (tibble: `roe_code`, `species`, mean `gain`
#'   in km), `upstream_km` (tibble: `roe_code`, `km`) and `model_weight`
#'   (tibble: `species`, `weight`, `probability`).
#' @export
predict_dam_removal_network <- function(sample, reach, reach_dam, visit) {
  upstream <- reach_dam |>
    dplyr::filter(position == "upstream") |>
    dplyr::distinct(reach_id = as.numeric(operation_id), roe_code)
  dams <- dplyr::count(upstream, reach_id, name = "dams")
  reach <- reach |>
    dplyr::left_join(dams, by = "reach_id") |>
    dplyr::mutate(dams = dplyr::coalesce(dams, 0L))
  dammed <- reach |> dplyr::filter(dams > 0)
  # Incidence dams x dammed reaches, weighted by reach length.
  code <- sort(unique(upstream$roe_code))
  pair <- upstream |> dplyr::filter(reach_id %in% dammed$reach_id)
  incidence <- Matrix::sparseMatrix(
    i = match(pair$roe_code, code), j = match(pair$reach_id, dammed$reach_id),
    x = dammed$length_km[match(pair$reach_id, dammed$reach_id)],
    dims = c(length(code), nrow(dammed))
  )
  # Design matrix of the fixed effects other than the dams, at the reaches.
  protocol <- names(which.max(table(visit$protocol)))
  design_data <- dammed |>
    dplyr::mutate(
      district = factor(district, levels(factor(visit$district))),
      protocol = factor(protocol, levels(factor(visit$protocol))),
      log_surface = stats::median(visit$log_surface)
    )
  design <- stats::model.matrix(
    ~ log_width + log_alt + I(log_alt^2) + log_dist + district + protocol + log_surface,
    design_data
  )
  species <- unique(vapply(sample, `[[`, character(1), "species"))
  model_weight <- purrr::map_dfr(sample, \(x) tibble::tibble(species = x$species, weight = x$weight, mlik = x$mlik)) |>
    dplyr::mutate(probability = exp(mlik - max(mlik)) / sum(exp(mlik - max(mlik))), .by = species)
  n_draws <- nrow(sample[[1]]$fixed)
  richness <- matrix(0, length(code), n_draws)
  gain_species <- list()
  set.seed(1)
  for (s in species) {
    fit <- sample[vapply(sample, \(x) x$species == s, logical(1))]
    probability <- model_weight$probability[model_weight$species == s][match(vapply(fit, `[[`, numeric(1), "weight"), model_weight$weight[model_weight$species == s])]
    pick <- sample(seq_along(fit), n_draws, replace = TRUE, prob = probability)
    change <- matrix(0, nrow(dammed), n_draws)
    for (k in unique(pick)) {
      d <- which(pick == k)
      beta <- fit[[k]]$fixed[d, , drop = FALSE]
      base <- design %*% t(beta[, colnames(design), drop = FALSE])
      catchment <- t(fit[[k]]$catchment[d, , drop = FALSE])[dammed$catchment, , drop = FALSE]
      catchment[is.na(catchment)] <- 0
      base <- base + catchment
      w <- fit[[k]]$weight
      now <- log1p(rank_dam_impact(dammed$dams, w))
      less <- log1p(rank_dam_impact(dammed$dams - 1, w))
      change[, d] <- stats::plogis(base + outer(less, beta[, "accumulation"])) -
        stats::plogis(base + outer(now, beta[, "accumulation"]))
    }
    gain <- as.matrix(incidence %*% change)
    richness <- richness + gain
    gain_species[[s]] <- tibble::tibble(roe_code = code, species = s, gain = rowMeans(gain))
  }
  upstream_km <- upstream |>
    dplyr::inner_join(dplyr::select(reach, reach_id, length_km), by = "reach_id") |>
    dplyr::summarise(km = sum(length_km), reaches = dplyr::n(), .by = roe_code)
  list(
    richness = tibble::tibble(
      roe_code = rep(code, n_draws), draw = rep(seq_len(n_draws), each = length(code)),
      gain = as.vector(richness)
    ),
    species = dplyr::bind_rows(gain_species),
    upstream_km = upstream_km,
    model_weight = dplyr::select(model_weight, species, weight, probability)
  )
}

#' Position and characteristics of the large dams in the river network.
#'
#' @param large_dam sf tibble, e.g. `snap_dams_network()`'s output.
#' @param dam_dam tibble of dams located relative to each other, e.g.
#'   `locate_operations_dams()` on the dams.
#' @param reach tibble, e.g. `prepare_network_reaches()`'s output.
#' @param reach_dam tibble with `operation_id` (reach), `roe_code`, e.g.
#'   `locate_operations_dams()` on the reaches.
#' @param rht sf tibble of RHT reaches (Strahler order, drained area).
#' @param network AMOBIO network with the dam nodes, e.g.
#'   `snap_dams_network()`'s output.
#' @param river_mouth tibble with `node_id`, `district`, `matched`.
#' @param basin sf tibble with `district`.
#' @return tibble, one row per dam placed on the network, with `rank` from
#'   the sea, `dams_upstream`, `dist_mouth_km`, `upstream_km` (reaches
#'   within the range of the visits), `catchment_km` and `upstream_share`
#'   (share of the river length of its catchment upstream of the dam), the
#'   `strahler` order and `drained_area` (km²) of the nearest RHT reach, and
#'   the ROE `kind`, `drop_height_min` and `altitude`.
#' @export
describe_dams <- function(large_dam, dam_dam, reach, reach_dam, rht, network, river_mouth, basin) {
  dam <- large_dam |>
    dplyr::filter(!is.na(node_id)) |>
    sf::st_transform(sf::st_crs(rht))
  dam <- dam |>
    sf::st_join(sf::st_simplify(sf::st_transform(basin, sf::st_crs(rht)), dTolerance = 200), join = sf::st_within) |>
    dplyr::filter(district %in% unique(reach$district))
  nearest <- sf::st_nearest_feature(dam, rht)
  coordinate <- sf::st_coordinates(sf::st_transform(dam, 4326))
  dam_point <- sf::st_drop_geometry(dam) |>
    dplyr::mutate(
      operation_id = roe_code, longitude = coordinate[, 1], latitude = coordinate[, 2],
      strahler = rht$strahler[nearest], drained_area = rht$surf_bv[nearest]
    )
  mouth <- compute_distance_to_mouth(network, dam_point, river_mouth, unique(dam_point$district)) |>
    assign_river_mouth(network, operation = _, river_mouth = river_mouth) |>
    dplyr::select(roe_code = operation_id, dist_mouth_m, mouth_node)
  reach_mouth <- reach |>
    dplyr::summarise(catchment_km = sum(length_km), .by = catchment)
  upstream <- reach_dam |>
    dplyr::distinct(reach_id = as.numeric(operation_id), roe_code) |>
    dplyr::inner_join(dplyr::select(reach, reach_id, length_km, catchment), by = "reach_id") |>
    dplyr::summarise(
      upstream_km = sum(length_km),
      catchment = dplyr::first(catchment[!is.na(catchment)], default = NA_integer_),
      .by = roe_code
    ) |>
    dplyr::left_join(reach_mouth, by = "catchment")
  rank <- dam_dam |>
    dplyr::filter(position == "upstream") |>
    dplyr::count(roe_code = operation_id, name = "dams_below")
  above <- dam_dam |>
    dplyr::filter(position == "downstream") |>
    dplyr::count(roe_code = operation_id, name = "dams_upstream")
  dam_point |>
    dplyr::select(roe_code, name, river, district, kind, drop_height_min, altitude, strahler, drained_area) |>
    dplyr::left_join(mouth, by = "roe_code") |>
    dplyr::left_join(rank, by = "roe_code") |>
    dplyr::left_join(above, by = "roe_code") |>
    dplyr::left_join(dplyr::select(upstream, roe_code, upstream_km, catchment_km), by = "roe_code") |>
    dplyr::mutate(
      rank = dplyr::coalesce(dams_below, 0L) + 1L,
      dams_upstream = dplyr::coalesce(dams_upstream, 0L),
      dist_mouth_km = dist_mouth_m / 1000,
      upstream_km = dplyr::coalesce(upstream_km, 0),
      upstream_share = upstream_km / catchment_km,
      altitude = suppressWarnings(as.numeric(altitude))
    ) |>
    dplyr::select(-dams_below, -dist_mouth_m)
}

#' Gain in the dendritic connectivity index from removing each dam.
#'
#' Cote et al. (2009), with the reaches as river length. Diadromous index:
#' a reach with D dams downstream counts its length times p^D, so removing
#' one of them adds length x p^(D - 1) (1 - p). Potamodromous index, with
#' impassable dams: removing a dam joins the free segment above it to the
#' one holding it, adding 2 x their product.
#'
#' @param reach tibble, e.g. `prepare_network_reaches()`'s output.
#' @param reach_dam tibble with `operation_id` (reach), `roe_code`,
#'   `river_dist_m`, e.g. `locate_operations_dams()` on the reaches.
#' @param passability numeric vector of dam passabilities p.
#' @return tibble, one row per dam, with `dci_d_<p>` (km) and `dci_p` (km²)
#'   gains and the river length `catchment_km` of its catchment.
#' @export
compute_dam_connectivity <- function(reach, reach_dam, passability) {
  reach <- reach |>
    dplyr::transmute(reach_id, length_km, catchment = paste(district, catchment))
  # Dams downstream of each reach, k = 1 for the dam nearest the sea.
  below <- reach_dam |>
    dplyr::distinct(reach_id = as.numeric(operation_id), roe_code, .keep_all = TRUE) |>
    dplyr::mutate(dams = dplyr::n(), k = rank(-river_dist_m, ties.method = "first"), .by = reach_id) |>
    dplyr::inner_join(reach, by = "reach_id")
  dci_d <- purrr::map(passability, \(p) {
    below |>
      dplyr::summarise(gain = sum(length_km * p^(dams - 1) * (1 - p)), .by = roe_code) |>
      dplyr::rename_with(\(x) paste0("dci_d_", p), "gain")
  }) |>
    purrr::reduce(dplyr::inner_join, by = "roe_code")
  # Free segment of each reach: the nearest dam downstream, or none.
  segment <- reach |>
    dplyr::left_join(
      below |> dplyr::filter(k == dams) |> dplyr::select(reach_id, segment = roe_code),
      by = "reach_id"
    ) |>
    dplyr::mutate(segment = dplyr::coalesce(segment, paste("none", catchment))) |>
    dplyr::summarise(segment_km = sum(length_km), .by = segment)
  next_below <- below |>
    dplyr::arrange(reach_id, k) |>
    dplyr::mutate(next_dam = dplyr::lag(roe_code), .by = reach_id) |>
    dplyr::distinct(roe_code, next_dam, catchment) |>
    dplyr::mutate(next_dam = dplyr::coalesce(next_dam, paste("none", catchment)))
  next_below |>
    dplyr::left_join(dplyr::rename(segment, above_km = segment_km), by = c("roe_code" = "segment")) |>
    dplyr::left_join(dplyr::rename(segment, holding_km = segment_km), by = c("next_dam" = "segment")) |>
    dplyr::left_join(
      dplyr::summarise(reach, catchment_km = sum(length_km), .by = catchment),
      by = "catchment"
    ) |>
    dplyr::transmute(
      roe_code, catchment_km,
      dci_p = 2 * dplyr::coalesce(above_km, 0) * dplyr::coalesce(holding_km, 0)
    ) |>
    dplyr::inner_join(dci_d, by = "roe_code")
}
