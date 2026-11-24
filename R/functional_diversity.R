#' Code the traits used for functional diversity.
#'
#' Effect traits of Dézerald et al. (diet, habitat use, size and life
#' history, reproduction, morphology), coded for [compute_functional_distance()]:
#' diet affinities as proportions (`diet.*`), rheophily as an ordered gradient
#' (1 limnophilic to 5 rheophilic), the vertical position as fuzzy affinities
#' to the bottom, the water column near the bottom and open water
#' (`living.*`), ordinal classes as integers and categories as characters.
#'
#' @param trait tibble, [read_dezerald_trait()]'s output.
#' @return tibble, one row per species, with `species` and the coded traits.
#' @export
code_functional_trait <- function(trait) {
  living <- c("BENT", "DEM", "PEL")
  diet <- trait |>
    dplyr::select(dplyr::starts_with("food.") & !dplyr::all_of("food.FFG")) |>
    as.matrix()
  living_affinity <- t(vapply(
    strsplit(trait$hab.LIVING, "..", fixed = TRUE),
    \(x) as.numeric(living %in% x) / length(x),
    numeric(length(living))
  ))
  trait |>
    dplyr::select(
      species, hab.TROPH, life.LONG, life.MAT, morp.LMAX, rep.ABSFEC,
      rep.REPMODE, rep.REPHAB, rep.CARE, rep.SPAWNMODE,
      morp.BLBD, morp.EHBD, morp.EHHD, morp.MOBD, morp.JLHD, morp.HDBD,
      morp.PFIBD, morp.PFLBD, morp.CFDCPD
    ) |>
    dplyr::mutate(
      hab.RHEO = match(trait$hab.RHEO, c("LIM", "EUR..LIM", "EUR", "EUR..RH", "RH")),
      .after = species
    ) |>
    dplyr::bind_cols(
      tibble::as_tibble(diet / rowSums(diet), .name_repair = \(x) sub("food", "diet", x)),
      tibble::as_tibble(living_affinity, .name_repair = \(x) paste0("living.", living))
    )
}

#' Functional distance between species.
#'
#' Gower-type distance: each variable gives a distance in \[0, 1\] (nominal:
#' 0 if same category, 1 otherwise; ordinal: absolute difference over the
#' range; fuzzy: half the sum of absolute differences of the affinities), the
#' distance of a group is the mean over its variables, and the distance
#' between species the mean over the groups, so that each group has the same
#' weight whatever its number of variables. Missing values are ignored pair
#' by pair.
#'
#' @param trait tibble, [code_functional_trait()]'s output.
#' @param variable tibble with `variable` (column of `trait`, or prefix of
#'   the columns of a fuzzy variable, e.g. "diet"), `group` and `type`
#'   ("nominal", "ordinal" or "fuzzy").
#' @return list with `distance` (species x species matrix) and `group`
#'   (array species x species x group of the group distances).
#' @export
compute_functional_distance <- function(trait, variable) {
  distance <- purrr::pmap(variable, \(variable, group, type) {
    if (type == "fuzzy") {
      x <- as.matrix(dplyr::select(trait, dplyr::starts_with(paste0(variable, "."))))
      return(as.matrix(stats::dist(x, method = "manhattan")) / 2)
    }
    x <- trait[[variable]]
    if (type == "ordinal") {
      return(abs(outer(x, x, "-")) / diff(range(x, na.rm = TRUE)))
    }
    outer(x, x, "!=") * 1
  })
  group <- vapply(
    unique(variable$group),
    \(g) apply(simplify2array(distance[variable$group == g]), c(1, 2), mean, na.rm = TRUE),
    matrix(0, nrow(trait), nrow(trait))
  )
  dimnames(group) <- list(trait$species, trait$species, unique(variable$group))
  list(distance = apply(group, c(1, 2), mean, na.rm = TRUE), group = group)
}

#' Functional richness and dispersion of each community.
#'
#' Functional richness is the volume of the convex hull of the species of the
#' community in the functional space, relative to the hull of all species
#' (NA with fewer species than axes + 1). Functional dispersion is the mean
#' distance of the species to the centroid of the community.
#'
#' @param coord matrix, species x axes, coordinates in the functional space
#'   (rownames: species).
#' @param presence tibble with `operation_id` and `species`.
#' @return tibble with `operation_id`, `richness`, `fric` and `fdis`.
#' @export
compute_functional_diversity <- function(coord, presence) {
  volume <- \(x) geometry::convhulln(x, options = "FA QJ")$vol
  pool <- volume(coord)
  presence |>
    dplyr::filter(species %in% rownames(coord)) |>
    dplyr::summarise(
      richness = dplyr::n(),
      fric = if (richness > ncol(coord)) volume(coord[species, , drop = FALSE]) / pool else NA_real_,
      fdis = mean(sqrt(rowSums(scale(coord[species, , drop = FALSE], scale = FALSE)^2))),
      .by = operation_id
    )
}

#' Functional diversity of communities with shuffled species traits.
#'
#' Null model keeping the richness of each community and the occupancy of
#' each species: the coordinates of the species in the functional space are
#' permuted among species.
#'
#' @inheritParams compute_functional_diversity
#' @param replicate integer, replicate number, also the random seed.
#' @return [compute_functional_diversity()]'s output with `replicate`.
#' @export
simulate_functional_diversity <- function(coord, presence, replicate) {
  set.seed(replicate)
  shuffled <- coord
  rownames(shuffled) <- sample(rownames(coord))
  compute_functional_diversity(shuffled, presence) |>
    dplyr::mutate(replicate = replicate)
}

#' Standardised effect size of functional diversity.
#'
#' @param observed tibble, [compute_functional_diversity()]'s output.
#' @param null tibble, [simulate_functional_diversity()]'s output for many
#'   replicates.
#' @return `observed` with, for `fric` and `fdis`, the standardised effect
#'   size (`*_ses`: observed minus null mean, over null sd) and the share of
#'   null values below the observed one (`*_rank`).
#' @export
compute_functional_ses <- function(observed, null) {
  null |>
    dplyr::summarise(
      dplyr::across(c(fric, fdis), list(mean = mean, sd = stats::sd)),
      fric_below = list(fric), fdis_below = list(fdis),
      .by = operation_id
    ) |>
    dplyr::inner_join(observed, by = "operation_id") |>
    dplyr::transmute(
      operation_id, richness, fric, fdis,
      fric_ses = (fric - fric_mean) / fric_sd,
      fdis_ses = (fdis - fdis_mean) / fdis_sd,
      fric_rank = purrr::map2_dbl(fric_below, fric, \(n, o) mean(n < o)),
      fdis_rank = purrr::map2_dbl(fdis_below, fdis, \(n, o) mean(n < o))
    )
}

#' Fit a model of a community index with dams and pressures.
#'
#' Same structure as [fit_pressure_presence_model()] (river position and
#' sampling controls, site and year iid effects, spatial field): Gaussian
#' for an index, binomial for a count out of `trials`. The year variance,
#' close to zero, gets a penalised complexity prior, and the optimisation is
#' restarted from its first mode: without these, repeated fits stopped at
#' different hyperparameters (marginal likelihood differing by ~20).
#'
#' @param visit tibble, e.g. the `pressure_visit_wide` target.
#' @param index tibble with `operation_id`, the column `response` and, for a
#'   binomial response, `trials`.
#' @param mesh INLA mesh, e.g. the `presence_mesh` target.
#' @param response character, column of `index`, e.g. "fric_ses".
#' @param term character, right-hand side of the fixed effects (columns of
#'   `visit` or `index`).
#' @param family character, "gaussian" or "binomial".
#' @param lincomb linear combinations of the fixed effects, e.g. from
#'   `INLA::inla.make.lincombs()`, or NULL.
#' @return list with `response`, `term`, `summary.fixed`,
#'   `summary.lincomb.derived` and `waic`.
#' @export
fit_functional_model <- function(visit, index, mesh, response, term, family = "gaussian", lincomb = NULL) {
  visit <- visit |>
    dplyr::inner_join(dplyr::rename(index, y = dplyr::all_of(response)), by = "operation_id") |>
    dplyr::filter(is.finite(y))
  trials <- if (family == "binomial") visit$trials else NULL
  data <- visit |>
    dplyr::transmute(
      y,
      dplyr::across(dplyr::all_of(all.vars(stats::as.formula(paste("~", term))))),
      district = factor(district), protocol = factor(protocol), log_surface,
      site, year
    )
  sampling <- stats::model.matrix(~ district + protocol, data)[, -1, drop = FALSE]
  colnames(sampling) <- make.names(colnames(sampling))
  data <- dplyr::bind_cols(dplyr::select(data, -district, -protocol), tibble::as_tibble(sampling))
  spde <- INLA::inla.spde2.pcmatern(mesh, prior.range = c(50, 0.5), prior.sigma = c(1, 0.01))
  stack <- INLA::inla.stack(
    data = list(y = data$y),
    A = list(INLA::inla.spde.make.A(mesh, loc = visit_coordinates(visit)), 1),
    effects = list(list(spatial = seq_len(spde$n.spde), Intercept = 1), dplyr::select(data, -y))
  )
  model <- INLA::inla(
    stats::as.formula(paste(
      "y ~ -1 + Intercept +", term, "+",
      paste(colnames(sampling), collapse = " + "), "+ log_surface +",
      "f(site, model = 'iid', hyper = list(prec = list(initial = 0))) +",
      "f(year, model = 'iid', hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01)))) +",
      "f(spatial, model = spde)"
    )),
    family = family, Ntrials = trials, data = INLA::inla.stack.data(stack), num.threads = 4,
    lincomb = lincomb,
    quantiles = c(0.045, 0.5, 0.955),
    control.predictor = list(A = INLA::inla.stack.A(stack)),
    control.compute = list(waic = TRUE),
    control.inla = list(tolerance = 1e-8), control.mode = list(restart = TRUE)
  )
  list(
    response = response, term = term, summary.fixed = model$summary.fixed,
    summary.lincomb.derived = model$summary.lincomb.derived, waic = model$waic$waic
  )
}

#' Size, life history and diet traits of the species, for community means.
#'
#' @param trait tibble, [code_functional_trait()]'s output.
#' @return tibble with `species` and `max_length`, `longevity`,
#'   `age_maturity`, `fecundity` (ordinal classes), `diet_breadth` (Shannon
#'   diversity of the diet affinities), and the diet proportions of
#'   `detritus`, `plants`, `small_prey` (< 1 mm), `medium_prey` (1-10 mm) and
#'   `large_prey` (> 10 mm).
#' @export
define_community_trait <- function(trait) {
  diet <- as.matrix(dplyr::select(trait, dplyr::starts_with("diet.")))
  trait |>
    dplyr::transmute(
      species,
      max_length = morp.LMAX, longevity = life.LONG, age_maturity = life.MAT, fecundity = rep.ABSFEC,
      diet_breadth = apply(diet, 1, \(p) -sum(p[p > 0] * log(p[p > 0]))),
      detritus = diet.DET1 + diet.DET2 + diet.DET3,
      plants = diet.PLANT4 + diet.PLANT5 + diet.PLANT6,
      small_prey = diet.ANI7, medium_prey = diet.ANI8, large_prey = diet.ANI9
    )
}

#' Community mean of each trait.
#'
#' @param trait tibble, [define_community_trait()]'s output.
#' @param presence tibble with `operation_id` and `species`.
#' @return tibble with `operation_id` and the mean of each trait over the
#'   species present (missing values ignored), centred and scaled.
#' @export
compute_community_trait_mean <- function(trait, presence) {
  presence |>
    dplyr::inner_join(trait, by = "species") |>
    dplyr::summarise(dplyr::across(-species, \(x) mean(x, na.rm = TRUE)), .by = operation_id) |>
    dplyr::mutate(dplyr::across(-operation_id, \(x) as.numeric(scale(x))))
}

#' Co-occurrence of functionally near-identical species along a pressure.
#'
#' Near-copies are the species pairs whose functional distance is in the
#' lowest `quantile` of all pairs. For each pair and each tercile of the
#' pressure, the observed share of visits where both species occur is
#' compared with the share expected if they occurred independently
#' (product of their frequencies in the tercile).
#'
#' @param distance list, [compute_functional_distance()]'s output.
#' @param presence tibble with `operation_id` and `species`.
#' @param visit tibble with `operation_id` and the column `pressure`.
#' @param pressure character, column of `visit`, e.g. "ammonium".
#' @param quantile numeric, share of the closest pairs kept, e.g. 0.05.
#' @return tibble with `species_1`, `species_2`, `distance`, `tercile`
#'   (1 low to 3 high), `observed`, `expected` and `visits`.
#' @export
summarise_near_copy <- function(distance, presence, visit, pressure, quantile) {
  d <- distance$distance
  pair <- which(upper.tri(d), arr.ind = TRUE)
  pair <- tibble::tibble(
    species_1 = rownames(d)[pair[, 1]], species_2 = rownames(d)[pair[, 2]], distance = d[pair]
  ) |>
    dplyr::filter(distance <= stats::quantile(distance, quantile))
  visit <- visit |>
    dplyr::transmute(operation_id, tercile = dplyr::ntile(.data[[pressure]], 3))
  occurrence <- presence |>
    dplyr::inner_join(visit, by = "operation_id") |>
    dplyr::mutate(visits = dplyr::n_distinct(operation_id), .by = tercile)
  frequency <- dplyr::summarise(occurrence, frequency = dplyr::n() / dplyr::first(visits), .by = c(species, tercile))
  pair |>
    dplyr::cross_join(dplyr::distinct(occurrence, tercile, visits)) |>
    dplyr::left_join(dplyr::rename(frequency, species_1 = species, frequency_1 = frequency), by = c("species_1", "tercile")) |>
    dplyr::left_join(dplyr::rename(frequency, species_2 = species, frequency_2 = frequency), by = c("species_2", "tercile")) |>
    dplyr::mutate(
      observed = purrr::pmap_dbl(list(species_1, species_2, tercile, visits), \(a, b, t, n) {
        site <- occurrence$operation_id[occurrence$tercile == t]
        sp <- occurrence$species[occurrence$tercile == t]
        length(intersect(site[sp == a], site[sp == b])) / n
      }),
      expected = dplyr::coalesce(frequency_1, 0) * dplyr::coalesce(frequency_2, 0)
    ) |>
    dplyr::select(species_1, species_2, distance, tercile, observed, expected, visits)
}

#' Number of functionally redundant species in each community.
#'
#' A species is redundant if another species of the community is a
#' near-copy: their functional distance is in the lowest `quantile` of the
#' distances between all species of the pool.
#'
#' @param distance list, [compute_functional_distance()]'s output.
#' @param presence tibble with `operation_id` and `species`.
#' @param quantile numeric, e.g. 0.05.
#' @return tibble with `operation_id`, `trials` (number of species) and
#'   `redundant` (number of species with a near-copy).
#' @export
compute_functional_redundancy <- function(distance, presence, quantile) {
  d <- distance$distance
  threshold <- stats::quantile(d[upper.tri(d)], quantile)
  close <- d <= threshold
  diag(close) <- FALSE
  presence |>
    dplyr::filter(species %in% rownames(d)) |>
    dplyr::summarise(
      trials = dplyr::n(),
      redundant = sum(rowSums(close[species, species, drop = FALSE]) > 0),
      .by = operation_id
    )
}

#' Dam visits with classes of a pressure.
#'
#' @param visit tibble, e.g. the `pressure_visit_wide` target.
#' @param pressure character, column of `visit`, e.g. "ammonium".
#' @param classes integer, number of classes (quantiles of the pressure).
#' @return `visit` with `pressure_class` and, for each class above the
#'   first, a sum-to-zero (effect) coded column `class_<k>`: 1 in class k,
#'   -1 in the first class, 0 otherwise, so that each coefficient is the
#'   deviation of a class from the average of the classes.
#' @export
add_pressure_class <- function(visit, pressure, classes) {
  visit <- dplyr::mutate(visit, pressure_class = dplyr::ntile(.data[[pressure]], classes))
  for (k in seq_len(classes)[-1]) {
    visit[[paste0("class_", k)]] <- (visit$pressure_class == k) - (visit$pressure_class == 1)
  }
  visit
}


#' One visit per site with its community and pressure class.
#'
#' @param visit tibble, e.g. the `pressure_visit_wide` target.
#' @param presence tibble with `operation_id` and `species`.
#' @param pressure character, column of `visit` used for the classes.
#' @param classes integer, number of classes (quantiles over sites).
#' @return list with `site` (the most recent visit of each site, with
#'   `class`) and `community` (site x species 0/1 matrix, same order).
#' @export
select_site_class <- function(visit, presence, pressure, classes) {
  site <- visit |>
    dplyr::slice_max(date, n = 1, with_ties = FALSE, by = site) |>
    dplyr::mutate(class = dplyr::ntile(.data[[pressure]], classes))
  community <- presence |>
    dplyr::filter(operation_id %in% site$operation_id) |>
    dplyr::mutate(present = 1) |>
    tidyr::pivot_wider(names_from = species, values_from = present, values_fill = 0)
  m <- as.matrix(community[match(site$operation_id, community$operation_id), -1])
  m[is.na(m)] <- 0
  list(site = site, community = m)
}

#' Dissimilarity between sites of the same pressure class.
#'
#' Pairs of sites in the same class and closer than `max_distance`
#' (straight line), with Sørensen dissimilarity, its turnover (Simpson) and
#' nestedness components (Baselga 2010), and the mean and absolute
#' difference of the covariates.
#'
#' @param site_class list, [select_site_class()]'s output.
#' @param covariate character, columns of the sites controlled, e.g.
#'   c("log_dist", "temperature").
#' @param max_distance numeric, km.
#' @return tibble, one row per pair: `site_1`, `site_2`, `class`,
#'   `sorensen`, `turnover`, `nestedness`, `log_distance`, `same_catchment`
#'   and `<covariate>_mean`, `<covariate>_diff`.
#' @export
compute_class_dissimilarity <- function(site_class, covariate, max_distance) {
  site <- site_class$site
  m <- site_class$community
  shared <- m %*% t(m)
  richness <- rowSums(m)
  geo <- as.matrix(stats::dist(visit_coordinates(site)))
  pair <- which(upper.tri(geo) & geo <= max_distance & outer(site$class, site$class, "=="), arr.ind = TRUE)
  i <- pair[, 1]
  j <- pair[, 2]
  a <- shared[pair]
  b <- richness[i] - a
  c <- richness[j] - a
  out <- tibble::tibble(
    site_1 = site$site[i], site_2 = site$site[j], class = site$class[i],
    sorensen = (b + c) / (2 * a + b + c),
    turnover = pmin(b, c) / (a + pmin(b, c)),
    log_distance = log(geo[pair]),
    same_catchment = as.numeric(site$catchment[i] == site$catchment[j])
  ) |>
    dplyr::mutate(nestedness = sorensen - turnover, .after = turnover)
  for (v in covariate) {
    out[[paste0(v, "_mean")]] <- (site[[v]][i] + site[[v]][j]) / 2
    out[[paste0(v, "_diff")]] <- abs(site[[v]][i] - site[[v]][j])
  }
  dplyr::filter(out, is.finite(sorensen), is.finite(turnover))
}

#' Distance-decay of dissimilarity per pressure class, with a bootstrap over sites.
#'
#' Linear model of each dissimilarity on class x log distance and the pair
#' covariates (centred and scaled), and its prediction for each class along
#' distance with the covariates at their mean. Replicate 0 uses all pairs;
#' other replicates resample sites with replacement and weight each pair by
#' the product of the number of times its sites are drawn.
#'
#' @param pair tibble, [compute_class_dissimilarity()]'s output.
#' @param replicate integer, bootstrap replicate (0: none), also the seed.
#' @return tibble with `replicate`, `response`, `class`, `distance` (km) and
#'   the predicted dissimilarity `predicted`.
#' @export
fit_class_dissimilarity_model <- function(pair, replicate) {
  covariate <- setdiff(names(pair), c("site_1", "site_2", "class", "sorensen", "turnover", "nestedness", "log_distance"))
  pair <- pair |>
    dplyr::mutate(class = factor(class), dplyr::across(dplyr::all_of(covariate), \(x) as.numeric(scale(x))))
  pair$weight <- 1
  if (replicate > 0) {
    set.seed(replicate)
    site <- unique(c(pair$site_1, pair$site_2))
    draw <- table(factor(sample(site, replace = TRUE), levels = site))
    pair$weight <- as.numeric(draw[as.character(pair$site_1)] * draw[as.character(pair$site_2)])
  }
  newdata <- tidyr::expand_grid(class = factor(levels(pair$class)), distance = c(5, 10, 25, 50, 100, 150, 200)) |>
    dplyr::mutate(log_distance = log(distance))
  for (v in covariate) newdata[[v]] <- 0
  purrr::map_dfr(c("sorensen", "turnover", "nestedness"), \(r) {
    fit <- stats::lm(
      stats::reformulate(c("class * log_distance", covariate), r),
      data = pair, weights = weight, subset = weight > 0
    )
    dplyr::mutate(newdata, replicate = replicate, response = r, predicted = stats::predict(fit, newdata)) |>
      dplyr::select(replicate, response, class, distance, predicted)
  })
}

#' Multivariate dispersion of communities per pressure class.
#'
#' Sites are placed by a principal coordinate analysis of the square root of
#' Sørensen dissimilarity; the covariates are partialled out (residuals of a linear model of each axis); the dispersion
#' of a class is the mean distance of its sites to the class centroid
#' (Anderson et al. 2006), with a bootstrap over sites.
#'
#' @param site_class list, [select_site_class()]'s output.
#' @param covariate character, columns of the sites partialled out.
#' @param replicate integer, number of bootstrap replicates.
#' @return tibble with `class`, `dispersion` and its 91% bootstrap interval
#'   (`lower`, `upper`).
#' @export
compute_class_dispersion <- function(site_class, covariate, replicate) {
  site <- site_class$site
  m <- site_class$community
  keep <- rowSums(m) > 0
  site <- site[keep, ]
  m <- m[keep, ]
  shared <- m %*% t(m)
  richness <- rowSums(m)
  sorensen <- 1 - 2 * shared / outer(richness, richness, "+")
  coord <- stats::cmdscale(stats::as.dist(sqrt(sorensen)), k = 20)
  residual <- stats::residuals(stats::lm(coord ~ ., data = as.data.frame(site[, covariate])))
  dispersion <- \(index) {
    tibble::tibble(class = site$class[index], row = index) |>
      dplyr::summarise(
        dispersion = mean(sqrt(rowSums(scale(residual[row, , drop = FALSE], scale = FALSE)^2))),
        .by = class
      )
  }
  set.seed(1)
  boot <- purrr::map_dfr(seq_len(replicate), \(r) {
    index <- unlist(lapply(split(seq_len(nrow(site)), site$class), \(x) sample(x, replace = TRUE)))
    dispersion(index)
  })
  dispersion(seq_len(nrow(site))) |>
    dplyr::left_join(
      dplyr::summarise(boot, lower = stats::quantile(dispersion, 0.045), upper = stats::quantile(dispersion, 0.955), .by = class),
      by = "class"
    ) |>
    dplyr::arrange(class)
}
