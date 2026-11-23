#' Response of each species to water quality and its occupancy.
#'
#' @param model list of presence model fits with a `PC1` term, e.g. the
#'   `species_additive_model` target.
#' @param visit tibble of the visits modelled, with `operation_id` and
#'   `site`, e.g. the `quality_visit` target.
#' @param presence tibble with `operation_id` and `species_valid`, one row
#'   per species present at a visit, e.g. the `river_presence` target.
#' @return tibble, one row per species: `beta`, `lower`, `upper` (effect of
#'   water quality on presence, logit scale, 91% credible interval) and
#'   `occupancy` (share of the sites where the species was caught at least
#'   once).
#' @export
summarise_species_response <- function(model, visit, presence) {
  occupancy <- presence |>
    dplyr::inner_join(dplyr::select(visit, operation_id, site), by = "operation_id") |>
    dplyr::distinct(species = species_valid, site) |>
    dplyr::count(species, name = "n_sites") |>
    dplyr::mutate(occupancy = n_sites / dplyr::n_distinct(visit$site))
  purrr::map_dfr(model, \(m) tibble::tibble(
    species = m$species,
    beta = m$summary.fixed["PC1", "mean"],
    lower = m$summary.fixed["PC1", "0.045quant"],
    upper = m$summary.fixed["PC1", "0.955quant"]
  )) |>
    dplyr::left_join(occupancy, by = "species")
}

#' Binary ecological traits of the fish species.
#'
#' @param trait tibble of the freshwaterecology.info traits, one row per
#'   species, e.g. the `fish_trait` target.
#' @return tibble with `species` and one logical column per trait (`NA` if
#'   unknown): `intolerant` and `tolerant` (water quality), `rheophilic`,
#'   `lithophilic` (spawning on gravel), `benthic` (feeding on the bottom)
#'   and `stenothermal` (narrow temperature range).
#' @export
define_fish_trait <- function(trait) {
  trait |>
    dplyr::transmute(
      species = esp_nom_latin,
      intolerant = `water quality tolerance` == "intol",
      tolerant = `water quality tolerance` == "tol",
      rheophilic = `rheophily (habitat)` %in% c("rheo", "rheo A", "rheo B"),
      lithophilic = `reproduction (fiBS)` == "lithophilic",
      benthic = `feeding habitat` == "ben",
      stenothermal = `temperature tolerance` == "steno"
    )
}

#' Number of species carrying each trait at each visit.
#'
#' @inheritParams summarise_species_response
#' @param trait tibble, `define_fish_trait()`'s output.
#' @return tibble with `operation_id`, `trait`, `n_trait` (species present
#'   with the trait) and `n_known` (species present whose trait is known).
#' @export
count_visit_trait <- function(visit, presence, trait) {
  presence |>
    dplyr::filter(operation_id %in% visit$operation_id) |>
    dplyr::distinct(operation_id, species = species_valid) |>
    dplyr::inner_join(trait, by = "species") |>
    tidyr::pivot_longer(-c(operation_id, species), names_to = "trait", values_to = "value") |>
    dplyr::filter(!is.na(value)) |>
    dplyr::summarise(n_trait = sum(value), n_known = dplyr::n(), .by = c(operation_id, trait))
}

#' Fit a model of the share of species carrying a trait.
#'
#' Binomial model of the number of species with the trait among the species
#' present whose trait is known, with the covariates, effort, site and year
#' effects and spatial field of [fit_stressor_structure_model()].
#'
#' @param visit tibble of the visits, e.g. the `quality_visit` target.
#' @param mesh INLA mesh, e.g. the `presence_mesh` target.
#' @param count tibble, `count_visit_trait()`'s output for one trait.
#' @param term character, fixed effects, e.g. "accumulation + PC1 + log_width".
#' @return list with `trait`, `term`, `summary.fixed`, `summary.hyperpar` and
#'   `waic`.
#' @export
fit_trait_model <- function(visit, mesh, count, term) {
  visit <- dplyr::inner_join(visit, count, by = "operation_id")
  data <- visit |>
    dplyr::transmute(
      n_trait, n_known,
      accumulation, PC1, log_width, log_alt, log_dist, log_slope,
      district = factor(district), protocol = factor(protocol), log_surface,
      site, year
    )
  sampling <- stats::model.matrix(~ district + protocol, data)[, -1, drop = FALSE]
  colnames(sampling) <- make.names(colnames(sampling))
  covariate <- dplyr::bind_cols(
    dplyr::select(data, -n_trait, -n_known, -district, -protocol), tibble::as_tibble(sampling)
  )
  spde <- INLA::inla.spde2.pcmatern(
    mesh, prior.range = c(50, 0.5), prior.sigma = c(1, 0.01)
  )
  stack <- INLA::inla.stack(
    data = list(n_trait = data$n_trait, n_known = data$n_known),
    A = list(INLA::inla.spde.make.A(mesh, loc = visit_coordinates(visit)), 1),
    effects = list(list(spatial = seq_len(spde$n.spde), Intercept = 1), covariate)
  )
  stack_data <- INLA::inla.stack.data(stack)
  model <- INLA::inla(
    stats::as.formula(paste(
      "n_trait ~ -1 + Intercept +", term, "+",
      paste(colnames(sampling), collapse = " + "), "+ log_surface +",
      "f(site, model = 'iid', hyper = list(prec = list(initial = 0))) +",
      "f(year, model = 'iid') + f(spatial, model = spde)"
    )),
    family = "binomial", Ntrials = stack_data$n_known, data = stack_data, num.threads = 4,
    quantiles = c(0.045, 0.5, 0.955),
    control.predictor = list(A = INLA::inla.stack.A(stack)),
    control.compute = list(config = FALSE, waic = TRUE)
  )
  out <- list(
    trait = unique(count$trait),
    term = term,
    summary.fixed = model$summary.fixed,
    summary.hyperpar = model$summary.hyperpar,
    waic = list(waic = model$waic$waic)
  )
  rm(model)
  gc()
  out
}

#' Read the fish trait table of Dézerald et al.
#'
#' Species names are validated against FishBase, as the ASPE species list
#' (`species_valid`, see [get_river_size()]), so that synonyms match (e.g.
#' *Aspius aspius* is *Leuciscus aspius*); ordinal traits coded with a prefix (e.g. "ls3") or as
#' a range (upper thermal tolerance) are converted to integers.
#'
#' @param file character, path to `dezerald_trait.csv`.
#' @return tibble, one row per species: `species`, `code` (ASPE code),
#'   `abundance` (share of the individuals in the source data, %), and the
#'   traits, named as in the file (e.g. `tol.OXY`, `food.ANI8`).
#' @export
read_dezerald_trait <- function(file) {
  ordinal <- c(
    "mig.DWC", "mig.DAC", "hab.RESP", "life.LONG", "life.MAT", "morp.LARVLENGTH",
    "morp.GILLRAK", "rep.REPPER", "rep.REPCYCIND"
  )
  utils::read.csv(
    file, sep = ";", fileEncoding = "latin1", check.names = FALSE,
    na.strings = c("NA", "")
  ) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      species = dplyr::coalesce(rfishbase::validate_names(sciname), sciname),
      dplyr::across(dplyr::all_of(ordinal), \(x) as.integer(sub("^[^0-9]*", "", x))),
      rep.UPTHERMTOL = as.integer(sub(".*\\.\\.", "", rep.UPTHERMTOL))
    ) |>
    dplyr::select(
      species, code = Code_TAXON, abundance = Percent_Eff1,
      dplyr::starts_with(c("mig.", "food.", "hab.", "life.", "morp.", "tol.", "rep.", "ice."))
    )
}

#' Read the trait descriptions of Dézerald et al.
#'
#' @param file character, path to `dezerald_metadata.xlsx` (Table S2, one
#'   sheet).
#' @return tibble, one row per trait category: `trait` (column name in
#'   [read_dezerald_trait()]'s output, e.g. "tol.OXY"), `name`, `source`,
#'   `category` and `explanation`.
#' @export
read_dezerald_metadata <- function(file) {
  dir <- tempfile()
  utils::unzip(file, files = c("xl/sharedStrings.xml", "xl/worksheets/sheet1.xml"), exdir = dir)
  on.exit(unlink(dir, recursive = TRUE))
  shared <- xml2::read_xml(file.path(dir, "xl/sharedStrings.xml"))
  xml2::xml_ns_strip(shared)
  string <- vapply(
    xml2::xml_find_all(shared, "//si"),
    \(si) paste(xml2::xml_text(xml2::xml_find_all(si, ".//t")), collapse = ""), ""
  )
  sheet <- xml2::read_xml(file.path(dir, "xl/worksheets/sheet1.xml"))
  xml2::xml_ns_strip(sheet)
  cell <- xml2::xml_find_all(sheet, "//sheetData/row/c")
  value <- xml2::xml_text(xml2::xml_find_first(cell, "v"))
  shared_cell <- xml2::xml_attr(cell, "t") %in% "s"
  value[shared_cell] <- string[as.integer(value[shared_cell]) + 1]
  reference <- xml2::xml_attr(cell, "r")
  tibble::tibble(
    row = as.integer(sub("^[A-Z]+", "", reference)),
    column = sub("[0-9]+$", "", reference),
    value = trimws(value)
  ) |>
    dplyr::filter(row > 2, column %in% c("A", "B", "C", "D", "E", "F")) |>
    tidyr::pivot_wider(names_from = column, values_from = value) |>
    dplyr::arrange(row) |>
    dplyr::transmute(
      source = A, name = B, trait = C, category = D,
      explanation = dplyr::if_else(is.na(F), E, paste(E, F, sep = " "))
    ) |>
    tidyr::fill(source, name, trait) |>
    dplyr::mutate(trait = dplyr::if_else(trait == "food", paste0("food.", category), trait))
}

#' Water quality variable predicted at the fish visits.
#'
#' Same model and prediction as [predict_water_quality_visit()], for one
#' variable of the station-year means instead of the water quality index.
#'
#' @param station_year tibble of station-year means, e.g. the
#'   `naiades_station_year` target.
#' @param water_quality tibble with `code_station`, `year`, `X`, `Y`, e.g.
#'   `compute_water_quality()`'s output (station coordinates).
#' @param mesh INLA mesh, e.g. `build_water_quality_mesh()`'s output.
#' @param visit tibble with `operation_id`, `sandre_code`, `year`.
#' @param variable character, column of `station_year`, e.g. "nitrate".
#' @return tibble with `operation_id`, `pressure` (= `variable`) and the
#'   posterior `mean` and `sd` of the variable, on its scale in
#'   `station_year` (log10 for nutrients).
#' @export
predict_pressure_visit <- function(station_year, water_quality, mesh, visit, variable) {
  station_year <- dplyr::filter(station_year, !is.na(.data[[variable]]))
  centre <- mean(station_year[[variable]])
  spread <- stats::sd(station_year[[variable]])
  data <- water_quality |>
    dplyr::select(code_station, year, X, Y) |>
    dplyr::inner_join(
      dplyr::transmute(station_year, code_station, year, PC1 = (.data[[variable]] - centre) / spread),
      by = c("code_station", "year")
    )
  predict_water_quality_visit(data, mesh, visit) |>
    dplyr::transmute(operation_id, pressure = variable, mean = centre + spread * mean, sd = spread * sd)
}

#' Dam visits with the water quality pressures.
#'
#' @param visit tibble of the dam visits, e.g. the `dam_visit` target.
#' @param pressure tibble, `predict_pressure_visit()`'s output for several
#'   variables.
#' @return `visit` restricted to the visits with all pressures, with one
#'   column per pressure (posterior mean), centred and scaled.
#' @export
prepare_pressure_visit <- function(visit, pressure) {
  wide <- pressure |>
    dplyr::select(operation_id, pressure, mean) |>
    tidyr::pivot_wider(names_from = pressure, values_from = mean) |>
    tidyr::drop_na()
  visit |>
    dplyr::inner_join(wide, by = "operation_id") |>
    dplyr::mutate(dplyr::across(dplyr::all_of(unique(pressure$pressure)), \(x) as.numeric(scale(x))))
}

#' Fit the presence model of a species with dams and water quality pressures.
#'
#' Same model as [fit_stressor_presence_model()], with any column of `visit`
#' available to `term` instead of the water quality index, e.g. the
#' pressures of [prepare_pressure_visit()] and the local dam proximities of
#' [add_local_dam()].
#'
#' @param visit tibble, e.g. `prepare_pressure_visit()`'s output.
#' @inheritParams fit_stressor_presence_model
#' @return list as [fit_stressor_presence_model()]'s output.
#' @export
fit_pressure_presence_model <- function(visit, presence, mesh, species, term) {
  data <- visit |>
    dplyr::transmute(
      present = as.numeric(
        operation_id %in% presence$operation_id[presence$species_valid == species]
      ),
      dplyr::across(dplyr::all_of(all.vars(stats::as.formula(paste("~", term))))),
      district = factor(district), protocol = factor(protocol), log_surface,
      site, year
    )
  sampling <- stats::model.matrix(~ district + protocol, data)[, -1, drop = FALSE]
  colnames(sampling) <- make.names(colnames(sampling))
  data <- dplyr::bind_cols(
    dplyr::select(data, -district, -protocol), tibble::as_tibble(sampling)
  )
  spde <- INLA::inla.spde2.pcmatern(
    mesh, prior.range = c(50, 0.5), prior.sigma = c(1, 0.01)
  )
  stack <- INLA::inla.stack(
    data = list(present = data$present),
    A = list(INLA::inla.spde.make.A(mesh, loc = visit_coordinates(visit)), 1),
    effects = list(list(spatial = seq_len(spde$n.spde), Intercept = 1), dplyr::select(data, -present))
  )
  model <- INLA::inla(
    stats::as.formula(paste(
      "present ~ -1 + Intercept +", term, "+",
      paste(colnames(sampling), collapse = " + "), "+ log_surface +",
      "f(site, model = 'iid', hyper = list(prec = list(initial = 0))) +",
      "f(year, model = 'iid') + f(spatial, model = spde)"
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
    summary.linear.predictor = model$summary.linear.predictor[seq_len(nrow(visit)), ],
    waic = list(waic = model$waic$waic)
  )
  rm(model)
  gc()
  out
}

#' Effects of the dams and pressures on the presence of each species.
#'
#' @param model list of presence model fits, e.g. the
#'   `species_pressure_model` target.
#' @param term character, terms kept, e.g. c("accumulation", "nitrate").
#' @return tibble with `species`, `term`, the posterior `mean` and `sd` of
#'   the effect (logit scale) and its 91% credible interval (`lower`,
#'   `upper`).
#' @export
summarise_pressure_response <- function(model, term) {
  purrr::map_dfr(model, \(m) {
    m$summary.fixed[term, ] |>
      tibble::as_tibble(rownames = "term") |>
      dplyr::transmute(
        species = m$species, term, mean, sd, lower = `0.045quant`, upper = `0.955quant`
      )
  })
}

#' Traits of the a priori trait-pressure hypotheses.
#'
#' @param trait tibble, [read_dezerald_trait()]'s output.
#' @return tibble with `species` and one numeric column per trait (binary
#'   0/1, or the ordinal class for `crossing`): `intolerant` (water
#'   quality), `omnivore`, `lithophilic` (spawns on gravel),
#'   `oxygen_intolerant`, `stenothermal`, `rheophilic`,
#'   `habitat_intolerant`, `diadromous` and `crossing` (largest obstacle
#'   height crossed, ICE class: 1 below 0.5 m, 2 below 1 m, 3 below 3 m).
#' @export
code_trait_predictor <- function(trait) {
  trait |>
    dplyr::transmute(
      species,
      intolerant = as.numeric(tol.GEN == "INTOL"),
      omnivore = as.numeric(food.FFG %in% c("OMNI", "OMNI..INSV")),
      lithophilic = as.numeric(rep.REPMODE == "LITH"),
      oxygen_intolerant = as.numeric(tol.OXY == "O2INTOL"),
      stenothermal = as.numeric(tol.THERM == "STTHER"),
      rheophilic = as.numeric(hab.RHEO == "RH"),
      habitat_intolerant = as.numeric(tol.HAB == "HINTOL"),
      diadromous = as.numeric(mig.MIGRATION %in% c("LONG..LMA", "LONG..LMC")),
      crossing = as.numeric(match(ice.CLASS0, c("0..0.5m", "0.5..1m", "1..3m")))
    )
}

#' Fit a meta-regression of species responses on a trait.
#'
#' Random-effects meta-regression: the effect of a pressure on the presence
#' of species s is beta_s = gamma_0 + gamma trait_s + u_s + e_s, with e_s of
#' known sd (posterior sd of beta_s) and u_s an iid species effect for the
#' variation not explained by the trait.
#'
#' @param response tibble, [summarise_pressure_response()]'s output.
#' @param predictor tibble, [code_trait_predictor()]'s output.
#' @param term character, pressure, e.g. "nitrate".
#' @param trait character, column of `predictor`, e.g. "omnivore".
#' @return tibble, one row: `term`, `trait`, `n` (species), the posterior
#'   `mean` of gamma and its 91% credible interval (`lower`, `upper`).
#' @export
fit_trait_response <- function(response, predictor, term, trait) {
  data <- response |>
    dplyr::filter(.data$term == !!term) |>
    dplyr::inner_join(dplyr::select(predictor, species, value = dplyr::all_of(trait)), by = "species") |>
    dplyr::filter(!is.na(value)) |>
    dplyr::mutate(species_index = dplyr::row_number())
  model <- INLA::inla(
    mean ~ value + f(species_index, model = "iid"),
    family = "gaussian", data = data, scale = 1 / data$sd^2,
    control.family = list(hyper = list(prec = list(initial = 0, fixed = TRUE))),
    quantiles = c(0.045, 0.5, 0.955), num.threads = 1
  )
  tibble::tibble(
    term = term, trait = trait, n = nrow(data),
    mean = model$summary.fixed["value", "mean"],
    lower = model$summary.fixed["value", "0.045quant"],
    upper = model$summary.fixed["value", "0.955quant"]
  )
}

#' Proximity of the nearest large dam downstream and upstream of each visit.
#'
#' Proximity is -log10(distance / 100 km) along the river, 0 beyond 100 km
#' or without dam: 1 at 10 km, 2 at 1 km. `dam_below` is the nearest dam
#' downstream of the site (the site may be in or above its reservoir),
#' `dam_above` the nearest dam upstream (the site is below it).
#'
#' @param visit tibble with `operation_id`, e.g. [prepare_pressure_visit()]'s
#'   output.
#' @param operation_dam tibble with `operation_id`, `position` and
#'   `river_dist_m`, e.g. `locate_operations_dams()`'s output.
#' @return `visit` with `dam_below` and `dam_above`, centred and scaled, and
#'   their raw values in `dam_below_raw` and `dam_above_raw`.
#' @export
add_local_dam <- function(visit, operation_dam) {
  proximity <- operation_dam |>
    dplyr::filter(operation_id %in% visit$operation_id) |>
    dplyr::summarise(distance = min(river_dist_m), .by = c(operation_id, position)) |>
    dplyr::mutate(
      proximity = pmax(0, -log10(distance / 1e5)),
      position = dplyr::if_else(position == "upstream", "dam_below_raw", "dam_above_raw")
    ) |>
    dplyr::select(-distance) |>
    tidyr::pivot_wider(names_from = position, values_from = proximity)
  visit |>
    dplyr::left_join(proximity, by = "operation_id") |>
    dplyr::mutate(
      dam_below_raw = dplyr::coalesce(dam_below_raw, 0),
      dam_above_raw = dplyr::coalesce(dam_above_raw, 0),
      dam_below = as.numeric(scale(dam_below_raw)),
      dam_above = as.numeric(scale(dam_above_raw))
    )
}

#' Summer summaries of the NAIADES parameters per station and year.
#'
#' Summer is June to September; a station-year is kept for a parameter with
#' at least 3 summer measurements. Nutrients are log10-transformed.
#'
#' @inheritParams compute_station_year_mean
#' @return tibble with `code_station`, `year` and the summer means
#'   `nitrate_summer`, `ammonium_summer`, `orthophosphate_summer`,
#'   `oxygen_saturation_summer` and `temperature_summer`, the summer minimum
#'   `oxygen_saturation_min` and maximum `temperature_max` (NA when too few
#'   measurements).
#' @export
compute_station_year_summer <- function(clean, code) {
  code <- code[c("nitrate", "ammonium", "orthophosphate", "oxygen_saturation", "temperature")]
  clean |>
    dplyr::filter(code_parametre %in% code, as.integer(format(date, "%m")) %in% 6:9) |>
    dplyr::mutate(
      parameter = names(code)[match(code_parametre, code)],
      year = as.integer(format(date, "%Y"))
    ) |>
    dplyr::summarise(
      mean = mean(value), min = min(value), max = max(value), n = dplyr::n(),
      .by = c(parameter, code_station, year)
    ) |>
    dplyr::filter(n >= 3) |>
    dplyr::transmute(
      code_station, year, parameter,
      summer = dplyr::if_else(parameter %in% c("oxygen_saturation", "temperature"), mean, log10(mean)),
      min = dplyr::if_else(parameter == "oxygen_saturation", min, NA),
      max = dplyr::if_else(parameter == "temperature", max, NA)
    ) |>
    tidyr::pivot_longer(c(summer, min, max), names_to = "statistic") |>
    dplyr::filter(!is.na(value)) |>
    tidyr::pivot_wider(names_from = c(parameter, statistic), values_from = value)
}

#' Status of the species in France from FishBase.
#'
#' @param species character, species names validated against FishBase.
#' @return tibble with `species`, `status` (FishBase status in France, e.g.
#'   "native", "introduced") and `non_native` (1 if introduced or not
#'   established, 0 otherwise).
#' @export
get_native_status <- function(species) {
  rfishbase::country(species) |>
    dplyr::filter(country == "France") |>
    dplyr::distinct(species = Species, status = Status) |>
    dplyr::mutate(non_native = as.numeric(status %in% c("introduced", "not established")))
}

#' Mean of a station-year variable over the year and the previous years.
#'
#' @param station_year tibble with `code_station`, `year` and `variable`,
#'   e.g. the `naiades_station_year` target.
#' @param variable character, column to average, e.g. "orthophosphate".
#' @param years integer, number of years averaged (the year and the
#'   `years - 1` previous ones); at least two of them must have a value.
#' @return tibble with `code_station`, `year` and the mean in
#'   `<variable>_<years>y`, for every year with at least two values in its
#'   window.
#' @export
average_previous_years <- function(station_year, variable, years) {
  value <- dplyr::select(station_year, code_station, year, value = dplyr::all_of(variable))
  value |>
    dplyr::distinct(code_station) |>
    dplyr::cross_join(tibble::tibble(lag = seq_len(years) - 1L)) |>
    dplyr::inner_join(value, by = "code_station", relationship = "many-to-many") |>
    dplyr::mutate(year = year + lag) |>
    dplyr::summarise(value = mean(value), n = dplyr::n(), .by = c(code_station, year)) |>
    dplyr::filter(n >= 2) |>
    dplyr::transmute(code_station, year, "{variable}_{years}y" := value)
}
