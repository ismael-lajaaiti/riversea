community_quality_targets <- list(
  tar_target(
    species_quality_response, # Water quality effect on presence (D16 additive models) and occupancy of each species.
    summarise_species_response(species_additive_model, quality_visit, river_presence)
  ),
  tar_target(
    species_trait, # Binary ecological traits of the fish species (freshwaterecology.info).
    define_fish_trait(fish_trait)
  ),
  tar_target(
    visit_trait, # Number of species with each trait at each visit with a water quality.
    count_visit_trait(quality_visit, river_presence, species_trait)
  ),
  tar_target(
    trait_name,
    c("intolerant", "tolerant", "rheophilic", "lithophilic", "benthic", "stenothermal")
  ),
  tar_target(
    trait_model, # Share of the species with each trait ~ dams + water quality.
    fit_trait_model(
      quality_visit, presence_mesh, dplyr::filter(visit_trait, trait == trait_name),
      paste(presence_term, "+ PC1")
    ),
    pattern = map(trait_name),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    dezerald_trait_file, # Fish traits of Dézerald et al. (72 species).
    "data/traits/dezerald_trait.csv",
    format = "file"
  ),
  tar_target(
    dezerald_metadata_file, # Description of the traits and their categories (Table S2).
    "data/traits/dezerald_metadata.xlsx",
    format = "file"
  ),
  tar_target(
    dezerald_trait,
    read_dezerald_trait(dezerald_trait_file)
  ),
  tar_target(
    dezerald_metadata,
    read_dezerald_metadata(dezerald_metadata_file)
  ),
  tar_target(
    pressure_name, # Water quality pressures modelled separately (NAIADES, log10 for nutrients).
    c("ammonium", "orthophosphate", "oxygen_saturation", "temperature")
  ),
  tar_target(
    naiades_station_year_summer, # Summer (June-September) means, oxygen minimum and temperature maximum.
    compute_station_year_summer(naiades_clean, naiades_parameter)
  ),
  tar_target(
    station_year_timing, # Each pressure averaged over the year and over the summer.
    dplyr::full_join(
      dplyr::select(naiades_station_year, code_station, year, dplyr::all_of(pressure_name)),
      dplyr::select(naiades_station_year_summer, code_station, year, dplyr::all_of(paste0(pressure_name, "_summer"))),
      by = c("code_station", "year")
    )
  ),
  tar_target(
    pressure_timing, # Averaging periods of each pressure (summer or year) and the matching column of `station_year_timing`.
    tidyr::expand_grid(
      pressure = pressure_name,
      duration = c("summer", "1y")
    ) |>
      dplyr::mutate(column = dplyr::if_else(duration == "1y", pressure, paste0(pressure, "_summer")))
  ),
  tar_target(
    pressure_visit, # Each pressure and duration predicted at the dam visits, as the water quality index (D14).
    predict_pressure_visit(station_year_timing, water_quality, water_quality_mesh, dam_visit, pressure_timing$column),
    pattern = map(pressure_timing),
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    pressure_visit_wide, # Dam visits with all pressures and durations and the local dam proximities, centred and scaled.
    prepare_pressure_visit(dam_visit, pressure_visit) |>
      add_local_dam(operation_dam)
  ),
  tar_target(
    species_pressure_model, # Presence models with dams (accumulation + local proximities) + the four annual pressures, for all species.
    fit_pressure_presence_model(
      pressure_visit_wide, river_presence, presence_mesh, dam_presence_species,
      paste(presence_term, "+ dam_below + dam_above +", paste(pressure_name, collapse = " + "))
    ),
    pattern = map(dam_presence_species),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    species_pressure_response, # Effects of dams and pressures on the presence of each species.
    summarise_pressure_response(species_pressure_model, c("accumulation", "dam_below", "dam_above", pressure_name))
  ),
  tar_target(
    native_status, # Status of the species in France (FishBase).
    get_native_status(dam_presence_species)
  ),
  tar_target(
    trait_predictor, # Traits of the a priori hypotheses (Dézerald et al.) and non-native status (FishBase).
    code_trait_predictor(dezerald_trait) |>
      dplyr::left_join(dplyr::select(native_status, species, non_native), by = "species")
  ),
  tar_target(
    trait_hypothesis, # A priori trait-pressure links and the expected sign of the trait effect on the response.
    tibble::tribble(
      ~term, ~trait, ~expected,
      "ammonium", "intolerant", -1,
      "ammonium", "omnivore", 1,
      "orthophosphate", "intolerant", -1,
      "orthophosphate", "omnivore", 1,
      "orthophosphate", "lithophilic", -1,
      "oxygen_saturation", "intolerant", 1,
      "oxygen_saturation", "oxygen_intolerant", 1,
      "temperature", "stenothermal", -1,
      "temperature", "rheophilic", -1,
      "accumulation", "diadromous", -1,
      "accumulation", "crossing", 1,
      "accumulation", "rheophilic", -1,
      "accumulation", "habitat_intolerant", -1,
      "dam_below", "rheophilic", -1,
      "dam_below", "lithophilic", -1,
      "dam_below", "crossing", 1,
      "dam_above", "rheophilic", -1,
      "dam_above", "lithophilic", -1,
      "accumulation", "non_native", 1,
      "dam_below", "non_native", 1,
      "dam_above", "non_native", 1
    )
  ),
  tar_target(
    trait_response, # Meta-regression of the species responses on the trait of each hypothesis.
    fit_trait_response(species_pressure_response, trait_predictor, trait_hypothesis$term, trait_hypothesis$trait),
    pattern = map(trait_hypothesis)
  ),
  tar_target(
    timing_setting, # Presence models with one pressure averaged over another duration, the others annual.
    dplyr::filter(pressure_timing, duration != "1y") |>
      dplyr::mutate(term = purrr::map2_chr(pressure, column, \(p, x) {
        paste(presence_term, "+ dam_below + dam_above +", paste(replace(pressure_name, pressure_name == p, x), collapse = " + "))
      }))
  ),
  tar_target(
    timing_model, # Presence models of each timing setting, for all species (fixed effects and WAIC only).
    {
      fit <- fit_pressure_presence_model(
        pressure_visit_wide, river_presence, presence_mesh, dam_presence_species, timing_setting$term
      )
      c(as.list(dplyr::select(timing_setting, pressure, duration, column)), fit[c("species", "summary.fixed", "waic")])
    },
    pattern = cross(timing_setting, dam_presence_species),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    timing_response, # Effect of each pressure on the presence of each species, per averaging duration.
    dplyr::bind_rows(
      dplyr::filter(species_pressure_response, term %in% pressure_name) |> dplyr::mutate(duration = "1y"),
      purrr::map_dfr(timing_model, \(m) {
        summarise_pressure_response(list(m), m$column) |> dplyr::mutate(term = m$pressure, duration = m$duration)
      })
    )
  ),
  tar_target(
    timing_hypothesis, # Trait hypotheses of the water quality pressures.
    dplyr::filter(trait_hypothesis, term %in% pressure_name)
  ),
  tar_target(
    timing_trait_response, # Trait meta-regressions of each water quality hypothesis, per averaging duration.
    purrr::map_dfr(unique(pressure_timing$duration), \(d) {
      fit_trait_response(dplyr::filter(timing_response, duration == d), trait_predictor, timing_hypothesis$term, timing_hypothesis$trait) |>
        dplyr::mutate(duration = d)
    }),
    pattern = map(timing_hypothesis)
  ),
  tar_target(
    timing_waic, # WAIC of the presence model of each species, per pressure and averaging duration (same visits).
    dplyr::bind_rows(
      tidyr::expand_grid(
        pressure = pressure_name, duration = "1y",
        purrr::map_dfr(species_pressure_model, \(m) tibble::tibble(species = m$species, waic = m$waic$waic))
      ),
      purrr::map_dfr(timing_model, \(m) {
        tibble::tibble(pressure = m$pressure, duration = m$duration, species = m$species, waic = m$waic$waic)
      })
    )
  )
)
