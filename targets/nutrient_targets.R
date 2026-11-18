nutrient_targets <- list(
  tar_target(
    naiades_parameter, # Sandre codes of the NAIADES physico-chemistry parameters used.
    c(
      ammonium = "1335", nitrite = "1339", nitrate = "1340", orthophosphate = "1433",
      total_phosphorus = "1350", kjeldahl_nitrogen = "1319", temperature = "1301",
      oxygen = "1311", oxygen_saturation = "1312"
    )
  ),
  tar_target(
    naiades_period, # Sampling period downloaded, covering the fish visits (2005-2021).
    c("2000-01-01", "2021-12-31")
  ),
  tar_target(
    naiades_station, # Stations measuring the parameters in the districts studied.
    get_naiades_station(
      dplyr::filter(basin, district %in% params$district_kept),
      naiades_parameter, naiades_period
    )
  ),
  tar_target(
    naiades_station_batch, # Stations queried together.
    split(
      naiades_station$code_station,
      ceiling(seq_along(naiades_station$code_station) / 10)
    ),
    iteration = "list"
  ),
  tar_target(
    naiades_pc, # Raw NAIADES measurements of the parameters, via Hub'Eau.
    download_naiades_pc(naiades_station_batch, naiades_parameter, naiades_period),
    pattern = map(naiades_station_batch)
  ),
  tar_target(
    naiades_clean, # Measurements on water, in a single unit, below-limit values at half the limit, without outliers.
    clean_naiades_pc(naiades_pc)
  ),
  tar_target(
    naiades_station_year, # Station-year means of the parameters (log10 nutrients), all parameters present.
    compute_station_year_mean(naiades_clean, naiades_parameter)
  ),
  tar_target(
    water_quality_pca, # PCA of the station-year means, PC1 increasing with water quality.
    fit_water_quality_pca(naiades_station_year)
  ),
  tar_target(
    water_quality, # Water quality index (PC1) per station-year, with coordinates.
    compute_water_quality(naiades_station_year, water_quality_pca, naiades_station)
  ),
  tar_target(
    water_quality_mesh,
    build_water_quality_mesh(water_quality)
  ),
  tar_target(
    water_quality_model, # Model of PC1 (spatial field + AR1 year + iid station) fitted on all station-years.
    fit_water_quality_model(water_quality, water_quality_mesh, rep(FALSE, nrow(water_quality))),
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    water_quality_fold, # Random folds of stations and of station-years.
    {
      set.seed(1)
      assign_water_quality_fold(water_quality, 10)
    }
  ),
  tar_target(
    water_quality_cv_setting, # Cross-validation runs: fold type and held-out fold.
    tidyr::expand_grid(type = c("station", "station_year"), fold = 1:10)
  ),
  tar_target(
    water_quality_cv, # Predictions of each fold by the model fitted on the other folds.
    {
      held_out <- water_quality_fold[[water_quality_cv_setting$type]] == water_quality_cv_setting$fold
      fit_water_quality_model(water_quality, water_quality_mesh, held_out)$prediction |>
        dplyr::filter(held_out) |>
        dplyr::mutate(type = water_quality_cv_setting$type, fold = water_quality_cv_setting$fold)
    },
    pattern = map(water_quality_cv_setting),
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    water_quality_visit, # PC1 predicted at the dam visits (D12) whose site is a NAIADES station.
    predict_water_quality_visit(water_quality, water_quality_mesh, dam_visit),
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    quality_visit, # Dam visits (D12) with a predicted PC1, centred and scaled.
    prepare_quality_visit(dam_visit, water_quality_visit)
  ),
  tar_target(
    stressor_setting, # Presence models of dams and water quality, with and without spatial field.
    tidyr::expand_grid(
      tibble::tibble(
        model = c("dam", "quality", "additive", "interaction"),
        stressor = c("accumulation", "PC1", "accumulation + PC1", "accumulation * PC1")
      ),
      field = c(TRUE, FALSE),
      species = "Salmo trutta"
    ) |>
      dplyr::mutate(term = paste0(stressor, sub("accumulation", "", presence_term, fixed = TRUE)))
  ),
  tar_target(
    stressor_model,
    fit_stressor_presence_model(
      quality_visit, river_presence, if (stressor_setting$field) presence_mesh else NULL,
      stressor_setting$species, stressor_setting$term
    ),
    pattern = map(stressor_setting),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    species_stressor_model, # Interaction model (dams x PC1, with field), best for the trout, for all species.
    fit_stressor_presence_model(
      quality_visit, river_presence, presence_mesh, dam_presence_species,
      paste0("accumulation * PC1", sub("accumulation", "", presence_term, fixed = TRUE))
    ),
    pattern = map(dam_presence_species),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    fish_trait_file, # Fish traits of freshwaterecology.info (expert-based tolerance classes), via the aspe package.
    "data/traits/freshwaterecology_fish_traits.csv",
    format = "file"
  ),
  tar_target(
    fish_trait, # One row per species (the file has one row per ASPE code, e.g. for carp varieties).
    readr::read_csv(fish_trait_file, comment = "#", show_col_types = FALSE) |>
      dplyr::select(-esp_code_alternatif, -esp_nom_commun) |>
      dplyr::distinct(esp_nom_latin, .keep_all = TRUE)
  )
)
