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
    water_quality_fold, # Spatial-block folds (50 km blocks) of the cross-validation.
    {
      set.seed(1)
      assign_block_fold(water_quality, 10)
    }
  ),
  tar_target(
    water_quality_cv, # Predictions of each fold by the model fitted on the other folds.
    fit_water_quality_model(water_quality, water_quality_mesh, water_quality_fold == water_quality_fold_index)$prediction |>
      dplyr::filter(held_out) |>
      dplyr::mutate(fold = water_quality_fold_index),
    pattern = map(water_quality_fold_index),
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    water_quality_fold_index, # Folds of the cross-validation.
    1:10
  )
)
