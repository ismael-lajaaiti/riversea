data_targets <- list(
  # Parameters.
  tar_target(
    params,
    list(
      occurence_min = 10,
      depth_max_rephy = 1, # Meter.
      year_min = 2000,
      distance_match_max = 10, # Kilometers.
      sampling_min = 20,
      river_network_max_dist = 100 # Kilometers.
    )
  ),
  tar_target(
    sea_data_raw,
    download_sea_raw_data(),
    format = "file"
  ),
  tar_target(
    solper_reftax,
    here(sea_data_raw, "Reftax_SIH.txt"),
    format = "file"
  ),
  # Sea survey data.
  tar_target(
    sea_data_tidy,
    preprocess_sea_data(sea_data_raw, solper_reftax)
  ),
  tar_target(
    sea_data_no_rare,
    remove_rare_sea_catch(sea_data_tidy, params$occurence_min)
  ),
  # Infer missing size.
  tar_target(
    sea_data_imputed,
    infer_missing_size(sea_data_no_rare)
  ),
  # Fish diet, validated against the literature.
  tar_target(
    diet_validated_file,
    "data/diet/marine_diet_validated.csv",
    format = "file"
  ),
  tar_target(
    diet_marine,
    get_diet_validated(diet_validated_file)
  ),
  # Fish diet for freshwater fish, validated against the literature.
  tar_target(
    diet_validated_river_file,
    "data/diet/river_diet_validated.csv",
    format = "file"
  ),
  tar_target(
    diet_river,
    get_diet_validated_river(diet_validated_river_file)
  ),
  # Merge freshwater into marine/estuarine; freshwater wins on overlap.
  tar_target(
    diet,
    merge_diet_river_marine(diet_marine, diet_river)
  ),
  tar_target(
    diet_coverage_check,
    assert_diet_coverage(
      diet,
      sea_data_tidy$catch,
      params$occurence_min
    )
  ),
  tar_target(
    size_extrema,
    get_size_extrema(sea_data_imputed$size)
  ),
  tar_target(
    foodweb_size_info,
    get_foodweb_size_info(sea_data_imputed$size, diet)
  ),
  tar_target(
    predation_window,
    get_predation_window(diet)
  ),
  tar_target(
    diet_resource_file,
    "data/diet/resource_diet.csv",
    format = "file"
  ),
  tar_target(
    diet_resource,
    read.csv(diet_resource_file),
  ),
  tar_target(
    resource_list,
    setdiff(
      names(diet_resource),
      c("light", "fish", "reference", "species")
    )
  ),
  tar_target(n_class_vals, seq(3, 9)),
  tar_target(
    metaweb_table,
    tibble(
      num_classes = n_class_vals,
      metaweb = list(get_metaweb(
        sea_data_imputed$size,
        diet,
        diet_resource,
        predation_window,
        num_classes = n_class_vals
      )$metaweb)
    ),
    pattern = map(n_class_vals),
  ),
  tar_target(num_classes, 5),
  tar_target(
    web_list,
    {
      stopifnot(diet_coverage_check)
      get_metaweb(
        sea_data_imputed$size,
        diet,
        diet_resource,
        predation_window,
        num_classes = num_classes,
        local = TRUE
      )
    }
  ),
  tar_target(
    metaweb_consistency_check,
    assert_metaweb_consistency(web_list$metaweb, resource_list)
  ),
  tar_target(
    local_foodwebs,
    {
      stopifnot(metaweb_consistency_check)
      prepare_local_foodwebs(web_list, sea_data_tidy, resource_list)
    }
  ),
  tar_target(dir_environment, "data/sea/raw/environment", format = "file"),
  tar_target(
    environment_data,
    read_environment_data(dir_environment)
  ),
  tar_target(
    fw_with_env,
    match_with_environment(local_foodwebs, environment_data)
  ),
  tar_target(
    hydrographic_area_dir,
    download_hydrographic_files(here("data", "geo", "hydrographic")),
    format = "file"
  ),
  tar_target(
    hydrographic_area,
    format_hydrographic_area(hydrographic_area_dir)
  ),
  tar_target(
    fw_with_district,
    match_foodweb_district(local_foodwebs, hydrographic_area)
  ),
  tar_target(
    amobio_paths,
    download_amobio_data(here("data", "amobio"), verbose = TRUE),
    format = "file"
  ),
  tar_target(average_time, "1y"),
  tar_target(
    amobio_metrics,
    extract_metrics_amobio(amobio_paths, average_time)
  ),
  tar_target(
    amobio_metrics_dedup,
    dedup_amobio_metrics(amobio_metrics)
  ),
  tar_target(
    amobio_nodes,
    extract_nodes_amobio(amobio_paths)
  ),
  tar_target(
    amobio_river_mouth,
    extract_river_mouth_amobio(amobio_paths)
  ),
  tar_target(
    river_mouth_district,
    match_mouth_district(
      amobio_river_mouth,
      hydrographic_area,
      params$distance_match_max
    )
  ),
  tar_target(
    amobio_network,
    extract_amobio_network(amobio_paths)
  ),
  tar_target(
    amobio_network_restricted,
    restrict_amobio_network(
      amobio_network,
      river_mouth_district,
      params$river_network_max_dist * 1000
    )
  ),
  tar_target(
    combined_amobio_data,
    combine_amobio_data(amobio_metrics_dedup, amobio_nodes)
  ),
  tar_target(
    aspe_file_foodweb,
    here::here("data", "river", "output_size2webs.rda"),
    format = "file"
  ),
  tar_target(
    aspe_file_code,
    here::here("data", "river", "output_individual_fish.rda"),
    format = "file"
  ),
  tar_target(aspe_data, get_aspe_data(aspe_file_foodweb, aspe_file_code)),
  tar_target(
    amobio_aspe_data,
    join_amobio_aspe(combined_amobio_data, aspe_data)
  ),
  tar_target(rephy_dir, download_rephy(), format = "file"),
  tar_target(
    rephy_files,
    list.files(rephy_dir, pattern = "^REPHY.*\\.csv$", full.names = TRUE),
    format = "file"
  ),
  tar_target(
    rephy_file_atlantic,
    stringr::str_subset(rephy_files, "Atlantique")
  ),
  tar_target(
    rephy_raw_data,
    purrr::map_dfr(rephy_file_atlantic, read_rephy)
  ),
  tar_target(
    rephy_data,
    extract_nutrients_rephy(rephy_raw_data)
  ),
  tar_target(
    rephy_district_data,
    match_rephy_district(rephy_data, hydrographic_area)
  ),
  tar_target(
    rephy_yearly,
    summarise_rephy_yearly(rephy_district_data)
  ),
  tar_target(
    rephy_data_surface,
    filter_surface_rephy(rephy_data, depth_max = params$depth_max_rephy),
  ),
  tar_target(
    rephy_data_recent,
    filter_year_rephy(rephy_data_surface, year_min = params$year_min)
  ),
  tar_target(
    rephy_data_sample,
    filter_sampling_rephy(rephy_data_recent, sampling_min = params$sampling_min)
  )
)
