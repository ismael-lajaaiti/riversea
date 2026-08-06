data_targets <- list(
  # Parameters.
  tar_target(
    params,
    list(
      occurence_min = 10, # Excluded.
      depth_max_rephy = 1, # Meter.
      year_min = 2005, # Included.
      year_max = 2022, # Included.
      distance_match_max = 10, # Kilometers.
      sampling_min = 20,
      river_station_network_max_dist = 5000, # Kilometers.
      n_size_class = 5, # Number of size classes for metaweb construction.
      amobio_average_window = "1y",
      district_kept = c("Loire", "Ardour-Garonne", "Seine", "Escaut-Somme")
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
  tar_target(
    sea_data_tidy,
    preprocess_sea_data(sea_data_raw, solper_reftax)
  ),
  tar_target(
    sea_data_no_rare,
    remove_rare_sea_catch(sea_data_tidy, params$occurence_min)
  ),
  tar_target(
    sea_data_imputed,
    infer_missing_size(sea_data_no_rare)
  ),
  tar_target(
    diet_validated_marine_file,
    "data/diet/marine_diet_validated.csv",
    format = "file"
  ),
  tar_target(
    diet_validated_river_file,
    "data/diet/river_diet_validated.csv",
    format = "file"
  ),
  tar_target(
    diet_missing_river_file,
    "data/diet/river_diet_missing.csv",
    format = "file"
  ),
  tar_target(
    diet_marine,
    get_diet_validated(diet_validated_marine_file)
  ),
  tar_target(
    diet_river,
    get_diet_validated_river(diet_validated_river_file)
  ),
  tar_target(
    diet_missing_river,
    get_diet_validated(diet_missing_river_file)
  ),
  tar_target(
    diet,
    merge_diet(merge_diet(diet_marine, diet_river), diet_missing_river)
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
    individual_fish_file,
    "data/river/output_individual_fish.rda",
    format = "file"
  ),
  tar_target(
    river_size,
    get_river_size(individual_fish_file)
  ),
  tar_target(
    river_size_no_rare,
    remove_rare_river_size(river_size, params$occurence_min)
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
  tar_target(
    hydrographic_basin_dir,
    download_hydrographic_basin(here("data", "geo", "hydrographic_basin")),
    format = "file"
  ),
  tar_target(
    hydrographic_basin,
    format_hydrographic_basin(hydrographic_basin_dir)
  ),
  tar_target(
    river_operation,
    get_river_operation(individual_fish_file, hydrographic_basin)
  ),
  tar_target(
    operation_id_check,
    assert_no_operation_id_collision(sea_data_tidy$trait, river_operation)
  ),
  tar_target(
    size_combined,
    {
      stopifnot(operation_id_check)
      merge_size(sea_data_imputed$size, river_size_no_rare)
    }
  ),
  tar_target(
    operation_combined,
    {
      stopifnot(operation_id_check)
      merge_operation(sea_data_tidy$trait, river_operation)
    }
  ),
  tar_target(
    operation_year_filtered,
    filter_year(operation_combined, params$year_min, params$year_max)
  ),
  tar_target(
    operation_location,
    classify_operation_location(operation_year_filtered, hydrographic_basin, params$district_kept)
  ),
  tar_target(
    size_year_filtered,
    size_combined |>
      dplyr::semi_join(operation_year_filtered, by = "operation_id")
  ),
  tar_target(
    species_traits,
    get_species_traits(sort(unique(size_year_filtered$species_valid)))
  ),
  tar_target(
    web_list,
    {
      stopifnot(diet_coverage_check)
      build_foodweb(
        size_year_filtered,
        diet,
        diet_resource,
        predation_window,
        num_classes = params$n_size_class,
        local = TRUE
      )
    }
  ),
  tar_target(
    metaweb_consistency_check,
    assert_metaweb_consistency(web_list$metaweb, resource_list)
  ),
  tar_target(
    foodweb_structure,
    {
      stopifnot(metaweb_consistency_check)
      measure_foodweb_structure(
        web_list,
        operation_year_filtered,
        resource_list
      )
    }
  ),
  tar_target(
    river_foodweb_structure,
    foodweb_structure |> dplyr::filter(survey == "river")
  ),
  tar_target(
    foodweb_size_info_combined,
    get_combined_foodweb_size_info(size_year_filtered, diet, operation_year_filtered)
  ),
  tar_target(dir_environment, "data/sea/raw/environment", format = "file"),
  tar_target(
    environment_data,
    read_environment_data(dir_environment)
  ),
  tar_target(
    fw_with_env,
    match_with_environment(
      foodweb_structure |> dplyr::filter(survey != "river"),
      environment_data
    )
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
    match_sea_district(
      foodweb_structure |> dplyr::filter(survey != "river"),
      hydrographic_area
    )
  ),
  tar_target(
    amobio_paths,
    download_amobio_data(here("data", "amobio"), verbose = TRUE),
    format = "file"
  ),
  tar_target(
    amobio_metrics,
    extract_metrics_amobio(amobio_paths, params$amobio_average_window)
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
      params$distance_match_max,
      params$district_kept
    )
  ),
  tar_target(
    river_network_full,
    extract_river_network(amobio_paths) |> patch_river_network()
  ),
  tar_target(
    river_network, # Restricted to the catchments of interest.
    restrict_river_network(
      river_network_full,
      river_mouth_district,
      params$river_station_network_max_dist * 1000
    )
  ),
  tar_target(
    river_network_geo,
    add_edge_geometry(river_network, amobio_paths)
  ),
  tar_target(
    river_station_snapped,
    snap_river_stations(river_operation, river_network_geo)
  ),
  tar_target(
    river_distance_to_mouth,
    compute_distance_to_mouth(
      river_network,
      river_station_snapped,
      river_mouth_district,
      params$district_kept
    )
  ),
  tar_target(
    combined_amobio_data,
    combine_amobio_data(amobio_metrics_dedup, amobio_nodes)
  ),
  tar_target(
    amobio_aspe_data,
    join_amobio_aspe(combined_amobio_data, river_foodweb_structure, river_operation)
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
    filter_year(
      from_date_to_year_month(rephy_data_surface),
      year_min = params$year_min,
      year_max = params$year_max
    )
  ),
  tar_target(
    rephy_data_sample,
    filter_sampling_rephy(rephy_data_recent, sampling_min = params$sampling_min)
  )
)
