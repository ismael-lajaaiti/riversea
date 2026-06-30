data_targets <- list(
  # Parameters.
  tar_target(
    params,
    list(
      occurence_min = 10
    )
  ),
  # Workshop river data.
  tar_target(
    workshop_zip,
    download_workshop_data(workshop_dir),
    format = "file"
  ),
  tar_target(
    solper_reftax,
    here(sea_data_raw, "Reftax_SIH.txt"),
    format = "file"
  ),
  tar_target(
    workshop_unzipped,
    unzip_workshop(workshop_zip, workshop_dir),
    format = "file"
  ),
  # Sea survey data.
  tar_target(
    sea_data_tidy,
    preprocess_sea_data(sea_data_raw, solper_reftax)
  ),
  # Infer missing size.
  tar_target(
    sea_data_imputed,
    infer_missing_size(sea_data_tidy)
  ),
  # Extract fish diet.
  tar_target(
    species_list,
    get_species_list(sea_data_imputed)
  ),
  tar_target(
    diet,
    get_diet_category(species_list)
  ),
  tar_target(
    diet_wide,
    widen_diet_category(diet)
  ),
  tar_target(
    diet_larvae,
    add_larvae(diet_wide),
  ),
  tar_target(
    list_no_rare,
    filter_out_rare(
      diet_larvae,
      sea_data_imputed,
      occurence_min = params$occurence_min
    )
  ),
  tar_target(diet_no_rare, list_no_rare$diet),
  tar_target(n_rare, list_no_rare$n_removed),
  # Life stage lengths.
  tar_target(
    species_diet,
    diet_no_rare |> pull(species) |> unique()
  ),
  tar_target(
    maturity_length,
    get_maturity_length(species_diet)
  ),
  tar_target(
    life_stage_length_file,
    "data/sea/raw/lifestage_size.csv",
    format = "file"
  ),
  tar_target(
    life_stage_length,
    get_life_stage_length(life_stage_length_file, maturity_length)
  ),
  tar_target(
    diet_and_length,
    merge_diet_length(diet_no_rare, life_stage_length)
  ),
  tar_target(
    diet_file,
    {
      path <- "data/diet/fishbase_sea.csv"
      readr::write_csv(diet_and_length, path)
      path
    },
    format = "file"
  ),
  tar_target(
    size_extrema,
    get_size_extrema(sea_data_imputed$size)
  ),
  tar_target(
    foodweb_size_info,
    get_foodweb_size_info(sea_data_imputed$size, diet_and_length)
  ),
  tar_target(
    predation_window,
    get_predation_window(diet_and_length)
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
        diet_and_length,
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
    get_metaweb(
      sea_data_imputed$size,
      diet_and_length,
      diet_resource,
      predation_window,
      num_classes = num_classes,
      local = TRUE
    )
  ),
  tar_target(
    local_foodwebs,
    prepare_local_foodwebs(web_list, sea_data_tidy, resource_list)
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
  tar_target(
    amobio_metrics,
    extract_metrics_amobio(amobio_paths)
  ),
  tar_target(
    amobio_nodes,
    extract_nodes_amobio(amobio_paths)
  ),
  tar_target(
    combined_amobio_data,
    combine_amobio_data(amobio_metrics, amobio_nodes)
  )
)
