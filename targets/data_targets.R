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
  tar_target(
    diet_size,
    merge_diet_size(diet_no_rare, maturity_length)
  ),
  tar_target(
    diet_file,
    {
      path <- "data/diet/fishbase_sea.csv"
      readr::write_csv(diet_size, path)
      path
    },
    format = "file"
  ),
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
    size_extrema,
    get_size_extrema(sea_data_imputed$size)
  ),
  tar_target(
    predation_window,
    get_predation_window(diet_size)
  ),
  tar_target(
    diet_resource,
    read.csv(diet_resource_file)
  ),
  tar_target(n_class_vals, seq(3, 9)),
  tar_target(
    metaweb_table,
    tibble(
      num_classes = n_class_vals,
      metaweb = list(get_metaweb(
        sea_data_imputed$size,
        diet_size,
        diet_resource,
        predation_window,
        num_classes = n_class_vals
      ))
    ),
    pattern = map(n_class_vals),
  )
)
