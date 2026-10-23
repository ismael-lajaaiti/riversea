dam_web_targets <- list(
  tar_target(
    dam_web_settings,
    list(
      min_sites = 30, # Species modelled: present at >= this number of sites.
      n_draws = 200 # Joint posterior draws per presence model.
    )
  ),
  tar_target(
    river_presence,
    size_year_filtered |>
      dplyr::filter(survey == "river") |>
      dplyr::distinct(operation_id, species_valid)
  ),
  tar_target(
    river_operation_mouth,
    inland_distance_to_mouth |>
      dplyr::semi_join(river_presence, by = "operation_id") |>
      assign_river_mouth(
        river_network,
        operation = _,
        river_mouth = sf::st_drop_geometry(river_mouth_district)
      )
  ),
  tar_target(
    dam_visit,
    prepare_dam_visits(
      continuum_operation, continuum_variables_operation, river_operation,
      river_effort, river_operation_mouth, operation_dam, large_dam,
      river_presence
    )
  ),
  tar_target(
    dam_presence_species,
    select_presence_species(river_presence, dam_visit, dam_web_settings$min_sites)
  ),
  tar_target(
    dam_presence_sample,
    sample_presence_dam(
      dam_visit, river_presence, dam_presence_species, dam_web_settings$n_draws
    ),
    pattern = map(dam_presence_species),
    iteration = "list"
  ),
  tar_target(
    dam_presence_probability,
    mean_presence_probability(dam_presence_sample)
  ),
  tar_target(
    dam_presence_effect, # Posterior of the dam effect of each species.
    purrr::map_dfr(dam_presence_sample, \(x) tibble::tibble(
      species = x$species, mean = mean(x$beta),
      lower = stats::quantile(x$beta, 0.025), upper = stats::quantile(x$beta, 0.975)
    ))
  ),
  tar_target(
    species_metaweb,
    aggregate_species_web(web_list$metaweb)
  ),
  tar_target(
    dam_web_simulated,
    simulate_dam_webs(dam_presence_sample, dam_visit, species_metaweb, resource_list)
  ),
  tar_target(
    dam_web_counterfactual,
    summarise_dam_counterfactual(dam_web_simulated, dam_visit)
  ),
  tar_target(
    dam_realisation, # Realisations of the generated webs refitted on the dams.
    seq_len(50)
  ),
  tar_target(
    dam_structure_generated,
    fit_structure_dam(
      dam_visit, extract_generated_webs(dam_web_simulated, dam_realisation), 200
    ) |>
      dplyr::mutate(realisation = dam_realisation),
    pattern = map(dam_realisation)
  ),
  tar_target(
    dam_structure_observed,
    fit_structure_dam(dam_visit, dam_web_observed, 200)
  ),
  tar_target(
    dam_weight, # Rank weights of the dams compared for the web structure.
    c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2)
  ),
  tar_target(
    dam_weight_structure,
    weight_structure_dams(dam_visit, dam_web_observed, dam_weight),
    pattern = map(dam_weight)
  ),
  tar_target(
    dam_weight_presence_grid, # Coarser grid for the 48 presence models.
    c(0, 0.5, 1, 1.5, 2)
  ),
  tar_target(
    dam_weight_presence,
    weight_presence_dams(dam_visit, river_presence, dam_presence_species, dam_weight_presence_grid),
    pattern = cross(dam_presence_species, dam_weight_presence_grid)
  ),
  tar_target(
    dam_web_one_less,
    simulate_one_dam_less(dam_presence_sample, dam_visit, species_metaweb, resource_list)
  ),
  tar_target(
    dam_removal,
    summarise_dam_removal(dam_web_simulated, dam_web_one_less, operation_dam, dam_visit)
  ),
  tar_target(
    dam_web_observed_size_class, # Size-class metrics (D11, D12), species-level richness.
    dam_web_observed |>
      dplyr::select(operation_id, richness) |>
      dplyr::inner_join(
        foodweb_structure |>
          dplyr::filter(survey == "river") |>
          dplyr::distinct(operation_id, .keep_all = TRUE) |>
          dplyr::transmute(
            operation_id = as.character(operation_id), trophic_length,
            piscivore = as.numeric(frac_piscivorous > 0), fish_diet_overlap
          ),
        by = "operation_id"
      )
  ),
  tar_target(
    dam_structure_size_class,
    fit_structure_dam(dam_visit, dam_web_observed_size_class, 200)
  ),
  tar_target(
    top_trophic_species,
    find_top_trophic_species(
      web_list$local, river_presence, dam_visit$operation_id,
      dam_presence_species, species_metaweb, resource_list
    )
  ),
  tar_target(
    network_reach,
    prepare_network_reaches(
      rht, basin, river_network_geo, river_network, river_mouth_district,
      dam_visit, continuum_variables_operation
    )
  ),
  tar_target(
    network_reach_dam, # Dams downstream of each reach.
    locate_operations_dams(
      dam_network$network,
      network_reach |> dplyr::transmute(operation_id = as.character(reach_id), node_id),
      large_dam,
      river_mouth_district
    ) |>
      dplyr::filter(position == "upstream")
  ),
  tar_target(
    dam_presence_weight_sample,
    sample_presence_weight(
      dam_visit, river_presence, dam_presence_species, dam_weight_presence_grid,
      dam_web_settings$n_draws
    ),
    pattern = cross(dam_presence_species, dam_weight_presence_grid),
    iteration = "list"
  ),
  tar_target(
    dam_removal_network,
    predict_dam_removal_network(dam_presence_weight_sample, network_reach, network_reach_dam, dam_visit)
  ),
  tar_target(
    dam_traits,
    describe_dams(
      large_dam, dam_dam, network_reach, network_reach_dam, rht, dam_network$network,
      river_mouth_district, basin
    )
  ),
  tar_target(
    dam_connectivity, # Gain in the dendritic connectivity index (Cote et al. 2009).
    compute_dam_connectivity(network_reach, network_reach_dam, c(0.25, 0.5, 0.75))
  ),
  tar_target(
    dam_web_observed,
    measure_observed_species_webs(
      river_presence, dam_visit$operation_id, dam_presence_species,
      species_metaweb, resource_list
    )
  )
)
