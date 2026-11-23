functional_diversity_targets <- list(
  tar_target(
    functional_variable, # Effect traits used for functional diversity, their group and type.
    dplyr::bind_rows(
      tibble::tribble(
        ~variable, ~group, ~type,
        "diet", "Diet", "fuzzy",
        "hab.RHEO", "Habitat", "ordinal",
        "hab.TROPH", "Habitat", "nominal",
        "living", "Habitat", "fuzzy",
        "morp.LMAX", "Size and life history", "ordinal",
        "life.LONG", "Size and life history", "ordinal",
        "life.MAT", "Size and life history", "ordinal",
        "rep.ABSFEC", "Size and life history", "ordinal",
        "rep.REPMODE", "Reproduction", "nominal",
        "rep.REPHAB", "Reproduction", "nominal",
        "rep.CARE", "Reproduction", "nominal",
        "rep.SPAWNMODE", "Reproduction", "nominal"
      ),
      tibble::tibble(
        variable = paste0("morp.", c("BLBD", "EHBD", "EHHD", "MOBD", "JLHD", "HDBD", "PFIBD", "PFLBD", "CFDCPD")),
        group = "Morphology", type = "ordinal"
      )
    )
  ),
  tar_target(
    functional_trait, # Coded effect traits of the river species.
    code_functional_trait(dplyr::filter(dezerald_trait, species %in% river_presence$species_valid))
  ),
  tar_target(
    functional_distance, # Functional distance between species, equal weight per trait group.
    compute_functional_distance(functional_trait, functional_variable)
  ),
  tar_target(
    functional_space, # PCoA of the functional distance and its quality (mean absolute deviation, Maire et al. 2015).
    mFD::quality.fspaces(
      sp_dist = stats::as.dist(functional_distance$distance), maxdim_pcoa = 8,
      deviation_weighting = c("absolute", "squared"), fdist_scaling = FALSE
    )[c("quality_fspaces", "details_fspaces")]
  ),
  tar_target(
    functional_axis_number, # Axes kept: mean absolute deviation 0.06, eigenvalues above the broken stick.
    4
  ),
  tar_target(
    functional_coord,
    functional_space$details_fspaces$sp_pc_coord[, seq_len(functional_axis_number)]
  ),
  tar_target(
    functional_presence, # Species caught at each dam visit.
    river_presence |>
      dplyr::filter(operation_id %in% dam_visit$operation_id) |>
      dplyr::distinct(operation_id, species = species_valid)
  ),
  tar_target(
    functional_diversity, # Functional richness and dispersion of each visit.
    compute_functional_diversity(functional_coord, functional_presence)
  ),
  tar_target(
    functional_replicate,
    seq_len(199)
  ),
  tar_target(
    functional_null, # Functional diversity with species coordinates shuffled among species.
    simulate_functional_diversity(functional_coord, functional_presence, functional_replicate),
    pattern = map(functional_replicate)
  ),
  tar_target(
    functional_ses, # Standardised effect sizes of functional richness and dispersion.
    compute_functional_ses(functional_diversity, functional_null)
  ),
  tar_target(
    functional_response,
    c("fric_ses", "fdis_ses", "fric", "fdis")
  ),
  tar_target(
    functional_model, # Functional diversity ~ dams + ammonium + temperature (phosphorus and oxygen left out: same pollution syndrome as ammonium).
    fit_functional_model(
      pressure_visit_wide, functional_ses, presence_mesh, functional_response,
      paste(presence_term, "+ dam_below + dam_above + ammonium + temperature")
    ),
    pattern = map(functional_response),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    functional_distinctiveness, # Distinctiveness (mean distance to the others) and uniqueness of each species.
    compute_functional_distinctiveness(functional_distance)
  ),
  tar_target(
    distinctiveness_response, # Meta-regression of the response of species to each pressure on their distinctiveness.
    fit_trait_response(
      species_pressure_response, functional_distinctiveness, distinctiveness_setting$term, distinctiveness_setting$trait
    ),
    pattern = map(distinctiveness_setting)
  ),
  tar_target(
    distinctiveness_setting,
    tidyr::expand_grid(
      term = c("ammonium", "orthophosphate", "oxygen_saturation", "temperature", "accumulation", "dam_below", "dam_above"),
      trait = c("distinctiveness_scaled", "uniqueness_scaled")
    )
  ),
  tar_target(
    pairwise_diversity, # Mean pairwise and nearest-neighbour distances, overall and per trait group.
    compute_pairwise_diversity(functional_distance, functional_presence)
  ),
  tar_target(
    pairwise_null, # Pairwise diversity with species shuffled.
    simulate_pairwise_diversity(functional_distance, functional_presence, functional_replicate),
    pattern = map(functional_replicate)
  ),
  tar_target(
    pairwise_ses,
    compute_pairwise_ses(pairwise_diversity, pairwise_null)
  ),
  tar_target(
    pairwise_response, # Indices and groups modelled.
    setdiff(names(pairwise_ses), "operation_id")
  ),
  tar_target(
    pairwise_model, # Pairwise diversity SES ~ dams + ammonium + temperature.
    fit_functional_model(
      pressure_visit_wide, pairwise_ses, presence_mesh, pairwise_response,
      paste(presence_term, "+ dam_below + dam_above + ammonium + temperature")
    ),
    pattern = map(pairwise_response),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    community_trait, # Size, life history and diet traits of the species.
    define_community_trait(functional_trait)
  ),
  tar_target(
    community_trait_mean, # Mean trait of the species present at each visit (scaled).
    compute_community_trait_mean(community_trait, functional_presence)
  ),
  tar_target(
    community_trait_name,
    setdiff(names(community_trait), "species")
  ),
  tar_target(
    community_trait_model, # Community mean trait ~ dams + ammonium + temperature.
    fit_functional_model(
      pressure_visit_wide, community_trait_mean, presence_mesh, community_trait_name,
      paste(presence_term, "+ dam_below + dam_above + ammonium + temperature")
    ),
    pattern = map(community_trait_name),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    functional_redundancy, # Number of species with a near-copy (5% closest pairs of the pool) at each visit.
    compute_functional_redundancy(functional_distance, functional_presence, 0.05) |>
      dplyr::mutate(log_richness = as.numeric(scale(log(trials))))
  ),
  tar_target(
    redundancy_setting, # Pressures in thirds and the other pressure controlled (phosphorus and oxygen left out: same pollution syndrome as ammonium).
    tibble::tibble(pressure = c("ammonium", "temperature"), covariate = c("temperature", "ammonium"))
  ),
  tar_target(
    redundancy_model, # Redundant species (out of the species present) ~ dams + other pressure + pressure in thirds.
    c(
      list(pressure = redundancy_setting$pressure),
      fit_functional_model(
        add_pressure_class(pressure_visit_wide, redundancy_setting$pressure, 3), functional_redundancy, presence_mesh,
        "redundant",
        paste(presence_term, "+ dam_below + dam_above +", redundancy_setting$covariate, "+ class_2 + class_3 + log_richness"),
        family = "binomial",
        lincomb = c( # Deviation of each third from the average of the thirds (logit).
          INLA::inla.make.lincomb(class_2 = -1, class_3 = -1) |> stats::setNames("low"),
          INLA::inla.make.lincomb(class_2 = 1) |> stats::setNames("medium"),
          INLA::inla.make.lincomb(class_3 = 1) |> stats::setNames("high")
        )
      )
    ),
    pattern = map(redundancy_setting),
    iteration = "list"
  ),
  tar_target(
    similarity_setting, # Pressures whose classes are compared, and the covariates controlled (river position, dams, other pressures; for ammonium, not phosphorus and oxygen).
    tibble::tibble(
      pressure = c("ammonium", "temperature"),
      covariate = list(
        c("log_dist", "log_width", "log_alt", "log_slope", "accumulation", "temperature"),
        c("log_dist", "log_width", "log_alt", "log_slope", "accumulation", "ammonium", "orthophosphate", "oxygen_saturation")
      )
    )
  ),
  tar_target(
    site_class, # Most recent visit of each site, its community and its third of the pressure.
    select_site_class(pressure_visit_wide, functional_presence, similarity_setting$pressure, 3),
    pattern = map(similarity_setting),
    iteration = "list"
  ),
  tar_target(
    class_dissimilarity, # Dissimilarity between sites of the same third of the pressure, less than 200 km apart.
    dplyr::mutate(
      compute_class_dissimilarity(site_class, similarity_setting$covariate[[1]], 200),
      pressure = similarity_setting$pressure
    ),
    pattern = map(site_class, similarity_setting),
    iteration = "list"
  ),
  tar_target(
    dissimilarity_replicate,
    0:199
  ),
  tar_target(
    class_dissimilarity_model, # Adjusted distance-decay per class, bootstrap over sites.
    dplyr::mutate(
      fit_class_dissimilarity_model(dplyr::select(class_dissimilarity, -pressure), dissimilarity_replicate),
      pressure = class_dissimilarity$pressure[1]
    ),
    pattern = cross(class_dissimilarity, dissimilarity_replicate)
  ),
  tar_target(
    class_dispersion, # Dispersion of communities around their class centroid, position and covariates partialled out.
    dplyr::mutate(
      compute_class_dispersion(site_class, similarity_setting$covariate[[1]], 199),
      pressure = similarity_setting$pressure
    ),
    pattern = map(site_class, similarity_setting)
  )
)
