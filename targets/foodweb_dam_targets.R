foodweb_dam_targets <- list(
  tar_target(
    node_biomass, # Individual weights of the ASPE file, summed by size class.
    compute_node_biomass(river_size_no_rare, web_list$size_class)
  ),
  tar_target(
    dam_foodweb, # Observed size-class webs of the dam visits (D12) and their structure.
    measure_dam_webs(river_foodweb_structure, dam_visit, node_biomass, resource_list)
  ),
  tar_target(
    species_removal, # Structure of each web with each fish species removed in turn.
    remove_web_species(dam_foodweb, resource_list)
  ),
  tar_target(
    species_removal_effect, # Standardised change of the first PCA axis per removed species.
    {
      set.seed(1)
      with_richness <- which(
        structure_axis_setting$axis == "PC1" & structure_axis_setting$richness
      )
      summarise_species_removal(
        compute_removal_pc1(
          species_removal, structure_pca,
          structure_axis_model[[with_richness]]$summary.fixed["log_richness", "mean"]
        ),
        min_sites = 30
      )
    }
  ),
  tar_target(
    structure_setting, # Metrics modelled, their likelihood and controls.
    tidyr::expand_grid(
      tibble::tibble(
        metric = c(
          "fish_richness", "trophic_length", "mean_trophic_level",
          "frac_piscivorous", "fish_diet_overlap", "trophic_redundancy",
          "log_biomass", "biomass_trophic_level", "piscivore_biomass",
          "biomass_pyramid"
        ),
        family = c("poisson", rep("gaussian", 9))
      ),
      richness = c(FALSE, TRUE)
    ) |>
      dplyr::filter(!(metric == "fish_richness" & richness)) |>
      dplyr::mutate(term = paste0(presence_term, ifelse(richness, " + log_richness", "")))
  ),
  tar_target(
    structure_model,
    fit_structure_model(
      dam_foodweb, presence_mesh, structure_setting$metric, structure_setting$family,
      structure_setting$term
    ),
    pattern = map(structure_setting),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    structure_pca, # PCA of the web metrics kept in D13, over the visits where all are defined.
    fit_structure_pca(
      dam_foodweb,
      c(
        "fish_richness", "mean_trophic_level", "fish_diet_overlap", "trophic_redundancy",
        "log_biomass", "biomass_trophic_level", "biomass_pyramid"
      )
    )
  ),
  tar_target(
    species_presence, # Presence of the fish species in the webs of the visits.
    compute_species_presence(dam_foodweb, resource_list)
  ),
  tar_target(
    structure_axis_setting, # PCA axes modelled, without and with richness.
    tidyr::expand_grid(axis = c("PC1", "PC2"), richness = c(FALSE, TRUE))
  ),
  tar_target(
    structure_axis_model, # Dam models of the PCA axes (scores scaled to unit variance).
    dam_foodweb |>
      dplyr::mutate(
        score = stats::predict(structure_pca, dam_foodweb)[, structure_axis_setting$axis],
        score = score / stats::sd(score, na.rm = TRUE)
      ) |>
      fit_structure_model(
        presence_mesh, "score", "gaussian",
        paste0(presence_term, ifelse(structure_axis_setting$richness, " + log_richness", ""))
      ),
    pattern = map(structure_axis_setting),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    structure_decomposition, # Structure change expected from the species lost to dams.
    {
      set.seed(1)
      decompose_dam_structure(presence_model, species_removal_effect, n_draws = 2000)
    }
  ),
  tar_target(
    structure_permutation, # Same with the dam effects shuffled among species.
    {
      set.seed(1)
      permute_dam_structure(presence_model, species_removal_effect, n_permutations = 2000)
    }
  )
)
