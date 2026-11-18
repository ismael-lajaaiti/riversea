foodweb_nutrient_targets <- list(
  tar_target(
    quality_foodweb, # Webs of the dam visits (D13) with a predicted water quality (D14), and their PCA scores.
    prepare_quality_foodweb(dam_foodweb, water_quality_visit, structure_pca)
  ),
  tar_target(
    quality_structure_setting, # Metrics modelled with dams only or dams + water quality, without and with richness.
    tidyr::expand_grid(
      tibble::tibble(
        metric = c(
          "fish_richness", "mean_trophic_level", "fish_diet_overlap", "trophic_redundancy",
          "log_biomass", "biomass_trophic_level", "biomass_pyramid", "PC1", "PC2"
        ),
        family = c("poisson", rep("gaussian", 8))
      ),
      model = c("dam", "additive"),
      richness = c(FALSE, TRUE)
    ) |>
      dplyr::filter(!(metric == "fish_richness" & richness)) |>
      dplyr::mutate(term = paste0(
        presence_term, ifelse(model == "additive", " + quality", ""),
        ifelse(richness, " + log_richness", "")
      ))
  ),
  tar_target(
    quality_structure_model,
    fit_stressor_structure_model(
      quality_foodweb, presence_mesh, quality_structure_setting$metric,
      quality_structure_setting$family, quality_structure_setting$term
    ),
    pattern = map(quality_structure_setting),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    species_additive_model, # Presence models with dams + water quality (D14 visits, with field), for all species.
    fit_stressor_presence_model(
      quality_visit, river_presence, presence_mesh, dam_presence_species,
      paste(presence_term, "+ PC1")
    ),
    pattern = map(dam_presence_species),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    mediation_stressor, # Stressor terms of the presence models: dams and water quality.
    c("accumulation", "PC1")
  ),
  tar_target(
    stressor_decomposition, # Structure change expected from the change of species presence, per stressor.
    {
      set.seed(1)
      decompose_dam_structure(
        species_additive_model, species_removal_effect, mediation_stressor, n_draws = 2000
      )
    },
    pattern = map(mediation_stressor),
    iteration = "list"
  ),
  tar_target(
    stressor_permutation, # Same with the stressor effects shuffled among species.
    {
      set.seed(1)
      permute_dam_structure(
        species_additive_model, species_removal_effect, mediation_stressor, n_permutations = 2000
      )
    },
    pattern = map(mediation_stressor),
    iteration = "list"
  )
)
