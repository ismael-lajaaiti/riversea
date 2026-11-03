presence_targets <- list(
  tar_target(
    presence_mesh,
    build_presence_mesh(dam_visit)
  ),
  tar_target(
    presence_term, # Dam and position terms of the presence models (D12).
    "accumulation + log_width + log_alt + I(log_alt^2) + log_dist + log_slope"
  ),
  tar_target(
    presence_model,
    fit_presence_model(
      dam_visit, river_presence, presence_mesh, dam_presence_species, presence_term
    ),
    pattern = map(dam_presence_species),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    presence_model_nofield, # Same models without the spatial field.
    fit_presence_model(
      dam_visit, river_presence, NULL, dam_presence_species, presence_term
    ),
    pattern = map(dam_presence_species),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  )
)
