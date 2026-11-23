sea_community_targets <- list(
  tar_target(
    pomet_file, # POMET catches, with gear, water body and estuarine guilds.
    here("data", "sea", "raw", "poisson_pomet.csv"),
    format = "file"
  ),
  tar_target(
    pomet_operation, # Gear, water body and depth of each POMET operation.
    read_pomet_operation(pomet_file)
  ),
  tar_target(
    estuary_guild, # Estuarine use guild of each species (POMET reference).
    read_estuary_guild(pomet_file)
  ),
  tar_target(
    sea_environment_dir,
    here("data", "sea", "raw", "environment"),
    format = "file"
  ),
  tar_target(
    sea_environment, # Yearly gridded environment of the Bay of Biscay coast (salinity, temperature, chlorophyll, discharge anomaly).
    read_sea_environment(sea_environment_dir)
  ),
  tar_target(
    sea_operation, # Estuary and coastal operations of the four districts (at sea and inland), with signed distance to the mouth (negative upstream), gear, water body and salinity (measured, else environment grid).
    operation_location |>
      dplyr::filter(survey %in% c("pomet", "nurse"), district %in% c("Loire", "Ardour-Garonne", "Seine", "Escaut-Somme")) |>
      dplyr::left_join(
        dplyr::bind_rows(
          dplyr::select(sea_offshore_distance_to_mouth, operation_id, dist_mouth_m),
          dplyr::transmute(dplyr::filter(inland_distance_to_mouth, survey != "river"), operation_id, dist_mouth_m = -dist_mouth_m)
        ),
        by = "operation_id"
      ) |>
      link_sea_operation(pomet_operation, coastal_water, sea_environment, 15)
  ),
  tar_target(
    estuary_guild_hand, # Guilds set by hand: species missing from the POMET reference, and corrections.
    tibble::tribble(
      ~species, ~guild, ~note,
      "Cottus gobio", "FW", "Freshwater species (estuarine resident in the POMET reference)",
      "Chelon labrosus", "MM", "Spawns at sea; juveniles grow in estuaries and lagoons, adults feed there seasonally and return to sea to spawn, so not resident; marine migrant as in most European classifications, though some class mullets as estuarine opportunists (diadromous in the POMET reference)",
      "Chelon auratus", "MM", "Spawns at sea; juveniles grow in estuaries and lagoons, adults feed there seasonally and return to sea to spawn, so not resident; marine migrant as in most European classifications, though some class mullets as estuarine opportunists (diadromous in the POMET reference)",
      "Neogobius melanostomus", "FW", "Euryhaline invasive goby, breeding in brackish and fresh water in its native range (could be an estuarine resident); classed freshwater because in our data it is caught almost only in the fresh tidal Seine (34 of 38 operations at salinity 0.5 or less, 2016-2018), where it arrived through rivers",
      "Squalius cephalus", "FW", "Freshwater species",
      "Barbatula barbatula", "FW", "Freshwater species",
      "Scophthalmus maximus", "MM", "Juveniles in shallow coastal nurseries",
      "Engraulis maeoticus", "MM", "Black Sea anchovy, probably a misidentified European anchovy (Engraulis encrasicolus, marine migrant)",
      "Mullus barbatus", "MS", "Marine species, occasional in estuaries",
      "Chelidonichthys lastoviza", "MS", "Marine species, occasional in estuaries",
      "Chelidonichthys cuculus", "MS", "Marine species, occasional in estuaries",
      "Raja microocellata", "MS", "Marine species, occasional in estuaries",
      "Lophius piscatorius", "MS", "Marine species, occasional in estuaries",
      "Trachinus draco", "MS", "Marine species, occasional in estuaries"
    )
  ),
  tar_target(
    species_guild, # Estuarine use guild of each species caught (POMET reference with validated names, completed and corrected by hand).
    assign_estuary_guild(
      dplyr::filter(sea_data_no_rare$catch, trait %in% sea_operation$operation_id)$species_valid,
      estuary_guild, estuary_guild_hand
    )
  ),
  tar_target(
    sea_catch, # Species caught at each sea operation, with their estuarine guild.
    sea_data_no_rare$catch |>
      tibble::as_tibble() |>
      dplyr::filter(trait %in% sea_operation$operation_id) |>
      dplyr::summarise(abundance = sum(abundance), .by = c(trait, species_valid)) |>
      dplyr::transmute(operation_id = trait, species = species_valid, abundance) |>
      dplyr::left_join(dplyr::select(species_guild, species, guild), by = "species")
  ),
  tar_target(
    sea_mouth, # River mouths of the four districts (AMOBIO), with the name of their river, pooled within 3 km (mouth_id).
    name_river_mouth(
      dplyr::filter(river_mouth_district, matched, district %in% c("Loire", "Ardour-Garonne", "Seine", "Escaut-Somme")),
      river_network, 50
    ) |>
      pool_river_mouth(3000)
  ),
  tar_target(
    mouth_size, # Catchment area of every (pooled) mouth, to select the main rivers: network routed to all mouths, grid outlets followed upstream up to 25 cells to the network, within 50 km of a mouth.
    find_mouth_outlet(river_network, label_network_mouth(river_network, sea_mouth), sea_mouth, flow_direction_file, 25, 50000) |>
      dplyr::summarise(area = sum(area), .by = mouth_id)
  ),
  tar_target(
    main_mouth, # Mouths of the main rivers (catchment of 1,000 km² or more).
    dplyr::filter(sea_mouth, mouth_id %in% mouth_size$mouth_id[mouth_size$area >= 1000])
  ),
  tar_target(
    network_mouth, # Main river mouth reached by each node of the river network (mouth_id, indexed by node id): small rivers joining an estuary join its main river.
    label_network_mouth(river_network, main_mouth)
  ),
  tar_target(
    mouth_outlet, # Grid outlets of the main rivers (same rule as mouth_size, network routed to the main mouths only).
    find_mouth_outlet(river_network, network_mouth, main_mouth, flow_direction_file, 25, 50000)
  ),
  tar_target(
    mouth_basin, # Catchment of each main river, as polygons (grid coarsened 3 times), for maps.
    build_mouth_basin(mouth_outlet, flow_direction_file, 3)
  ),
  tar_target(
    sea_mouth_area, # Main rivers: name and position (mouth nearest the largest outlet), number of mouths and outlets, catchment area (km²).
    summarise_mouth_pool(main_mouth, mouth_outlet) |>
      dplyr::rows_update(
        tibble::tribble( # Names set by hand: pools named after a small member, canalised rivers.
          ~mouth_id, ~river,
          191442, "la Gironde (Garonne, Dordogne)",
          826832, "la Somme",
          769167, "la Rance",
          790194, "le Trieux",
          791271, "le Jaudy",
          770997, "la Sée",
          741250, "l'Aulne",
          666825, "la Laïta"
        ),
        by = "mouth_id", unmatched = "ignore"
      )
  ),
  tar_target(
    sea_operation_mouth, # Main river of each sea operation: along the river network inland, the nearest main mouth along the sea within 15 km at sea; others dropped. By hand: the Gironde estuary (FRFT09 lower, FRFT35 upper), as the AMOBIO Gironde ends at its head (Bec d'Ambès), 20-70 km from its operations.
    link_operation_mouth(sea_operation, inland_operations_snapped, network_mouth, main_mouth, sea_graph, 15000) |>
      dplyr::bind_rows(
        sea_operation |>
          dplyr::filter(zone_id %in% c("FRFT09", "FRFT35")) |>
          dplyr::transmute(operation_id, mouth_id = 191442, link = "Gironde estuary (by hand)")
      ) |>
      dplyr::distinct(operation_id, .keep_all = TRUE)
  ),
  tar_target(
    mouth_landcover, # Share of each CORINE level-2 class in the catchment of the outlets.
    compute_station_landcover(corine_catchment_file, dplyr::distinct(mouth_outlet, code_station, cell)),
    pattern = map(corine_catchment_file)
  ),
  tar_target(
    mouth_weather, # Yearly weather of the catchment of the outlets.
    compute_station_weather(worldclim_catchment_file, dplyr::distinct(mouth_outlet, code_station, cell))
  ),
  tar_target(
    mouth_catchment, # Catchment land cover, climate and yearly weather anomalies of each main river.
    summarise_mouth_catchment(mouth_outlet, mouth_landcover, mouth_weather)
  ),
  tar_target(
    guild_visit, # Guild counts at the operations linked to a main river, with scaled drivers.
    prepare_guild_visit(sea_operation, sea_operation_mouth, mouth_catchment, sea_catch)
  ),
  tar_target(
    guild_mesh, # Fine mesh of the spatial field of the guild models.
    build_guild_mesh(dplyr::distinct(guild_visit, operation_id, longitude, latitude))
  ),
  tar_target(
    guild_term, # Fixed effects of the guild models (besides gear and month).
    "salinity + salinity2 + intensive + log_area + precipitation_anomaly + temperature_anomaly"
  ),
  tar_target(
    guild_code, # Guilds modelled.
    c("MS", "MM", "ER", "CA", "FW")
  ),
  tar_target(
    guild_model, # Model of the share of each guild among the species of an operation (MM timed first: 35 s).
    fit_guild_model(guild_visit, guild_mesh, guild_code, guild_term),
    pattern = map(guild_code),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  )
)
