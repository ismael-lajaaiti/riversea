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
    pomet_size_file, # POMET catches by size, with salinity, temperature and oxygen of each operation.
    here("data", "sea", "raw", "poisson_tailles_pomet.csv"),
    format = "file"
  ),
  tar_target(
    pomet_water, # Water temperature and oxygen of each POMET operation (suspect oxygen campaigns set to NA).
    read_pomet_water(pomet_size_file)
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
    sea_operation, # Estuary and coastal operations of the four districts (at sea and inland; two misplaced ones dropped), with signed distance to the mouth (negative upstream), gear, water body and salinity (measured, else environment grid).
    operation_location |>
      dplyr::filter(survey %in% c("pomet", "nurse"), district %in% c("Loire", "Ardour-Garonne", "Seine", "Escaut-Somme")) |>
      dplyr::left_join(
        dplyr::bind_rows(
          dplyr::select(sea_offshore_distance_to_mouth, operation_id, dist_mouth_m),
          dplyr::transmute(dplyr::filter(inland_distance_to_mouth, survey != "river"), operation_id, dist_mouth_m = -dist_mouth_m)
        ),
        by = "operation_id"
      ) |>
      link_sea_operation(pomet_operation, coastal_water, sea_environment, 15) |>
      dplyr::filter(!operation_id %in% c(61831282, 61831526)) # 2012 "Baie de Morlaix" (salinity 34) placed in central Brittany: wrong coordinates.
  ),
  tar_target(
    sea_district_extent, # Box (longitude, latitude) of the sea operations of each district, plus 0.1 degree.
    sea_operation |>
      dplyr::summarise(
        xmin = min(longitude) - 0.1, ymin = min(latitude) - 0.1,
        xmax = max(longitude) + 0.1, ymax = max(latitude) + 0.1,
        .by = district
      )
  ),
  tar_target(
    emodnet_substrate_file, # EMODnet Geology seabed substrate (1:250,000, Folk classes) of each district.
    download_emodnet_substrate(
      here("data", "geo", "emodnet"),
      unlist(sea_district_extent[c("xmin", "ymin", "xmax", "ymax")]),
      gsub("\\W", "", sea_district_extent$district)
    ),
    pattern = map(sea_district_extent),
    format = "file"
  ),
  tar_target(
    emodnet_bathymetry_file, # EMODnet bathymetry (2022, 1/16 arc-minute) of each district.
    download_emodnet_bathymetry(
      here("data", "geo", "emodnet"),
      unlist(sea_district_extent[c("xmin", "ymin", "xmax", "ymax")]),
      gsub("\\W", "", sea_district_extent$district)
    ),
    pattern = map(sea_district_extent),
    format = "file"
  ),
  tar_target(
    sea_operation_seabed, # Seabed substrate (4 classes + unmapped) and EMODnet depth of each sea operation.
    link_sea_seabed(sea_operation, emodnet_substrate_file, emodnet_bathymetry_file)
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
    all_mouth_outlet, # Grid outlets of every (pooled) mouth: network routed to all mouths, grid outlets followed upstream up to 25 cells to the network, within 50 km of a mouth.
    find_mouth_outlet(river_network, label_network_mouth(river_network, sea_mouth), sea_mouth, flow_direction_file, 25, 50000)
  ),
  tar_target(
    mouth_size, # Catchment area of every (pooled) mouth, to select the main rivers.
    dplyr::summarise(all_mouth_outlet, area = sum(area), .by = mouth_id)
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
  ),
  tar_target(
    sea_presence_visit, # Operations of the continuum presence models (salinity, distance, depth, substrate).
    prepare_sea_presence_visit(sea_operation, sea_operation_seabed)
  ),
  tar_target(
    pomet_presence_visit, # POMET operations with measured salinity, temperature and depth (oxygen dropped: no effect on sole or flounder).
    sea_presence_visit |>
      dplyr::filter(survey == "pomet") |>
      dplyr::inner_join(pomet_water, by = "operation_id") |>
      dplyr::filter(!is.na(temperature)) |>
      dplyr::select(-depth_nurse, -oxygen)
  ),
  tar_target(
    pomet_estuary, # Estuary of each POMET zone (and the lower Seine water body), set by hand: the Gironde and Seine zones pooled; Blavet-Scorff and Aven-Belon pooled as they share one mouth; type: bay (coastal water body) or estuary.
    tibble::tribble(
      ~zone_id, ~estuary, ~type,
      "FRFT09", "Gironde", "estuary", "FRFT35", "Gironde", "estuary", "FRFT31", "Gironde", "estuary",
      "FRFT32", "Gironde", "estuary", "FRFT33", "Gironde", "estuary",
      "FRHT01", "Seine", "estuary", "FRHT02", "Seine", "estuary", "FRHT03", "Seine", "estuary", # FRHT03: lower Seine, no POMET zone
      "FRHT07", "Seine", "estuary", "FRHC16", "Seine", "estuary",
      "FRHT08", "Dives", "estuary", "FRHC14", "Baie de Caen (Orne)", "bay", "FRHC10", "Baie des Veys", "bay",
      "FRGC01", "Mont-Saint-Michel", "bay", "FRAC05", "Somme", "bay", "FRFC02", "Pertuis Charentais", "bay",
      "FRGC11", "Baie de Morlaix", "bay", "FRGT02", "Rance", "estuary", "FRGT03", "Trieux", "estuary",
      "FRGT04", "Jaudy", "estuary", "FRGT05", "Léguer", "estuary", "FRGT08", "Aber Wrac'h", "estuary",
      "FRGT09", "Aber Benoit", "estuary", "FRGT10", "Élorn", "estuary", "FRGT12", "Aulne", "estuary",
      "FRGT13", "Goyen", "estuary", "FRGT14", "Pont-l'Abbé", "estuary", "FRGT15", "Odet", "estuary",
      "FRGT16", "Aven-Belon", "estuary", "FRGT17", "Aven-Belon", "estuary", "FRGT18", "Laïta", "estuary",
      "FRGT19", "Blavet-Scorff", "estuary", "FRGT20", "Blavet-Scorff", "estuary", "FRGT21", "Etel", "estuary",
      "FRGT23", "Auray", "estuary", "FRGT27", "Vilaine", "estuary", "FRGT28", "Loire", "estuary",
      "FRGT30", "Lay", "estuary", "FRGT31", "Sèvre Niortaise", "estuary"
    )
  ),
  tar_target(
    pomet_estuary_mouth, # River mouths (sea_mouth mouth_id) of each estuary, set by hand.
    tibble::tribble(
      ~estuary, ~mouth_id,
      "Gironde", 191442, "Seine", 815337, "Dives", 809653, "Baie de Caen (Orne)", 807792,
      "Baie des Veys", 813905, "Baie des Veys", 813711, # Douve (with the Taute), Vire
      "Mont-Saint-Michel", 763892, "Mont-Saint-Michel", 770997, # Couesnon, Sée (with the Sélune)
      "Somme", 826832, "Pertuis Charentais", 252754, "Pertuis Charentais", 241976, # Charente, Seudre
      "Baie de Morlaix", 778912, "Baie de Morlaix", 782210, # Dourduff (Morlaix river), Penzé
      "Rance", 769167, "Trieux", 790194, "Jaudy", 791271, "Léguer", 785819, "Aber Wrac'h", 781312,
      "Aber Benoit", 778408, "Élorn", 753292, "Aulne", 741250, "Goyen", 708292, "Pont-l'Abbé", 682347,
      "Odet", 685607, "Aven-Belon", 675067, "Laïta", 666825, "Blavet-Scorff", 663195, "Etel", 655616,
      "Auray", 648078, "Vilaine", 640183, "Loire", 623205, "Lay", 550316, "Sèvre Niortaise", 549691
    )
  ),
  tar_target(
    transitional_water_file, # WFD transitional water bodies (Sandre 2022), with the Seine ones missing from the 2019 polygons.
    download_transitional_water(here("data", "geo", "transitional")),
    format = "file"
  ),
  tar_target(
    pomet_estuary_descriptor, # Catchment area and area (transitional water bodies, else coastal water body) of each estuary.
    describe_pomet_estuary(pomet_estuary, pomet_estuary_mouth, mouth_size, transitional_water_file, coastal_water)
  ),
  tar_target(
    pomet_presence_mesh, # SPDE mesh over the POMET operations.
    build_guild_mesh(pomet_presence_visit)
  ),
  tar_target(
    pomet_presence_variant, # POMET-only model variants of the test species: one change to the reference each (substrate and oxygen tested, not in the reference: no gain).
    {
      substrate <- "substratesand + substratecoarse + substraterock + substrateunmapped"
      reference <- c(
        salinity = "log_salinity + log_salinity2", distance = "distance", depth = "depth_pomet",
        temperature = "temperature", seabed = ""
      )
      change <- \(...) {
        x <- utils::modifyList(as.list(reference), list(...))
        paste(unlist(x)[nzchar(unlist(x))], collapse = " + ")
      }
      term <- c(
        reference = change(),
        `linear log salinity` = change(salinity = "log_salinity"),
        `quadratic temperature` = change(temperature = "temperature + temperature2"),
        `no temperature` = change(temperature = ""),
        `no distance` = change(distance = ""),
        `no depth` = change(depth = ""),
        `with substrate` = change(seabed = substrate),
        `no field` = change()
      )
      tidyr::expand_grid(species = c("Solea solea", "Platichthys flesus"), variant = names(term)) |>
        dplyr::mutate(term = term[variant], field = variant != "no field", prior_range = list(c(20, 0.5)))
    }
  ),
  tar_target(
    pomet_presence_model, # POMET-only presence model of each test species and variant.
    c(
      list(variant = pomet_presence_variant$variant),
      fit_sea_presence_model(
        pomet_presence_visit, sea_catch, if (pomet_presence_variant$field) pomet_presence_mesh else NULL,
        pomet_presence_variant$species, pomet_presence_variant$term, pomet_presence_variant$prior_range[[1]],
        character(0)
      )
    ),
    pattern = map(pomet_presence_variant),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    pomet_oxygen_model, # Oxygen test: POMET reference model with and without oxygen, on the operations with oxygen.
    {
      visit <- sea_presence_visit |>
        dplyr::filter(survey == "pomet") |>
        dplyr::inner_join(pomet_water, by = "operation_id") |>
        dplyr::filter(!is.na(temperature), !is.na(oxygen))
      term <- unname(pomet_presence_variant$term[pomet_presence_variant$variant == "reference"][1])
      term <- c(`without oxygen` = term, `with oxygen` = paste(term, "+ oxygen"), `quadratic oxygen` = paste(term, "+ oxygen + oxygen2"))
      purrr::map(names(term), \(v) c(
        list(variant = v),
        fit_sea_presence_model(visit, sea_catch, pomet_presence_mesh, pomet_species, term[[v]], c(20, 0.5), character(0))
      ))
    },
    pattern = map(pomet_species),
    iteration = "list"
  ),
  tar_target(pomet_species, c("Solea solea", "Platichthys flesus")),
  tar_target(
    pomet_estuary_visit, # POMET operations with the descriptors of their estuary: log catchment area, log estuary area, bay or estuary.
    pomet_presence_visit |>
      dplyr::inner_join(dplyr::select(pomet_estuary, zone_id, estuary, type), by = "zone_id") |>
      dplyr::left_join(pomet_estuary_descriptor, by = "estuary") |>
      dplyr::mutate(log_catchment = log(catchment), log_estuary_area = log(estuary_area), bay = as.numeric(type == "bay"))
  ),
  tar_target(
    pomet_estuary_variant, # Estuary descriptors added to the POMET reference model, with and without spatial field.
    {
      reference <- unname(pomet_presence_variant$term[pomet_presence_variant$variant == "reference"][1])
      added <- c(
        `no descriptor` = "", catchment = "log_catchment", area = "log_estuary_area", type = "bay",
        `all three` = "log_catchment + log_estuary_area + bay"
      )
      tidyr::expand_grid(species = pomet_species, descriptor = names(added), field = c(TRUE, FALSE)) |>
        dplyr::mutate(
          variant = paste0(descriptor, dplyr::if_else(field, "", ", no field")),
          term = dplyr::if_else(nzchar(added[descriptor]), paste(reference, "+", added[descriptor]), reference)
        )
    }
  ),
  tar_target(
    pomet_estuary_model, # POMET presence models with estuary descriptors.
    c(
      list(variant = pomet_estuary_variant$variant, descriptor = pomet_estuary_variant$descriptor),
      fit_sea_presence_model(
        pomet_estuary_visit, sea_catch, if (pomet_estuary_variant$field) pomet_presence_mesh else NULL,
        pomet_estuary_variant$species, pomet_estuary_variant$term, c(20, 0.5), character(0)
      )
    ),
    pattern = map(pomet_estuary_variant),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    hydrometry_station, # Gauging stations (HYDRO, Hub'Eau) in the districts studied.
    get_hydrometry_station(dplyr::filter(basin, district %in% params$district_kept))
  ),
  tar_target(
    hydrometry_station_flow, # Cell of the flow-direction grid of each gauging station.
    snap_station_flow(hydrometry_station, rht, flow_direction_file, 1000)
  ),
  tar_target(
    hydrometry_station_mouth, # River mouth reached by each gauging station along the flow directions, and drained area.
    link_station_mouth(hydrometry_station_flow, all_mouth_outlet, flow_direction_file)
  ),
  tar_target(
    hydrometry_batch, # Gauging stations draining to the POMET estuaries, queried together.
    {
      code <- sort(hydrometry_station_mouth$code_station[hydrometry_station_mouth$mouth_id %in% pomet_estuary_mouth$mouth_id])
      split(code, ceiling(seq_along(code) / 50))
    },
    iteration = "list"
  ),
  tar_target(
    hydrometry_discharge, # Monthly mean discharge of these stations, 2005-2019 (m3/s).
    download_monthly_discharge(hydrometry_batch, c("2005-01-01", "2019-12-31")),
    pattern = map(hydrometry_batch)
  ),
  tar_target(
    estuary_discharge, # Yearly discharge entering each estuary: runoff of the most downstream gauging station of each mouth measured that year (10+ complete years), times the catchment.
    summarise_estuary_discharge(hydrometry_discharge, hydrometry_station_mouth, pomet_estuary_mouth, mouth_size, c(2005, 2019), 10)
  ),
  tar_target(
    estuary_nutrient_mix, # Yearly river inputs of each estuary, mixing the non-nested stations of its branches in proportion to their discharge.
    mix_estuary_nutrient(
      estuary_nutrient$station, naiades_station_year, station_flow,
      dplyr::inner_join(hydrometry_station_mouth, dplyr::select(hydrometry_station_flow, code_station, cell), by = "code_station"),
      estuary_discharge, flow_direction_file, c(2005, 2019), 0.8
    )
  ),
  tar_target(
    pomet_nutrient_visit, # POMET operations with the yearly river inputs of their estuary (mixed branches: log10 mg/l, oxygen saturation %; log10 discharge, m3/s), each split into the estuary mean (`_between`) and the deviation of the year (`_within`); operations without all of them dropped.
    {
      river <- c("nitrate", "ammonium", "orthophosphate", "river_oxygen", "log_discharge")
      estuary_year <- estuary_nutrient_mix$estuary_year |>
        dplyr::select(-cover) |>
        dplyr::inner_join(
          dplyr::transmute(estuary_discharge$estuary_year, estuary, year, log_discharge = log10(discharge)),
          by = c("estuary", "year")
        ) |>
        dplyr::mutate(
          dplyr::across(dplyr::all_of(river), list(between = mean, within = \(x) x - mean(x)), .names = "{.col}_{.fn}"),
          .by = estuary
        )
      pomet_estuary_visit |>
        dplyr::select(-catchment, -estuary_area, -log_catchment, -log_estuary_area, -bay) |>
        dplyr::inner_join(estuary_year, by = c("estuary", "year"))
    }
  ),
  tar_target(
    pomet_nutrient_variant, # Each river input added to the POMET reference model, alone, with and without spatial field.
    {
      reference <- unname(pomet_presence_variant$term[pomet_presence_variant$variant == "reference"][1])
      tidyr::expand_grid(species = pomet_species, nutrient = c("none", "nitrate", "ammonium", "orthophosphate"), field = c(TRUE, FALSE)) |>
        dplyr::mutate(
          variant = paste0(nutrient, dplyr::if_else(field, "", ", no field")),
          term = dplyr::if_else(nutrient == "none", reference, paste(reference, "+", nutrient))
        )
    }
  ),
  tar_target(
    pomet_nutrient_model, # POMET presence models with river nutrients.
    c(
      list(variant = pomet_nutrient_variant$variant, nutrient = pomet_nutrient_variant$nutrient),
      fit_sea_presence_model(
        pomet_nutrient_visit, sea_catch, if (pomet_nutrient_variant$field) pomet_presence_mesh else NULL,
        pomet_nutrient_variant$species, pomet_nutrient_variant$term, c(20, 0.5), character(0)
      )
    ),
    pattern = map(pomet_nutrient_variant),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    pomet_split_variant, # River inputs split into estuary mean and yearly deviation, alone or with discharge, with an iid estuary effect, with and without spatial field.
    {
      reference <- unname(pomet_presence_variant$term[pomet_presence_variant$variant == "reference"][1])
      split <- \(x) paste0(x, "_between + ", x, "_within")
      input <- c("nitrate", "ammonium", "orthophosphate", "river_oxygen")
      added <- c(
        none = "", log_discharge = split("log_discharge"),
        stats::setNames(purrr::map_chr(input, split), input),
        stats::setNames(paste(purrr::map_chr(input, split), "+", split("log_discharge")), paste0(input, " + log_discharge"))
      )
      tidyr::expand_grid(species = pomet_species, input = names(added), field = c(TRUE, FALSE)) |>
        dplyr::mutate(
          variant = paste0(input, dplyr::if_else(field, "", ", no field")),
          term = paste(reference, "+ estuary", dplyr::if_else(nzchar(added[input]), paste("+", added[input]), ""))
        )
    }
  ),
  tar_target(
    pomet_split_model, # POMET presence models with the between- and within-estuary parts of each river input.
    c(
      list(variant = pomet_split_variant$variant, input = pomet_split_variant$input),
      fit_sea_presence_model(
        pomet_nutrient_visit, sea_catch, if (pomet_split_variant$field) pomet_presence_mesh else NULL,
        pomet_split_variant$species, pomet_split_variant$term, c(20, 0.5), character(0)
      )
    ),
    pattern = map(pomet_split_variant),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    station_mouth, # River mouth reached by each NAIADES station along the flow directions.
    link_station_mouth(station_flow, all_mouth_outlet, flow_direction_file)
  ),
  tar_target(
    estuary_nutrient, # Yearly river nutrients of each estuary: most downstream NAIADES station of each mouth, outside the transitional and coastal waters, with 10+ years in 2005-2019, catchment-weighted.
    summarise_estuary_nutrient(
      naiades_station_year, station_mouth, pomet_estuary_mouth, mouth_size, naiades_station,
      dplyr::bind_rows(
        sf::st_transform(sf::st_read(transitional_water_file, quiet = TRUE)["CdEuMasseD"], 4326),
        sf::st_transform(coastal_water["zone_id"], 4326)
      ),
      c(2005, 2019), 10
    )
  ),
  tar_target(
    pomet_presence_block, # Estuary systems of the POMET operations for cross-validation (steps under 10 km).
    block_sea_operations(pomet_presence_visit, 10000)
  ),
  tar_target(
    pomet_presence_fold, # Cross-validation folds: each block left out, reference model with and without field.
    dplyr::filter(pomet_presence_variant, variant %in% c("reference", "no field")) |>
      tidyr::expand_grid(block = sort(unique(pomet_presence_block$block)))
  ),
  tar_target(
    pomet_presence_cv, # POMET reference model fitted without one block; predictions at the left-out operations and fixed effects.
    {
      test <- pomet_presence_block$operation_id[pomet_presence_block$block == pomet_presence_fold$block]
      fit <- fit_sea_presence_model(
        pomet_presence_visit, sea_catch, if (pomet_presence_fold$field) pomet_presence_mesh else NULL,
        pomet_presence_fold$species, pomet_presence_fold$term, pomet_presence_fold$prior_range[[1]], test
      )
      list(
        species = pomet_presence_fold$species, variant = pomet_presence_fold$variant, block = pomet_presence_fold$block,
        summary.fixed = fit$summary.fixed, prediction = dplyr::filter(fit$prediction, operation_id %in% test)
      )
    },
    pattern = map(pomet_presence_fold),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  )
)
