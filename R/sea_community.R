#' Gear, water body and depth of the POMET operations.
#'
#' @param file path to `poisson_pomet.csv`.
#' @return tibble, one row per operation: `operation_id`, `gear` (beam
#'   trawl width), `zone_id` (Water Framework Directive water body code,
#'   e.g. "FRGT20"), `zone_name` and `depth` (m).
#' @export
read_pomet_operation <- function(file) {
  utils::read.csv(file, sep = ";", fileEncoding = "latin1") |>
    tibble::as_tibble() |>
    dplyr::distinct(
      operation_id = as.character(ID_interne_prelevement),
      gear = sub("chalut a perche_", "beam trawl ", Engin_peche),
      zone = ID_Zone, depth = as.numeric(sub(",", ".", Prof))
    ) |>
    dplyr::mutate(
      zone_id = trimws(sub(" - .*", "", zone)),
      zone_name = trimws(sub("^[^-]* - ", "", zone)),
      .keep = "unused"
    ) |>
    dplyr::distinct(operation_id, .keep_all = TRUE)
}

#' Estuarine use guild of the fish species (POMET reference).
#'
#' @param file path to `poisson_pomet.csv`.
#' @return tibble with `species` (name validated against FishBase, as the
#'   catch data, e.g. *Liza ramada* becomes *Chelon ramada*) and `guild`:
#'   "MS" marine straggler, "MM" marine migrant (juveniles using estuaries),
#'   "ER" estuarine resident, "FW" freshwater, "CA" diadromous
#'   (catadromous and anadromous).
#' @export
read_estuary_guild <- function(file) {
  guild <- utils::read.csv(file, sep = ";", fileEncoding = "latin1") |>
    dplyr::filter(!is.na(Ecologique), Ecologique != "") |>
    dplyr::distinct(species = NomScient, guild = Ecologique)
  guild |>
    dplyr::mutate(species = dplyr::coalesce(suppressMessages(rfishbase::validate_names(species)), species)) |>
    dplyr::distinct(species, .keep_all = TRUE)
}

#' Estuarine use guild of the species caught, with corrections by hand.
#'
#' @param catch_species character, species of the catches.
#' @param guild tibble, [read_estuary_guild()]'s output.
#' @param hand tibble with `species`, `guild` and `note`: guilds set by hand,
#'   for species missing from the POMET reference or whose guild is
#'   corrected.
#' @return tibble with `species`, `guild`, `source` ("POMET reference",
#'   "corrected by hand" or "added by hand") and `note`.
#' @export
assign_estuary_guild <- function(catch_species, guild, hand) {
  tibble::tibble(species = sort(unique(catch_species))) |>
    dplyr::left_join(guild, by = "species") |>
    dplyr::left_join(dplyr::rename(hand, guild_hand = guild), by = "species") |>
    dplyr::mutate(
      source = dplyr::case_when(
        is.na(guild_hand) ~ "POMET reference",
        is.na(guild) ~ "added by hand",
        .default = "corrected by hand"
      ),
      guild = dplyr::coalesce(guild_hand, guild)
    ) |>
    dplyr::select(species, guild, source, note)
}

#' Yearly gridded environment of the Bay of Biscay coast.
#'
#' @param dir path to the folder of `Environment_QM_<year>.csv` files.
#' @return tibble with `year`, `longitude`, `latitude`, `salinity`
#'   (`sali_mean`), `temperature` (`temp_mean`), `chlorophyll`
#'   (`chl_mean`) and `discharge_anomaly` (`QmMAnomalie_jan_avr`).
#' @export
read_sea_environment <- function(dir) {
  file <- list.files(dir, pattern = "^Environment_QM_\\d{4}\\.csv$", full.names = TRUE)
  purrr::map_dfr(file, \(f) {
    utils::read.csv(f) |>
      tibble::as_tibble() |>
      dplyr::transmute(
        year = as.integer(gsub("\\D", "", basename(f))),
        longitude = x, latitude = y, salinity = sali_mean, temperature = temp_mean,
        chlorophyll = chl_mean, discharge_anomaly = QmMAnomalie_jan_avr
      )
  })
}

#' Sea operations with gear, water body and salinity.
#'
#' Water body: the POMET zone for POMET operations; for the nursery survey,
#' the water body containing the operation, or the nearest within
#' `dist_max` km. Salinity: measured when available, otherwise the yearly
#' mean of the nearest cell of the environment grid of the same year
#' (within 5 km).
#'
#' @param operation tibble with `operation_id`, `survey`, `year`, `month`,
#'   `longitude`, `latitude` and `salinity`.
#' @param pomet tibble, [read_pomet_operation()]'s output.
#' @param coastal_water sf, water body polygons with `zone_id`.
#' @param environment tibble, [read_sea_environment()]'s output.
#' @param dist_max numeric, km.
#' @return `operation` with `gear`, `zone_id`, `zone_type` ("T"
#'   transitional, "C" coastal), `depth`, `salinity` (completed) and
#'   `salinity_source` ("measured", "grid" or NA).
#' @export
link_sea_operation <- function(operation, pomet, coastal_water, environment, dist_max) {
  point <- operation |>
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) |>
    sf::st_transform(2154)
  water <- sf::st_transform(coastal_water["zone_id"], 2154)
  nearest <- sf::st_nearest_feature(point, water)
  distance <- as.numeric(sf::st_distance(point, water[nearest, ], by_element = TRUE)) / 1000
  operation <- operation |>
    dplyr::mutate(zone_polygon = dplyr::if_else(distance <= dist_max, water$zone_id[nearest], NA)) |>
    dplyr::left_join(pomet, by = "operation_id") |>
    dplyr::mutate(
      gear = dplyr::if_else(survey == "nurse", "beam trawl nursery", gear),
      zone_id = dplyr::if_else(survey == "pomet", dplyr::coalesce(zone_id, zone_polygon), zone_polygon),
      zone_type = substr(zone_id, 4, 4)
    ) |>
    dplyr::select(-zone_polygon)
  grid <- sf::st_as_sf(environment, coords = c("longitude", "latitude"), crs = 4326) |>
    sf::st_transform(2154)
  grid_salinity <- purrr::map_dbl(seq_len(nrow(operation)), \(i) {
    cell <- grid[environment$year == operation$year[i], ]
    if (nrow(cell) == 0) return(NA_real_)
    k <- sf::st_nearest_feature(point[i, ], cell)
    if (as.numeric(sf::st_distance(point[i, ], cell[k, ])) > 5000) return(NA_real_)
    cell$salinity[k]
  })
  operation |>
    dplyr::mutate(
      salinity_source = dplyr::case_when(
        !is.na(salinity) ~ "measured", !is.na(grid_salinity) ~ "grid", .default = NA
      ),
      salinity = dplyr::coalesce(salinity, grid_salinity)
    )
}

#' Name of the river of each mouth.
#'
#' The river of a mouth is the name of the nearest named reach (the reach
#' at the mouth is often unnamed, e.g. the Loire), searched up to `order`
#' nodes away.
#'
#' @param mouth sf, river mouths with `node_id`, e.g. the
#'   `river_mouth_district` target.
#' @param network list(nodes, edges), e.g. the `river_network` target.
#' @param order integer, maximum number of nodes searched from the mouth.
#' @return `mouth` as a tibble with `longitude`, `latitude` and `river` (NA
#'   if no named reach).
#' @export
name_river_mouth <- function(mouth, network, order) {
  edge <- sf::st_drop_geometry(network$edges)
  g <- igraph::graph_from_edgelist(as.matrix(edge[c("from", "to")]), directed = FALSE)
  river <- vapply(mouth$node_id, \(id) {
    for (k in seq_len(order)) {
      near <- igraph::ego(g, k, id)[[1]]
      reach <- stats::na.omit(edge$name[igraph::E(g)[.inc(near)]])
      if (length(reach) > 0) return(names(which.max(table(reach))))
    }
    NA_character_
  }, character(1))
  coordinate <- sf::st_coordinates(sf::st_transform(mouth, 4326))
  sf::st_drop_geometry(mouth) |>
    tibble::as_tibble() |>
    dplyr::mutate(longitude = coordinate[, 1], latitude = coordinate[, 2], river = river)
}

#' Pool river mouths close to each other.
#'
#' Mouths closer than `dist` m (single linkage) share an estuary or a bay
#' (e.g. the Orne and the Dan at Caen) and are pooled: their catchments
#' are summed and their operations share them.
#'
#' @param mouth tibble with `node_id`, `longitude`, `latitude`, e.g.
#'   [name_river_mouth()]'s output.
#' @param dist numeric, m.
#' @return `mouth` with `mouth_id`: node id of the first mouth of the pool.
#' @export
pool_river_mouth <- function(mouth, dist) {
  point <- sf::st_transform(sf::st_as_sf(mouth, coords = c("longitude", "latitude"), crs = 4326), 2154)
  near <- sf::st_is_within_distance(point, point, dist = dist)
  pool <- seq_len(nrow(mouth))
  repeat {
    updated <- vapply(near, \(i) min(pool[i]), integer(1))
    if (identical(updated, pool)) break
    pool <- updated
  }
  dplyr::mutate(mouth, mouth_id = node_id[pool])
}

#' Pooled river mouths with their catchment size.
#'
#' A pool of mouths is named and placed after the mouth nearest to its
#' largest outlet (e.g. the Garonne for the Gironde).
#'
#' @param mouth tibble, [pool_river_mouth()]'s output.
#' @param mouth_outlet tibble, [find_mouth_outlet()]'s output.
#' @return tibble, one row per pool with an outlet: `mouth_id`,
#'   `district`, `river`, `longitude`, `latitude`, `mouths`, `outlets` and
#'   `area` (km²).
#' @export
summarise_mouth_pool <- function(mouth, mouth_outlet) {
  largest <- dplyr::slice_max(mouth_outlet, area, by = mouth_id, with_ties = FALSE)
  mouth |>
    dplyr::inner_join(dplyr::select(largest, mouth_id, x = longitude, y = latitude), by = "mouth_id") |>
    dplyr::mutate(mouths = dplyr::n(), .by = mouth_id) |>
    dplyr::slice_min((longitude - x)^2 + (cos(y * pi / 180) * (latitude - y))^2, by = mouth_id, with_ties = FALSE) |>
    dplyr::select(mouth_id, district, river, longitude, latitude, mouths) |>
    dplyr::inner_join(
      dplyr::summarise(mouth_outlet, outlets = dplyr::n(), area = sum(area), .by = mouth_id),
      by = "mouth_id"
    )
}

#' River mouth reached by each node of the river network.
#'
#' Each node follows its cheapest route to the nearest mouth along the
#' network (`route_to_mouth()`, avoiding artificial channels).
#'
#' @param network list(nodes, edges), e.g. the `river_network` target.
#' @param mouth tibble with `node_id` and `mouth_id`, e.g.
#'   [pool_river_mouth()]'s output.
#' @return vector indexed by node id: `mouth_id` of the mouth reached (NA
#'   if none).
#' @export
label_network_mouth <- function(network, mouth) {
  route <- route_to_mouth(network, mouth$node_id)
  up <- route$parent
  up[mouth$node_id] <- mouth$node_id
  repeat { # Pointer jumping: each node points twice as far up at each pass, until the mouth.
    jump <- up[up]
    if (identical(jump, up)) break
    up <- jump
  }
  mouth$mouth_id[match(up[-route$source], mouth$node_id)]
}

#' Downstream cell of each cell of the flow-direction grid.
#'
#' @param direction SpatRaster of flow directions (ESRI D8 codes, 0 for
#'   outlets, NA for the sea).
#' @return integer vector, one per cell: the downstream cell, NA for the
#'   outlets (code 0, or flowing into the sea or out of the grid) and the
#'   sea.
#' @noRd
flow_downstream <- function(direction) {
  code <- terra::values(direction, mat = FALSE)
  ncol <- terra::ncol(direction)
  cell <- seq_along(code)
  down <- cell + c(1L, ncol + 1L, ncol, ncol - 1L, -1L, -ncol - 1L, -ncol, -ncol + 1L)[match(code, 2^(0:7))]
  down[!is.na(down) & (down < 1 | down > length(code))] <- NA
  down[!is.na(down) & abs((cell - 1) %% ncol - (down - 1) %% ncol) > 1] <- NA
  down[!is.na(down) & is.na(code[down])] <- NA
  down
}

#' Outlets of the flow-direction grid reaching each river mouth.
#'
#' The outlets are the cells of the flow-direction grid flowing into the
#' sea (or a sink); their catchments are disjoint and cover the land. From
#' each outlet within `dist_mouth` m of a mouth, the main river (upstream
#' cell draining the most) is followed for up to `step` cells; the outlet
#' takes the mouth of the first network nodes met (within one cell of the
#' river). Following the river matters where the network and the grid end
#' at different places: the AMOBIO Vilaine ends at the Arzal dam, upstream
#' of the grid outlet; on the grid, the Gironde is sea, so the Garonne and
#' the Dordogne end upstream of the mouth. The distance to the mouths
#' excludes the outlets of rivers joined to the network by canals but
#' reaching another sea (e.g. the Rhône).
#'
#' @param network list(nodes, edges), e.g. the `river_network` target.
#' @param network_mouth vector, [label_network_mouth()]'s output.
#' @param mouth tibble with `longitude` and `latitude`, the mouths.
#' @param direction_file path, e.g. `crop_flow_direction()`'s output.
#' @param step integer, maximum number of cells followed upstream.
#' @param dist_mouth numeric, m.
#' @return tibble with `mouth_id`, `code_station` (outlet id, named as the
#'   stations so that the catchment functions of D15 apply), `cell`,
#'   `longitude`, `latitude`, `area` (km² drained) and `step` (cells
#'   followed upstream to meet the network).
#' @export
find_mouth_outlet <- function(network, network_mouth, mouth, direction_file, step, dist_mouth) {
  direction <- terra::rast(direction_file)
  code <- terra::values(direction, mat = FALSE)
  down <- flow_downstream(direction)
  upstream <- terra::flowAccumulation(direction, weight = terra::cellSize(direction, unit = "km"))
  drained <- terra::values(upstream, mat = FALSE)
  # Mouth of the network nodes in each cell (most common), widened by one cell.
  node <- sf::st_transform(sf::st_as_sf(network$nodes), 4326)
  node_mouth <- tibble::tibble(
    cell = terra::cellFromXY(direction, sf::st_coordinates(node)),
    mouth_id = network_mouth[node$node_id]
  ) |>
    dplyr::filter(!is.na(cell), !is.na(mouth_id)) |>
    dplyr::count(cell, mouth_id) |>
    dplyr::slice_max(n, by = cell, with_ties = FALSE)
  near <- terra::rast(direction)
  near[node_mouth$cell] <- node_mouth$mouth_id
  near <- terra::values(terra::focal(near, 3, "modal", na.policy = "only", na.rm = TRUE), mat = FALSE)
  near[node_mouth$cell] <- node_mouth$mouth_id
  # Outlets near a mouth.
  cell <- which(!is.na(code) & is.na(down))
  point <- sf::st_as_sf(data.frame(cell = cell, terra::xyFromCell(direction, cell)), coords = c("x", "y"), crs = 4326)
  mouth <- sf::st_as_sf(mouth, coords = c("longitude", "latitude"), crs = 4326)
  cell <- cell[lengths(sf::st_is_within_distance(point, mouth, dist = dist_mouth)) > 0]
  # Main upstream cell of each cell: the upstream neighbour draining most.
  from <- which(!is.na(down))
  from <- from[order(drained[from])]
  main <- rep(NA_integer_, length(down))
  main[down[from]] <- from
  stem <- matrix(NA_integer_, length(cell), step + 1)
  stem[, 1] <- cell
  for (k in seq_len(step)) stem[, k + 1] <- main[stem[, k]]
  met <- matrix(near[stem], nrow(stem))
  first <- apply(!is.na(met), 1, \(x) match(TRUE, x))
  coordinate <- terra::xyFromCell(direction, cell)
  tibble::tibble(
    mouth_id = met[cbind(seq_along(cell), first)],
    code_station = paste0("cell_", cell),
    cell = cell,
    longitude = coordinate[, 1],
    latitude = coordinate[, 2],
    area = drained[cell],
    step = first - 1L
  ) |>
    dplyr::filter(!is.na(mouth_id))
}

#' Catchment of each river mouth, as polygons.
#'
#' Each cell of the flow-direction grid takes the outlet it drains to
#' (pointer jumping along the flow directions), then the mouth of that
#' outlet. The grid is coarsened by `fact` (most common mouth) before
#' conversion to polygons, for maps.
#'
#' @param mouth_outlet tibble, [find_mouth_outlet()]'s output.
#' @param direction_file path, e.g. `crop_flow_direction()`'s output.
#' @param fact integer, aggregation factor.
#' @return sf polygons (WGS84) with `mouth_id`.
#' @export
build_mouth_basin <- function(mouth_outlet, direction_file, fact) {
  direction <- terra::rast(direction_file)
  up <- flow_downstream(direction)
  outlet <- which(!is.na(terra::values(direction, mat = FALSE)) & is.na(up))
  up[outlet] <- outlet
  repeat {
    jump <- up[up]
    if (identical(jump, up)) break
    up <- jump
  }
  basin <- terra::rast(direction)
  terra::values(basin) <- mouth_outlet$mouth_id[match(up, mouth_outlet$cell)]
  basin <- terra::aggregate(basin, fact, "modal", na.rm = TRUE)
  names(basin) <- "mouth_id"
  sf::st_as_sf(terra::as.polygons(basin)) |>
    sf::st_make_valid()
}

#' Sea operations linked to their river mouth.
#'
#' Operations inland (upper estuaries) take the mouth reached along the
#' river network; operations at sea take the nearest mouth along the sea
#' (coast-aware path of D8), among `mouth`, if it is within `dist_max` m.
#' Operations without a mouth are dropped.
#'
#' @param operation tibble with `operation_id`, `longitude`, `latitude`,
#'   `inland`, e.g. the `sea_operation` target.
#' @param snapped tibble with `operation_id` and `node_id` (network node),
#'   e.g. the `inland_operations_snapped` target.
#' @param network_mouth integer vector, [label_network_mouth()]'s output.
#' @param mouth tibble with `node_id`, `mouth_id`, `longitude`,
#'   `latitude`: the mouths operations at sea can be linked to.
#' @param sea_graph igraph, e.g. the `sea_graph` target.
#' @param dist_max numeric, m.
#' @return tibble with `operation_id`, `mouth_id` and `link` ("river
#'   network" or "sea path").
#' @export
link_operation_mouth <- function(operation, snapped, network_mouth, mouth, sea_graph, dist_max) {
  inland <- operation |>
    dplyr::filter(inland) |>
    dplyr::inner_join(dplyr::select(snapped, operation_id, node_id), by = "operation_id") |>
    dplyr::transmute(operation_id, mouth_id = network_mouth[node_id], link = "river network")
  sea <- dplyr::filter(operation, !inland)
  sea_node <- snap_to_sea_graph(sf::st_as_sf(sea, coords = c("longitude", "latitude"), crs = 4326), sea_graph)$node_id
  mouth_node <- snap_to_sea_graph(sf::st_as_sf(mouth, coords = c("longitude", "latitude"), crs = 4326), sea_graph)$node_id
  distance <- igraph::distances(sea_graph, v = unique(mouth_node), to = unique(sea_node), weights = igraph::E(sea_graph)$weight)
  nearest <- mouth$mouth_id[match(unique(mouth_node), mouth_node)][apply(distance, 2, which.min)]
  nearest[apply(distance, 2, min) > dist_max] <- NA
  dplyr::bind_rows(
    inland,
    tibble::tibble(operation_id = sea$operation_id, mouth_id = nearest[match(sea_node, unique(sea_node))], link = "sea path")
  ) |>
    dplyr::filter(!is.na(mouth_id))
}

#' Catchment drivers of each river mouth.
#'
#' The catchment of a mouth is the union of the catchments of its outlets;
#' land cover shares and climate are their area-weighted means. Land cover:
#' arable land (CORINE 21), other farmland (22-24) and artificial surfaces
#' (11-14), mean of the vintages. Climate: mean precipitation and air
#' temperature over 2000-2021; weather of each year: log ratio of
#' precipitation to the mean (and of the previous year), and difference of
#' temperature to the mean.
#'
#' @param mouth_outlet tibble, [find_mouth_outlet()]'s output.
#' @param landcover tibble, `compute_station_landcover()`'s output for the
#'   outlets.
#' @param weather tibble, `compute_station_weather()`'s output for the
#'   outlets.
#' @return list with `mouth` (one row per mouth: `mouth_id`, `outlets`,
#'   `area`, `arable`, `other_farmland`, `artificial`, `precipitation`,
#'   `air_temperature`) and `mouth_year` (`mouth_id`, `year`,
#'   `precipitation_anomaly`, `precipitation_previous`,
#'   `temperature_anomaly`).
#' @export
summarise_mouth_catchment <- function(mouth_outlet, landcover, weather) {
  cover <- landcover |>
    dplyr::mutate(group = dplyr::case_when(
      class == 21 ~ "arable", class %in% 22:24 ~ "other_farmland", class %in% 11:14 ~ "artificial"
    )) |>
    dplyr::filter(!is.na(group)) |>
    dplyr::summarise(share = sum(share), .by = c(code_station, year, group)) |>
    dplyr::summarise(share = mean(share), .by = c(code_station, group)) |>
    tidyr::pivot_wider(names_from = group, values_from = share)
  climate <- weather |>
    dplyr::summarise(precipitation = mean(precipitation), air_temperature = mean(temperature), .by = code_station)
  mouth <- mouth_outlet |>
    dplyr::inner_join(cover, by = "code_station") |>
    dplyr::inner_join(climate, by = "code_station") |>
    dplyr::summarise(
      dplyr::across(c(arable, other_farmland, artificial, precipitation, air_temperature), \(x) stats::weighted.mean(x, area)),
      outlets = dplyr::n(),
      area = sum(area),
      .by = mouth_id
    )
  mouth_year <- mouth_outlet |>
    dplyr::inner_join(weather, by = "code_station", relationship = "many-to-many") |>
    dplyr::summarise(
      dplyr::across(c(precipitation, temperature), \(x) stats::weighted.mean(x, area)),
      .by = c(mouth_id, year)
    ) |>
    dplyr::arrange(mouth_id, year) |>
    dplyr::mutate(
      precipitation_anomaly = log(precipitation / mean(precipitation)),
      precipitation_previous = dplyr::lag(precipitation_anomaly),
      temperature_anomaly = temperature - mean(temperature),
      .by = mouth_id
    ) |>
    dplyr::select(mouth_id, year, precipitation_anomaly, precipitation_previous, temperature_anomaly)
  list(mouth = mouth, mouth_year = mouth_year)
}

#' Guild counts at the sea operations, with their drivers.
#'
#' One row per operation and guild: number of species of the guild (`n`)
#' among the species with a guild caught (`total`), for the operations
#' linked to a main river, with salinity, a catch and the weather of the
#' year. Drivers are centred and scaled: `salinity` (and `salinity2`, square
#' of the scaled salinity), `arable`, `artificial`, `intensive` (arable +
#' artificial: the two are correlated between the catchments), `log_area`,
#' `precipitation_anomaly`, `temperature_anomaly`.
#'
#' @param operation tibble, e.g. the `sea_operation` target.
#' @param operation_mouth tibble, [link_operation_mouth()]'s output.
#' @param catchment list, [summarise_mouth_catchment()]'s output.
#' @param catch tibble with `operation_id`, `species`, `guild`, e.g. the
#'   `sea_catch` target.
#' @return tibble.
#' @export
prepare_guild_visit <- function(operation, operation_mouth, catchment, catch) {
  count <- catch |>
    dplyr::filter(!is.na(guild)) |>
    dplyr::count(operation_id, guild) |>
    tidyr::complete(operation_id, guild, fill = list(n = 0)) |>
    dplyr::mutate(total = sum(n), .by = operation_id)
  operation |>
    dplyr::select(operation_id, year, month, longitude, latitude, survey, gear, salinity) |>
    dplyr::inner_join(dplyr::select(operation_mouth, operation_id, mouth_id), by = "operation_id") |>
    dplyr::inner_join(catchment$mouth, by = "mouth_id") |>
    dplyr::inner_join(catchment$mouth_year, by = c("mouth_id", "year")) |>
    dplyr::filter(!is.na(salinity), !is.na(precipitation_anomaly)) |>
    dplyr::mutate(
      log_area = log(area),
      intensive = arable + artificial,
      dplyr::across(
        c(salinity, arable, artificial, intensive, log_area, precipitation_anomaly, temperature_anomaly),
        \(x) as.numeric(scale(x))
      ),
      salinity2 = salinity^2
    ) |>
    dplyr::inner_join(count, by = "operation_id")
}

#' Mesh of the spatial field of the guild models.
#'
#' Fine triangles (edges up to 10 km) over the operations, for structure
#' within estuaries and bays; larger ones (up to 50 km) in the outer buffer.
#'
#' @param visit tibble with `longitude`, `latitude`.
#' @return INLA mesh, coordinates in km (Lambert-93).
#' @export
build_guild_mesh <- function(visit) {
  INLA::inla.mesh.2d(loc = unique(visit_coordinates(visit)), max.edge = c(10, 50), cutoff = 2)
}

#' Fit the model of the share of a guild among the species of an operation.
#'
#' Binomial model with logit link of the number of species of `guild`
#' among the species caught, on `term`, sampling (gear and month as fixed
#' effects), an iid catchment (main river) effect, an AR1 year effect and a
#' Matern spatial field (SPDE on `mesh`; PC priors: range below 20 km with
#' probability 0.5, so that the field captures structure within estuaries
#' rather than between catchments).
#'
#' @param visit tibble, [prepare_guild_visit()]'s output.
#' @param mesh INLA mesh, [build_guild_mesh()]'s output.
#' @param guild character, guild code.
#' @param term character, right-hand side of the fixed effects.
#' @return list with `guild`, `term`, `n` (operations), `summary.fixed`,
#'   `summary.hyperpar`, `summary.mouth` (catchment effects), `summary.year`,
#'   `waic` and `seconds` (fitting time).
#' @export
fit_guild_model <- function(visit, mesh, guild, term) {
  start <- Sys.time()
  data <- visit[visit$guild == guild, ] |>
    dplyr::mutate(gear = factor(gear), month = factor(month), year_index = year - min(year) + 1L)
  sampling <- stats::model.matrix(~ gear + month, data)[, -1, drop = FALSE]
  colnames(sampling) <- make.names(colnames(sampling))
  covariate <- dplyr::bind_cols(
    dplyr::select(data, dplyr::all_of(all.vars(stats::as.formula(paste("~", term)))), mouth_id, year_index),
    tibble::as_tibble(sampling)
  )
  spde <- INLA::inla.spde2.pcmatern(mesh, prior.range = c(20, 0.5), prior.sigma = c(1, 0.01))
  stack <- INLA::inla.stack(
    data = list(n = data$n, total = data$total),
    A = list(INLA::inla.spde.make.A(mesh, loc = visit_coordinates(data)), 1),
    effects = list(list(spatial = seq_len(spde$n.spde), Intercept = 1), covariate)
  )
  stack_data <- INLA::inla.stack.data(stack)
  model <- INLA::inla(
    stats::as.formula(paste(
      "n ~ -1 + Intercept +", term, "+", paste(colnames(sampling), collapse = " + "), "+",
      "f(mouth_id, model = 'iid', hyper = list(prec = list(initial = 0, prior = 'pc.prec', param = c(1, 0.01)))) +",
      "f(year_index, model = 'ar1', hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01)))) +",
      "f(spatial, model = spde)"
    )),
    family = "binomial", Ntrials = stack_data$total, data = stack_data, num.threads = 4,
    quantiles = c(0.045, 0.5, 0.955),
    control.predictor = list(A = INLA::inla.stack.A(stack)),
    control.compute = list(waic = TRUE)
  )
  list(
    guild = guild, term = term, n = nrow(data),
    summary.fixed = model$summary.fixed, summary.hyperpar = model$summary.hyperpar,
    summary.mouth = model$summary.random$mouth_id, summary.year = model$summary.random$year_index,
    waic = model$waic$waic, seconds = as.numeric(difftime(Sys.time(), start, units = "secs"))
  )
}
