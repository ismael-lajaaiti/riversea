#' Download the AMOBIO data base
#'
#' From data.gouv at the link
#' https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/EQYVLP
#'
#' @param dir where the data is downloaded.
#'
#' @return nothing
#' @export
download_amobio_data <- function(dir, verbose = FALSE) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  Sys.setenv("DATAVERSE_SERVER" = "entrepot.recherche.data.gouv.fr")
  dataset <- dataverse::get_dataset("doi:10.57745/EQYVLP")
  files <- dataset$files
  paths <- purrr::map2_chr(files$id, files$label, \(id, filename) {
    if (verbose) {
      print(paste("Downloading", filename))
    }
    path <- file.path(dir, filename)
    writeBin(dataverse::get_file(id, progress = verbose), path)
    path
  })
  paths
}

extract_metrics_amobio <- function(paths, average_time) {
  metric_file <- paths[grepl("METRICS_AMOBIO.Rdata", paths)]
  df_metric <- get(base::load(metric_file))
  pc_x_cols <- names(df_metric) |> stringr::str_subset("^PC_X")
  df_metric <- df_metric |>
    dplyr::select(-dplyr::any_of(pc_x_cols)) |>
    dplyr::rename_with(tolower) |>
    dplyr::select(
      node_id, id, date_operation,
      dplyr::any_of(amobio_selected_vars(average_time))
    ) |>
    dplyr::rename(dplyr::any_of(rename_amobio(average_time))) |>
    dplyr::mutate(
      node_id = as.integer(as.character(node_id)),
      sandre_code = stringr::str_pad(as.character(sandre_code), 8, pad = "0")
    )
}

extract_nodes_amobio <- function(paths) {
  network_file <- paths[grepl("NETWORK_AMOBIO.Rdata", paths)]
  network <- get(load(network_file))
  network |>
    sfnetworks::activate("nodes") |>
    tibble::as_tibble() |>
    tibble::as_tibble() |>
    dplyr::filter(IS_FISH == "FISH") |>
    dplyr::select(node_id, geometry, ELEVATION_IGN25M) |>
    dplyr::rename_with(tolower) |>
    dplyr::rename(dplyr::any_of(rename_amobio()))
}

extract_river_mouth_amobio <- function(paths) {
  network_file <- paths[grepl("NETWORK_AMOBIO.Rdata", paths)]
  network <- get(base::load(network_file))
  network |>
    sfnetworks::activate("nodes") |>
    tibble::as_tibble() |>
    tibble::as_tibble() |>
    dplyr::filter(IS_RIVERMOUTH_LEVEL1 == "RIVERMOUTH_LEVEL1") |>
    dplyr::select(node_id, geometry)
}

match_mouth_district <- function(river_mouth, district, dist_max, district_kept) {
  hydro_kept <- district |>
    filter(district %in% district_kept)

  river_mouth_sf <- river_mouth |>
    sf::st_as_sf() |>
    sf::st_transform(4326)

  joined <- river_mouth_sf |>
    sf::st_join(hydro_kept["district"], join = sf::st_within)
  matched <- joined |>
    dplyr::filter(!is.na(district)) |>
    dplyr::mutate(dist_km = 0)
  unmatched <- joined |> dplyr::filter(is.na(district))

  nearest_idx <- sf::st_nearest_feature(unmatched, hydro_kept)
  unmatched <- unmatched |>
    dplyr::mutate(
      district = hydro_kept$district[nearest_idx],
      dist_km = as.numeric(units::set_units(
        sf::st_distance(
          unmatched,
          hydro_kept[nearest_idx, ],
          by_element = TRUE
        ),
        "km"
      ))
    )

  river_mouth_district <- dplyr::bind_rows(matched, unmatched) |>
    dplyr::mutate(matched = dist_km <= dist_max)
}

#' Load the AMOBIO river network topology, stripped down to the minimum
#'
#' The raw network is ~1GB on disk / ~2.9GB in memory (1M nodes, 1.1M edges,
#' ~80 columns each - mostly links to unrelated monitoring networks, plus
#' full edge geometries). Keeps only what a shortest-path distance
#' computation needs: node geometry and the river-mouth flag, and edge
#' topology with the precomputed length (dropping edge geometry, which is
#' most of the memory).
#'
#' @param paths character vector of AMOBIO file paths.
#'
#' @return list(nodes, edges): `nodes` has `node_id`, `geometry`,
#'   `is_river_mouth`; `edges` has `from`, `to`, `length` (meters).
#' @export
extract_river_network <- function(paths) {
  network_file <- paths[grepl("NETWORK_AMOBIO.Rdata", paths)]
  network <- get(base::load(network_file))

  nodes <- network |>
    sfnetworks::activate("nodes") |>
    tibble::as_tibble() |>
    tibble::as_tibble() |>
    dplyr::transmute(
      node_id,
      geometry,
      is_river_mouth = IS_RIVERMOUTH_LEVEL1 %in% "RIVERMOUTH_LEVEL1"
    )

  edges <- network |>
    sfnetworks::activate("edges") |>
    tibble::as_tibble() |>
    tibble::as_tibble() |>
    dplyr::select(from, to, length = length_meters_num)

  list(nodes = nodes, edges = edges)
}

#' Patch a known missing link in the AMOBIO river network.
#'
#' A ~6,000-node component around the Vendée/Deux-Sèvres area is
#' topologically disconnected from the rest of the network - no node in it
#' is flagged a river mouth, and it has no path to any node outside it -
#' despite the closest pair of nodes across the gap being only 162m apart,
#' right at the Sèvre Nantaise/Loire confluence in Nantes. That's a
#' missing/broken digitized link, not a real hydrological separation (see
#' the D8 notebook for the full investigation), and adding this one edge
#' reconnects the entire component (verified: all of its nodes end up in
#' the same component as the rest of the network afterwards).
#'
#' @param network list(nodes, edges), see `extract_river_network()`.
#'
#' @return `network` with one additional edge added to `edges`, its
#' `length` the straight-line distance between the two nodes.
#' @export
patch_river_network <- function(network) {
  from_id <- 614458L
  to_id <- 614552L
  endpoints <- network$nodes |>
    dplyr::filter(node_id %in% c(from_id, to_id)) |>
    sf::st_as_sf()
  patch_length <- as.numeric(sf::st_distance(
    endpoints[endpoints$node_id == from_id, ],
    endpoints[endpoints$node_id == to_id, ]
  ))

  patch_edge <- tibble::tibble(from = from_id, to = to_id, length = patch_length)
  network$edges <- dplyr::bind_rows(network$edges, patch_edge)
  network
}

#' Restrict the AMOBIO network to within `max_dist` of the kept river mouths
#'
#' Grows the network outward (channel length, not straight-line) from every
#' matched river mouth at once, using a single Dijkstra run from a virtual
#' source connected to all of them with zero-weight edges - much cheaper than
#' one shortest-path search per mouth.
#'
#' @param network list(nodes, edges), see [extract_river_network()].
#' @param river_mouth tibble with `node_id`, `matched`, e.g. from
#'   [match_mouth_district()].
#' @param max_dist maximum channel-length distance to keep, in meters.
#'
#' @return list(nodes, edges) restricted to nodes within `max_dist`.
#' @export
restrict_river_network <- function(network, river_mouth, max_dist) {
  g <- igraph::graph_from_edgelist(
    as.matrix(network$edges[c("from", "to")]),
    directed = FALSE
  )
  igraph::E(g)$weight <- network$edges$length

  seeds <- river_mouth$node_id[river_mouth$matched]
  source <- igraph::vcount(g) + 1L
  g <- igraph::add_vertices(g, 1)
  g <- igraph::add_edges(g, as.vector(rbind(source, as.integer(seeds))), weight = 0)

  dist <- igraph::distances(g, v = source, weights = igraph::E(g)$weight)[1, ]
  kept_id <- setdiff(which(dist <= max_dist), source)

  list(
    nodes = dplyr::filter(network$nodes, node_id %in% kept_id),
    edges = dplyr::filter(network$edges, from %in% kept_id & to %in% kept_id)
  )
}

#' Snap operations to the nearest AMOBIO river network edge.
#'
#' Matches against `restrict_river_network()`'s kept-catchment-only
#' output. `snap_dist_m` is the straight-line distance to the nearest point
#' on the nearest edge linestring, not to the nearest node - an operation
#' sitting midway along a straight stretch between two nodes isn't
#' penalized for it. `node_id` is still the nearest *node* (needed
#' downstream for network routing), so it can be a step or two farther
#' than `snap_dist_m` alone would suggest.
#'
#' @param operation tibble with `operation_id`, `longitude`, `latitude`
#'   (WGS84), `survey`, `district`, e.g. `classify_operation_location()`'s
#'   output, restricted to inland operations.
#' @param network list(nodes, edges) with edge linestring geometry
#'   (Lambert-93), e.g. `add_edge_geometry()`'s output.
#'
#' @return one row per distinct operation (`operation_id`, `longitude`,
#' `latitude`, `survey`, `district`), with `node_id` (nearest network
#' node) and `snap_dist_m` (straight-line distance to the nearest edge,
#' in meters) added.
#' @export
snap_operations <- function(operation, network) {
  nodes <- sf::st_as_sf(network$nodes)
  edges <- sf::st_as_sf(network$edges)
  op_sf <- operation |>
    dplyr::distinct(operation_id, longitude, latitude, survey, district) |>
    sf::st_as_sf(
      coords = c("longitude", "latitude"), crs = 4326, remove = FALSE
    ) |>
    sf::st_transform(sf::st_crs(nodes))

  nearest_node <- sf::st_nearest_feature(op_sf, nodes)
  nearest_edge <- sf::st_nearest_feature(op_sf, edges)
  op_sf |>
    dplyr::mutate(
      node_id = nodes$node_id[nearest_node],
      snap_dist_m = as.numeric(
        sf::st_distance(op_sf, edges[nearest_edge, ], by_element = TRUE)
      )
    ) |>
    sf::st_drop_geometry()
}

#' Compute the channel-length distance from each operation to the nearest
#' river mouth of its own catchment.
#'
#' Runs a separate multi-source Dijkstra per district (seeded from that
#' district's own matched river mouths only), so an operation's distance
#' is always to a mouth of its own catchment. `dist_mouth_m` includes the
#' operation's own offset from its snapped node.
#'
#' @param network list(nodes, edges) with edge `length`, e.g.
#'   `restrict_river_network()`'s output.
#' @param operation tibble with `node_id`, `longitude`, `latitude`,
#'   `district`, e.g. `snap_operations()`'s output.
#' @param river_mouth tibble with `node_id`, `district`, `matched`, e.g.
#'   `match_mouth_district()`'s output.
#' @param district_kept character vector of districts to compute for.
#'
#' @return `operation`, restricted to `district_kept`, with `dist_mouth_m`
#' added (channel-length distance to the nearest same-district river
#' mouth, in meters).
#' @export
compute_distance_to_mouth <- function(network,
                                       operation,
                                       river_mouth,
                                       district_kept) {
  operation_kept <- operation |> dplyr::filter(district %in% district_kept)
  nodes <- sf::st_as_sf(network$nodes)
  operation_sf <- operation_kept |>
    sf::st_as_sf(
      coords = c("longitude", "latitude"), crs = 4326, remove = FALSE
    ) |>
    sf::st_transform(sf::st_crs(nodes))
  node_dist_m <- as.numeric(sf::st_distance(
    operation_sf, nodes[match(operation_kept$node_id, nodes$node_id), ],
    by_element = TRUE
  ))

  purrr::map_dfr(district_kept, function(d) {
    seeds <- river_mouth$node_id[
      river_mouth$matched & river_mouth$district == d
    ]

    g <- igraph::graph_from_edgelist(
      as.matrix(network$edges[c("from", "to")]),
      directed = FALSE
    )
    igraph::E(g)$weight <- network$edges$length
    source <- igraph::vcount(g) + 1L
    g <- igraph::add_vertices(g, 1)
    g <- igraph::add_edges(
      g, as.vector(rbind(source, as.integer(seeds))),
      weight = 0
    )

    dist <- igraph::distances(g, v = source, weights = igraph::E(g)$weight)[1, ]

    is_d <- operation_kept$district == d
    operation_kept[is_d, ] |>
      dplyr::mutate(dist_mouth_m = dist[node_id] + node_dist_m[is_d])
  })
}

#' Channel-network path from one operation to its nearest same-district
#' river mouth.
#'
#' Same idea as `sea_mouth_path()`, but for the AMOBIO river network: the
#' path is returned as the actual channel geometry of the edges it walks
#' (not straight lines between nodes), since river channels wind.
#'
#' @param op_row single-row tibble with `node_id`, `district`, e.g. one
#'   row of `compute_distance_to_mouth()`'s output.
#' @param network list(nodes, edges) with edge `length` and `geometry`,
#'   e.g. `add_edge_geometry()`'s output.
#' @param river_mouth tibble with `node_id`, `district`, `matched`, e.g.
#'   `match_mouth_district()`'s output.
#'
#' @return list(path, mouth_coords): `path` an sf linestring (WGS84),
#' `mouth_coords` the reached mouth's coordinates (WGS84).
#' @export
river_mouth_path <- function(op_row, network, river_mouth) {
  seeds <- river_mouth$node_id[
    river_mouth$matched & river_mouth$district == op_row$district
  ]
  g <- igraph::graph_from_edgelist(
    as.matrix(sf::st_drop_geometry(network$edges)[c("from", "to")]),
    directed = FALSE
  )
  igraph::E(g)$weight <- network$edges$length
  source <- igraph::vcount(g) + 1L
  g <- igraph::add_vertices(g, 1)
  g <- igraph::add_edges(
    g, as.vector(rbind(source, as.integer(seeds))),
    weight = 0
  )
  sp <- igraph::shortest_paths(
    g,
    from = op_row$node_id, to = source,
    weights = igraph::E(g)$weight, output = "both"
  )

  path_edges <- as.integer(sp$epath[[1]])
  real_edges <- path_edges[path_edges <= nrow(network$edges)]
  path_geom <- sf::st_transform(
    sf::st_union(sf::st_geometry(network$edges)[real_edges]), 4326
  )

  path_nodes <- as.integer(sp$vpath[[1]])
  path_nodes <- path_nodes[path_nodes != source]
  mouth_node_id <- utils::tail(path_nodes, 1)
  nodes_sf <- sf::st_as_sf(network$nodes) |> sf::st_transform(4326)
  mouth_coords <- sf::st_coordinates(
    nodes_sf[nodes_sf$node_id == mouth_node_id, ]
  )[1, ]

  list(
    path = sf::st_sf(geometry = path_geom),
    mouth_coords = mouth_coords
  )
}

#' Plot one river-network path checkpoint: the shortest channel-network
#' path highlighted against the rest of the local hydrographic network.
#'
#' @param op_row single-row tibble with `longitude`, `latitude`.
#' @param path list(path, mouth_coords), e.g. `river_mouth_path()`'s
#'   output.
#' @param network list(nodes, edges) with edge geometry (any CRS), e.g.
#'   `add_edge_geometry()`'s output - used to draw the rest of the local
#'   network as context.
#' @param basin_geom sf polygon(s), the catchment(s) to draw as
#'   geographic context, e.g. `basin` filtered to `op_row`'s district.
#' @param label panel label (e.g. `"A"`).
#'
#' @return ggplot.
#' @export
plot_river_path_checkpoint <- function(op_row, path, network, basin_geom, label) {
  op_pt <- sf::st_as_sf(
    op_row,
    coords = c("longitude", "latitude"), crs = 4326, remove = FALSE
  )
  mouth_pt <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_point(path$mouth_coords), crs = 4326
  ))

  bbox <- sf::st_bbox(path$path)
  pad_x <- max((bbox["xmax"] - bbox["xmin"]) * 0.6, 0.03)
  pad_y <- max((bbox["ymax"] - bbox["ymin"]) * 0.6, 0.03)
  xlim <- c(bbox["xmin"] - pad_x, bbox["xmax"] + pad_x)
  ylim <- c(bbox["ymin"] - pad_y, bbox["ymax"] + pad_y)

  zoom_box <- sf::st_transform(
    sf::st_as_sfc(sf::st_bbox(
      c(
        xmin = unname(xlim[1]), xmax = unname(xlim[2]),
        ymin = unname(ylim[1]), ymax = unname(ylim[2])
      ),
      crs = 4326
    )),
    sf::st_crs(network$edges)
  )
  context_edges <- sf::st_transform(
    sf::st_filter(network$edges, zoom_box), 4326
  )

  ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = basin_geom, fill = "grey97", color = "grey40", linewidth = 0.4
    ) +
    ggplot2::geom_sf(data = context_edges, color = "grey65", linewidth = 0.2) +
    ggplot2::geom_sf(data = path$path, color = "steelblue", linewidth = 1.4) +
    ggplot2::geom_sf(data = op_pt, color = "black", size = 3.5, shape = 17) +
    ggplot2::geom_sf(data = mouth_pt, color = "black", size = 3.5, shape = 8) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
    ggplot2::coord_sf(xlim = xlim, ylim = ylim) +
    ggplot2::labs(title = label, x = NULL, y = NULL) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 5, r = 14, b = 5, l = 5)
    )
}

#' Build a water-only triangular mesh over the sea, excluding land.
#'
#' Domain is a buffered bounding box around the sea-side operations and
#' matched river mouths, minus `land`. `fmesher::fm_mesh_2d()` auto-pads
#' the mesh beyond this boundary regardless of land - `in_water` flags
#' which triangles actually lie in water.
#'
#' @param operation tibble with `longitude`, `latitude`, `survey`,
#'   `inland`, `district`, e.g. `classify_operation_location()`'s output.
#' @param mouth sf tibble with `district`, `matched`, `geometry`
#'   (WGS84), e.g. `match_mouth_district()`'s output.
#' @param land sf polygons (WGS84), e.g. `land_polygon()`'s output.
#' @param district_kept character vector of districts to build the mesh
#'   for.
#' @param buffer_deg degrees of padding around the sea-side operations
#'   and matched mouths.
#'
#' @return list(mesh, in_water, land): `mesh` is an `fm_mesh_2d` object,
#' `in_water` a logical vector (one per mesh triangle), `land` the
#' (cropped, simplified) polygon used as the boundary.
#' @export
build_sea_mesh <- function(operation, mouth, land, district_kept, buffer_deg) {
  op_sea <- operation |>
    dplyr::filter(survey != "river", !inland, district %in% district_kept)
  mouth_kept <- mouth |>
    dplyr::filter(matched, district %in% district_kept) |>
    sf::st_as_sf() |>
    sf::st_transform(4326)
  mouth_coords <- sf::st_coordinates(mouth_kept)

  domain_bbox <- sf::st_bbox(c(
    xmin = min(op_sea$longitude, mouth_coords[, 1]) - buffer_deg,
    xmax = max(op_sea$longitude, mouth_coords[, 1]) + buffer_deg,
    ymin = min(op_sea$latitude, mouth_coords[, 2]) - buffer_deg,
    ymax = max(op_sea$latitude, mouth_coords[, 2]) + buffer_deg
  ), crs = 4326)
  domain_box_sfc <- sf::st_as_sfc(domain_bbox)

  land_crop <- sf::st_intersection(sf::st_make_valid(land), domain_box_sfc)
  land_union <- sf::st_make_valid(sf::st_union(land_crop))
  land_simple <- sf::st_simplify(
    land_union,
    dTolerance = 0.005, preserveTopology = TRUE
  )
  land_parts <- sf::st_cast(land_simple, "POLYGON", warn = FALSE)
  land_parts <- land_parts[as.numeric(sf::st_area(land_parts)) > 1e5]

  water_poly <- sf::st_difference(domain_box_sfc, sf::st_union(land_parts))
  water_seg <- fmesher::fm_as_segm(water_poly)
  mesh <- fmesher::fm_mesh_2d(
    boundary = water_seg,
    max.edge = c(0.05, 0.3),
    cutoff = 0.01
  )

  tv <- mesh$graph$tv
  loc <- mesh$loc[, 1:2]
  centroids <- cbind(
    (loc[tv[, 1], 1] + loc[tv[, 2], 1] + loc[tv[, 3], 1]) / 3,
    (loc[tv[, 1], 2] + loc[tv[, 2], 2] + loc[tv[, 3], 2]) / 3
  )
  centroids_sf <- sf::st_as_sf(
    as.data.frame(centroids),
    coords = c("V1", "V2"), crs = 4326
  )
  in_water <- lengths(sf::st_intersects(centroids_sf, water_poly)) > 0

  list(mesh = mesh, in_water = in_water, land = land_simple)
}

#' Build a graph from a sea mesh's water triangles.
#'
#' Vertices keep `x`/`y` (WGS84) as attributes. Edges are water-triangle
#' sides, weighted by geodesic distance. Vertices only used by
#' land-extension triangles (see `build_sea_mesh()`) stay isolated.
#'
#' @param sea_mesh list(mesh, in_water, land), e.g. `build_sea_mesh()`'s
#'   output.
#'
#' @return igraph graph, undirected, weighted (`weight`, meters).
#' @export
build_sea_graph <- function(sea_mesh) {
  tv <- sea_mesh$mesh$graph$tv[sea_mesh$in_water, , drop = FALSE]
  loc <- sea_mesh$mesh$loc[, 1:2]

  edges <- rbind(tv[, c(1, 2)], tv[, c(2, 3)], tv[, c(3, 1)])
  edges <- t(apply(edges, 1, sort))
  edges <- unique(edges)

  p1 <- sf::st_sfc(
    lapply(seq_len(nrow(edges)), \(i) sf::st_point(loc[edges[i, 1], ])),
    crs = 4326
  )
  p2 <- sf::st_sfc(
    lapply(seq_len(nrow(edges)), \(i) sf::st_point(loc[edges[i, 2], ])),
    crs = 4326
  )
  w <- as.numeric(sf::st_distance(p1, p2, by_element = TRUE))

  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  igraph::E(g)$weight <- w
  igraph::V(g)$x <- loc[seq_len(igraph::vcount(g)), 1]
  igraph::V(g)$y <- loc[seq_len(igraph::vcount(g)), 2]
  g
}

#' Nearest reachable node of a sea graph, for each point.
#'
#' @param points sf points (WGS84).
#' @param sea_graph igraph graph with vertex attributes `x`/`y`, e.g.
#'   `build_sea_graph()`'s output.
#'
#' @return list(node_id, offset_m): `node_id` the nearest reachable
#' vertex id, `offset_m` the distance from each point to it (meters).
#' @export
snap_to_sea_graph <- function(points, sea_graph) {
  reachable <- as.integer(igraph::V(sea_graph)[igraph::degree(sea_graph) > 0])
  reachable_sf <- sf::st_as_sf(
    data.frame(
      node_id = reachable,
      x = igraph::V(sea_graph)$x[reachable],
      y = igraph::V(sea_graph)$y[reachable]
    ),
    coords = c("x", "y"), crs = 4326
  )
  nearest <- sf::st_nearest_feature(points, reachable_sf)
  list(
    node_id = reachable_sf$node_id[nearest],
    offset_m = as.numeric(sf::st_distance(
      points, reachable_sf[nearest, ],
      by_element = TRUE
    ))
  )
}

#' Compute the coast-aware distance from each sea-side operation to the
#' nearest river mouth of its own district.
#'
#' Snaps each operation and matched mouth to the nearest reachable node
#' of `sea_graph` (see `snap_to_sea_graph()`), then runs a separate
#' multi-source Dijkstra per district (seeded from that district's own
#' matched mouths only) - see `compute_distance_to_mouth()` for the same
#' pattern on the river side.
#'
#' @param operation tibble with `longitude`, `latitude`, `survey`,
#'   `inland`, `district`, e.g. `classify_operation_location()`'s output.
#' @param mouth sf tibble with `district`, `matched`, `geometry`
#'   (WGS84), e.g. `match_mouth_district()`'s output.
#' @param sea_graph igraph graph with vertex attributes `x`/`y`, e.g.
#'   `build_sea_graph()`'s output.
#' @param district_kept character vector of districts to compute for.
#'
#' @return sea-side rows of `operation`, restricted to `district_kept`,
#' with `dist_mouth_m` added (meters).
#' @export
compute_sea_distance_to_mouth <- function(operation,
                                           mouth,
                                           sea_graph,
                                           district_kept) {
  op_sea <- operation |>
    dplyr::filter(survey != "river", !inland, district %in% district_kept)
  op_sf <- op_sea |>
    sf::st_as_sf(
      coords = c("longitude", "latitude"), crs = 4326, remove = FALSE
    )
  snap_op <- snap_to_sea_graph(op_sf, sea_graph)
  op_sea$node_id <- snap_op$node_id
  node_dist_m <- snap_op$offset_m

  mouth_kept <- mouth |>
    dplyr::filter(matched, district %in% district_kept) |>
    sf::st_as_sf() |>
    sf::st_transform(4326)
  mouth_kept$node_id <- snap_to_sea_graph(mouth_kept, sea_graph)$node_id

  purrr::map_dfr(district_kept, function(d) {
    seeds <- mouth_kept$node_id[mouth_kept$district == d]
    source <- igraph::vcount(sea_graph) + 1L
    g2 <- igraph::add_vertices(sea_graph, 1)
    g2 <- igraph::add_edges(
      g2, as.vector(rbind(source, as.integer(seeds))),
      weight = 0
    )
    dist <- igraph::distances(
      g2, v = source, weights = igraph::E(g2)$weight
    )[1, ]

    is_d <- op_sea$district == d
    op_sea[is_d, ] |>
      dplyr::mutate(dist_mouth_m = dist[node_id] + node_dist_m[is_d])
  })
}

#' Combine river-network and sea-mesh distances into a single dataframe.
#'
#' @param inland_distance tibble from `compute_distance_to_mouth()`.
#' @param sea_offshore_distance tibble from
#'   `compute_sea_distance_to_mouth()`.
#'
#' @return tibble with `operation_id`, `longitude`, `latitude`, `survey`,
#' `district`, `dist_mouth_m`, `inland` (TRUE if `dist_mouth_m` was
#' computed via the river network, FALSE if via the sea mesh).
#' @export
combine_distance_to_mouth <- function(inland_distance, sea_offshore_distance) {
  cols <- c(
    "operation_id", "longitude", "latitude", "survey", "district",
    "dist_mouth_m", "inland"
  )
  dplyr::bind_rows(
    inland_distance |>
      dplyr::mutate(inland = TRUE) |>
      dplyr::select(dplyr::all_of(cols)),
    sea_offshore_distance |>
      dplyr::mutate(inland = FALSE) |>
      dplyr::select(dplyr::all_of(cols))
  )
}

#' Water-triangle edges of a sea mesh, as sf linestrings.
#'
#' @param sea_mesh list(mesh, in_water, land), e.g. `build_sea_mesh()`'s
#'   output.
#'
#' @return sf linestrings (WGS84), one per water-triangle edge.
#' @export
sea_mesh_edges <- function(sea_mesh) {
  tv <- sea_mesh$mesh$graph$tv[sea_mesh$in_water, , drop = FALSE]
  loc <- sea_mesh$mesh$loc[, 1:2]
  edges <- rbind(tv[, c(1, 2)], tv[, c(2, 3)], tv[, c(3, 1)])
  edges <- t(apply(edges, 1, sort))
  edges <- unique(edges)
  seg <- lapply(
    seq_len(nrow(edges)), \(i) sf::st_linestring(loc[edges[i, ], ])
  )
  sf::st_sf(geometry = sf::st_sfc(seg, crs = 4326))
}

#' Plot a sea mesh, zoomed to a bounding box.
#'
#' @param sea_mesh list(mesh, in_water, land), e.g. `build_sea_mesh()`'s
#'   output.
#' @param mesh_edges sf linestrings, e.g. `sea_mesh_edges()`'s output.
#' @param xlim,ylim numeric length-2, plot extent (degrees).
#' @param title plot title.
#'
#' @return ggplot.
#' @export
plot_sea_mesh <- function(sea_mesh, mesh_edges, xlim, ylim, title = NULL) {
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = sea_mesh$land, fill = "grey75", color = NA) +
    ggplot2::geom_sf(data = mesh_edges, color = "steelblue", linewidth = 0.15) +
    ggplot2::coord_sf(xlim = xlim, ylim = ylim) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' Coast-aware path from one operation to its nearest same-district
#' river mouth.
#'
#' @param op_row single-row tibble with `node_id`, `district`, e.g. one
#'   row of `compute_sea_distance_to_mouth()`'s output.
#' @param mouth sf tibble with `district`, `node_id` (already snapped to
#'   `sea_graph`, see `snap_to_sea_graph()`).
#' @param sea_graph igraph graph with vertex attributes `x`/`y`, e.g.
#'   `build_sea_graph()`'s output.
#'
#' @return list(path, mouth_coords): `path` an sf linestring (WGS84),
#' `mouth_coords` the reached mouth's coordinates.
#' @export
sea_mouth_path <- function(op_row, mouth, sea_graph) {
  seeds <- mouth$node_id[mouth$district == op_row$district]
  source <- igraph::vcount(sea_graph) + 1L
  g2 <- igraph::add_vertices(sea_graph, 1)
  g2 <- igraph::add_edges(
    g2, as.vector(rbind(source, as.integer(seeds))),
    weight = 0
  )
  sp <- igraph::shortest_paths(
    g2,
    from = op_row$node_id, to = source,
    weights = igraph::E(g2)$weight, output = "vpath"
  )
  path_nodes <- as.integer(sp$vpath[[1]])
  path_nodes <- path_nodes[path_nodes != source]
  loc <- cbind(igraph::V(sea_graph)$x, igraph::V(sea_graph)$y)
  list(
    path = sf::st_sf(geometry = sf::st_sfc(
      sf::st_linestring(loc[path_nodes, , drop = FALSE]), crs = 4326
    )),
    mouth_coords = loc[utils::tail(path_nodes, 1), ]
  )
}

#' Plot one path checkpoint: straight line vs. coast-aware mesh path.
#'
#' @param op_row single-row tibble with `longitude`, `latitude`,
#'   `dist_mouth_m`.
#' @param path list(path, mouth_coords), e.g. `sea_mouth_path()`'s
#'   output.
#' @param mesh_edges sf linestrings, e.g. `sea_mesh_edges()`'s output.
#' @param land sf polygons, e.g. `build_sea_mesh()`'s `land` output.
#' @param title plot title.
#' @param lang `"en"` or `"fr"`, language of the subtitle.
#' @param show_subtitle whether to add a subtitle with the straight-line
#'   and mesh-path distances.
#'
#' @return ggplot.
#' @export
plot_sea_path_checkpoint <- function(op_row, path, mesh_edges, land, title,
                                      lang = c("en", "fr"), show_subtitle = TRUE) {
  lang <- match.arg(lang)
  op_pt <- sf::st_as_sf(
    op_row,
    coords = c("longitude", "latitude"), crs = 4326, remove = FALSE
  )
  mouth_pt <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_point(path$mouth_coords), crs = 4326
  ))
  straight <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_linestring(rbind(
      c(op_row$longitude, op_row$latitude), path$mouth_coords
    )),
    crs = 4326
  ))
  straight_dist_m <- as.numeric(sf::st_distance(op_pt, mouth_pt))

  bbox <- sf::st_bbox(path$path)
  pad_x <- max((bbox["xmax"] - bbox["xmin"]) * 0.6, 0.03)
  pad_y <- max((bbox["ymax"] - bbox["ymin"]) * 0.6, 0.03)

  subtitle <- if (!show_subtitle) {
    NULL
  } else if (lang == "fr") {
    sprintf(
      "ligne droite : %.0fm | maillage : %.0fm",
      straight_dist_m, op_row$dist_mouth_m
    )
  } else {
    sprintf(
      "straight-line: %.0fm | mesh path: %.0fm",
      straight_dist_m, op_row$dist_mouth_m
    )
  }

  ggplot2::ggplot() +
    ggplot2::geom_sf(data = land, fill = "grey80", color = NA) +
    ggplot2::geom_sf(data = mesh_edges, color = "grey75", linewidth = 0.15) +
    ggplot2::geom_sf(
      data = straight, color = "firebrick",
      linewidth = 0.8, linetype = "dashed"
    ) +
    ggplot2::geom_sf(data = path$path, color = "steelblue", linewidth = 1.4) +
    ggplot2::geom_sf(data = op_pt, color = "black", size = 3.5, shape = 17) +
    ggplot2::geom_sf(data = mouth_pt, color = "black", size = 3.5, shape = 8) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
    ggplot2::coord_sf(
      xlim = c(bbox["xmin"] - pad_x, bbox["xmax"] + pad_x),
      ylim = c(bbox["ymin"] - pad_y, bbox["ymax"] + pad_y)
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = 9),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 5, r = 14, b = 5, l = 5)
    )
}

#' Length of a straight line's overlap with land, per operation.
#'
#' For each operation, draws a straight line to the nearest mouth of its
#' own district and measures how much of it overlaps `land`. Only counts
#' past a minimal overlap length - `st_intersects()` alone also flags a
#' line that merely touches a boundary vertex, not a real crossing.
#'
#' @param operation sf tibble with `operation_id`, `district`, `geometry`
#'   (WGS84).
#' @param mouth sf tibble with `district`, `geometry` (WGS84).
#' @param land sf polygons (WGS84).
#' @param district_kept character vector of districts to compute for.
#'
#' @return tibble with `operation_id`, `crossing_len_m` (meters, 0 if no
#' meaningful overlap).
#' @export
straight_line_land_crossing <- function(operation, mouth, land, district_kept) {
  land_union <- sf::st_union(land)
  mouth_by_district <- split(mouth, mouth$district)

  purrr::map_dfr(district_kept, function(d) {
    op_d <- operation[operation$district == d, ]
    mouth_d <- mouth_by_district[[d]]
    idx <- sf::st_nearest_feature(op_d, mouth_d)
    lines <- mapply(
      \(p, q) sf::st_linestring(rbind(p, q)),
      sf::st_coordinates(op_d) |> split(seq_len(nrow(op_d))),
      sf::st_coordinates(mouth_d[idx, ]) |> split(seq_len(nrow(op_d))),
      SIMPLIFY = FALSE
    )
    lines_sfc <- sf::st_sfc(lines, crs = 4326)
    has_inter <- lengths(sf::st_intersects(lines_sfc, land)) > 0
    len_m <- sapply(seq_along(lines_sfc), function(i) {
      if (!has_inter[i]) {
        return(0)
      }
      ii <- suppressWarnings(sf::st_intersection(lines_sfc[i], land_union))
      if (length(ii) == 0) {
        return(0)
      }
      sum(as.numeric(sf::st_length(sf::st_cast(ii, "LINESTRING"))))
    })
    tibble::tibble(operation_id = op_d$operation_id, crossing_len_m = len_m)
  })
}

#' Attach edge linestring geometry to an already-restricted network.
#'
#' `extract_river_network()` drops edge geometry to keep the
#' Dijkstra restriction cheap - `restrict_river_network()`'s node-to-node
#' distance only needs the precomputed `length`, not the shape of the
#' channel between them. Snapping to the nearest point on the true channel
#' (rather than the nearest node) needs that shape back, so this re-reads
#' just the geometry for the edges already kept by the restriction,
#' instead of redoing it with geometry attached throughout.
#'
#' @param network list(nodes, edges) already restricted, e.g.
#'   `restrict_river_network()`'s output.
#' @param paths character vector of AMOBIO file paths.
#'
#' @return `network` with `edges` gaining a `geometry` column
#' (linestring, same CRS as `nodes`).
#' @export
add_edge_geometry <- function(network, paths) {
  network_file <- paths[grepl("NETWORK_AMOBIO.Rdata", paths)]
  raw_edges <- get(base::load(network_file)) |>
    sfnetworks::activate("edges") |>
    sf::st_as_sf() |>
    dplyr::select(from, to, geometry) |>
    dplyr::distinct(from, to, .keep_all = TRUE)

  edges <- network$edges |>
    dplyr::inner_join(raw_edges, by = c("from", "to")) |>
    sf::st_as_sf()

  list(nodes = network$nodes, edges = edges)
}

#' Map of France, operations colored by distance to river mouth.
#'
#' @param distance tibble with `longitude`, `latitude`, `dist_mouth_m`,
#'   e.g. `combine_distance_to_mouth()`'s output.
#' @param district_kept character vector of districts to keep.
#' @param basin sf tibble with `district`, `geometry`, e.g.
#'   `format_basin()`'s output.
#'
#' @return ggplot.
#' @export
plot_distance_to_mouth_map <- function(distance, district_kept, basin) {
  france <- rnaturalearth::ne_countries(
    geounit = "france", type = "map_units", scale = "medium",
    returnclass = "sf"
  )
  basin_kept <- basin |> dplyr::filter(district %in% district_kept)

  ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = france, fill = NA, color = "grey80", linewidth = 0.25
    ) +
    ggplot2::geom_sf(
      data = basin_kept, fill = "grey94", color = "grey20", linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = distance,
      ggplot2::aes(longitude, latitude, color = dist_mouth_m / 1000),
      size = 0.5, alpha = 0.8
    ) +
    ggplot2::coord_sf(xlim = c(-5.3, 8.3), ylim = c(41.2, 51.2), expand = FALSE) +
    ggplot2::scale_color_viridis_c(name = "Distance à l'embouchure (km)") +
    ggplot2::labs(x = "Longitude", y = "Latitude") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "top"
    )
}

#' Three-panel figure illustrating the distance-to-mouth computation: one
#' example path on the river network, one example path on the coastal
#' mesh, and a France-wide map of the result for every operation.
#'
#' @param river_network_geo list(nodes, edges) with edge geometry, e.g.
#'   `add_edge_geometry()`'s output.
#' @param river_mouth_district tibble with `node_id`, `district`,
#'   `matched`, e.g. `match_mouth_district()`'s output.
#' @param inland_distance_to_mouth tibble, e.g.
#'   `compute_distance_to_mouth()`'s output.
#' @param sea_graph igraph graph, e.g. `build_sea_graph()`'s output.
#' @param sea_mesh list(mesh, in_water, land), e.g. `build_sea_mesh()`'s
#'   output.
#' @param sea_offshore_distance_to_mouth tibble, e.g.
#'   `compute_sea_distance_to_mouth()`'s output.
#' @param distance_to_mouth tibble, e.g.
#'   `combine_distance_to_mouth()`'s output.
#' @param district_kept character vector of districts to keep.
#' @param basin sf tibble with `district`, `geometry`, e.g.
#'   `format_basin()`'s output.
#' @param river_operation_id operation_id of the river-network example.
#' @param sea_operation_id operation_id of the coastal-mesh example.
#'
#' @return patchwork object, 1 row x 3 columns.
#' @export
plot_distance_to_mouth_overview <- function(river_network_geo,
                                             river_mouth_district,
                                             inland_distance_to_mouth,
                                             sea_graph,
                                             sea_mesh,
                                             sea_offshore_distance_to_mouth,
                                             distance_to_mouth,
                                             district_kept,
                                             basin,
                                             river_operation_id,
                                             sea_operation_id) {
  nice_theme()

  op_river <- inland_distance_to_mouth |>
    dplyr::filter(operation_id == river_operation_id)
  path_river <- river_mouth_path(op_river, river_network_geo, river_mouth_district)
  basin_river <- basin |> dplyr::filter(district == op_river$district)
  panel_a <- plot_river_path_checkpoint(
    op_river, path_river, river_network_geo, basin_river, "A"
  )

  mouth_sea <- river_mouth_district |>
    dplyr::filter(matched, district %in% district_kept) |>
    sf::st_as_sf() |>
    sf::st_transform(4326)
  mouth_sea$node_id <- snap_to_sea_graph(mouth_sea, sea_graph)$node_id
  op_sea <- sea_offshore_distance_to_mouth |>
    dplyr::filter(operation_id == sea_operation_id)
  path_sea <- sea_mouth_path(op_sea, mouth_sea, sea_graph)
  mesh_edges <- sea_mesh_edges(sea_mesh)
  panel_b <- plot_sea_path_checkpoint(
    op_sea, path_sea, mesh_edges, sea_mesh$land,
    "B", show_subtitle = FALSE
  )

  panel_c <- plot_distance_to_mouth_map(distance_to_mouth, district_kept, basin) +
    ggplot2::labs(title = "C")

  ((panel_a / panel_b) | panel_c) +
    patchwork::plot_layout(widths = c(1, 1.5)) &
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 15, face = "bold"),
      axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 11)
    )
}

#' Deduplicate AMOBIO metrics on (sandre_code, date).
#'
#' Some nodes have several candidate physicochemical stations
#' (`METADATA_PC_CdStation`) for the same date, which is dropped upstream.
#' This leaves rows that share `(sandre_code, date)` but disagree on the
#' `pc_n_*` nutrient columns, so those are averaged across candidates while
#' every other column (identical across candidates) is kept as is.
#'
#' @param metrics tibble
#'
#' @return tibble with one row per (sandre_code, date)
dedup_amobio_metrics <- function(metrics) {
  pc_cols <- names(metrics) |> stringr::str_subset("^pc_n_")
  other_cols <- setdiff(names(metrics), c(pc_cols, "sandre_code", "date"))

  metrics |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(pc_cols),
        \(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
      ),
      dplyr::across(dplyr::all_of(other_cols), dplyr::first),
      .by = c(sandre_code, date)
    )
}

#' Combine the AMOBIO tibble in a single one.
#'
#' Also select variables of interest.
#'
#' @param metrics tibble
#' @param nodes tibble
#'
#' @return sf object CRS Lambert-93
#' @export
combine_amobio_data <- function(metrics, nodes) {
  inner_join(nodes, metrics, by = join_by(node_id)) |>
    sf::st_as_sf()
}

#' Utility function to visualise AMOBIO data.
#'
#' @param amobio_data dataframe
#' @param metric string specifying the metric name to plot
#' @param title string giving the title of the plot
#' @param subtitle string giving the subtitle of the plot
#' @param legend string giving the legend (color caption) of the plot
#' @param log_transform whether or not to log transform the metric
#' @param year string to select point of a given year
#'
#' @return plot
#' @export
plot_amobio_metric <- function(
  amobio_data,
  metric,
  title = NULL,
  subtitle = NULL,
  legend = NULL,
  f_transform = \(x) x,
  year = NA
) {
  data <- amobio_data |>
    select(all_of(metric), date) |>
    filter(!is.na(.data[[metric]])) |>
    distinct() |>
    mutate(to_plot = f_transform(.data[[metric]]))

  if (!is.na(year)) {
    data <- data |>
      filter(lubridate::year(date) == year)
  }

  coast <- rnaturalearth::ne_coastline(scale = "large", returnclass = "sf") |>
    sf::st_transform(sf::st_crs(data))
  boundaries <- sf::st_bbox(sf::st_as_sf(data))

  ggplot2::ggplot() +
    ggplot2::geom_sf(data = coast) +
    ggplot2::geom_sf(
      data = data,
      ggplot2::aes(color = to_plot),
      size = 0.5
    ) +
    ggplot2::coord_sf(
      xlim = boundaries[c("xmin", "xmax")],
      ylim = boundaries[c("ymin", "ymax")]
    ) +
    hrbrthemes::scale_color_flexoki_continuous() +
    ggplot2::labs(
      x = "Latitude",
      y = "Longitude",
      color = legend,
      title = title,
      subtitle = subtitle
    )
}

#' Join AMOBIO and river food web data.
#'
#' @param amobio dataframe.
#' @param river_foodweb river food web structure, as returned by
#' `measure_foodweb_structure()` (filtered to the river survey).
#' @param river_operation river operation metadata, as returned by
#' `get_river_operation()`, providing the `sandre_code`/`date` join keys
#' that `river_foodweb` itself doesn't carry.
#'
#' @return dataframe.
#' @export
join_amobio_aspe <- function(amobio, river_foodweb, river_operation) {
  aspe <- river_foodweb |>
    left_join(
      river_operation |> select(operation_id, sandre_code, date),
      by = join_by(operation_id)
    )
  # dplyr::inner_join() on an sf object currently errors (sf 1.1.1 / dplyr
  # 1.2.1 incompatibility), so geometry is carried through as plain x/y
  # columns across the join and reattached afterwards instead.
  amobio_df <- amobio |>
    mutate(
      geom_x = sf::st_coordinates(geometry)[, 1],
      geom_y = sf::st_coordinates(geometry)[, 2]
    ) |>
    sf::st_drop_geometry()
  inner_join(amobio_df, aspe, by = join_by(sandre_code, date)) |>
    sf::st_as_sf(coords = c("geom_x", "geom_y"), crs = sf::st_crs(amobio))
}

#' Environmental variables from the AMOBIO databsed selected for the analysis
#'
#' @param average_time string indicating the duration of the average for
#' chemical variabels. Possible values are "1y" (one year), "3m" (three months),
#' "5y" (five years) and "all" (all available points).
#'
#' @return vector of string
#' @export
amobio_selected_vars <- function(average_time) {
  c(
    "elevation",
    "poe_ddam_l3_dam",
    "poe_ddam_l2_dam",
    "poe_bh5_l6",
    "poe_bh5_l5",
    paste("pc_n_ptotal_s1", average_time, sep = "_"),
    paste("pc_n_orthophosp_s1", average_time, sep = "_"),
    paste("pc_n_no3_s1", average_time, sep = "_"),
    paste("pc_n_no2_s1", average_time, sep = "_"),
    paste("pc_n_nh4_s1", average_time, sep = "_"),
    paste("pc_n_nkj_s1", average_time, sep = "_")
  )
}

#' Rename selected environmental variables of the AMOBIO database.
#'
#' @return vector of strings.
#' @export
rename_amobio <- function(average_time = "1y") {
  c(
    "date" = "date_operation",
    "sandre_code" = "id",
    "elevation" = "elevation_ign25m",
    "dam_distance_upstream" = "poe_ddam_l2_dam",
    "dam_distance_downstream" = "poe_ddam_l3_dam",
    "barrier_downstream" = "poe_bh5_l5",
    "barrier_upstream" = "poe_bh5_l6",
    "total_phosphorus" = paste("pc_n_ptotal_s1", average_time, sep = "_"),
    "phosphate" = paste("pc_n_orthophosp_s1", average_time, sep = "_"),
    "ammonium" = paste("pc_n_nh4_s1", average_time, sep = "_"),
    "nitrite" = paste("pc_n_no2_s1", average_time, sep = "_"),
    "nitrate" = paste("pc_n_no3_s1", average_time, sep = "_"),
    "kjeldahl_nitrogen" = paste("pc_n_nkj_s1", average_time, sep = "_")
  )
}

#' Summarise completeness of environmental variables
#'
#' For each environmental variable, report the number of observations, the
#' number and share of missing values, and the range of observed (non-NA)
#' values. Useful to spot variables with poor coverage before modelling.
#'
#' @param data tibble containing the environmental variables as columns.
#' @param vars character vector of variable names to summarise. Defaults to
#' the AMOBIO environmental variables (see [rename_amobio()]).
#'
#' @return tibble with one row per variable.
#' @export
summarise_env_na <- function(
  data,
  vars = setdiff(names(rename_amobio()), c("date", "sandre_code"))
) {
  n_total <- nrow(data)
  purrr::map(vars, \(v) {
    x <- data[[v]]
    n_na <- sum(is.na(x))
    tibble::tibble(
      variable = v,
      n = n_total,
      n_na = n_na,
      pct_na = round(100 * n_na / n_total, 1),
      min = suppressWarnings(min(x, na.rm = TRUE)),
      median = suppressWarnings(stats::median(x, na.rm = TRUE)),
      max = suppressWarnings(max(x, na.rm = TRUE))
    )
  }) |>
    purrr::list_rbind()
}

#' Summarise temporal coverage of environmental variables
#'
#' For each year and environmental variable, count the number of
#' non-missing observations, out of the number of points sampled that year.
#' Useful to spot years or variables with poor sampling coverage.
#'
#' @param data tibble containing a `date` column and the environmental
#' variables as columns.
#' @param vars character vector of variable names to summarise. Defaults to
#' the AMOBIO environmental variables (see [rename_amobio()]).
#'
#' @return tibble with one row per (year, variable).
#' @export
summarise_env_temporal_coverage <- function(
  data,
  vars = setdiff(names(rename_amobio()), c("date", "sandre_code"))
) {
  data |>
    sf::st_drop_geometry() |>
    dplyr::mutate(year = lubridate::year(date)) |>
    dplyr::reframe(
      dplyr::across(dplyr::all_of(vars), \(x) sum(!is.na(x))),
      n = dplyr::n(),
      .by = year
    ) |>
    tidyr::pivot_longer(
      dplyr::all_of(vars),
      names_to = "variable",
      values_to = "n_obs"
    ) |>
    dplyr::mutate(
      pct_obs = 100 * n_obs / n,
      variable = factor(variable, levels = rev(vars))
    ) |>
    dplyr::arrange(year, variable)
}

#' Plot temporal coverage of environmental variables
#'
#' Heatmap of the share of non-missing observations for each environmental
#' variable, across years.
#'
#' @param data tibble containing a `date` column and the environmental
#' variables as columns.
#' @param vars character vector of variable names to summarise. Defaults to
#' the AMOBIO environmental variables (see [rename_amobio()]).
#'
#' @return ggplot
#' @export
plot_env_temporal_coverage <- function(
  data,
  vars = setdiff(names(rename_amobio()), c("date", "sandre_code"))
) {
  coverage <- summarise_env_temporal_coverage(data, vars)

  ggplot2::ggplot(
    coverage,
    ggplot2::aes(x = year, y = variable, fill = n_obs)
  ) +
    ggplot2::geom_tile(color = "#FAFAF8", linewidth = 0.6) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty()) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

#' Summarise extreme values in environmental variables
#'
#' Flags extreme values using Tukey's fences (values below
#' `Q1 - 1.5 * IQR` or above `Q3 + 1.5 * IQR` are flagged), computed on the
#' log10 scale by default since nutrient concentrations are typically
#' right-skewed. Non-positive values are excluded from the log-scale fences
#' (reported separately) since they cannot be log-transformed.
#'
#' @param data tibble containing the environmental variables as columns.
#' @param vars character vector of variable names to check. Defaults to the
#' nutrient variables (phosphorus and nitrogen derivatives).
#' @param log_transform whether to compute the fences on the log10 scale.
#'
#' @return tibble with one row per variable.
#' @export
summarise_env_outliers <- function(
  data,
  vars = c(
    "total_phosphorus", "phosphate", "ammonium",
    "nitrite", "nitrate", "kjeldahl_nitrogen"
  ),
  log_transform = TRUE
) {
  purrr::map(vars, \(v) {
    x <- data[[v]]
    x <- x[!is.na(x)]
    x_pos <- if (log_transform) x[x > 0] else x
    x_t <- if (log_transform) log10(x_pos) else x_pos

    q <- stats::quantile(x_t, c(0.25, 0.75))
    iqr <- q[[2]] - q[[1]]
    lower <- q[[1]] - 1.5 * iqr
    upper <- q[[2]] + 1.5 * iqr
    is_outlier <- x_t < lower | x_t > upper

    tibble::tibble(
      variable = v,
      n = length(x),
      n_outlier = sum(is_outlier),
      pct_outlier = round(100 * sum(is_outlier) / length(x_t), 2),
      fence_low = if (log_transform) 10^lower else lower,
      fence_high = if (log_transform) 10^upper else upper,
      min = min(x),
      max = max(x)
    )
  }) |>
    purrr::list_rbind()
}

#' Boxplot of environmental variable distributions
#'
#' Boxplots on a log10 scale to visually spot extreme values: points beyond
#' the whiskers are flagged as outliers using Tukey's fences (1.5 * IQR),
#' the same convention as [summarise_env_outliers()] and
#' [ggplot2::geom_boxplot()]. Non-positive values are dropped since they
#' cannot be shown on a log scale. Optionally overlays a per-variable
#' reference threshold, e.g. a regulatory quality-class boundary.
#'
#' @param data tibble containing the environmental variables as columns.
#' @param vars character vector of variable names to plot. Defaults to the
#' nutrient variables (phosphorus and nitrogen derivatives).
#' @param thresholds optional tibble with a `variable` column (matching
#' `vars`) and a `threshold_col` column giving the reference value to
#' overlay for each variable, e.g. the SEQ-Eau grid (see
#' `data/river/seq_eau_nutrient_thresholds.csv`).
#' @param threshold_col name of the column in `thresholds` holding the
#' reference value.
#' @param threshold_label legend label for the threshold reference line.
#'
#' @return ggplot
#' @export
plot_env_boxplot <- function(
  data,
  vars = c(
    "total_phosphorus", "phosphate", "ammonium",
    "nitrite", "nitrate", "kjeldahl_nitrogen"
  ),
  thresholds = NULL,
  threshold_col = "mauvais_min",
  threshold_label = "SEQ-Eau \"mauvais\" threshold"
) {
  p <- data |>
    sf::st_drop_geometry() |>
    dplyr::select(dplyr::all_of(vars)) |>
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to = "variable", values_to = "value"
    ) |>
    dplyr::filter(!is.na(value), value > 0) |>
    ggplot2::ggplot(ggplot2::aes(x = variable, y = value)) +
    ggplot2::geom_boxplot(
      fill = "#cde2fb", color = "#104281",
      outlier.color = "#db8474", outlier.alpha = 0.5
    )

  if (!is.null(thresholds)) {
    thresholds_plot <- thresholds |>
      dplyr::filter(variable %in% vars) |>
      dplyr::rename(
        threshold = dplyr::all_of(threshold_col)
      )

    p <- p +
      ggplot2::geom_errorbar(
        data = thresholds_plot,
        ggplot2::aes(
          x = variable, ymin = threshold, ymax = threshold,
          color = threshold_label
        ),
        inherit.aes = FALSE,
        linewidth = 0.8, width = 0.7
      ) +
      ggplot2::scale_color_manual(
        name = NULL,
        values = stats::setNames("#9f66b3", threshold_label)
      )
  }

  p +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = NULL, y = "Concentration [log10 scale]") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}
