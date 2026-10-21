#' Download the obstacles to river flow referential (ROE).
#'
#' The ROE (Référentiel des Obstacles à l'Écoulement, OFB) lists the
#' artificial obstacles of French rivers (dams, weirs, sills, ...). The
#' mainland layer (`ObstEcoul_FXX`) is downloaded from the SANDRE WFS
#' service, as a zipped GeoPackage.
#'
#' @param dir directory where to store the zip file.
#' @return Path of `roe_fxx.zip`.
#' @export
download_roe <- function(dir) {
  url <- paste0(
    "https://services.sandre.eaufrance.fr/geo/obs?SERVICE=WFS&VERSION=2.0.0",
    "&REQUEST=GetFeature&TYPENAMES=sa:ObstEcoul_FXX&OUTPUTFORMAT=GPKG"
  )
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  path <- file.path(dir, "roe_fxx.zip")
  curl::curl_download(url, path)
  path
}

#' Read the obstacles to river flow referential (ROE).
#'
#' @param zip_file path of the zipped ROE GeoPackage, e.g. `download_roe()`'s
#'   output.
#' @return sf tibble, one row per obstacle.
#' @export
read_roe <- function(zip_file) {
  dir <- tempfile("roe")
  gpkg <- utils::unzip(zip_file, exdir = dir)
  sf::read_sf(gpkg)
}

#' Extract the ROE obstacles placed on the AMOBIO river network.
#'
#' AMOBIO nodes include the ROE obstacles (`OE_ID`), so obstacles are
#' positioned on the network without snapping.
#'
#' @param paths character vector of AMOBIO file paths.
#' @return tibble with `node_id`, `roe_code`.
#' @export
extract_obstacle_nodes_amobio <- function(paths) {
  network_file <- paths[grepl("NETWORK_AMOBIO.Rdata", paths)]
  network <- get(base::load(network_file))
  network |>
    sfnetworks::activate("nodes") |>
    tibble::as_tibble() |>
    tibble::as_tibble() |>
    dplyr::filter(!is.na(OE_ID)) |>
    dplyr::distinct(roe_code = OE_ID, .keep_all = TRUE) |>
    dplyr::select(node_id, roe_code)
}

#' Select the large dams of the ROE.
#'
#' Large dams are the existing dams and weirs (ROE types "barrage" and
#' "seuil en rivière") whose drop height class is 5-10 m or >= 10 m. Other
#' types (culverts, bridge sills, levees, groynes, fish-farm grilles) are
#' left out. The class is used instead of the measured drop height, which
#' is missing for most obstacles and holds obvious errors (weirs of more
#' than 100 m), while the class agrees with the measured height for
#' well-known dams (see D12).
#'
#' @param roe sf tibble, e.g. `read_roe()`'s output.
#' @param obstacle_node tibble with `node_id`, `roe_code`, e.g.
#'   `extract_obstacle_nodes_amobio()`'s output.
#' @return sf tibble, one row per large dam, with `kind` ("dam" or "weir"),
#'   `drop_height_min` (5 or 10 m), the measured `drop_height` (m) and the
#'   AMOBIO `node_id` (`NA` if the dam is not on the network).
#' @export
select_large_dams <- function(roe, obstacle_node) {
  roe |>
    dplyr::filter(
      LbEtOuvrage == "Existant",
      CdHautChutClObstEcoul %in% c("7", "8"),
      grepl("barrage|seuil", LbTypeOuvrage, ignore.case = TRUE)
    ) |>
    dplyr::transmute(
      roe_code = CdObstEcoul,
      name = NomPrincipalObstEcoul,
      type = LbTypeOuvrage,
      kind = dplyr::if_else(grepl("seuil", type, ignore.case = TRUE), "weir", "dam"),
      river = NomEntiteHydrographique,
      use = LbUsageObstEcoul1,
      altitude = AltiPointCarOuvrage,
      drop_height_min = dplyr::if_else(CdHautChutClObstEcoul == "8", 10, 5),
      drop_height = suppressWarnings(as.numeric(HautChutEtObstEcoul))
    ) |>
    dplyr::left_join(obstacle_node, by = "roe_code")
}

#' Apply the manual review of the large dams.
#'
#' The review table lists, after a check on aerial imagery, the dams left
#' out by the snapping that should be placed on the network (`decision`
#' "snap") and the ROE codes that duplicate another code of the same
#' structure (`drop_roe_code`, e.g. a power plant next to its weir).
#'
#' @param dam sf tibble, e.g. `select_large_dams()`'s output.
#' @param review tibble with `roe_code`, `decision`, `drop_roe_code`, e.g.
#'   read from `data/dam_review.csv`.
#' @return `dam` without the duplicates, with `force_snap` (TRUE for the
#'   dams to place on the network whatever their distance).
#' @export
review_large_dams <- function(dam, review) {
  dam |>
    dplyr::filter(!roe_code %in% review$drop_roe_code) |>
    dplyr::mutate(
      force_snap = roe_code %in% review$roe_code[review$decision == "snap"]
    )
}

#' Snap the large dams that are not on the network.
#'
#' A large dam without a node (its ROE code is on no AMOBIO node) is placed
#' on the nearest network edge if it lies within `max_dist_m` of it, or if
#' its `force_snap` is TRUE (see `review_large_dams()`): the
#' edge is split at the projection of the dam, with a new node, and its
#' length is shared between the two new edges in proportion. Several dams
#' on the same edge split it in turn.
#'
#' @param network list(nodes, edges) with edge `length`, e.g.
#'   `restrict_river_network()`'s output.
#' @param network_geo the same network with edge geometry (from `from` to
#'   `to`), e.g. `add_edge_geometry()`'s output.
#' @param dam sf tibble with `node_id`, `force_snap`, e.g.
#'   `review_large_dams()`'s output.
#' @param max_dist_m maximum distance between a dam and the edge (m).
#' @return list(network, dam): `network` with the new nodes and split edges;
#'   `dam` with `node_id` of the snapped dams filled, `match` ("code",
#'   "snapped" or `NA`, left out) and `snap_dist_m`, the distance to the
#'   nearest edge of the dams without a node.
#' @export
snap_dams_network <- function(network, network_geo, dam, max_dist_m) {
  edges <- sf::st_as_sf(network_geo$edges)
  dam <- sf::st_transform(dam, sf::st_crs(edges))
  off <- which(is.na(dam$node_id))
  nearest <- sf::st_nearest_feature(dam[off, ], edges)
  dam$snap_dist_m <- NA_real_
  dam$snap_dist_m[off] <- as.numeric(
    sf::st_distance(dam[off, ], edges[nearest, ], by_element = TRUE)
  )
  dam$match <- dplyr::if_else(is.na(dam$node_id), NA_character_, "code")
  near <- dam$snap_dist_m[off] <= max_dist_m | dam$force_snap[off]
  split <- tibble::tibble(
    dam_row = off[near],
    from = edges$from[nearest[near]],
    to = edges$to[nearest[near]],
    fraction = sf::st_line_project(
      sf::st_geometry(edges)[nearest[near]], sf::st_geometry(dam)[off[near]],
      normalized = TRUE
    )
  ) |>
    dplyr::arrange(from, to, fraction) |>
    dplyr::mutate(node_id = max(network$nodes$node_id) + dplyr::row_number())
  split_edge <- dplyr::semi_join(network$edges, split, by = c("from", "to"))
  new_edge <- split |>
    dplyr::inner_join(split_edge, by = c("from", "to")) |>
    dplyr::reframe(
      new_from = c(from[1], node_id),
      new_to = c(node_id, to[1]),
      # At least 1 cm: a zero-length edge would make two nodes each
      # other's parent towards the mouth.
      length = pmax(diff(c(0, fraction, 1)) * length[1], 0.01),
      .by = c(from, to)
    ) |>
    dplyr::select(from = new_from, to = new_to, length)
  new_node <- tibble::tibble(
    node_id = split$node_id,
    geometry = sf::st_geometry(dam)[split$dam_row],
    is_river_mouth = FALSE
  )
  dam$node_id[split$dam_row] <- split$node_id
  dam$match[split$dam_row] <- "snapped"
  list(
    network = list(
      nodes = dplyr::bind_rows(network$nodes, new_node),
      edges = dplyr::bind_rows(
        dplyr::anti_join(network$edges, split, by = c("from", "to")), new_edge
      )
    ),
    dam = dam
  )
}

#' Locate operations upstream and downstream of dams along the network.
#'
#' The river network is oriented by the shortest-path tree from the river
#' mouths: each node flows to its parent towards the nearest mouth. An
#' operation is upstream of a dam if the dam lies on its path to the mouth,
#' and downstream if it lies on the dam's path to the mouth. Operations
#' and dams on the same node are left out.
#'
#' @param network list(nodes, edges) with edge `length`, e.g.
#'   `restrict_river_network()`'s output.
#' @param operation tibble with `operation_id`, `node_id`, e.g.
#'   `snap_operations()`'s output.
#' @param dam tibble with `roe_code`, `node_id`, e.g. `select_large_dams()`'s
#'   output.
#' @param river_mouth tibble with `node_id`, `matched`, e.g.
#'   `match_mouth_district()`'s output.
#' @return tibble, one row per (operation, dam) pair on the same flow path,
#'   with `position` of the operation ("upstream" or "downstream" of the
#'   dam) and the `river_dist_m` between them.
#' @export
locate_operations_dams <- function(network, operation, dam, river_mouth) {
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

  # Parent of each node: the neighbour on its shortest path to the mouth.
  edge <- igraph::as_data_frame(g, what = "edges")
  edge <- data.frame(
    node = c(edge$to, edge$from),
    parent = c(edge$from, edge$to),
    weight = c(edge$weight, edge$weight)
  )
  edge <- edge[abs(dist[edge$parent] + edge$weight - dist[edge$node]) < 1e-6, ]
  edge <- edge[!duplicated(edge$node), ]
  parent <- rep(NA_integer_, length(dist))
  parent[edge$node] <- edge$parent
  parent[source] <- NA_integer_

  # Pairs (start, target) where target is on the path from start to the mouth.
  walk <- function(start, target) {
    current <- parent[start]
    i <- seq_along(start)
    pairs <- list()
    while (length(i) > 0) {
      hit <- current %in% target
      pairs[[length(pairs) + 1]] <- data.frame(start = start[i[hit]], target = current[hit])
      keep <- !is.na(current) & current != source
      i <- i[keep]
      current <- parent[current[keep]]
    }
    pairs <- unique(do.call(rbind, pairs))
    pairs$river_dist_m <- dist[pairs$start] - dist[pairs$target]
    pairs
  }

  # Nodes outside the restricted network are absent from the graph.
  operation_node <- intersect(operation$node_id, network$nodes$node_id)
  dam_node <- intersect(dam$node_id, network$nodes$node_id)
  dplyr::bind_rows(
    walk(operation_node, dam_node) |>
      dplyr::transmute(op_node = start, dam_node = target, position = "upstream", river_dist_m),
    walk(dam_node, operation_node) |>
      dplyr::transmute(op_node = target, dam_node = start, position = "downstream", river_dist_m)
  ) |>
    dplyr::inner_join(
      dplyr::select(operation, operation_id, op_node = node_id),
      by = "op_node", relationship = "many-to-many"
    ) |>
    dplyr::inner_join(
      dplyr::select(sf::st_drop_geometry(dam), roe_code, dam_node = node_id),
      by = "dam_node", relationship = "many-to-many"
    ) |>
    dplyr::select(operation_id, roe_code, position, river_dist_m)
}
