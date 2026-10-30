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

#' Select the large dams of the ROE.
#'
#' Large dams are the dams and weirs (ROE types "barrage" and "seuil en
#' rivière") with a drop height of 5 m or more, existing or destroyed (the
#' destroyed large dams with sampled sites upstream were removed at the end
#' of the sampling or after). Other types (culverts, bridge sills, levees,
#' groynes, fish-farm grilles) are left out. The drop height class is
#' unknown for many obstacles, so the drop height is taken from (see D12):
#'
#' 1. the ROE drop height class (5-10 m or >= 10 m);
#' 2. when the class is unknown, the structure height of dams (not weirs):
#'    >= 7.5 m for a drop of 5 m or more and >= 15 m for 10 m or more. Among
#'    the ROE dams with both heights known, nine in ten or more above these
#'    thresholds have such a drop.
#'
#' The measured drop height and the structure height of weirs are not used:
#' they hold impossible values (weirs of tens to hundreds of metres, often
#' their length).
#'
#' @param roe sf tibble, e.g. `read_roe()`'s output.
#' @return sf tibble, one row per large dam, with `kind` ("dam" or "weir"),
#'   `status`, the ROE `river` name, `drop_height_min` (5 or 10 m), the
#'   measured `drop_height` (m), the `structure_height` (m) and the
#'   `evidence` of its drop height ("class" or "structure height").
#' @export
select_large_dams <- function(roe) {
  obstacle <- roe |>
    dplyr::filter(
      LbEtOuvrage %in% c("Existant", "Détruit entièrement", "Détruit partiellement"),
      grepl("barrage|seuil", LbTypeOuvrage, ignore.case = TRUE)
    ) |>
    dplyr::transmute(
      roe_code = CdObstEcoul,
      name = NomPrincipalObstEcoul,
      type = LbTypeOuvrage,
      kind = dplyr::if_else(grepl("seuil", type, ignore.case = TRUE), "weir", "dam"),
      status = LbEtOuvrage,
      river = NomEntiteHydrographique,
      use = LbUsageObstEcoul1,
      altitude = AltiPointCarOuvrage,
      class = CdHautChutClObstEcoul,
      drop_height = suppressWarnings(as.numeric(HautChutEtObstEcoul)),
      structure_height = suppressWarnings(as.numeric(HautMaxTer))
    )
  unknown <- !obstacle$class %in% as.character(1:8) # Class "0": unknown.
  obstacle |>
    dplyr::mutate(
      evidence = dplyr::case_when(
        class %in% c("7", "8") ~ "class",
        unknown & kind == "dam" & dplyr::coalesce(structure_height >= 7.5, FALSE) ~ "structure height"
      ),
      drop_height_min = dplyr::case_when(
        class == "8" ~ 10,
        evidence == "structure height" & structure_height >= 15 ~ 10,
        TRUE ~ 5
      )
    ) |>
    dplyr::filter(!is.na(evidence)) |>
    dplyr::select(-class)
}

#' Simplified river name, for comparing names across sources.
#'
#' @param name character vector of river names.
#' @return lower-case names without accents, generic words ("rivière",
#'   "ruisseau", articles, ...) and punctuation.
#' @export
simplify_river_name <- function(name) {
  name <- tolower(iconv(name, to = "ASCII//TRANSLIT"))
  name <- gsub("['-]", " ", name)
  name <- gsub(
    "\\b(riviere|fleuve|ruisseau|torrent|gave|la|le|les|l|d|de|du|des|et|a)\\b",
    " ", name
  )
  trimws(gsub("\\s+", " ", name))
}

#' Place the large dams on the river network.
#'
#' Each dam is placed on the nearest network edge within `max_dist_m`,
#' preferring an edge with the dam's ROE river name: at a confluence or
#' beside a canal, the nearest edge can be the wrong river. The edge is split at
#' the projection of the dam, with a new node, and its length is shared
#' between the two new edges in proportion. Several dams on the same edge
#' split it in turn.
#'
#' @param network list(nodes, edges) with edge `length`, e.g.
#'   `restrict_river_network()`'s output.
#' @param network_geo the same network with edge geometry (from `from` to
#'   `to`) and `name`, e.g. `add_edge_geometry()`'s output.
#' @param dam sf tibble with `river`, e.g. `select_large_dams()`'s output.
#' @param max_dist_m maximum distance between a dam and the edge (m).
#' @return list(network, dam): `network` with the new nodes and split edges;
#'   `dam` with the `node_id` of the placed dams (`NA` if left out),
#'   `match` ("nearest", "name" when an edge with its river name was
#'   preferred to the nearest, `NA` if left out) and `snap_dist_m`, the
#'   distance to the edge.
#' @export
snap_dams_network <- function(network, network_geo, dam, max_dist_m) {
  edges <- sf::st_as_sf(network_geo$edges)
  dam <- sf::st_transform(dam, sf::st_crs(edges))
  nearest <- sf::st_nearest_feature(dam, edges)
  # Nearest edge with the dam's river name, within `max_dist_m`.
  river <- simplify_river_name(dam$river)
  # Buffer intersection: indexed, unlike st_is_within_distance() here.
  candidate <- sf::st_intersects(sf::st_buffer(dam, max_dist_m), edges)
  named <- mapply(function(edge, r, i) {
    if (is.na(r) || r == "") return(NA_integer_)
    name <- simplify_river_name(edges$name[edge])
    same <- !is.na(name) & name != "" &
      (vapply(name, grepl, logical(1), x = r, fixed = TRUE) | grepl(r, name, fixed = TRUE))
    if (!any(same)) return(NA_integer_)
    edge <- edge[same]
    edge[which.min(as.numeric(sf::st_distance(dam[i, ], edges[edge, ])))]
  }, candidate, river, seq_len(nrow(dam)))
  chosen <- dplyr::coalesce(named, nearest)
  dam$snap_dist_m <- as.numeric(sf::st_distance(dam, edges[chosen, ], by_element = TRUE))
  placed <- dam$snap_dist_m <= max_dist_m
  dam$match <- dplyr::case_when(
    !placed ~ NA_character_,
    !is.na(named) & named != nearest ~ "name",
    TRUE ~ "nearest"
  )
  split <- tibble::tibble(
    dam_row = which(placed),
    from = edges$from[chosen[placed]],
    to = edges$to[chosen[placed]],
    fraction = sf::st_line_project(
      sf::st_geometry(edges)[chosen[placed]], sf::st_geometry(dam)[placed],
      normalized = TRUE
    )
  ) |>
    dplyr::arrange(from, to, fraction) |>
    dplyr::mutate(node_id = max(network$nodes$node_id) + dplyr::row_number())
  # One edge per pair of nodes: parallel duplicates are replaced by the split edge.
  split_edge <- network$edges |>
    dplyr::semi_join(split, by = c("from", "to")) |>
    dplyr::distinct(from, to, .keep_all = TRUE)
  new_edge <- split |>
    dplyr::inner_join(split_edge, by = c("from", "to")) |>
    dplyr::reframe(
      new_from = c(from[1], node_id),
      new_to = c(node_id, to[1]),
      # At least 1 cm: a zero-length edge would make two nodes each
      # other's parent towards the mouth.
      length = pmax(diff(c(0, fraction, 1)) * length[1], 0.01),
      name = name[1],
      artificial = artificial[1],
      .by = c(from, to)
    ) |>
    dplyr::select(from = new_from, to = new_to, length, name, artificial)
  new_node <- tibble::tibble(
    node_id = split$node_id,
    geometry = sf::st_geometry(dam)[split$dam_row],
    is_river_mouth = FALSE
  )
  dam$node_id <- NA_integer_
  dam$node_id[split$dam_row] <- split$node_id
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

#' Count the dams close to each other along a river as one.
#'
#' Two dams less than `max_dist_m` apart on the same flow path act as one
#' barrier for fish moving along the river (and are often two ROE codes of
#' one structure: a weir and its power plant, a dam and its spillway).
#' Going up from the sea, a dam is merged into the nearest kept dam
#' downstream if it lies within `max_dist_m` of it, so that a series of
#' dams is not merged into one beyond that distance.
#'
#' @param network list(nodes, edges), e.g. `snap_dams_network()`'s output.
#' @param dam sf tibble with `roe_code`, `node_id`, `match`, e.g.
#'   `snap_dams_network()`'s output.
#' @param river_mouth tibble with `node_id`, `matched`, e.g.
#'   `match_mouth_district()`'s output.
#' @param max_dist_m maximum river distance between two dams counted as one
#'   (m).
#' @return `dam` with `merged_into` (ROE code of the kept dam, `NA` if
#'   kept); the merged dams get `match` "merged" and no `node_id`, so that
#'   they are not counted.
#' @export
merge_close_dams <- function(network, dam, river_mouth, max_dist_m) {
  placed <- dam |>
    sf::st_drop_geometry() |>
    dplyr::filter(!is.na(node_id)) |>
    dplyr::select(roe_code, node_id)
  close <- locate_operations_dams(
    network, dplyr::rename(placed, operation_id = roe_code), dam, river_mouth
  ) |>
    dplyr::filter(position == "upstream", river_dist_m <= max_dist_m, operation_id != roe_code)
  route <- route_to_mouth(network, river_mouth$node_id[river_mouth$matched])
  placed <- placed[order(route$dist[placed$node_id]), ]
  kept <- character(0)
  dam$merged_into <- NA_character_
  for (code in placed$roe_code) {
    below <- close[close$operation_id == code & close$roe_code %in% kept, ]
    if (nrow(below) > 0) {
      dam$merged_into[dam$roe_code == code] <- below$roe_code[which.min(below$river_dist_m)]
    } else {
      kept <- c(kept, code)
    }
  }
  merged <- !is.na(dam$merged_into)
  dam$match[merged] <- "merged"
  dam$node_id[merged] <- NA_integer_
  dam
}

#' Locate operations upstream and downstream of dams along the network.
#'
#' The river network is oriented by the shortest-path tree from the river
#' mouths: each node flows to its parent towards the nearest mouth, avoiding
#' artificial channels (see `route_to_mouth()`). An
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
  route <- route_to_mouth(network, river_mouth$node_id[river_mouth$matched])
  parent <- route$parent
  dist <- route$dist
  source <- route$source

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
