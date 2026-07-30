#' Stop with a combined message if any `mclapply()` worker failed.
#'
#' `mclapply()` catches worker errors into `try-error` objects instead of
#' propagating them, so a failed worker would otherwise end up silently
#' mixed into the result list. `conditionMessage()` doesn't work directly
#' on a `try-error` (it isn't itself a condition) - the actual condition is
#' stashed in its "condition" attribute.
#'
#' @param results output of `parallel::mclapply()`.
#' @param caller name of the calling function, used in the error message.
#'
#' @return nothing; called for its side effect of stopping on failure.
#' @noRd
.stop_on_mclapply_errors <- function(results, caller) {
  failed <- vapply(results, function(x) inherits(x, "try-error"), logical(1))
  if (any(failed)) {
    msgs <- vapply(results[failed], function(x) {
      conditionMessage(attr(x, "condition"))
    }, character(1))
    stop(
      caller, "(): ", sum(failed), " worker(s) failed:\n",
      paste(msgs, collapse = "\n"),
      call. = FALSE
    )
  }
}

#' Build the metaweb representing all possible interactions.
#'
#' @param size_table data frame of each individual fish size, with an
#' `operation_id` column identifying the sampling operation (see
#' `merge_size()`).
#' @param diet_fish table of fish diet.
#' @param diet_resource table of resource diet.
#' @param predation_window table of fish predation window, defining the minimal
#' and maximal size at which they can eat (for piscivorous).
#' @param num_classes number of size classes used to build the metaweb, default
#' set to 5.
#' @param selected_resources which resources are considered at the base of the
#' web, by default we consider: zooplankton, phytoplankton, biofilm, zoobenthos,
#' macrophyte and detritus.
#' @param local whether or not to build local food webs.
#'
#' @return metaweb as a matrix.
#' @export
build_foodweb <- function(
  size_table,
  diet_fish,
  diet_resource,
  predation_window,
  local = FALSE,
  num_classes = 5,
  selected_resources = c(
    "biofilm",
    "detritus",
    "phytoplankton",
    "macrophyte",
    "worm",
    "crustacean",
    "mollusk",
    "echinoderm",
    "insect",
    "zooplankton"
  )
) {
  # Format input data.
  size_table <- size_table |>
    rename(species_code = species_valid, size = length)
  diet_fish <- diet_fish |>
    rename(
      species_code = species,
      size_min = length_min,
      size_max = length_max,
    )
  diet_resource <- diet_resource |>
    rename(species_code = species) |>
    filter_out(species_code == "fish")
  predation_window <- predation_window |> rename(species_code = species)
  # Clean data.
  size_clean <- foodwebbuilder::remove_missing_species(
    ind_measure = size_table,
    fish_diet_shift = diet_fish,
    pred_win = predation_window
  )
  size_classes <- foodwebbuilder::compute_size_classes(
    ind_measure = size_clean,
    num_classes = num_classes
  )
  metaweb <- foodwebbuilder::build_metaweb(
    tab_size_classes = size_classes,
    pred_win = predation_window,
    fish_diet_shift = diet_fish,
    resource_diet_shift = diet_resource,
    num_classes = num_classes,
    selected_resources = selected_resources
  )
  res <- list(metaweb = metaweb, size_class = size_classes)
  if (local) {
    res$local <- build_local_foodweb_parallel(
      ind_measure = size_clean,
      local_id = "operation_id",
      metaweb = metaweb,
      tab_size_classes = size_classes,
      selected_resources = selected_resources
    )
  }
  res
}

#' Build local food webs in parallel.
#'
#' `foodwebbuilder::build_local_foodweb()` re-scans the entire `ind_measure`
#' table once per operation, so its cost grows with
#' `n_operations * n_individuals` rather than `n_individuals`. This splits
#' `ind_measure` into `n_workers` operation-disjoint chunks up front (so
#' each worker only ever scans its own small slice) and processes the
#' chunks in parallel forked processes, which share `metaweb` and
#' `tab_size_classes` via copy-on-write rather than each copying them —
#' memory use stays close to that of a single call, dominated by
#' `ind_measure` itself. Forking is Linux/macOS only; `mc.cores` silently
#' falls back to serial on Windows.
#'
#' @param ind_measure output of `foodwebbuilder::remove_missing_species()`.
#' @param local_id column in `ind_measure` identifying the sampling operation.
#' @param metaweb metaweb matrix, as returned by
#' `foodwebbuilder::build_metaweb()`.
#' @param tab_size_classes size class table, as returned by
#' `foodwebbuilder::compute_size_classes()`.
#' @param selected_resources which resources are considered at the base of
#' the web.
#' @param n_workers number of parallel workers, capped at the number of
#' distinct operations.
#'
#' @return named list of local food web matrices, one per operation id.
#' @export
build_local_foodweb_parallel <- function(
  ind_measure,
  local_id,
  metaweb,
  tab_size_classes,
  selected_resources,
  n_workers = parallel::detectCores()
) {
  ids <- unique(ind_measure[[local_id]])
  n_workers <- min(n_workers, length(ids))
  batch_of_id <- stats::setNames(
    as.integer(cut(seq_along(ids), n_workers, labels = FALSE)), ids
  )
  batch <- batch_of_id[as.character(ind_measure[[local_id]])]
  chunks <- split(ind_measure, batch)

  results <- parallel::mclapply(
    chunks,
    function(chunk) {
      foodwebbuilder::build_local_foodweb(
        ind_measure = chunk,
        local_id = local_id,
        metaweb = metaweb,
        tab_size_classes = tab_size_classes,
        selected_resources = selected_resources
      )
    },
    mc.cores = n_workers
  )
  .stop_on_mclapply_errors(results, "build_local_foodweb_parallel")
  do.call(c, unname(results))
}

#' Combine sea and river individual fish size data into one table.
#'
#' Relies on `assert_no_operation_id_collision()` having been run first, so
#' `operation_id` values are carried through unprefixed.
#'
#' @param sea_size sea individual size table, as in `sea_data_imputed$size`.
#' @param river_size river individual size table, as returned by
#' `remove_rare_river_size()`.
#'
#' @return combined tibble with columns operation_id, species_valid, length,
#' survey. Can be passed directly to `build_foodweb()`.
#' @export
merge_size <- function(sea_size, river_size) {
  sea <- sea_size |>
    transmute(
      operation_id = trait,
      species_valid,
      length,
      survey
    )
  river <- river_size |>
    transmute(
      operation_id = as.character(operation_id),
      species_valid,
      length,
      survey = "river"
    )
  bind_rows(sea, river)
}

#' Combine sea and river operation metadata into one table.
#'
#' @param sea_operation sea operation metadata, as in `sea_data_tidy$trait`.
#' @param river_operation river operation metadata, as returned by
#' `get_river_operation()`.
#'
#' @return combined tibble with columns operation_id, year, month, longitude,
#' latitude, survey.
#' @export
merge_operation <- function(sea_operation, river_operation) {
  sea <- sea_operation |>
    transmute(
      operation_id = trait,
      year,
      month,
      longitude,
      latitude,
      survey
    )
  river <- river_operation |>
    transmute(
      operation_id,
      year,
      month,
      longitude,
      latitude,
      survey = "river"
    )
  bind_rows(sea, river)
}

#' Assert that operation ids never collide between sea and river surveys.
#'
#' `foodwebbuilder::build_local_foodweb()` groups individuals into local
#' food webs by `operation_id` alone, with no survey qualifier, so a
#' collision between a sea and a river id would silently pool two
#' unrelated samples into one fabricated local food web instead of
#' erroring. `merge_size()`/`merge_operation()` rely on ids staying unique
#' across surveys without prefixing them; run this first so a future data
#' refresh that introduces a collision fails loudly instead.
#'
#' @param sea_operation sea operation metadata, as in `sea_data_tidy$trait`.
#' @param river_operation river operation metadata, as returned by
#' `get_river_operation()`.
#'
#' @return TRUE, invisibly, if no operation id is shared between surveys.
#' @export
assert_no_operation_id_collision <- function(sea_operation, river_operation) {
  collisions <- intersect(
    as.character(sea_operation$trait),
    as.character(river_operation$operation_id)
  )
  if (length(collisions) > 0) {
    stop(
      "Operation id(s) collide between sea and river surveys: ",
      paste(collisions, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Compute the connectance of the web.
#'
#' @param web
#'
#' @return numeric
#' @export
get_connectance <- function(web) {
  sum(web) / prod(dim(web))
}

#' Compute the trophic richness of food web.
#'
#' Trophic richness count class sizes of the same species as different trophic
#' species, so it is not always equal to the species richness.
#'
#' @param web
#'
#' @return numeric
#' @export
get_trophic_richness <- function(web) {
  dimension <- dim(web)
  if (dimension[1] != dimension[2]) {
    stop("Dimension of the adjacency matrix are different.")
  }
  dimension[1]
}

#' Compute the species richness of the food web
#'
#' @param web matrix
#'
#' @return numeric
#' @export
get_species_richness <- function(web) {
  dimension <- dim(web)
  if (dimension[1] != dimension[2]) {
    stop("Dimension of the adjacency matrix are different.")
  }
  sub("_.*", "", colnames(web)) |>
    unique() |>
    length()
}

#' Fraction of piscivorous fishes.
#'
#' @param web adjacency matix
#' @param resource list of resources so we exclude them when computing the
#' fraction.
#'
#' @return numeric
#' @export
get_frac_piscivorous <- function(web, resource) {
  S_cor <- ncol(web) - length(resource)
  if (S_cor == 0) {
    return(0)
  }
  web |>
    as.data.frame() |>
    tibble::rownames_to_column("prey") |>
    filter_out(prey %in% resource) |>
    summarise(across(where(is.numeric), \(x) if_else(sum(x) >= 1, 1, 0))) |>
    rowwise() |>
    mutate(frac_piscivorous = sum(c_across(where(is.numeric))) / S_cor) |>
    pull(frac_piscivorous)
}

#' Get trophic length from a food web adjacency matrix
#'
#' @param web adjacency matrix of type matrix.
#' @param tl_method string, cheddar method to compute trophic levels.
#'
#' @return numeric, trophic length, that is the maximal trophic level of the
#' community.
#' @export
get_trophic_length <- function(web, tl_method = "ChainAveragedTL") {
  nodes <- data.frame(node = colnames(web))
  links <- web |>
    as.data.frame() |>
    mutate(resource = rownames(web)) |>
    tidyr::pivot_longer(-one_of("resource"), names_to = "consumer") |>
    filter(value == 1) |>
    select(-value)
  # trophic_length <-
  community <- cheddar::Community(
    nodes = nodes,
    properties = list(title = "Community"),
    trophic.links = links
  )
  cheddar::TrophicLevels(community) |>
    as.data.frame() |>
    pull("ChainAveragedTL") |>
    max(na.rm = TRUE)
}

#' Mean diet overlap between predators in a food web.
#'
#' For every pair of predators (nodes with >=1 prey), the Jaccard similarity
#' of their prey sets is computed and averaged across all pairs - low values
#' mean predators partition prey between them (complementarity), high values
#' mean they eat largely the same things (redundancy). Unlike connectance,
#' this is close to independent of species/trophic richness.
#'
#' @param web adjacency matrix, rows are prey and columns are predators.
#'
#' @return numeric, or `NA` if fewer than 2 predators.
#' @export
get_diet_overlap <- function(web) {
  predators <- which(colSums(web) > 0)
  if (length(predators) < 2) {
    return(NA_real_)
  }
  M <- web[, predators, drop = FALSE]
  inter <- t(M) %*% M
  sizes <- colSums(M)
  union_mat <- outer(sizes, sizes, "+") - inter
  jaccard <- inter / union_mat
  diag(jaccard) <- NA
  mean(jaccard[upper.tri(jaccard)], na.rm = TRUE)
}

#' Collapse a trophic-species web to a species-level adjacency matrix.
#'
#' `build_metaweb()`/`build_local_foodweb()` split each species into size
#' classes (e.g. `"Anguilla anguilla_1"`, `"_2"`, ...), which makes a
#' network plot unreadable at the scale of the full metaweb. This unions
#' interactions across size classes, keeping an edge whenever any size
#' class of the prey is eaten by any size class of the predator.
#'
#' @param web adjacency matrix, rows are prey and columns are predators,
#' with trophic-species names of the form `"species_class"`.
#'
#' @return 0/1 adjacency matrix at the species level.
#' @export
aggregate_species_web <- function(web) {
  species <- sub("_.*", "", rownames(web))
  species_levels <- unique(species)
  # M: trophic-species x species membership indicator. t(M) %*% web %*% M
  # sums prey (rows) and predators (columns) into species groups in one step,
  # keeping row = prey / column = predator.
  M <- matrix(
    0, nrow(web), length(species_levels),
    dimnames = list(NULL, species_levels)
  )
  M[cbind(seq_along(species), match(species, species_levels))] <- 1
  agg <- t(M) %*% web %*% M
  (agg > 0) * 1
}

#' Classify species in a metaweb as resource, piscivorous or non-piscivorous.
#'
#' Uses the diet table's per-stage `fish` flag rather than the metaweb's own
#' links, which union all size classes and would mark nearly every species
#' piscivorous.
#'
#' @param node character vector of species/node names, e.g. `colnames(web)`.
#' @param resource character vector of basal resource node names.
#' @param diet table of fish diet, as used by `build_foodweb()`.
#'
#' @return character vector, one of "resource"/"piscivorous"/
#' "non-piscivorous", in the same order as `node`.
#' @export
classify_metaweb_species <- function(node, resource, diet) {
  piscivorous_species <- diet |>
    filter(fish == 1) |>
    pull(species) |>
    unique()
  case_when(
    node %in% resource ~ "resource",
    node %in% piscivorous_species ~ "piscivorous",
    TRUE ~ "non-piscivorous"
  )
}

#' Trophic level via iterative prey-averaging, anchored at basal resources.
#'
#' `TL = 1` (the classical convention: producers = 1), fixed, for resources
#' with no outgoing edges (true autotrophs/inputs - e.g. detritus,
#' phytoplankton). Every other node, including resource categories that
#' themselves have prey (e.g. `worm` eats `mollusk`/`crustacean`, per
#' `diet_resource`), is `1 + mean(TL(prey))`, solved by fixed-point
#' iteration instead of bottom-up recursion because the metaweb has cycles
#' (species preying across size classes) - the same reason
#' `cheddar::TrophicLevels()` returns all-`NA` on the full metaweb.
#'
#' @param g directed igraph, edges prey -> predator (see `agg` in
#' `aggregate_species_web()`).
#' @param resource character vector of basal resource node names.
#' @param max_iter maximum number of fixed-point iterations.
#' @param tol convergence tolerance on the largest per-node change.
#'
#' @return named numeric vector of trophic levels, one per vertex of `g`.
#' @noRd
compute_prey_averaged_trophic_level <- function(
  g, resource, max_iter = 200, tol = 1e-8
) {
  nodes <- igraph::V(g)$name
  adj <- igraph::as_adjacency_matrix(g, sparse = FALSE)
  is_basal <- nodes %in% resource & colSums(adj) == 0
  tl <- stats::setNames(rep(1, length(nodes)), nodes)
  for (i in seq_len(max_iter)) {
    prey_counts <- colSums(adj)
    prey_tl_sum <- as.vector(tl %*% adj)
    new_tl <- ifelse(prey_counts > 0, 1 + prey_tl_sum / prey_counts, 1)
    new_tl[is_basal] <- 1
    names(new_tl) <- nodes
    if (max(abs(new_tl - tl)) < tol) {
      return(new_tl)
    }
    tl <- new_tl
  }
  tl
}

#' Build a trophic-level network layout for a food web.
#'
#' Aggregates `web` to species level and lays it out with `ggraph`'s
#' Sugiyama algorithm: nodes in rows by trophic level (rounded to the
#' nearest integer), spread horizontally within each row to minimise edge
#' crossings; `y` is then set to the actual (continuous) trophic level.
#' Deterministic - no seed needed. Shared by `plot_metaweb()`,
#' `plot_metaweb_interactive()` and `plot_local_foodweb()`.
#'
#' @param web adjacency matrix, as `aggregate_species_web()` expects.
#' @param resource character vector of basal resource node names.
#'
#' @return list(graph = igraph object, layout = `ggraph` layout data frame
#' with `x`/`y`/`name` columns).
#' @noRd
build_trophic_layout <- function(web, resource) {
  agg <- aggregate_species_web(web)
  g <- igraph::graph_from_adjacency_matrix(agg, mode = "directed")
  tl <- compute_prey_averaged_trophic_level(g, resource)[igraph::V(g)$name]
  layer <- round(tl)
  layout <- ggraph::create_layout(g, layout = "sugiyama", layers = layer)
  # Sugiyama spaces nodes within a row to minimise edge crossings, which
  # gives densely-connected rows (e.g. apex predators preying on each
  # other) far more spread than sparse rows (e.g. species that only eat
  # resources) regardless of node count - one row can end up 7x wider than
  # another with more nodes in it. Rank-based x keeps Sugiyama's
  # crossing-minimising order but gives every row the same per-node gap.
  # (Note: dplyr verbs would strip the layout object's class/attributes
  # that ggraph::ggraph() needs, so this reassigns the column in place.)
  layout$x <- ave(
    layout$x, layer[layout$name],
    FUN = function(x) rank(x, ties.method = "first")
  )
  # Only Sugiyama's x is wanted - its own y treats layer 0 as a root at the
  # top (org-chart style), and layers are rounded anyway, so y is set from
  # the actual (continuous) trophic level.
  layout$y <- tl[layout$name]
  list(graph = g, layout = layout)
}

trophic_group_colors <- c(
  resource = "grey60",
  "non-piscivorous" = "#4fae94",
  piscivorous = "#c26b51"
)

# theme_void() leaves the background transparent; force it white so the
# plot doesn't depend on whatever it's composited against.
theme_blank_white <- function() {
  theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

#' Attach a body-length column to a trophic layout, for node sizing.
#'
#' Resource nodes (no length data) are set to half the smallest fish size.
#'
#' @param layout `ggraph` layout data frame with a `name` column.
#' @param size individual fish size table with `species_valid`/`length`
#' columns, e.g. `size_year_filtered`, already restricted to whatever scope
#' the node sizes should reflect (e.g. one operation for a local food web).
#' @param metric "mean" or "max" length per species.
#'
#' @return `layout` with a `body_length` column added.
#' @noRd
attach_body_length <- function(layout, size, metric = c("mean", "max")) {
  metric <- match.arg(metric)
  agg_fn <- if (metric == "max") max else mean
  species_size <- size |>
    group_by(species_valid) |>
    summarise(body_length = agg_fn(length, na.rm = TRUE), .groups = "drop")

  layout |>
    left_join(species_size, by = c(name = "species_valid")) |>
    mutate(
      body_length = if_else(
        group == "resource", min(body_length, na.rm = TRUE) / 2, body_length
      )
    )
}

#' Abbreviate a binomial species name to "G. species" (e.g. for crowded
#' node labels). Names with no genus/species space (resource categories
#' like "detritus") are left untouched.
#'
#' @param x character vector of names.
#'
#' @return character vector.
#' @noRd
abbreviate_species_name <- function(x) {
  sub("^([A-Z])[a-z]+ ", "\\1. ", x)
}

#' Plot a food web as a trophic network, aggregated at the species level.
#'
#' Nodes are laid out with `build_trophic_layout()`'s Sugiyama layers
#' (rows = trophic level, minimising edge crossings within each row),
#' colored by `classify_metaweb_species()` and sized by body length (see
#' `attach_body_length()`).
#'
#' @param web adjacency matrix, as `aggregate_species_web()` expects.
#' @param resource character vector of basal resource node names.
#' @param diet table of fish diet, passed to `classify_metaweb_species()`.
#' @param size individual fish size table with `species_valid`/`length`
#' columns, e.g. `size_year_filtered`.
#' @param metric "mean" or "max" length per species.
#'
#' @return ggplot.
#' @export
plot_metaweb <- function(web, resource, diet, size, metric = c("mean", "max")) {
  metric <- match.arg(metric)
  ld <- build_trophic_layout(web, resource)
  ld$layout$group <- classify_metaweb_species(ld$layout$name, resource, diet)
  ld$layout <- attach_body_length(ld$layout, size, metric)

  ggraph::ggraph(ld$layout) +
    ggraph::geom_edge_link(alpha = 0.05, width = 0.2, color = "grey40") +
    ggraph::geom_node_point(aes(color = group, size = body_length)) +
    scale_color_manual(values = trophic_group_colors) +
    scale_size_continuous(
      range = c(1, 12), name = paste(metric, "length (cm)")
    ) +
    labs(color = NULL) +
    theme_blank_white()
}

#' Plot a single local food web as a trophic network, with species labels.
#'
#' Same layout, coloring and sizing as `plot_metaweb()`, but meant for a
#' single operation's local food web (tens of species rather than ~150),
#' small enough that species names can be shown directly with repelled
#' labels.
#'
#' @inheritParams plot_metaweb
#' @param size individual fish size table, restricted to this operation's
#' own catch (e.g. `size_year_filtered |> filter(operation_id == id)`), so
#' node size reflects the fish actually caught there, not a global average.
#' @param title optional plot title, e.g. the operation's location.
#' @param size_limits fixed `c(min, max)` for the size scale, so several
#' `plot_local_foodweb()` plots share one legend (e.g. via
#' `patchwork::plot_layout(guides = "collect")`) instead of each getting
#' its own. Default `NULL` auto-scales to this plot's own data.
#' @param label_size species label font size.
#'
#' @return ggplot.
#' @export
plot_local_foodweb <- function(
  web, resource, diet, size, title = NULL, metric = c("mean", "max"),
  size_limits = NULL, label_size = 4
) {
  metric <- match.arg(metric)
  ld <- build_trophic_layout(web, resource)
  ld$layout$group <- classify_metaweb_species(ld$layout$name, resource, diet)
  ld$layout <- attach_body_length(ld$layout, size, metric)
  ld$layout$label <- abbreviate_species_name(ld$layout$name)
  # Resources fall back to half this plot's own local min fish length,
  # which can be below a shared `size_limits` floor from another plot's
  # (larger) fish - clamp so they don't silently drop off the scale.
  if (!is.null(size_limits)) {
    ld$layout$body_length <- pmax(ld$layout$body_length, size_limits[1])
  }

  ggraph::ggraph(ld$layout) +
    ggraph::geom_edge_link(alpha = 0.15, width = 0.3, color = "grey40") +
    ggraph::geom_node_point(aes(color = group, size = body_length)) +
    ggraph::geom_node_text(
      aes(label = label),
      repel = TRUE, size = label_size, max.overlaps = 30,
      segment.size = 0.2, segment.alpha = 0.5, color = "grey20", seed = 1
    ) +
    scale_color_manual(values = trophic_group_colors) +
    scale_size_continuous(
      range = c(2, 14), limits = size_limits,
      name = paste(metric, "length (cm)")
    ) +
    labs(color = NULL, title = title) +
    theme_blank_white() +
    theme(plot.title = element_text(face = "bold"))
}

#' Interactive version of `plot_metaweb()`, with species names on hover.
#'
#' Rebuilds the same Sugiyama layout as `plot_metaweb()` but as plain
#' `ggplot2` point/segment geoms (rather than `ggraph`'s specialised network
#' geoms), since those translate reliably through `plotly::ggplotly()`. Text
#' labels are mapped to `plotly`'s hover tooltip instead of being drawn on
#' the plot, so species names stay legible on zoom without cluttering the
#' static view.
#'
#' @inheritParams plot_metaweb
#'
#' @return a `plotly` htmlwidget.
#' @export
plot_metaweb_interactive <- function(
  web, resource, diet, size, metric = c("mean", "max")
) {
  metric <- match.arg(metric)
  ld <- build_trophic_layout(web, resource)
  ld$layout$group <- classify_metaweb_species(ld$layout$name, resource, diet)
  ld$layout <- attach_body_length(ld$layout, size, metric)
  nodes <- ld$layout |> transmute(node = name, x, y, group, body_length)

  edges <- igraph::as_data_frame(ld$graph, what = "edges") |>
    rename(node = from) |>
    left_join(nodes, by = "node") |>
    rename(x_start = x, y_start = y, node_from = node) |>
    rename(node = to) |>
    left_join(nodes, by = "node") |>
    rename(x_end = x, y_end = y)

  p <- ggplot() +
    geom_segment(
      data = edges,
      aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
      color = "grey40", alpha = 0.05, linewidth = 0.2
    ) +
    geom_point(
      data = nodes,
      aes(x = x, y = y, color = group, size = body_length, text = node)
    ) +
    scale_color_manual(values = trophic_group_colors) +
    scale_size_continuous(
      range = c(1, 10), name = paste(metric, "length (cm)")
    ) +
    labs(color = NULL) +
    theme_blank_white()

  plotly::ggplotly(p, tooltip = "text")
}

#' Share of basal-resource feeding links going to each resource category.
#'
#' For one local food web, computes each resource's out-degree (number of
#' predators feeding on it) as a share of all resource-based links, giving
#' the resource composition of the community's diet. `NA` for every
#' resource not present in `web` (rather than 0), so it can't be confused
#' with a resource that was present but unused.
#'
#' @param web adjacency matrix of one local food web.
#' @param resource character vector of basal resource node names.
#'
#' @return named numeric vector, one share per resource in `resource`.
#' @export
summarise_resource_diet <- function(web, resource) {
  resource_rows <- intersect(resource, rownames(web))
  links <- rowSums(web[resource_rows, , drop = FALSE])
  total <- sum(links)
  share <- stats::setNames(rep(NA_real_, length(resource)), resource)
  if (total > 0) {
    share[] <- 0
    share[names(links)] <- links / total
  }
  share
}

plot_sizeclass_connectance <- function(metaweb_table) {
  metaweb_table |>
    mutate(connectance = purrr::map_dbl(metaweb, get_connectance)) |>
    ggplot2::ggplot(ggplot2::aes(num_classes, connectance)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Number of size classes", y = "Metaweb connectance")
}

get_trophic_breadth <- function(web) {
  breadth_values <- apply(web, 2, mean)
  list(
    median = stats::median(breadth_values),
    mean = mean(breadth_values),
    max = max(breadth_values),
    q90 = stats::quantile(breadth_values, 0.9)[[1]]
  )
}

#' Compute structural metrics for one local food web.
#'
#' `get_trophic_length()` (via `cheddar`) dominates the cost of this by
#' roughly two orders of magnitude over the other metrics combined, so
#' bundling every metric into one call keeps each food web's work in a
#' single unit for `measure_foodweb_structure()` to parallelize over.
#'
#' @param foodweb adjacency matrix for one local food web.
#' @param resource resource species list.
#'
#' @return one-row tibble of structural metrics, including the input
#' `foodweb` itself as a list-column.
compute_foodweb_metrics <- function(foodweb, resource) {
  # Metrics are computed into locals before assembling the tibble: naming a
  # tibble() column "foodweb" while a later column expression still refers
  # to the original `foodweb` argument would silently shadow it, since
  # tibble() columns are evaluated sequentially and can see prior columns.
  breadth <- get_trophic_breadth(foodweb)
  log_trophic_richness <- log(get_trophic_richness(foodweb))
  log_species_richness <- log(get_species_richness(foodweb))
  connectance <- get_connectance(foodweb)
  trophic_length <- get_trophic_length(foodweb)
  frac_piscivorous <- get_frac_piscivorous(foodweb, resource = resource)
  diet_overlap <- get_diet_overlap(foodweb)
  tibble::tibble(
    foodweb = list(foodweb),
    log_trophic_richness = log_trophic_richness,
    log_species_richness = log_species_richness,
    connectance = connectance,
    trophic_length = trophic_length,
    trophic_breadth_q90 = breadth$q90,
    trophic_breadth_median = breadth$median,
    frac_piscivorous = frac_piscivorous,
    diet_overlap = diet_overlap
  )
}

#' Compute local food web metrics and attach operation metadata, in parallel.
#'
#' Splits the local food webs into `n_workers` chunks and computes
#' `compute_foodweb_metrics()` for each chunk in a forked worker (see
#' `build_local_foodweb_parallel()` for the memory-sharing rationale).
#' `get_trophic_length()` dominates the per-web cost, so this is worth
#' parallelizing even though each web's computation itself is cheap.
#'
#' @param web_list output of `build_foodweb()` with `local = TRUE`.
#' @param operation tibble with an `operation_id` column identifying each
#' sampling operation, plus whatever metadata to attach (e.g. year,
#' longitude, latitude).
#' @param resource resource species list.
#' @param n_workers number of parallel workers, capped at the number of
#' local food webs.
#'
#' @return tibble of local food webs with metrics and operation metadata.
#' @export
measure_foodweb_structure <- function(
  web_list,
  operation,
  resource,
  n_workers = parallel::detectCores()
) {
  webs <- web_list$local
  n_workers <- min(n_workers, length(webs))
  batch <- as.integer(cut(seq_along(webs), n_workers, labels = FALSE))
  chunks <- split(webs, batch)

  results <- parallel::mclapply(
    chunks,
    function(chunk) {
      purrr::map_dfr(
        chunk, compute_foodweb_metrics,
        resource = resource, .id = "operation_id"
      )
    },
    mc.cores = n_workers
  )
  .stop_on_mclapply_errors(results, "measure_foodweb_structure")
  foodweb <- dplyr::bind_rows(results)
  foodweb |> left_join(operation, by = join_by(operation_id))
}

match_with_environment <- function(foodweb, environment) {
  foodweb |>
    filter(year %in% unique(environment$year)) |>
    group_by(year) |>
    group_modify(~ {
      env_year <- environment |>
        filter(year == .y$year) |>
        sf::st_as_sf(coords = c("x", "y"), crs = 4326)
      fw_year <- .x |>
        sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
      idx <- st_nearest_feature(fw_year, env_year)
      env_matched <- env_year |>
        slice(idx) |>
        select(-year)
      dist_km <- st_distance(fw_year, env_matched, by_element = TRUE) / 1000
      bind_cols(
        fw_year |>
          mutate(
            longitude = st_coordinates(geometry)[, 1],
            latitude = st_coordinates(geometry)[, 2]
          ) |>
          st_drop_geometry(),
        env_matched |> st_drop_geometry()
      ) |>
        mutate(dist_to_env = as.numeric(dist_km))
    }) |>
    ungroup()
}

get_foodweb_size_info <- function(size, diet) {
  size <- size |>
    select(trait, species_valid, length) |>
    rename(species = species_valid)
  diet <- diet |>
    filter(fish == 1) |>
    select(species, length_min, length_max)
  size |>
    inner_join(diet, by = join_by(species), relationship = "many-to-many") |>
    filter(length >= length_min & length <= length_max) |>
    group_by(trait) |>
    summarise(
      max_fish_size = max(length),
      mean_fish_size = mean(length),
    ) |>
    ungroup()
}

#' Get max/mean fish size info for the combined sea/river dataset.
#'
#' Wraps `get_foodweb_size_info()`, which expects a `trait` id column (the
#' sea survey's naming), and reattaches operation metadata so the result can
#' be mapped or compared across surveys directly.
#'
#' @param size combined individual size table, as returned by `merge_size()`
#' and filtered (see `size_year_filtered`).
#' @param diet table of fish diet.
#' @param operation combined operation metadata, as returned by
#' `merge_operation()` (see `operation_year_filtered`).
#'
#' @return tibble of `get_foodweb_size_info()` output joined with `operation`.
#' @export
get_combined_foodweb_size_info <- function(size, diet, operation) {
  get_foodweb_size_info(size |> rename(trait = operation_id), diet) |>
    rename(operation_id = trait) |>
    left_join(operation, by = "operation_id")
}

#' Get metadata for each ASPE fishing operation.
#'
#' `sandre_code` and `date` are kept (rather than only the derived `year`/
#' `month`) because `join_amobio_aspe()` needs them to match river
#' operations against AMOBIO environmental data.
#'
#' @param individual_fish_file path to `output_individual_fish.rda`.
#'
#' @return tibble with columns operation_id, year, month, date, sandre_code,
#' longitude, latitude.
#' @export
get_river_operation <- function(individual_fish_file) {
  out <- get(base::load(individual_fish_file))

  out$fishing_operation |>
    dplyr::left_join(out$station, by = dplyr::join_by(site_id)) |>
    dplyr::transmute(
      operation_id = as.character(operation_id),
      year = lubridate::year(date),
      month = lubridate::month(date),
      date,
      sandre_code,
      longitude = x,
      latitude = y
    )
}

#' Get individual fish sizes from the ASPE river survey data.
#'
#' @param individual_fish_file path to `output_individual_fish.rda`.
#'
#' @return tibble with columns operation_id, batch_id, length (cm),
#' latin_name, species_valid, is_valid. Can be passed directly to
#' `build_foodweb()`.
#' @export
get_river_size <- function(individual_fish_file) {
  out <- get(base::load(individual_fish_file))

  size <- out$fish_individuals |>
    dplyr::mutate(length = size_mm / 10) |>
    dplyr::select(operation_id, batch_id, species_code, length)

  code_to_name <- out$species_ref_aspe |>
    dplyr::select(species_code, latin_name) |>
    dplyr::distinct() |>
    dplyr::mutate(
      # rfishbase can't disambiguate the bare "Salmo trutta" string (4
      # candidate SpecCodes in FishBase's synonym table), but resolves the
      # subspecies form cleanly to the same canonical name.
      latin_name = dplyr::if_else(
        latin_name == "Salmo trutta", "Salmo trutta fario", latin_name
      ),
      species_valid = rfishbase::validate_names(latin_name),
      is_valid = !is.na(species_valid)
    )

  size |>
    dplyr::left_join(code_to_name, by = dplyr::join_by(species_code)) |>
    dplyr::select(-species_code)
}

#' Determine which species clear an occurrence threshold.
#'
#' A species is considered rare if its occurrence — the number of distinct
#' traits (fishing operations / tows) it was recorded in — is equal to or
#' below `occurence_min`. Rows failing name validation (`is_valid == FALSE`)
#' are ignored when computing occurrence.
#'
#' @param data tibble with `species_valid`, `trait` and `is_valid` columns.
#' @param occurence_min numeric threshold defining a rare species.
#'
#' @return character vector of species_valid values that are not rare.
#' @export
get_no_rare_species <- function(data, occurence_min) {
  data |>
    dplyr::filter(is_valid) |>
    dplyr::group_by(species_valid) |>
    dplyr::summarise(n = dplyr::n_distinct(trait), .groups = "drop") |>
    dplyr::filter(n > occurence_min) |>
    dplyr::pull(species_valid)
}

#' Remove rare species from the river size data frame.
#'
#' @param river_size tibble as returned by `get_river_size()`.
#' @param occurence_min numeric threshold defining a rare species.
#'
#' @return river_size, trimmed of rare and unvalidated species.
#' @export
remove_rare_river_size <- function(river_size, occurence_min) {
  no_rare_species <- river_size |>
    dplyr::rename(trait = operation_id) |>
    get_no_rare_species(occurence_min)
  river_size |>
    dplyr::filter(species_valid %in% no_rare_species)
}
