#' Remove a fish species from a food web, with secondary extinctions.
#'
#' All size classes of `species` are removed, then any consumer (non-resource
#' node) left without prey is removed too, repeatedly, until every consumer
#' has at least one prey. Resources are never removed.
#'
#' @param web adjacency matrix, rows are prey and columns are predators.
#' @param species species name, as in the node names before `_<size class>`.
#' @param resource resource species list.
#'
#' @return the reduced adjacency matrix, with attribute `n_secondary`: the
#' number of nodes lost to secondary extinctions.
#' @export
remove_species <- function(web, species, resource) {
  keep <- sub("_[0-9]+$", "", colnames(web)) != species
  web <- web[keep, keep, drop = FALSE]
  n_secondary <- 0
  repeat {
    starving <- colSums(web) == 0 & !colnames(web) %in% resource
    if (!any(starving)) break
    n_secondary <- n_secondary + sum(starving)
    web <- web[!starving, !starving, drop = FALSE]
  }
  attr(web, "n_secondary") <- n_secondary
  web
}

#' Prey-averaged trophic level of every node of a food web.
#'
#' Same computation as [get_trophic_length()], which keeps only the maximum.
#'
#' @param web adjacency matrix, rows are prey and columns are predators.
#'
#' @return named numeric vector, one value per node.
#' @export
get_trophic_levels <- function(web) {
  links <- web |>
    as.data.frame() |>
    dplyr::mutate(resource = rownames(web)) |>
    tidyr::pivot_longer(-"resource", names_to = "consumer") |>
    dplyr::filter(value == 1) |>
    dplyr::select(-value)
  community <- cheddar::Community(
    nodes = data.frame(node = colnames(web)),
    properties = list(title = "Community"),
    trophic.links = links
  )
  cheddar::PreyAveragedTrophicLevel(community)
}

#' Food web metrics used to measure species contributions.
#'
#' @param web adjacency matrix, rows are prey and columns are predators.
#' @param resource resource species list.
#'
#' @return named numeric vector.
#' @export
contribution_metrics <- function(web, resource) {
  c(
    trophic_length = get_trophic_length(web),
    connectance = get_connectance(web),
    fish_diet_overlap = get_fish_diet_overlap(web, resource)
  )
}

#' Contribution of each fish species to the structure of many food webs, in
#' parallel.
#'
#' Runs [species_contribution()] on each web, in forked workers, leaving 4
#' cores free so the machine stays responsive.
#'
#' @param foodweb tibble with `operation_id` and the `foodweb` list-column,
#'   e.g. `foodweb_structure` restricted to the river operations.
#' @param resource resource species list.
#'
#' @return tibble, one row per (operation, fish species).
#' @export
measure_species_contribution <- function(foodweb, resource) {
  n_workers <- max(1, parallel::detectCores() - 4)
  n_workers <- min(n_workers, nrow(foodweb))
  batch <- as.integer(cut(seq_len(nrow(foodweb)), n_workers, labels = FALSE))
  results <- parallel::mclapply(
    split(seq_len(nrow(foodweb)), batch),
    \(rows) {
      purrr::map(rows, \(i) {
        species_contribution(foodweb$foodweb[[i]], resource) |>
          dplyr::mutate(operation_id = foodweb$operation_id[i], .before = 1)
      }) |>
        dplyr::bind_rows()
    },
    mc.cores = n_workers
  )
  .stop_on_mclapply_errors(results, "measure_species_contribution")
  dplyr::bind_rows(results)
}

#' Contribution of each fish species to the structure of a food web.
#'
#' Each fish species is removed in turn (see [remove_species()]) and the
#' metrics are recomputed. `delta_*` is the metric of the full web minus the
#' metric without the species (positive: the species raises the metric).
#' `relative_*` is `delta_*` minus its mean over the fish species of the web:
#' every removal lowers richness by one, so this removes the effect of losing
#' any one species and keeps the effect of the species' identity. It is `NA`
#' in webs with a single fish species. `trophic_level` is the species' mean
#' prey-averaged trophic level over its size classes, in the full web.
#'
#' @param web adjacency matrix, rows are prey and columns are predators.
#' @param resource resource species list.
#'
#' @return tibble, one row per fish species of the web.
#' @export
species_contribution <- function(web, resource) {
  node_species <- sub("_[0-9]+$", "", colnames(web))
  species <- setdiff(unique(node_species), resource)
  full <- contribution_metrics(web, resource)
  trophic_level <- get_trophic_levels(web)[colnames(web)]
  rows <- lapply(species, \(sp) {
    reduced <- remove_species(web, sp, resource)
    delta <- full - contribution_metrics(reduced, resource)
    tibble::tibble(
      species = sp,
      trophic_level = mean(trophic_level[node_species == sp]),
      n_secondary = attr(reduced, "n_secondary"),
      delta_trophic_length = delta[["trophic_length"]],
      delta_connectance = delta[["connectance"]],
      delta_fish_diet_overlap = delta[["fish_diet_overlap"]]
    )
  })
  dplyr::bind_rows(rows) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::starts_with("delta_"),
        \(x) if (length(x) > 1) x - mean(x, na.rm = TRUE) else NA_real_,
        .names = "{sub('delta_', 'relative_', .col)}"
      )
    )
}
