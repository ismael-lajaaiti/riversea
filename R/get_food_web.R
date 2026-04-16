#' Build the metaweb representing all possible interactions.
#'
#' @param size_table data frame of each individual fish size and the operation
#' id.
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
#' @param local_id column name indicating the location.
#'
#' @return metaweb as a matrix.
#' @export
get_metaweb <- function(
  size_table,
  diet_fish,
  diet_resource,
  predation_window,
  local = FALSE,
  num_classes = 5,
  local_id = "trait",
  selected_resources = c(
    "zooplankton",
    "phytoplankton",
    "biofilm",
    "zoobenthos",
    "macrophyte",
    "detritus"
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
  predation_window <- predation_window |> rename(species_code = species)
  diet_resource <- diet_resource |>
    mutate(species = if_else(
      species == "biofilms",
      "biofilm",
      species
    ))
  diet_resource$species_code <- diet_resource$species
  diet_resource <- diet_resource |> rename(
    detritus = det,
    biofilm = biof,
    macrophyte = macroph,
    phytoplankton = phytopl,
    zooplankton = zoopl,
    zoobenthos = zoob
  )
  # Clean data.
  size_clean <- remove_missing_species(
    ind_measure = size_table,
    fish_diet_shift = diet_fish,
    pred_win = predation_window
  )
  size_classes <- compute_size_classes(
    ind_measure = size_clean,
    num_classes = num_classes
  )
  metaweb <- build_metaweb(
    tab_size_classes = size_classes,
    pred_win = predation_window,
    fish_diet_shift = diet_fish,
    resource_diet_shift = diet_resource,
    num_classes = num_classes,
    selected_resources = selected_resources
  )
  res <- list(metaweb = metaweb, size_class = size_classes)
  if (local) {
    res$local <- build_local_foodweb(
      ind_measure = size_clean,
      local_id = local_id,
      metaweb = metaweb,
      tab_size_classes = size_classes,
      selected_resources = selected_resources
    )
  }
  res
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

plot_sizeclass_connectance <- function(metaweb_table) {
  metaweb_table |>
    mutate(connectance = purrr::map_dbl(metaweb, get_connectance)) |>
    ggplot(aes(num_classes, connectance)) +
    geom_point() +
    geom_line() +
    labs(x = "Number of size classes", y = "Metaweb connectance")
}
