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
    pivot_longer(-one_of("resource"), names_to = "consumer") |>
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

plot_sizeclass_connectance <- function(metaweb_table) {
  metaweb_table |>
    mutate(connectance = purrr::map_dbl(metaweb, get_connectance)) |>
    ggplot(aes(num_classes, connectance)) +
    geom_point() +
    geom_line() +
    labs(x = "Number of size classes", y = "Metaweb connectance")
}

prepare_local_foodwebs <- function(web_list, sea_data_tidy) {
  foodweb <- enframe(web_list$local, name = "trait", value = "foodweb") |>
    mutate(
      log_trophic_richness = log(purrr::map_dbl(foodweb, get_trophic_richness)),
      log_species_richness = log(purrr::map_dbl(foodweb, get_species_richness)),
      connectance = purrr::map_dbl(foodweb, get_connectance),
      trophic_length = purrr::map_dbl(foodweb, get_trophic_length)
    )
  trait <- sea_data_tidy$trait
  foodweb |> left_join(trait, by = join_by(trait))
}

match_with_environment <- function(foodweb, environment) {
  foodweb |>
    filter(year %in% unique(environment$year)) |>
    group_by(year) |>
    group_modify(~ {
      env_year <- environment |>
        filter(year == .y$year) |>
        st_as_sf(coords = c("x", "y"), crs = 4326)
      fw_year <- .x |>
        st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
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
