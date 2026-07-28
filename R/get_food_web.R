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
#' @param local_id column name in `size_table` identifying the sampling
#' operation (fishing operation / tow), used as the grouping unit for local
#' food webs. Ignored if `size_table` already has an `operation_id` column.
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
  if (!"operation_id" %in% names(size_table)) {
    size_table <- size_table |> rename(operation_id = all_of(local_id))
  }
  if (!"batch_id" %in% names(size_table)) {
    # foodwebbuilder requires a batch_id column but never uses its value
    # (only carries it through) - synthesize a unique one for data with no
    # natural batching concept, e.g. sea surveys where every fish is
    # individually measured.
    size_table <- size_table |> mutate(batch_id = dplyr::row_number())
  }
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
    res$local <- foodwebbuilder::build_local_foodweb(
      ind_measure = size_clean,
      local_id = "operation_id",
      metaweb = metaweb,
      tab_size_classes = size_classes,
      selected_resources = selected_resources
    )
  }
  res
}

get_metaweb_river <- function(size_table, diet_fish, predation_window) {}

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

#' Compute local food web metrics and attach operation metadata.
#'
#' @param web_list output of `get_metaweb()` with `local = TRUE`.
#' @param operation tibble with an `operation_id` column identifying each
#' sampling operation, plus whatever metadata to attach (e.g. year,
#' longitude, latitude).
#' @param resource resource species list.
#'
#' @return tibble of local food webs with metrics and operation metadata.
#' @export
prepare_local_foodwebs <- function(web_list, operation, resource) {
  foodweb <- tibble::enframe(
    web_list$local, name = "operation_id", value = "foodweb"
  ) |>
    mutate(
      log_trophic_richness = log(purrr::map_dbl(foodweb, get_trophic_richness)),
      log_species_richness = log(purrr::map_dbl(foodweb, get_species_richness)),
      connectance = purrr::map_dbl(foodweb, get_connectance),
      trophic_length = purrr::map_dbl(foodweb, get_trophic_length),
      trophic_breadth_q90 = purrr::map_dbl(
        foodweb, \(x) get_trophic_breadth(x)$q90
      ),
      trophic_breadth_median = purrr::map_dbl(
        foodweb, \(x) get_trophic_breadth(x)$median
      ),
      frac_piscivorous = purrr::map_dbl(
        foodweb,
        \(x) get_frac_piscivorous(x, resource = resource)
      )
    )
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
    select(-c(measured, year)) |>
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

#' Get metadata for each ASPE fishing operation.
#'
#' @param individual_fish_file path to `output_individual_fish.rda`.
#'
#' @return tibble with columns operation_id, year, longitude, latitude.
#' @export
get_river_operation <- function(individual_fish_file) {
  out <- get(base::load(individual_fish_file))

  out$fishing_operation |>
    dplyr::left_join(out$station, by = dplyr::join_by(site_id)) |>
    dplyr::transmute(
      operation_id = as.character(operation_id),
      year = lubridate::year(date),
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
#' `get_metaweb()`.
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
