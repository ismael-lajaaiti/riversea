#' Get fish diet with standardized categories.
#'
#' @param species A character vector of species names.
#' @return A data frame with species name ("species"), predator stage ("stage")
#' and prey category ("prey_category").
#' @import dplyr
#' @export
get_diet_category <- function(species) {
  items <- rfishbase::fooditems(species)
  diet_table <- items |>
    select(
      Species,
      PredatorStage,
      FoodI,
      FoodII,
      FoodIII,
      CommonessII
    ) |>
    filter_out(grepl("rare", CommonessII)) |>
    rename(species = Species, stage = PredatorStage) |>
    mutate_prey_category() |>
    select(species, stage, prey_category) |>
    distinct() |>
    arrange(species, desc(stage))
  missing <- tibble(species = species) |>
    filter_out(species %in% diet_table$species)
  diet_table |>
    bind_rows(missing) |>
    arrange(species, stage)
}

get_missing_diet <- function(species_init, diet, sea_data) {
  species_with_info <- diet |>
    filter_out(is.na(prey_category)) |>
    pull(species) |>
    unique()
  species_no_info <- sort(setdiff(species_init, species_with_info))
  sea_data$catch |>
    filter(species %in% species_no_info) |>
    group_by(species) |>
    summarise(occurence = n(), .groups = "drop")
}

mutate_prey_category <- function(df) {
  df |>
    mutate(
      prey_category = case_when(
        FoodII == "phytoplankton" ~ "phytoplankton",
        FoodII == "other plants" & grepl("algae", tolower(FoodIII)) ~ "biofilm",
        FoodII == "other plants" & FoodIII == "periphyton" ~ "biofilm",
        FoodII == "other plants" ~ "macrophyte",
        FoodI == "nekton" ~ "fish",
        FoodI == "zoobenthos" & grepl("crust", FoodII) ~ "crustacean",
        FoodI == "zoobenthos" & FoodII == "insects" ~ "insect",
        FoodI == "zoobenthos" & FoodII == "mollusks" ~ "mollusk",
        FoodI == "zoobenthos" & FoodII == "echinoderms" ~ "echinoderm",
        FoodI == "zoobenthos" & FoodII == "worms" ~ "worm",
        FoodI == "zoobenthos" ~ NA_character_, # Ignore rare zoobenthos groups.
        TRUE ~ FoodI
      )
    )
}

#' Widen the dataframe of species diets.
#'
#' @param diet_category Output of get_diet_category.
#' @return Same as get_diet_categery, but in wide format.
#' @import dplyr
#' @importFrom tidyr pivot_wider
#' @export
widen_diet_category <- function(diet_category) {
  diet_category_wide <- diet_category |>
    mutate(eat = 1) |>
    pivot_wider(
      names_from = prey_category,
      values_from = eat,
      values_fill = list(eat = 0)
    ) |>
    select(-c(others, "NA"))
  cols <- names(diet_category_wide)
  diet_category_wide |>
    select(species, stage, sort(setdiff(cols, c("species", "stage")))) |>
    arrange(species, desc(stage))
}

#' Get fish diet with standardized categories in a wide format.
#'
#' @param species List of species names.
#' @return Same as get_diet_categery, but in wide format.
#' @export
get_diet_category_wide <- function(species) {
  species |>
    get_diet_category() |>
    widen_diet_category()
}

#' Get the species list of the data set.
#'
#' Assume that species names are given by 'species_valid' column.
#'
#' @param data_list
#'
#' @return species list
#' @export
get_species_list <- function(data_list) {
  data_list$size |>
    pull(species_valid) |>
    unique()
}

#' Add missing larvae stages to the diet table.
#'
#' @param diet_table wide diet table.
#'
#' @return data.frame diet table with larvae stage.
#' @export
add_larvae <- function(diet_table) {
  species_no_larvae <- diet_table |>
    group_by(species) |>
    summarise(has_larvae = "larvae" %in% stage) |>
    filter(!has_larvae) |>
    pull(species)
  larvae_rows <- diet_table |>
    filter(
      species %in% species_no_larvae,
    ) |>
    mutate(
      stage = "larvae",
      across(-c(species, stage), ~0)
    ) |>
    distinct() |>
    mutate(zooplankton = 1)
  diet_table |>
    rbind(larvae_rows) |>
    arrange(species, desc(stage))
}

#' Filter rare species from diet table.
#'
#' Rare species are defined on occurences.
#' A species is rare if its occurence is lower than `occurence_min`.
#'
#' @param diet_table diet table wide
#' @param sea_data data list of sea surveys
#' @param occurence_min set to 10 by default, define when a species is rare
#'
#' @return data.frame diet table without rare species
#' @export
filter_out_rare <- function(diet_table, sea_data, occurence_min = 10) {
  species_to_keep <- sea_data$catch |>
    group_by(species_valid) |>
    summarise(occurence = n(), .groups = "drop") |>
    filter(occurence >= occurence_min) |>
    pull(species_valid)
  diet <- diet_table |> filter(species %in% species_to_keep)
  n_removed <- length(setdiff(unique(diet_table$species), species_to_keep))
  list(diet = diet, n_removed = n_removed)
}

#' Get the list of piscivorous species from the diet table.
#'
#' @param diet_table
#'
#' @return vector of piscivorous species
get_piscivorous <- function(diet_table) {
  diet_table |>
    filter(fish == 1) |>
    pull(species) |>
    unique()
}

#' Create predation window table.
#'
#' Beta_min and beta_max corresponds respectively to the minimal and maximal
#' size on which species can feed on, relative to its own body size.
#' For now, and for simplicity, we take constant values for piscivorous species.
#'
#' @param diet_table
#'
#' @return data frame
#' @export
get_predation_window <- function(diet_table) {
  sp_piscivorous <- get_piscivorous(diet_table)
  sp_all <- diet_table |>
    pull(species) |>
    unique()
  sp_remaining <- setdiff(sp_all, sp_piscivorous)
  bind_rows(
    tibble(
      species = sp_piscivorous,
      beta_min = 0.03,
      beta_max = 0.45
    ),
    tibble(
      species = sp_remaining,
      beta_min = 0,
      beta_max = 0
    )
  )
}

#' Import the literature-validated diet table for estuarine and coastal fish.
#'
#' Replaces the fishbase-scraped diet table. Treats a blank `length_max` as
#' an open-ended upper bound (`Inf`), consistent with the adult stage
#' convention used elsewhere in the diet tables.
#'
#' @param file path to the validated diet csv.
#'
#' @return data.frame with columns species, stage, length_min, length_max
#' and prey category columns.
#' @export
get_diet_validated <- function(file) {
  categories <- c(
    "crustacean", "insect", "worm", "zooplankton", "detritus", "macrophyte",
    "mollusk", "fish", "biofilm", "phytoplankton", "echinoderm"
  )
  read.csv(file, stringsAsFactors = FALSE) |>
    mutate(
      length_max = if_else(is.na(length_max), Inf, length_max),
      across(all_of(categories), as.numeric)
    ) |>
    select(species, stage, length_min, length_max, all_of(categories))
}

#' Merge two diet tables.
#'
#' For a species present in both tables, `diet_new` wins entirely: all of
#' its rows are kept and `diet`'s rows for that species are dropped, rather
#' than mixing rows from both sources for the same species. Small
#' inter-stage boundary buffers are closed (see `close_diet_size_gaps()`)
#' so no species has a size dead zone that would make `build_metaweb()`
#' error out.
#'
#' @param diet diet table, as returned by `get_diet_validated()`.
#' @param diet_new diet table to merge in, taking precedence on any species
#' overlap.
#'
#' @return combined diet table, same schema as `diet`.
#' @export
merge_diet <- function(diet, diet_new) {
  diet |>
    filter(!species %in% diet_new$species) |>
    bind_rows(diet_new) |>
    close_diet_size_gaps()
}
