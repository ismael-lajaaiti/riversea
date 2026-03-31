#' Get fish diet with standardized categories.
#'
#' @param species A character vector of species names.
#' @return A data frame with species name ("species"), predator stage ("stage") and prey category ("prey_category").
#' @import dplyr
#' @export
get_diet_category <- function(species) {
  items <- rfishbase::fb_tbl("diet_items")
  diet_species <- rfishbase::diet(species)
  diet_category <- diet_species |>
    left_join(items, by = "DietCode") |>
    select(
      Species,
      SampleStage,
      FoodI,
      FoodII,
      FoodIII,
      DietPercent
    ) |>
    filter(DietPercent > 10) |>
    rename(species = Species, stage = SampleStage) |>
    mutate_prey_category() |>
    select(species, stage, prey_category) |>
    distinct() |>
    arrange(species, desc(stage))
  missing_species <- tibble(species = species) |>
    filter(!species %in% diet_category$species) |>
    pull(species)
  diet_fooditems <- get_missing_diet(missing_species)
  still_missing <- tibble(species = missing_species) |>
    filter_out(species %in% diet_fooditems$species)
  diet_category |>
    bind_rows(diet_fooditems) |>
    bind_rows(still_missing) |>
    arrange(species, stage)
}

get_missing_diet <- function(species) {
  items <- rfishbase::fooditems(species)
  diet_category <- items |>
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
  diet_category
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
    select(-c(others, "NA")) |>
    arrange(species, desc(stage))
  diet_category_wide
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
  diet_table |> filter(species %in% species_to_keep)
}
