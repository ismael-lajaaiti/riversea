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
      DietCode,
      FoodI,
      FoodII,
      FoodIII,
      ItemName,
      Stage,
      DietPercent
    ) |>
    filter(DietPercent > 10) |>
    rename(PredatorStage = SampleStage, PreyStage = Stage) |>
    mutate(
      prey_category = case_when(
        FoodII == "phytoplankton" ~ "phytoplankton",
        FoodII == "other plants" &
          tolower(ItemName) == "macrophyte" ~ "macrophyte",
        FoodII == "other plants" &
          FoodIII == "terrestrial plants" ~ "macrophyte",
        FoodII == "other plants" &
          grepl("macrophyte", tolower(ItemName)) ~ "macrophyte",
        FoodII == "other plants" & grepl("algae", tolower(FoodIII)) ~ "biofilm",
        FoodI == "nekton" ~ "fish",
        TRUE ~ FoodI
      )
    ) |>
    select(Species, PredatorStage, prey_category) |>
    distinct(Species, PredatorStage, prey_category) |>
    arrange(Species, desc(PredatorStage)) |>
    rename(species = Species, stage = PredatorStage)
  diet_category
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
    select(-others) |>
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
