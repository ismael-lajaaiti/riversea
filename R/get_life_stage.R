#' Get maturity length using fishbase
#'
#' @param species list of species
#'
#' @return data frame of species and their maturity length.
#' @export
get_maturity_length <- function(species) {
  lm <- rfishbase::maturity(species) |>
    select(Species, Lm) |>
    filter(!is.na(Lm)) |>
    group_by(Species) |>
    summarise(Lm = mean(Lm), .groups = "drop") |>
    rename(species = Species, maturity_length = Lm)
  missing <- data.frame(
    species = setdiff(species, lm$species),
    maturity_length = NA
  )
  rbind(lm, missing)
}

#' Fill diet table with corresponding life stage size
#'
#' @param diet diet table
#' @param maturity_length table of species maturity length
#'
#' @return data.frame
#' @export
merge_diet_size <- function(diet, maturity_length) {
  species_no_lm <- maturity_length |>
    filter(is.na(maturity_length)) |>
    pull(species)
  diet <- diet |>
    filter_out(is.na(stage)) |>
    filter_out(stage == "recruits/juv.") |>
    left_join(maturity_length, by = "species") |>
    mutate(
      length_min = case_when(
        stage == "larvae" ~ 0,
        stage == "juv./adults" ~ 2.01,
        stage == "adults" ~ maturity_length + 0.01
      ),
      length_max = case_when(
        stage == "larvae" ~ 2,
        stage == "juv./adults" ~ maturity_length,
        stage == "adults" ~ Inf
      )
    ) |>
    select(species, length_min, length_max, everything()) |>
    group_by(species) |>
    mutate(
      has_adult = any(stage == "adults"),
      length_max = if_else(stage == "juv./adults" & !has_adult, Inf, length_max)
    ) |>
    ungroup() |>
    select(-has_adult)
  merged_rows <- diet |>
    filter(species %in% species_no_lm) |>
    group_by(species) |>
    summarise(
      length_min = 2.01,
      length_max = Inf,
      stage = "combined", # or "combined" if you prefer
      across(zoobenthos:macrophyte, max),
      .groups = "drop"
    )
  diet |>
    select(-maturity_length) |>
    filter_out(species %in% species_no_lm & stage != "larvae") |>
    rbind(merged_rows) |>
    arrange(species, desc(stage))
}
