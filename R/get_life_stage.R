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

filter_out_only_larvae <- function(diet) {
  sp_only_larvae <- diet |>
    group_by(species) |>
    summarise(has_larvae = "larvae" %in% stage, n = n(), .groups = "drop") |>
    filter(has_larvae & n == 1) |>
    pull(species)
  diet |> filter_out(species %in% sp_only_larvae)
}

#' Fill diet table with corresponding life stage size
#'
#' @param diet diet table
#' @param maturity_length table of species maturity length
#'
#' @return data.frame
#' @export
merge_diet_size <- function(diet, maturity_length, epsilon = 1e-6) {
  species_no_lm <- maturity_length |>
    filter(is.na(maturity_length)) |>
    pull(species)
  species_no_juv <- diet |> # Species with no juvenile stage.
    group_by(species) |>
    summarise(has_juv = "juv./adults" %in% stage, .groups = "drop") |>
    pull(species)
  diet <- diet |>
    filter_out(is.na(stage)) |>
    filter_out(stage == "recruits/juv.") |>
    filter_out_only_larvae() |>
    left_join(maturity_length, by = "species") |>
    mutate(
      maturity_length =
        if_else(species %in% species_no_juv, 2 + epsilon, maturity_length)
    ) |>
    mutate(
      length_min = case_when(
        stage == "larvae" ~ 0,
        stage == "juv./adults" ~ 2 + 2 * epsilon,
        stage == "adults" ~ maturity_length + epsilon
      ),
      length_max = case_when(
        stage == "larvae" ~ 2 + epsilon,
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
      length_min = 2 + 2 * epsilon,
      length_max = Inf,
      stage = "combined", # or "combined" if you prefer
      across(biofilm:zooplankton, max),
      .groups = "drop"
    )
  diet |>
    select(-maturity_length) |>
    filter_out(species %in% species_no_lm & stage != "larvae") |>
    rbind(merged_rows) |>
    arrange(species, desc(stage))
}
