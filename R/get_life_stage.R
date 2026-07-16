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

#' Get life stage length for species.
#'
#' To get life stage length we use data provided by Anik to get lg0max, lg1max
#' and length at maturity.
#' When the data is missing, we use fishbase which can provide length at
#' maturity.
#'
#'
#' @param life_stage_file path to the file
#' @param maturity_length data frame with maturity length for each species
#' extracted from fishbase
#'
#' @return dataframe
#' @export
get_life_stage_length <- function(
  life_stage_file,
  maturity_length,
  epsilon = 1e-6
) {
  df <- utils::read.csv2(life_stage_file) |>
    rename(
      species = Espece,
      lg0_max = LG0max,
      maturity_length = L50
    ) |>
    select(species, lg0_max, maturity_length) |>
    mutate(species = rfishbase::validate_names(species))
  df_join <- left_join(maturity_length, df, by = join_by(species)) |>
    mutate(
      juvenile_length_max = coalesce(
        maturity_length.y,
        maturity_length.x,
        lg0_max
      )
    ) |>
    select(-c(maturity_length.x, maturity_length.y, lg0_max)) |>
    transmute(
      species,
      larvae_length_min = 0,
      larvae_length_max = 2 + epsilon,
      juvenile_length_min = 2 + 2 * epsilon,
      juvenile_length_max =
        pmax(juvenile_length_max + epsilon, 2 + 2 * epsilon, na.rm = TRUE),
      adult_length_min = juvenile_length_max + epsilon,
      adult_length_max = Inf
    ) |>
    tidyr::pivot_longer(
      -species,
      names_to = c("stage", ".value"),
      names_sep = "_length_"
    ) |>
    rename(length_max = max, length_min = min)
}

merge_diet_length <- function(diet, length, epsilon = 1e-6) {
  diet_and_length <- diet |>
    filter_out(is.na(stage)) |>
    mutate(
      stage = case_when(
        stage == "recruits/juv." ~ "juvenile",
        stage == "juv./adults" ~ "adult",
        stage == "adults" ~ "adult",
        TRUE ~ stage
      )
    ) |>
    tidyr::pivot_longer(!c(species, stage), names_to = "category", values_to = "eat") |>
    filter_out(eat == 0) |>
    distinct() |>
    pivot_wider(
      names_from = category,
      values_from = eat,
      values_fill = list(eat = 0)
    ) |>
    arrange(species, desc(stage)) |>
    inner_join(length, by = join_by(species, stage)) |>
    select(species, stage, length_min, length_max, everything()) |>
    filter_out_only_larvae()
  without_juvenile <- diet_and_length |>
    group_by(species) |>
    summarise(with_juvenile = "juvenile" %in% stage, .groups = "drop") |>
    filter_out(with_juvenile) |>
    pull(species)
  diet_and_length |> mutate(
    length_min = if_else(
      species %in% without_juvenile & stage == "adult",
      2 + 2 * epsilon,
      length_min
    )
  )
}
