#' Load raw sea survey data.
#'
#' @param data_folder Path to data folder.
#' @return A nested list with one element per survey.
#' Each survey contains two data frames: `catch` and `size`.
#' @export
load_raw_sea_data <- function(data_folder) {
  out <- list(nurse = list(), pomet = list(), solper = list())

  read_file <- function(file_name) {
    read.csv(file.path(data_folder, file_name), sep = ";")
  }

  out$pomet$catch <- read_file("poisson_pomet.csv")
  out$pomet$size <- read_file("poisson_tailles_pomet.csv")
  out$nurse$catch <- read_file("Captures_Nurse_1980_2023.csv")
  out$nurse$size <- read_file("Tailles_Nurse_1980_2023.csv")
  out$solper$catch <- read_file("captures_SOLPER_2005_2011.csv")
  out$solper$size <- read_file("tailles_SOLPER_2005_2011.csv")

  # TODO: Add reftax.
  # TODO: Add trait.

  out
}

#' Download sea survey data.
#'
#' @param data_folder Where to store data files.
#' @return Nothing.
#' @export
download_sea_data <- function(data_folder) {} # TODO: Write the function.

#' Clean sea data.
#'
#' @param raw_data output of `load_raw_sea_data`.
#'
#' @return List containing cleaned data frames for each sea campagne.
#' @export
clean_sea_data <- function(raw_data) {
  clean_data <- list(
    pomet = list(),
    nurse = list(),
    solper = list()
  )

  # Pomet - Catch.
  clean_data$pomet$catch <- raw_data$pomet$catch |>
    select(c(
      ID_interne_prelevement,
      Annee,
      NomScient,
      Nind_esp
    )) |>
    rename(
      trait = ID_interne_prelevement,
      year = Annee,
      species = NomScient,
      abundance = Nind_esp,
    )

  # Pomet - Size.
  clean_data$pomet$size <- raw_data$pomet$size |>
    select(c(
      ID_interne_prelevement,
      NomScient,
      Longueur_fourche_mm,
      Nind_esp_taille,
      Annee
    )) |>
    rename(
      trait = ID_interne_prelevement,
      year = Annee,
      species = NomScient,
      length = Longueur_fourche_mm,
      batch_size = Nind_esp_taille
    ) |>
    filter(!is.na(length)) |>
    uncount(batch_size) |>
    mutate(length = length |> str_replace(",", ".") |> as.numeric()) |>
    mutate(length = length / 10) # Convert to cm.

  # Nurse - Catch.
  clean_data$nurse$catch <- raw_data$nurse$catch |>
    select(
      -c(
        Campagne,
        Poids
      )
    ) |>
    rename(
      trait = Trait,
      year = Annee,
      species = Espece,
      abundance = Nombre
    )

  # Nurse - Size.
  clean_data$nurse$size <- raw_data$nurse$size |>
    select(
      -c(
        Campagne,
        Sexe,
        Maturite,
        Poids,
        Age
      )
    ) |>
    rename(
      trait = Trait,
      year = Annee,
      species = Espece,
      length = Longueur,
      batch_size = Nombre
    ) |>
    filter(!is.na(length)) |>
    uncount(batch_size) |>
    mutate(length = length |> str_replace(",", ".") |> as.numeric())

  # Solper - Catch.
  clean_data$solper$catch <- raw_data$solper$catch |>
    select(
      -c(
        Campagne,
        Poids
      )
    ) |>
    rename(
      trait = Trait,
      year = Annee,
      species = Espece,
      abundance = Nombre
    )

  # Solper - Size.
  clean_data$solper$size <- raw_data$solper$size |>
    select(
      -c(
        Campagne,
        Poids,
        Sexe,
        Maturite,
        Age
      )
    ) |>
    rename(
      trait = Trait,
      year = Annee,
      species = Espece,
      length = Longueur,
      batch_size = Nombre
    ) |>
    filter(!is.na(length)) |>
    uncount(batch_size) |>
    mutate(length = length |> str_replace(",", ".") |> as.numeric())

  clean_data
}

#' Combine data frames of different surveys.
#'
#' @param clean_data that is the output of `clean_sea_data`.
#' @return combined a list of combined data frame for $catch and $size.
#' @export
download_sea_data <- function(data_folder) {} # TODO: Write the function.

#' Clean sea data.
#'
#' @param raw_data output of `load_raw_sea_data`.
#'
#' @return List containing cleaned data frames for each sea campagne.
#' @export
combine_sea_data <- function(clean_data) {
  combined <- list(catch = c(), size = c())

  combined$catch <- rbind(
    clean_data$pomet$catch |> mutate(survey = "pomet"),
    clean_data$nurse$catch |> mutate(survey = "nurse"),
    clean_data$solper$catch |> mutate(survey = "solper")
  )

  combined$size <- rbind(
    clean_data$pomet$size |> mutate(survey = "pomet"),
    clean_data$nurse$size |> mutate(survey = "nurse"),
    clean_data$solper$size |> mutate(survey = "solper")
  )

  combined
}

#' Entire pipeline to prepare the sea data survey for analysis.
#'
#' Pipe utility functions to process raw sea data.
#'
#' @param dir Directory of raw data.
#' @return list of data.frame $catch and $size.
#' @export
preprocess_sea_data <- function(dir) {
  tidy <- load_raw_sea_data(dir) |>
    clean_sea_data() |>
    combine_sea_data()

  tidy
}
