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
  out$pomet$trait <- read_file("poisson_pomet.csv")
  out$nurse$catch <- read_file("Captures_Nurse_1980_2023.csv")
  out$nurse$size <- read_file("Tailles_Nurse_1980_2023.csv")
  out$nurse$trait <- read_file("Traits_Nurse_1980_2023.csv")
  out$solper$catch <- read_file("captures_SOLPER_2005_2011.csv")
  out$solper$size <- read_file("tailles_SOLPER_2005_2011.csv")
  out$solper$trait <- NA # For the moment we ignore this survey.

  # TODO: Add reftax.

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
  # Pomet - Trait.
  clean_data$pomet$trait <- raw_data$pomet$trait |>
    select(
      ID_interne_prelevement,
      Annee,
      Mois,
      Coord_Deb_xmin,
      Coord_Deb_ymin,
      Coord_Fin_xmin,
      Coord_Fin_ymin
    ) |>
    rename(
      trait = ID_interne_prelevement,
      year = Annee,
      month = Mois,
      x_start = Coord_Deb_xmin,
      y_start = Coord_Deb_ymin,
      x_end = Coord_Fin_xmin,
      y_end = Coord_Fin_ymin,
    ) |>
    mutate(
      x_start = x_start |> str_replace(",", ".") |> as.numeric(),
      x_end = x_end |> str_replace(",", ".") |> as.numeric(),
      y_start = y_start |> str_replace(",", ".") |> as.numeric(),
      y_end = y_end |> str_replace(",", ".") |> as.numeric(),
      longitude = (x_start + x_end) / 2,
      latitude = (y_start + y_end) / 2
    ) |>
    select(-c(
      x_start,
      y_start,
      x_end,
      y_end
    ))
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
  # Nurse - Trait
  clean_data$nurse$trait <- raw_data$nurse$trait |>
    select(
      Trait,
      Annee,
      Mois,
      Lat,
      Long
    ) |>
    rename(
      trait = Trait,
      year = Annee,
      month = Mois,
      latitude = Lat,
      longitude = Long
    )
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
  # Solper - Trait.
  clean_data$solper$trait <- NA # For now, we exclude SOLPER survey.

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

  combined$trait <- rbind(
    clean_data$pomet$trait |> mutate(survey = "pomet"),
    clean_data$nurse$trait |> mutate(survey = "nurse")
  )

  combined
}

#' Validate species scientific names.
#'
#' For now, we remove all species who have not a name validated by fishbase.
#'
#' @param data_list
#'
#' @return List of cleaned data frame.
#' @export
validate_species_names <- function(data_list, solper_reftax) {
  data <- data_list[c("catch", "size")] |> purrr::map(clean_name, solper_reftax)
  data$trait <- data_list$trait
  data
}

#' Clean species name of a dataframe using fishbase.
#'
#' @param df
#'
#' @return df with cleaned name.
#' @export
clean_name <- function(df, solper_reftax) {
  df |>
    transform_solper_name(solper_reftax) |>
    mutate(
      species_valid = rfishbase::validate_names(species),
      is_valid = !is.na(species_valid)
    )
}

#' Transform solper code to scientific species name
#'
#' @param df data frame with species code names
#' @param fname path to solper reftax file
#' @param solper whether or not to include solper survey data
#'
#' @return data frame with species scientifc names
#' @export
transform_solper_name <- function(df, fname, solper = FALSE) {
  if (solper) {
    solper_reftax <- read.csv2(fname) |>
      select(C_VALIDE, L_VALIDE) |>
      rename(code = C_VALIDE, scientific_name = L_VALIDE)
    df_solper <- df |> filter(survey == "solper")
    df_other <- df |> filter(survey != "solper")
    df_solper <- df_solper |>
      left_join(solper_reftax, by = join_by(species == code)) |>
      select(-species) |>
      rename(species = scientific_name)
    res <- bind_rows(df_other, df_solper)
  } else {
    res <- df |> filter(survey != "solper")
  }
  res
}

#' Keep only species with valid name in data list.
#'
#' @param data_list
#'
#' @return data_list where non-valid species have been removed.
#' @export
filter_valid <- function(data_list) {
  data_list[c("catch", "size")] |> purrr::map(function(df) {
    df |>
      filter(is_valid) |>
      select(-c(is_valid, species))
  })
}

#' Add number of measured individuals to the catch data frame.
#'
#' @param valid_data
#'
#' @return List of data frames.
#' @export
add_observation_number <- function(valid_data) {
  d_measured <- valid_data$size |>
    group_by(trait, year, survey, species_valid) |>
    summarise(n_measured = n(), .groups = "drop")
  valid_data$catch <- valid_data$catch |>
    left_join(d_measured, by = c("trait", "year", "survey", "species_valid")) |>
    mutate(
      abundance_cor = if_else(
        n_measured > abundance | is.na(abundance),
        n_measured,
        abundance
      )
    ) |>
    filter(!is.na(abundance_cor))
  valid_data
}

#' Impute missing fish length.
#'
#' Infer missing missing fish size using a truncated normal distribution.
#' Mean and variance are given by the observed length distribution
#' of the species in the sample.
#' If there three or less observations the variance is set to zero.
#' Extrema of the truncated distribution are given by the species length
#' extrema observed the entire dataset.
#'
#' @param data
#'
#' @return list of data.frame.
#' @export
impute_size <- function(data) {
  extrema <- data$size |>
    filter(!is.na(length)) |>
    group_by(species_valid) |>
    summarise(size_min = min(length), size_max = max(length), .groups = "drop")
  size_params <- data$size |>
    filter(!is.na(length)) |>
    group_by(trait, year, species_valid, survey) |>
    summarise(
      mu = mean(length),
      sigma = sd(length),
      .groups = "drop"
    ) |>
    left_join(data$catch) |>
    filter(n_measured < abundance_cor) |>
    mutate(sigma = if_else(is.na(sigma), 0, sigma)) |>
    select(-c(abundance, n_measured, abundance_cor)) |>
    left_join(extrema, by = "species_valid")
  missing_imputed <- data$catch |>
    mutate(n_missing = abundance_cor - n_measured) |>
    select(-c(abundance, abundance_cor, n_measured)) |>
    filter(n_missing > 0) |>
    left_join(size_params, by = c("trait", "year", "survey", "species_valid")) |>
    rowwise() |>
    mutate(
      length = list(
        truncdist::rtrunc(
          n_missing,
          spec = "norm",
          a = size_min,
          b = size_max,
          mean = mu,
          sd = sigma
        )
      )
    ) |>
    select(-c(n_missing, mu, sigma, size_min, size_max)) |>
    unnest(length) |>
    mutate(measured = FALSE)
  data$size <- data$size |>
    mutate(measured = TRUE)
  data$size <- bind_rows(data$size, missing_imputed)
  data
}

#' Get minimal and maximal size for each fish species
#'
#' @param data_size
#'
#' @return data frame of size extrema
#' @export
get_size_extrema <- function(data_size) {
  data_size |>
    filter(!is.na(length)) |>
    filter(measured) |>
    group_by(species_valid) |>
    summarise(
      length_min = min(length),
      length_max = max(length),
      .groups = "drop"
    )
}

#' Entire pipeline to infer missing sizes from the tidy data list.
#'
#' @param tidy_data
#'
#' @return list of data.frame.
#' @export
infer_missing_size <- function(tidy_data) {
  tidy_data |>
    filter_valid() |>
    add_observation_number() |>
    impute_size()
}

#' Entire pipeline to prepare the sea data survey for analysis.
#'
#' Pipe utility functions to process raw sea data.
#'
#' @param dir Directory of raw data.
#' @return list of data.frame $catch and $size.
#' @export
preprocess_sea_data <- function(dir, solper_reftax) {
  tidy <- load_raw_sea_data(dir) |>
    clean_sea_data() |>
    combine_sea_data() |>
    validate_species_names(solper_reftax)
  tidy$trait <- tidy$trait |> distinct()
  tidy
}
