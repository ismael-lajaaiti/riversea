#' Load raw sea survey data.
#'
#' @param data_folder Path to data folder.
#' @return A nested list with one element per survey ("pomet", "nurse", "solper").
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
#' @return
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

    clean_data$pomet$catch <- raw_data$pomet$catch |>
        select(-c(
            Ecologique,
            Trophique,
            Position,
            Commentaire,
            Date,
            Heure,
            Prof,
            Engin_peche,
            Zone
        )) |>
        rename(
            id_passage = ID_interne_passage,
            id_prelevement = ID_interne_prelevement,
            zone = ID_Zone,
            year = Annee,
            month = Mois,
            trait = Trait,
            x_start = Coord_Deb_xmin,
            x_end = Coord_Fin_xmin,
            y_start = Coord_Deb_ymin,
            y_end = Coord_Fin_ymin,
            species = NomScient,
            abundance = Nind_esp,
            weight = Pds_esp
        )

    clean_data$pomet$size <- raw_data$pomet$size |>
        select(c(
            ID_interne_passage,
            ID_interne_prelevement,
            NomScient,
            Longueur_fourche_mm,
            Nind_esp_taille
        )) |>
        rename(
            id_passage = ID_interne_passage,
            id_prelevement = ID_interne_prelevement,
            species = NomScient,
            size = Longueur_fourche_mm,
            batch_size = Nind_esp_taille
        )

    clean_data
}
