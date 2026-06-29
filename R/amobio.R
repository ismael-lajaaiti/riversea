#' Download the AMOBIO data base
#'
#' From data.gouv at the link
#' https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/EQYVLP
#'
#' @param dir where the data is downloaded.
#'
#' @return nothing
#' @export
download_amobio_data <- function(dir, verbose = FALSE) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  Sys.setenv("DATAVERSE_SERVER" = "entrepot.recherche.data.gouv.fr")
  dataset <- dataverse::get_dataset("doi:10.57745/EQYVLP")
  files <- dataset$files
  paths <- purrr::map2_chr(files$id, files$label, \(id, filename) {
    if (verbose) {
      print(paste("Downloading", filename))
    }
    path <- file.path(dir, filename)
    writeBin(dataverse::get_file(id, progress = verbose), path)
    path
  })
  paths
}

extract_metrics_amobio <- function(paths) {
  metric_file <- paths[grepl("METRICS_AMOBIO.Rdata", paths)]
  df_metric <- get(load(metric_file))
  df_metric |>
    dplyr::select(
      -matches(
        "B_FISH_|B_DIA_|B_INV_|NC_|H_|HM_|T_|P_|IS_|COMPARTMENT|META"
      ),
      -matches("^PC_(?!N_)", perl = TRUE)
    ) |>
    dplyr::rename_with(tolower) |>
    dplyr::rename(dplyr::any_of(rename_amobio())) |>
    dplyr::mutate(
      node_id = as.integer(as.character(node_id)),
      sandre_code = stringr::str_pad(as.character(sandre_code), 8, pad = "0")
    )
}

extract_nodes_amobio <- function(paths) {
  network_file <- paths[grepl("NETWORK_AMOBIO.Rdata", paths)]
  network <- get(load(network_file))
  network |>
    sfnetworks::activate("nodes") |>
    tibble::as_tibble() |>
    tibble::as_tibble() |>
    dplyr::filter(IS_FISH == "FISH") |>
    dplyr::select(node_id, geometry, ELEVATION_IGN25M) |>
    dplyr::rename_with(tolower) |>
    dplyr::rename(dplyr::any_of(rename_amobio()))
}

#' Dictionnary of renamed columns in the AMOBIO database.
#'
#' To keep track of the renaming and ensure consistency.
#'
#' @return named vector
#' @export
rename_amobio <- function() {
  c(
    "date" = "date_operation",
    "sandre_code" = "id",
    "elevation" = "elevation_ign25m"
  )
}
