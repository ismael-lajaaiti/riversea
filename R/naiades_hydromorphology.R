#' Average bankfull width per station, from Naiades hydromorphology data.
#'
#' Averages across all available survey dates per station (`sandre_code`);
#' the two known sentinel/invalid values (0 and 999) are dropped.
#'
#' @param path path to the Naiades hydromorphology `operation.csv` file.
#'
#' @return tibble with `sandre_code` and `river_width` (m).
#' @export
get_river_width <- function(path) {
  readr::read_delim(
    path,
    delim = ";",
    locale = readr::locale(encoding = "latin1", decimal_mark = "."),
    show_col_types = FALSE
  ) |>
    dplyr::filter(
      !is.na(LargeurPleinBord), !(LargeurPleinBord %in% c(0, 999))
    ) |>
    dplyr::group_by(sandre_code = CdStationMesureEauxSurface) |>
    dplyr::summarise(river_width = mean(LargeurPleinBord), .groups = "drop")
}
