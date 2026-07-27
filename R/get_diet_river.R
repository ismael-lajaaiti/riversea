#' Clean and format the validated freshwater diet table.
#'
#' Freshwater diet categories predate the marine/estuarine schema and use a
#' different, coarser resolution: zoobenthos has no direct equivalent and is
#' spread across mollusk, crustacean and insect (worm and echinoderm have no
#' equivalent category and are set to 0); phytobenthos has no equivalent in
#' the marine categories and is dropped. Sizes are recorded in millimeters
#' in the source file and are converted to centimeters, consistent with the
#' marine diet table.
#'
#' @param file path to the validated freshwater diet csv.
#'
#' @return data.frame with the same schema as `get_diet_validated()`:
#' species, stage, length_min, length_max and prey category columns.
#' @export
get_diet_validated_river <- function(file) {
  read.csv(file, stringsAsFactors = FALSE) |>
    transmute(
      species = gsub("_", " ", species_name),
      stage = as.character(stage),
      length_min = size_min / 10,
      length_max = size_max / 10,
      crustacean = zoob,
      insect = zoob,
      worm = 0,
      zooplankton = zoopl,
      detritus = det,
      macrophyte = macroph,
      mollusk = zoob,
      fish,
      biofilm = biof,
      phytoplankton = phytopl,
      echinoderm = 0
    ) |>
    arrange(species, length_min)
}
