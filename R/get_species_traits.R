#' Get species-level ecological traits from FishBase.
#'
#' Retrieves, for each species, its demersal/pelagic position, whether it is
#' diadromous, whether it is native to `country_name`, plus a
#' freshwater/marine/both environment derived from FishBase's Fresh/Saltwater
#' flags.
#'
#' `is_diadromous` is `TRUE` for FishBase `AnaCat` categories that cross the
#' fresh/salt water boundary (`amphidromous`, `anadromous`, `catadromous`),
#' `FALSE` for every other documented category (`non-migratory`,
#' `potamodromous`, `oceanodromous`), and `NA` only when FishBase leaves
#' `AnaCat` genuinely blank rather than coding it "non-migratory" - see the
#' notebook for how that's interpreted downstream.
#'
#' `is_native` is `FALSE` ("introduced") or `TRUE` ("native"): FishBase's
#' "endemic" and "not established" are folded into `TRUE` for simplicity, and
#' the three Ponto-Caspian gobies missing a France record in FishBase's
#' `country()` table (documented recent invaders in the literature) are
#' hardcoded to `FALSE`.
#'
#' @param species character vector of FishBase-validated species names.
#' @param country_name country to look up native status for.
#'
#' @return data frame with columns species, demers_pelag, is_diadromous,
#' environment, is_native.
#' @export
get_species_traits <- function(species, country_name = "France") {
  diadromous_cats <- c("amphidromous", "anadromous", "catadromous")

  habitat <- rfishbase::species(species) |>
    select(
      species = Species,
      fresh = Fresh,
      saltwater = Saltwater,
      demers_pelag = DemersPelag,
      migratory_category = AnaCat
    ) |>
    mutate(
      migratory_category = dplyr::na_if(trimws(migratory_category), ""),
      is_diadromous = case_when(
        is.na(migratory_category) ~ NA,
        migratory_category %in% diadromous_cats ~ TRUE,
        TRUE ~ FALSE
      ),
      environment = case_when(
        fresh == 1 & saltwater == 1 ~ "both",
        fresh == 1 & saltwater == 0 ~ "freshwater",
        fresh == 0 & saltwater == 1 ~ "marine",
        TRUE ~ NA_character_
      )
    ) |>
    select(species, demers_pelag, is_diadromous, environment)

  invasive_undocumented <- c(
    "Neogobius melanostomus", "Ponticola kessleri", "Proterorhinus semilunaris"
  )

  status <- rfishbase::country(species) |>
    filter(country == country_name) |>
    select(species = Species, native_status = Status) |>
    distinct() |>
    mutate(is_native = native_status != "introduced") |>
    select(species, is_native)

  left_join(habitat, status, by = "species") |>
    mutate(
      is_native = dplyr::if_else(
        species %in% invasive_undocumented, FALSE, is_native
      )
    )
}
