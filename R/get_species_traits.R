#' Get species-level ecological traits from FishBase.
#'
#' Retrieves, for each species, its demersal/pelagic position, migratory
#' category and native/introduced status in `country_name`, plus a
#' freshwater/estuarine/marine/diadromous environment derived from FishBase's
#' Fresh/Brack/Saltwater flags.
#'
#' FishBase leaves `AnaCat` (migratory category) blank for most non-diadromous
#' species rather than coding it "non-migratory"; this is left as `NA` rather
#' than imputed - see the notebook for how it is interpreted downstream.
#'
#' `native_status` is collapsed to just "native"/"introduced": "endemic" and
#' "not established" are folded into "native" for simplicity, and the three
#' Ponto-Caspian gobies missing a France record in FishBase's `country()`
#' table (documented recent invaders in the literature) are hardcoded to
#' "introduced".
#'
#' @param species character vector of FishBase-validated species names.
#' @param country_name country to look up native/introduced status for.
#'
#' @return data frame with columns species, demers_pelag, migratory_category,
#' environment, native_status.
#' @export
get_species_traits <- function(species, country_name = "France") {
  habitat <- rfishbase::species(species) |>
    select(
      species = Species,
      fresh = Fresh,
      brack = Brack,
      saltwater = Saltwater,
      demers_pelag = DemersPelag,
      migratory_category = AnaCat
    ) |>
    mutate(
      migratory_category = dplyr::na_if(migratory_category, ""),
      # `brack` is a brackish-tolerance flag, not a habitat category on its
      # own: it is set for most freshwater and marine species too, so it
      # only distinguishes "estuarine" from the fresh/saltwater axis below.
      environment = case_when(
        fresh == 1 & saltwater == 1 ~ "diadromous",
        fresh == 1 & saltwater == 0 ~ "freshwater",
        fresh == 0 & saltwater == 1 ~ "marine",
        TRUE ~ "estuarine"
      )
    ) |>
    select(species, demers_pelag, migratory_category, environment)

  invasive_undocumented <- c(
    "Neogobius melanostomus", "Ponticola kessleri", "Proterorhinus semilunaris"
  )

  status <- rfishbase::country(species) |>
    filter(country == country_name) |>
    select(species = Species, native_status = Status) |>
    distinct() |>
    mutate(
      native_status = case_when(
        native_status %in% c("endemic", "not established") ~ "native",
        TRUE ~ native_status
      )
    )

  left_join(habitat, status, by = "species") |>
    mutate(
      native_status = dplyr::if_else(
        species %in% invasive_undocumented, "introduced", native_status
      )
    )
}
