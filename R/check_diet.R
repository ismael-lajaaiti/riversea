#' Check that the diet table covers every commonly observed species.
#'
#' A species is considered commonly observed if its occurrence is strictly
#' greater than `occurence_min`, where occurrence is the number of distinct
#' traits (fishing operations) in which the species was caught, not the raw
#' number of catch rows (a species can have several rows within the same
#' trait, e.g. one per size class).
#'
#' @param diet diet table with a `species` column.
#' @param catch sea survey catch data with `species_valid` and `trait`
#' columns.
#' @param occurence_min occurrence threshold.
#'
#' @return character vector of commonly observed species missing from
#' `diet`. Empty if the diet table covers every commonly observed species.
#' @export
check_diet_species_coverage <- function(diet, catch, occurence_min) {
  species_observed <- catch |>
    filter(!is.na(species_valid)) |>
    group_by(species_valid) |>
    summarise(occurence = n_distinct(trait), .groups = "drop") |>
    filter(occurence > occurence_min) |>
    pull(species_valid)
  setdiff(species_observed, unique(diet$species))
}

#' Tolerance (in centimeters) for the small boundary buffer literature
#' sources use between diet stages, e.g. a stage ending at 3.0 cm and the
#' next one starting at 3.1 cm. Below this, `close_diet_size_gaps()`
#' treats a gap as formatting rather than a real missing-diet interval.
#' @noRd
diet_gap_tolerance <- 0.15

#' Close small boundary buffers between consecutive diet stages.
#'
#' `foodwebbuilder::build_metaweb()` matches a fish's size to a diet
#' interval with `size >= length_min & size < length_max`. A tiny buffer
#' between stages (e.g. `length_max` of 3.0 followed by the next stage's
#' `length_min` of 3.1) is common in the source literature tables and looks
#' harmless, but a size-class midpoint that lands exactly in that 0.1 cm gap
#' matches neither interval and makes `build_metaweb()` error out. This
#' closes any gap up to `diet_gap_tolerance` by snapping a stage's
#' `length_max` to the next stage's `length_min`, so stages touch exactly
#' with no dead zone. Gaps larger than the tolerance are left untouched, so
#' they still surface as real coverage problems in
#' `check_diet_size_coverage()`.
#'
#' @param diet diet table with species, stage, length_min, length_max
#' columns.
#'
#' @return diet table with small inter-stage gaps closed.
#' @export
close_diet_size_gaps <- function(diet) {
  diet |>
    arrange(species, length_min) |>
    group_by(species) |>
    mutate(
      next_min = lead(length_min),
      length_max = if_else(
        !is.na(next_min) & (next_min - length_max) <= diet_gap_tolerance,
        next_min,
        length_max
      )
    ) |>
    ungroup() |>
    select(-next_min)
}

#' Check that each species' diet is defined over a contiguous size range.
#'
#' Flags three failure modes: a gap between consecutive stages large enough
#' (greater than `diet_gap_tolerance`, to accommodate the small boundary
#' buffer literature sources use between stages) that a fish could fall
#' into it, a species whose first stage does not start at size 0, and a
#' species whose last stage does not reach an unbounded size.
#' `length_min`/`length_max` are assumed to be in centimeters. Run
#' `close_diet_size_gaps()` first to avoid flagging (and to stop
#' `build_metaweb()` erroring on) inter-stage buffers that are just
#' formatting.
#'
#' @param diet diet table with species, stage, length_min, length_max
#' columns.
#'
#' @return tibble with columns species, issue, detail. Empty if every
#' species' diet is defined over its whole size range.
#' @export
check_diet_size_coverage <- function(diet) {
  by_species <- diet |>
    arrange(species, length_min) |>
    group_by(species)

  gaps <- by_species |>
    mutate(gap = length_min - lag(length_max)) |>
    ungroup() |>
    filter(!is.na(gap) & gap > diet_gap_tolerance) |>
    transmute(
      species,
      issue = "size_gap",
      detail = paste0("gap of ", round(gap, 3), " before stage '", stage, "'")
    )

  not_zero <- diet |>
    group_by(species) |>
    summarise(min_length = min(length_min), .groups = "drop") |>
    filter(min_length > 1e-3) |>
    transmute(
      species,
      issue = "does_not_start_at_zero",
      detail = paste0("starts at ", min_length)
    )

  not_infinite <- diet |>
    group_by(species) |>
    summarise(max_length = max(length_max), .groups = "drop") |>
    filter(is.finite(max_length)) |>
    transmute(
      species,
      issue = "does_not_reach_infinity",
      detail = paste0("stops at ", max_length)
    )

  bind_rows(gaps, not_zero, not_infinite)
}

#' Assert that the diet table passes all coverage checks.
#'
#' Wraps `check_diet_species_coverage()` and `check_diet_size_coverage()`,
#' stopping with a combined report if either finds an error.
#'
#' @param diet diet table with species, stage, length_min, length_max
#' columns.
#' @param catch sea survey catch data with `species_valid` and `trait`
#' columns.
#' @param occurence_min occurrence threshold.
#'
#' @return TRUE, invisibly, if the diet table passes every check.
#' @export
assert_diet_coverage <- function(diet, catch, occurence_min) {
  missing_species <- check_diet_species_coverage(diet, catch, occurence_min)
  size_issues <- check_diet_size_coverage(diet)

  problems <- c(
    if (length(missing_species) > 0) {
      paste0("missing species: ", paste(missing_species, collapse = ", "))
    },
    if (nrow(size_issues) > 0) {
      paste0(
        "size coverage: ",
        paste(
          size_issues$species,
          size_issues$issue,
          sep = " - ",
          collapse = "; "
        )
      )
    }
  )
  if (length(problems) > 0) {
    stop(
      "Diet table failed coverage checks:\n- ",
      paste(problems, collapse = "\n- "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
