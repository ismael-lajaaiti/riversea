#' Check the internal consistency of a metaweb adjacency matrix.
#'
#' Flags two failure modes: a trophic species that is fully disconnected, and a
#' non-resource (fish) trophic species with no prey at all, which would starve
#' in any local food web built from this metaweb.
#'
#' @param metaweb square 0/1 adjacency matrix, rows are prey and columns are
#' predators.
#' @param resource_list character vector of basal resource node names.
#'
#' @return tibble with columns species, issue. Empty if the metaweb passes
#' both checks.
#' @export
check_metaweb_consistency <- function(metaweb, resource_list) {
  stopifnot(nrow(metaweb) == ncol(metaweb), all(metaweb %in% c(0, 1)))

  row_sums <- rowSums(metaweb)
  col_sums <- colSums(metaweb)

  isolated <- rownames(metaweb)[row_sums == 0 & col_sums == 0]
  no_prey <- setdiff(
    colnames(metaweb)[col_sums == 0], c(resource_list, isolated)
  )

  bind_rows(
    tibble(species = isolated, issue = "isolated_node"),
    tibble(species = no_prey, issue = "consumer_with_no_prey")
  )
}

#' Assert that a metaweb passes the consistency check.
#'
#' Wraps `check_metaweb_consistency()`, stopping with a report if any
#' violation is found. Intended as a single gate call in the targets
#' pipeline.
#'
#' @param metaweb square 0/1 adjacency matrix, rows are prey and columns are
#' predators.
#' @param resource_list character vector of basal resource node names.
#'
#' @return TRUE, invisibly, if the metaweb passes the check.
#' @export
assert_metaweb_consistency <- function(metaweb, resource_list) {
  violations <- check_metaweb_consistency(metaweb, resource_list)
  if (nrow(violations) > 0) {
    stop(
      "Metaweb failed consistency check:\n- ",
      paste(
        violations$species, violations$issue,
        sep = " - ", collapse = "\n- "
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
