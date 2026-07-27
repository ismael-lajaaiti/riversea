test_that("check_diet_species_coverage flags a commonly observed species missing from the diet table", {
  diet <- tibble::tibble(species = c("Abramis brama", "Alosa alosa"))
  catch <- tibble::tibble(
    species_valid = c(
      rep("Abramis brama", 15),
      rep("Alosa alosa", 15),
      rep("Cyprinus carpio", 15)
    ),
    trait = rep(1:15, 3)
  )
  missing <- check_diet_species_coverage(diet, catch, occurence_min = 10)
  expect_equal(missing, "Cyprinus carpio")
})

test_that("check_diet_species_coverage passes when every commonly observed species is covered", {
  diet <- tibble::tibble(species = c("Abramis brama", "Alosa alosa"))
  catch <- tibble::tibble(
    species_valid = c(rep("Abramis brama", 15), rep("Alosa alosa", 15)),
    trait = rep(1:15, 2)
  )
  expect_length(check_diet_species_coverage(diet, catch, occurence_min = 10), 0)
})

test_that("check_diet_species_coverage uses a strict threshold", {
  diet <- tibble::tibble(species = character())
  catch <- tibble::tibble(species_valid = rep("Rare species", 10), trait = 1:10) # exactly at threshold
  expect_length(check_diet_species_coverage(diet, catch, occurence_min = 10), 0)
})

test_that("check_diet_species_coverage ignores unmatched (NA) species names", {
  diet <- tibble::tibble(species = character())
  catch <- tibble::tibble(species_valid = rep(NA_character_, 15), trait = 1:15)
  expect_length(check_diet_species_coverage(diet, catch, occurence_min = 10), 0)
})

test_that("check_diet_species_coverage counts distinct traits, not raw catch rows", {
  diet <- tibble::tibble(species = character())
  # 20 catch rows (e.g. one per size class) but only 3 distinct traits:
  # occurrence should be 3, not 20, and stay below the threshold.
  catch <- tibble::tibble(
    species_valid = rep("Common in one trait", 20),
    trait = rep(1:3, length.out = 20)
  )
  expect_length(check_diet_species_coverage(diet, catch, occurence_min = 10), 0)
})

test_that("close_diet_size_gaps closes a small boundary buffer exactly", {
  diet <- tibble::tribble(
    ~species, ~stage, ~length_min, ~length_max,
    "Sp A", "larvae", 0, 3.0,
    "Sp A", "adult", 3.1, Inf
  )
  closed <- close_diet_size_gaps(diet)
  expect_equal(closed$length_max[closed$stage == "larvae"], 3.1)
  # a size-class midpoint sitting exactly on the old gap is now covered
  expect_true(any(3.05 >= closed$length_min & 3.05 < closed$length_max))
})

test_that("close_diet_size_gaps leaves a real, larger gap untouched", {
  diet <- tibble::tribble(
    ~species, ~stage, ~length_min, ~length_max,
    "Sp A", "larvae", 0, 2,
    "Sp A", "adult", 10, Inf
  )
  closed <- close_diet_size_gaps(diet)
  expect_equal(closed$length_max[closed$stage == "larvae"], 2)
})

test_that("check_diet_size_coverage flags a gap between stages", {
  diet <- tibble::tribble(
    ~species, ~stage, ~length_min, ~length_max,
    "Sp A", "larvae", 0, 2,
    "Sp A", "adult", 10, Inf
  )
  issues <- check_diet_size_coverage(diet)
  expect_equal(issues$issue, "size_gap")
})

test_that("check_diet_size_coverage tolerates a small literature boundary buffer", {
  diet <- tibble::tribble(
    ~species, ~stage, ~length_min, ~length_max,
    "Sp A", "larvae", 0, 2,
    "Sp A", "adult", 2.1, Inf
  )
  expect_equal(nrow(check_diet_size_coverage(diet)), 0)
})

test_that("check_diet_size_coverage flags a species whose last stage does not reach infinity", {
  diet <- tibble::tribble(
    ~species, ~stage, ~length_min, ~length_max,
    "Sp A", "larvae", 0, 30
  )
  issues <- check_diet_size_coverage(diet)
  expect_true("does_not_reach_infinity" %in% issues$issue)
})

test_that("check_diet_size_coverage flags a species not starting at zero", {
  diet <- tibble::tribble(
    ~species, ~stage, ~length_min, ~length_max,
    "Sp A", "larvae", 2, Inf
  )
  issues <- check_diet_size_coverage(diet)
  expect_true("does_not_start_at_zero" %in% issues$issue)
})

test_that("assert_diet_coverage stops with a combined message when checks fail", {
  diet <- tibble::tribble(
    ~species, ~stage, ~length_min, ~length_max,
    "Sp A", "larvae", 0, 2,
    "Sp A", "adult", 10, Inf
  )
  catch <- tibble::tibble(species_valid = rep("Sp B", 15), trait = 1:15)
  expect_error(assert_diet_coverage(diet, catch, occurence_min = 10))
})

test_that("assert_diet_coverage returns TRUE invisibly when checks pass", {
  diet <- tibble::tribble(
    ~species, ~stage, ~length_min, ~length_max,
    "Sp A", "larvae", 0, 2.000001,
    "Sp A", "adult", 2.000002, Inf
  )
  catch <- tibble::tibble(species_valid = rep("Sp A", 15), trait = 1:15)
  expect_true(assert_diet_coverage(diet, catch, occurence_min = 10))
})

test_that("check_diet_size_coverage passes on a contiguous, well-formed diet table", {
  diet <- tibble::tribble(
    ~species, ~stage, ~length_min, ~length_max,
    "Sp A", "larvae", 0, 2.000001,
    "Sp A", "juvenile", 2.000002, 9.000001,
    "Sp A", "adult", 9.000002, Inf
  )
  expect_equal(nrow(check_diet_size_coverage(diet)), 0)
})
