test_that("merge_diet_river_marine lets freshwater win on overlapping species", {
  diet_marine <- tibble::tribble(
    ~species, ~stage, ~length_min, ~length_max, ~fish,
    "Sp A", "larvae", 0, 2, 0,
    "Sp A", "adult", 2, Inf, 1,
    "Sp B", "larvae", 0, Inf, 0
  )
  diet_river <- tibble::tribble(
    ~species, ~stage, ~length_min, ~length_max, ~fish,
    "Sp A", "larvae", 0, 3, 0,
    "Sp A", "adult", 3, Inf, 0,
    "Sp C", "larvae", 0, Inf, 0
  )
  merged <- merge_diet_river_marine(diet_marine, diet_river)

  sp_a <- merged[merged$species == "Sp A", ]
  expect_equal(nrow(sp_a), 2)
  expect_equal(sp_a$fish[sp_a$stage == "adult"], 0) # river value, not marine's 1

  expect_true("Sp B" %in% merged$species) # marine-only, kept
  expect_true("Sp C" %in% merged$species) # river-only, kept
})

test_that("merge_diet_river_marine closes small boundary buffers in the merged table", {
  diet_marine <- tibble::tibble(species = "Sp A", stage = "larvae", length_min = 0, length_max = Inf)
  diet_river <- tibble::tribble(
    ~species, ~stage, ~length_min, ~length_max,
    "Sp B", "larvae", 0, 3.0,
    "Sp B", "adult", 3.1, Inf
  )
  merged <- merge_diet_river_marine(diet_marine, diet_river)
  expect_equal(merged$length_max[merged$species == "Sp B" & merged$stage == "larvae"], 3.1)
})

test_that("merge_diet_river_marine keeps all species when there is no overlap", {
  diet_marine <- tibble::tibble(species = "Sp A", stage = "larvae", length_min = 0, length_max = Inf)
  diet_river <- tibble::tibble(species = "Sp B", stage = "larvae", length_min = 0, length_max = Inf)
  merged <- merge_diet_river_marine(diet_marine, diet_river)
  expect_setequal(merged$species, c("Sp A", "Sp B"))
})
