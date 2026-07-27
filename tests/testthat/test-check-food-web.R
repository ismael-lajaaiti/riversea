build_test_metaweb <- function() {
  nodes <- c("phytoplankton", "fishA", "fishB", "fishC", "fishD")
  m <- matrix(0, nrow = 5, ncol = 5, dimnames = list(nodes, nodes))
  m["phytoplankton", "fishA"] <- 1 # fishA eats phytoplankton
  m["fishA", "fishB"] <- 1 # fishB eats fishA
  m["fishC", "fishB"] <- 1 # fishB also eats fishC, but fishC itself eats nothing
  # fishD is left fully disconnected: no prey, no predators
  m
}

test_that("check_metaweb_consistency flags a consumer with no prey", {
  violations <- check_metaweb_consistency(build_test_metaweb(), resource_list = "phytoplankton")
  expect_true(any(violations$species == "fishC" & violations$issue == "consumer_with_no_prey"))
})

test_that("check_metaweb_consistency flags a fully isolated node", {
  violations <- check_metaweb_consistency(build_test_metaweb(), resource_list = "phytoplankton")
  expect_true(any(violations$species == "fishD" & violations$issue == "isolated_node"))
})

test_that("check_metaweb_consistency does not double-flag an isolated node as a consumer with no prey", {
  violations <- check_metaweb_consistency(build_test_metaweb(), resource_list = "phytoplankton")
  expect_false(any(violations$species == "fishD" & violations$issue == "consumer_with_no_prey"))
})

test_that("check_metaweb_consistency does not flag resources, apex predators or normal consumers", {
  violations <- check_metaweb_consistency(build_test_metaweb(), resource_list = "phytoplankton")
  expect_false(any(violations$species %in% c("phytoplankton", "fishA", "fishB")))
})

test_that("check_metaweb_consistency passes on a fully consistent metaweb", {
  nodes <- c("phytoplankton", "fishA", "fishB")
  m <- matrix(0, nrow = 3, ncol = 3, dimnames = list(nodes, nodes))
  m["phytoplankton", "fishA"] <- 1
  m["fishA", "fishB"] <- 1
  expect_equal(nrow(check_metaweb_consistency(m, resource_list = "phytoplankton")), 0)
})

test_that("check_metaweb_consistency rejects a non-square matrix", {
  m <- matrix(0, nrow = 2, ncol = 3)
  expect_error(check_metaweb_consistency(m, resource_list = character()))
})

test_that("check_metaweb_consistency rejects non-binary values", {
  nodes <- c("a", "b")
  m <- matrix(c(0, 2, 0, 0), nrow = 2, dimnames = list(nodes, nodes))
  expect_error(check_metaweb_consistency(m, resource_list = character()))
})

test_that("assert_metaweb_consistency stops with a report when checks fail", {
  expect_error(assert_metaweb_consistency(build_test_metaweb(), resource_list = "phytoplankton"))
})

test_that("assert_metaweb_consistency returns TRUE invisibly when checks pass", {
  nodes <- c("phytoplankton", "fishA", "fishB")
  m <- matrix(0, nrow = 3, ncol = 3, dimnames = list(nodes, nodes))
  m["phytoplankton", "fishA"] <- 1
  m["fishA", "fishB"] <- 1
  expect_true(assert_metaweb_consistency(m, resource_list = "phytoplankton"))
})
