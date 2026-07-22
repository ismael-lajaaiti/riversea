options(
  renv.config.repos.override = c(
    CRAN = "https://packagemanager.posit.co/cran/latest",
    INLA = "https://inla.r-inla-download.org/R/stable",
    rOpenSci = "https://ropensci.r-universe.dev"
  ),
  renv.config.pak.enabled = TRUE
)
source("renv/activate.R")

# Must run after renv/activate.R: when the renv autoloader is disabled (as in
# CI), activate.R re-sources ~/.Rprofile, which can reset options(repos) and
# drop the INLA repo needed to resolve DESCRIPTION's Imports via pak.
options(
  repos = c(
    CRAN = "https://packagemanager.posit.co/cran/latest",
    INLA = "https://inla.r-inla-download.org/R/stable",
    rOpenSci = "https://ropensci.r-universe.dev"
  )
)
