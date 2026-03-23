#' Download the data from the spatiotemporal workshop from Zenodo.
#'
#' @param dir Directory where to store the zip file.
#' @return Path of the downloaded zip file.
#' @export
download_workshop_data <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  zen4R::download_zenodo("10.5281/zenodo.17962542", path = dir)
  file.path(dir, "miste_data.zip")
}

#' Unzip data from the spatiotemporal workshop.
#' Data is stored within zenodo folder.
#'
#' @param zip_file path of the zip file.
#' @param dir directory where to store the unzipped files.
#' @return Path of the unzipped directory.
#' @export
unzip_workshop <- function(zip_file, dir) {
  unzip(zip_file, exdir = dir)
  file.path(dir, "zenodo")
}

#' Create the plot of the DAG using ggplot.
#'
#' This represents the DAG we assume for our study.
#' That is, how environmental variable can shape food web structure.
#'
#' @return ggplot
#' @export
create_plot_dag <- function() {
  dag <- dagify(
    C ~ S + Comp. + Env.,
    S ~ Env.,
    Comp. ~ Env.,
    latent = "Comp.",
    exposure = "Env.",
    outcome = "C"
  )
  ggdag(dag) + theme_dag()
}
