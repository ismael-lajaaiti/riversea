download_workshop_data <- function(dir) {
    if (!dir.exists(dir)) {
        dir.create(dir)
    }
    zen4R::download_zenodo("10.5281/zenodo.17962542", path = dir)
    zip_path <- file.path(dir, "miste_data.zip")
    return(zip_path)
}

unzip_workshop <- function(zip_file, dir) {
    unzip(zip_file, exdir = dir)
    return(file.path(dir, "zenodo"))
}

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
