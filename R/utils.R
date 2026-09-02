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
  utils::unzip(zip_file, exdir = dir)
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
  dag <- ggdag::dagify(
    TL ~ S + Comp + Env,
    S ~ Env,
    Comp ~ Env,
    latent = "Comp",
    exposure = "Env",
    outcome = "TL",
    labels = c(
      TL = "Food web structure",
      S = "Richness",
      Comp = "Composition",
      Env = "Environment"
    )
  )
  ggdag::ggdag(dag, use_labels = "label") + ggdag::theme_dag()
}

#' Plot venn diagram of common species between surveys
#'
#' Styled to match `nice_theme()`: thin ink-grey outlines and labels, a
#' restrained blue fill, no legend (redundant with the in-plot counts).
#'
#' @param df data frame
#'
#' @return plot
#' @export
plot_venn_diagram <- function(df) {
  ink <- "grey20"
  # split() names each group after its own factor level - avoids mismatches
  # that group_split() + setNames(unique(df$survey)) can silently produce.
  distinct_df <- df |> distinct(species_valid, survey)
  sets <- split(distinct_df$species_valid, distinct_df$survey)

  ggVennDiagram::ggVennDiagram(
    sets,
    label = "count",
    label_geom = "text",
    label_color = ink,
    label_alpha = 1,
    label_size = 3.2,
    set_color = ink,
    set_size = 4,
    edge_size = 0.4
  ) +
    ggplot2::scale_fill_gradient(low = "white", high = "#a9c9e0") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.15)) +
    ggplot2::theme(legend.position = "none")
}

#' Get relative path.
#'
#' @param x path relative to the project root
#'
#' @return path relative to `output/`
#' @export
to_rel <- function(x) {
  file.path("..", x)
}

#' Plot sampling station sea and river survey.
#'
#' @param sea_data tidy data frame of sea data.
#' @param station_file file containing information about river station.
#'
#' @return plot.
#' @export
plot_sampling <- function(sea_data, station_file) {
  d_trait <- sea_data$trait
  d_sf <- d_trait |>
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
  france <- rnaturalearth::ne_countries(
    geounit = "france",
    type = "map_units",
    scale = "medium",
    returnclass = "sf"
  )
  rivers <- rnaturalearth::ne_download(
    scale = 10, type = "rivers_lake_centerlines", category = "physical",
    returnclass = "sf"
  )
  rivers_fr <- sf::st_intersection(rivers, france)
  d_station <- utils::read.csv(station_file, sep = ";") |>
    sf::st_as_sf(coords = c("lat", "long"), crs = 4326)
  d_sampling <- rbind(
    d_sf |> select(survey, geometry),
    d_station |> mutate(survey = "onema") |> select(survey, geometry)
  )
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = france) +
    ggplot2::geom_sf(data = rivers_fr, color = "grey") +
    ggplot2::geom_sf(data = d_sampling, ggplot2::aes(color = survey))
}

#' Map and temporal coverage of fishing operations by data source.
#'
#' Progress-report figure: a national map of stations colored by survey
#' (river operations relabelled `ASPE`), next to yearly operation counts per
#' survey. Both panels are restricted to `district_kept`, so the map matches
#' the catchments actually used in the analysis rather than the full,
#' national extent of each raw survey. The four kept catchments are drawn
#' as polygons on the map for spatial context.
#'
#' @param operation tibble with `longitude`, `latitude`, `survey`, `year`,
#'   `district`, e.g. `classify_operation_location()`'s output.
#' @param district_kept character vector of districts to keep.
#' @param basin sf tibble with `district`, `geometry`, e.g.
#'   `format_basin()`'s output.
#' @param lang "en" or "fr" axis labels.
#'
#' @return patchwork object combining the map and temporal-coverage panels.
#' @export
plot_station_overview <- function(
  operation, district_kept, basin, lang = c("en", "fr")
) {
  lang <- match.arg(lang)
  nice_theme() # self-contained: don't depend on caller-side theme_set()
  pal <- c(ASPE = "#2a78d6", Pomet = "#eb6834", Nurse = "#1baf7a")

  d <- operation |>
    filter(district %in% district_kept) |>
    mutate(
      source = recode(survey, pomet = "Pomet", nurse = "Nurse", river = "ASPE"),
      source = factor(source, levels = c("ASPE", "Pomet", "Nurse"))
    )
  d_sf <- d |> sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

  france <- rnaturalearth::ne_countries(
    geounit = "france", type = "map_units", scale = "medium",
    returnclass = "sf"
  )
  basin_kept <- basin |> filter(district %in% district_kept)

  big_text <- ggplot2::theme(
    axis.text = ggplot2::element_text(size = 13),
    axis.title = ggplot2::element_text(size = 15),
    legend.text = ggplot2::element_text(size = 14)
  )

  map_panel <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = france, fill = NA, color = "grey80", linewidth = 0.25) +
    ggplot2::geom_sf(data = basin_kept, fill = "grey97", color = "grey40", linewidth = 0.5) +
    ggplot2::geom_sf(
      data = d_sf |> arrange(source),
      ggplot2::aes(color = source), size = 0.3, alpha = 0.5
    ) +
    ggplot2::coord_sf(xlim = c(-5.3, 8.3), ylim = c(41.2, 51.2), expand = FALSE) +
    ggplot2::scale_color_manual(values = pal, name = NULL) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 3.5, alpha = 1))
    ) +
    ggplot2::labs(x = "Longitude", y = "Latitude") +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    big_text

  temporal_labs <- if (lang == "fr") {
    list(x = "Année", y = "Nombre d'opérations")
  } else {
    list(x = "Year", y = "Number of operations")
  }

  temporal_panel <- d |>
    count(source, year) |>
    ggplot2::ggplot(ggplot2::aes(year, n, fill = source)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~source, ncol = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::labs(x = temporal_labs$x, y = temporal_labs$y) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank()
    ) +
    big_text

  ((map_panel + temporal_panel) +
    patchwork::plot_layout(widths = c(1.3, 1), guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A")) &
    ggplot2::theme(
      legend.position = "bottom",
      plot.tag = ggplot2::element_text(size = 15, face = "bold")
    )
}

#' Map of food web structure metrics along the continuum.
#'
#' Four small maps - species richness, trophic chain length, connectance,
#' diet overlap - one point per local food web (operation), restricted to
#' `district_kept`. Each metric gets its own panel/color scale (very
#' different units and ranges), same basin/coastline context as
#' `plot_station_overview()`.
#'
#' @param foodweb_structure tibble with `operation_id`, `longitude`,
#' `latitude`, `log_species_richness`, `trophic_length`, `connectance`,
#' `diet_overlap`, e.g. `measure_foodweb_structure()`'s output.
#' @param operation_location tibble with `operation_id`, `district`, e.g.
#' `classify_operation_location()`'s output.
#' @param district_kept character vector of districts to keep.
#' @param basin sf tibble with `district`, `geometry`, e.g.
#' `format_basin()`'s output.
#'
#' @return patchwork object, four panels.
#' @export
plot_foodweb_structure_map <- function(
  foodweb_structure, operation_location, district_kept, basin
) {
  nice_theme()

  d <- foodweb_structure |>
    left_join(
      operation_location |> select(operation_id, district),
      by = "operation_id"
    ) |>
    filter(district %in% district_kept) |>
    mutate(richness = exp(log_species_richness)) |>
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

  france <- rnaturalearth::ne_countries(
    geounit = "france", type = "map_units", scale = "medium",
    returnclass = "sf"
  )
  basin_kept <- basin |> filter(district %in% district_kept)

  metrics <- c(
    richness = "Richesse spécifique",
    trophic_length = "Longueur de chaîne trophique",
    connectance = "Connectance",
    diet_overlap = "Similarité des régimes alimentaires"
  )

  panels <- lapply(names(metrics), function(var) {
    ggplot2::ggplot() +
      ggplot2::geom_sf(
        data = france, fill = NA, color = "grey80", linewidth = 0.25
      ) +
      ggplot2::geom_sf(
        data = basin_kept, fill = "grey97", color = "grey40", linewidth = 0.4
      ) +
      ggplot2::geom_sf(
        data = d, ggplot2::aes(color = .data[[var]]), size = 0.3, alpha = 0.6
      ) +
      ggplot2::coord_sf(xlim = c(-5.3, 8.3), ylim = c(41.2, 51.2), expand = FALSE) +
      ggplot2::scale_color_viridis_c(
        name = NULL, breaks = scales::breaks_pretty(n = 3)
      ) +
      ggplot2::labs(x = "Longitude", y = "Latitude", title = metrics[[var]]) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 13, face = "bold"),
        axis.text = ggplot2::element_text(size = 9),
        axis.title = ggplot2::element_text(size = 10),
        legend.text = ggplot2::element_text(size = 9)
      )
  })

  patchwork::wrap_plots(panels, ncol = 2)
}

#' Diet-category composition of each survey's fish community.
#'
#' Stacked bar of prey categories represented in the diet of the species
#' caught in each survey (river operations relabelled `ASPE`, matching the
#' other report figures). A species counts toward a category if any of its
#' life stages consume it.
#'
#' @param diet tibble with `species`, `stage`, and one 0/1 column per prey
#'   category, e.g. `merge_diet()`'s output.
#' @param size tibble with `species_valid`, `survey`, e.g.
#'   `size_year_filtered`.
#'
#' @return ggplot.
#' @export
plot_diet_categories <- function(diet, size) {
  nice_theme()

  category_labels <- c(
    detritus = "Détritus", biofilm = "Biofilm", macrophyte = "Macrophytes",
    phytoplankton = "Phytoplancton", zooplankton = "Zooplancton",
    worm = "Vers", mollusk = "Mollusques", crustacean = "Crustacés",
    insect = "Insectes", echinoderm = "Échinodermes", fish = "Poissons"
  )
  category_cols <- names(category_labels)

  species_diet <- diet |>
    group_by(species) |>
    summarise(
      across(all_of(category_cols), ~ as.integer(any(. == 1))),
      .groups = "drop"
    )

  d <- size |>
    distinct(species_valid, survey) |>
    mutate(
      source = recode(survey, pomet = "Pomet", nurse = "Nurse", river = "ASPE"),
      source = factor(source, levels = c("ASPE", "Pomet", "Nurse"))
    ) |>
    inner_join(species_diet, by = c("species_valid" = "species")) |>
    tidyr::pivot_longer(
      all_of(category_cols), names_to = "category", values_to = "present"
    ) |>
    filter(present == 1) |>
    count(source, category) |>
    group_by(source) |>
    mutate(prop = n / sum(n)) |>
    ungroup() |>
    mutate(
      category = factor(
        unname(category_labels[category]),
        levels = category_labels[category_cols]
      )
    )

  ggplot2::ggplot(d, ggplot2::aes(source, prop, fill = category)) +
    ggplot2::geom_col() +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::scale_fill_brewer(palette = "Paired", name = "Catégorie de proie") +
    ggplot2::labs(x = NULL, y = "Proportion d'espèces") +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 13),
      legend.text = ggplot2::element_text(size = 11),
      legend.title = ggplot2::element_text(size = 12)
    )
}

plot_network_groups <- function(diet_resource) {
  adj <- t(as.matrix(diet_resource |> select(-c(light, species, reference))))
  colnames(adj) <- diet_resource$species
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "directed")
  metanetwork::build_metanet(metaweb = g) |>
    metanetwork::compute_TL() |>
    metanetwork::ggmetanet()
}

read_environment_data <- function(dir, year_start = 1993, year_end = 2023) {
  year_vals <- seq(year_start, year_end)
  purrr::map(year_vals, function(year) {
    utils::read.csv(here::here(dir, paste0("Environment_QM_", year, ".csv"))) |>
      mutate(year = year)
  }) |>
    purrr::list_rbind()
}

download_coastal_water <- function(dir) {
  url_sea_polygones <- c(
    estuary = "https://services.sandre.eaufrance.fr/telechargement/geo/MDO/MasseDEau/vedl_2019/PolygMasseDEauTransition/PolygMasseDEauTransition_FXX-shp.zip",
    coast = "https://services.sandre.eaufrance.fr/telechargement/geo/MDO/MasseDEau/vrap_2022/MasseDEauCotiere/MasseDEauCotiere_FRA-shp.zip"
  )
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  purrr::walk2(
    url_sea_polygones, names(url_sea_polygones),
    \(url, name) {
      zip_file <- here::here(dir, paste0(name, ".zip"))
      utils::download.file(url, destfile = zip_file)
      utils::unzip(zipfile = zip_file, exdir = here::here(dir))
      file.remove(zip_file)
    }
  )
  dir
}

format_coastal_water <- function(dir) {
  # District classification from
  # https://www.donnees.statistiques.developpement-durable.gouv.fr/lesessentiels/essentiels/dce-district.html
  file_list <- list.files(path = dir, pattern = ".shp")
  purrr::map(
    here::here(dir, file_list),
    \(x) sf::read_sf(x) |> dplyr::select(CdEuMasseD)
  ) |>
    purrr::map(\(x) sf::st_transform(x, crs = 4326)) |>
    bind_rows() |>
    rename(zone_id = CdEuMasseD) |>
    sf::st_make_valid() |>
    mutate(district = case_when(
      str_detect(zone_id, "^FRA") ~ "Escaut-Somme",
      str_detect(zone_id, "^FRB") ~ "Meuse",
      str_detect(zone_id, "^FRC") ~ "Rhin",
      str_detect(zone_id, "^FRD") ~ "Rhone",
      str_detect(zone_id, "^FRE") ~ "Corse",
      str_detect(zone_id, "^FRF") ~ "Ardour-Garonne",
      str_detect(zone_id, "^FRG") ~ "Loire",
      str_detect(zone_id, "^FRH") ~ "Seine",
      TRUE ~ NA_character_
    ))
}

download_basin <- function(dir) {
  url <- "https://services.sandre.eaufrance.fr/telechargement/geo/ETH/BDTopage/2025/BassinHydrographique/BassinHydrographique_FXX-shp.zip"
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  zip_file <- here::here(dir, "bassin.zip")
  utils::download.file(url, destfile = zip_file)
  utils::unzip(zipfile = zip_file, exdir = dir)
  file.remove(zip_file)
  dir
}

format_basin <- function(dir) {
  # LbBH is the basin's official name (BD Topage, SANDRE) - mapped to this
  # project's district scheme. Rhin-Meuse is a single combined polygon
  # upstream (unlike the FRB/FRC split used elsewhere) - not an issue since
  # neither Rhin nor Meuse is a kept district.
  file_list <- list.files(path = dir, pattern = ".shp")
  sf::read_sf(here::here(dir, file_list)) |>
    sf::st_transform(crs = 4326) |>
    sf::st_make_valid() |>
    mutate(district = case_when(
      LbBH == "Artois-Picardie" ~ "Escaut-Somme",
      LbBH == "Seine-Normandie" ~ "Seine",
      LbBH == "Loire-Bretagne" ~ "Loire",
      LbBH == "Adour-Garonne" ~ "Ardour-Garonne",
      LbBH == "Rhône-Méditerranée" ~ "Rhone",
      LbBH == "Corse" ~ "Corse",
      LbBH == "Rhin-Meuse" ~ "Rhin-Meuse",
      TRUE ~ NA_character_
    )) |>
    select(district, geometry)
}

#' Assign districts and whether the operation is inland or sea-side.
#' and restrict `river` operations to the catchments of interest.
#'
#' Inland test if operation is in a basin. If they miss, they are sea-side.
#' Inland operations catchment come from the basin they fall in.
#' Sea-side operations catchment come from the nearest coastal water they are
#' to.
#'
#' @param operation tibble with `longitude`, `latitude`, `survey`.
#' @param basin sf tibble with `district`, `geometry`.
#' @param coastal_water sf tibble with `district`, `geometry`.
#' @param district_kept character vector of districts to keep.
#' @param dist_max_sea maximum distance, in km, to fall back to the
#'   nearest `coastal_water` polygon for sea-side operations.
#'
#' @return `operation` with `district` and `inland` added (TRUE if inside a
#' basin polygon), `river` rows restricted to `district_kept`.
#' @export
classify_operation_location <- function(
  operation,
  basin,
  coastal_water,
  district_kept,
  dist_max_sea
) {
  basin_valid <- basin |> filter(!is.na(district))
  coastal_water_valid <- coastal_water |> filter(!is.na(district))
  op_geom <- operation |>
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
    sf::st_geometry()

  district_basin <- in_district(op_geom, basin_valid)
  # Operation 61837633 matches the basin polygon at the coast but is
  # sea-side - an implausible ~8.2km network snap distance gives it away.
  district_basin[operation$operation_id == "61837633"] <- NA_character_
  inland <- !is.na(district_basin)
  district_coastal <- rep(NA_character_, length(op_geom))
  district_coastal[!inland] <- nearest_district(
    op_geom[!inland],
    coastal_water_valid,
    dist_max = dist_max_sea
  )

  operation |>
    mutate(
      inland = inland,
      district = dplyr::coalesce(district_basin, district_coastal)
    ) |>
    filter(survey != "river" | district %in% district_kept)
}

#' Get district of the polygon each point falls within.
#'
#' @param points sf or sfc points.
#' @param hydro_zone sf polygons with `district`.
#'
#' @return character vector of `district`, one per row of `points`, `NA`
#' outside every polygon.
in_district <- function(points, hydro_zone) {
  sf::st_join(
    sf::st_as_sf(points),
    hydro_zone["district"],
    join = sf::st_within
  )$district
}

#' District of the nearest polygon, within a distance threshold.
#'
#' For a pure containment check (no fallback), use `in_district()`
#' instead - cheaper on complex polygons, see its docstring for why.
#'
#' @param points sf points.
#' @param hydro_zone sf polygons with `district`.
#' @param dist_max maximum distance to accept a match, in km.
#'
#' @return character vector of `district`, one per row of `points`, `NA`
#' beyond `dist_max`.
nearest_district <- function(points, hydro_zone, dist_max) {
  nearest_idx <- sf::st_nearest_feature(points, hydro_zone)
  dist_km <- as.numeric(sf::st_distance(
    points,
    hydro_zone[nearest_idx, ],
    by_element = TRUE
  )) / 1000
  ifelse(dist_km <= dist_max, hydro_zone$district[nearest_idx], NA_character_)
}

match_sea_district <- function(foodweb, hydro_zone, dist_max = 5) {
  fw_sf <- foodweb |>
    select(-c(
      foodweb,
      log_trophic_richness,
      log_species_richness,
      connectance,
      trophic_length
    )) |>
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
    sf::st_join(hydro_zone, join = sf::st_within)

  # Food web within zone.
  fw_within <- fw_sf |>
    filter(!is.na(zone_id))
  fw_outside <- fw_sf |>
    filter(is.na(zone_id))

  # Match outside point to nearest district.
  fw_nearest <- fw_outside |>
    mutate(district = nearest_district(geometry, hydro_zone, dist_max)) |>
    filter(!is.na(district))

  fw_within |>
    rbind(fw_nearest) |>
    mutate(
      longitude = sf::st_coordinates(geometry)[, 1],
      latitude  = sf::st_coordinates(geometry)[, 2]
    ) |>
    sf::st_drop_geometry()
}

#' Convert a dataframe to sf object
#'
#' Assumes coordinates are given by the longitude and latitude columns.
#'
#' @param x data.frame
#'
#' @return sf object
#' @export
to_sf <- function(x) {
  x |> sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
}

#' Plot data points on the France map
#'
#' @param data data frame
#'
#' @return plot by ggplot
#' @export
plot_data_on_map <- function(data) {
  coastline <- rnaturalearth::ne_coastline(scale = "large", returnclass = "sf")
  ggplot2::ggplot(data, ggplot2::aes(x = longitude, y = latitude)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::geom_sf(data = coastline, inherit.aes = FALSE, linewidth = 0.3) +
    ggplot2::coord_sf(
      xlim = range(data$longitude),
      ylim = range(data$latitude)
    ) +
    ggplot2::labs(x = NULL, y = NULL)
}

#' Theme for plots of the RIVERSEA project
#'
#' A clean, minimalist theme suitable for scientific publication: sober
#' typography, thin axis lines instead of a panel border, faint horizontal
#' gridlines only, and unobtrusive facet strips and legends.
#'
#' @param base_size base font size, in points.
#'
#' @return nothing, sets the theme globally via [ggplot2::theme_set()].
#' @export
nice_theme <- function(base_size = 11) {
  ink <- "grey20"
  muted <- "grey45"
  faint <- "grey85"
  ggplot2::theme_set(
    ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(
        # Text
        plot.title = ggplot2::element_text(
          face = "bold", size = ggplot2::rel(1.15), color = ink,
          margin = ggplot2::margin(b = 4)
        ),
        plot.subtitle = ggplot2::element_text(
          color = muted, size = ggplot2::rel(0.95),
          margin = ggplot2::margin(b = 10)
        ),
        plot.caption = ggplot2::element_text(
          color = muted, size = ggplot2::rel(0.75), hjust = 0
        ),
        axis.title = ggplot2::element_text(color = ink, size = ggplot2::rel(0.95)),
        axis.text = ggplot2::element_text(color = muted, size = ggplot2::rel(0.85)),
        # Axes and grid: thin axis lines instead of a panel border, only
        # faint horizontal gridlines to guide the eye across the plot.
        axis.line = ggplot2::element_line(color = ink, linewidth = 0.3),
        axis.ticks = ggplot2::element_line(color = ink, linewidth = 0.3),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(color = faint, linewidth = 0.3),
        # Facets
        strip.text = ggplot2::element_text(
          face = "bold", color = ink, size = ggplot2::rel(0.9)
        ),
        strip.background = ggplot2::element_blank(),
        # Legend
        legend.position = "bottom",
        legend.title = ggplot2::element_text(color = ink, size = ggplot2::rel(0.9)),
        legend.text = ggplot2::element_text(color = muted, size = ggplot2::rel(0.85)),
        legend.key = ggplot2::element_blank(),
        # Backgrounds and margins
        plot.margin = ggplot2::margin(t = 10, r = 15, b = 10, l = 15)
      )
  )
}

download_sea_raw_data <- function(dest_dir = "data/sea/raw") {
  parent_dir <- dirname(dest_dir)
  dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)

  zip_path <- file.path(parent_dir, "sea.zip")

  zen4R::download_zenodo(
    doi = "10.5281/zenodo.21395668",
    path = parent_dir,
    files = "sea.zip"
  )

  utils::unzip(zip_path, exdir = parent_dir)
  file.remove(zip_path)

  dest_dir
}

#' Get coastline from rnaturalearth package.
#'
#'
#' @return sf object
#' @export
coastline <- function() {
  rnaturalearth::ne_coastline(scale = "large", returnclass = "sf")
}

#' Get land polygons (islands included) from rnaturalearth package.
#'
#' @return sf polygons, WGS84.
#' @export
land_polygon <- function() {
  rnaturalearth::ne_download(scale = 10, type = "land", category = "physical", returnclass = "sf")
}

#' Filter data such that year is in the specified time window
#'
#' @param data data tibble with `year` column.
#' @param year_min minimum year to be kept.
#' @param year_max maximal year to be kept.
#'
#' @return data tibble filtered.
#' @export
filter_year <- function(data, year_min, year_max) {
  data |>
    filter(year >= year_min, year <= year_max)
}

#' Convert date to year and month
#'
#' @param data tibble with `date` column.
#' @param keep logical, whether or not to keep the original `date` column.
#'
#' @return tibble
#' @export
from_date_to_year_month <- function(data, keep = FALSE) {
  data <- data |>
    mutate(
      year = lubridate::year(date),
      month = lubridate::month(date)
    )
  if (keep) data else data |> select(-date)
}
