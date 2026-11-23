#' Download the HydroSHEDS flow directions of Europe.
#'
#' HydroSHEDS v1 (Lehner et al. 2008) derives, from the SRTM elevation model,
#' the direction in which water flows out of each 15 arc-second cell (~450 m
#' x 310 m in France), towards one of its 8 neighbours (codes 1, 2, 4, ...,
#' 128 for east, south-east, south, ..., north-east). The download is
#' skipped if the file exists.
#'
#' @param dir directory where to store the zip file.
#' @return Path of `hyd_eu_dir_15s.zip`.
#' @export
download_hydrosheds <- function(dir) {
  path <- file.path(dir, "hyd_eu_dir_15s.zip")
  if (!file.exists(path)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    curl::curl_download("https://data.hydrosheds.org/file/hydrosheds-v1-dir/hyd_eu_dir_15s.zip", path)
  }
  path
}

#' Crop the HydroSHEDS flow directions.
#'
#' @param zip_file path, e.g. `download_hydrosheds()`'s output.
#' @param extent numeric, `c(xmin, xmax, ymin, ymax)` in longitude-latitude;
#'   must contain the whole catchments of the stations.
#' @param dir directory where to write the raster.
#' @return Path of `flow_direction_study_area.tif`.
#' @export
crop_flow_direction <- function(zip_file, extent, dir) {
  direction <- terra::rast(paste0("/vsizip/", zip_file, "/hyd_eu_dir_15s.tif"))
  path <- file.path(dir, "flow_direction_study_area.tif")
  terra::writeRaster(terra::crop(direction, terra::ext(extent)), path, overwrite = TRUE, datatype = "INT1U")
  path
}

#' Write the catchment mean of rasters at every cell.
#'
#' Flow accumulation sums a value over all the cells upstream of each cell
#' (the cell itself included). The catchment mean of a variable is the
#' accumulation of value x cell area divided by the accumulation of the cell
#' area, both over the cells where the variable is known. All layers must
#' have the same unknown cells (checked), so the known area is accumulated
#' once. The means are written to disk, so that the catchment mean at any
#' cell is then read directly (see `read_catchment_mean()`).
#'
#' @param direction SpatRaster of flow directions.
#' @param value SpatRaster on the grid of `direction`, one layer per
#'   variable.
#' @param path path of the GeoTIFF written.
#' @return `path`, with one layer per layer of `value` (same names) and a
#'   last layer `coverage`: the fraction of the catchment area where the
#'   layers are known.
#' @noRd
write_catchment_mean <- function(direction, value, path) {
  unknown <- terra::global(is.na(value), "sum")[, 1]
  stopifnot("all layers must have the same unknown cells" = length(unique(unknown)) == 1)
  area <- terra::cellSize(direction, unit = "km")
  known <- terra::flowAccumulation(direction, weight = terra::ifel(is.na(value[[1]]), 0, area))
  mean <- terra::rast(lapply(seq_len(terra::nlyr(value)), \(k) {
    terra::flowAccumulation(direction, weight = terra::ifel(is.na(value[[k]]), 0, value[[k]] * area)) / known
  }))
  coverage <- known / terra::flowAccumulation(direction, weight = area)
  out <- c(mean, coverage)
  names(out) <- c(names(value), "coverage")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(out, path, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
  path
}

#' Catchment mean of rasters at given cells, read from disk.
#'
#' @param path path, `write_catchment_mean()`'s output.
#' @param cell integer, cells of the flow-direction grid.
#' @return matrix, one row per cell, one column per variable, with the
#'   attribute `coverage`.
#' @noRd
read_catchment_mean <- function(path, cell) {
  value <- as.matrix(terra::rast(path)[cell])
  mean <- value[, colnames(value) != "coverage", drop = FALSE]
  attr(mean, "coverage") <- value[, "coverage"]
  mean
}

#' Snap the NAIADES stations to the flow-direction grid.
#'
#' The cell of a station is, among the cells within `radius`, the one whose
#' upstream area is closest to the drained area of the nearest RHT reach
#' (Réseau Hydrographique Théorique, computed on a finer network): the
#' coarse grid often puts the river one or two cells away from the station.
#'
#' @param station sf points with `code_station`, e.g. the `naiades_station`
#'   target.
#' @param rht sf lines with `surf_bv` (km²), e.g. `read_rht()`'s output.
#' @param direction_file path, e.g. `crop_flow_direction()`'s output.
#' @param radius numeric, search radius (m).
#' @return tibble with `code_station`, the `cell` of the station, its
#'   `longitude`, `latitude` and upstream `area` (km²), and `rht_area` (km²)
#'   and `rht_distance` (m) of the nearest RHT reach.
#' @export
snap_station_flow <- function(station, rht, direction_file, radius) {
  direction <- terra::rast(direction_file)
  upstream <- terra::flowAccumulation(direction, weight = terra::cellSize(direction, unit = "km"))
  lambert <- sf::st_transform(station["code_station"], 2154)
  reach <- sf::st_nearest_feature(lambert, rht)
  rht_area <- rht$surf_bv[reach]
  zone <- terra::vect(sf::st_transform(sf::st_buffer(lambert, radius), 4326))
  candidate <- terra::extract(upstream, zone, cells = TRUE, ID = TRUE)
  names(candidate) <- c("ID", "area", "cell")
  candidate <- candidate |>
    dplyr::filter(!is.na(area)) |>
    dplyr::slice_min(abs(log(area / rht_area[ID])), by = ID, with_ties = FALSE)
  coordinate <- terra::xyFromCell(direction, candidate$cell)
  tibble::tibble(
    code_station = station$code_station[candidate$ID],
    cell = candidate$cell,
    longitude = coordinate[, 1],
    latitude = coordinate[, 2],
    area = candidate$area,
    rht_area = rht_area[candidate$ID],
    rht_distance = as.numeric(sf::st_distance(lambert, rht[reach, ], by_element = TRUE))[candidate$ID]
  )
}

#' Download a vintage of CORINE Land Cover for mainland France.
#'
#' CORINE Land Cover (Copernicus, produced in France by the SDES and IGN)
#' maps land cover in 44 classes (3-level nomenclature) with a minimum
#' mapping unit of 25 ha, by photo-interpretation of satellite images.
#' Vector data in Lambert-93 from the SDES. The 2018 vintage is distributed
#' as a 7z archive inside a zip: it is extracted with 7-Zip (`7z` command)
#' and zipped again with the layout of the other vintages
#' (`CLC18_FR_RGF_SHP/CLC18_FR_RGF.*`). The download is skipped if the file
#' exists.
#'
#' @param dir directory where to store the zip file.
#' @param year integer, vintage: 2000, 2006, 2012 or 2018.
#' @return Path of `clc_<year>.zip`.
#' @export
download_corine <- function(dir, year) {
  millesime <- "https://www.donnees.statistiques.developpement-durable.gouv.fr/donneesCLC/CLC/millesime/"
  url <- c(
    "2000" = paste0(millesime, "CLC00_FR_RGF_SHP.zip"),
    "2006" = paste0(millesime, "CLC06_FR_RGF_SHP.zip"),
    "2012" = paste0(millesime, "CLC12_FR_RGF_SHP.zip"),
    "2018" = "https://www.statistiques.developpement-durable.gouv.fr/media/9687/download?inline"
  )[[as.character(year)]]
  path <- file.path(dir, sprintf("clc_%d.zip", year))
  if (file.exists(path)) {
    return(path)
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  if (year != 2018) {
    curl::curl_download(url, path)
    return(path)
  }
  temp <- tempfile("corine")
  shp_dir <- file.path(temp, "CLC18_FR_RGF_SHP")
  dir.create(shp_dir, recursive = TRUE)
  on.exit(unlink(temp, recursive = TRUE))
  curl::curl_download(url, file.path(temp, "download.zip"))
  archive <- utils::unzip(file.path(temp, "download.zip"), exdir = temp)
  system2("7z", c("e", "-y", shQuote(paste0("-o", shp_dir)), shQuote(archive), shQuote("CLC18_FR_RGF.*")), stdout = FALSE)
  system2("7z", c("a", "-tzip", shQuote(normalizePath(path, mustWork = FALSE)), shQuote(shp_dir)), stdout = FALSE)
  path
}

#' Read a vintage of CORINE Land Cover.
#'
#' @param zip_file path, e.g. `download_corine()`'s output.
#' @param frame optional sf object; only the polygons intersecting its
#'   bounding box are read.
#' @return sf tibble with the 3-digit CORINE `code` of each polygon,
#'   Lambert-93.
#' @export
read_corine <- function(zip_file, frame = NULL) {
  year <- as.integer(gsub("\\D", "", basename(zip_file)))
  layer <- sprintf("CLC%02d_FR_RGF", year %% 100)
  shp <- paste0("/vsizip/", zip_file, "/", sub("\\.zip$", "", layer), "_SHP/", layer, ".shp")
  wkt <- if (is.null(frame)) character(0) else {
    sf::st_as_text(sf::st_as_sfc(sf::st_bbox(sf::st_transform(frame, 2154))))
  }
  corine <- sf::read_sf(shp, wkt_filter = wkt)
  corine |>
    dplyr::transmute(code = as.integer(.data[[grep("^CODE_", names(corine), value = TRUE)[1]]])) |>
    sf::st_set_crs(2154)
}

#' Polygons of a vintage of CORINE Land Cover.
#'
#' @param zip_file path, e.g. `download_corine()`'s output.
#' @return tibble with the `year` (vintage), 3-digit `code` and `area` (ha)
#'   of each polygon.
#' @export
read_corine_polygon <- function(zip_file) {
  year <- as.integer(gsub("\\D", "", basename(zip_file)))
  layer <- sprintf("CLC%02d_FR_RGF", year %% 100)
  sf::st_read(
    paste0("/vsizip/", zip_file, "/", layer, "_SHP/", layer, ".dbf"),
    query = sprintf("SELECT CODE_%02d AS code, AREA_HA AS area FROM %s", year %% 100, layer),
    quiet = TRUE
  ) |>
    tibble::as_tibble() |>
    dplyr::transmute(year = year, code = as.integer(code), area)
}

#' Rasterise a vintage of CORINE Land Cover.
#'
#' The polygons are reprojected to longitude-latitude and rasterised with
#' GDAL, streaming from disk without loading them into R; each cell gets
#' the code of the polygon containing its centre (0 outside the map).
#'
#' @param zip_file path, e.g. `download_corine()`'s output.
#' @param extent numeric, `c(xmin, xmax, ymin, ymax)` in longitude-latitude.
#' @param resolution numeric, cell size (degrees).
#' @param path path of the raster to write.
#' @return `path`, a raster of 3-digit CORINE codes.
#' @export
rasterize_corine <- function(zip_file, extent, resolution, path) {
  year <- as.integer(gsub("\\D", "", basename(zip_file)))
  layer <- sprintf("CLC%02d_FR_RGF", year %% 100)
  field <- sprintf("CODE_%02d", year %% 100)
  polygon <- tempfile(fileext = ".gpkg")
  on.exit(unlink(polygon))
  sf::gdal_utils(
    "vectortranslate", paste0("/vsizip/", zip_file, "/", layer, "_SHP/", layer, ".shp"), polygon,
    options = c("-t_srs", "EPSG:4326", "-select", field)
  )
  sf::gdal_utils(
    "rasterize", polygon, path,
    options = c(
      "-a", field, "-te", extent[c(1, 3, 2, 4)], "-tr", resolution, resolution,
      "-ot", "UInt16", "-init", 0, "-co", "COMPRESS=DEFLATE"
    )
  )
  path
}

#' Catchment land cover at every cell of the flow-direction grid.
#'
#' The CORINE polygons are reprojected and rasterised with GDAL (streaming
#' from disk, without loading them into R) on a grid 5 times finer than the
#' flow-direction grid (3 arc-seconds, ~90 m). The fraction of each
#' flow-direction cell covered by each level-2 class (e.g. 21 arable land,
#' 23 pastures) is then the mean over its 5 x 5 fine cells; cells outside the
#' CORINE map (outside France, sea) are unknown. The share of each class in
#' the catchment of every cell is its catchment mean, written to disk (see
#' `write_catchment_mean()`).
#'
#' @param zip_file path, e.g. `download_corine()`'s output.
#' @param direction_file path, e.g. `crop_flow_direction()`'s output.
#' @param dir directory where to write the raster.
#' @return path of the raster `catchment_corine_<year>.tif`, one layer per
#'   level-2 class (share 0-1) and `coverage`.
#' @export
build_catchment_landcover <- function(zip_file, direction_file, dir) {
  option <- terra::terraOptions(print = FALSE)
  terra::terraOptions(memmax = 1, progress = 0) # process the fine grid by blocks
  on.exit(terra::terraOptions(memmax = option$memmax, progress = option$progress))
  direction <- terra::rast(direction_file)
  year <- as.integer(gsub("\\D", "", basename(zip_file)))
  fine <- tempfile(fileext = ".tif")
  on.exit(unlink(fine), add = TRUE)
  rasterize_corine(zip_file, as.vector(terra::ext(direction)), terra::res(direction)[1] / 5, fine)
  class <- c(11:14, 21:24, 31:33, 41:42, 51:52)
  level <- terra::rast(fine) %/% 10
  cover <- terra::rast(lapply(class, \(k) terra::aggregate(level == k, 5, mean)))
  mapped <- sum(cover)
  cover <- terra::ifel(mapped > 0, cover / mapped, NA)
  names(cover) <- class
  write_catchment_mean(direction, cover, file.path(dir, sprintf("catchment_corine_%d.tif", year)))
}

#' Land cover of the upstream catchment of the stations.
#'
#' @param catchment_file path, `build_catchment_landcover()`'s output.
#' @param station_flow tibble with `code_station` and `cell`, e.g.
#'   `snap_station_flow()`'s output.
#' @return tibble with `code_station`, `year` (CORINE vintage), `class`
#'   (level-2 code), `share` (% of the mapped catchment area) and `coverage`
#'   (% of the catchment area mapped by CORINE).
#' @export
compute_station_landcover <- function(catchment_file, station_flow) {
  share <- read_catchment_mean(catchment_file, station_flow$cell)
  tibble::as_tibble(share) |>
    dplyr::mutate(
      code_station = station_flow$code_station,
      year = as.integer(gsub("\\D", "", basename(catchment_file))),
      coverage = 100 * attr(share, "coverage")
    ) |>
    tidyr::pivot_longer(-c(code_station, year, coverage), names_to = "class", values_to = "share") |>
    dplyr::transmute(code_station, year, class = as.integer(class), share = 100 * share, coverage)
}

#' Crop the WorldClim monthly weather of 2000-2021.
#'
#' WorldClim 2.1 historical monthly weather: CRU-TS 4.06 monthly anomalies
#' (interpolated from weather stations, 0.5°) downscaled with the WorldClim
#' 2.1 climatology, at 2.5 arc-minutes (~4.5 km x 3 km in France). The
#' monthly rasters are read directly from the remote zip files, cropped and
#' stacked. Skipped if the file exists.
#'
#' @param dir directory where to write the raster.
#' @param variable character, `"prec"` (total precipitation, mm), `"tmin"`
#'   or `"tmax"` (mean daily minimum or maximum temperature, °C).
#' @param extent numeric, `c(xmin, xmax, ymin, ymax)` in longitude-latitude.
#' @return Path of the raster, one layer per month (named `YYYY-MM`).
#' @export
download_worldclim <- function(dir, variable, extent) {
  path <- file.path(dir, sprintf("worldclim_%s_2000-2021.tif", variable))
  if (file.exists(path)) {
    return(path)
  }
  url <- "/vsizip//vsicurl/https://geodata.ucdavis.edu/climate/worldclim/2_1/hist/cts4.06/2.5m"
  month <- format(seq(as.Date("2000-01-01"), as.Date("2021-12-01"), by = "month"), "%Y-%m")
  decade <- c("2000-2009", "2010-2019", "2020-2021")[findInterval(as.integer(substr(month, 1, 4)), c(2000, 2010, 2020))]
  weather <- terra::rast(lapply(seq_along(month), \(i) {
    terra::crop(terra::rast(sprintf(
      "%s/wc2.1_cruts4.06_2.5m_%s_%s.zip/wc2.1_2.5m_%s_%s.tif",
      url, variable, decade[i], variable, month[i]
    )), terra::ext(extent))
  }))
  names(weather) <- month
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(weather, path, datatype = "FLT4S")
  path
}

#' Yearly weather of the WorldClim cells.
#'
#' @param worldclim_file paths of the precipitation, minimum and maximum
#'   temperature rasters, e.g. `download_worldclim()`'s outputs.
#' @return SpatRaster with, for each year, `precipitation_<year>` (annual
#'   total, mm), `temperature_<year>` (annual mean of the daily mean, °C)
#'   and `summer_tmax_<year>` (June-August mean of the daily maximum, °C).
#' @export
summarise_worldclim <- function(worldclim_file) {
  read <- \(v) terra::rast(grep(paste0("_", v, "_"), worldclim_file, value = TRUE))
  precipitation <- read("prec")
  mean_temperature <- (read("tmin") + read("tmax")) / 2
  tmax <- read("tmax")
  year <- substr(names(precipitation), 1, 4)
  summer <- substr(names(precipitation), 6, 7) %in% c("06", "07", "08")
  yearly <- c(
    terra::tapp(precipitation, year, sum),
    terra::tapp(mean_temperature, year, mean),
    terra::tapp(tmax[[summer]], year[summer], mean)
  )
  names(yearly) <- paste(rep(c("precipitation", "temperature", "summer_tmax"), each = length(unique(year))), unique(year), sep = "_")
  yearly
}

#' Catchment weather at every cell of the flow-direction grid.
#'
#' The yearly WorldClim weather (see `summarise_worldclim()`) is split into
#' the flow-direction cells (10 x 10 per WorldClim cell), and the weather of
#' the catchment of every cell is its catchment mean, written to disk (see
#' `write_catchment_mean()`).
#'
#' @param worldclim_file paths, e.g. `download_worldclim()`'s outputs.
#' @param direction_file path, e.g. `crop_flow_direction()`'s output.
#' @param dir directory where to write the raster.
#' @return path of the raster `catchment_worldclim.tif`, with
#'   `precipitation_<year>`, `temperature_<year>`, `summer_tmax_<year>` and
#'   `coverage`.
#' @export
build_catchment_weather <- function(worldclim_file, direction_file, dir) {
  direction <- terra::rast(direction_file)
  weather <- terra::resample(summarise_worldclim(worldclim_file), direction, method = "near")
  write_catchment_mean(direction, weather, file.path(dir, "catchment_worldclim.tif"))
}

#' Yearly weather of the upstream catchment of the stations.
#'
#' @param catchment_file path, `build_catchment_weather()`'s output.
#' @param station_flow tibble with `code_station` and `cell`, e.g.
#'   `snap_station_flow()`'s output.
#' @return tibble with `code_station`, `year`, `precipitation`,
#'   `temperature` and `summer_tmax` (see `summarise_worldclim()`).
#' @export
compute_station_weather <- function(catchment_file, station_flow) {
  mean <- read_catchment_mean(catchment_file, station_flow$cell)
  tibble::as_tibble(mean) |>
    dplyr::mutate(code_station = station_flow$code_station) |>
    tidyr::pivot_longer(-code_station, names_to = c(".value", "year"), names_pattern = "(.*)_(\\d{4})") |>
    dplyr::mutate(year = as.integer(year))
}

#' Station-years of the driver models, with standardised drivers.
#'
#' Keeps the stations whose catchment matches the RHT drained area within a
#' factor 1.5. Land cover is the mean of the CORINE vintages (between
#' stations only). Climate is the mean weather of the catchment over
#' 2000-2021, and the weather of the year is its anomaly: log ratio of the
#' precipitation of the year (and of the previous year) to the mean, and
#' difference of the air temperature to the mean.
#'
#' @param station_year tibble, e.g. `compute_station_year_mean()`'s output.
#' @param water_quality tibble, e.g. `compute_water_quality()`'s output.
#' @param station_flow tibble, e.g. `snap_station_flow()`'s output.
#' @param landcover tibble, e.g. `compute_station_landcover()`'s outputs.
#' @param weather tibble, e.g. `compute_station_weather()`'s output.
#' @return tibble, one row per station-year: `code_station`, `year`, `X`,
#'   `Y`, the parameters and `PC1`, and the drivers scaled to mean 0 and SD
#'   1: `arable` (CORINE 21), `other_farmland` (22-24: permanent crops,
#'   pastures, heterogeneous agricultural areas), `artificial` (11-14) (shares
#'   of the catchment), `precipitation`, `air_temperature` (climate), `area` (log10
#'   catchment area), `precipitation_anomaly`, `precipitation_previous`,
#'   `temperature_anomaly`.
#' @export
prepare_driver_data <- function(station_year, water_quality, station_flow, landcover, weather) {
  station <- dplyr::filter(station_flow, abs(log(area / rht_area)) < log(1.5))
  cover <- landcover |>
    dplyr::mutate(type = dplyr::case_when(
      class == 21 ~ "arable", class %in% 22:24 ~ "other_farmland", class %in% 11:14 ~ "artificial"
    )) |>
    dplyr::filter(!is.na(type)) |>
    dplyr::summarise(share = sum(share), .by = c(code_station, year, type)) |>
    dplyr::summarise(share = mean(share), .by = c(code_station, type)) |>
    tidyr::pivot_wider(names_from = type, values_from = share)
  weather <- weather |>
    dplyr::arrange(code_station, year) |>
    dplyr::mutate(
      precipitation_anomaly = log(precipitation / mean(precipitation)),
      precipitation_previous = dplyr::lag(precipitation_anomaly),
      temperature_anomaly = temperature - mean(temperature),
      air_temperature = mean(temperature),
      precipitation = mean(precipitation),
      .by = code_station
    ) |>
    dplyr::select(code_station, year, precipitation, air_temperature, precipitation_anomaly, precipitation_previous, temperature_anomaly)
  driver <- c(
    "arable", "other_farmland", "artificial", "precipitation", "air_temperature", "area",
    "precipitation_anomaly", "precipitation_previous", "temperature_anomaly"
  )
  station_year |>
    dplyr::inner_join(water_quality, by = c("code_station", "year")) |>
    dplyr::inner_join(cover, by = "code_station") |>
    dplyr::inner_join(dplyr::transmute(station, code_station, area = log10(area)), by = "code_station") |>
    dplyr::inner_join(weather, by = c("code_station", "year")) |>
    tidyr::drop_na(dplyr::all_of(driver)) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(driver), \(x) as.numeric(scale(x))))
}

#' Fit the driver model of a parameter.
#'
#' Gaussian model of the standardised `response` on the drivers in `term`,
#' with a Matern spatial field (SPDE on `mesh`), an AR1 year effect and an
#' iid station effect, as the water quality model of D14.
#'
#' @param data tibble, e.g. `prepare_driver_data()`'s output.
#' @param mesh INLA mesh, e.g. `build_water_quality_mesh()`'s output.
#' @param response character, column of `data`.
#' @param term character, right-hand side of the driver effects.
#' @return list with `response`, `term`, `n` (station-years),
#'   `summary.fixed`, `summary.hyperpar` and `waic`.
#' @export
fit_driver_model <- function(data, mesh, response, term) {
  data$y <- as.numeric(scale(data[[response]]))
  spde <- INLA::inla.spde2.pcmatern(mesh, prior.range = c(50, 0.5), prior.sigma = c(1, 0.01))
  covariate <- dplyr::select(data, -dplyr::any_of(c("y", "X", "Y")))
  stack <- INLA::inla.stack(
    data = list(y = data$y),
    A = list(INLA::inla.spde.make.A(mesh, loc = as.matrix(data[c("X", "Y")])), 1),
    effects = list(
      list(field = seq_len(spde$n.spde), Intercept = 1),
      dplyr::rename(covariate, station = code_station)
    )
  )
  model <- INLA::inla(
    stats::as.formula(paste(
      "y ~ -1 + Intercept +", term,
      "+ f(field, model = spde) + f(year, model = 'ar1') + f(station, model = 'iid')"
    )),
    family = "gaussian", data = INLA::inla.stack.data(stack), num.threads = 4,
    quantiles = c(0.045, 0.5, 0.955),
    control.predictor = list(A = INLA::inla.stack.A(stack)),
    control.compute = list(waic = TRUE),
    control.inla = list(int.strategy = "eb")
  )
  out <- list(
    response = response,
    term = term,
    n = nrow(data),
    summary.fixed = model$summary.fixed,
    summary.hyperpar = model$summary.hyperpar,
    waic = model$waic$waic
  )
  rm(model)
  gc()
  out
}
