nutrient_driver_targets <- list(
  tar_target(
    study_extent, # Longitude-latitude box (xmin, xmax, ymin, ymax) containing the catchments of the stations.
    c(-5.5, 7.5, 42, 51.5)
  ),
  tar_target(
    hydrosheds_zip, # HydroSHEDS v1 flow directions of Europe (15 arc-seconds).
    download_hydrosheds(here("data", "geo", "hydrosheds")),
    format = "file"
  ),
  tar_target(
    flow_direction_file, # Flow directions cropped to the study extent.
    crop_flow_direction(hydrosheds_zip, study_extent, here("data", "geo", "hydrosheds")),
    format = "file"
  ),
  tar_target(
    station_flow, # Cell of the flow-direction grid of each NAIADES station, matching the RHT drained area within 1 km.
    snap_station_flow(naiades_station, rht, flow_direction_file, 1000)
  ),
  tar_target(
    corine_year, # CORINE Land Cover vintages (2018 needs 7-Zip).
    c(2000, 2006, 2012, 2018)
  ),
  tar_target(
    corine_zip,
    download_corine(here("data", "geo", "corine"), corine_year),
    pattern = map(corine_year),
    format = "file"
  ),
  tar_target(
    corine_polygon, # Code and area of the CORINE polygons, per vintage.
    read_corine_polygon(corine_zip),
    pattern = map(corine_zip)
  ),
  tar_target(
    corine_map_file, # CORINE 2018 rasterised at 30 arc-seconds (~0.9 km x 0.6 km), for maps.
    rasterize_corine(
      grep("2018", corine_zip, value = TRUE), study_extent, 1 / 120,
      here("data", "geo", "corine", "clc_2018_map.tif")
    ),
    format = "file"
  ),
  tar_target(
    station_landcover, # Share of each CORINE level-2 class in the upstream catchment of the stations.
    compute_station_landcover(corine_zip, flow_direction_file, station_flow),
    pattern = map(corine_zip),
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    worldclim_variable, # WorldClim monthly weather variables: precipitation, minimum and maximum temperature.
    c("prec", "tmin", "tmax")
  ),
  tar_target(
    worldclim_file, # Monthly weather of 2000-2021 cropped to the study extent.
    download_worldclim(here("data", "geo", "worldclim"), worldclim_variable, study_extent),
    pattern = map(worldclim_variable),
    format = "file"
  ),
  tar_target(
    station_weather, # Yearly weather of the upstream catchment of the stations.
    compute_station_weather(worldclim_file, flow_direction_file, station_flow),
    memory = "transient",
    garbage_collection = TRUE
  ),
  tar_target(
    driver_data, # Station-years of the stations with fish visits, with standardised drivers (land cover, climate, weather anomalies, catchment area).
    prepare_driver_data(
      naiades_station_year, water_quality,
      dplyr::filter(station_flow, code_station %in% water_quality_visit$code_station), # stations with fish visits
      station_landcover, station_weather
    )
  ),
  tar_target(
    driver_setting, # One driver model per parameter.
    tibble::tibble(
      response = c("nitrate", "ammonium", "orthophosphate", "total_phosphorus", "PC1"),
      term = paste(
        "arable + other_farmland + artificial + precipitation + air_temperature",
        "+ area + precipitation_anomaly + precipitation_previous + temperature_anomaly"
      )
    )
  ),
  tar_target(
    driver_model, # Driver models: drivers + spatial field + AR1 year + iid station.
    fit_driver_model(driver_data, water_quality_mesh, driver_setting$response, driver_setting$term),
    pattern = map(driver_setting),
    iteration = "list",
    memory = "transient",
    garbage_collection = TRUE
  )
)
