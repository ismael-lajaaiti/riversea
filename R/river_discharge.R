#' Query the Hub'Eau hydrometry API.
#'
#' @param endpoint character, e.g. "referentiel/stations" or "obs_elab".
#' @param query list of query parameters.
#' @return list, the parsed JSON response.
#' @noRd
query_hydrometry <- function(endpoint, query) {
  httr2::request("https://hubeau.eaufrance.fr/api/v2/hydrometrie") |>
    httr2::req_url_path_append(endpoint) |>
    httr2::req_url_query(!!!query) |>
    httr2::req_retry(max_tries = 5) |>
    httr2::req_throttle(rate = 5) |>
    httr2::req_perform() |>
    httr2::resp_body_json(simplifyVector = TRUE)
}

#' Gauging stations of the French hydrometric network in the area.
#'
#' Stations of the HYDRO database (Hub'Eau hydrometry API), open or closed.
#'
#' @param area sf polygons of the study area.
#' @return sf points (WGS84) with `code_station` and `libelle_station`.
#' @export
get_hydrometry_station <- function(area) {
  area <- sf::st_transform(area, 4326)
  query <- list(
    bbox = paste(sf::st_bbox(area), collapse = ","), format = "json", size = 10000,
    fields = "code_station,libelle_station,longitude_station,latitude_station"
  )
  first <- query_hydrometry("referentiel/stations", c(query, page = 1))
  pages <- seq_len(ceiling(first$count / 10000))[-1]
  dplyr::bind_rows(
    first$data,
    purrr::map(pages, \(page) query_hydrometry("referentiel/stations", c(query, page = page))$data)
  ) |>
    dplyr::filter(!is.na(longitude_station), !is.na(latitude_station)) |>
    sf::st_as_sf(coords = c("longitude_station", "latitude_station"), crs = 4326) |>
    sf::st_filter(area)
}

#' Monthly mean discharge of gauging stations.
#'
#' Validated or not, as served by Hub'Eau (`obs_elab`, QmM).
#'
#' @param station character, station codes (up to ~100 per call: the API
#'   serves 20,000 results per page).
#' @param period character, first and last dates ("YYYY-MM-DD").
#' @return tibble with `code_station`, `date` (first day of the month) and
#'   `discharge` (m³/s).
#' @export
download_monthly_discharge <- function(station, period) {
  response <- query_hydrometry("obs_elab", list(
    code_entite = paste(station, collapse = ","), grandeur_hydro_elab = "QmM",
    date_debut_obs_elab = period[1], date_fin_obs_elab = period[2], size = 20000,
    fields = "code_station,date_obs_elab,resultat_obs_elab"
  ))
  stopifnot(response$count <= 20000)
  if (response$count == 0) {
    return(tibble::tibble(code_station = character(0), date = as.Date(character(0)), discharge = numeric(0)))
  }
  tibble::as_tibble(response$data) |>
    dplyr::transmute(code_station, date = as.Date(date_obs_elab), discharge = resultat_obs_elab / 1000)
}

#' Yearly discharge entering each estuary.
#'
#' The gauging stations of each river mouth of an estuary with all 12
#' monthly means (missing values dropped) in at least `years_min` years of
#' `period` are eligible. Each year, a mouth takes the runoff (discharge
#' per km²) of its most downstream eligible station measured that year
#' (largest drained area on the flow-direction grid); the discharge of the
#' estuary is the mean runoff of its mouths (weighted by their catchment
#' area) times the catchment area of the estuary: a mouth missing a year
#' takes the runoff of the other mouths.
#'
#' @param discharge tibble, [download_monthly_discharge()]'s output.
#' @param station_mouth tibble with `code_station`, `mouth_id` and `area`
#'   (km² drained), [link_station_mouth()]'s output for the gauging
#'   stations.
#' @param mouth tibble with `estuary` and `mouth_id`.
#' @param mouth_size tibble with `mouth_id` and `area`.
#' @param period integer, first and last years.
#' @param years_min integer.
#' @return list: `station` (tibble `estuary`, `mouth_id`, `code_station`,
#'   `area`, `share` of the mouth catchment drained, `years`: years used
#'   for the mouth), `station_year` (tibble `code_station`, `year`,
#'   `discharge`, m³/s, and `runoff`, l/s/km², for all the gauging stations
#'   of the estuaries), and `estuary_year` (tibble `estuary`, `year`,
#'   `discharge`, m³/s).
#' @export
summarise_estuary_discharge <- function(discharge, station_mouth, mouth, mouth_size, period, years_min) {
  station_year <- discharge |>
    dplyr::filter(!is.na(discharge), discharge >= 0) |>
    dplyr::mutate(year = as.integer(format(date, "%Y"))) |>
    dplyr::filter(dplyr::between(year, period[1], period[2])) |>
    dplyr::summarise(months = dplyr::n_distinct(date), discharge = mean(discharge), .by = c(code_station, year)) |>
    dplyr::filter(months == 12) |>
    dplyr::select(-months) |>
    dplyr::inner_join(dplyr::select(station_mouth, code_station, area), by = "code_station") |>
    dplyr::mutate(runoff = 1000 * discharge / area) |>
    dplyr::select(-area)
  mouth_year <- mouth |>
    dplyr::left_join(dplyr::rename(mouth_size, catchment = area), by = "mouth_id") |>
    dplyr::inner_join(station_mouth, by = "mouth_id") |>
    dplyr::inner_join(dplyr::count(station_year, code_station, name = "years"), by = "code_station") |>
    dplyr::filter(years >= years_min) |>
    dplyr::inner_join(station_year, by = "code_station") |>
    dplyr::slice_max(area, by = c(estuary, mouth_id, year), with_ties = FALSE)
  catchment <- mouth |>
    dplyr::left_join(dplyr::rename(mouth_size, catchment = area), by = "mouth_id") |>
    dplyr::summarise(catchment = sum(catchment), .by = estuary)
  estuary_year <- mouth_year |>
    dplyr::summarise(runoff = stats::weighted.mean(runoff, catchment), .by = c(estuary, year)) |>
    dplyr::inner_join(catchment, by = "estuary") |>
    dplyr::transmute(estuary, year, discharge = runoff * catchment / 1000)
  list(
    station = mouth_year |>
      dplyr::summarise(area = area[1], share = area[1] / catchment[1], years = dplyr::n(), .by = c(estuary, mouth_id, code_station)),
    station_year = dplyr::filter(station_year, code_station %in% station_mouth$code_station[station_mouth$mouth_id %in% mouth$mouth_id]),
    estuary_year = estuary_year
  )
}

#' Stations downstream of each station along the flow directions.
#'
#' @param cell integer, cells of the stations on the flow-direction grid.
#' @param direction_file path, e.g. `crop_flow_direction()`'s output.
#' @return tibble with `from` and `to`, indices in `cell`: station `to` is
#'   downstream of station `from` (another cell on its path to the sea).
#' @noRd
find_downstream_station <- function(cell, direction_file) {
  down <- flow_downstream(terra::rast(direction_file))
  current <- cell
  owner <- seq_along(cell)
  pair <- list()
  repeat {
    current <- down[current]
    keep <- !is.na(current)
    current <- current[keep]
    owner <- owner[keep]
    if (length(current) == 0) break
    hit <- match(current, cell)
    if (any(!is.na(hit))) {
      # Every station sharing the hit cell is downstream.
      to <- purrr::map(current[!is.na(hit)], \(x) which(cell == x))
      pair[[length(pair) + 1]] <- tibble::tibble(from = rep(owner[!is.na(hit)], lengths(to)), to = unlist(to))
    }
  }
  dplyr::distinct(dplyr::bind_rows(pair, tibble::tibble(from = integer(0), to = integer(0))))
}

#' Yearly river inputs of each estuary, mixing the stations of its branches.
#'
#' For each estuary, the eligible NAIADES stations (outside the water
#' bodies, enough years: `status` "chosen" or "upstream" in `candidate`)
#' are taken from the largest drained area down, skipping the stations
#' upstream of one already taken: the stations kept are not nested, each
#' samples a different branch reaching the estuary. Their yearly values are
#' mixed in proportion to their discharge: drained area times the yearly
#' runoff (discharge per km²) of the gauging station on the same branch
#' (nested with the station, closest drained area); without one, the runoff
#' of the gauging stations of the estuary that year, else the mean runoff of
#' the other stations of the estuary that year; without any, drained area
#' alone. Nutrients are mixed as concentrations (then log10), oxygen
#' saturation as is. A year is kept when the stations measured that year
#' drain at least `cover` of the drained area of all the stations kept.
#'
#' @param candidate tibble, `station` of [summarise_estuary_nutrient()].
#' @param station_year tibble with `code_station`, `year`, the log10 yearly
#'   means `ammonium`, `nitrate`, `orthophosphate` and `oxygen_saturation`.
#' @param station_flow tibble with `code_station` and `cell` of the NAIADES
#'   stations, [snap_station_flow()]'s output.
#' @param gauge tibble with `code_station`, `cell`, `mouth_id` and `area` of
#'   the gauging stations.
#' @param discharge list, [summarise_estuary_discharge()]'s output.
#' @param direction_file path, e.g. `crop_flow_direction()`'s output.
#' @param period integer, first and last years.
#' @param cover numeric, share of the drained area.
#' @return list: `station` (tibble `estuary`, `mouth_id`, `code_station`,
#'   `area`, `share` of the estuary catchment (all its mouths), `gauge`: gauging station
#'   giving its runoff, NA if none), and `estuary_year` (tibble `estuary`,
#'   `year`, `ammonium`, `nitrate`, `orthophosphate` (log10 mg/l),
#'   `river_oxygen` (%), `cover`: share of the estuary catchment drained
#'   by the stations measured that year).
#' @export
mix_estuary_nutrient <- function(candidate, station_year, station_flow, gauge, discharge, direction_file, period, cover) {
  eligible <- candidate |>
    dplyr::filter(status %in% c("chosen", "upstream")) |>
    dplyr::mutate(catchment = area / share) |>
    dplyr::inner_join(dplyr::select(station_flow, code_station, cell), by = "code_station")
  gauge <- dplyr::filter(gauge, mouth_id %in% eligible$mouth_id)
  node <- dplyr::bind_rows(
    dplyr::transmute(eligible, kind = "station", code_station, mouth_id, cell, area),
    dplyr::transmute(gauge, kind = "gauge", code_station, mouth_id, cell, area)
  ) |>
    dplyr::distinct(kind, code_station, mouth_id, .keep_all = TRUE)
  pair <- find_downstream_station(node$cell, direction_file) |>
    dplyr::filter(node$mouth_id[from] == node$mouth_id[to])
  # Stations kept: largest first, none downstream of another kept one.
  kept <- purrr::map_dfr(split(eligible, eligible$estuary), \(e) {
    e <- dplyr::arrange(e, dplyr::desc(area))
    index <- match(paste("station", e$code_station, e$mouth_id), paste(node$kind, node$code_station, node$mouth_id))
    taken <- integer(0)
    for (i in index) if (!any(pair$to[pair$from == i] %in% taken)) taken <- c(taken, i)
    e[index %in% taken, ]
  })
  # Gauging station on the same branch: nested with the station, closest area.
  gauge_year <- discharge$station_year
  kept$gauge <- purrr::map_chr(seq_len(nrow(kept)), \(k) {
    i <- match(paste("station", kept$code_station[k], kept$mouth_id[k]), paste(node$kind, node$code_station, node$mouth_id))
    nested <- unique(c(pair$to[pair$from == i], pair$from[pair$to == i]))
    nested <- nested[node$kind[nested] == "gauge" & node$code_station[nested] %in% gauge_year$code_station]
    if (length(nested) == 0) return(NA_character_)
    node$code_station[nested[which.min(abs(log(node$area[nested] / kept$area[k])))]]
  })
  estuary_runoff <- discharge$station |>
    dplyr::inner_join(gauge_year, by = "code_station") |>
    dplyr::summarise(estuary_runoff = stats::weighted.mean(runoff, area), .by = c(estuary, year))
  total <- kept |>
    dplyr::summarise(total = sum(area), catchment = catchment[1], .by = c(estuary, mouth_id)) |>
    dplyr::summarise(total = sum(total), catchment = sum(catchment), .by = estuary)
  estuary_year <- kept |>
    dplyr::inner_join(dplyr::filter(station_year, dplyr::between(year, period[1], period[2])), by = "code_station") |>
    dplyr::left_join(dplyr::select(gauge_year, gauge = code_station, year, runoff), by = c("gauge", "year")) |>
    dplyr::left_join(estuary_runoff, by = c("estuary", "year")) |>
    dplyr::mutate(runoff = dplyr::coalesce(runoff, estuary_runoff)) |>
    dplyr::mutate(runoff = dplyr::if_else(is.na(runoff), mean(runoff, na.rm = TRUE), runoff), .by = c(estuary, year)) |>
    dplyr::mutate(runoff = dplyr::if_else(is.na(runoff), 1, runoff), weight = area * runoff) |>
    dplyr::summarise(
      measured = sum(area),
      dplyr::across(c(ammonium, nitrate, orthophosphate), \(x) log10(sum(weight * 10^x) / sum(weight))),
      river_oxygen = sum(weight * oxygen_saturation) / sum(weight),
      .by = c(estuary, year)
    ) |>
    dplyr::inner_join(total, by = "estuary") |>
    dplyr::filter(measured >= cover * total) |>
    dplyr::mutate(cover = measured / catchment) |>
    dplyr::select(estuary, year, ammonium, nitrate, orthophosphate, river_oxygen, cover)
  list(
    station = kept |>
      dplyr::left_join(dplyr::select(total, estuary, estuary_catchment = catchment), by = "estuary") |>
      dplyr::mutate(share = area / estuary_catchment) |>
      dplyr::select(estuary, mouth_id, code_station, area, share, gauge),
    estuary_year = estuary_year
  )
}
