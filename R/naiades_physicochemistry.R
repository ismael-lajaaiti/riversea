#' Query the Hub'Eau river quality API.
#'
#' Hub'Eau serves the NAIADES database. Transient errors are retried.
#'
#' @param endpoint character, e.g. "analyse_pc" or "station_pc".
#' @param query named list of query parameters.
#' @return parsed JSON response, with `count` and `data`.
#' @export
query_hubeau <- function(endpoint, query) {
  httr2::request("https://hubeau.eaufrance.fr/api/v2/qualite_rivieres") |>
    httr2::req_url_path_append(endpoint) |>
    httr2::req_url_query(!!!query) |>
    httr2::req_retry(max_tries = 5) |>
    httr2::req_throttle(rate = 5) |>
    httr2::req_perform() |>
    httr2::resp_body_json(simplifyVector = TRUE)
}

#' NAIADES physico-chemistry stations measuring the parameters in the area.
#'
#' @param area sf polygons of the study area.
#' @param parameter character, Sandre codes of the parameters.
#' @param period character, first and last sampling dates ("YYYY-MM-DD").
#' @return sf points (WGS84) with `code_station`, `libelle_station` and
#'   `type_entite_hydro`, for the stations within `area` with at least one
#'   measurement of the parameters in the period.
#' @export
get_naiades_station <- function(area, parameter, period) {
  area <- sf::st_transform(area, 4326)
  query <- list(
    bbox = paste(sf::st_bbox(area), collapse = ","),
    code_parametre = paste(parameter, collapse = ","),
    date_debut_prelevement = period[1],
    date_fin_prelevement = period[2],
    fields = "code_station,libelle_station,type_entite_hydro,longitude,latitude",
    exact_count = "true",
    size = 5000
  )
  first <- query_hubeau("station_pc", c(query, page = 1))
  # The API does not serve results beyond the 20,000th.
  stopifnot(first$count <= 20000)
  pages <- seq_len(ceiling(first$count / 5000))[-1]
  station <- dplyr::bind_rows(
    first$data,
    purrr::map(pages, \(page) query_hubeau("station_pc", c(query, page = page))$data)
  )
  station |>
    dplyr::filter(!is.na(longitude), !is.na(latitude)) |>
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
    sf::st_filter(area)
}

#' NAIADES physico-chemistry measurements of stations.
#'
#' The API does not serve results beyond the 20,000th of a query, so the
#' query is split by station, then by period, until each part fits.
#'
#' @param station character, station codes.
#' @param parameter character, Sandre codes of the parameters.
#' @param period character, first and last sampling dates ("YYYY-MM-DD").
#' @return tibble, one row per measurement, with the station, the date and
#'   time of sampling, the parameter, the fraction analysed, the result and
#'   its unit, the remark code (e.g. 10 for a result below the
#'   quantification limit), the quantification limit, and the
#'   qualification and status codes of the result.
#' @export
download_naiades_pc <- function(station, parameter, period) {
  field <- c(
    "code_station", "date_prelevement", "heure_prelevement", "code_parametre",
    "libelle_parametre", "code_fraction", "resultat", "symbole_unite",
    "code_remarque", "limite_quantification", "code_qualification", "code_statut"
  )
  query <- list(
    code_station = paste(station, collapse = ","),
    code_parametre = paste(parameter, collapse = ","),
    date_debut_prelevement = period[1],
    date_fin_prelevement = period[2],
    fields = paste(field, collapse = ",")
  )
  count <- query_hubeau("analyse_pc", c(query, size = 1))$count
  if (count == 0) {
    return(NULL)
  }
  if (count > 20000) {
    if (length(station) > 1) {
      half <- seq_along(station) <= length(station) / 2
      return(dplyr::bind_rows(
        download_naiades_pc(station[half], parameter, period),
        download_naiades_pc(station[!half], parameter, period)
      ))
    }
    middle <- as.Date(period[1]) + floor(diff(as.Date(period)) / 2)
    return(dplyr::bind_rows(
      download_naiades_pc(station, parameter, c(period[1], format(middle))),
      download_naiades_pc(station, parameter, c(format(middle + 1), period[2]))
    ))
  }
  # Pages are read in CSV, several times faster to serve than JSON.
  purrr::map(seq_len(ceiling(count / 10000)), \(page) {
    httr2::request("https://hubeau.eaufrance.fr/api/v2/qualite_rivieres/analyse_pc.csv") |>
      httr2::req_url_query(!!!c(query, size = 10000, page = page)) |>
      httr2::req_retry(max_tries = 5) |>
      httr2::req_throttle(rate = 5) |>
      httr2::req_perform() |>
      httr2::resp_body_string() |>
      I() |>
      readr::read_delim(delim = ";", col_types = readr::cols(.default = "c"))
  }) |>
    dplyr::bind_rows() |>
    dplyr::transmute(
      code_station = as.character(code_station),
      date = as.Date(date_prelevement),
      time = as.character(heure_prelevement),
      code_parametre = as.character(code_parametre),
      libelle_parametre = as.character(libelle_parametre),
      code_fraction = as.character(code_fraction),
      result = as.numeric(resultat),
      unit = as.character(symbole_unite),
      code_remarque = as.character(code_remarque),
      quantification_limit = as.numeric(limite_quantification),
      code_qualification = as.character(code_qualification),
      code_statut = as.character(code_statut)
    )
}

#' Clean the NAIADES physico-chemistry measurements.
#'
#' Keeps the measurements on water (not on sediments), drops those
#' qualified as incorrect or uncertain and those without a result, puts the
#' nutrients expressed as N or P back in the unit of the ion (the generic
#' "mg/L" label is the same unit as the explicit one), sets the results
#' below the detection or quantification limit to half the limit (the
#' result field then holds the limit), and removes outliers: values outside
#' physical bounds for temperature and oxygen, and, for nutrients, values more
#' than 5 median absolute deviations away from the median of the station on
#' the log scale (stations with at least 10 measurements).
#'
#' @param pc tibble, `download_naiades_pc()`'s output.
#' @return `pc` cleaned, with `value` (the result used) and `below_limit`.
#' @export
clean_naiades_pc <- function(pc) {
  unit <- tibble::tribble(
    ~code_parametre, ~unit, ~factor,
    "1335", "mg(NH4)/L", 1, "1335", "mg/L", 1, "1335", "mg(N)/L", 18.04 / 14.01,
    "1339", "mg(NO2)/L", 1, "1339", "mg/L", 1, "1339", "mg(N)/L", 46.01 / 14.01,
    "1340", "mg(NO3)/L", 1, "1340", "mg/L", 1, "1340", "mg(N)/L", 62.00 / 14.01,
    "1340", "mg(CO3)/L", 1, # mislabelled nitrate, values in mg(NO3)/L
    "1433", "mg(PO4)/L", 1, "1433", "mg/L", 1, "1433", "mg(P)/L", 94.97 / 30.97,
    "1350", "mg(P)/L", 1, "1350", "mg/L", 1,
    "1319", "mg(N)/L", 1, "1319", "mg/L", 1,
    "1301", "°C", 1,
    "1311", "mg(O2)/L", 1, "1311", "mg/L", 1,
    "1312", "%", 1
  )
  bound <- tibble::tribble(
    ~code_parametre, ~lower, ~upper,
    "1301", 0, 35, "1311", 0, 25, "1312", 0, 250
  )
  pc |>
    dplyr::filter(
      code_fraction %in% c("3", "22", "23"), # filtered, unknown and raw water
      !code_qualification %in% c("2", "3"), # incorrect, uncertain
      code_remarque %in% c("1", "2", "7", "10"),
      code_remarque == "1" | !code_parametre %in% c("1301", "1311", "1312"), # limit flags are errors
      !is.na(result)
    ) |>
    dplyr::inner_join(unit, by = c("code_parametre", "unit")) |>
    dplyr::mutate(
      below_limit = code_remarque != "1",
      value = result * factor / dplyr::if_else(below_limit, 2, 1)
    ) |>
    dplyr::filter(value > 0) |>
    dplyr::left_join(bound, by = "code_parametre") |>
    dplyr::filter(is.na(lower) | (value >= lower & value <= upper)) |>
    dplyr::mutate(
      log_value = log(value),
      deviation = abs(log_value - stats::median(log_value)) / stats::mad(log_value),
      outlier = is.na(lower) & dplyr::n() >= 10 & !below_limit & is.finite(deviation) & deviation > 5,
      .by = c(code_station, code_parametre)
    ) |>
    dplyr::filter(!outlier) |>
    dplyr::select(-factor, -lower, -upper, -log_value, -deviation, -outlier)
}
