#' Download the theoretical hydrographic network (RHT).
#'
#' The RHT has no official public download. A copy is stored, uncompressed,
#' inside the MiSTE workshop archive on Zenodo (10.5281/zenodo.17962542,
#' CC-BY 4.0). Only its bytes are downloaded (~52 MB), with HTTP range
#' requests, instead of the whole archive (~874 MB): the archive's central
#' directory gives the position of the RHT file.
#'
#' @param dir directory where to store the zip file.
#' @return Path of `reseau_hydro_rht_2020.zip`.
#' @export
download_rht <- function(dir) {
  url <- "https://zenodo.org/api/records/17962542/files/miste_data.zip/content"
  file_name <- "reseau_hydro_rht_2020.zip"
  range_handle <- function(from, to) {
    curl::new_handle(range = sprintf("%.0f-%.0f", from, to))
  }
  fetch <- function(from, to) {
    curl::curl_fetch_memory(url, range_handle(from, to))
  }
  # Little-endian unsigned integer of `n` bytes, at 0-based offset `at`.
  uint <- function(bytes, at, n) {
    sum(as.numeric(bytes[at + seq_len(n)]) * 256^(seq_len(n) - 1))
  }

  headers <- curl::parse_headers(fetch(0, 0)$headers)
  content_range <- grep("^content-range", headers, ignore.case = TRUE, value = TRUE)
  size <- as.numeric(sub(".*/", "", content_range))

  # End of central directory record, within the archive's last 64 KB.
  tail_start <- max(0, size - 65536)
  tail <- fetch(tail_start, size - 1)$content
  eocd <- max(grepRaw(as.raw(c(0x50, 0x4b, 0x05, 0x06)), tail, all = TRUE)) - 1
  cd_size <- uint(tail, eocd + 12, 4)
  cd_offset <- uint(tail, eocd + 16, 4)
  stopifnot("Zip64 archives are not supported" = cd_offset < 2^32 - 1)
  cd <- fetch(cd_offset, cd_offset + cd_size - 1)$content

  # Walk the central directory entries until the RHT file.
  pos <- 0
  repeat {
    stopifnot("RHT file not found in the archive" = pos < length(cd))
    name_len <- uint(cd, pos + 28, 2)
    name <- rawToChar(cd[pos + 46 + seq_len(name_len)])
    if (basename(name) == file_name) break
    pos <- pos + 46 + name_len + uint(cd, pos + 30, 2) + uint(cd, pos + 32, 2)
  }
  stopifnot("RHT file is compressed in the archive" = uint(cd, pos + 10, 2) == 0)
  data_size <- uint(cd, pos + 20, 4)
  header_offset <- uint(cd, pos + 42, 4)

  # The file data starts after its local header.
  header <- fetch(header_offset, header_offset + 29)$content
  data_offset <- header_offset + 30 + uint(header, 26, 2) + uint(header, 28, 2)

  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  path <- file.path(dir, file_name)
  curl::curl_download(
    url, path,
    handle = range_handle(data_offset, data_offset + data_size - 1)
  )
  path
}

#' Read the theoretical hydrographic network (RHT) with its attributes.
#'
#' The RHT (Pella et al., 2012) is a river network for mainland France
#' derived from the IGN digital elevation model. Reach geometries
#' (`rht_lbt93.shp`) are joined to their environmental attributes
#' (`Attr_RHT_Avril2020.csv`) by `id_drain`. A Strahler order of 0 marks
#' the few reaches where the RHT models were not computed (discharge at its
#' 0.001 floor, no width or depth), so it is recoded to `NA`.
#'
#' @param zip_file path of `reseau_hydro_rht_2020.zip`.
#'
#' @return sf tibble, one row per reach, CRS Lambert-93.
#' @export
read_rht <- function(zip_file) {
  dir <- tempfile("rht")
  utils::unzip(zip_file, exdir = dir)
  attributes <- readr::read_delim(
    file.path(dir, "Attr_RHT_Avril2020.csv"),
    delim = ";", show_col_types = FALSE
  )
  sf::read_sf(file.path(dir, "rht_lbt93.shp")) |>
    dplyr::select(id_drain = ID_DRAIN) |>
    dplyr::inner_join(attributes, by = "id_drain") |>
    dplyr::mutate(strahler = dplyr::na_if(strahler, 0))
}

#' Flag operations correctly matched to their RHT reach.
#'
#' Snapping distance alone misses wrong matches: stations on small streams
#' absent from the RHT are snapped to a nearby larger river, even at short
#' distance. A match is kept if the RHT upstream catchment area (`surf_bv`)
#' and the AMOBIO drained area of the station agree within a factor
#' `area_ratio_max`, or, when the AMOBIO area is unknown, if the snapping
#' distance is at most `snap_dist_max` (see D11).
#'
#' @param operation tibble, e.g. `snap_operations_rht()`'s output.
#' @param river_operation tibble with `operation_id`, `sandre_code`.
#' @param amobio_metrics tibble with `sandre_code`, `drained_area`, e.g.
#'   `dedup_amobio_metrics()`'s output.
#' @param area_ratio_max maximum ratio between the two catchment areas.
#' @param snap_dist_max maximum snapping distance (m).
#'
#' @return `operation` with `drained_area`, `area_ratio` and the logical
#' `rht_match` added.
#' @export
check_rht_match <- function(operation,
                            river_operation,
                            amobio_metrics,
                            area_ratio_max,
                            snap_dist_max) {
  drained_area <- amobio_metrics |>
    dplyr::filter(!is.na(drained_area)) |>
    dplyr::distinct(sandre_code, .keep_all = TRUE) |>
    dplyr::select(sandre_code, drained_area)
  station <- river_operation |>
    dplyr::distinct(operation_id, sandre_code)

  operation |>
    dplyr::left_join(station, by = "operation_id") |>
    dplyr::left_join(drained_area, by = "sandre_code") |>
    dplyr::mutate(
      area_ratio = surf_bv / drained_area,
      rht_match = dplyr::if_else(
        is.na(area_ratio),
        snap_dist_m <= snap_dist_max,
        area_ratio <= area_ratio_max & area_ratio >= 1 / area_ratio_max
      )
    ) |>
    dplyr::select(-sandre_code)
}

#' Snap operations to the nearest RHT reach.
#'
#' @param operation tibble with `operation_id`, `longitude`, `latitude`,
#'   `district`, e.g. `classify_operation_location()`'s output.
#' @param rht sf tibble, e.g. `read_rht()`'s output.
#'
#' @return tibble with one row per operation, the attributes of its nearest
#' reach, and `snap_dist_m` the straight-line distance to that reach.
#' @export
snap_operations_rht <- function(operation, rht) {
  op_sf <- operation |>
    dplyr::distinct(operation_id, longitude, latitude, district) |>
    sf::st_as_sf(
      coords = c("longitude", "latitude"), crs = 4326, remove = FALSE
    ) |>
    sf::st_transform(sf::st_crs(rht))

  nearest <- sf::st_nearest_feature(op_sf, rht)
  op_sf |>
    dplyr::mutate(
      snap_dist_m = as.numeric(
        sf::st_distance(op_sf, rht[nearest, ], by_element = TRUE)
      )
    ) |>
    sf::st_drop_geometry() |>
    dplyr::bind_cols(sf::st_drop_geometry(rht)[nearest, ])
}
