# =============================================================================
# R/embedding_tl/gee_embeddings.R
# Extract GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL at global core locations
# and download 64-band embedding raster over the local AOI.
# =============================================================================

.EMB_COLLECTION <- "GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL"
.EMB_N_BANDS    <- 64L
.EMB_SCALE      <- 10L   # native resolution (metres)


# Internal: build mean embedding image averaged over a range of years.
.build_emb_mean <- function(years) {
  suppressPackageStartupMessages(library(rgee))
  start <- paste0(min(years), "-01-01")
  end   <- paste0(max(years) + 1L, "-01-01")
  ee$ImageCollection(.EMB_COLLECTION)$
    filterDate(start, end)$
    mean()
}


# Internal: rename bands to emb_1 … emb_N using an ee$Image expression list.
# GEE band names for this collection are unnamed by default; we rename in R
# after download rather than in GEE to avoid EE string-expression complexity.
.rename_emb_bands <- function(df, n = .EMB_N_BANDS) {
  emb_cols_old <- grep("^(emb_|embedding_|b)", names(df), value = TRUE,
                       ignore.case = TRUE)
  emb_cols_old <- emb_cols_old[seq_len(min(n, length(emb_cols_old)))]
  if (length(emb_cols_old) == n) {
    names(df)[match(emb_cols_old, names(df))] <- paste0("emb_", seq_len(n))
  }
  df
}


# -----------------------------------------------------------------------------
# extract_global_embeddings()
# -----------------------------------------------------------------------------
# For each unique profile in profiles_df (or a path to the covariates CSV),
# extract the mean 64-d embedding vector averaged over `years`.
#
# Returns a data.frame: profile_id, dataset, emb_1 … emb_64.
# -----------------------------------------------------------------------------
extract_global_embeddings <- function(profiles_df, gee_project = NULL,
                                      years = 2023:2025) {
  suppressPackageStartupMessages({
    library(rgee); library(dplyr); library(readr)
  })

  if (is.character(profiles_df) && length(profiles_df) == 1L &&
      file.exists(profiles_df)) {
    profiles_df <- read_csv(profiles_df, show_col_types = FALSE)
  }

  initialize_gee(gee_project)

  locs <- profiles_df |>
    mutate(profile_id = as.character(profile_id)) |>
    distinct(profile_id, dataset, latitude, longitude) |>
    filter(!is.na(latitude), !is.na(longitude))

  message(sprintf(
    "[EMB] Extracting embeddings at %d global core locations (years: %d–%d)...",
    nrow(locs), min(years), max(years)
  ))

  emb_img <- .build_emb_mean(years)

  # Build FeatureCollection from lat/lon pairs
  features <- lapply(seq_len(nrow(locs)), function(i) {
    row  <- locs[i, ]
    geom <- ee$Geometry$Point(c(as.numeric(row$longitude),
                                as.numeric(row$latitude)))
    ee$Feature(geom, list(profile_id = as.character(row$profile_id),
                          dataset    = as.character(row$dataset)))
  })
  fc <- ee$FeatureCollection(features)

  # reduceRegions: consistent with rest of extraction pipeline
  result_fc <- emb_img$reduceRegions(
    collection = fc,
    reducer    = ee$Reducer$first(),
    scale      = .EMB_SCALE,
    tileScale  = 4L
  )

  data <- result_fc$getInfo()$features
  if (length(data) == 0L)
    stop("[EMB] No embedding values returned — check collection availability and years.")

  rows <- lapply(data, function(f) {
    props <- f$properties
    props <- props[setdiff(names(props), c("system:index", ".geo", "first"))]
    props <- lapply(props, function(v) if (is.null(v)) NA else v)
    as.data.frame(props, stringsAsFactors = FALSE)
  })
  result <- dplyr::bind_rows(rows)

  # Standardise embedding band names to emb_1 … emb_64
  result <- .rename_emb_bands(result)

  n_emb <- sum(grepl("^emb_", names(result)))
  message(sprintf("[EMB] Extracted %d profiles × %d embedding bands",
                  nrow(result), n_emb))
  result
}


# -----------------------------------------------------------------------------
# extract_aoi_embedding_raster()
# -----------------------------------------------------------------------------
# Downloads the mean 64-band embedding image for the local AOI extent and
# writes it directly to outputs/embedding/aoi_embedding_raster.tif.
# No Google Drive authentication required — uses GEE's getDownloadURL path
# (via = "getInfo"), which is suitable for AOIs up to ~32 MB uncompressed.
#
# covar_file : path to local covariate raster (defines AOI extent + CRS)
# years      : integer vector of years to average
# -----------------------------------------------------------------------------
extract_aoi_embedding_raster <- function(covar_file, gee_project = NULL,
                                          years = 2023:2025) {
  suppressPackageStartupMessages({ library(rgee); library(terra) })

  initialize_gee(gee_project)   # no drive needed

  # Derive AOI bounding box in WGS84 from the covariate raster
  r    <- rast(covar_file)
  r_ll <- project(r[[1]], "EPSG:4326")
  e    <- ext(r_ll)
  aoi  <- ee$Geometry$BBox(e$xmin, e$ymin, e$xmax, e$ymax)

  message(sprintf(
    "[EMB] AOI: lon [%.4f, %.4f]  lat [%.4f, %.4f]",
    e$xmin, e$xmax, e$ymin, e$ymax
  ))
  message(sprintf("[EMB] Downloading AOI embedding raster (%d bands, %d–%d avg)...",
                  .EMB_N_BANDS, min(years), max(years)))

  emb_img <- .build_emb_mean(years)$clip(aoi)

  out_dir  <- file.path("outputs", "embedding")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(out_dir, "aoi_embedding_raster.tif")

  # via = "getInfo": direct download to local path, no Drive required.
  # If the AOI is very large (> ~32 MB uncompressed) this will error with a
  # size limit message — reduce scale or switch to via = "drive".
  result <- rgee::ee_as_raster(
    image      = emb_img,
    region     = aoi,
    dsn        = out_path,
    scale      = .EMB_SCALE,
    via        = "getInfo",
    lazy       = FALSE,
    quiet      = TRUE
  )

  result <- if (inherits(result, c("RasterBrick", "RasterStack", "RasterLayer"))) {
    terra::rast(result)
  } else {
    terra::rast(out_path)
  }

  if (nlyr(result) == .EMB_N_BANDS) {
    names(result) <- paste0("emb_", seq_len(.EMB_N_BANDS))
  } else {
    message(sprintf("[EMB] Warning: expected %d bands, got %d — check collection.",
                    .EMB_N_BANDS, nlyr(result)))
  }

  message(sprintf("[EMB] Written to %s  (%d × %d px × %d bands)",
                  out_path, nrow(result), ncol(result), nlyr(result)))
  result
}
