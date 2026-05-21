# =============================================================================
# R/preanalysis/gee_covariates.R
# GEE covariate extraction — exact R/rgee port of
#   CoastalBlueCarbon_GlobalCoreCovariate_Extraction.ipynb
#
# Produces the canonical 27-band stack:
#   Group 1 — Topography & Channels (7):
#     elevation_m, slope, elevationRelMHW, twi, dist_to_channel_m,
#     tidal_flat_prob, coastal_dist_m
#   Group 2 — Sentinel-1 SAR (3):
#     VV_mean, VH_mean, VVVH_ratio
#   Group 3 — Sentinel-2 Optical & Phenology (15):
#     B, G, R, B5, B6, B7, NIR, SWIR1, SWIR2
#     NDVI_median, LSWI_median, mNDWI_median, NDVI_stdDev, SAVI_median, tidal_wetness
#   Group 4 — Climate (2):
#     MAT_C, MAP_mm
#
# All extraction parameters (date ranges, cloud thresholds, batch sizes,
# scale, season filter) match the Python notebook exactly.
# =============================================================================

# Canonical band order — must match GoogleEarthEngineAOICovariateAnalysis.js
# and CoastalBlueCarbon_LargeScaleCovariateExtraction.ipynb
CANONICAL_BANDS <- c(
  "elevation_m", "slope", "elevationRelMHW", "twi", "dist_to_channel_m",
  "tidal_flat_prob", "coastal_dist_m",
  "VV_mean", "VH_mean", "VVVH_ratio",
  "B", "G", "R", "B5", "B6", "B7", "NIR", "SWIR1", "SWIR2",
  "NDVI_median", "LSWI_median", "mNDWI_median",
  "NDVI_stdDev", "SAVI_median", "tidal_wetness",
  "MAT_C", "MAP_mm"
)
stopifnot(length(CANONICAL_BANDS) == 27L)

# GEE system columns emitted by reduceRegions — drop from all results
.GEE_SYSTEM_COLS <- c("system:index", ".geo", "first")

# ── Date ranges (must match GEE JS script) ───────────────────────────────────
.S2_START  <- "2020-01-01"
.S2_END    <- "2023-12-31"
.SAR_START <- "2020-01-01"
.SAR_END   <- "2023-12-31"
.TC_START  <- "2000-01-01"
.TC_END    <- "2022-12-31"

# ── S2 extraction parameters ─────────────────────────────────────────────────
.S2_MAX_CLOUD   <- 20L   # % cloud cover threshold (collection-level pre-filter)
.S2_IMAGE_LIMIT <- 20L   # max scenes per batch area (least cloudy first)
.S2_BUFFER_M    <- 5000L # metres buffer around each batch of points
# scale = 30 rather than 10: for point extraction at core locations the 9×
# reduction in pixel count is the dominant speedup; sub-pixel accuracy gain
# from 10 m is negligible for ecosystem-scale carbon stock prediction.
.S2_SCALE       <- 30L
.TILE_SCALE     <- 8L    # GEE tileScale — 8 gives more memory tiles for large reductions

# Bands needed inside the cloud-mask / index functions — select early so GEE
# only loads these tiles, not all 13 S2 bands.
.S2_BANDS_FULL  <- c("B2", "B3", "B4", "B5", "B6", "B7", "B8", "B11", "B12", "QA60")
.S2_BANDS_NDVI  <- c("B4", "B8", "QA60")   # only NDVI (B8-B4) needed; no tidal mask


# =============================================================================
# INTERNAL HELPERS — image stack builders
# =============================================================================

.build_topo_stack <- function() {
  dem         <- ee$Image("NASA/NASADEM_HGT/001")$select("elevation")
  elevation_m <- dem$rename("elevation_m")
  slope       <- ee$Terrain$slope(dem)$rename("slope")

  # Elevation relative to Mean High Water: DEM minus 0.5 m MHW-above-MSL proxy.
  # Consistent global approximation; matches Python notebook exactly.
  elevRelMHW <- dem$subtract(0.5)$rename("elevationRelMHW")

  # TWI: ln(upslope_area / tan(slope))
  # upslope_area proxied by counting valid DEM pixels within a 20-pixel circle
  slope_rad <- ee$Terrain$slope(dem)$multiply(pi / 180)
  tan_slope <- slope_rad$tan()$max(0.001)
  contrib   <- dem$gte(-9999)$unmask(0L)$
    reduceNeighborhood(
      reducer = ee$Reducer$sum(),
      kernel  = ee$Kernel$circle(radius = 20, units = "pixels")
    )$max(1)
  twi <- contrib$divide(tan_slope)$log()$rename("twi")

  # dist_to_channel_m: distance to JRC surface water occurrence > 30%
  # focalMax(30 m) closes canopy-gap artefacts; fastDistanceTransform × 30 → metres
  channel_mask <- ee$Image("JRC/GSW1_4/GlobalSurfaceWater")$
    select("occurrence")$gt(30L)$unmask(0L)$
    focalMax(30, "circle", "meters")
  dist_to_channel <- channel_mask$
    fastDistanceTransform(500L, "pixels", "squared_euclidean")$
    sqrt()$multiply(30)$rename("dist_to_channel_m")$float()

  # tidal_flat_prob: Murray et al. 2019 intertidal classification (ImageCollection)
  tidal_flat_prob <- tryCatch({
    ee$ImageCollection("UQ/murray/Intertidal/v1_1/global_intertidal")$
      filterBounds(ee$Geometry$BBox(-180, -90, 180, 90))$
      mosaic()$
      select("classification")$eq(1L)$unmask(0L)$
      rename("tidal_flat_prob")$float()
  }, error = function(e) {
    message("[GEE] Murray intertidal unavailable — tidal_flat_prob set to 0")
    ee$Image(0)$rename("tidal_flat_prob")$float()
  })

  # coastal_dist_m: distance to JRC occurrence > 50% (open water)
  water_mask   <- ee$Image("JRC/GSW1_4/GlobalSurfaceWater")$
    select("occurrence")$gt(50L)$unmask(0L)
  coastal_dist <- water_mask$
    fastDistanceTransform(500L, "pixels", "squared_euclidean")$
    sqrt()$multiply(30)$rename("coastal_dist_m")$float()

  elevation_m$addBands(slope)$addBands(elevRelMHW)$
    addBands(twi)$addBands(dist_to_channel)$
    addBands(tidal_flat_prob)$addBands(coastal_dist)
}


.build_sar_stack <- function() {
  s1_col <- ee$ImageCollection("COPERNICUS/S1_GRD")$
    filterDate(.SAR_START, .SAR_END)$
    filter(ee$Filter$listContains("transmitterReceiverPolarisation", "VV"))$
    filter(ee$Filter$listContains("transmitterReceiverPolarisation", "VH"))$
    filter(ee$Filter$eq("instrumentMode", "IW"))$
    map(function(img) img$updateMask(img$select("VV")$gt(-30)))  # pixel-level noise mask

  s1_mean <- s1_col$mean()
  vv      <- s1_mean$select("VV")$rename("VV_mean")
  vh      <- s1_mean$select("VH")$rename("VH_mean")
  vvvh    <- s1_mean$select("VV")$subtract(s1_mean$select("VH"))$rename("VVVH_ratio")

  vv$addBands(vh)$addBands(vvvh)
}


.build_climate_stack <- function() {
  terra_mean <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")$
    filterDate(.TC_START, .TC_END)$
    select(c("tmmn", "tmmx", "pr"))$
    mean()

  # MAT (°C): raw units are °C × 10
  mat_img <- terra_mean$expression(
    "((tmmn + tmmx) / 2.0) / 10.0",
    list("tmmn" = terra_mean$select("tmmn"),
         "tmmx" = terra_mean$select("tmmx"))
  )$rename("MAT_C")

  # MAP (mm/year): raw is mm/month → × 12
  map_img <- terra_mean$select("pr")$multiply(12)$rename("MAP_mm")

  mat_img$addBands(map_img)
}


# Per-batch NDVI_stdDev: spatially filtered to each batch's bounding box.
# stdDev uses the FULL unranked cloud-masked collection (no image limit) to
# capture the complete range of phenological variability across the season.
#
# NOTE: No tidal mask — same reasoning as .build_s2_median(). Cloud mask only.
# NOTE: filterBounds(region) is REQUIRED here. Without it GEE must compute
#       stdDev across the global S2 archive → computation timeout for every batch.
.build_ndvi_stddev_local <- function(region) {
  process <- function(image) {
    qa         <- image$select("QA60")
    cloud_mask <- qa$bitwiseAnd(1024L)$eq(0L)$And(qa$bitwiseAnd(2048L)$eq(0L))
    image$divide(10000)$
      updateMask(cloud_mask)$
      normalizedDifference(c("B8", "B4"))$rename("NDVI")  # → stdDev gives "NDVI_stdDev"
  }

  ee$ImageCollection("COPERNICUS/S2_SR_HARMONIZED")$
    filterDate(.S2_START, .S2_END)$
    filterBounds(region)$                                              # spatial filter first
    filter(ee$Filter$lt("CLOUDY_PIXEL_PERCENTAGE", .S2_MAX_CLOUD))$
    filter(ee$Filter$calendarRange(5, 9, "month"))$
    select(.S2_BANDS_NDVI)$
    map(process)$
    reduce(ee$Reducer$stdDev())$
    rename("NDVI_stdDev")
}


# Per-batch S2 median spatially filtered to the batch's bounding box + buffer.
# This prevents GEE from computing a full global mosaic per call.
#
# NOTE: No tidal mask applied — see .build_ndvi_stddev_img() for rationale.
# Tidal wetland cores (EM/SG) must include inundated observations; cloud-only
# masking produces an unbiased median over the full inundation cycle.
#
# Performance notes applied:
#   ① filterBounds(region) BEFORE .map() — spatial filter before per-image processing
#   ② CLOUDY_PIXEL_PERCENTAGE < 20 at collection level
#   ③ .select(.S2_BANDS_FULL) BEFORE .map() — load only needed bands per tile
#   ④ QA60 bitmask + ÷10000 scaling inside process()
.build_s2_median <- function(region) {
  process <- function(image) {
    qa         <- image$select("QA60")
    cloud_mask <- qa$bitwiseAnd(1024L)$eq(0L)$And(qa$bitwiseAnd(2048L)$eq(0L))
    image$divide(10000)$
      updateMask(cloud_mask)
  }

  ee$ImageCollection("COPERNICUS/S2_SR_HARMONIZED")$
    filterDate(.S2_START, .S2_END)$
    filterBounds(region)$                                              # ① spatial filter first
    filter(ee$Filter$lt("CLOUDY_PIXEL_PERCENTAGE", .S2_MAX_CLOUD))$   # ② collection pre-filter
    filter(ee$Filter$calendarRange(5, 9, "month"))$
    limit(.S2_IMAGE_LIMIT, "CLOUDY_PIXEL_PERCENTAGE")$
    select(.S2_BANDS_FULL)$                                           # ③ band selection before map
    map(process)$                                                      # ④ cloud mask + scale inside
    median()$clip(region)
}


# ── Band selector functions for S2 ──────────────────────────────────────────

# 9 raw reflectance bands (includes Red-Edge B5/B6/B7 for seagrass/saltmarsh)
.s2_select_raw <- function(s2) {
  s2$select("B2")$rename("B")$
    addBands(s2$select("B3")$rename("G"))$
    addBands(s2$select("B4")$rename("R"))$
    addBands(s2$select("B5"))$   # Red-Edge 705 nm — already named B5
    addBands(s2$select("B6"))$   # Red-Edge 740 nm — already named B6
    addBands(s2$select("B7"))$   # Red-Edge 783 nm — already named B7
    addBands(s2$select("B8")$rename("NIR"))$
    addBands(s2$select("B11")$rename("SWIR1"))$
    addBands(s2$select("B12")$rename("SWIR2"))
}

# 5 derived indices (EVI removed; tidal_wetness = Nedkov 2017 TC Wetness)
.s2_select_derived <- function(s2) {
  ndvi <- s2$normalizedDifference(c("B8", "B4"))$rename("NDVI_median")
  lswi <- s2$normalizedDifference(c("B8", "B11"))$rename("LSWI_median")
  mndwi <- s2$normalizedDifference(c("B3", "B11"))$rename("mNDWI_median")
  savi <- s2$expression(
    "((NIR - RED) / (NIR + RED + 0.5)) * 1.5",
    list("NIR" = s2$select("B8"), "RED" = s2$select("B4"))
  )$rename("SAVI_median")
  # Tasseled Cap Wetness — Nedkov 2017 Sentinel-2 SR coefficients
  tidal_wetness <- s2$expression(
    "0.1511*B + 0.1973*G + 0.3283*R + 0.3407*NIR + (-0.7117)*SWIR1 + (-0.4559)*SWIR2",
    list(
      "B"     = s2$select("B2"),
      "G"     = s2$select("B3"),
      "R"     = s2$select("B4"),
      "NIR"   = s2$select("B8"),
      "SWIR1" = s2$select("B11"),
      "SWIR2" = s2$select("B12")
    )
  )$rename("tidal_wetness")

  ndvi$addBands(lswi)$addBands(mndwi)$addBands(savi)$addBands(tidal_wetness)
}


# =============================================================================
# INTERNAL HELPERS — batch extraction engine
# =============================================================================

# Build an ee.FeatureCollection directly from a data.frame's lat/lon columns.
# Mirrors the Python notebook's loop exactly — no sf or geojsonio required.
.df_to_ee_fc <- function(df) {
  features <- lapply(seq_len(nrow(df)), function(i) {
    row  <- df[i, ]
    geom <- ee$Geometry$Point(c(as.numeric(row$longitude), as.numeric(row$latitude)))
    ee$Feature(geom, list(
      profile_id = as.character(row$profile_id),
      dataset    = as.character(row$dataset)
    ))
  })
  ee$FeatureCollection(features)
}


# Convert a GEE FeatureCollection result to a data.frame.
# use_drive = TRUE: async Drive export (for large result sets or slow connections).
# use_drive = FALSE (default): synchronous getInfo() (fine for ≤ ~200 features/call).
#
# Drive export requires ee_Initialize(drive = TRUE) to have been called first.
.ee_fc_result_to_df <- function(result_fc, name, use_drive = FALSE) {
  suppressPackageStartupMessages(library(dplyr))

  if (use_drive) {
    task_name <- paste0("rgee_", gsub("[^a-zA-Z0-9]", "_", substr(name, 1L, 40L)))
    task <- rgee::ee_table_to_drive(
      collection  = result_fc,
      description = task_name,
      fileFormat  = "CSV",
      folder      = "rgee_exports"
    )
    task$start()
    message(sprintf("  [Drive] task '%s' submitted — waiting...", task_name))
    rgee::ee_monitoring(task, max_attempts = 200L, quiet = TRUE)

    local_file <- tempfile(fileext = ".csv")
    rgee::ee_drive_to_local(task, dsn = local_file, quiet = TRUE)
    df <- readr::read_csv(local_file, show_col_types = FALSE)
    df <- df[, setdiff(names(df), .GEE_SYSTEM_COLS), drop = FALSE]
    return(df)
  }

  # Synchronous path
  data <- result_fc$getInfo()$features
  if (length(data) == 0L) return(data.frame())
  rows <- lapply(data, function(f) {
    props <- f$properties
    props <- props[setdiff(names(props), .GEE_SYSTEM_COLS)]
    # GEE returns NULL for masked / out-of-extent pixels; coerce to NA so
    # as.data.frame() doesn't see mixed lengths (0 vs 1) and error out.
    props <- lapply(props, function(v) if (is.null(v)) NA else v)
    as.data.frame(props, stringsAsFactors = FALSE)
  })
  dplyr::bind_rows(rows)
}


# Generic batched reduceRegions for a pre-built ee.Image.
# Returns data.frame with profile_id + extracted band columns.
#
# use_drive: passed to .ee_fc_result_to_df() — set TRUE for large point sets.
.extract_batch <- function(image, profiles_df, name,
                            batch_size = 100L, scale = 30L, use_drive = FALSE) {
  suppressPackageStartupMessages(library(dplyr))

  n         <- nrow(profiles_df)
  n_batches <- ceiling(n / batch_size)
  all_rows  <- list()
  n_failed  <- 0L

  message(sprintf("[GEE] Extracting %s (%d pts, batch=%d, scale=%dm, tileScale=%d)",
                  name, n, batch_size, scale, .TILE_SCALE))

  for (i in seq(1L, n, by = batch_size)) {
    end_idx   <- min(i + batch_size - 1L, n)
    batch_df  <- profiles_df[i:end_idx, ]
    fc        <- .df_to_ee_fc(batch_df)
    batch_num <- ceiling(i / batch_size)

    tryCatch({
      result_fc <- image$reduceRegions(
        collection = fc,
        reducer    = ee$Reducer$first(),
        scale      = scale,
        tileScale  = .TILE_SCALE   # ④ tileScale = 4 (4× more GEE memory tiles)
      )
      batch_df_out <- .ee_fc_result_to_df(result_fc, name, use_drive)
      all_rows <- c(all_rows, list(batch_df_out))

      if (batch_num %% 10L == 0L || batch_num == n_batches)
        message(sprintf("  Batch %d/%d OK  (%d rows so far)",
                        batch_num, n_batches, sum(vapply(all_rows, nrow, 0L))))
    }, error = function(e) {
      n_failed <<- n_failed + 1L
      message(sprintf("  Batch %d/%d FAILED: %s", batch_num, n_batches, conditionMessage(e)))
    })
  }

  message(sprintf("[GEE] %s complete — %d rows, %d batches failed",
                  name, sum(vapply(all_rows, nrow, 0L)), n_failed))

  if (length(all_rows) == 0L) return(data.frame())
  dplyr::bind_rows(all_rows)
}


# Combined S2 batch extraction: builds ONE spatially filtered S2 median per
# batch, then extracts all 14 bands (raw reflectance + derived indices) in a
# single reduceRegions() call.  Previously, raw and derived were separate
# targets → 2 × n_batches S2 composites computed; this halves that to n_batches.
#
# use_drive : passed to .ee_fc_result_to_df()
.extract_batch_s2_all <- function(profiles_df, batch_size = 25L,
                                   scale = .S2_SCALE, use_drive = FALSE) {
  suppressPackageStartupMessages(library(dplyr))

  n         <- nrow(profiles_df)
  n_batches <- ceiling(n / batch_size)
  all_rows  <- list()
  n_failed  <- 0L

  message(sprintf("[GEE] Extracting S2 all bands — raw (9) + derived (5) combined",))
  message(sprintf("      %d pts | batch=%d | scale=%dm | tileScale=%d | buffer=%dm",
                  n, batch_size, scale, .TILE_SCALE, .S2_BUFFER_M))

  for (i in seq(1L, n, by = batch_size)) {
    end_idx   <- min(i + batch_size - 1L, n)
    batch_df  <- profiles_df[i:end_idx, ]
    fc        <- .df_to_ee_fc(batch_df)
    batch_num <- ceiling(i / batch_size)

    tryCatch({
      region   <- fc$geometry()$bounds()$buffer(.S2_BUFFER_M)
      s2_local <- .build_s2_median(region)   # ONE composite per batch

      # Merge raw + derived into a single 14-band image before reduceRegions
      img <- .s2_select_raw(s2_local)$addBands(.s2_select_derived(s2_local))

      result_fc <- img$reduceRegions(
        collection = fc,
        reducer    = ee$Reducer$first(),
        scale      = scale,
        tileScale  = .TILE_SCALE
      )
      batch_df_out <- .ee_fc_result_to_df(result_fc, "S2 all bands", use_drive)
      all_rows <- c(all_rows, list(batch_df_out))

      if (batch_num %% 10L == 0L || batch_num == n_batches)
        message(sprintf("  Batch %d/%d OK  (%d rows so far)",
                        batch_num, n_batches, sum(vapply(all_rows, nrow, 0L))))
    }, error = function(e) {
      n_failed <<- n_failed + 1L
      message(sprintf("  Batch %d/%d FAILED: %s", batch_num, n_batches, conditionMessage(e)))
    })
  }

  message(sprintf("[GEE] S2 all bands complete — %d rows, %d batches failed",
                  sum(vapply(all_rows, nrow, 0L)), n_failed))

  if (length(all_rows) == 0L) return(data.frame())
  dplyr::bind_rows(all_rows)
}


# =============================================================================
# PUBLIC API — one target per extraction group
# =============================================================================

extract_topo <- function(profiles_df, gee_project = NULL, use_drive = FALSE) {
  suppressPackageStartupMessages(library(rgee))
  initialize_gee(gee_project)
  message("[GEE] Building topography + channels stack...")
  stack <- .build_topo_stack()
  .extract_batch(stack, profiles_df, "Topography & Channels (7 bands)",
                 batch_size = 500L, scale = 30L, use_drive = use_drive)
}

extract_sar <- function(profiles_df, gee_project = NULL, use_drive = FALSE) {
  suppressPackageStartupMessages(library(rgee))
  initialize_gee(gee_project)
  message("[GEE] Building Sentinel-1 SAR stack...")
  stack <- .build_sar_stack()
  .extract_batch(stack, profiles_df, "Sentinel-1 SAR (3 bands)",
                 batch_size = 100L, scale = 30L, use_drive = use_drive)
}

extract_ndvi_stddev <- function(profiles_df, gee_project = NULL,
                                 batch_size = 25L, use_drive = FALSE) {
  suppressPackageStartupMessages({ library(rgee); library(dplyr) })
  initialize_gee(gee_project)

  n         <- nrow(profiles_df)
  n_batches <- ceiling(n / batch_size)
  all_rows  <- list()
  n_failed  <- 0L

  message(sprintf("[GEE] Extracting NDVI_stdDev (%d pts, batch=%d, scale=30m, tileScale=%d, buffer=%dm)",
                  n, batch_size, .TILE_SCALE, .S2_BUFFER_M))

  for (i in seq(1L, n, by = batch_size)) {
    end_idx   <- min(i + batch_size - 1L, n)
    batch_df  <- profiles_df[i:end_idx, ]
    fc        <- .df_to_ee_fc(batch_df)
    batch_num <- ceiling(i / batch_size)

    tryCatch({
      region <- fc$geometry()$bounds()$buffer(.S2_BUFFER_M)
      img    <- .build_ndvi_stddev_local(region)

      result_fc <- img$reduceRegions(
        collection = fc,
        reducer    = ee$Reducer$first(),
        scale      = 30L,
        tileScale  = .TILE_SCALE
      )
      batch_df_out <- .ee_fc_result_to_df(result_fc, "NDVI_stdDev", use_drive)
      all_rows <- c(all_rows, list(batch_df_out))

      if (batch_num %% 10L == 0L || batch_num == n_batches)
        message(sprintf("  Batch %d/%d OK  (%d rows so far)",
                        batch_num, n_batches, sum(vapply(all_rows, nrow, 0L))))
    }, error = function(e) {
      n_failed <<- n_failed + 1L
      message(sprintf("  Batch %d/%d FAILED: %s", batch_num, n_batches, conditionMessage(e)))
    })
  }

  message(sprintf("[GEE] NDVI_stdDev complete — %d rows, %d batches failed",
                  sum(vapply(all_rows, nrow, 0L)), n_failed))

  if (length(all_rows) == 0L) return(data.frame())
  dplyr::bind_rows(all_rows)
}

extract_s2_all <- function(profiles_df, gee_project = NULL, use_drive = FALSE) {
  suppressPackageStartupMessages(library(rgee))
  initialize_gee(gee_project)
  .extract_batch_s2_all(profiles_df, batch_size = 25L,
                         scale = .S2_SCALE, use_drive = use_drive)
}

extract_climate <- function(profiles_df, gee_project = NULL, use_drive = FALSE) {
  suppressPackageStartupMessages(library(rgee))
  initialize_gee(gee_project)
  message(sprintf("[GEE] Building TerraClimate stack (%s–%s)...",
                  substr(.TC_START, 1, 4), substr(.TC_END, 1, 4)))
  stack <- .build_climate_stack()
  .extract_batch(stack, profiles_df,
                 sprintf("TerraClimate MAT/MAP (%s–%s)",
                         substr(.TC_START, 1, 4), substr(.TC_END, 1, 4)),
                 batch_size = 500L, scale = 4000L, use_drive = use_drive)
}


# =============================================================================
# combine_covariates()
# =============================================================================
# Left-joins all extraction results onto profiles_df, enforces canonical
# column order, and fills missing bands with NA (with a warning).
#
# Returns a data.frame ready for write_covariates_csv().
# =============================================================================
combine_covariates <- function(profiles_df, topo, sar, ndvi_sd, s2, climate) {
  suppressPackageStartupMessages(library(dplyr))

  .merge_gee <- function(main, sub) {
    if (is.null(sub) || !is.data.frame(sub) || nrow(sub) == 0L) return(main)
    sub <- sub |>
      mutate(profile_id = as.character(profile_id)) |>
      select(-any_of(c(.GEE_SYSTEM_COLS, "dataset")))
    dplyr::left_join(main, sub, by = "profile_id")
  }

  result <- profiles_df |>
    mutate(profile_id = as.character(profile_id))

  for (df in list(climate, topo, sar, s2, ndvi_sd)) {
    result <- .merge_gee(result, df)
  }

  # Flag and fill any missing canonical bands
  missing <- setdiff(CANONICAL_BANDS, names(result))
  if (length(missing) > 0L) {
    warning(sprintf("[covariates] %d canonical band(s) missing (extraction failed): %s",
                    length(missing), paste(missing, collapse = ", ")))
    for (b in missing) result[[b]] <- NA_real_
  }

  # Drop any leaked removed-band columns
  removed_bands <- c("tpi", "EVI_median", "NDBI_median", "brightness", "greenness",
                     "sg_soc_0_30cm", "sg_soc_30_100cm", "sg_soc_0_100cm", "first")
  result <- select(result, -any_of(removed_bands))

  # Enforce column order: meta columns first, then canonical bands
  meta_cols <- setdiff(names(result), CANONICAL_BANDS)
  result[, c(meta_cols, CANONICAL_BANDS)]
}


# =============================================================================
# write_covariates_csv()
# =============================================================================
# Writes the final covariate data.frame to disk and returns the path
# (format = "file" target).
# =============================================================================
write_covariates_csv <- function(global_covariates, path) {
  suppressPackageStartupMessages(library(readr))

  n_profiles <- nrow(global_covariates)
  n_cols     <- ncol(global_covariates)
  n_complete <- sum(complete.cases(global_covariates[, CANONICAL_BANDS]))
  pct_complete <- round(100 * n_complete / n_profiles, 1)

  readr::write_csv(global_covariates, path)

  message(sprintf("[covariates] Saved: %s", path))
  message(sprintf("[covariates] %d profiles × %d cols | %d/%d (%.1f%%) with complete covariates",
                  n_profiles, n_cols, n_complete, n_profiles, pct_complete))

  # Warn if NA rate is high
  for (band in CANONICAL_BANDS) {
    na_rate <- mean(is.na(global_covariates[[band]]))
    if (na_rate > 0.05)
      message(sprintf("  ⚠ %s: %.1f%% NA (check GEE extraction)", band, 100 * na_rate))
  }

  path
}
