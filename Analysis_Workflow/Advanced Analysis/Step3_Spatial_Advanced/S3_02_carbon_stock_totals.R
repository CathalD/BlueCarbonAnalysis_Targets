# ============================================================================
# STEP 3-02: CARBON STOCK TOTALS
# ============================================================================
# PURPOSE: Aggregate per-depth kriging and RF prediction rasters into total
#          carbon stocks per stratum and AOI-wide, with VM0033-compliant
#          conservative estimates. Also reads Stage 2a outputs if available.
#
# INPUTS:
#   - outputs/predictions/kriging/aggregated_stocks/*.tif  (Mg C/ha, from S2b kriging)
#   - outputs/predictions/rf/stock_layer*.tif              (Mg C/ha or kg/m², from S3_01)
#   - AOI_FILE (GeoJSON/SHP with stratum column)
#   - outputs/Stage2a_SimpleExtrapolation/tables/ (if available)
#
# OUTPUTS (all in outputs/carbon_stocks/):
#   - carbon_stocks_by_stratum_kriging.csv
#   - carbon_stocks_overall_kriging.csv
#   - carbon_stocks_conservative_vm0033_kriging.csv
#   - carbon_stocks_by_stratum_rf.csv
#   - carbon_stocks_overall_rf.csv
#   - carbon_stocks_conservative_vm0033_rf.csv
#   - carbon_stocks_method_comparison.csv
#   - outputs/Stage3_AdvancedSpatial/Carbon_Stock_Totals/carbon_stock_totals_all_methods.csv
# ============================================================================

# ── PATH RESOLVER ────────────────────────────────────────────────────────────
local({
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) {
    setwd(dirname(dirname(dirname(dirname(this_file)))))
    return()
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NULL)
    if (!is.null(active) && nchar(active) > 0) {
      setwd(dirname(dirname(dirname(dirname(normalizePath(active))))))
      return()
    }
  }
  message("⚠ Could not auto-detect project root. Please run:\n",
          "  setwd('/path/to/BlueCarbon_Workflow_V1.0')\nbefore sourcing this script.")
})
# ─────────────────────────────────────────────────────────────────────────────

if (file.exists("Analysis_Workflow/blue_carbon_config.R")) {
  source("Analysis_Workflow/blue_carbon_config.R")
} else {
  stop("Configuration file not found.")
}

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  STEP 3-02 — CARBON STOCK TOTALS                        ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

# ── Logging ──────────────────────────────────────────────────────────────────
log_file <- file.path("logs", paste0("carbon_stock_totals_", Sys.Date(), ".log"))
if (!dir.exists("logs")) dir.create("logs", recursive = TRUE)

log_msg <- function(msg, level = "INFO") {
  entry <- sprintf("[%s] %s: %s", format(Sys.time(), "%H:%M:%S"), level, msg)
  cat(entry, "\n")
  cat(entry, "\n", file = log_file, append = TRUE)
}
log_msg("=== STEP 3-02: CARBON STOCK TOTALS STARTED ===")

# ── Packages ──────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(readr)
})

# ── Output directories ────────────────────────────────────────────────────────
dir.create("outputs/carbon_stocks", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Advanced_analysis/Step3_Spatial_Advanced/Carbon_Stock_Totals",
           recursive = TRUE, showWarnings = FALSE)

# ── Config helpers ────────────────────────────────────────────────────────────
processing_crs     <- if (exists("PROCESSING_CRS")) PROCESSING_CRS else 3347
aoi_stratum_field  <- if (exists("AOI_STRATUM_FIELD")) AOI_STRATUM_FIELD else "stratum"
aoi_available      <- exists("AOI_FILE") && !is.null(AOI_FILE) && file.exists(AOI_FILE)

# ============================================================================
# HELPER: Extract zonal statistics from a raster + AOI polygons
# ============================================================================
#' @param rast_path  Path to GeoTIFF (Mg C/ha or kg/m² — specify units)
#' @param aoi_sf     sf object with 'stratum' column, projected
#' @param units_in   "MgCha" or "kgm2" — auto-converts to Mg C/ha internally
#' @return data.frame with stratum, area_ha, mean_MgCha, total_Mg, se_total_Mg
extract_zonal_stats <- function(rast_path, aoi_sf, units_in = "MgCha") {
  if (!file.exists(rast_path)) return(NULL)

  tryCatch({
    r <- terra::rast(rast_path)
    r <- terra::project(r, paste0("EPSG:", processing_crs))

    aoi_vect <- terra::vect(aoi_sf)

    # Zonal mean per stratum polygon
    zonal_mean <- terra::extract(r, aoi_vect, fun = "mean", na.rm = TRUE, ID = TRUE)
    zonal_sd   <- terra::extract(r, aoi_vect, fun = "sd",   na.rm = TRUE, ID = TRUE)
    zonal_n    <- terra::extract(r, aoi_vect, fun = function(x) sum(!is.na(x)), ID = TRUE)

    layer_col <- names(r)[1]

    result <- aoi_sf %>%
      sf::st_drop_geometry() %>%
      mutate(ID = seq_len(n())) %>%
      left_join(zonal_mean %>% rename(mean_val = !!layer_col), by = "ID") %>%
      left_join(zonal_sd   %>% rename(sd_val   = !!layer_col), by = "ID") %>%
      left_join(zonal_n    %>% rename(n_cells  = !!layer_col), by = "ID") %>%
      select(-ID) %>%
      mutate(
        # Convert to Mg C/ha if input is kg/m²
        mean_MgCha = if (units_in == "kgm2") mean_val * 10 else mean_val,
        sd_MgCha   = if (units_in == "kgm2") sd_val   * 10 else sd_val,
        area_ha    = as.numeric(sf::st_area(aoi_sf)) / 10000,
        # Total carbon stock
        total_Mg   = mean_MgCha * area_ha,
        # Standard error of total
        se_total_Mg = (sd_MgCha / sqrt(pmax(n_cells, 1))) * area_ha
      ) %>%
      select(stratum, area_ha, mean_MgCha, sd_MgCha, n_cells, total_Mg, se_total_Mg)

    return(result)
  }, error = function(e) {
    log_msg(sprintf("  extract_zonal_stats failed for %s: %s", basename(rast_path), e$message), "WARNING")
    return(NULL)
  })
}

# ============================================================================
# HELPER: Build the 3 standard CSV files for a method
# ============================================================================
build_stock_csvs <- function(by_stratum_df, method_name) {

  method_safe <- tolower(gsub("[^A-Za-z0-9]", "_", method_name))

  # ── 1. By stratum ──────────────────────────────────────────────────────────
  by_stratum_out <- by_stratum_df %>%
    mutate(method = method_name) %>%
    arrange(desc(total_Mg))

  write_csv(by_stratum_out,
            sprintf("outputs/carbon_stocks/carbon_stocks_by_stratum_%s.csv", method_safe))

  # ── 2. Overall (sum across strata) ────────────────────────────────────────
  overall <- by_stratum_df %>%
    summarise(
      method         = method_name,
      n_strata       = n(),
      total_area_ha  = sum(area_ha, na.rm = TRUE),
      total_Mg       = sum(total_Mg, na.rm = TRUE),
      se_total_Mg    = sqrt(sum(se_total_Mg^2, na.rm = TRUE)),
      .groups        = "drop"
    ) %>%
    mutate(
      ci95_lower_Mg = total_Mg - 1.96 * se_total_Mg,
      ci95_upper_Mg = total_Mg + 1.96 * se_total_Mg,
      Site_Wide_Pooled = total_Mg   # alias used by S3_03_reporting.R
    )

  write_csv(overall,
            sprintf("outputs/carbon_stocks/carbon_stocks_overall_%s.csv", method_safe))

  # ── 3. VM0033 conservative (90% one-sided lower bound) ────────────────────
  conservative <- overall %>%
    mutate(
      conservative_Mg = total_Mg - 1.645 * se_total_Mg,
      note = "VM0033 conservative estimate: 90% one-sided lower confidence bound"
    ) %>%
    select(method, total_Mg, se_total_Mg, conservative_Mg, note)

  write_csv(conservative,
            sprintf("outputs/carbon_stocks/carbon_stocks_conservative_vm0033_%s.csv", method_safe))

  log_msg(sprintf("Saved 3 CSVs for method: %s (total = %.1f Mg C)", method_name, overall$total_Mg))

  return(list(
    by_stratum   = by_stratum_out,
    overall      = overall,
    conservative = conservative
  ))
}

# ============================================================================
# 1. LOAD AOI (required for zonal stats)
# ============================================================================

if (!aoi_available) {
  cat("⚠ AOI_FILE not set — carbon stock totals will be estimated from raster summaries only.\n\n")
  aoi_proj <- NULL
} else {
  aoi_raw  <- sf::st_read(AOI_FILE, quiet = TRUE)
  aoi_raw  <- aoi_raw %>% rename(stratum = !!aoi_stratum_field)
  aoi_proj <- sf::st_transform(aoi_raw, processing_crs)
  aoi_proj$area_ha <- as.numeric(sf::st_area(aoi_proj)) / 10000
  log_msg(sprintf("AOI loaded: %d features", nrow(aoi_proj)))
}

# ============================================================================
# 2. KRIGING TOTALS
# ============================================================================

kriging_agg_dir <- "outputs/predictions/kriging/aggregated_stocks"
kriging_tifs    <- list.files(kriging_agg_dir, pattern = "\\.tif$", full.names = TRUE)

method_results <- list()

if (length(kriging_tifs) > 0 && !is.null(aoi_proj)) {
  log_msg(sprintf("Processing %d kriging aggregated TIFs...", length(kriging_tifs)))

  # Sum all depth layers into one total stock raster
  tryCatch({
    stack_list <- lapply(kriging_tifs, terra::rast)
    # Resample to common grid if needed
    ref <- stack_list[[1]]
    stack_aligned <- lapply(stack_list, function(r) {
      if (!terra::compareGeom(r, ref, stopOnError = FALSE)) {
        terra::resample(r, ref, method = "bilinear")
      } else r
    })
    total_rast <- Reduce(`+`, stack_aligned)
    names(total_rast) <- "total_carbon_MgCha"

    # Zonal stats per stratum polygon
    all_stratum_stats <- lapply(seq_len(nrow(aoi_proj)), function(i) {
      poly_i <- aoi_proj[i, ]
      extract_zonal_stats(
        rast_path = tempfile(fileext = ".tif"),  # use in-memory rast instead
        aoi_sf    = poly_i,
        units_in  = "MgCha"
      )
    })

    # Better approach: extract directly
    aoi_vect  <- terra::vect(aoi_proj)
    z_mean    <- terra::extract(total_rast, aoi_vect, fun = "mean", na.rm = TRUE)
    z_sd      <- terra::extract(total_rast, aoi_vect, fun = "sd",   na.rm = TRUE)
    z_n       <- terra::extract(total_rast, aoi_vect,
                                fun = function(x) sum(!is.na(x)))

    kriging_by_stratum <- aoi_proj %>%
      sf::st_drop_geometry() %>%
      mutate(
        mean_MgCha  = z_mean[[2]],
        sd_MgCha    = z_sd[[2]],
        n_cells     = z_n[[2]],
        total_Mg    = mean_MgCha * area_ha,
        se_total_Mg = (sd_MgCha / sqrt(pmax(n_cells, 1))) * area_ha
      ) %>%
      select(stratum, area_ha, mean_MgCha, sd_MgCha, n_cells, total_Mg, se_total_Mg)

    method_results[["kriging"]] <- build_stock_csvs(kriging_by_stratum, "Kriging")

  }, error = function(e) {
    log_msg(sprintf("Kriging total processing failed: %s", e$message), "WARNING")
  })

} else {
  if (length(kriging_tifs) == 0) log_msg("No kriging aggregated TIFs found — skipping kriging totals", "WARNING")
  if (is.null(aoi_proj))         log_msg("No AOI — skipping kriging totals", "WARNING")
}

# ============================================================================
# 3. RANDOM FOREST TOTALS
# ============================================================================

rf_tif_dir <- "outputs/predictions/rf"
rf_tifs    <- list.files(rf_tif_dir, pattern = "stock_layer.*\\.tif$", full.names = TRUE)

if (length(rf_tifs) == 0) {
  # Try all TIFs in the RF directory
  rf_tifs <- list.files(rf_tif_dir, pattern = "\\.tif$", full.names = TRUE)
  rf_tifs <- rf_tifs[!grepl("(aoa|se_|variance)", basename(rf_tifs), ignore.case = TRUE)]
}

if (length(rf_tifs) > 0 && !is.null(aoi_proj)) {
  log_msg(sprintf("Processing %d RF TIFs...", length(rf_tifs)))

  tryCatch({
    # RF TIFs are in kg/m² — convert to Mg C/ha (×10) before summing
    stack_list <- lapply(rf_tifs, function(p) terra::rast(p) * 10)
    ref <- stack_list[[1]]
    stack_aligned <- lapply(stack_list, function(r) {
      if (!terra::compareGeom(r, ref, stopOnError = FALSE)) {
        terra::resample(r, ref, method = "bilinear")
      } else r
    })
    total_rf_rast <- Reduce(`+`, stack_aligned)
    names(total_rf_rast) <- "total_carbon_MgCha"

    aoi_vect <- terra::vect(aoi_proj)
    z_mean   <- terra::extract(total_rf_rast, aoi_vect, fun = "mean", na.rm = TRUE)
    z_sd     <- terra::extract(total_rf_rast, aoi_vect, fun = "sd",   na.rm = TRUE)
    z_n      <- terra::extract(total_rf_rast, aoi_vect,
                               fun = function(x) sum(!is.na(x)))

    rf_by_stratum <- aoi_proj %>%
      sf::st_drop_geometry() %>%
      mutate(
        mean_MgCha  = z_mean[[2]],
        sd_MgCha    = z_sd[[2]],
        n_cells     = z_n[[2]],
        total_Mg    = mean_MgCha * area_ha,
        se_total_Mg = (sd_MgCha / sqrt(pmax(n_cells, 1))) * area_ha
      ) %>%
      select(stratum, area_ha, mean_MgCha, sd_MgCha, n_cells, total_Mg, se_total_Mg)

    method_results[["rf"]] <- build_stock_csvs(rf_by_stratum, "Random Forest")

  }, error = function(e) {
    log_msg(sprintf("RF total processing failed: %s", e$message), "WARNING")
  })

} else {
  if (length(rf_tifs) == 0) log_msg("No RF prediction TIFs found — skipping RF totals", "WARNING")
}

# ============================================================================
# 4. SIMPLE EXTRAPOLATION (Stage 2a) — read if available
# ============================================================================

s2a_path <- "outputs/Basic_analysis/Step2_Spatial_Basic/SimpleExtrapolation/tables/simple_extrapolation_aoi_total.csv"
simple_extrap_result <- NULL

if (file.exists(s2a_path)) {
  tryCatch({
    s2a <- read_csv(s2a_path, show_col_types = FALSE)
    simple_extrap_result <- s2a %>%
      mutate(method = "Simple Extrapolation")
    log_msg("Loaded Stage 2a simple extrapolation results")
  }, error = function(e) {
    log_msg(sprintf("Could not load Stage 2a results: %s", e$message), "WARNING")
  })
}

# ============================================================================
# 5. METHOD COMPARISON TABLE
# ============================================================================

comparison_rows <- list()

# Kriging
if (!is.null(method_results[["kriging"]])) {
  ov <- method_results[["kriging"]]$overall
  cv <- method_results[["kriging"]]$conservative
  comparison_rows[["kriging"]] <- data.frame(
    method              = "Kriging",
    total_area_ha       = ov$total_area_ha,
    total_carbon_Mg     = ov$total_Mg,
    se_total_Mg         = ov$se_total_Mg,
    ci95_lower_Mg       = ov$ci95_lower_Mg,
    ci95_upper_Mg       = ov$ci95_upper_Mg,
    conservative_Mg     = cv$conservative_Mg,
    stringsAsFactors    = FALSE
  )
}

# RF
if (!is.null(method_results[["rf"]])) {
  ov <- method_results[["rf"]]$overall
  cv <- method_results[["rf"]]$conservative
  comparison_rows[["rf"]] <- data.frame(
    method              = "Random Forest",
    total_area_ha       = ov$total_area_ha,
    total_carbon_Mg     = ov$total_Mg,
    se_total_Mg         = ov$se_total_Mg,
    ci95_lower_Mg       = ov$ci95_lower_Mg,
    ci95_upper_Mg       = ov$ci95_upper_Mg,
    conservative_Mg     = cv$conservative_Mg,
    stringsAsFactors    = FALSE
  )
}

# Simple Extrapolation
if (!is.null(simple_extrap_result)) {
  comparison_rows[["simple"]] <- data.frame(
    method           = "Simple Extrapolation",
    total_area_ha    = simple_extrap_result$total_area_ha,
    total_carbon_Mg  = simple_extrap_result$total_carbon_Mg,
    se_total_Mg      = simple_extrap_result$se_total_kg / 1000,
    ci95_lower_Mg    = simple_extrap_result$total_carbon_kg_lower / 1000,
    ci95_upper_Mg    = simple_extrap_result$total_carbon_kg_upper / 1000,
    conservative_Mg  = simple_extrap_result$conservative_kg / 1000,
    stringsAsFactors = FALSE
  )
}

if (length(comparison_rows) > 0) {
  comparison_df <- bind_rows(comparison_rows)
  write_csv(comparison_df, "outputs/carbon_stocks/carbon_stocks_method_comparison.csv")

  # Copy to Stage3 folder
  write_csv(comparison_df,
            "outputs/Advanced_analysis/Step3_Spatial_Advanced/Carbon_Stock_Totals/carbon_stock_totals_all_methods.csv")

  log_msg(sprintf("Saved method comparison: %d methods", nrow(comparison_df)))

  # Print comparison
  cat("\n========================================\n")
  cat("CARBON STOCK METHOD COMPARISON\n")
  cat("========================================\n\n")
  for (i in seq_len(nrow(comparison_df))) {
    r <- comparison_df[i, ]
    cat(sprintf("%-22s  %.1f Mg C  (95%% CI: %.1f – %.1f)  Conservative: %.1f Mg C\n",
                r$method, r$total_carbon_Mg, r$ci95_lower_Mg, r$ci95_upper_Mg, r$conservative_Mg))
  }
  cat("\n")

} else {
  log_msg("No method results available — check inputs", "WARNING")
  cat("⚠ No carbon stock results were generated. Verify:\n")
  cat("  - Kriging TIFs exist in outputs/predictions/kriging/aggregated_stocks/\n")
  cat("  - RF TIFs exist in outputs/predictions/rf/\n")
  cat("  - AOI_FILE is set in blue_carbon_config.R\n\n")
}

# ============================================================================
# SUMMARY
# ============================================================================

cat("Outputs saved to:\n")
cat("  outputs/carbon_stocks/\n")
cat("  outputs/Stage3_AdvancedSpatial/Carbon_Stock_Totals/\n\n")

log_msg("=== STEP 3-02 COMPLETE ===")
