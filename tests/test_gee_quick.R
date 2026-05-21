# =============================================================================
# tests/test_gee_quick.R
# Quick smoke test for GEE covariate extraction.
#
# Run this BEFORE submitting the full Janousek extraction to confirm that
# every extraction function returns non-NA values for known wetland sites.
#
# Usage:
#   source("tests/test_gee_quick.R")
#
# Prerequisites: GEE authenticated (source the setup script once if not done):
#   source("Pre-Analysis Data Preparation/GEE_R/00_setup_rgee.R")
# =============================================================================

cat("\n=== GEE extraction quick test ===\n\n")

# ── Setup ─────────────────────────────────────────────────────────────────────
if (!file.exists("_targets.R")) setwd("..")  # run from repo root

lapply(
  list.files("R/preanalysis", pattern = "\\.R$", full.names = TRUE),
  source
)

GEE_PROJECT <- "north-star-project-470316"

# ── Test points ───────────────────────────────────────────────────────────────
# Four globally distributed coastal wetland sites covering both ecosystem types.
# Chosen to be well inside known tidal marsh / seagrass / mangrove areas so
# S2 should return real values even at 30 m scale.
test_profiles <- data.frame(
  profile_id = c("test_EM_USA",  "test_EM_EUR",  "test_SG_AUS",  "test_EM_ASI"),
  dataset    = c("test",         "test",          "test",          "test"),
  latitude   = c( 31.45,          53.18,          -33.88,           1.35),
  longitude  = c(-81.28,           5.40,          151.23,         103.82),
  ecosystem  = c("EM",           "EM",            "SG",            "EM"),
  stringsAsFactors = FALSE
)

cat("Test points:\n")
print(test_profiles[, c("profile_id", "latitude", "longitude", "ecosystem")])
cat("\n")

# ── Helper: report results ────────────────────────────────────────────────────
check_result <- function(name, df, expected_cols) {
  if (!is.data.frame(df) || nrow(df) == 0) {
    cat(sprintf("  FAIL  %s — returned empty data frame\n", name))
    return(invisible(NULL))
  }
  missing_cols <- setdiff(expected_cols, names(df))
  na_cols      <- names(df)[sapply(df[, expected_cols[expected_cols %in% names(df)],
                                       drop = FALSE], function(x) all(is.na(x)))]
  if (length(missing_cols) > 0)
    cat(sprintf("  WARN  %s — missing cols: %s\n", name, paste(missing_cols, collapse = ", ")))
  if (length(na_cols) > 0)
    cat(sprintf("  WARN  %s — all-NA cols: %s\n", name, paste(na_cols, collapse = ", ")))
  if (length(missing_cols) == 0 && length(na_cols) == 0)
    cat(sprintf("  OK    %s — %d rows, all expected cols present and non-NA\n", name, nrow(df)))
  invisible(df)
}

# ── Run each extractor ────────────────────────────────────────────────────────

cat("--- Climate (TerraClimate) ---\n")
res_climate <- tryCatch(
  extract_climate(test_profiles, gee_project = GEE_PROJECT),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)
check_result("Climate", res_climate, c("MAT_C", "MAP_mm"))

cat("\n--- Topography & Channels ---\n")
res_topo <- tryCatch(
  extract_topo(test_profiles, gee_project = GEE_PROJECT),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)
check_result("Topo", res_topo,
             c("elevation_m", "slope", "elevationRelMHW", "twi",
               "dist_to_channel_m", "tidal_flat_prob", "coastal_dist_m"))

cat("\n--- Sentinel-1 SAR ---\n")
res_sar <- tryCatch(
  extract_sar(test_profiles, gee_project = GEE_PROJECT),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)
check_result("SAR", res_sar, c("VV_mean", "VH_mean", "VVVH_ratio"))

cat("\n--- NDVI_stdDev ---\n")
res_ndvi_sd <- tryCatch(
  extract_ndvi_stddev(test_profiles, gee_project = GEE_PROJECT),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)
check_result("NDVI_stdDev", res_ndvi_sd, "NDVI_stdDev")

cat("\n--- Sentinel-2 (raw + derived, combined) ---\n")
res_s2 <- tryCatch(
  extract_s2_all(test_profiles, gee_project = GEE_PROJECT),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)
check_result("S2 all", res_s2,
             c("B", "G", "R", "B5", "B6", "B7", "NIR", "SWIR1", "SWIR2",
               "NDVI_median", "LSWI_median", "mNDWI_median", "SAVI_median",
               "tidal_wetness"))

# ── Full combine + canonical check ────────────────────────────────────────────
cat("\n--- combine_covariates() ---\n")
if (!any(sapply(list(res_climate, res_topo, res_sar, res_ndvi_sd, res_s2), is.null))) {
  combined <- tryCatch(
    combine_covariates(test_profiles,
                       topo    = res_topo,
                       sar     = res_sar,
                       ndvi_sd = res_ndvi_sd,
                       s2      = res_s2,
                       climate = res_climate),
    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(combined)) {
    n_complete <- sum(complete.cases(combined[, CANONICAL_BANDS]))
    cat(sprintf("  %d/%d rows with all 27 canonical bands non-NA\n",
                n_complete, nrow(combined)))
    cat("\nFull result:\n")
    print(t(combined[, CANONICAL_BANDS]))
  }
} else {
  cat("  Skipped — one or more extraction steps failed above\n")
}

cat("\n=== Test complete ===\n\n")
