# R/data_prep.R
# ============================================================================
# PURPOSE: Load and merge raw core location and sample CSVs.
#
# INPUTS:
#   locations_path — file path to core_locations.csv
#   samples_path   — file path to core_samples.csv
#
# OUTPUT:
#   Merged data frame: one row per sample, all location columns joined.
#   No QC flags. No carbon stock calculations. Just clean merged data.
#
# WHY SEPARATE FROM QC?
#   - You can inspect raw merged data (tar_load(cores_raw)) before QC
#   - If QC thresholds change, only the qc target re-runs, not loading
#   - Errors in loading vs errors in QC are immediately distinguishable
# ============================================================================
load_raw_data <- function(locations_path, samples_path) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
  })
  message("[data_prep] Reading raw CSVs...")
  locations <- read_csv(locations_path, show_col_types = FALSE)
  samples   <- read_csv(samples_path,   show_col_types = FALSE)

  # ── Validate required columns ──────────────────────────────────────────────
  required_loc <- c("core_id", "longitude", "latitude", "stratum")
  required_smp <- c("core_id", "depth_top_cm", "depth_bottom_cm", "soc_g_kg")
  missing_loc <- setdiff(required_loc, names(locations))
  missing_smp <- setdiff(required_smp, names(samples))
  if (length(missing_loc) > 0)
    stop("core_locations.csv is missing required columns: ",
         paste(missing_loc, collapse = ", "))
  if (length(missing_smp) > 0)
    stop("core_samples.csv is missing required columns: ",
         paste(missing_smp, collapse = ", "))

  # ── Merge and coerce types ─────────────────────────────────────────────────
  # suppressWarnings on as.numeric() allows "N/A" strings to become NA
  # without flooding the console. The QC step will flag these.
  cores <- samples %>%
    left_join(locations, by = "core_id") %>%
    mutate(
      depth_top_cm       = suppressWarnings(as.numeric(depth_top_cm)),
      depth_bottom_cm    = suppressWarnings(as.numeric(depth_bottom_cm)),
      soc_g_kg           = suppressWarnings(as.numeric(soc_g_kg)),
      bulk_density_g_cm3 = suppressWarnings(as.numeric(bulk_density_g_cm3)),
      depth_cm           = (depth_top_cm + depth_bottom_cm) / 2,
      layer_thickness_cm = depth_bottom_cm - depth_top_cm
    )

  message(sprintf("[data_prep] Loaded: %d samples from %d cores.",
                  nrow(cores), n_distinct(cores$core_id)))
  cores
}
