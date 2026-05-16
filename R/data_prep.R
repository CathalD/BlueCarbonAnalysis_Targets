# R/data_prep.R
# PURPOSE: Load, merge, and prepare core data for analysis.
#   1. Reads core_locations.csv and core_samples.csv
#   2. Merges on core_id
#   3. Fills missing bulk density from cfg$BD_DEFAULTS (by stratum)
#   4. Computes carbon_stock_kg_m2 = SOC x BD x thickness / 100
load_raw_data <- function(locations_path, samples_path, cfg = NULL) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
  })
  message("[data_prep] Reading raw CSVs...")
  locations <- read_csv(locations_path, show_col_types = FALSE)
  samples   <- read_csv(samples_path,   show_col_types = FALSE)

  required_loc <- c("core_id", "longitude", "latitude", "stratum")
  required_smp <- c("core_id", "depth_top_cm", "depth_bottom_cm", "soc_g_kg")
  missing_loc  <- setdiff(required_loc, names(locations))
  missing_smp  <- setdiff(required_smp, names(samples))
  if (length(missing_loc) > 0)
    stop("core_locations.csv is missing columns: ", paste(missing_loc, collapse = ", "))
  if (length(missing_smp) > 0)
    stop("core_samples.csv is missing columns: ", paste(missing_smp, collapse = ", "))

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

  # Apply BD defaults where bulk_density_g_cm3 is missing
  bd_defaults <- cfg$BD_DEFAULTS
  if (!is.null(bd_defaults) && length(bd_defaults) > 0) {
    cores <- cores %>%
      mutate(
        bd_estimated = is.na(bulk_density_g_cm3),
        bulk_density_g_cm3 = if_else(
          bd_estimated,
          vapply(stratum, function(s) {
            val <- bd_defaults[[s]]
            if (is.null(val)) NA_real_ else as.numeric(val)
          }, numeric(1)),
          bulk_density_g_cm3
        )
      )
  } else {
    cores <- cores %>% mutate(bd_estimated = is.na(bulk_density_g_cm3))
  }

  cores <- cores %>%
    mutate(
      carbon_stock_kg_m2 = (soc_g_kg * bulk_density_g_cm3 * layer_thickness_cm) / 100
    )

  message(sprintf("[data_prep] %d samples from %d cores. BD defaults applied: %d.",
                  nrow(cores), n_distinct(cores$core_id), sum(cores$bd_estimated, na.rm = TRUE)))
  cores
}
