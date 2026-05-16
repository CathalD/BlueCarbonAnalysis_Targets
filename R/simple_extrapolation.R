# R/simple_extrapolation.R
# PURPOSE: Step 2a — extrapolate core carbon stocks to the full AOI.
#
# METHOD: Mean stock per stratum (from harmonized depth profiles) x stratum area.
#
# INPUTS:
#   stratum_summary — output of stratum_summary target (one row per stratum x depth)
#   cfg             — named list from load_config()
#
# OUTPUT — data frame with one row per stratum:
#   stratum, n_cores, total_stock_kg_m2, total_stock_MgC_ha
#   + area_m2, total_stock_kg, total_stock_MgC  (only if AOI file is available)
simple_extrapolation <- function(stratum_summary, cfg) {
  suppressPackageStartupMessages({ library(dplyr) })

  # Sum harmonized stocks across depth intervals to get total per-stratum mean
  per_stratum <- stratum_summary |>
    group_by(stratum) |>
    summarise(
      n_cores           = max(n_cores),
      total_stock_kg_m2 = sum(mean_stock, na.rm = TRUE),
      total_stock_MgC_ha = sum(mean_stock, na.rm = TRUE) * 10,
      .groups = "drop"
    )

  # If AOI file exists, compute total stocks by stratum area
  aoi_path <- cfg$AOI_FILE
  if (!is.null(aoi_path) && file.exists(aoi_path)) {
    suppressPackageStartupMessages({ library(sf) })
    message("[step2] Reading AOI: ", aoi_path)
    aoi <- st_read(aoi_path, quiet = TRUE)

    stratum_field <- cfg$AOI_STRATUM_FIELD
    if (is.null(stratum_field)) {
      message("[step2] AOI_STRATUM_FIELD is NULL — computing total AOI area only.")
      total_area_m2 <- sum(st_area(aoi))
      per_stratum <- per_stratum |>
        mutate(
          area_m2        = as.numeric(total_area_m2) / n(),
          total_stock_kg = total_stock_kg_m2 * area_m2,
          total_stock_MgC = total_stock_kg / 1000
        )
    } else {
      aoi_areas <- aoi |>
        mutate(area_m2 = as.numeric(st_area(geometry))) |>
        st_drop_geometry() |>
        group_by(stratum = .data[[stratum_field]]) |>
        summarise(area_m2 = sum(area_m2), .groups = "drop")

      per_stratum <- per_stratum |>
        left_join(aoi_areas, by = "stratum") |>
        mutate(
          total_stock_kg  = total_stock_kg_m2 * area_m2,
          total_stock_MgC = total_stock_kg / 1000
        )
    }
    message(sprintf("[step2] Extrapolation complete. Total AOI area: %.1f ha.",
                    sum(per_stratum$area_m2, na.rm = TRUE) / 10000))
  } else {
    message("[step2] No AOI file found — returning per-stratum density only.")
    message("[step2] Set AOI_FILE in blue_carbon_config.R to compute total stocks.")
  }

  per_stratum
}
