# R/summarise.R
# PURPOSE: Compute per-stratum, per-depth summary statistics from harmonized cores.
summarise_strata <- function(cores_harmonized) {
  suppressPackageStartupMessages({ library(dplyr) })
  cores_harmonized |>
    group_by(stratum, depth_cm_midpoint) |>
    summarise(
      n_cores          = n_distinct(core_id),
      mean_stock       = mean(carbon_stock_kg_m2, na.rm = TRUE),
      sd_stock         = sd(carbon_stock_kg_m2,   na.rm = TRUE),
      mean_soc         = mean(soc_harmonized,      na.rm = TRUE),
      mean_bd          = mean(bd_harmonized,        na.rm = TRUE),
      pct_extrapolated = mean(is_extrapolated) * 100,
      .groups = "drop"
    )
}
