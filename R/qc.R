# R/qc.R
# ============================================================================
# PURPOSE: Apply QA/QC flags, bulk density defaults, and calculate carbon stocks.
#
# INPUTS:
#   cores_raw — merged data frame from load_raw_data()
#   cfg       — named list from load_config()
#
# OUTPUT:
#   Same data frame with additional columns:
#     qa_depth_valid     — depth_top and depth_bottom are numeric and ordered
#     qa_soc_valid       — SOC within [QC_SOC_MIN, QC_SOC_MAX]
#     qa_bd_valid        — BD within [QC_BD_MIN, QC_BD_MAX] (NA BD passes this flag)
#     qa_has_location    — latitude and longitude are present
#     qa_pass            — all four flags are TRUE
#     bd_estimated       — TRUE where BD was NA and a default was applied
#     carbon_stock_kg_m2 — SOC x BD x thickness / 100 (NA if qa_pass = FALSE)
#
# QC THRESHOLDS — from blue_carbon_config.R:
#   QC_SOC_MIN = 0      g/kg
#   QC_SOC_MAX = 500    g/kg
#   QC_BD_MIN  = 0.1    g/cm³
#   QC_BD_MAX  = 3.0    g/cm³
#
# BD DEFAULTS — from blue_carbon_config.R:
#   BD_DEFAULTS = list(IM = 0.8, NM = 0.8, MF = 0.8)
#   Applied only where bulk_density_g_cm3 is NA.
#   bd_estimated = TRUE flags these rows for transparency in downstream reporting.
# ============================================================================
run_qc <- function(cores_raw, cfg) {
  suppressPackageStartupMessages({
    library(dplyr)
  })
  message("[qc] Applying QA/QC flags...")

  # ── Pull thresholds from config ────────────────────────────────────────────
  qc_soc_min  <- cfg$QC_SOC_MIN   # 0
  qc_soc_max  <- cfg$QC_SOC_MAX   # 500
  qc_bd_min   <- cfg$QC_BD_MIN    # 0.1
  qc_bd_max   <- cfg$QC_BD_MAX    # 3.0
  bd_defaults <- cfg$BD_DEFAULTS  # list(IM = 0.8, NM = 0.8, MF = 0.8)

  # ── Validate strata ────────────────────────────────────────────────────────
  if (!is.null(cfg$VALID_STRATA)) {
    invalid_strata <- setdiff(unique(cores_raw$stratum), cfg$VALID_STRATA)
    if (length(invalid_strata) > 0) {
      warning(sprintf(
        "[qc] Stratum names not in VALID_STRATA: %s",
        paste(invalid_strata, collapse = ", ")
      ))
    }
  }

  # ── QA flags ──────────────────────────────────────────────────────────────
  cores <- cores_raw %>%
    mutate(
      # Both depth bounds must be numeric and bottom > top
      qa_depth_valid = !is.na(depth_top_cm) &
                       !is.na(depth_bottom_cm) &
                       depth_top_cm >= 0 &
                       depth_bottom_cm > depth_top_cm,
      # SOC must be present and within range
      qa_soc_valid = !is.na(soc_g_kg) &
                     soc_g_kg >= qc_soc_min &
                     soc_g_kg <= qc_soc_max,
      # BD range check — NA BD is NOT flagged here; missing BD is handled
      # separately below via defaults (bd_estimated tracks this)
      qa_bd_valid = is.na(bulk_density_g_cm3) |
                    (bulk_density_g_cm3 >= qc_bd_min &
                     bulk_density_g_cm3 <= qc_bd_max),
      # Both coordinates must be present
      qa_has_location = !is.na(latitude) & !is.na(longitude),
      # Overall pass: all four must be TRUE
      qa_pass = qa_depth_valid & qa_soc_valid & qa_bd_valid & qa_has_location
    )

  # ── Bulk density defaults ─────────────────────────────────────────────────
  # BD_DEFAULTS is a named list: list(IM = 0.8, NM = 0.8, MF = 0.8)
  #
  # We use vapply() to look up each row's stratum key in the list.
  # vapply() guarantees a numeric(1) result per element and errors clearly
  # if a stratum is not found — safer than sapply() or map() here.
  #
  # IMPORTANT: Use [[ ]] not [ ] for list lookup.
  #   bd_defaults["IM"]   returns list(IM = 0.8)  — a list, not a number
  #   bd_defaults[["IM"]] returns 0.8             — correct
  if (!is.null(bd_defaults) && length(bd_defaults) > 0) {
    cores <- cores %>%
      mutate(
        bd_estimated = is.na(bulk_density_g_cm3),
        bulk_density_g_cm3 = if_else(
          bd_estimated,
          vapply(
            stratum,
            function(s) {
              val <- bd_defaults[[s]]
              if (is.null(val)) {
                warning(sprintf(
                  "[qc] No BD default for stratum '%s'. Returning NA.", s
                ))
                NA_real_
              } else {
                as.numeric(val)
              }
            },
            numeric(1)
          ),
          bulk_density_g_cm3
        )
      )
  } else {
    cores <- cores %>% mutate(bd_estimated = is.na(bulk_density_g_cm3))
  }

  # ── Carbon stock calculation ───────────────────────────────────────────────
  # Formula: SOC (g/kg) x BD (g/cm³) x thickness (cm) / 100 = kg C/m²
  # Only calculated for QA-passed rows — NA returned otherwise.
  cores <- cores %>%
    mutate(
      carbon_stock_kg_m2 = if_else(
        qa_pass,
        (soc_g_kg * bulk_density_g_cm3 * layer_thickness_cm) / 100,
        NA_real_
      )
    )

  # ── Summary to console ────────────────────────────────────────────────────
  n_total  <- nrow(cores)
  n_pass   <- sum(cores$qa_pass, na.rm = TRUE)
  n_bd_est <- sum(cores$bd_estimated, na.rm = TRUE)
  message(sprintf("[qc] Total: %d | Pass: %d (%.1f%%) | Fail: %d",
                  n_total, n_pass, round(100 * n_pass / n_total, 1), n_total - n_pass))
  message(sprintf("[qc] BD defaults applied: %d samples (bd_estimated = TRUE)", n_bd_est))
  if ((n_total - n_pass) > 0) {
    message(sprintf("[qc]   Depth invalid:    %d", sum(!cores$qa_depth_valid, na.rm = TRUE)))
    message(sprintf("[qc]   SOC out of range: %d", sum(!cores$qa_soc_valid,   na.rm = TRUE)))
    message(sprintf("[qc]   BD out of range:  %d", sum(!cores$qa_bd_valid,    na.rm = TRUE)))
    message(sprintf("[qc]   Missing location: %d", sum(!cores$qa_has_location, na.rm = TRUE)))
  }
  cores
}
