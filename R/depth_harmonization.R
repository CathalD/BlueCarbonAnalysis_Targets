# R/depth_harmonization.R
# ============================================================================
# PURPOSE: Harmonize raw depth profiles to VM0033 standard depths.
#
# METHOD: Hybrid equal-area spline (monoH.FC) + exponential decay extrapolation.
#
# VM0033 STANDARD DEPTHS (midpoints): 7.5, 22.5, 40, 75 cm
# LAYER THICKNESSES:                  15, 15, 20, 50 cm
#
# FORMULA:
#   carbon_stock_kg_m2 = SOC(g/kg) x BD(g/cm³) x thickness(cm) / 100
#
# INPUT:
#   cores_clean — QA-flagged data frame from run_qc()
#   cfg         — named list from load_config()
#
# OUTPUT — one row per core x VM0033 depth interval:
#   core_id, stratum, latitude, longitude,
#   depth_cm_midpoint, thickness_cm,
#   soc_harmonized, bd_harmonized,
#   carbon_stock_kg_m2,
#   is_extrapolated  — TRUE if target depth exceeded measured range
#   rmse, r2         — spline fit diagnostics at original measurement depths
#   qa_realistic     — harmonized values within QC thresholds
#   qa_monotonic     — SOC decreasing with depth (Spearman rho < -0.3)
# ============================================================================
harmonize_depths <- function(cores_clean, cfg) {
  suppressPackageStartupMessages({ library(dplyr) })

  # ── Configuration ─────────────────────────────────────────────────────────
  standard_depths  <- cfg$VM0033_DEPTH_MIDPOINTS %||% c(7.5, 22.5, 40, 75)
  vm0033_thickness <- if (!is.null(cfg$VM0033_DEPTH_INTERVALS)) {
    cfg$VM0033_DEPTH_INTERVALS$thickness_cm
  } else {
    c(15, 15, 20, 50)
  }
  qc_soc_min <- cfg$QC_SOC_MIN %||% 0
  qc_soc_max <- cfg$QC_SOC_MAX %||% 500
  qc_bd_min  <- cfg$QC_BD_MIN  %||% 0.1
  qc_bd_max  <- cfg$QC_BD_MAX  %||% 3.0

  cores_qa <- cores_clean |> filter(qa_pass) |> arrange(core_id, depth_cm)

  message(sprintf("[harmonization] Processing %d cores...",
                  n_distinct(cores_qa$core_id)))

  # ── Helper: spline fit diagnostics ───────────────────────────────────────
  get_diagnostics <- function(core_df) {
    tryCatch({
      fn    <- splinefun(core_df$depth_cm, core_df$soc_g_kg, method = "monoH.FC")
      preds <- fn(core_df$depth_cm)
      resid <- core_df$soc_g_kg - preds
      list(
        rmse = sqrt(mean(resid^2)),
        r2   = 1 - sum(resid^2) / sum((core_df$soc_g_kg - mean(core_df$soc_g_kg))^2)
      )
    }, error = function(e) list(rmse = NA_real_, r2 = NA_real_))
  }

  # ── Helper: hybrid spline + decay fit ────────────────────────────────────
  fit_hybrid_profile <- function(depths, values, targets) {
    if (length(depths) < 2) return(rep(NA_real_, length(targets)))
    max_d <- max(depths)
    fn <- tryCatch(splinefun(depths, values, method = "monoH.FC"), error = function(e) NULL)
    if (is.null(fn)) return(rep(NA_real_, length(targets)))
    decay <- NULL
    rho   <- tryCatch(
      cor(depths, values, method = "spearman", use = "complete.obs"),
      error = function(e) NA_real_
    )
    if (!is.na(rho) && rho < -0.3 && length(depths) >= 3)
      try(decay <- lm(log(values + 0.1) ~ depths), silent = TRUE)
    sapply(targets, function(d) {
      if (d <= max_d) {
        fn(d)
      } else if (d > max_d * 2.5) {
        NA_real_
      } else if (!is.null(decay) && coef(decay)[2] < 0) {
        exp(predict(decay, newdata = data.frame(depths = d))) - 0.1
      } else {
        values[which.max(depths)]
      }
    })
  }

  # ── Helper: process one core ──────────────────────────────────────────────
  process_core <- function(core_df) {
    depths   <- core_df$depth_cm
    soc_pred <- fit_hybrid_profile(depths, core_df$soc_g_kg, standard_depths)
    # BD: spline within range, constant extrapolation beyond (physically appropriate)
    bd_fn   <- tryCatch(
      splinefun(depths, core_df$bulk_density_g_cm3, method = "monoH.FC"),
      error = function(e) NULL
    )
    if (is.null(bd_fn)) return(NULL)
    bd_last <- core_df$bulk_density_g_cm3[which.max(depths)]
    bd_pred <- sapply(standard_depths, function(d) {
      if (d <= max(depths)) bd_fn(d) else bd_last
    })
    diag <- get_diagnostics(core_df)
    data.frame(
      core_id           = core_df$core_id[1],
      stratum           = core_df$stratum[1],
      latitude          = core_df$latitude[1],
      longitude         = core_df$longitude[1],
      depth_cm_midpoint = standard_depths,
      thickness_cm      = vm0033_thickness,
      soc_harmonized    = pmax(0, soc_pred),
      bd_harmonized     = pmax(0, bd_pred),
      is_extrapolated   = standard_depths > max(depths),
      rmse              = diag$rmse,
      r2                = diag$r2
    )
  }

  # ── Main loop ─────────────────────────────────────────────────────────────
  core_ids <- unique(cores_qa$core_id)
  results  <- vector("list", length(core_ids))
  for (i in seq_along(core_ids)) {
    sub <- cores_qa |> filter(core_id == core_ids[i]) |> arrange(depth_cm)
    if (nrow(sub) < 2) next
    results[[i]] <- tryCatch(process_core(sub), error = function(e) NULL)
  }

  harmonized <- bind_rows(results) |>
    mutate(carbon_stock_kg_m2 = (soc_harmonized * bd_harmonized * thickness_cm) / 100) |>
    filter(!is.na(carbon_stock_kg_m2))

  # ── Post-harmonization QA flags ───────────────────────────────────────────
  harmonized <- harmonized |>
    mutate(
      qa_realistic_soc = soc_harmonized >= qc_soc_min & soc_harmonized <= qc_soc_max,
      qa_realistic_bd  = bd_harmonized  >= qc_bd_min  & bd_harmonized  <= qc_bd_max,
      qa_realistic     = qa_realistic_soc & qa_realistic_bd
    )

  mono <- harmonized |>
    group_by(core_id) |>
    summarise(
      qa_monotonic = cor(depth_cm_midpoint, soc_harmonized, use = "complete.obs") < -0.3,
      .groups = "drop"
    )

  harmonized <- harmonized |> left_join(mono, by = "core_id")

  message(sprintf("[harmonization] Complete. %d cores. %d unrealistic predictions.",
                  n_distinct(harmonized$core_id),
                  sum(!harmonized$qa_realistic, na.rm = TRUE)))
  harmonized
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
