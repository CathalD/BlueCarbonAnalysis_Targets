library(targets)
library(tarchetypes)

tar_option_set(
  packages = c("dplyr", "readr", "tidyr", "ggplot2", "sf")
)

tar_source("R/")

list(
  # ── INPUT FILE TRACKING ───────────────────────────────────────────────────
  tar_target(
    locations_file,
    "Pre-Analysis Data Preparation/data_raw/core_locations.csv",
    format = "file"
  ),
  tar_target(
    samples_file,
    "Pre-Analysis Data Preparation/data_raw/core_samples.csv",
    format = "file"
  ),

  # ── CONFIGURATION ─────────────────────────────────────────────────────────
  tar_target(cfg, load_config("blue_carbon_config.R")),

  # ── STEP 1a: LOAD & PREPARE DATA ──────────────────────────────────────────
  tar_target(
    cores_raw,
    load_raw_data(locations_path = locations_file, samples_path = samples_file, cfg = cfg)
  ),

  # ── STEP 1b: EDA AND HARMONIZATION (parallel) ─────────────────────────────
  tar_target(eda_plots,        run_eda(cores_raw, cfg)),
  tar_target(cores_harmonized, harmonize_depths(cores_raw, cfg)),

  # ── STEP 1c: SUMMARY STATISTICS ───────────────────────────────────────────
  tar_target(
    stratum_summary,
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
  ),

  # ── STEP 2a: SIMPLE EXTRAPOLATION ─────────────────────────────────────────
  tar_target(step2_extrapolation, simple_extrapolation(stratum_summary, cfg)),

  # ── REPORT ────────────────────────────────────────────────────────────────
  tar_quarto(report_step1, path = "reports/step1_nonspatial.qmd")
)
