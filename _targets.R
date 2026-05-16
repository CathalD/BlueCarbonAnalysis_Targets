library(targets)
library(tarchetypes)

# Packages available in every target's environment.
tar_option_set(
  packages = c("dplyr", "readr", "tidyr", "ggplot2", "sf")
)

# Sources all .R files in R/ — functions become available to all targets.
tar_source("R/")

# ============================================================================
# PIPELINE PLAN
# Each tar_target(name, command) declares one step:
#   name    — the reference name (used by downstream targets as argument names)
#   command — R expression that computes the target's value
#
# targets infers execution order by matching argument names to target names.
# ============================================================================
list(
  # ── INPUT FILE TRACKING ───────────────────────────────────────────────────
  # format = "file" tracks the path AND a hash of the file's contents.
  # If the CSV changes on disk, all downstream targets become outdated.
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

  # ── STEP 1a: LOAD RAW DATA ────────────────────────────────────────────────
  # Depends on: locations_file, samples_file
  # If either CSV changes → cores_raw re-runs → everything downstream re-runs
  tar_target(
    cores_raw,
    load_raw_data(locations_path = locations_file, samples_path = samples_file)
  ),

  # ── STEP 1b: QA/QC ────────────────────────────────────────────────────────
  # Depends on: cores_raw, cfg
  # If QC_SOC_MAX changes in config → cfg re-runs → cores_clean re-runs
  # cores_raw stays cached (raw files did not change)
  tar_target(cores_clean, run_qc(cores_raw, cfg)),

  # ── STEP 1c: EDA AND HARMONIZATION (parallel branches) ───────────────────
  # These two do not depend on each other — they are parallel branches.
  # With the crew package installed, they can run simultaneously.
  tar_target(eda_plots,        run_eda(cores_clean, cfg)),
  tar_target(cores_harmonized, harmonize_depths(cores_clean, cfg)),

  # ── STEP 1d: SUMMARY STATISTICS ──────────────────────────────────────────
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

  # ── STEP 1e: QUARTO REPORT ────────────────────────────────────────────────
  # tar_quarto() renders the .qmd and detects its tar_read() dependencies
  # automatically. The report re-renders if any upstream target changes.
  tar_quarto(report_step1, path = "reports/step1_nonspatial.qmd")
)
