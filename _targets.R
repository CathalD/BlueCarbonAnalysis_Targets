library(targets)
library(tarchetypes)

tar_option_set(
  packages = c("dplyr", "readr", "tidyr", "ggplot2", "sf", "terra", "randomForest")
)

tar_source("R/")

# Usage:
#   targets::tar_make()          — or tm() shortcut in .Rprofile
#   Then run _targets_rf.R if a covariate raster is available.

list(
  # ── INPUT FILE TRACKING ───────────────────────────────────────────────────
  tar_target(locations_file, "Pre-Analysis Data Preparation/data_raw/core_locations.csv", format = "file"),
  tar_target(samples_file,   "Pre-Analysis Data Preparation/data_raw/core_samples.csv",   format = "file"),
  tar_target(config_file,    "blue_carbon_config.R",                                       format = "file"),

  # ── CONFIGURATION ─────────────────────────────────────────────────────────
  tar_target(cfg, load_config(config_file)),

  # ── STEP 1: LOAD, PREPARE, HARMONIZE ──────────────────────────────────────
  tar_target(cores_raw,        load_raw_data(locations_file, samples_file, cfg)),
  tar_target(eda_plots,        run_eda(cores_raw, cfg)),
  tar_target(cores_harmonized, harmonize_depths(cores_raw, cfg)),
  tar_target(stratum_summary,  summarise_strata(cores_harmonized)),

  # ── STEP 2: SIMPLE EXTRAPOLATION ──────────────────────────────────────────
  tar_target(step2_extrapolation, simple_extrapolation(stratum_summary, cfg)),

  # ── REPORT ────────────────────────────────────────────────────────────────
  tar_quarto(report_nonspatial, path = "reports/step1_nonspatial.qmd")
)
