library(targets)
library(tarchetypes)
library(geotargets)

tar_option_set(
  packages = c("dplyr", "readr", "tidyr", "ggplot2", "sf", "terra", "ranger")
)

tar_source("R/")

# Usage:
#   Main pipeline first: targets::tar_make()                     (uses _targets.yaml "main" config)
#   Transfer pipeline:   targets::tar_make(script = "_targets_transfer.R", store = "_targets_transfer")
#   Or set config name:  targets::tar_config_set(config = "transfer"); targets::tar_make()

list(
  # ── INPUT FILE TRACKING ───────────────────────────────────────────────────
  tar_target(global_layers_file,
             "Pre-Analysis Data Preparation/data_global/combined_layers_filtered.csv",
             format = "file"),
  tar_target(global_covar_file,
             "Pre-Analysis Data Preparation/data_raw/CorePoints_Covariates_BC_Canada.csv",
             format = "file"),
  tar_target(config_file_tl, "blue_carbon_config.R", format = "file"),

  # ── CONFIGURATION ─────────────────────────────────────────────────────────
  tar_target(cfg_tl, load_config(config_file_tl)),

  # covar_file_tl path comes from config, so it depends on cfg_tl
  tar_target(covar_file_tl, cfg_tl$COVARIATE_RASTER, format = "file"),

  # ── LOCAL DATA FROM MAIN PIPELINE ────────────────────────────────────────
  # Reads cores_harmonized from the main pipeline store (_targets/).
  # Note: targets cannot auto-detect cross-store dependencies — re-run the
  # transfer pipeline manually after updating the main pipeline.
  tar_target(cores_harmonized_tl,
             targets::tar_read(cores_harmonized, store = "_targets")),

  # ── GLOBAL DATA HARMONIZATION ────────────────────────────────────────────
  tar_target(global_harmonized,
             harmonize_global_layers(global_layers_file, cfg_tl)),

  # ── COMBINED TL TRAINING DATA ────────────────────────────────────────────
  tar_target(tl_data,
             prepare_tl_data(cores_harmonized_tl, global_harmonized,
                             global_covar_file, covar_file_tl)),

  # ── TRANSFER LEARNING MODELS (per VM0033 depth) ──────────────────────────
  tar_target(tl_models, train_tl(tl_data, cfg_tl)),

  # ── SPATIAL PREDICTION RASTERS (4 bands per depth) ───────────────────────
  tar_terra_rast(tl_rasters, predict_tl_rasters(tl_models, covar_file_tl)),

  # ── MAPS & VALIDATION SUMMARY ────────────────────────────────────────────
  tar_target(tl_maps, plot_tl_maps(tl_rasters, tl_models, cfg_tl)),

  # ── REPORT ───────────────────────────────────────────────────────────────
  tar_quarto(report_tl, path = "reports/step4_transfer_learning.qmd")
)
