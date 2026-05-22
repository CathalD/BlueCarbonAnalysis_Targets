library(targets)
library(tarchetypes)
library(geotargets)

tar_option_set(
  packages = c("dplyr", "readr", "tidyr", "ggplot2", "sf", "terra", "randomForest")
)

tar_source("R/")

# Usage:
#   Main pipeline first: targets::tar_make()
#   Then RF pipeline:    targets::tar_make(script = "_targets_rf.R", store = "_targets_rf")
#   Or use shortcut:     tmrf()

list(
  # ── CONFIGURATION ─────────────────────────────────────────────────────────
  tar_target(config_file_rf, "blue_carbon_config.R", format = "file"),
  tar_target(cfg_rf, load_config(config_file_rf)),

  # covar_file_rf path comes from config, so it depends on cfg_rf
  tar_target(covar_file_rf, cfg_rf$COVARIATE_RASTER, format = "file"),

  # ── LOCAL DATA FROM MAIN PIPELINE ─────────────────────────────────────────
  # Reads cores_harmonized from the main pipeline store (_targets/).
  # Re-run this pipeline manually after updating the main pipeline.
  tar_target(cores_harmonized_rf,
             targets::tar_read(cores_harmonized, store = "_targets")),

  # ── STEP 3: RANDOM FOREST ─────────────────────────────────────────────────
  tar_target(rf_data,            prepare_rf_data(cores_harmonized_rf, covar_file_rf)),
  tar_target(rf_models,          train_rf(rf_data)),
  tar_terra_rast(rf_rasters,     predict_rf_rasters(rf_models, covar_file_rf)),
  tar_target(rf_importance_plot, plot_rf_importance(rf_models, cfg_rf)),
  tar_target(rf_maps,            plot_rf_maps(rf_rasters, cfg_rf)),

  # ── REPORT ────────────────────────────────────────────────────────────────
  tar_quarto(report_rf, path = "reports/step3_random_forest.qmd")
)
