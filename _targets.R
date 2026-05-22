library(targets)
library(tarchetypes)
library(geotargets)

tar_option_set(
  packages = c("dplyr", "readr", "tidyr", "ggplot2", "sf", "terra", "randomForest")
)

tar_source("R/")

list(
  # ── INPUT FILE TRACKING ───────────────────────────────────────────────────
  tar_target(locations_file, "Pre-Analysis Data Preparation/data_raw/core_locations.csv", format = "file"),
  tar_target(samples_file,   "Pre-Analysis Data Preparation/data_raw/core_samples.csv",   format = "file"),
  tar_target(config_file,    "blue_carbon_config.R",                                       format = "file"),

  # ── CONFIGURATION ─────────────────────────────────────────────────────────
  tar_target(cfg, load_config(config_file)),

  # Covariate raster path — NULL when not configured or file is absent.
  # Downstream RF targets check for NULL and skip gracefully so the
  # non-spatial pipeline (Steps 1–2 + report_nonspatial) always completes.
  tar_target(covar_file, {
    path <- cfg$COVARIATE_RASTER
    if (is.null(path) || nchar(trimws(path)) == 0) {
      message("[covar_file] No raster path configured — spatial analysis will be skipped.")
      return(NULL)
    }
    if (!file.exists(path)) {
      message("[covar_file] Raster not found at: ", path, " — spatial analysis will be skipped.")
      return(NULL)
    }
    path
  }),

  # ── STEP 1: LOAD, PREPARE, HARMONIZE ──────────────────────────────────────
  tar_target(cores_raw,        load_raw_data(locations_file, samples_file, cfg)),
  tar_target(eda_plots,        run_eda(cores_raw, cfg)),
  tar_target(cores_harmonized, harmonize_depths(cores_raw, cfg)),
  tar_target(stratum_summary,  summarise_strata(cores_harmonized)),

  # ── STEP 2: SIMPLE EXTRAPOLATION ──────────────────────────────────────────
  tar_target(step2_extrapolation, simple_extrapolation(stratum_summary, cfg)),

  # ── STEP 3: RANDOM FOREST ─────────────────────────────────────────────────
  # Each target returns NULL early when covar_file is NULL, so the full
  # pipeline succeeds even without a raster. tar_terra_rast is replaced with
  # tar_target here so that a NULL return is serialisable.
  tar_target(rf_data, {
    if (is.null(covar_file)) return(NULL)
    prepare_rf_data(cores_harmonized, covar_file)
  }),
  tar_target(rf_models, {
    if (is.null(rf_data)) return(NULL)
    train_rf(rf_data)
  }),
  tar_target(rf_rasters, {
    if (is.null(rf_models)) return(NULL)
    predict_rf_rasters(rf_models, covar_file)
  }),
  tar_target(rf_importance_plot, {
    if (is.null(rf_models)) return(NULL)
    plot_rf_importance(rf_models, cfg)
  }),
  tar_target(rf_maps, {
    if (is.null(rf_rasters)) return(NULL)
    plot_rf_maps(rf_rasters, cfg)
  }),

  # ── REPORTS ───────────────────────────────────────────────────────────────
  tar_quarto(report_nonspatial, path = "reports/step1_nonspatial.qmd"),
  tar_quarto(report_rf,         path = "reports/step3_random_forest.qmd")
)
