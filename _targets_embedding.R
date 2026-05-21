# =============================================================================
# _targets_embedding.R
# Embedding-weighted transfer learning pipeline (Model 2).
#
# Uses GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL (64-band foundation model
# embeddings) to compute cosine similarity between each global Janousek
# core and the local AOI, replacing the Wadoux RF domain classifier used
# in _targets_transfer.R.
#
# Run with:
#   targets::tar_make(script = "_targets_embedding.R",
#                     store  = "_targets_embedding")
#
# Prerequisites:
#   1. Main pipeline completed (_targets store exists):
#        targets::tar_make()
#   2. Preanalysis pipeline completed (global covariates CSV written):
#        targets::tar_make(script="_targets_preanalysis.R",
#                          store="_targets_preanalysis")
#   3. rgee authenticated with Drive enabled (run once):
#        library(rgee); ee_Initialize(user="your@email.com", drive=TRUE)
#
# Pipeline structure:
#   Phase 1 — GEE embedding extraction:
#     global_embeddings     : 64-d embedding at each global core location
#     aoi_embedding_raster  : 64-band embedding image over local AOI
#   Phase 2 — Similarity weights:
#     emb_weights           : cosine similarity → instance weights per core
#   Phase 3 — Embedding-weighted transfer learning:
#     emb_tl_data           : harmonized global + local training data
#     emb_tl_models         : per-depth weighted RF + bias correction
#     emb_tl_rasters        : 4-band GeoTIFF per depth (same format as TL)
#     emb_tl_maps           : comparison maps + similarity plots
# =============================================================================

library(targets)
library(tarchetypes)
library(geotargets)

tar_option_set(
  packages = c("dplyr", "readr", "tidyr", "ggplot2", "sf",
               "terra", "ranger", "rgee")
)

# Shared R modules (config, data_prep, depth_harmonization, transfer_learning,
# random_forest, etc.) plus embedding-specific modules.
tar_source("R/")

GEE_PROJECT <- "north-star-project-470316"
EMB_YEARS   <- 2023:2025   # years to average for the embedding

list(

  # ── Configuration ──────────────────────────────────────────────────────────
  tar_target(config_file_emb, "blue_carbon_config.R", format = "file"),
  tar_target(cfg_emb, load_config(config_file_emb)),
  tar_target(covar_file_emb, cfg_emb$COVARIATE_RASTER, format = "file"),

  # ── Local cores from main pipeline ─────────────────────────────────────────
  tar_target(
    cores_harmonized_emb,
    targets::tar_read(cores_harmonized, store = "_targets")
  ),

  # ── Global data ────────────────────────────────────────────────────────────
  tar_target(
    global_layers_file_emb,
    "Pre-Analysis Data Preparation/data_global/combined_layers_filtered.csv",
    format = "file"
  ),
  tar_target(
    global_harmonized_emb,
    harmonize_global_layers(global_layers_file_emb, cfg_emb)
  ),
  tar_target(
    global_covar_file_emb,
    "Pre-Analysis Data Preparation/data_raw/CorePoints_Covariates_BC_Canada.csv",
    format = "file"
  ),

  # ── Phase 1: GEE embedding extraction ──────────────────────────────────────

  # 64-d embedding vector at each of the ~952 global core locations.
  # Uses GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL averaged over EMB_YEARS.
  tar_target(
    global_embeddings,
    extract_global_embeddings(
      global_covar_file_emb,
      gee_project = GEE_PROJECT,
      years       = EMB_YEARS
    )
  ),

  # 64-band embedding raster covering the local AOI extent.
  # Requires Drive authentication (ee_Initialize with drive = TRUE).
  tar_terra_rast(
    aoi_embedding_raster,
    extract_aoi_embedding_raster(
      covar_file_emb,
      gee_project = GEE_PROJECT,
      years       = EMB_YEARS
    )
  ),

  # ── Phase 2: Cosine similarity → instance weights ──────────────────────────
  # For each global core: weight = (cosine_sim)^alpha, normalised so mean = 1.
  # alpha = 5 strongly favours cores whose embedding closely matches the AOI.
  tar_target(
    emb_weights,
    compute_embedding_weights(
      aoi_emb_rast  = aoi_embedding_raster,
      global_emb_df = global_embeddings,
      alpha         = 5
    )
  ),

  # ── Phase 3: Embedding-weighted transfer learning ──────────────────────────

  tar_target(
    emb_tl_data,
    prepare_emb_tl_data(
      cores_harmonized = cores_harmonized_emb,
      global_harmonized = global_harmonized_emb,
      global_covar_file = global_covar_file_emb,
      local_covar_file  = covar_file_emb
    )
  ),

  tar_target(
    emb_tl_models,
    train_emb_tl(emb_tl_data, emb_weights, cfg_emb)
  ),

  tar_terra_rast(
    emb_tl_rasters,
    predict_emb_tl_rasters(emb_tl_models, covar_file_emb)
  ),

  tar_target(
    emb_tl_maps,
    plot_emb_tl_maps(emb_tl_rasters, emb_tl_models, cfg_emb)
  ),

  # ── Report ─────────────────────────────────────────────────────────────────
  tar_quarto(report_emb_tl, path = "reports/step5_embedding_tl.qmd")

)
