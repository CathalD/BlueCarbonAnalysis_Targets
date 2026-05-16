# ============================================================================
# STEP 4-01: TRANSFER LEARNING FOR BLUE CARBON — 3-MODEL COMPARISON
# ============================================================================
# PURPOSE: Leverage global blue carbon data (Janousek/WoSIS) to improve local
#          carbon stock predictions using three complementary approaches:
#
#   MODEL 1 — GEE Covariates:
#     Uses spectral/SAR/topographic covariates extracted via Google Earth Engine
#     (NDVI, LSWI, mNDWI, VV, elevation, etc.) as bridge variables between
#     global and local cores. Standard Wadoux-weighted domain adaptation.
#
#   MODEL 2 — Landscape Embeddings:
#     Uses deep-learning embedding similarity scores (5 cluster affinities)
#     derived from satellite imagery as the bridge variable set. These
#     capture holistic landscape character rather than individual band values.
#
#   MODEL 3 — Combined (GEE + Embeddings):
#     Concatenates both covariate sets for maximum information. Best when
#     embeddings and GEE covariates capture complementary aspects of site
#     similarity.
#
# METHODOLOGY (all 3 models):
#   1. Wadoux Instance Weighting — domain classifier identifies which global
#      cores are most similar to local environment
#   2. Weighted Global RF Prior — trained on similarity-weighted global data
#   3. Simple Bias Correction — mean residual from global→local prediction
#      (robust with tiny N; avoids overfitting a second RF)
#   4. Leave-One-Core-Out CV — spatially independent validation
#   5. Bootstrap Uncertainty — 500 resamplings of local residuals for PI
#   6. Spatial Predictions — 4-band rasters per depth if covariate rasters
#      are available
#
# PREREQUISITES:
#   1. Run P4_3b_Depth_Harmonization_Global.R to produce:
#        data_processed/global_cores_harmonized_VM0033.csv
#   2. Run S1_01_data_prep.R to produce:
#        data_processed/cores_harmonized_bluecarbon.csv
#   3. (For Models 1+3 spatial predictions) Rename covariate TIF to match
#        COVARIATE_RASTER path in blue_carbon_config.R
#   4. (For Models 2+3 spatial predictions) Place embedding similarity
#        GeoTIFF at path EMBEDDING_RASTER in blue_carbon_config.R
#
# INPUTS:
#   - data_processed/global_cores_harmonized_VM0033.csv
#       (from P4_3b; columns: profile_id, depth_cm_midpoint, carbon_stock_kg_m2,
#        latitude, longitude, + GEE covariates)
#   - Pre-Analysis Data Preparation/data_raw/CorePoints_EmbeddingSimilarity_BC_Canada.csv
#       (columns: profile_id, latitude, longitude, embedding_sim_to_cluster_0..4,
#        embedding_max_sim, embedding_best_cluster)
#   - data_processed/cores_harmonized_bluecarbon.csv  (local field cores)
#   - COVARIATE_RASTER (optional, for spatial predictions)
#   - EMBEDDING_RASTER  (optional, for Models 2+3 spatial predictions)
#
# OUTPUTS (all under outputs/Transfer_learning/):
#   ├── Model1_GEE/
#   │   ├── tables/   transfer_validation_gee.csv, variable_importance_gee_depth_*.csv
#   │   ├── models/   transfer_model_gee_depth_*.rds
#   │   └── maps/     depth_*_predictions_gee.tif  (4-band: Global_Prior, Local_Only,
#   │                                                Transfer_Final, Difference)
#   ├── Model2_Embeddings/
#   │   ├── tables/   transfer_validation_embeddings.csv
#   │   ├── models/   transfer_model_embeddings_depth_*.rds
#   │   └── maps/     depth_*_predictions_embeddings.tif
#   ├── Model3_Combined/
#   │   ├── tables/   transfer_validation_combined.csv
#   │   ├── models/   transfer_model_combined_depth_*.rds
#   │   └── maps/     depth_*_predictions_combined.tif
#   └── Comparison/
#       └── model_comparison_all_depths.csv
# ============================================================================

# ── PATH RESOLVER ────────────────────────────────────────────────────────────
# Script is 3 levels deep: BlueCarbon_Workflow_V1.0/Analysis_Workflow/
#   Transfer Learning/S4_01_transfer_learning.R
local({
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) {
    setwd(dirname(dirname(dirname(this_file))))
    return()
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NULL)
    if (!is.null(active) && nchar(active) > 0) {
      setwd(dirname(dirname(dirname(normalizePath(active)))))
      return()
    }
  }
  message("⚠ Could not auto-detect project root. Please run:\n",
          "  setwd('/path/to/BlueCarbon_Workflow_V1.0')\nbefore sourcing this script.")
})
# ─────────────────────────────────────────────────────────────────────────────

if (file.exists("Analysis_Workflow/blue_carbon_config.R")) {
  source("Analysis_Workflow/blue_carbon_config.R")
} else {
  stop("Configuration file not found. Ensure blue_carbon_config.R exists.")
}

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  STEP 4-01 — TRANSFER LEARNING (3-MODEL COMPARISON)     ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

# ── Logging ──────────────────────────────────────────────────────────────────
log_dir  <- "logs"
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
log_file <- file.path(log_dir, paste0("transfer_learning_", Sys.Date(), ".log"))

log_msg <- function(msg, level = "INFO") {
  entry <- sprintf("[%s] %s: %s", format(Sys.time(), "%H:%M:%S"), level, msg)
  cat(entry, "\n")
  cat(entry, "\n", file = log_file, append = TRUE)
}

log_msg("=== STEP 4-01: TRANSFER LEARNING STARTED ===")

# ── Packages ─────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(ranger)
  library(terra)
  library(sf)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# ── Output directories ────────────────────────────────────────────────────────
tl_root <- "outputs/Transfer_learning"
model_dirs <- list(
  gee       = file.path(tl_root, "Model1_GEE"),
  embed     = file.path(tl_root, "Model2_Embeddings"),
  combined  = file.path(tl_root, "Model3_Combined"),
  compare   = file.path(tl_root, "Comparison")
)

for (md in model_dirs) {
  for (sub in c("tables", "models", "maps")) {
    dir.create(file.path(md, sub), recursive = TRUE, showWarnings = FALSE)
  }
}
dir.create(model_dirs$compare, recursive = TRUE, showWarnings = FALSE)

# ── Configuration ─────────────────────────────────────────────────────────────
SEED               <- 42
N_BOOT             <- 500
N_THRESHOLD_REDUCED <- 15
VM0033_DEPTHS      <- c(7.5, 22.5, 40, 75)

set.seed(SEED)

# ── Covariate sets ────────────────────────────────────────────────────────────

# GEE bridge variables (spectral + SAR + topographic)
# Reduced set used when N_local < N_THRESHOLD_REDUCED
GEE_VARS_FULL <- c(
  "NDVI_median", "LSWI_median", "mNDWI_median", "SAVI_median",
  "EVI_median",  "NDBI_median",
  "B", "G", "R", "NIR", "SWIR1", "SWIR2",
  "brightness", "greenness", "wetness",
  "VV_mean", "VH_mean", "VVVH_ratio",
  "elevation_m", "elevationRelMHW", "slope"
)

GEE_VARS_REDUCED <- c(
  "NDVI_median", "LSWI_median", "mNDWI_median",
  "VV_mean", "elevation_m", "elevationRelMHW"
)

# Embedding bridge variables (cluster similarity scores)
EMBED_VARS <- c(
  "embedding_sim_to_cluster_0",
  "embedding_sim_to_cluster_1",
  "embedding_sim_to_cluster_2",
  "embedding_sim_to_cluster_3",
  "embedding_sim_to_cluster_4"
)

# ============================================================================
# SECTION 1: LOAD AND VALIDATE ALL INPUT DATA
# ============================================================================

cat("\n─────────────────────────────────────────────────\n")
cat(" SECTION 1: LOADING DATA\n")
cat("─────────────────────────────────────────────────\n\n")

# ── 1a. Global harmonized data (from P4_3b) ──────────────────────────────────
global_harmonized_path <- "data_processed/global_cores_harmonized_VM0033.csv"
if (!file.exists(global_harmonized_path)) {
  stop(sprintf(
    "Global harmonized data not found: %s\n  → Run P4_3b_Depth_Harmonization_Global.R first.",
    global_harmonized_path
  ))
}

global_harm <- read_csv(global_harmonized_path, show_col_types = FALSE) %>%
  mutate(data_source = "global")

log_msg(sprintf("Global harmonized data: %d rows, %d profiles",
                nrow(global_harm),
                n_distinct(global_harm$profile_id %||% global_harm$core_id)))

# ── 1b. GEE covariates for global cores ──────────────────────────────────────
gee_path <- file.path(DATA_RAW_DIR, "CorePoints_Covariates_BC_Canada.csv")
if (!file.exists(gee_path)) {
  stop(sprintf(
    "GEE covariate CSV not found: %s\n  → Ensure CorePoints_Covariates_BC_Canada.csv is in %s",
    gee_path, DATA_RAW_DIR
  ))
}

gee_covs <- read_csv(gee_path, show_col_types = FALSE) %>%
  mutate(profile_id = as.character(profile_id))

log_msg(sprintf("GEE covariate CSV: %d profiles", nrow(gee_covs)))

# ── 1c. Embedding similarity scores for global cores ─────────────────────────
embed_path <- file.path(DATA_RAW_DIR, "CorePoints_EmbeddingSimilarity_BC_Canada.csv")
if (!file.exists(embed_path)) {
  stop(sprintf(
    "Embedding similarity CSV not found: %s\n  → Ensure CorePoints_EmbeddingSimilarity_BC_Canada.csv is in %s",
    embed_path, DATA_RAW_DIR
  ))
}

embed_covs <- read_csv(embed_path, show_col_types = FALSE) %>%
  mutate(profile_id = as.character(profile_id))

log_msg(sprintf("Embedding similarity CSV: %d profiles", nrow(embed_covs)))

# ── 1d. Local harmonized data ─────────────────────────────────────────────────
local_path <- "data_processed/cores_harmonized_bluecarbon.csv"
if (!file.exists(local_path)) {
  stop(sprintf("Local harmonized data not found: %s\n  → Run S1_01_data_prep.R first.", local_path))
}

local_df <- read_csv(local_path, show_col_types = FALSE) %>%
  mutate(data_source = "local")

n_local_cores <- n_distinct(local_df$core_id)
log_msg(sprintf("Local harmonized data: %d rows, %d cores", nrow(local_df), n_local_cores))

# ── 1e. Join covariates to global harmonized data ────────────────────────────
# global_harm needs profile_id as character for join
global_harm <- global_harm %>%
  mutate(profile_id = as.character(
    if ("profile_id" %in% names(.)) profile_id else core_id
  ))

# Join GEE covariates (profile-level → depth-level: same covariates at all depths)
global_gee <- global_harm %>%
  left_join(
    gee_covs %>% select(profile_id, any_of(GEE_VARS_FULL)),
    by = "profile_id"
  )

# Trim GEE_VARS to those actually present after join
GEE_VARS_FULL    <- intersect(GEE_VARS_FULL,    names(global_gee))
GEE_VARS_REDUCED <- intersect(GEE_VARS_REDUCED, names(global_gee))

log_msg(sprintf("GEE vars available after join: %d full, %d reduced",
                length(GEE_VARS_FULL), length(GEE_VARS_REDUCED)))

if (length(GEE_VARS_REDUCED) == 0) {
  stop("No GEE bridge variables found after joining covariates to global data. Check profile_id join.")
}

# Join embedding covariates
global_embed <- global_harm %>%
  left_join(
    embed_covs %>% select(profile_id, any_of(EMBED_VARS)),
    by = "profile_id"
  )

EMBED_VARS <- intersect(EMBED_VARS, names(global_embed))
log_msg(sprintf("Embedding vars available after join: %d", length(EMBED_VARS)))

if (length(EMBED_VARS) == 0) {
  stop("No embedding variables found after joining. Check profile_id join between embedding CSV and harmonized global data.")
}

# Combined: join both to one global dataframe
global_combined <- global_gee %>%
  left_join(
    embed_covs %>% select(profile_id, any_of(EMBED_VARS)),
    by = "profile_id"
  )

# ── 1f. Check covariate rasters (optional) ────────────────────────────────────
gee_raster_available    <- FALSE
embed_raster_available  <- FALSE
gee_stack    <- NULL
embed_stack  <- NULL

# GEE covariate raster (same as used by RF step)
if (exists("COVARIATE_RASTER") && !is.null(COVARIATE_RASTER) && file.exists(COVARIATE_RASTER)) {
  tryCatch({
    gee_stack <- terra::rast(COVARIATE_RASTER)
    GEE_VARS_FULL    <- intersect(GEE_VARS_FULL,    names(gee_stack))
    GEE_VARS_REDUCED <- intersect(GEE_VARS_REDUCED, names(gee_stack))
    gee_raster_available <- TRUE
    log_msg(sprintf("GEE covariate raster loaded: %d bands available", nlyr(gee_stack)))
  }, error = function(e) {
    log_msg(sprintf("Could not load GEE covariate raster: %s — spatial predictions will be skipped for Models 1+3", e$message), "WARNING")
  })
} else {
  log_msg("COVARIATE_RASTER not set or file not found — spatial predictions skipped for Models 1+3", "WARNING")
}

# Embedding GeoTIFF (user-provided 5-band raster of embedding similarity scores)
if (exists("EMBEDDING_RASTER") && !is.null(EMBEDDING_RASTER) && file.exists(EMBEDDING_RASTER)) {
  tryCatch({
    embed_stack <- terra::rast(EMBEDDING_RASTER)
    # Ensure bands are named to match EMBED_VARS; try to rename if not
    if (!all(EMBED_VARS %in% names(embed_stack))) {
      n_bands <- nlyr(embed_stack)
      n_expected <- length(EMBED_VARS)
      if (n_bands == n_expected) {
        names(embed_stack) <- EMBED_VARS
        log_msg("Renamed embedding raster bands to match EMBED_VARS")
      } else {
        log_msg(sprintf(
          "Embedding raster has %d bands but %d expected (%s). Spatial predictions skipped for Models 2+3.",
          n_bands, n_expected, paste(EMBED_VARS, collapse = ", ")
        ), "WARNING")
        embed_stack <- NULL
      }
    }
    if (!is.null(embed_stack)) {
      embed_raster_available <- TRUE
      log_msg(sprintf("Embedding raster loaded: %d bands", nlyr(embed_stack)))
    }
  }, error = function(e) {
    log_msg(sprintf("Could not load embedding raster: %s — spatial predictions will be skipped for Models 2+3", e$message), "WARNING")
  })
} else {
  log_msg("EMBEDDING_RASTER not set or file not found — spatial predictions skipped for Models 2+3", "WARNING")
}

# ── 1g. Extract local covariate values from rasters ───────────────────────────
# Build a single local_covs dataframe with GEE + embedding values for each local core

local_coords <- local_df %>%
  select(core_id, longitude, latitude) %>%
  distinct()

local_gee_covs    <- local_coords
local_embed_covs  <- local_coords

if (gee_raster_available && !is.null(gee_stack)) {
  tryCatch({
    local_vect <- terra::vect(local_coords, geom = c("longitude", "latitude"), crs = "EPSG:4326")
    local_vect <- terra::project(local_vect, terra::crs(gee_stack))
    needed_gee <- intersect(c(GEE_VARS_FULL, GEE_VARS_REDUCED), names(gee_stack))
    extracted  <- terra::extract(gee_stack[[needed_gee]], local_vect, ID = FALSE)
    local_gee_covs <- bind_cols(local_coords, extracted)
    log_msg(sprintf("Extracted %d GEE covariates for %d local cores from raster",
                    ncol(extracted), nrow(local_coords)))
  }, error = function(e) {
    log_msg(sprintf("GEE raster extraction failed: %s — falling back to nearest-neighbour", e$message), "WARNING")
    gee_raster_available <<- FALSE
  })
}

if (!gee_raster_available) {
  # Fallback: nearest-neighbour from global GEE CSV
  log_msg("Using nearest-neighbour matching to assign GEE covariates to local cores")
  global_gee_pts <- gee_covs %>%
    select(profile_id, latitude, longitude, any_of(GEE_VARS_FULL)) %>%
    drop_na(latitude, longitude)

  nn_gee <- local_coords %>%
    rowwise() %>%
    mutate({
      dists  <- sqrt((global_gee_pts$longitude - longitude)^2 +
                       (global_gee_pts$latitude  - latitude)^2)
      nn_idx <- which.min(dists)
      nn_row <- global_gee_pts[nn_idx, ] %>% select(any_of(GEE_VARS_FULL))
      nn_row
    }) %>%
    ungroup()

  local_gee_covs <- bind_cols(local_coords, nn_gee %>% select(any_of(GEE_VARS_FULL)))
}

if (embed_raster_available && !is.null(embed_stack)) {
  tryCatch({
    local_vect_e <- terra::vect(local_coords, geom = c("longitude", "latitude"), crs = "EPSG:4326")
    local_vect_e <- terra::project(local_vect_e, terra::crs(embed_stack))
    extracted_e  <- terra::extract(embed_stack[[EMBED_VARS]], local_vect_e, ID = FALSE)
    local_embed_covs <- bind_cols(local_coords, extracted_e)
    log_msg(sprintf("Extracted %d embedding scores for %d local cores from raster",
                    ncol(extracted_e), nrow(local_coords)))
  }, error = function(e) {
    log_msg(sprintf("Embedding raster extraction failed: %s — falling back to nearest-neighbour", e$message), "WARNING")
    embed_raster_available <<- FALSE
  })
}

if (!embed_raster_available) {
  # Fallback: nearest-neighbour from embedding CSV
  log_msg("Using nearest-neighbour matching to assign embedding scores to local cores")
  embed_pts <- embed_covs %>%
    select(profile_id, latitude, longitude, any_of(EMBED_VARS)) %>%
    drop_na(latitude, longitude)

  nn_embed <- local_coords %>%
    rowwise() %>%
    mutate({
      dists  <- sqrt((embed_pts$longitude - longitude)^2 +
                       (embed_pts$latitude  - latitude)^2)
      nn_idx <- which.min(dists)
      nn_row <- embed_pts[nn_idx, ] %>% select(any_of(EMBED_VARS))
      nn_row
    }) %>%
    ungroup()

  local_embed_covs <- bind_cols(local_coords, nn_embed %>% select(any_of(EMBED_VARS)))
}

# Join all covariate sets to local depth data
local_model_gee <- local_df %>%
  left_join(local_gee_covs %>% select(core_id, any_of(GEE_VARS_FULL)), by = "core_id")

local_model_embed <- local_df %>%
  left_join(local_embed_covs %>% select(core_id, any_of(EMBED_VARS)), by = "core_id")

local_model_combined <- local_df %>%
  left_join(local_gee_covs   %>% select(core_id, any_of(GEE_VARS_FULL)), by = "core_id") %>%
  left_join(local_embed_covs %>% select(core_id, any_of(EMBED_VARS)),    by = "core_id")

cat(sprintf("Local data ready: %d cores × %d depths\n", n_local_cores, length(VM0033_DEPTHS)))
cat(sprintf("  Model 1 (GEE):      %d bridge vars available\n", length(GEE_VARS_FULL)))
cat(sprintf("  Model 2 (Embed):    %d bridge vars available\n", length(EMBED_VARS)))
cat(sprintf("  Model 3 (Combined): %d bridge vars available\n",
            length(GEE_VARS_FULL) + length(EMBED_VARS)))

# ============================================================================
# SECTION 2: TRANSFER LEARNING ENGINE
# ============================================================================
# Shared function — runs all 5 stages (Wadoux weighting, global prior,
# bias correction, LOOCV, bootstrap) for one model type across all depths.
# Returns a list: results_df, models_list (by depth), raster_predict_fn
# ============================================================================

run_transfer_learning <- function(
  model_name,          # e.g. "GEE", "Embeddings", "Combined"
  global_df,           # global data with carbon_stock_kg_m2 + covariates
  local_df,            # local data with carbon_stock_kg_m2 + covariates
  bridge_vars_full,    # full covariate vector
  bridge_vars_reduced, # reduced covariate vector (small N)
  out_dir              # output root (model_dirs$gee etc.)
) {

  cat(sprintf("\n═══════════════════════════════════════════════════\n"))
  cat(sprintf(" MODEL: %s\n", model_name))
  cat(sprintf("═══════════════════════════════════════════════════\n\n"))
  log_msg(sprintf("=== Starting Transfer Learning: %s ===", model_name))

  # ── Stage 0: Wadoux domain-classifier weights ─────────────────────────────
  # Train once using the reduced var set (robust with small N domain classifier)
  vars_for_weighting <- bridge_vars_reduced

  global_x <- global_df %>%
    select(all_of(vars_for_weighting)) %>%
    drop_na() %>%
    mutate(is_target = 0L)

  local_x <- local_df %>%
    select(all_of(vars_for_weighting)) %>%
    drop_na() %>%
    mutate(is_target = 1L)

  n_global_w <- nrow(global_x)
  n_local_w  <- nrow(local_x)

  cat(sprintf("  Domain classifier: %d global + %d local observations\n",
              n_global_w, n_local_w))

  if (n_local_w < 2) {
    log_msg(sprintf("%s: insufficient local data for domain classifier (N=%d)", model_name, n_local_w), "WARNING")
    return(NULL)
  }

  combined_x <- bind_rows(global_x, local_x)

  set.seed(SEED)
  rf_domain <- tryCatch(
    ranger(
      is_target ~ .,
      data          = combined_x,
      num.trees     = 500,
      probability   = TRUE,
      min.node.size = 5,
      seed          = SEED
    ),
    error = function(e) {
      log_msg(sprintf("Domain classifier failed: %s", e$message), "WARNING")
      NULL
    }
  )

  if (is.null(rf_domain)) {
    # Fallback: equal weights
    log_msg("Using equal weights (domain classifier failed)", "WARNING")
    global_df_w <- global_df %>% drop_na(all_of(vars_for_weighting)) %>% mutate(wadoux_weight = 1.0)
  } else {
    pred_global  <- predict(rf_domain, data = global_x)$predictions
    p_target <- if (is.matrix(pred_global)) {
      col_idx <- which(colnames(pred_global) == "1")
      if (length(col_idx) == 0) col_idx <- 2
      pred_global[, col_idx]
    } else pred_global

    p_target       <- pmin(pmax(p_target, 0.01), 0.99)
    wadoux_weights <- p_target / (1 - p_target)
    wadoux_weights <- wadoux_weights / mean(wadoux_weights)

    cat(sprintf("  Wadoux weights: min=%.2f  median=%.2f  max=%.2f\n",
                min(wadoux_weights), median(wadoux_weights), max(wadoux_weights)))
    eff_n <- sum(wadoux_weights)^2 / sum(wadoux_weights^2)
    cat(sprintf("  Effective sample size: %.0f / %d (%.1f%%)\n",
                eff_n, length(wadoux_weights), 100 * eff_n / length(wadoux_weights)))

    global_df_w <- global_df %>%
      drop_na(all_of(vars_for_weighting)) %>%
      mutate(wadoux_weight = wadoux_weights)

    write_csv(
      data.frame(weight = wadoux_weights, p_target = p_target),
      file.path(out_dir, "tables", "wadoux_weights.csv")
    )
  }

  # ── Stage 1–5: Per-depth loop ─────────────────────────────────────────────
  results_list <- list()
  models_list  <- list()

  for (d in VM0033_DEPTHS) {
    cat(sprintf("\n  ── Depth: %.1f cm ──\n", d))

    g_data <- global_df_w %>%
      filter(depth_cm_midpoint == d) %>%
      drop_na(carbon_stock_kg_m2, all_of(bridge_vars_reduced))

    l_data <- local_df %>%
      filter(depth_cm_midpoint == d) %>%
      drop_na(carbon_stock_kg_m2, all_of(bridge_vars_reduced))

    n_global <- nrow(g_data)
    n_local  <- nrow(l_data)
    n_cores  <- n_distinct(l_data$core_id)

    cat(sprintf("    Global N: %d | Local N: %d (%d cores)\n", n_global, n_local, n_cores))

    if (n_local < 2) {
      cat(sprintf("    SKIPPING depth %.1f — insufficient local data (N=%d)\n", d, n_local))
      next
    }
    if (n_global < 10) {
      cat(sprintf("    WARNING: very few global observations at this depth (N=%d)\n", n_global))
    }

    # Select covariate set
    covars <- if (n_local < N_THRESHOLD_REDUCED) {
      cat(sprintf("    Using REDUCED covariate set (%d vars)\n", length(bridge_vars_reduced)))
      bridge_vars_reduced
    } else {
      cat(sprintf("    Using FULL covariate set (%d vars)\n", length(bridge_vars_full)))
      bridge_vars_full
    }

    # ── STAGE 1: Weighted global RF prior ──────────────────────────────────
    formula_str <- paste("carbon_stock_kg_m2 ~", paste(covars, collapse = " + "))

    set.seed(SEED)
    rf_global <- tryCatch(
      ranger(
        formula      = as.formula(formula_str),
        data         = g_data,
        case.weights = g_data$wadoux_weight,
        num.trees    = 1000,
        mtry         = max(2, floor(length(covars) / 3)),
        min.node.size = 5,
        importance   = "permutation",
        seed         = SEED
      ),
      error = function(e) {
        log_msg(sprintf("  Global RF failed at depth %.1f: %s", d, e$message), "WARNING")
        NULL
      }
    )

    if (is.null(rf_global)) next

    var_imp <- data.frame(
      variable   = names(rf_global$variable.importance),
      importance = rf_global$variable.importance,
      depth_cm   = d,
      model      = model_name
    ) %>% arrange(desc(importance))

    cat(sprintf("    Top 3 vars: %s\n", paste(head(var_imp$variable, 3), collapse = ", ")))
    write_csv(var_imp,
              file.path(out_dir, "tables",
                        sprintf("variable_importance_depth_%.0f.csv", d * 10)))

    # ── STAGE 2: Local bias estimation ─────────────────────────────────────
    l_data$global_pred <- predict(rf_global, data = l_data)$predictions
    l_data$residual    <- l_data$carbon_stock_kg_m2 - l_data$global_pred

    bias_mean <- mean(l_data$residual)
    bias_sd   <- sd(l_data$residual)

    cat(sprintf("    Global→Local bias: %.3f ± %.3f kg/m²\n", bias_mean, bias_sd))

    # ── STAGE 3: Leave-one-core-out CV ─────────────────────────────────────
    unique_cores <- unique(l_data$core_id)
    cv_results   <- data.frame()

    for (held_core in unique_cores) {
      train_d <- l_data[l_data$core_id != held_core, ]
      test_d  <- l_data[l_data$core_id == held_core, ]
      if (nrow(train_d) < 2) next
      cv_bias          <- mean(train_d$residual)
      test_d$cv_pred   <- test_d$global_pred + cv_bias
      cv_results <- bind_rows(cv_results,
                              test_d %>% select(core_id, depth_cm_midpoint,
                                                carbon_stock_kg_m2, global_pred,
                                                cv_pred, residual))
    }

    if (nrow(cv_results) >= 2) {
      ss_tot       <- sum((cv_results$carbon_stock_kg_m2 - mean(cv_results$carbon_stock_kg_m2))^2)
      ss_res_cv    <- sum((cv_results$carbon_stock_kg_m2 - cv_results$cv_pred)^2)
      ss_res_glob  <- sum((cv_results$carbon_stock_kg_m2 - cv_results$global_pred)^2)
      r2_cv        <- 1 - ss_res_cv   / ss_tot
      r2_global    <- 1 - ss_res_glob / ss_tot
      rmse_cv      <- sqrt(mean((cv_results$carbon_stock_kg_m2 - cv_results$cv_pred)^2))
      rmse_global  <- sqrt(mean((cv_results$carbon_stock_kg_m2 - cv_results$global_pred)^2))
      mae_cv       <- mean(abs(cv_results$carbon_stock_kg_m2 - cv_results$cv_pred))
      cat(sprintf("    CV (bias-corr): R²=%.3f  RMSE=%.3f kg/m²  MAE=%.3f\n",
                  r2_cv, rmse_cv, mae_cv))
      cat(sprintf("    CV (global-only): R²=%.3f  RMSE=%.3f kg/m²\n", r2_global, rmse_global))
    } else {
      r2_cv <- r2_global <- rmse_cv <- rmse_global <- mae_cv <- NA
      cat("    CV not possible (insufficient splits)\n")
    }

    # ── STAGE 4: Bootstrap uncertainty ─────────────────────────────────────
    set.seed(SEED)
    boot_biases <- numeric(N_BOOT)
    for (b in seq_len(N_BOOT)) {
      boot_idx       <- sample(seq_len(n_local), replace = TRUE)
      boot_biases[b] <- mean(l_data$residual[boot_idx])
    }
    bias_se       <- sd(boot_biases)
    residual_var  <- var(l_data$residual)
    pred_se_total <- sqrt(bias_se^2 + residual_var)

    cat(sprintf("    Bootstrap bias SE: %.3f  Residual SD: %.3f\n",
                bias_se, sqrt(residual_var)))

    # ── STAGE 5: Save model object ──────────────────────────────────────────
    model_obj <- list(
      model_name       = model_name,
      depth_cm         = d,
      global_rf        = rf_global,
      bias_correction  = bias_mean,
      bias_se          = bias_se,
      residual_sd      = sqrt(residual_var),
      pred_se_total    = pred_se_total,
      predictors       = covars,
      n_global         = n_global,
      n_local          = n_local,
      cv_r2            = r2_cv,
      cv_rmse          = rmse_cv,
      cv_r2_global     = r2_global,
      cv_rmse_global   = rmse_global,
      method           = sprintf("Wadoux_weighted_%s_plus_bias_correction", model_name),
      created          = Sys.time()
    )

    # Attach prediction function
    model_obj$predict_with_uncertainty <- local({
      .rf   <- rf_global
      .bias <- bias_mean
      .se   <- pred_se_total
      function(newdata, ci_level = 0.90) {
        pred      <- predict(.rf, data = newdata)$predictions
        pred_corr <- pred + .bias
        z         <- qnorm(1 - (1 - ci_level) / 2)
        data.frame(
          prediction = pred_corr,
          se         = .se,
          lower      = pred_corr - z * .se,
          upper      = pred_corr + z * .se,
          ci_level   = ci_level
        )
      }
    })

    model_path <- file.path(out_dir, "models",
                            sprintf("transfer_model_depth_%.0f.rds", d * 10))
    saveRDS(model_obj, model_path)

    models_list[[as.character(d)]] <- model_obj

    results_list[[as.character(d)]] <- data.frame(
      model            = model_name,
      depth_cm         = d,
      n_global         = n_global,
      n_local          = n_local,
      n_covariates     = length(covars),
      bias_correction  = bias_mean,
      bias_se          = bias_se,
      residual_sd      = sqrt(residual_var),
      pred_se_total    = pred_se_total,
      cv_r2            = r2_cv,
      cv_rmse          = rmse_cv,
      cv_mae           = mae_cv,
      cv_r2_global_only = r2_global,
      cv_rmse_global_only = rmse_global
    )
  } # end depth loop

  if (length(results_list) == 0) {
    log_msg(sprintf("%s: no depths were successfully modelled", model_name), "WARNING")
    return(NULL)
  }

  results_df <- bind_rows(results_list)
  write_csv(results_df,
            file.path(out_dir, "tables",
                      sprintf("transfer_validation_%s.csv",
                              tolower(gsub(" ", "_", model_name)))))

  log_msg(sprintf("%s: %d depths modelled successfully", model_name, nrow(results_df)))

  return(list(results = results_df, models = models_list))
}

# ============================================================================
# SECTION 3: RUN ALL THREE MODELS
# ============================================================================

cat("\n─────────────────────────────────────────────────\n")
cat(" SECTION 3: RUNNING TRANSFER LEARNING MODELS\n")
cat("─────────────────────────────────────────────────\n")

tl_gee      <- run_transfer_learning(
  model_name          = "GEE",
  global_df           = global_gee,
  local_df            = local_model_gee,
  bridge_vars_full    = GEE_VARS_FULL,
  bridge_vars_reduced = GEE_VARS_REDUCED,
  out_dir             = model_dirs$gee
)

tl_embed    <- run_transfer_learning(
  model_name          = "Embeddings",
  global_df           = global_embed,
  local_df            = local_model_embed,
  bridge_vars_full    = EMBED_VARS,
  bridge_vars_reduced = EMBED_VARS,     # same set regardless of N (only 5 vars)
  out_dir             = model_dirs$embed
)

tl_combined <- run_transfer_learning(
  model_name          = "Combined",
  global_df           = global_combined,
  local_df            = local_model_combined,
  bridge_vars_full    = c(GEE_VARS_FULL,    EMBED_VARS),
  bridge_vars_reduced = c(GEE_VARS_REDUCED, EMBED_VARS),
  out_dir             = model_dirs$combined
)

# ============================================================================
# SECTION 4: SPATIAL PREDICTION RASTERS
# ============================================================================
# For each model that has covariate rasters available, generate 4-band
# per-depth prediction rasters:
#   Band 1 — Global_Prior   (Wadoux-weighted RF, uncorrected)
#   Band 2 — Local_Only     (constant raster = local mean at that depth)
#   Band 3 — Transfer_Final (bias-corrected)
#   Band 4 — Difference     (Transfer_Final − Local_Only)
# ============================================================================

cat("\n─────────────────────────────────────────────────\n")
cat(" SECTION 4: SPATIAL PREDICTION RASTERS\n")
cat("─────────────────────────────────────────────────\n\n")

generate_prediction_rasters <- function(
  model_result,   # output of run_transfer_learning()
  covar_stack,    # terra SpatRaster with bridge vars as band names
  local_df,       # for computing local mean at each depth
  out_dir,
  model_name
) {
  if (is.null(model_result) || is.null(covar_stack)) return(invisible(NULL))

  models_list <- model_result$models

  for (depth_char in names(models_list)) {
    d   <- as.numeric(depth_char)
    obj <- models_list[[depth_char]]

    covars_needed <- obj$predictors
    missing       <- setdiff(covars_needed, names(covar_stack))

    if (length(missing) > 0) {
      log_msg(sprintf("%s depth %.1fcm: missing raster bands %s — skipping",
                      model_name, d, paste(missing, collapse = ", ")), "WARNING")
      next
    }

    cat(sprintf("  Generating rasters: %s depth %.1f cm\n", model_name, d))

    pred_stack <- covar_stack[[covars_needed]]

    # Band 1: Global prior (uncorrected)
    global_prior <- tryCatch(
      terra::predict(
        pred_stack,
        obj$global_rf,
        fun = function(model, newdata) predict(model, data = newdata)$predictions,
        na.rm = TRUE
      ),
      error = function(e) {
        log_msg(sprintf("  Raster prediction failed: %s", e$message), "WARNING")
        NULL
      }
    )
    if (is.null(global_prior)) next
    names(global_prior) <- "Global_Prior"

    # Band 2: Local-only constant
    local_mean  <- mean(local_df %>% filter(depth_cm_midpoint == d) %>%
                          pull(carbon_stock_kg_m2), na.rm = TRUE)
    local_only  <- global_prior * 0 + local_mean
    names(local_only) <- "Local_Only"

    # Band 3: Transfer Final (bias-corrected)
    transfer_final <- global_prior + obj$bias_correction
    names(transfer_final) <- "Transfer_Final"

    # Band 4: Difference
    difference <- transfer_final - local_only
    names(difference) <- "Difference"

    output_stack <- c(global_prior, local_only, transfer_final, difference)
    out_path     <- file.path(out_dir, "maps",
                              sprintf("depth_%.0f_predictions_%s.tif",
                                      d * 10,
                                      tolower(model_name)))

    terra::writeRaster(output_stack, out_path, overwrite = TRUE,
                       gdal = c("COMPRESS=LZW", "TILED=YES"))
    log_msg(sprintf("Saved: %s", basename(out_path)))

    cat(sprintf("    Global Prior mean: %.3f  Transfer Final mean: %.3f kg/m²\n",
                as.numeric(terra::global(global_prior,  "mean", na.rm = TRUE)),
                as.numeric(terra::global(transfer_final, "mean", na.rm = TRUE))))
  }
}

# Model 1: GEE rasters
if (gee_raster_available) {
  generate_prediction_rasters(tl_gee, gee_stack, local_model_gee,
                               model_dirs$gee, "GEE")
} else {
  cat("  Model 1 (GEE): spatial predictions skipped — no covariate raster\n")
}

# Model 2: Embedding rasters
if (embed_raster_available) {
  generate_prediction_rasters(tl_embed, embed_stack, local_model_embed,
                               model_dirs$embed, "Embeddings")
} else {
  cat("  Model 2 (Embeddings): spatial predictions skipped — no embedding raster\n")
}

# Model 3: Combined — needs both raster types
if (gee_raster_available && embed_raster_available) {
  combined_stack <- c(gee_stack[[intersect(c(GEE_VARS_FULL, GEE_VARS_REDUCED), names(gee_stack))]],
                      embed_stack[[EMBED_VARS]])
  generate_prediction_rasters(tl_combined, combined_stack, local_model_combined,
                               model_dirs$combined, "Combined")
} else {
  cat("  Model 3 (Combined): spatial predictions skipped — needs both GEE + embedding rasters\n")
}

# ============================================================================
# SECTION 5: CROSS-MODEL COMPARISON TABLE
# ============================================================================

cat("\n─────────────────────────────────────────────────\n")
cat(" SECTION 5: MODEL COMPARISON\n")
cat("─────────────────────────────────────────────────\n\n")

comparison_rows <- list()
for (res in list(tl_gee, tl_embed, tl_combined)) {
  if (!is.null(res)) comparison_rows <- c(comparison_rows, list(res$results))
}

if (length(comparison_rows) > 0) {
  comparison_df <- bind_rows(comparison_rows) %>%
    arrange(depth_cm, model)

  write_csv(comparison_df,
            file.path(model_dirs$compare, "model_comparison_all_depths.csv"))
  log_msg("Saved cross-model comparison table")

  # Console summary
  cat("Cross-model CV performance (Leave-One-Core-Out):\n\n")
  cat(sprintf("  %-12s  %-8s  %-8s  %-8s  %-8s\n",
              "Model", "Depth", "R²(TL)", "RMSE(TL)", "RMSE(Global)"))
  cat(sprintf("  %s\n", strrep("-", 54)))
  for (i in seq_len(nrow(comparison_df))) {
    r <- comparison_df[i, ]
    cat(sprintf("  %-12s  %-8.1f  %-8.3f  %-8.3f  %-8.3f\n",
                r$model, r$depth_cm,
                ifelse(is.na(r$cv_r2),   NA_real_, r$cv_r2),
                ifelse(is.na(r$cv_rmse), NA_real_, r$cv_rmse),
                ifelse(is.na(r$cv_rmse_global_only), NA_real_, r$cv_rmse_global_only)))
  }
  cat("\n")
} else {
  cat("⚠ No model results available for comparison.\n")
}

# ============================================================================
# SECTION 6: PRESENTATION COPY
# ============================================================================

cat("\n─────────────────────────────────────────────────\n")
cat(" SECTION 6: PRESENTATION COPY\n")
cat("─────────────────────────────────────────────────\n\n")

presentation_dir <- file.path(tl_root, "Global_vs_Local_Comparisons")
dir.create(presentation_dir, recursive = TRUE, showWarnings = FALSE)

# Copy comparison table
comp_table <- file.path(model_dirs$compare, "model_comparison_all_depths.csv")
if (file.exists(comp_table)) {
  file.copy(comp_table,
            file.path(presentation_dir, "Transfer_Learning_Model_Comparison.csv"),
            overwrite = TRUE)
}

# Copy per-model validation summaries
for (nm in c("GEE", "Embeddings", "Combined")) {
  src_dir <- switch(nm, GEE = model_dirs$gee, Embeddings = model_dirs$embed, Combined = model_dirs$combined)
  val_file <- file.path(src_dir, "tables",
                        sprintf("transfer_validation_%s.csv", tolower(nm)))
  if (file.exists(val_file)) {
    file.copy(val_file,
              file.path(presentation_dir,
                        sprintf("Transfer_Validation_%s.csv", nm)),
              overwrite = TRUE)
    cat(sprintf("  Copied: Transfer_Validation_%s.csv\n", nm))
  }
}

# Copy best prediction rasters (Transfer_Final band from deepest depth = 75cm)
for (nm in c("GEE", "Embeddings", "Combined")) {
  src_dir <- switch(nm, GEE = model_dirs$gee, Embeddings = model_dirs$embed, Combined = model_dirs$combined)
  tif_files <- list.files(file.path(src_dir, "maps"), pattern = "\\.tif$", full.names = TRUE)
  if (length(tif_files) > 0) {
    deepest_tif <- tif_files[which.max(as.numeric(
      gsub(".*depth_([0-9]+)_predictions.*", "\\1", basename(tif_files))
    ))]
    dst <- file.path(presentation_dir,
                     sprintf("Transfer_Learning_%s_Final_Map_0_to_100cm.tif", nm))
    file.copy(deepest_tif, dst, overwrite = TRUE)
    cat(sprintf("  Copied: Transfer_Learning_%s_Final_Map_0_to_100cm.tif\n", nm))
  }
}

# ============================================================================
# SECTION 7: FINAL CONSOLE SUMMARY
# ============================================================================

cat("\n╔══════════════════════════════════════════════════════════╗\n")
cat("║  STEP 4-01 COMPLETE — TRANSFER LEARNING                 ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

cat("Models trained:\n")
cat(sprintf("  Model 1 — GEE Covariates:     %s\n",
            if (!is.null(tl_gee))      sprintf("%d depths", nrow(tl_gee$results))      else "FAILED"))
cat(sprintf("  Model 2 — Embeddings:         %s\n",
            if (!is.null(tl_embed))    sprintf("%d depths", nrow(tl_embed$results))    else "FAILED"))
cat(sprintf("  Model 3 — Combined:           %s\n",
            if (!is.null(tl_combined)) sprintf("%d depths", nrow(tl_combined$results)) else "FAILED"))

cat("\nSpatial predictions:\n")
cat(sprintf("  GEE raster:       %s\n", if (gee_raster_available)   "✓ generated" else "⚠ skipped (no COVARIATE_RASTER)"))
cat(sprintf("  Embedding raster: %s\n", if (embed_raster_available) "✓ generated" else "⚠ skipped (no EMBEDDING_RASTER)"))

cat("\nOutputs:\n")
cat(sprintf("  %s\n", tl_root))
cat("  ├── Model1_GEE/\n")
cat("  ├── Model2_Embeddings/\n")
cat("  ├── Model3_Combined/\n")
cat("  └── Comparison/model_comparison_all_depths.csv\n\n")

cat("To use a model for predictions:\n")
cat("  m <- readRDS('outputs/Transfer_learning/Model1_GEE/models/transfer_model_depth_75.rds')\n")
cat("  p <- m$predict_with_uncertainty(newdata, ci_level = 0.90)\n\n")

log_msg("=== STEP 4-01: TRANSFER LEARNING COMPLETE ===")
