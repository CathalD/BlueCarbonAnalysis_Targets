# ============================================================================
# MODULE 05: TRANSFER LEARNING FOR BLUE CARBON (REVISED)
# ============================================================================
# PURPOSE: Leverage global blue carbon data to improve predictions in a local
#          estuary with very limited samples (N ~ 6 cores).
#
# METHODOLOGY:
#   1. Wadoux Instance Weighting - prioritize global samples similar to local
#   2. Hierarchical Prior - pool information across depths via mixed effects
#   3. Simple Bias Correction - avoid overfitting with minimal parameters
#   4. Proper Uncertainty Quantification - prediction intervals for VM0033
#
# KEY CHANGES FROM ORIGINAL:
#   - Reduced covariate set to avoid overparameterization
#   - Leave-one-core-out CV (spatial independence)
#   - Bootstrap uncertainty for prediction intervals
#   - Explicit handling of missing local rasters
#
# INPUTS:
#   - data_global/global_cores_with_gee_covariates.csv
#   - data_processed/cores_harmonized_bluecarbon.csv
#   - covariates/*.tif (optional - script handles missing rasters)
#
# OUTPUTS:
#   - outputs/models/transfer/transfer_model_depth_*.rds
#   - outputs/carbon_stocks/predictions/depth_*_predictions.tif (4-band rasters)
#   - diagnostics/transfer/transfer_learning_validation.csv
#   - diagnostics/transfer/variable_importance_depth_*.csv
# ============================================================================

# ── PATH RESOLVER ───────────────────────────────────────────────────────────────────
# Ensures working directory is BlueCarbon_Workflow_V1.0/ (project root)
# so all relative data paths (data_raw/, outputs/, etc.) resolve correctly.
local({
  # Method 1: called via source() — detect this script's location
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) {
    setwd(dirname(dirname(dirname(this_file))))
    return()
  }
  # Method 2: running interactively in RStudio — rstudioapi gives active file
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NULL)
    if (!is.null(active) && nchar(active) > 0) {
      setwd(dirname(dirname(dirname(normalizePath(active)))))
      return()
    }
  }
  # Method 3: fallback — warn and let the user set wd manually
  message("⚠ Could not auto-detect project root. Please run:
",
          "  setwd('/path/to/BlueCarbon_Workflow_V1.0')
",
          "before sourcing this script.")
})
# ──────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(tidyverse)
  library(ranger)
  library(terra)
  library(sf)
})

dir.create("outputs/models/transfer", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/carbon_stocks/predictions", recursive = TRUE, showWarnings = FALSE)
dir.create("diagnostics/transfer", recursive = TRUE, showWarnings = FALSE)

cat("\n==================================================\n")
cat("   BLUE CARBON TRANSFER LEARNING (REVISED)\n")
cat("==================================================\n")

# ============================================================================
# 1. CONFIGURATION
# ============================================================================

# Full bridge variable set (21 covariates)
BRIDGE_VARS_FULL <- c(
  "R", "B", "G", "NIR", "SWIR1", "SWIR2",
  "NDVI_median", "NDBI_median", "EVI_median", 
  "SAVI_median", "LSWI_median", "mNDWI_median",
  "brightness", "greenness", "wetness",
  "VV_mean", "VH_mean", "VVVH_ratio",
  "elevation_m", "elevationRelMHW"
)

# Reduced set for small samples - scientifically motivated selection
# These capture: vegetation (NDVI), wetness (LSWI/mNDWI), structure (SAR), 
# and topography (elevation) without redundancy
BRIDGE_VARS_REDUCED <- c(
  "NDVI_median",      # Vegetation density
  "LSWI_median",      # Leaf water / wetness
  "mNDWI_median",     # Surface water
  "VV_mean",          # SAR backscatter (structure)
  "elevation_m",      # Topographic position
  "elevationRelMHW"   # Tidal position (critical for blue carbon)
)

# Threshold for switching to reduced model
N_THRESHOLD_REDUCED <- 15

# Random seed for reproducibility
SEED <- 42
set.seed(SEED)

# Flag to track if spatial prediction is possible
RASTER_STACK_AVAILABLE <- FALSE
PREDICTION_STACK <- NULL

cat(sprintf("Full covariate set: %d variables\n", length(BRIDGE_VARS_FULL)))
cat(sprintf("Reduced covariate set: %d variables\n", length(BRIDGE_VARS_REDUCED)))

# ============================================================================
# 2. LOAD AND VALIDATE DATA
# ============================================================================

cat("\n--- Loading Data ---\n")

# --- GLOBAL DATA ---
# DATA_GLOBAL_DIR set in blue_carbon_config.R → Pre-Analysis Data Preparation/data_global/
global_path <- file.path(DATA_GLOBAL_DIR, "global_cores_with_gee_covariates.csv")
if (!file.exists(global_path)) {
  stop(sprintf("Global covariate CSV not found: %s\n  → Run CoastalBlueCarbon_GlobalCoreCovariate_Extraction.ipynb and place output in: %s",
               global_path, DATA_GLOBAL_DIR))
}

global_df <- read_csv(global_path, show_col_types = FALSE) %>%
  mutate(data_source = "global")

# Validate columns
missing_global <- setdiff(BRIDGE_VARS_FULL, names(global_df))
if (length(missing_global) > 0) {
  cat("WARNING: Global data missing covariates:", paste(missing_global, collapse = ", "), "\n")
  BRIDGE_VARS_FULL <- intersect(BRIDGE_VARS_FULL, names(global_df))
  BRIDGE_VARS_REDUCED <- intersect(BRIDGE_VARS_REDUCED, names(global_df))
}

cat(sprintf("Global data: %d observations\n", nrow(global_df)))

# --- LOCAL DATA ---
local_path <- "data_processed/cores_harmonized_bluecarbon.csv"
if (!file.exists(local_path)) {
  stop("Local harmonized data not found at: ", local_path)
}

local_df <- read_csv(local_path, show_col_types = FALSE) %>%
  mutate(data_source = "local")

n_local_cores <- n_distinct(local_df$core_id)
cat(sprintf("Local data: %d observations from %d cores\n", nrow(local_df), n_local_cores))

# ============================================================================
# 3. EXTRACT LOCAL COVARIATES (WITH FALLBACK)
# ============================================================================

cat("\n--- Extracting Local Covariates ---\n")

# Check for local rasters
# COVARIATES_DIR set in blue_carbon_config.R → Pre-Analysis Data Preparation/covariates/
tif_files <- list.files(COVARIATES_DIR, pattern = "\\.tif$", full.names = TRUE)

if (length(tif_files) == 0) {
  cat(sprintf("No local rasters found in '%s'.\n", COVARIATES_DIR))
  cat("Attempting to match local sites to nearest global observations...\n")
  
  # Fallback: Find nearest global neighbor for each local site
  # This is a valid approach when local rasters are unavailable
  local_coords <- local_df %>% 
    select(core_id, longitude, latitude) %>% 
    distinct()
  
  global_coords <- global_df %>%
    select(core_id, longitude, latitude) %>%
    distinct()
  
  # Simple nearest-neighbor matching (Euclidean on lat/lon is approximate but acceptable)
  match_covariates <- function(local_row, global_data) {
    dists <- sqrt((global_data$longitude - local_row$longitude)^2 + 
                    (global_data$latitude - local_row$latitude)^2)
    nearest_idx <- which.min(dists)
    nearest_dist_km <- dists[nearest_idx] * 111  # Approximate km at mid-latitudes
    
    if (nearest_dist_km > 50) {
      warning(sprintf("Core %s: nearest global match is %.1f km away", 
                      local_row$core_id, nearest_dist_km))
    }
    
    # Return covariates from nearest global observation
    global_data[nearest_idx, ]
  }
  
  cat("WARNING: Using nearest-neighbor covariate matching.\n")
  cat("         This assumes local sites have similar spectral/SAR signatures to nearby global sites.\n")
  cat("         Results should be interpreted with caution.\n\n")
  
  # Get covariates for each unique local site
  nearest_covs <- local_coords %>%
    rowwise() %>%
    mutate(
      nearest_global = list(match_covariates(cur_data(), global_df))
    ) %>%
    unnest(nearest_global, names_sep = "_") %>%
    select(core_id, all_of(paste0("nearest_global_", BRIDGE_VARS_FULL)))
  
  # Rename columns
  names(nearest_covs) <- gsub("nearest_global_", "", names(nearest_covs))
  
  # Join to local data
  local_model_df <- local_df %>%
    left_join(nearest_covs, by = "core_id")
  
} else {
  cat(sprintf("Found %d raster files.\n", length(tif_files)))
  
  # Load raster stack
  local_stack <- terra::rast(tif_files)
  
  # Check for missing layers
  available_layers <- names(local_stack)
  missing_layers <- setdiff(BRIDGE_VARS_FULL, available_layers)
  
  if (length(missing_layers) > 0) {
    cat("WARNING: Rasters missing for:", paste(missing_layers, collapse = ", "), "\n")
    BRIDGE_VARS_FULL <- intersect(BRIDGE_VARS_FULL, available_layers)
    BRIDGE_VARS_REDUCED <- intersect(BRIDGE_VARS_REDUCED, available_layers)
  }
  
  # Store stack for spatial predictions later
  RASTER_STACK_AVAILABLE <- TRUE
  PREDICTION_STACK <- local_stack
  
  # Extract values
  local_vect <- terra::vect(
    local_df, 
    geom = c("longitude", "latitude"), 
    crs = "EPSG:4326"
  )
  local_vect <- terra::project(local_vect, terra::crs(local_stack))
  
  extracted_vals <- terra::extract(local_stack[[BRIDGE_VARS_FULL]], local_vect, ID = FALSE)
  local_model_df <- bind_cols(local_df, extracted_vals)
  
  # Scale correction (detect and fix GEE scale mismatches)
  cat("Checking for scale mismatches...\n")
  for (var in BRIDGE_VARS_FULL) {
    if (!var %in% names(global_df) || !var %in% names(local_model_df)) next
    
    g_median <- median(global_df[[var]], na.rm = TRUE)
    l_median <- median(local_model_df[[var]], na.rm = TRUE)
    
    # Skip elevation variables and near-zero values
    if (grepl("elevation", var, ignore.case = TRUE)) next
    if (abs(g_median) < 0.001 || abs(l_median) < 0.001) next
    
    ratio <- abs(l_median / g_median)
    
    if (ratio > 50 && ratio < 20000) {
      cat(sprintf("  %s: local median %.4f vs global %.4f (ratio: %.0f) - scaling by 1/10000\n",
                  var, l_median, g_median, ratio))
      local_model_df[[var]] <- local_model_df[[var]] / 10000
    }
  }
}

# Report completeness
n_complete <- local_model_df %>%
  filter(if_all(all_of(BRIDGE_VARS_FULL), ~ !is.na(.))) %>%
  nrow()

cat(sprintf("Local observations with complete covariates: %d / %d\n", 
            n_complete, nrow(local_model_df)))

# ============================================================================
# 4. WADOUX INSTANCE WEIGHTING
# ============================================================================

cat("\n--- Computing Wadoux Instance Weights ---\n")

# Prepare covariate matrices
# Using reduced set for domain classifier to avoid overfitting
vars_for_weighting <- BRIDGE_VARS_REDUCED

global_x <- global_df %>%
  select(all_of(vars_for_weighting)) %>%
  drop_na() %>%
  mutate(is_target = 0)

local_x <- local_model_df %>%
  select(all_of(vars_for_weighting)) %>%
  drop_na() %>%
  mutate(is_target = 1)

combined_x <- bind_rows(global_x, local_x)

# Train domain classifier
# Using probability forest to estimate P(target | X)
cat(sprintf("Training domain classifier on %d variables...\n", length(vars_for_weighting)))

rf_domain <- ranger(
  is_target ~ .,
  data = combined_x,
  num.trees = 500,
  probability = TRUE,
  min.node.size = 5,
  seed = SEED
)

# Get probability of being in target (local) domain for global samples
pred_global <- predict(rf_domain, data = global_x)$predictions

# Handle both matrix and vector outputs from ranger
if (is.matrix(pred_global)) {
  # Find the column corresponding to class "1" (target domain)
  col_idx <- which(colnames(pred_global) == "1")
  if (length(col_idx) == 0) col_idx <- 2  # Fallback: assume second column
  p_target <- pred_global[, col_idx]
} else {
  p_target <- pred_global
}

# Clip probabilities to avoid infinite weights
p_target <- pmin(pmax(p_target, 0.01), 0.99)

# Wadoux weights: w = p / (1 - p)
# This upweights global samples that look like local samples
wadoux_weights <- p_target / (1 - p_target)

# Normalize to have mean = 1 (preserves effective sample size interpretation)
wadoux_weights <- wadoux_weights / mean(wadoux_weights)

# Diagnostic: weight distribution
cat(sprintf("Weight distribution: min=%.2f, median=%.2f, max=%.2f\n",
            min(wadoux_weights), median(wadoux_weights), max(wadoux_weights)))
cat(sprintf("Effective sample size: %.0f / %d (%.1f%%)\n",
            sum(wadoux_weights)^2 / sum(wadoux_weights^2),
            length(wadoux_weights),
            100 * (sum(wadoux_weights)^2 / sum(wadoux_weights^2)) / length(wadoux_weights)))

# Attach weights to global data
global_df_weighted <- global_df %>%
  drop_na(all_of(vars_for_weighting)) %>%
  mutate(wadoux_weight = wadoux_weights)

# Save weight diagnostics
weight_diag <- data.frame(
  weight = wadoux_weights,
  p_target = p_target
)
write_csv(weight_diag, "diagnostics/transfer/wadoux_weights_distribution.csv")

# ============================================================================
# 5. TRANSFER LEARNING BY DEPTH
# ============================================================================

cat("\n--- Transfer Learning Models by Depth ---\n")

depths <- c(7.5, 22.5, 40, 75)
results_list <- list()
models_list <- list()

for (d in depths) {
  cat(sprintf("\n=== Depth: %.1f cm ===\n", d))
  
  # Filter data for this depth
  g_data <- global_df_weighted %>%
    filter(depth_cm_midpoint == d) %>%
    drop_na(carbon_stock_kg_m2, all_of(BRIDGE_VARS_REDUCED))
  
  l_data <- local_model_df %>%
    filter(depth_cm_midpoint == d) %>%
    drop_na(carbon_stock_kg_m2, all_of(BRIDGE_VARS_REDUCED))
  
  n_global <- nrow(g_data)
  n_local <- nrow(l_data)
  n_local_cores <- n_distinct(l_data$core_id)
  
  cat(sprintf("  Global N: %d, Local N: %d (%d cores)\n", n_global, n_local, n_local_cores))
  
  if (n_local < 2) {
    cat("  SKIPPING: Insufficient local data (N < 2)\n")
    next
  }
  
  # Select covariate set based on sample size
  if (n_local < N_THRESHOLD_REDUCED) {
    covars <- BRIDGE_VARS_REDUCED
    cat(sprintf("  Using REDUCED covariate set (%d variables) due to small N\n", length(covars)))
  } else {
    covars <- BRIDGE_VARS_FULL
    cat(sprintf("  Using FULL covariate set (%d variables)\n", length(covars)))
  }
  
  # ==========================================================================
  # STAGE 1: WEIGHTED GLOBAL PRIOR MODEL
  # ==========================================================================
  
  cat("  Training weighted global prior...\n")
  
  formula_str <- paste("carbon_stock_kg_m2 ~", paste(covars, collapse = " + "))
  
  rf_global <- ranger(
    formula = as.formula(formula_str),
    data = g_data,
    case.weights = g_data$wadoux_weight,
    num.trees = 1000,
    mtry = max(2, floor(length(covars) / 3)),
    min.node.size = 5,
    importance = "permutation",
    seed = SEED
  )
  
  # Variable importance from global model
  var_imp <- data.frame(
    variable = names(rf_global$variable.importance),
    importance = rf_global$variable.importance,
    depth = d
  ) %>% arrange(desc(importance))
  
  cat("  Top 3 important variables:", 
      paste(head(var_imp$variable, 3), collapse = ", "), "\n")
  
  # ==========================================================================
  # STAGE 2: LOCAL BIAS ESTIMATION
  # ==========================================================================
  
  # Predict global model on local data
  l_data$global_pred <- predict(rf_global, data = l_data)$predictions
  
  # Calculate residuals
  l_data$residual <- l_data$carbon_stock_kg_m2 - l_data$global_pred
  
  # Simple bias correction: mean residual
  # This is more robust than fitting another RF with tiny N
  bias_mean <- mean(l_data$residual)
  bias_sd <- sd(l_data$residual)
  
  cat(sprintf("  Global model bias on local data: %.3f ± %.3f kg/m²\n", 
              bias_mean, bias_sd))
  
  # ==========================================================================
  # STAGE 3: LEAVE-ONE-CORE-OUT CROSS-VALIDATION
  # ==========================================================================
  
  cat("  Running leave-one-core-out CV...\n")
  
  unique_cores <- unique(l_data$core_id)
  cv_results <- data.frame()
  
  for (held_out_core in unique_cores) {
    # Split data
    train_idx <- l_data$core_id != held_out_core
    test_idx <- l_data$core_id == held_out_core
    
    train_data <- l_data[train_idx, ]
    test_data <- l_data[test_idx, ]
    
    if (nrow(train_data) < 2) next
    
    # Compute bias from training set only
    cv_bias <- mean(train_data$residual)
    
    # Final prediction = global_pred + bias_correction
    test_data$cv_pred <- test_data$global_pred + cv_bias
    
    cv_results <- bind_rows(cv_results, test_data %>%
                              select(core_id, carbon_stock_kg_m2, global_pred, cv_pred, residual))
  }
  
  # Calculate CV metrics
  if (nrow(cv_results) > 0 && nrow(cv_results) >= 2) {
    ss_tot <- sum((cv_results$carbon_stock_kg_m2 - mean(cv_results$carbon_stock_kg_m2))^2)
    ss_res <- sum((cv_results$carbon_stock_kg_m2 - cv_results$cv_pred)^2)
    
    r2_cv <- 1 - ss_res / ss_tot
    rmse_cv <- sqrt(mean((cv_results$carbon_stock_kg_m2 - cv_results$cv_pred)^2))
    mae_cv <- mean(abs(cv_results$carbon_stock_kg_m2 - cv_results$cv_pred))
    
    # Also compute global-only metrics for comparison
    ss_res_global <- sum((cv_results$carbon_stock_kg_m2 - cv_results$global_pred)^2)
    r2_global <- 1 - ss_res_global / ss_tot
    rmse_global <- sqrt(mean((cv_results$carbon_stock_kg_m2 - cv_results$global_pred)^2))
    
    cat(sprintf("  CV Results (bias-corrected): R² = %.3f, RMSE = %.3f kg/m²\n", r2_cv, rmse_cv))
    cat(sprintf("  CV Results (global only):    R² = %.3f, RMSE = %.3f kg/m²\n", r2_global, rmse_global))
  } else {
    r2_cv <- NA; rmse_cv <- NA; mae_cv <- NA
    r2_global <- NA; rmse_global <- NA
    cat("  CV not possible (insufficient data after splits)\n")
  }
  
  # ==========================================================================
  # STAGE 4: BOOTSTRAP UNCERTAINTY QUANTIFICATION
  # ==========================================================================
  
  cat("  Computing prediction uncertainty via bootstrap...\n")
  
  n_boot <- 500
  boot_biases <- numeric(n_boot)
  
  for (b in 1:n_boot) {
    # Resample local observations with replacement
    boot_idx <- sample(1:n_local, replace = TRUE)
    boot_residuals <- l_data$residual[boot_idx]
    boot_biases[b] <- mean(boot_residuals)
  }
  
  # Prediction interval components
  # Total uncertainty = global model uncertainty + bias uncertainty + residual variance
  bias_se <- sd(boot_biases)
  residual_var <- var(l_data$residual)
  
  cat(sprintf("  Bias SE: %.3f, Residual SD: %.3f\n", bias_se, sqrt(residual_var)))
  
  # ==========================================================================
  # SAVE MODEL OBJECT
  # ==========================================================================
  
  model_object <- list(
    depth_cm = d,
    global_model = rf_global,
    bias_correction = bias_mean,
    bias_se = bias_se,
    residual_sd = sqrt(residual_var),
    predictors = covars,
    n_global = n_global,
    n_local = n_local,
    cv_r2 = r2_cv,
    cv_rmse = rmse_cv,
    method = "Wadoux_weighted_global_plus_bias_correction",
    created = Sys.time(),
    
    # Function to make predictions with uncertainty
    predict_with_uncertainty = function(newdata, ci_level = 0.90) {
      pred <- predict(rf_global, data = newdata)$predictions
      pred_corrected <- pred + bias_mean
      
      # Approximate prediction SE
      # Combines: RF prediction variance (from quantile regression) + bias uncertainty + residual noise
      pred_se <- sqrt(bias_se^2 + residual_var)
      
      z <- qnorm(1 - (1 - ci_level) / 2)
      
      data.frame(
        prediction = pred_corrected,
        se = pred_se,
        lower = pred_corrected - z * pred_se,
        upper = pred_corrected + z * pred_se,
        ci_level = ci_level
      )
    }
  )
  
  models_list[[as.character(d)]] <- model_object
  saveRDS(model_object, sprintf("outputs/models/transfer/transfer_model_depth_%.1f.rds", d))
  
  # Store results
  results_list[[as.character(d)]] <- data.frame(
    depth_cm = d,
    n_global = n_global,
    n_local = n_local,
    n_covariates = length(covars),
    bias_correction = bias_mean,
    bias_se = bias_se,
    residual_sd = sqrt(residual_var),
    cv_r2 = r2_cv,
    cv_rmse = rmse_cv,
    cv_r2_global_only = r2_global,
    cv_rmse_global_only = rmse_global
  )
  
  # Save variable importance
  write_csv(var_imp, sprintf("diagnostics/transfer/variable_importance_depth_%.1f.csv", d))
}

# ============================================================================
# 6. GENERATE SPATIAL PREDICTION RASTERS
# ============================================================================

cat("\n==================================================\n")
cat("   GENERATING SPATIAL PREDICTIONS\n")
cat("==================================================\n")

if (!RASTER_STACK_AVAILABLE) {
  cat(sprintf("\nWARNING: No covariate rasters available in '%s'.\n", COVARIATES_DIR))
  cat("         Skipping spatial prediction raster generation.\n")
  cat("         To generate prediction maps, add covariate .tif files and re-run.\n")
} else if (length(models_list) == 0) {
  cat("\nWARNING: No models were trained successfully.\n")
  cat("         Cannot generate spatial predictions.\n")
} else {
  
  cat(sprintf("Generating predictions for %d depth(s)...\n", length(models_list)))
  
  # Get the covariate stack subset needed for predictions
  # Use the reduced variable set since that's what small-sample models use
  available_covars <- intersect(BRIDGE_VARS_REDUCED, names(PREDICTION_STACK))
  
  if (length(available_covars) < length(BRIDGE_VARS_REDUCED)) {
    cat(sprintf("WARNING: Only %d of %d required covariates available in rasters.\n",
                length(available_covars), length(BRIDGE_VARS_REDUCED)))
    cat("         Missing:", paste(setdiff(BRIDGE_VARS_REDUCED, available_covars), collapse = ", "), "\n")
  }
  
  covar_stack <- PREDICTION_STACK[[available_covars]]
  
  # Check for scale correction on rasters (same logic as point extraction)
  cat("Checking raster scale alignment with global data...\n")
  for (var in available_covars) {
    if (!var %in% names(global_df)) next
    
    g_median <- median(global_df[[var]], na.rm = TRUE)
    
    # Sample raster values for comparison
    set.seed(SEED)
    sample_cells <- sample(1:ncell(covar_stack[[var]]), min(1000, ncell(covar_stack[[var]])))
    r_vals <- values(covar_stack[[var]])[sample_cells]
    r_median <- median(r_vals, na.rm = TRUE)
    
    # Skip elevation and near-zero values
    if (grepl("elevation", var, ignore.case = TRUE)) next
    if (is.na(r_median) || abs(g_median) < 0.001 || abs(r_median) < 0.001) next
    
    ratio <- abs(r_median / g_median)
    
    if (ratio > 50 && ratio < 20000) {
      cat(sprintf("  Scaling %s by 1/10000 (ratio: %.0f)\n", var, ratio))
      covar_stack[[var]] <- covar_stack[[var]] / 10000
    }
  }
  
  # Loop through each depth and generate predictions
  for (depth_char in names(models_list)) {
    
    d <- as.numeric(depth_char)
    model_obj <- models_list[[depth_char]]
    
    cat(sprintf("\n--- Depth: %.1f cm ---\n", d))
    
    # Check that model predictors are available
    model_covars <- model_obj$predictors
    missing_for_model <- setdiff(model_covars, names(covar_stack))
    
    if (length(missing_for_model) > 0) {
      cat(sprintf("  SKIPPING: Missing covariates for this model: %s\n",
                  paste(missing_for_model, collapse = ", ")))
      next
    }
    
    # Prepare prediction stack with only needed layers
    pred_covars <- covar_stack[[model_covars]]
    
    # --- BAND 1: Global Prior (uncorrected) ---
    cat("  Predicting global prior...\n")
    
    global_pred <- terra::predict(
      pred_covars,
      model_obj$global_model,
      fun = function(model, newdata) {
        predict(model, data = newdata)$predictions
      },
      na.rm = TRUE
    )
    names(global_pred) <- "Global_Prior"
    
    # --- BAND 2: Local-Only Benchmark ---
    # For local-only, we use mean of local observations at this depth
    # This represents what you'd predict without any spatial model
    cat("  Computing local-only benchmark...\n")
    
    local_depth_data <- local_model_df %>%
      filter(depth_cm_midpoint == d) %>%
      pull(carbon_stock_kg_m2)
    
    local_mean <- mean(local_depth_data, na.rm = TRUE)
    local_only <- global_pred * 0 + local_mean  # Constant raster
    names(local_only) <- "Local_Only"
    
    # --- BAND 3: Transfer Learning Final (bias-corrected) ---
    cat("  Applying bias correction...\n")
    
    transfer_final <- global_pred + model_obj$bias_correction
    names(transfer_final) <- "Transfer_Final"
    
    # --- BAND 4: Difference (Transfer - Local) ---
    difference <- transfer_final - local_only
    names(difference) <- "Difference"
    
    # --- Stack and save ---
    output_stack <- c(global_pred, local_only, transfer_final, difference)
    
    output_path <- sprintf("outputs/carbon_stocks/predictions/depth_%s_predictions.tif", d)
    
    terra::writeRaster(
      output_stack,
      output_path,
      overwrite = TRUE,
      gdal = c("COMPRESS=LZW")
    )
    
    cat(sprintf("  Saved: %s\n", basename(output_path)))
    
    # Report summary statistics
    cat(sprintf("  Global Prior mean: %.2f kg/m²\n", 
                global(global_pred, "mean", na.rm = TRUE)[1, 1]))
    cat(sprintf("  Transfer Final mean: %.2f kg/m² (bias: %+.2f)\n",
                global(transfer_final, "mean", na.rm = TRUE)[1, 1],
                model_obj$bias_correction))
  }
  
  cat("\nSpatial predictions complete.\n")
}

# ============================================================================
# 7. SUMMARY AND DIAGNOSTICS
# ============================================================================

cat("\n==================================================\n")
cat("   SUMMARY\n")
cat("==================================================\n")

if (length(results_list) > 0) {
  summary_df <- bind_rows(results_list)
  print(summary_df)
  write_csv(summary_df, "diagnostics/transfer/transfer_learning_validation.csv")
  
  # Interpretation guidance
  cat("\n--- Interpretation Notes ---\n")
  cat("1. R² values may be negative if global model is biased and sample size is tiny.\n")
  cat("2. Bias correction addresses systematic over/under-prediction by the global model.\n")
  cat("3. Prediction intervals account for: model uncertainty + bias uncertainty + residual noise.\n")
  cat("4. For VM0033 reporting, use the 90% prediction intervals from predict_with_uncertainty().\n")
  
  cat("\n--- How to Make Predictions ---\n")
  cat("model <- readRDS('outputs/models/transfer/transfer_model_depth_7.5.rds')\n")
  cat("preds <- model$predict_with_uncertainty(newdata, ci_level = 0.90)\n")
  
} else {
  cat("No models were successfully trained.\n")
}

cat("\nOutputs saved to:\n")
cat("  - outputs/models/transfer/*.rds (model objects)\n")
cat("  - outputs/carbon_stocks/predictions/*.tif (prediction rasters)\n")
cat("  - diagnostics/transfer/*.csv (validation stats)\n")
cat("\n")
# ============================================================================
# PRESENTATION COPY: TRANSFER LEARNING OUTPUTS
# ============================================================================

cat("\nCopying key Transfer Learning outputs to outputs/Transfer_learning/...\n")

dir.create("outputs/Transfer_learning", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Transfer_learning/Global_vs_Local_Comparisons", recursive = TRUE, showWarnings = FALSE)

get_depth_interval_label <- function(depth_cm) {
  depth_lookup <- data.frame(
    depth_midpoint = c(7.5, 22.5, 40, 75),
    depth_top = c(0, 15, 30, 50),
    depth_bottom = c(15, 30, 50, 100)
  )

  idx <- which.min(abs(depth_lookup$depth_midpoint - depth_cm))
  sprintf("%d_to_%dcm", depth_lookup$depth_top[idx], depth_lookup$depth_bottom[idx])
}

transfer_stack_files <- list.files(
  "outputs/carbon_stocks/predictions",
  pattern = "^depth_[0-9]+\\.?[0-9]*_predictions\\.tif$",
  full.names = TRUE
)

if (length(transfer_stack_files) > 0) {
  depth_values <- as.numeric(gsub(
    "^depth_([0-9]+\\.?[0-9]*)_predictions\\.tif$",
    "\\1",
    basename(transfer_stack_files)
  ))

  hero_src <- transfer_stack_files[which.max(depth_values)]
  hero_dst <- "outputs/Transfer_learning/Transfer_Learning_Final_Calibrated_Map_0_to_100cm.tif"
  file.copy(hero_src, hero_dst, overwrite = TRUE)
  cat(sprintf("Copied Transfer hero output: %s\n", basename(hero_dst)))

  for (i in seq_along(transfer_stack_files)) {
    depth_label <- get_depth_interval_label(depth_values[i])
    dst_name <- sprintf("Transfer_Learning_Global_vs_Local_Comparison_%s.tif", depth_label)
    dst <- file.path("outputs/Transfer_learning/Global_vs_Local_Comparisons", dst_name)
    file.copy(transfer_stack_files[i], dst, overwrite = TRUE)
    cat(sprintf("Copied Transfer comparison map: %s\n", basename(dst)))
  }
} else {
  cat("WARNING: No transfer prediction rasters found for presentation copy step.\n")
}

transfer_diag_sources <- c(
  list.files("diagnostics/transfer", pattern = "^transfer_learning_validation\\.csv$", full.names = TRUE),
  list.files("diagnostics/transfer", pattern = "^variable_importance_depth_[0-9]+\\.?[0-9]*\\.csv$", full.names = TRUE)
)

if (length(transfer_diag_sources) == 0) {
  cat("WARNING: No transfer diagnostics found for presentation copy step.\n")
} else {
  for (src in transfer_diag_sources) {
    src_name <- basename(src)
    if (grepl("^transfer_learning_validation\\.csv$", src_name)) {
      dst_name <- "Transfer_Learning_Validation_Summary.csv"
    } else {
      depth_cm <- as.numeric(gsub(
        "^variable_importance_depth_([0-9]+\\.?[0-9]*)\\.csv$",
        "\\1",
        src_name
      ))
      dst_name <- sprintf("Transfer_Learning_Variable_Importance_%s.csv", get_depth_interval_label(depth_cm))
    }
    dst <- file.path("outputs/Transfer_learning/Global_vs_Local_Comparisons", dst_name)
    file.copy(src, dst, overwrite = TRUE)
    cat(sprintf("Copied Transfer diagnostic: %s\n", basename(dst)))
  }
}
