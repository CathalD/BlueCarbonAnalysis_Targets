# ============================================================================
# MODULE 04: BLUE CARBON STRATIFIED KRIGING PREDICTIONS
# ============================================================================
# PURPOSE: Spatial interpolation of carbon stocks (kg/m²) using stratified
#          kriging by ecosystem type. Models harmonized carbon stocks directly
#          rather than SOC concentration.
# INPUTS:
#   - data_processed/cores_harmonized_bluecarbon.rds (from Module 03)
#   - covariates/*.tif (optional, for covariate-assisted kriging)
# OUTPUTS:
#   - outputs/predictions/kriging/carbon_stock_maps/*.tif    (kg/m² predictions by stratum/depth)
#   - outputs/predictions/kriging/standard_error_maps/*.tif  (combined SE, kg/m²)
#   - outputs/predictions/kriging/variance_maps/*.tif        (kriging + harmonization variance)
#   - outputs/predictions/kriging/aggregated_stocks/*.tif    (Mg C/ha per depth layer + totals)
#   - outputs/models/kriging/*.rds                           (variogram models, CV results)
#   - outputs/Basic_analysis/Spatial_Maps_by_Depth/          (PNG summary maps)
# ============================================================================

# ============================================================================
# SETUP
# ============================================================================

# Load configuration
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
if (file.exists("Analysis_Workflow/blue_carbon_config.R")) {
  source("Analysis_Workflow/blue_carbon_config.R")
} else {
  stop("Configuration file not found. Run 00b_setup_directories.R first.")
}

# Initialize logging
log_file <- file.path("logs", paste0("kriging_", Sys.Date(), ".log"))

log_message <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- sprintf("[%s] %s: %s", timestamp, level, msg)
  cat(log_entry, "\n")
  cat(log_entry, "\n", file = log_file, append = TRUE)
}

log_message("=== MODULE 04: STRATIFIED KRIGING ===")

# Load packages
suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(terra)
  library(gstat)
})

# Check for automap (optional)
has_automap <- requireNamespace("automap", quietly = TRUE)
if (has_automap) {
  library(automap)
  log_message("automap package available - using automatic fitting")
} else {
  log_message("automap not available - using manual fitting", "WARNING")
}

# Check for CAST (optional, for spatial CV and AOA)
has_cast <- requireNamespace("CAST", quietly = TRUE)
if (has_cast) {
  library(CAST)
  log_message("CAST package available - using spatial CV and AOA")
} else {
  log_message("CAST not available - using standard k-fold CV", "WARNING")
}

log_message("Packages loaded successfully")

# ============================================================================
# KRIGING METHOD OPTIONS
# ============================================================================

# Compare kriging methods (Ordinary vs Simple)
# Set to TRUE to perform both and compare performance
COMPARE_KRIGING_METHODS <- TRUE

# Default kriging method if not comparing: "ordinary" or "simple"
DEFAULT_KRIGING_METHOD <- "ordinary"

log_message(sprintf("Kriging method comparison: %s",
                    ifelse(COMPARE_KRIGING_METHODS, "Enabled (OK vs SK)",
                           paste0("Disabled (", DEFAULT_KRIGING_METHOD, " only)"))))

# Create output directories (subdirectory structure created per-stratum during processing)
dir.create("outputs/predictions/kriging",                   recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/predictions/kriging/carbon_stock_maps", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/predictions/kriging/standard_error_maps", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/predictions/kriging/variance_maps",     recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/predictions/kriging/aggregated_stocks", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/models/kriging",                        recursive = TRUE, showWarnings = FALSE)
dir.create("diagnostics/variograms",                        recursive = TRUE, showWarnings = FALSE)
dir.create("diagnostics/crossvalidation",                   recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD DATA
# ============================================================================

log_message("Loading harmonized data...")

if (!file.exists("data_processed/cores_harmonized_bluecarbon.rds")) {
  stop("ERROR: Harmonized data missing. Run Module 03 first.")
}

cores_harmonized <- readRDS("data_processed/cores_harmonized_bluecarbon.rds") %>%
  filter(depth_cm_midpoint %in% STANDARD_DEPTHS) %>%
  filter(qa_realistic) %>%
  mutate(depth_cm = depth_cm_midpoint)  # Create depth_cm column for compatibility

# DATA QUALITY CHECKS
log_message(sprintf("Loaded: %d predictions from %d cores", 
                    nrow(cores_harmonized), n_distinct(cores_harmonized$core_id)))

# Check for required columns
required_cols <- c("core_id", "longitude", "latitude", "depth_cm_midpoint", "carbon_stock_kg_m2", "stratum")
missing_cols <- setdiff(required_cols, names(cores_harmonized))
if (length(missing_cols) > 0) {
  stop(sprintf("ERROR: Missing required columns: %s", paste(missing_cols, collapse=", ")))
}

# Load harmonization metadata
harmonization_metadata <- NULL
if (file.exists("data_processed/harmonization_metadata.rds")) {
  harmonization_metadata <- readRDS("data_processed/harmonization_metadata.rds")
  log_message(sprintf("Harmonization method: %s", harmonization_metadata$method))
}

# Load Module 01 QA data
vm0033_compliance <- NULL
if (file.exists("data_processed/vm0033_compliance.rds")) {
  vm0033_compliance <- readRDS("data_processed/vm0033_compliance.rds")
  log_message("Loaded sampling adequacy data")
}

# Filter to standard depths only (VM0033 midpoints)
cores_standard <- cores_harmonized %>%
  filter(depth_cm %in% STANDARD_DEPTHS)

log_message(sprintf(vm_label("VM0033 target depths: %s cm", "Standard target depths: %s cm"), paste(STANDARD_DEPTHS, collapse = ", ")))
log_message(sprintf("Standard depth predictions: %d", nrow(cores_standard)))
depth_cm = cores_harmonized$depth_cm_midpoint
# Check for required columns
required_cols <- c("core_id", "longitude", "latitude", "stratum",
                   "depth_cm_midpoint", "carbon_stock_kg_m2")
missing <- setdiff(required_cols, names(cores_standard))
if (length(missing) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")))
}

log_message(sprintf("Modeling carbon stocks in kg/m² (standard units throughout workflow)"))

# Standardize core type names if present
if ("core_type" %in% names(cores_standard) || "core_type_clean" %in% names(cores_standard)) {
  if (!"core_type_clean" %in% names(cores_standard)) {
    cores_standard <- cores_standard %>%
      mutate(
        core_type_clean = case_when(
          tolower(core_type) %in% c("hr", "high-res", "high resolution", "high res") ~ "HR",
          tolower(core_type) %in% c("paired composite", "paired comp", "paired") ~ "Paired Composite",
          tolower(core_type) %in% c("unpaired composite", "unpaired comp", "unpaired", "composite", "comp") ~ "Unpaired Composite",
          TRUE ~ ifelse(is.na(core_type), "Unknown", core_type)
        )
      )
  }
  
  log_message("Core type distribution:")
  core_type_summary <- cores_standard %>%
    distinct(core_id, core_type_clean) %>%
    count(core_type_clean)
  for (i in 1:nrow(core_type_summary)) {
    log_message(sprintf("  %s: %d cores",
                        core_type_summary$core_type_clean[i],
                        core_type_summary$n[i]))
  }
}

# ============================================================================
# SAMPLE SIZE ASSESSMENT & RECOMMENDATIONS
# ============================================================================

log_message("\n=== SAMPLE SIZE ASSESSMENT ===")

sample_assessment <- cores_standard %>%
  group_by(stratum) %>%
  summarise(
    n_cores = n_distinct(core_id),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  mutate(
    # Assessment categories
    status = case_when(
      n_cores >= 10 ~ "SUFFICIENT",
      n_cores >= 5 ~ "MINIMUM MET",
      n_cores >= 3 ~ "BELOW RECOMMENDED",
      TRUE ~ "INSUFFICIENT"
    ),
    # Calculate additional cores needed
    additional_for_minimum = pmax(0, 5 - n_cores),
    additional_for_recommended = pmax(0, 10 - n_cores),
    # Determine approach
    analysis_approach = case_when(
      n_cores >= 5 ~ "Stratified kriging",
      n_cores >= 3 ~ "Stratified kriging (high uncertainty)",
      TRUE ~ "Pool with other strata"
    )
  )

cat("\n========================================\n")
cat("SAMPLE SIZE STATUS BY STRATUM\n")
cat("========================================\n\n")

for (i in 1:nrow(sample_assessment)) {
  cat(sprintf("Stratum: %s\n", sample_assessment$stratum[i]))
  cat(sprintf("  Current: %d cores (%d samples)\n", 
              sample_assessment$n_cores[i], 
              sample_assessment$n_samples[i]))
  cat(sprintf("  Status: %s\n", sample_assessment$status[i]))
  cat(sprintf("  Approach: %s\n", sample_assessment$analysis_approach[i]))
  
  if (sample_assessment$additional_for_minimum[i] > 0) {
    cat(sprintf("  ⚠️  ACTION NEEDED: Collect %d more cores to meet the statistical best-practice minimum (n=5)\n",
                sample_assessment$additional_for_minimum[i]))
  }
  
  if (sample_assessment$additional_for_recommended[i] > 0) {
    cat(sprintf("  📊 RECOMMENDED: Collect %d more cores for robust estimates (n=10)\n",
                sample_assessment$additional_for_recommended[i]))
  }
  
  cat("\n")
}

# Overall site assessment
total_cores <- sum(sample_assessment$n_cores)
n_sufficient <- sum(sample_assessment$status == "SUFFICIENT")
n_strata <- nrow(sample_assessment)

cat("========================================\n")
cat("OVERALL SITE ASSESSMENT\n")
cat("========================================\n\n")
cat(sprintf("Total cores collected: %d\n", total_cores))
cat(sprintf("Strata with sufficient data: %d/%d\n", n_sufficient, n_strata))

if (n_sufficient == n_strata) {
  cat("✅ STATUS: Ready for high-confidence carbon stock estimation\n\n")
} else if (total_cores >= 3 * n_strata) {
  cat("⚠️  STATUS: Analysis possible but with increased uncertainty\n")
  cat("   RECOMMENDATION: Collect additional samples before crediting\n\n")
} else {
  cat("❌ STATUS: Insufficient samples for stratified analysis\n")
  cat("   RECOMMENDATION: Pool strata or collect more samples\n\n")
}

# Determine analysis strategy
if (all(sample_assessment$n_cores < 5)) {
  log_message("ANALYSIS STRATEGY: Pooling all strata due to small sample sizes", "WARNING")
  log_message("This provides site-wide estimates but loses stratum-specific detail", "INFO")
  
  POOL_STRATA <- TRUE
  
  cores_standard <- cores_standard %>%
    mutate(stratum_original = stratum,
           stratum = "Site_Wide_Pooled")
  
  strata <- c("Site_Wide_Pooled")
  
} else if (any(sample_assessment$n_cores < 5)) {
  log_message("ANALYSIS STRATEGY: Mixed - some strata will be pooled", "WARNING")
  
  POOL_STRATA <- FALSE
  
  # Pool only insufficient strata
  insufficient_strata <- sample_assessment$stratum[sample_assessment$n_cores < 3]
  
  if (length(insufficient_strata) > 0) {
    log_message(sprintf("Pooling insufficient strata: %s", 
                        paste(insufficient_strata, collapse = ", ")))
    
    cores_standard <- cores_standard %>%
      mutate(stratum_original = stratum,
             stratum = ifelse(stratum %in% insufficient_strata, 
                              "Pooled_Insufficient", 
                              stratum))
  }
  
  strata <- unique(cores_standard$stratum)
  
} else {
  log_message("ANALYSIS STRATEGY: Full stratified analysis", "INFO")
  POOL_STRATA <- FALSE
  strata <- unique(cores_standard$stratum)
}

# Save assessment for reporting
saveRDS(sample_assessment, "diagnostics/sample_size_assessment.rds")
write.csv(sample_assessment, "diagnostics/sample_size_assessment.csv", row.names = FALSE)

log_message("Sample size assessment complete")


# ============================================================================
# LOAD GEE STRATUM BOUNDARY MASKS
# ============================================================================

log_message("Loading GEE stratum boundary masks...")

# Directory for GEE exported strata
gee_strata_dir <- file.path(DATA_RAW_DIR, "gee_strata")   # DATA_RAW_DIR set in config

# Load stratum masks if available
stratum_masks <- list()

if (dir.exists(gee_strata_dir)) {
  # Expected GEE export file patterns
  # Normalize stratum names to match file naming conventions
  stratum_file_map <- c(
    "Upper Marsh" = "upper_marsh.tif",
    "Mid Marsh" = "mid_marsh.tif",
    "Lower Marsh" = "lower_marsh.tif",
    "Underwater Vegetation" = "underwater_vegetation.tif",
    "Open Water" = "open_water.tif"
  )
  
  for (stratum_name in names(stratum_file_map)) {
    file_path <- file.path(gee_strata_dir, stratum_file_map[stratum_name])
    
    if (file.exists(file_path)) {
      mask_rast <- tryCatch({
        r <- rast(file_path)
        
        # Reproject to processing CRS if needed
        if (st_crs(r)$input != paste0("EPSG:", PROCESSING_CRS)) {
          r <- project(r, paste0("EPSG:", PROCESSING_CRS), method = "near")
        }
        
        # Binary mask: 1 = inside stratum, 0/NA = outside
        r[r == 0] <- NA
        
        r
      }, error = function(e) {
        log_message(sprintf("  Failed to load %s: %s", stratum_name, e$message), "WARNING")
        NULL
      })
      
      if (!is.null(mask_rast)) {
        stratum_masks[[stratum_name]] <- mask_rast
        log_message(sprintf("  Loaded: %s (%d x %d cells)",
                            stratum_name,
                            ncol(mask_rast),
                            nrow(mask_rast)))
      }
    }
  }
  
  if (length(stratum_masks) > 0) {
    log_message(sprintf("Loaded %d stratum boundary masks from GEE", length(stratum_masks)))
  } else {
    log_message("No GEE stratum masks found - predictions will use buffered extent", "WARNING")
  }
} else {
  log_message(sprintf("GEE strata directory not found: %s", gee_strata_dir), "WARNING")
  log_message("Predictions will use buffered extent without stratum masks", "WARNING")
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Fit variogram with automatic or manual model selection
#' @param sp_data Spatial data
#' @param formula Formula for kriging
#' @param cutoff Maximum distance for variogram
#' @return Best variogram model
fit_variogram_auto <- function(sp_data, formula, cutoff = NULL) {
  
  tryCatch({
    # Calculate empirical variogram
    vgm_emp <- variogram(formula, sp_data, cutoff = cutoff, width = KRIGING_WIDTH)
    
    # Use automap if available, otherwise manual fitting
    if (has_automap) {
      # Automatic model fitting with automap
      vgm_fit <- automap::autofitVariogram(
        formula = formula,
        input_data = sp_data,
        model = c("Sph", "Exp", "Gau"),
        cutoff = cutoff,
        width = KRIGING_WIDTH,
        verbose = FALSE
      )
      
      return(list(
        empirical = vgm_emp,
        model = vgm_fit$var_model,
        sserr = vgm_fit$sserr
      ))
      
    } else {
      # Manual variogram fitting (fallback)
      log_message("  Using manual variogram fitting", "INFO")
      
      # Get initial parameter estimates
      max_semivar <- max(vgm_emp$gamma, na.rm = TRUE)
      max_dist <- max(vgm_emp$dist, na.rm = TRUE)
      
      # Try multiple models
      models_to_try <- c("Sph", "Exp", "Gau")
      best_model <- NULL
      best_sse <- Inf
      
      for (model_name in models_to_try) {
        # Initial parameters
        init_model <- vgm(
          psill = max_semivar * 0.8,
          model = model_name,
          range = max_dist / 3,
          nugget = max_semivar * 0.1
        )
        
        # Fit model
        fitted_model <- tryCatch({
          fit.variogram(vgm_emp, init_model, fit.method = 7)
        }, error = function(e) NULL)
        
        if (!is.null(fitted_model)) {
          # Calculate SSE
          sse <- attr(fitted_model, "SSErr")
          if (!is.null(sse) && sse < best_sse) {
            best_sse <- sse
            best_model <- fitted_model
          }
        }
      }
      
      if (is.null(best_model)) {
        # Last resort: simple spherical model
        best_model <- vgm(
          psill = max_semivar * 0.8,
          model = "Sph",
          range = max_dist / 3,
          nugget = max_semivar * 0.1
        )
        best_sse <- NA
      }
      
      return(list(
        empirical = vgm_emp,
        model = best_model,
        sserr = best_sse
      ))
    }
    
  }, error = function(e) {
    log_message(sprintf("Variogram fitting error: %s", e$message), "ERROR")
    return(NULL)
  })
}

#' Test for spatial autocorrelation in residuals
#'
#' Uses Moran's I test to detect spatial autocorrelation in kriging residuals.
#' Significant autocorrelation suggests the variogram model may need adjustment.
#'
#' @param residuals Vector of kriging residuals
#' @param coords Coordinate matrix (x, y)
#' @param k Number of nearest neighbors for spatial weights (default 5)
#' @return Moran's I test result or NULL if test fails
#'
#' @details
#' Moran's I ranges from -1 (perfect dispersion) to +1 (perfect clustering).
#' Values near 0 indicate no spatial autocorrelation (desired for residuals).
#'
#' For kriging validation:
#' - Non-significant Moran's I → Model adequately captures spatial structure
#' - Significant positive I → Residuals still spatially correlated (underfitting)
#' - Significant negative I → Over-smoothing (less common)
#'
#' Requires spdep package. Returns NULL gracefully if unavailable.
#'
#' Reference: Moran, P.A.P. (1950). Notes on continuous stochastic phenomena.
#' Biometrika, 37(1/2), 17-23.
#'
#' @examples
#' moran_test <- test_spatial_autocorrelation(residuals, coords)
#' if (!is.null(moran_test) && moran_test$p.value < 0.05) {
#'   # Residuals show significant spatial autocorrelation
#' }
test_spatial_autocorrelation <- function(residuals, coords, k = 5) {
  
  # Check if spdep package is available
  if (!requireNamespace("spdep", quietly = TRUE)) {
    log_message("spdep package not available - skipping Moran's I test",
                level = "WARNING")
    log_message("Install with: install.packages('spdep')", level = "INFO")
    return(NULL)
  }
  
  # Need at least 10 points for meaningful test
  if (length(residuals) < 10) {
    log_message("Too few points for Moran's I test (need >= 10)",
                level = "INFO")
    return(NULL)
  }
  
  tryCatch({
    # Create spatial neighbors (k nearest neighbors)
    k_actual <- min(k, floor(length(residuals) / 2))  # Ensure k is reasonable
    nb <- spdep::knn2nb(spdep::knearneigh(coords, k = k_actual))
    
    # Convert to spatial weights list
    listw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
    
    # Perform Moran's I test
    moran_result <- spdep::moran.test(residuals, listw, zero.policy = TRUE)
    
    # Interpret results
    moran_i <- moran_result$estimate["Moran I statistic"]
    p_value <- moran_result$p.value
    
    if (p_value < 0.05) {
      if (moran_i > 0) {
        log_message(
          sprintf("Significant positive spatial autocorrelation detected (I = %.3f, p = %.4f)",
                  moran_i, p_value),
          level = "WARNING"
        )
        log_message("  → Residuals still spatially correlated - variogram may need adjustment",
                    level = "WARNING")
      } else {
        log_message(
          sprintf("Significant negative spatial autocorrelation detected (I = %.3f, p = %.4f)",
                  moran_i, p_value),
          level = "WARNING"
        )
        log_message("  → Possible over-smoothing", level = "WARNING")
      }
    } else {
      log_message(
        sprintf("No significant spatial autocorrelation in residuals (I = %.3f, p = %.4f)",
                moran_i, p_value),
        level = "INFO"
      )
      log_message("  ✓ Variogram model adequately captures spatial structure",
                  level = "INFO")
    }
    
    return(moran_result)
    
  }, error = function(e) {
    log_message(sprintf("Moran's I test failed: %s", e$message),
                level = "WARNING")
    return(NULL)
  })
}

#' Perform cross-validation for kriging
#' @param sp_data Spatial data
#' @param formula Formula for kriging
#' @param vgm_model Variogram model
#' @param use_spatial_cv Use CAST spatial CV (if available)
#' @return CV results
crossvalidate_kriging <- function(sp_data, formula, vgm_model, use_spatial_cv = TRUE) {
  
  tryCatch({
    
    # Use CAST spatial CV if available
    if (use_spatial_cv && has_cast && length(sp_data) >= 10) {
      log_message("    Using CAST spatial cross-validation...")
      
      # Convert sp to sf for CAST
      sp_sf <- st_as_sf(sp_data)
      
      # Create spatial folds using CAST
      # CreateSpacetimeFolds separates training/test sets spatially
      spatial_folds <- tryCatch({
        CreateSpacetimeFolds(
          x = sp_sf,
          spacevar = "geometry",
          k = min(CV_FOLDS, floor(length(sp_data) / 3))  # Adjust folds for small datasets
        )
      }, error = function(e) {
        log_message(sprintf("    Spatial fold creation failed: %s", e$message), "WARNING")
        NULL
      })
      
      if (!is.null(spatial_folds)) {
        # Manual spatial CV with CAST folds
        predictions <- rep(NA, length(sp_data))
        observations <- sp_data[[all.vars(formula)[1]]]
        
        for (fold in 1:length(spatial_folds$index)) {
          train_idx <- spatial_folds$index[[fold]]
          test_idx <- spatial_folds$indexOut[[fold]]
          
          if (length(train_idx) >= 5 && length(test_idx) > 0) {
            # Train on fold
            train_data <- sp_data[train_idx, ]
            
            # Predict test set
            test_data <- sp_data[test_idx, ]
            
            pred_result <- tryCatch({
              krige(
                formula = formula,
                locations = train_data,
                newdata = test_data,
                model = vgm_model
              )
            }, error = function(e) {
              NULL
            })
            
            if (!is.null(pred_result)) {
              predictions[test_idx] <- pred_result$var1.pred
            }
          }
        }
        
        # Calculate metrics
        valid_idx <- !is.na(predictions) & !is.na(observations)
        if (sum(valid_idx) > 0) {
          residuals <- observations[valid_idx] - predictions[valid_idx]
          
          rmse <- sqrt(mean(residuals^2))
          mae <- mean(abs(residuals))
          me <- mean(residuals)
          
          ss_res <- sum(residuals^2)
          ss_tot <- sum((observations[valid_idx] - mean(observations[valid_idx]))^2)
          r2 <- 1 - (ss_res / ss_tot)
          
          return(list(
            rmse = rmse,
            mae = mae,
            me = me,
            r2 = r2,
            method = "spatial_cv",
            n_folds = length(spatial_folds$index),
            n_validated = sum(valid_idx)
          ))
        }
      }
    }
    
    # Fallback: standard k-fold CV
    log_message("    Using standard k-fold cross-validation...")
    
    cv_result <- krige.cv(
      formula = formula,
      locations = sp_data,
      model = vgm_model,
      nfold = CV_FOLDS,
      verbose = FALSE
    )
    
    # Calculate metrics
    rmse <- sqrt(mean(cv_result$residual^2, na.rm = TRUE))
    mae <- mean(abs(cv_result$residual), na.rm = TRUE)
    me <- mean(cv_result$residual, na.rm = TRUE)
    
    # R² calculation
    ss_res <- sum(cv_result$residual^2, na.rm = TRUE)
    ss_tot <- sum((cv_result$observed - mean(cv_result$observed, na.rm = TRUE))^2, na.rm = TRUE)
    r2 <- 1 - (ss_res / ss_tot)
    
    # Test for spatial autocorrelation in residuals
    log_message("    Testing spatial autocorrelation in residuals...", "INFO")
    coords <- st_coordinates(cv_result)
    moran_test <- test_spatial_autocorrelation(cv_result$residual, coords, k = 5)
    
    return(list(
      rmse = rmse,
      mae = mae,
      me = me,
      r2 = r2,
      method = "kfold_cv",
      cv_data = cv_result,
      moran_test = moran_test
    ))
    
  }, error = function(e) {
    log_message(sprintf("Cross-validation error: %s", e$message), "ERROR")
    return(NULL)
  })
}

#' Calculate Area of Applicability based on distance to samples
#' @param pred_raster Prediction raster
#' @param sample_points Sample locations (sf or sp)
#' @param max_distance Maximum acceptable distance (meters)
#' @return AOA raster (1 = inside AOA, NA = outside AOA)
calculate_spatial_aoa <- function(pred_raster, sample_points, max_distance = NULL) {
  
  tryCatch({
    # Convert sample points to sf if needed
    if (inherits(sample_points, "Spatial")) {
      sample_points <- st_as_sf(sample_points)
    }
    
    # Get prediction grid centroids
    pred_points <- as.points(pred_raster, na.rm = FALSE)
    pred_sf <- st_as_sf(pred_points)
    
    # Calculate distance to nearest sample for each prediction point
    dist_matrix <- st_distance(pred_sf, sample_points)
    min_dist <- apply(dist_matrix, 1, min)
    
    # Determine threshold
    if (is.null(max_distance)) {
      # Use median nearest-neighbor distance * 3 as threshold
      nn_dists <- apply(st_distance(sample_points, sample_points) + diag(Inf, nrow(sample_points)),
                        1, min)
      max_distance <- median(nn_dists, na.rm = TRUE) * 3
    }
    
    # Create AOA mask: 1 = inside AOA, NA = outside
    aoa_values <- ifelse(min_dist <= max_distance, 1, NA)
    
    # Convert to raster
    aoa_raster <- pred_raster
    values(aoa_raster) <- aoa_values
    
    # Calculate % area in AOA
    pct_in_aoa <- sum(!is.na(aoa_values)) / length(aoa_values) * 100
    
    return(list(
      aoa_raster = aoa_raster,
      max_distance = max_distance,
      pct_in_aoa = pct_in_aoa,
      min_distances = min_dist
    ))
    
  }, error = function(e) {
    log_message(sprintf("AOA calculation error: %s", e$message), "ERROR")
    return(NULL)
  })
}

#' Plot variogram with model fit
#' @param vgm_emp Empirical variogram
#' @param vgm_model Fitted model
#' @param stratum_name Stratum name
#' @param depth_cm Depth
plot_variogram <- function(vgm_emp, vgm_model, stratum_name, depth_cm) {
  
  tryCatch({
    png(file.path("diagnostics/variograms", 
                  sprintf("variogram_%s_%.0fcm.png", 
                          gsub(" ", "_", stratum_name), depth_cm)),
        width = 8, height = 6, units = "in", res = 300)
    
    plot(vgm_emp, vgm_model,
         main = sprintf("Carbon Stock Variogram: %s at %.0f cm", stratum_name, depth_cm),
         xlab = "Distance (m)",
         ylab = "Semivariance (kg/m²)²")
    
    dev.off()
    
  }, error = function(e) {
    log_message(sprintf("Variogram plot error: %s", e$message), "WARNING")
  })
}

# ============================================================================
# STRATIFIED KRIGING BY DEPTH
# ============================================================================

log_message("Starting stratified kriging by depth and stratum...")

# Initialize results storage
kriging_results <- list()
cv_results_all <- data.frame()
variogram_models <- list()

# Get unique strata
strata <- unique(cores_standard$stratum)
log_message(sprintf("Processing %d strata: %s", 
                    length(strata), paste(strata, collapse = ", ")))

# Process each depth
for (depth in STANDARD_DEPTHS) {
  
  log_message(sprintf("\n=== Processing depth: %.1f cm ===", depth))
  
  # Filter to this depth
  cores_depth <- cores_standard %>%
    filter(depth_cm == depth)
  
  log_message(sprintf("Samples at %.1f cm: %d from %d cores",
                      depth, nrow(cores_depth), n_distinct(cores_depth$core_id)))
  
  # Process each stratum
  for (stratum_name in strata) {
    
    log_message(sprintf("Processing: %s at %.1f cm", stratum_name, depth))
    
    # Filter to this stratum
    cores_stratum <- cores_depth %>%
      filter(stratum == stratum_name)
    
    n_samples <- nrow(cores_stratum)
    
    if (n_samples < 5) {
      log_message(sprintf("Skipping %s (n=%d, need ≥5)", stratum_name, n_samples), "WARNING")
      next
    }
    
    log_message(sprintf("  Samples: %d", n_samples))
    
    # Convert to spatial object
    cores_sf <- st_as_sf(cores_stratum,
                         coords = c("longitude", "latitude"),
                         crs = INPUT_CRS)
    
    # Transform to processing CRS
    cores_sf <- st_transform(cores_sf, PROCESSING_CRS)
    
    # Convert to sp object (required by gstat)
    cores_sp <- as(cores_sf, "Spatial")
    
    # ========================================================================
    # EXTRACT MODULE 03 UNCERTAINTIES
    # ========================================================================
    
    # Check if uncertainty fields are available from Module 03
    # Note: Carbon stock uncertainties would need to be propagated from SOC and BD uncertainties
    # For now, we use kriging variance as the primary uncertainty measure
    has_uncertainties <- "is_interpolated" %in% names(cores_stratum)
    
    harmonization_var_mean <- 0  # Default if no uncertainty data
    
    if (has_uncertainties) {
      # Check interpolation status
      n_interpolated <- sum(cores_stratum$is_interpolated, na.rm = TRUE)
      n_measured <- sum(!cores_stratum$is_interpolated, na.rm = TRUE)
      
      log_message(sprintf("  Depths: %d measured, %d interpolated",
                          n_measured, n_interpolated))
      
      # Optional: Include measurement CV if available
      if ("measurement_cv" %in% names(cores_stratum)) {
        mean_cv <- mean(cores_stratum$measurement_cv, na.rm = TRUE)
        if (!is.na(mean_cv) && mean_cv > 0) {
          log_message(sprintf("  Mean measurement CV: %.1f%%", mean_cv * 100))
        }
      }
      
      log_message("  Using kriging variance for uncertainty quantification")
    } else {
      log_message("  No Module 03 metadata - using kriging variance only", "WARNING")
    }
    
    # ========================================================================
    # FIT VARIOGRAM
    # ========================================================================
    
    log_message("  Fitting variogram...")
    
    # Determine cutoff (max 1/3 of study area extent)
    bbox <- st_bbox(cores_sf)
    max_dist <- sqrt((bbox$xmax - bbox$xmin)^2 + (bbox$ymax - bbox$ymin)^2)
    cutoff <- min(max_dist / 3, KRIGING_MAX_DISTANCE)
    
    # Fit variogram
    vgm_result <- fit_variogram_auto(
      sp_data = cores_sp,
      formula = carbon_stock_kg_m2 ~ 1,
      cutoff = cutoff
    )
    
    if (is.null(vgm_result)) {
      log_message("  Variogram fitting failed - skipping", "ERROR")
      next
    }
    
    log_message(sprintf("  Model: %s, SSErr: %.2f", 
                        vgm_result$model$model[2], vgm_result$sserr))
    
    # Plot variogram
    plot_variogram(vgm_result$empirical, vgm_result$model, 
                   stratum_name, depth)
    
    # Store variogram model
    model_key <- sprintf("%s_%.0fcm", gsub(" ", "_", stratum_name), depth)
    variogram_models[[model_key]] <- vgm_result$model
    
    # ========================================================================
    # CROSS-VALIDATION
    # ========================================================================
    
    log_message("  Performing cross-validation...")
    
    cv_result <- crossvalidate_kriging(
      sp_data = cores_sp,
      formula = carbon_stock_kg_m2 ~ 1,
      vgm_model = vgm_result$model
    )
    
    if (!is.null(cv_result)) {
      cv_method_label <- ifelse(!is.null(cv_result$method),
                                ifelse(cv_result$method == "spatial_cv", "Spatial CV", "K-fold CV"),
                                "K-fold CV")
      log_message(sprintf("  CV (%s) RMSE: %.2f, R²: %.3f",
                          cv_method_label, cv_result$rmse, cv_result$r2))
      
      # Store CV results for ordinary kriging
      cv_results_all <- rbind(cv_results_all, data.frame(
        stratum = stratum_name,
        depth_cm = depth,
        n_samples = n_samples,
        cv_rmse = cv_result$rmse,
        cv_mae = cv_result$mae,
        cv_me = cv_result$me,
        cv_r2 = cv_result$r2,
        cv_method = ifelse(!is.null(cv_result$method), cv_result$method, "kfold_cv"),
        kriging_type = "ordinary",
        model_type = vgm_result$model$model[2],
        sserr = vgm_result$sserr
      ))
      
      # ======================================================================
      # OPTIONAL: COMPARE WITH SIMPLE KRIGING
      # ======================================================================
      
      if (COMPARE_KRIGING_METHODS) {
        log_message("  Comparing with Simple Kriging...")
        
        # For simple kriging, we need to specify a known mean (beta parameter)
        # Use the sample mean for this stratum/depth
        global_mean <- mean(cores_sp$carbon_stock_kg_m2, na.rm = TRUE)
        
        # Simple kriging CV using beta parameter
        cv_result_sk <- tryCatch({
          # Note: gstat's krige.cv doesn't directly support beta parameter
          # We implement it by de-meaning the data, doing OK, then adding mean back
          cores_sp_demeaned <- cores_sp
          cores_sp_demeaned$carbon_stock_demeaned <- cores_sp$carbon_stock_kg_m2 - global_mean
          
          # Perform CV on demeaned data
          cv_sk <- krige.cv(
            formula = carbon_stock_demeaned ~ 1,
            locations = cores_sp_demeaned,
            model = vgm_result$model,
            nfold = CV_FOLDS,
            verbose = FALSE
          )
          
          # Add mean back to predictions
          cv_sk$var1.pred <- cv_sk$var1.pred + global_mean
          cv_sk$observed <- cv_sk$observed + global_mean
          cv_sk$residual <- cv_sk$observed - cv_sk$var1.pred
          
          # Calculate metrics
          rmse_sk <- sqrt(mean(cv_sk$residual^2, na.rm = TRUE))
          mae_sk <- mean(abs(cv_sk$residual), na.rm = TRUE)
          me_sk <- mean(cv_sk$residual, na.rm = TRUE)
          ss_res_sk <- sum(cv_sk$residual^2, na.rm = TRUE)
          ss_tot_sk <- sum((cv_sk$observed - mean(cv_sk$observed, na.rm = TRUE))^2, na.rm = TRUE)
          r2_sk <- 1 - (ss_res_sk / ss_tot_sk)
          
          list(rmse = rmse_sk, mae = mae_sk, me = me_sk, r2 = r2_sk)
        }, error = function(e) {
          log_message(sprintf("    Simple kriging CV failed: %s", e$message), "WARNING")
          NULL
        })
        
        if (!is.null(cv_result_sk)) {
          log_message(sprintf("  Simple Kriging CV RMSE: %.2f, R²: %.3f",
                              cv_result_sk$rmse, cv_result_sk$r2))
          
          # Compare methods
          rmse_diff <- cv_result$rmse - cv_result_sk$rmse
          better_method <- ifelse(cv_result_sk$rmse < cv_result$rmse, "Simple", "Ordinary")
          
          log_message(sprintf("  Method comparison: %s kriging performs better (ΔRMSE: %.2f)",
                              better_method, abs(rmse_diff)))
          
          # Store SK results
          cv_results_all <- rbind(cv_results_all, data.frame(
            stratum = stratum_name,
            depth_cm = depth,
            n_samples = n_samples,
            cv_rmse = cv_result_sk$rmse,
            cv_mae = cv_result_sk$mae,
            cv_me = cv_result_sk$me,
            cv_r2 = cv_result_sk$r2,
            cv_method = ifelse(!is.null(cv_result$method), cv_result$method, "kfold_cv"),
            kriging_type = "simple",
            model_type = vgm_result$model$model[2],
            sserr = vgm_result$sserr
          ))
        }
      }
    }
    
    # ========================================================================
    # CREATE PREDICTION GRID
    # ========================================================================
    
    log_message("  Creating prediction grid...")
    
    # Determine prediction extent
    # Use GEE mask extent if available, otherwise use buffered sample extent
    if (stratum_name %in% names(stratum_masks)) {
      # Use mask extent for more accurate prediction area
      extent_vec <- st_bbox(stratum_masks[[stratum_name]])
      log_message("  Using GEE mask extent for prediction grid")
    } else {
      # Fallback: buffer around sample locations
      buffer_dist <- 500  # 500m buffer
      bbox_buffered <- st_buffer(st_as_sfc(bbox), buffer_dist)
      extent_vec <- st_bbox(bbox_buffered)
      log_message("  Using buffered sample extent for prediction grid")
    }
    
    # Calculate grid dimensions
    x_range <- extent_vec["xmax"] - extent_vec["xmin"]
    y_range <- extent_vec["ymax"] - extent_vec["ymin"]
    ncols <- ceiling(x_range / KRIGING_CELL_SIZE)
    nrows <- ceiling(y_range / KRIGING_CELL_SIZE)
    
    # Create prediction grid
    pred_grid <- st_as_sf(
      expand.grid(
        x = seq(extent_vec["xmin"], extent_vec["xmax"], length.out = ncols),
        y = seq(extent_vec["ymin"], extent_vec["ymax"], length.out = nrows)
      ),
      coords = c("x", "y"),
      crs = PROCESSING_CRS
    )
    
    # Convert to sp
    pred_grid_sp <- as(pred_grid, "Spatial")
    
    log_message(sprintf("  Prediction grid: %d cells (%d x %d)",
                        length(pred_grid_sp), ncols, nrows))
    
    # ========================================================================
    # PERFORM KRIGING
    # ========================================================================
    
    log_message("  Performing kriging...")
    
    krige_result <- tryCatch({
      krige(
        formula = carbon_stock_kg_m2 ~ 1,
        locations = cores_sp,
        newdata = pred_grid_sp,
        model = vgm_result$model
      )
    }, error = function(e) {
      log_message(sprintf("  Kriging failed: %s", e$message), "ERROR")
      NULL
    })
    
    if (is.null(krige_result)) {
      next
    }
    
    # ========================================================================
    # SAVE RESULTS
    # ========================================================================
    
    log_message("  Saving predictions...")
    
    # Convert kriging results to raster
    # First convert to SpatVector, then rasterize
    krige_vect <- vect(krige_result)
    
    # Create empty raster template
    template_rast <- rast(
      extent = ext(krige_vect),
      resolution = KRIGING_CELL_SIZE,
      crs = sprintf("EPSG:%d", PROCESSING_CRS)
    )
    
    # Rasterize predictions
    pred_raster <- rasterize(
      x = krige_vect,
      y = template_rast,
      field = "var1.pred",
      fun = "mean"
    )
    
    # Rasterize variance
    var_raster <- rasterize(
      x = krige_vect,
      y = template_rast,
      field = "var1.var",
      fun = "mean"
    )
    
    # Apply stratum mask if available
    if (stratum_name %in% names(stratum_masks)) {
      log_message("  Applying GEE stratum boundary mask...")
      
      # Resample mask to match prediction raster resolution
      mask_resampled <- resample(stratum_masks[[stratum_name]],
                                 pred_raster,
                                 method = "near")
      
      # Apply mask (keep only cells within stratum boundary)
      pred_raster <- mask(pred_raster, mask_resampled)
      var_raster <- mask(var_raster, mask_resampled)
      
      log_message(sprintf("  Masked to stratum boundary (%d valid cells)",
                          sum(!is.na(values(pred_raster)))))
    }
    
    # ========================================================================
    # COMBINE UNCERTAINTIES
    # ========================================================================
    
    # Combine kriging variance with harmonization variance
    # Total variance = kriging variance + harmonization variance
    var_combined <- var_raster + harmonization_var_mean
    
    if (has_uncertainties) {
      log_message(sprintf("  Combined uncertainty: kriging var = %.2f, harmonization var = %.2f",
                          mean(values(var_raster), na.rm = TRUE),
                          harmonization_var_mean))
    }
    
    # Calculate standard errors
    se_kriging <- sqrt(var_raster)
    se_combined <- sqrt(var_combined)
    
    # Create organised output subdirectories
    dir.create("outputs/predictions/kriging/carbon_stock_maps",  recursive = TRUE, showWarnings = FALSE)
    dir.create("outputs/predictions/kriging/standard_error_maps", recursive = TRUE, showWarnings = FALSE)
    dir.create("outputs/predictions/kriging/variance_maps",       recursive = TRUE, showWarnings = FALSE)
    dir.create("outputs/predictions/kriging/aggregated_stocks",   recursive = TRUE, showWarnings = FALSE)

    # Save prediction raster
    pred_file <- file.path("outputs/predictions/kriging/carbon_stock_maps",
                           sprintf("carbon_stock_%s_%.0fcm.tif",
                                   gsub(" ", "_", stratum_name), depth))
    writeRaster(pred_raster, pred_file, overwrite = TRUE)

    # Save kriging variance raster (kriging uncertainty only)
    var_kriging_file <- file.path("outputs/predictions/kriging/variance_maps",
                                  sprintf("variance_kriging_%s_%.0fcm.tif",
                                          gsub(" ", "_", stratum_name), depth))
    writeRaster(var_raster, var_kriging_file, overwrite = TRUE)

    # Save combined variance raster (kriging + harmonization)
    var_combined_file <- file.path("outputs/predictions/kriging/variance_maps",
                                   sprintf("variance_combined_%s_%.0fcm.tif",
                                           gsub(" ", "_", stratum_name), depth))
    writeRaster(var_combined, var_combined_file, overwrite = TRUE)

    # Save kriging standard error
    se_kriging_file <- file.path("outputs/predictions/kriging/standard_error_maps",
                                 sprintf("se_kriging_%s_%.0fcm.tif",
                                         gsub(" ", "_", stratum_name), depth))
    writeRaster(se_kriging, se_kriging_file, overwrite = TRUE)

    # Save combined standard error (recommended for uncertainty reporting)
    se_combined_file <- file.path("outputs/predictions/kriging/standard_error_maps",
                                  sprintf("se_combined_%s_%.0fcm.tif",
                                          gsub(" ", "_", stratum_name), depth))
    writeRaster(se_combined, se_combined_file, overwrite = TRUE)
    
    # ========================================================================
    # CALCULATE AREA OF APPLICABILITY (AOA)
    # ========================================================================
    
    aoa_file <- NULL
    aoa_max_dist <- NA
    aoa_pct <- NA
    
    if (ENABLE_AOA && has_cast) {
      log_message("  Calculating Area of Applicability...")
      
      # Use variogram range as max distance if available
      vgm_range <- vgm_result$model$range[2]
      max_aoa_dist <- ifelse(!is.na(vgm_range) && vgm_range > 0,
                             vgm_range * 1.5,  # 1.5x variogram range
                             NULL)  # NULL = auto-calculate
      
      aoa_result <- calculate_spatial_aoa(
        pred_raster = pred_raster,
        sample_points = cores_sp,
        max_distance = max_aoa_dist
      )
      
      if (!is.null(aoa_result)) {
        # Save AOA raster
        aoa_file <- file.path("outputs/predictions/aoa",
                              sprintf("aoa_%s_%.0fcm.tif",
                                      gsub(" ", "_", stratum_name), depth))
        dir.create("outputs/predictions/aoa", recursive = TRUE, showWarnings = FALSE)
        writeRaster(aoa_result$aoa_raster, aoa_file, overwrite = TRUE)
        
        aoa_max_dist <- aoa_result$max_distance
        aoa_pct <- aoa_result$pct_in_aoa
        
        log_message(sprintf("  AOA: %.1f%% of area (max dist: %.0f m)",
                            aoa_pct, aoa_max_dist))
      }
    }
    
    # Store result
    result_key <- sprintf("%s_%.0fcm", gsub(" ", "_", stratum_name), depth)
    kriging_results[[result_key]] <- list(
      stratum = stratum_name,
      depth_cm = depth,
      n_samples = n_samples,
      prediction_file = pred_file,
      variance_kriging_file = var_kriging_file,
      variance_combined_file = var_combined_file,
      se_kriging_file = se_kriging_file,
      se_combined_file = se_combined_file,
      aoa_file = aoa_file,
      mean_prediction = mean(values(pred_raster), na.rm = TRUE),
      sd_prediction = sd(values(pred_raster), na.rm = TRUE),
      mean_se_combined = mean(values(se_combined), na.rm = TRUE),
      has_harmonization_uncertainty = has_uncertainties,
      harmonization_var_mean = harmonization_var_mean,
      aoa_max_distance = aoa_max_dist,
      aoa_pct_in = aoa_pct
    )
    
    log_message(sprintf("  Mean prediction: %.2f ± %.2f kg/m²",
                        kriging_results[[result_key]]$mean_prediction,
                        kriging_results[[result_key]]$sd_prediction))
  }
}

# ============================================================================
# DEPTH LAYER AGGREGATION: CARBON STOCK SUMMARIES
# ============================================================================

log_message(vm_label("\n=== VM0033 Layer Aggregation ===", "\n=== Depth Layer Aggregation ==="))
log_message("NOTE: Carbon stocks already calculated in Module 03 (kg/m² per depth interval)")
log_message(vm_label("      This section aggregates layers and converts to VM0033 reporting units (Mg C/ha)", "      This section aggregates layers and converts to reporting units (Mg C/ha)"))

# Create organised output directories for aggregated depth layer stocks
dir.create("outputs/predictions/kriging/aggregated_stocks", recursive = TRUE, showWarnings = FALSE)

# Process each stratum
for (stratum_name in unique(cores_standard$stratum)) {
  
  log_message(sprintf("\nAggregating stocks for: %s", stratum_name))
  
  # Initialize stack rasters for this stratum
  stock_layers <- list()
  
  # Process each VM0033 layer
  for (i in 1:nrow(VM0033_DEPTH_INTERVALS)) {
    
    depth_top <- VM0033_DEPTH_INTERVALS$depth_top[i]
    depth_bottom <- VM0033_DEPTH_INTERVALS$depth_bottom[i]
    depth_midpoint <- VM0033_DEPTH_INTERVALS$depth_midpoint[i]
    thickness_cm <- VM0033_DEPTH_INTERVALS$thickness_cm[i]
    
    log_message(sprintf("  Layer %d: %d-%d cm (midpoint %.1f cm, thickness %d cm)",
                        i, depth_top, depth_bottom, depth_midpoint, thickness_cm))
    
    # Load carbon stock prediction raster for this depth midpoint
    # Note: Carbon stocks are already calculated in Module 03 as kg/m² per depth increment
    stock_file_kgm2 <- file.path("outputs/predictions/kriging/carbon_stock_maps",
                                 sprintf("carbon_stock_%s_%.0fcm.tif",
                                         gsub(" ", "_", stratum_name),
                                         depth_midpoint))
    
    if (!file.exists(stock_file_kgm2)) {
      log_message(sprintf("    Carbon stock file not found: %s - skipping",
                          basename(stock_file_kgm2)), "WARNING")
      next
    }
    
    stock_raster_kgm2 <- rast(stock_file_kgm2)
    
    # Convert carbon stock from kg/m² to Mg C/ha for VM0033 reporting
    # Conversion: 1 kg/m² = 10 Mg/ha
    stock_raster <- stock_raster_kgm2 * 10
    
    # Save layer stock
    stock_file <- file.path("outputs/predictions/kriging/aggregated_stocks",
                            sprintf("stock_layer%d_%s_%d-%dcm.tif",
                                    i, gsub(" ", "_", stratum_name),
                                    depth_top, depth_bottom))
    writeRaster(stock_raster, stock_file, overwrite = TRUE)
    
    stock_layers[[paste0("layer_", i)]] <- stock_raster
    
    mean_stock <- mean(values(stock_raster), na.rm = TRUE)
    log_message(sprintf("    Mean stock: %.2f Mg C/ha", mean_stock))
    
    # Propagate uncertainty to stock (Mg C/ha)
    # Carbon stock SE in kg/m² from kriging
    se_file_kgm2 <- file.path("outputs/predictions/kriging/standard_error_maps",
                              sprintf("se_combined_%s_%.0fcm.tif",
                                      gsub(" ", "_", stratum_name),
                                      depth_midpoint))
    
    if (file.exists(se_file_kgm2)) {
      se_raster_kgm2 <- rast(se_file_kgm2)

      # Convert SE from kg/m² to Mg C/ha: multiply by 10
      stock_se_raster <- se_raster_kgm2 * 10

      # Save stock SE
      stock_se_file <- file.path("outputs/predictions/kriging/aggregated_stocks",
                                 sprintf("stock_se_layer%d_%s_%d-%dcm.tif",
                                         i, gsub(" ", "_", stratum_name),
                                         depth_top, depth_bottom))
      writeRaster(stock_se_raster, stock_se_file, overwrite = TRUE)
      
      mean_stock_se <- mean(values(stock_se_raster), na.rm = TRUE)
      log_message(sprintf("    Mean stock SE: %.2f Mg C/ha", mean_stock_se))
    }
  }
  
  # Calculate cumulative stocks if we have all layers
  if (length(stock_layers) > 0) {
    
    # Total stock to 1m depth (sum of all layers)
    if (length(stock_layers) == 4) {
      total_stock <- Reduce(`+`, stock_layers)
      
      total_file <- file.path("outputs/predictions/kriging/aggregated_stocks",
                              sprintf("stock_total_0-100cm_%s.tif",
                                      gsub(" ", "_", stratum_name)))
      writeRaster(total_stock, total_file, overwrite = TRUE)
      
      mean_total <- mean(values(total_stock), na.rm = TRUE)
      log_message(sprintf("  Total stock (0-100 cm): %.2f Mg C/ha", mean_total))
    }
    
    # Cumulative stocks for common reporting depths
    # 0-30 cm (top 2 layers)
    if (length(stock_layers) >= 2) {
      stock_0_30 <- stock_layers[[1]] + stock_layers[[2]]
      
      stock_file <- file.path("outputs/predictions/kriging/aggregated_stocks",
                              sprintf("stock_cumulative_0-30cm_%s.tif",
                                      gsub(" ", "_", stratum_name)))
      writeRaster(stock_0_30, stock_file, overwrite = TRUE)
      
      mean_030 <- mean(values(stock_0_30), na.rm = TRUE)
      log_message(sprintf("  Cumulative stock (0-30 cm): %.2f Mg C/ha", mean_030))
    }
    
    # 0-50 cm (top 3 layers)
    if (length(stock_layers) >= 3) {
      stock_0_50 <- stock_layers[[1]] + stock_layers[[2]] + stock_layers[[3]]
      
      stock_file <- file.path("outputs/predictions/kriging/aggregated_stocks",
                              sprintf("stock_cumulative_0-50cm_%s.tif",
                                      gsub(" ", "_", stratum_name)))
      writeRaster(stock_0_50, stock_file, overwrite = TRUE)
      
      mean_050 <- mean(values(stock_0_50), na.rm = TRUE)
      log_message(sprintf("  Cumulative stock (0-50 cm): %.2f Mg C/ha", mean_050))
    }
  }
}

log_message(vm_label("\nVM0033 layer aggregation complete", "\nDepth layer aggregation complete"))

# ============================================================================
# KRIGING CARBON STOCK SUMMARY MAPS
# ============================================================================
# Generates per-stratum 4-panel maps (one panel per depth layer) and a
# total 0-100 cm map, matching the Transfer Learning visualisation style.
# Saved to: outputs/Basic_analysis/Spatial_Maps_by_Depth/
# ============================================================================

log_message("Generating kriging carbon stock summary maps...")

suppressPackageStartupMessages({
  if (!requireNamespace("viridis", quietly = TRUE)) {
    message("  Package 'viridis' not available — skipping summary maps. Install with: install.packages('viridis')")
  }
})

if (requireNamespace("viridis", quietly = TRUE)) {
  library(viridis)

  dir.create("outputs/Basic_analysis/Spatial_Maps_by_Depth", recursive = TRUE, showWarnings = FALSE)

  # Depth interval labels matching VM0033 standard intervals
  depth_labels <- c("A. 0-15 cm", "B. 15-30 cm", "C. 30-50 cm", "D. 50-100 cm")

  for (stratum_name in unique(cores_standard$stratum)) {

    stratum_token <- gsub(" ", "_", stratum_name)

    tryCatch({

      # ── Load all four depth layer rasters ──────────────────────────────────
      stock_rasters <- list()
      for (i in 1:nrow(VM0033_DEPTH_INTERVALS)) {
        depth_top    <- VM0033_DEPTH_INTERVALS$depth_top[i]
        depth_bottom <- VM0033_DEPTH_INTERVALS$depth_bottom[i]
        stock_path <- file.path("outputs/predictions/kriging/aggregated_stocks",
                                sprintf("stock_layer%d_%s_%d-%dcm.tif",
                                        i, stratum_token, depth_top, depth_bottom))
        if (file.exists(stock_path)) {
          stock_rasters[[i]] <- rast(stock_path)
        } else {
          stock_rasters[[i]] <- NULL
        }
      }

      available <- !sapply(stock_rasters, is.null)
      if (sum(available) == 0) {
        log_message(sprintf("  Skipping map for %s — no aggregated stock rasters found", stratum_name), "WARNING")
        next
      }

      # ── Shared colour scale (2nd–98th percentile across all layers) ────────
      all_vals <- unlist(lapply(stock_rasters[available], function(r) {
        v <- terra::values(r, na.rm = TRUE)
        if (length(v) > 5e5) sample(v, 5e5) else v
      }))
      stock_limits <- quantile(all_vals, probs = c(0.02, 0.98), na.rm = TRUE)

      # ── 4-panel depth map ──────────────────────────────────────────────────
      map4_path <- file.path("outputs/Basic_analysis/Spatial_Maps_by_Depth",
                             sprintf("kriging_carbon_stock_%s_by_depth.png", stratum_token))

      png(map4_path, width = 1400, height = 1200, res = 150)
      par(mfrow = c(2, 2), mar = c(2, 2, 3, 5.5), oma = c(0, 0, 3, 0))

      for (i in 1:4) {
        panel_label <- if (i <= length(depth_labels)) depth_labels[i] else sprintf("Layer %d", i)
        if (!is.null(stock_rasters[[i]])) {
          terra::plot(stock_rasters[[i]],
                      main  = panel_label,
                      col   = viridis(100),
                      range = stock_limits,
                      axes  = FALSE, box = FALSE,
                      plg   = list(title = "Mg C/ha", cex = 0.85))
        } else {
          plot.new()
          title(main = sprintf("%s\n(no data)", panel_label), cex.main = 0.9)
        }
      }

      mtext(sprintf("Kriging Carbon Stock — %s", stratum_name),
            outer = TRUE, cex = 1.2, font = 2)
      dev.off()
      log_message(sprintf("  Saved 4-panel depth map: %s", basename(map4_path)))

      # ── Total 0-100 cm map ─────────────────────────────────────────────────
      total_path <- file.path("outputs/predictions/kriging/aggregated_stocks",
                              sprintf("stock_total_0-100cm_%s.tif", stratum_token))

      if (file.exists(total_path)) {
        total_rast <- rast(total_path)
        total_vals <- terra::values(total_rast, na.rm = TRUE)
        total_limits <- quantile(total_vals, probs = c(0.02, 0.98), na.rm = TRUE)

        total_map_path <- file.path("outputs/Basic_analysis/Spatial_Maps_by_Depth",
                                    sprintf("kriging_carbon_stock_%s_total_0-100cm.png", stratum_token))

        png(total_map_path, width = 900, height = 800, res = 150)
        par(mar = c(2, 2, 3.5, 6))
        terra::plot(total_rast,
                    main  = sprintf("Total Carbon Stock (0-100 cm)\n%s", stratum_name),
                    col   = viridis(100),
                    range = total_limits,
                    axes  = FALSE, box = FALSE,
                    plg   = list(title = "Mg C/ha"))
        dev.off()
        log_message(sprintf("  Saved total stock map:  %s", basename(total_map_path)))
      }

    }, error = function(e) {
      try(dev.off(), silent = TRUE)
      log_message(sprintf("  WARNING: Could not create map for %s — %s", stratum_name, e$message), "WARNING")
    })
  }

  log_message("Carbon stock summary maps complete")

} else {
  log_message("Skipping carbon stock maps — install 'viridis' package to enable", "WARNING")
}

# ============================================================================
# AOI OVERVIEW SUMMARY: RAW CORES / KRIGING MAP / DEPTH PROFILES
# ============================================================================
# Three-panel summary per stratum:
#   Panel A — Raw core locations coloured by total carbon stock (0-100 cm)
#   Panel B — Kriging interpolation (total 0-100 cm, Mg C/ha)
#   Panel C — Depth profile bar chart (harmonised mean per depth layer)
#
# Requires 'viridis'. Skipped gracefully if rasters are missing.
# Saved to: outputs/Basic_analysis/Spatial_Maps_by_Depth/
# ============================================================================

log_message("Generating AOI overview summary figures...")

if (requireNamespace("viridis", quietly = TRUE)) {

  # Compute per-core total carbon stock (sum across depth layers)
  core_totals_summary <- cores_standard %>%
    dplyr::group_by(core_id, stratum) %>%
    dplyr::summarise(
      longitude          = mean(longitude, na.rm = TRUE),
      latitude           = mean(latitude,  na.rm = TRUE),
      total_stock_kgm2   = sum(carbon_stock_kg_m2, na.rm = TRUE),
      n_depths           = dplyr::n(),
      .groups = "drop"
    )

  # Depth layer labels for bar chart
  depth_bar_labels <- if (exists("VM0033_DEPTH_INTERVALS")) {
    sprintf("%d-%d cm", VM0033_DEPTH_INTERVALS$depth_top, VM0033_DEPTH_INTERVALS$depth_bottom)
  } else {
    paste0(STANDARD_DEPTHS, " cm")
  }

  for (stratum_name in unique(cores_standard$stratum)) {

    stratum_token  <- gsub(" ", "_", stratum_name)
    stratum_cores  <- dplyr::filter(core_totals_summary, stratum == stratum_name)
    stratum_depths <- cores_standard %>%
      dplyr::filter(stratum == stratum_name) %>%
      dplyr::group_by(depth_cm_midpoint) %>%
      dplyr::summarise(
        mean_stock_mgha = mean(carbon_stock_kg_m2, na.rm = TRUE) * 10,  # kg/m² → Mg C/ha
        .groups = "drop"
      ) %>%
      dplyr::arrange(depth_cm_midpoint)

    if (nrow(stratum_cores) == 0) next

    tryCatch({

      out_path <- file.path("outputs/Basic_analysis/Spatial_Maps_by_Depth",
                            sprintf("aoi_overview_%s.png", stratum_token))

      # Load total kriging stock raster (if available)
      total_rast_path <- file.path("outputs/predictions/kriging/aggregated_stocks",
                                   sprintf("stock_total_0-100cm_%s.tif", stratum_token))
      has_kriging_rast <- file.exists(total_rast_path)
      if (has_kriging_rast) {
        total_rast <- rast(total_rast_path)
      }

      png(out_path, width = 1600, height = 600, res = 150)

      # Layout: Panel A and B are square spatial panels; Panel C is a bar chart
      layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 0.85))
      par(oma = c(0, 0, 2.5, 0))

      # ── Panel A: Raw core locations ─────────────────────────────────────
      par(mar = c(3, 3, 3, 1))

      # Colour-scale for raw point stocks
      pt_vals <- stratum_cores$total_stock_kgm2
      pt_cols_n <- 50
      pt_breaks <- seq(min(pt_vals, na.rm=TRUE), max(pt_vals, na.rm=TRUE) + 1e-9, length.out = pt_cols_n + 1)
      pt_pal <- viridis(pt_cols_n)
      pt_col_idx <- findInterval(pt_vals, pt_breaks, rightmost.closed = TRUE)
      pt_col_idx <- pmax(1, pmin(pt_col_idx, pt_cols_n))
      pt_colors  <- pt_pal[pt_col_idx]

      plot(stratum_cores$longitude, stratum_cores$latitude,
           col = pt_colors, pch = 21, bg = pt_colors, cex = 1.8,
           xlab = "Longitude", ylab = "Latitude",
           main = "A. Raw Core Locations\n(total stock, kg/m²)",
           cex.main = 0.9, cex.axis = 0.7, cex.lab = 0.75)
      # Colour-bar legend (vertical strip)
      legend_x <- grconvertX(0.85, "npc", "user")
      legend_y_bot <- grconvertY(0.1, "npc", "user")
      legend_y_top <- grconvertY(0.9, "npc", "user")
      legend_strip_x <- c(legend_x, legend_x + diff(range(stratum_cores$longitude, na.rm=TRUE)) * 0.04)
      n_leg <- 10
      leg_breaks <- seq(min(pt_vals), max(pt_vals), length.out = n_leg + 1)
      for (j in seq_len(n_leg)) {
        rect(legend_strip_x[1], legend_y_bot + (j-1)/n_leg * (legend_y_top - legend_y_bot),
             legend_strip_x[2], legend_y_bot + j/n_leg     * (legend_y_top - legend_y_bot),
             col = viridis(n_leg)[j], border = NA)
      }

      # ── Panel B: Kriging prediction raster ──────────────────────────────
      par(mar = c(3, 1, 3, 5.5))
      if (has_kriging_rast) {
        k_vals    <- terra::values(total_rast, na.rm = TRUE)
        k_limits  <- quantile(k_vals, probs = c(0.02, 0.98), na.rm = TRUE)
        terra::plot(total_rast,
                    main  = "B. Kriging Interpolation\n(0-100 cm, Mg C/ha)",
                    col   = viridis(100),
                    range = k_limits,
                    axes  = FALSE, box = FALSE,
                    plg   = list(title = "Mg C/ha", cex = 0.8),
                    cex.main = 0.9)
      } else {
        plot.new()
        title(main = "B. Kriging Interpolation\n(raster not available)", cex.main = 0.9)
        text(0.5, 0.5, "Run kriging module\nto generate raster", cex = 0.9, col = "grey40")
      }

      # ── Panel C: Harmonised depth profile bar chart ──────────────────────
      par(mar = c(5, 3.5, 3, 1.5))
      if (nrow(stratum_depths) > 0) {
        bar_vals   <- stratum_depths$mean_stock_mgha
        bar_labels <- if (length(depth_bar_labels) == nrow(stratum_depths)) {
          depth_bar_labels
        } else {
          paste0(stratum_depths$depth_cm_midpoint, " cm")
        }
        bar_colors <- viridis(nrow(stratum_depths), option = "D", direction = -1)

        bp <- barplot(bar_vals,
                      names.arg  = bar_labels,
                      col        = bar_colors,
                      horiz      = FALSE,
                      las        = 2,
                      cex.names  = 0.75,
                      cex.axis   = 0.75,
                      ylab       = "Mean Carbon Stock (Mg C/ha)",
                      main       = "C. Harmonised Depth Profile",
                      cex.main   = 0.9,
                      cex.lab    = 0.8,
                      border     = "white")
        # Value labels above bars
        text(bp, bar_vals + max(bar_vals, na.rm=TRUE) * 0.03,
             labels = sprintf("%.1f", bar_vals), cex = 0.7, col = "grey20")
      } else {
        plot.new()
        title(main = "C. Harmonised Depth Profile\n(no data)", cex.main = 0.9)
      }

      mtext(sprintf("Site Overview — %s", stratum_name), outer = TRUE, cex = 1.2, font = 2)
      dev.off()
      log_message(sprintf("  Saved AOI overview: %s", basename(out_path)))

    }, error = function(e) {
      try(dev.off(), silent = TRUE)
      log_message(sprintf("  WARNING: Could not create overview for %s — %s",
                          stratum_name, e$message), "WARNING")
    })
  }

  log_message("AOI overview summary figures complete")

} else {
  log_message("Skipping AOI overview figures — install 'viridis' package to enable", "WARNING")
}

# ============================================================================
# SAVE VARIOGRAM MODELS AND CV RESULTS
# ============================================================================

log_message("Saving variogram models and CV results...")

# Save variogram models
saveRDS(variogram_models, "outputs/models/kriging/variogram_models.rds")

# Save CV results
write.csv(cv_results_all, "diagnostics/crossvalidation/kriging_cv_results.csv", 
          row.names = FALSE)
saveRDS(cv_results_all, "diagnostics/crossvalidation/kriging_cv_results.rds")

# Save kriging results summary
saveRDS(kriging_results, "outputs/models/kriging/kriging_results_summary.rds")

log_message("Saved models and diagnostics")

# ============================================================================
# CREATE CV SUMMARY PLOTS
# ============================================================================

if (nrow(cv_results_all) > 0) {
  
  log_message("Creating CV summary plots...")
  
  suppressPackageStartupMessages({
    library(ggplot2)
    library(tidyr)
  })
  
  # CV RMSE by stratum and depth
  p_cv_rmse <- ggplot(cv_results_all, aes(x = factor(depth_cm), y = cv_rmse, 
                                          fill = stratum)) +
    geom_col(position = "dodge", alpha = 0.7) +
    scale_fill_manual(values = STRATUM_COLORS) +
    labs(
      title = "Kriging Cross-Validation RMSE",
      subtitle = "By stratum and depth (carbon stock predictions)",
      x = "Depth (cm)",
      y = "CV RMSE (kg/m²)",
      fill = "Stratum"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("diagnostics/crossvalidation/cv_rmse_by_stratum_depth.png",
         p_cv_rmse, width = 10, height = 6, dpi = 300)
  
  # CV R² by stratum and depth
  p_cv_r2 <- ggplot(cv_results_all, aes(x = factor(depth_cm), y = cv_r2, 
                                        fill = stratum)) +
    geom_col(position = "dodge", alpha = 0.7) +
    geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
    scale_fill_manual(values = STRATUM_COLORS) +
    labs(
      title = "Kriging Cross-Validation R²",
      subtitle = "By stratum and depth (red line = 0.7 threshold)",
      x = "Depth (cm)",
      y = "CV R²",
      fill = "Stratum"
    ) +
    ylim(0, 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("diagnostics/crossvalidation/cv_r2_by_stratum_depth.png",
         p_cv_r2, width = 10, height = 6, dpi = 300)
  
  # ==========================================================================
  # MODEL COMPARISON TABLE
  # ==========================================================================
  
  log_message("Creating model comparison table...")
  
  # Summary table by stratum and depth
  comparison_table <- cv_results_all %>%
    group_by(stratum, depth_cm, kriging_type) %>%
    summarise(
      n = mean(n_samples),
      RMSE = mean(cv_rmse, na.rm = TRUE),
      MAE = mean(cv_mae, na.rm = TRUE),
      R2 = mean(cv_r2, na.rm = TRUE),
      Model = first(model_type),
      .groups = "drop"
    ) %>%
    arrange(stratum, depth_cm, kriging_type)
  
  # Save as CSV
  write.csv(comparison_table,
            "diagnostics/crossvalidation/model_comparison_table.csv",
            row.names = FALSE)
  
  # Print to console
  cat("\n=== MODEL COMPARISON TABLE ===\n")
  print(comparison_table, n = 100)
  
  # If comparing OK vs SK, create comparison plot
  if (COMPARE_KRIGING_METHODS && "simple" %in% cv_results_all$kriging_type) {
    
    # Comparison plot: OK vs SK performance
    comparison_wide <- cv_results_all %>%
      select(stratum, depth_cm, kriging_type, cv_rmse, cv_r2) %>%
      pivot_wider(names_from = kriging_type,
                  values_from = c(cv_rmse, cv_r2),
                  names_sep = "_")
    
    if (all(c("cv_rmse_ordinary", "cv_rmse_simple") %in% names(comparison_wide))) {
      
      # RMSE comparison
      p_method_comp <- ggplot(comparison_wide,
                              aes(x = cv_rmse_ordinary, y = cv_rmse_simple)) +
        geom_point(aes(color = stratum), size = 3) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        geom_text(aes(label = paste0(depth_cm, "cm")),
                  hjust = -0.2, size = 3) +
        scale_color_manual(values = STRATUM_COLORS) +
        labs(
          title = "Kriging Method Comparison",
          subtitle = "Ordinary vs Simple Kriging CV RMSE",
          x = "Ordinary Kriging RMSE (kg/m²)",
          y = "Simple Kriging RMSE (kg/m²)",
          color = "Stratum",
          caption = "Points below line indicate Simple Kriging performs better"
        ) +
        theme_minimal() +
        theme(legend.position = "right")
      
      ggsave("diagnostics/crossvalidation/method_comparison_ok_vs_sk.png",
             p_method_comp, width = 8, height = 6, dpi = 300)
      
      log_message("Saved method comparison plot")
    }
  }
  
  log_message("Saved CV plots and comparison tables")
}

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("MODULE 04 COMPLETE\n")
cat("========================================\n\n")

cat("Kriging Summary:\n")
cat("----------------------------------------\n")
cat(sprintf("Depths processed: %d\n", length(STANDARD_DEPTHS)))
cat(sprintf("Strata processed: %d\n", length(strata)))
cat(sprintf("Total predictions created: %d\n", length(kriging_results)))
cat(sprintf("Variogram models saved: %d\n", length(variogram_models)))

if (nrow(cv_results_all) > 0) {
  cat("\nCross-Validation Performance:\n")
  cat("----------------------------------------\n")
  
  cv_summary <- cv_results_all %>%
    group_by(stratum) %>%
    summarise(
      mean_rmse = mean(cv_rmse),
      mean_r2 = mean(cv_r2),
      .groups = "drop"
    )
  
  for (i in 1:nrow(cv_summary)) {
    cat(sprintf("%s: RMSE=%.2f, R²=%.3f\n",
                cv_summary$stratum[i],
                cv_summary$mean_rmse[i],
                cv_summary$mean_r2[i]))
  }
}

cat("\nOutputs:\n")
cat("  Carbon Stock Maps (kg/m²):        outputs/predictions/kriging/carbon_stock_maps/\n")
cat("  Standard Error Maps:              outputs/predictions/kriging/standard_error_maps/\n")
cat("  Variance Maps:                    outputs/predictions/kriging/variance_maps/\n")
cat(sprintf("  %s (Mg C/ha): outputs/predictions/kriging/aggregated_stocks/\n", vm_label("VM0033 Aggregated Stocks", "Depth Layer Stocks")))
cat("  Area of Applicability:            outputs/predictions/aoa/\n")
cat("  Variogram Models:                 outputs/models/kriging/\n")
cat("  Diagnostics:                      diagnostics/variograms/ and diagnostics/crossvalidation/\n")

cat("\nKey Changes:\n")
cat("  ✓ Modeling carbon stocks (kg/m²) directly from Module 03 harmonization\n")
cat("  ✓ No SOC→stock conversion needed (stocks pre-calculated from SOC + BD)\n")
cat(vm_label("  ✓ VM0033 reporting stocks converted to Mg C/ha (multiply by 10)\n", "  ✓ Depth layer stocks converted to Mg C/ha (multiply by 10)\n"))

cat("\nNext steps:\n")
cat("  1. Review variogram plots in diagnostics/variograms/\n")
cat("  2. Check CV results in diagnostics/crossvalidation/\n")
cat("  3. Review carbon stock prediction maps by stratum and depth\n")
cat("  4. Run: source('05_raster_predictions_rf_bluecarbon.R')\n\n")

# ============================================================================
# QUICK-ACCESS COPY: BASIC ANALYSIS OUTPUTS (SPATIAL MAPS BY DEPTH)
# ============================================================================

log_message("Copying key Module 04 outputs to outputs/Basic_analysis/Spatial_Maps_by_Depth...")
dir.create("outputs/Basic_analysis", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Basic_analysis/Spatial_Maps_by_Depth", recursive = TRUE, showWarnings = FALSE)

format_stratum_label <- function(stratum_token) {
  gsub("_", " ", stratum_token)
}

get_depth_interval_label <- function(depth_cm) {
  if (exists("VM0033_DEPTH_INTERVALS") && all(c("depth_midpoint", "depth_top", "depth_bottom") %in% names(VM0033_DEPTH_INTERVALS))) {
    depth_index <- which.min(abs(VM0033_DEPTH_INTERVALS$depth_midpoint - depth_cm))
    depth_top <- VM0033_DEPTH_INTERVALS$depth_top[depth_index]
    depth_bottom <- VM0033_DEPTH_INTERVALS$depth_bottom[depth_index]
    return(sprintf("%d_to_%dcm", depth_top, depth_bottom))
  }
  sprintf("%.0fcm", depth_cm)
}

build_basic_spatial_name <- function(src_file) {
  src_name <- tools::file_path_sans_ext(basename(src_file))
  src_ext <- tools::file_ext(src_file)

  parsed <- regexec("^(carbon_stock|se_combined|se_kriging|variance_combined|variance_kriging)_(.+)_([0-9]+)cm$", src_name)
  parsed_parts <- regmatches(src_name, parsed)[[1]]

  if (length(parsed_parts) == 4) {
    metric_key <- parsed_parts[2]
    stratum_key <- parsed_parts[3]
    depth_cm <- as.numeric(parsed_parts[4])
    depth_label <- get_depth_interval_label(depth_cm)
    stratum_label <- gsub(" ", "_", format_stratum_label(stratum_key))

    metric_label <- switch(
      metric_key,
      carbon_stock = "Kriging_Carbon_Stock_Map",
      se_combined = "Combined_Standard_Error_Uncertainty_Map",
      se_kriging = "Kriging_Standard_Error_Uncertainty_Map",
      variance_combined = "Combined_Variance_Uncertainty_Map",
      variance_kriging = "Kriging_Variance_Uncertainty_Map",
      metric_key
    )

    return(sprintf("%s_%s_%s.%s", metric_label, stratum_label, depth_label, src_ext))
  }

  sprintf("Spatial_Output_%s.%s", src_name, src_ext)
}

basic_spatial_sources <- c(
  list.files("outputs/predictions/kriging/carbon_stock_maps",  pattern = "\\.tif$", full.names = TRUE),
  list.files("outputs/predictions/kriging/standard_error_maps", pattern = "\\.tif$", full.names = TRUE),
  list.files("outputs/predictions/kriging/variance_maps",       pattern = "\\.tif$", full.names = TRUE)
)

if (length(basic_spatial_sources) == 0) {
  log_message("No kriging spatial outputs found for Basic_analysis copy step", "WARNING")
} else {
  for (src in basic_spatial_sources) {
    dst_name <- build_basic_spatial_name(src)
    dst <- file.path("outputs/Basic_analysis/Spatial_Maps_by_Depth", dst_name)
    file.copy(src, dst, overwrite = TRUE)
    log_message(sprintf("Copied to Basic_analysis/Spatial_Maps_by_Depth: %s", basename(dst)))
  }
}

log_message("=== MODULE 04 COMPLETE ===")
