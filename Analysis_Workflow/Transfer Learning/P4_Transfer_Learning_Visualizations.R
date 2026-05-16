# ============================================================================
# MODULE 06: VISUALIZE TRANSFER LEARNING PREDICTIONS (REVISED)
# ============================================================================
# PURPOSE:   Load and visualize pre-generated prediction rasters from Module 05.
#            Creates comparison plots showing Global Prior, Local-Only,
#            Transfer Learning Final, and Difference maps.
#
# INPUTS:    
#   - outputs/carbon_stocks/predictions/depth_*_predictions.tif
#     (4-band GeoTIFFs: Global_Prior, Local_Only, Transfer_Final, Difference)
#   - diagnostics/transfer/transfer_learning_validation.csv (optional)
#
# OUTPUTS:   
#   - outputs/carbon_stocks/visualizations/comparison_plot_depth_*.png
#   - outputs/carbon_stocks/visualizations/all_depths_overview.png
#   - outputs/carbon_stocks/visualizations/prediction_statistics.csv
#
# USAGE:     Source this script from your project root directory, or set
#            PROJECT_ROOT below to an absolute path.
#
# NOTE:      This script does NO modeling. Run Module 05 first.
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
  library(terra)
  library(viridis)
  library(scico)
})

# ============================================================================
# CONFIGURATION
# ============================================================================

# Set project root (use relative path "." if running from project directory)
# Change this to an absolute path if needed for your environment
PROJECT_ROOT <- "."

# Depth intervals to visualize (must match Module 05 outputs)
DEPTHS <- c(7.5, 22.5, 40, 75)

# Expected band names in prediction rasters
EXPECTED_BANDS <- c("Global_Prior", "Local_Only", "Transfer_Final", "Difference")

# Quantile thresholds for color scaling
QUANTILE_LOW <- 0.02
QUANTILE_HIGH <- 0.98

# Interactive mode (set FALSE for batch processing)
INTERACTIVE <- interactive()

# Resolve namespace conflicts explicitly
select <- dplyr::select
filter <- dplyr::filter

# ============================================================================
# SETUP PATHS
# ============================================================================

INPUT_DIR <- file.path(PROJECT_ROOT, "outputs/carbon_stocks/predictions")
OUTPUT_DIR <- file.path(PROJECT_ROOT, "outputs/carbon_stocks/visualizations")
STATS_PATH <- file.path(PROJECT_ROOT, "diagnostics/transfer/transfer_learning_validation.csv")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("\n==================================================\n")
cat("   VISUALIZING TRANSFER LEARNING PREDICTIONS\n")
cat("==================================================\n")
cat(sprintf("Project root: %s\n", normalizePath(PROJECT_ROOT, mustWork = FALSE)))
cat(sprintf("Input directory: %s\n", INPUT_DIR))
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))

# ============================================================================
# 1. VALIDATE INPUT RASTERS
# ============================================================================

cat("\n--- Checking Input Rasters ---\n")

available_depths <- c()
raster_info <- list()

for (d in DEPTHS) {
  input_path <- file.path(INPUT_DIR, sprintf("depth_%s_predictions.tif", d))
  
  if (!file.exists(input_path)) {
    cat(sprintf("  [MISSING] depth_%s_predictions.tif\n", d))
    next
  }
  
  # Load and validate
  tryCatch({
    r <- rast(input_path)
    n_bands <- nlyr(r)
    band_names <- names(r)
    
    if (n_bands < 4) {
      cat(sprintf("  [WARNING] depth_%s: only %d bands (expected 4)\n", d, n_bands))
    }
    
    available_depths <- c(available_depths, d)
    raster_info[[as.character(d)]] <- list(
      path = input_path,
      bands = n_bands,
      names = band_names,
      extent = ext(r),
      crs = crs(r, describe = TRUE)$name
    )
    
    cat(sprintf("  [OK] depth_%s: %d bands, CRS: %s\n", 
                d, n_bands, raster_info[[as.character(d)]]$crs))
    
  }, error = function(e) {
    cat(sprintf("  [ERROR] depth_%s: %s\n", d, e$message))
  })
}

if (length(available_depths) == 0) {
  stop("\nNo valid prediction rasters found in: ", INPUT_DIR, "\n\n",
       "To generate these rasters, you must:\n",
       "  1. Have covariate .tif files in your 'covariates/' folder\n",
       "  2. Run Module 05 (transfer learning) successfully\n\n",
       "Module 05 will generate prediction rasters only if covariate rasters exist.\n",
       "If you only have point data (no rasters), spatial visualization is not possible.\n")
}

cat(sprintf("\nFound %d of %d expected depths: %s\n", 
            length(available_depths), 
            length(DEPTHS),
            paste(available_depths, collapse = ", ")))

# ============================================================================
# 2. LOAD VALIDATION STATISTICS (OPTIONAL)
# ============================================================================

cat("\n--- Loading Validation Statistics ---\n")

stats_df <- NULL

if (file.exists(STATS_PATH)) {
  tryCatch({
    stats_df <- read_csv(STATS_PATH, show_col_types = FALSE)
    
    # Standardize column names (handle both old and new naming conventions)
    if ("depth" %in% names(stats_df) && !"depth_cm" %in% names(stats_df)) {
      stats_df <- stats_df %>% rename(depth_cm = depth)
    }
    
    cat(sprintf("  Loaded statistics for %d depths\n", nrow(stats_df)))
    cat(sprintf("  Available columns: %s\n", paste(names(stats_df), collapse = ", ")))
    
  }, error = function(e) {
    cat(sprintf("  Warning: Could not load stats file: %s\n", e$message))
    stats_df <- NULL
  })
} else {
  cat("  No validation statistics file found (plots will omit R²/RMSE)\n")
}

# ============================================================================
# 3. HELPER FUNCTIONS
# ============================================================================

#' Safely compute quantile range from raster using terra::global
#' Avoids loading full raster into memory
compute_raster_quantiles <- function(r, probs = c(0.02, 0.98)) {
  # terra::global with quantile function
  # Fall back to minmax if quantile fails
  tryCatch({
    vals <- values(r, na.rm = TRUE)
    if (length(vals) > 1e6) {
      # Sample for very large rasters
      vals <- sample(vals, 1e6)
    }
    quantile(vals, probs = probs, na.rm = TRUE)
  }, error = function(e) {
    mm <- minmax(r)
    c(mm[1], mm[2])
  })
}

#' Create a 4-panel comparison plot for one depth
#' Returns TRUE on success, FALSE on failure
create_depth_comparison <- function(depth, input_dir, output_dir, stats_df = NULL,
                                    save_png = TRUE, display = FALSE) {
  
  input_path <- file.path(input_dir, sprintf("depth_%s_predictions.tif", depth))
  
  if (!file.exists(input_path)) {
    warning(sprintf("Raster not found for depth %s", depth))
    return(FALSE)
  }
  
  # Load raster stack
  pred_stack <- rast(input_path)
  
  # Extract layers by name if possible, otherwise by position
  layer_names <- names(pred_stack)
  
  get_layer <- function(pattern, position) {
    idx <- grep(pattern, layer_names, ignore.case = TRUE)
    if (length(idx) > 0) {
      return(pred_stack[[idx[1]]])
    } else if (position <= nlyr(pred_stack)) {
      return(pred_stack[[position]])
    } else {
      return(NULL)
    }
  }
  
  map_global <- get_layer("global", 1)
  map_local <- get_layer("local", 2)
  map_transfer <- get_layer("transfer|final", 3)
  map_difference <- get_layer("diff", 4)
  
  if (is.null(map_transfer)) {
    warning(sprintf("Could not identify transfer layer for depth %s", depth))
    return(FALSE)
  }
  
  # Compute color scales
  # Carbon stock layers share a common scale
  stock_layers <- list(map_global, map_local, map_transfer)
  stock_layers <- stock_layers[!sapply(stock_layers, is.null)]
  
  all_stock_vals <- unlist(lapply(stock_layers, function(r) {
    v <- values(r, na.rm = TRUE)
    if (length(v) > 5e5) sample(v, 5e5) else v
  }))
  
  stock_limits <- quantile(all_stock_vals, probs = c(QUANTILE_LOW, QUANTILE_HIGH), na.rm = TRUE)
  
  # Difference layer uses symmetric diverging scale
  diff_limits <- NULL
  if (!is.null(map_difference)) {
    diff_vals <- values(map_difference, na.rm = TRUE)
    if (length(diff_vals) > 5e5) diff_vals <- sample(diff_vals, 5e5)
    diff_max <- max(abs(quantile(diff_vals, probs = c(0.05, 0.95), na.rm = TRUE)))
    diff_limits <- c(-diff_max, diff_max)
  }
  
  # Get validation stats for title
  title_stats <- ""
  if (!is.null(stats_df) && depth %in% stats_df$depth_cm) {
    ds <- stats_df %>% filter(depth_cm == depth)
    
    # Handle different column naming conventions
    r2_val <- ds$cv_r2[1] %||% ds$transfer_r2[1] %||% NA
    rmse_val <- ds$cv_rmse[1] %||% ds$transfer_rmse[1] %||% NA
    
    if (!is.na(r2_val) && !is.na(rmse_val)) {
      title_stats <- sprintf(" | CV R² = %.2f | RMSE = %.2f kg/m²", r2_val, rmse_val)
    }
  }
  
  main_title <- sprintf("Carbon Stock Predictions: %.1f cm depth%s", depth, title_stats)
  
  # Calculate summary statistics
  summary_stats <- data.frame(
    depth_cm = depth,
    global_mean = if (!is.null(map_global)) global(map_global, "mean", na.rm = TRUE)[1, 1] else NA,
    local_mean = if (!is.null(map_local)) global(map_local, "mean", na.rm = TRUE)[1, 1] else NA,
    transfer_mean = global(map_transfer, "mean", na.rm = TRUE)[1, 1],
    transfer_sd = global(map_transfer, "sd", na.rm = TRUE)[1, 1],
    diff_mean = if (!is.null(map_difference)) global(map_difference, "mean", na.rm = TRUE)[1, 1] else NA,
    diff_sd = if (!is.null(map_difference)) global(map_difference, "sd", na.rm = TRUE)[1, 1] else NA
  )
  
  # Create plot
  plot_func <- function() {
    # Determine layout based on available layers
    n_panels <- sum(!sapply(list(map_global, map_local, map_transfer, map_difference), is.null))
    
    if (n_panels == 4) {
      par(mfrow = c(2, 2), mar = c(2, 2, 3, 5), oma = c(0, 0, 3, 0))
    } else if (n_panels == 3) {
      par(mfrow = c(1, 3), mar = c(2, 2, 3, 5), oma = c(0, 0, 3, 0))
    } else {
      par(mfrow = c(1, n_panels), mar = c(2, 2, 3, 5), oma = c(0, 0, 3, 0))
    }
    
    # Panel A: Global Prior
    if (!is.null(map_global)) {
      plot(map_global,
           main = "A. Global Prior",
           col = viridis(100),
           range = stock_limits,
           axes = FALSE, box = FALSE,
           plg = list(title = "kg/m²"))
    }
    
    # Panel B: Local-Only
    if (!is.null(map_local)) {
      plot(map_local,
           main = "B. Local-Only Benchmark",
           col = viridis(100),
           range = stock_limits,
           axes = FALSE, box = FALSE,
           plg = list(title = "kg/m²"))
    }
    
    # Panel C: Transfer Final
    plot(map_transfer,
         main = "C. Transfer Learning (Final)",
         col = viridis(100),
         range = stock_limits,
         axes = FALSE, box = FALSE,
         plg = list(title = "kg/m²"))
    
    # Panel D: Difference
    if (!is.null(map_difference) && !is.null(diff_limits)) {
      plot(map_difference,
           main = "D. Difference (Transfer - Local)",
           col = scico(100, palette = "vik"),
           range = diff_limits,
           axes = FALSE, box = FALSE,
           plg = list(title = "kg/m²"))
    }
    
    mtext(main_title, outer = TRUE, cex = 1.2, font = 2)
  }
  
  # Save to PNG
  if (save_png) {
    png_path <- file.path(output_dir, sprintf("comparison_plot_depth_%s.png", depth))
    
    tryCatch({
      png(png_path, width = 1400, height = 1200, res = 150)
      plot_func()
      dev.off()
      cat(sprintf("  Saved: %s\n", basename(png_path)))
    }, error = function(e) {
      try(dev.off(), silent = TRUE)
      warning(sprintf("Failed to save PNG for depth %s: %s", depth, e$message))
    })
  }
  
  # Display in R graphics device
  if (display && INTERACTIVE) {
    plot_func()
  }
  
  return(summary_stats)
}

# Null coalescing operator (for R < 4.4 compatibility)
`%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

# ============================================================================
# 4. GENERATE INDIVIDUAL DEPTH PLOTS
# ============================================================================

cat("\n--- Generating Comparison Plots ---\n")

summary_list <- list()

for (d in available_depths) {
  cat(sprintf("\nProcessing depth: %.1f cm\n", d))
  
  result <- create_depth_comparison(
    depth = d,
    input_dir = INPUT_DIR,
    output_dir = OUTPUT_DIR,
    stats_df = stats_df,
    save_png = TRUE,
    display = FALSE
  )
  
  if (is.data.frame(result)) {
    summary_list[[as.character(d)]] <- result
  }
}

# ============================================================================
# 5. CREATE MULTI-DEPTH OVERVIEW
# ============================================================================

cat("\n--- Creating Multi-Depth Overview ---\n")

# Load all transfer final predictions
transfer_rasters <- list()

for (d in available_depths) {
  input_path <- file.path(INPUT_DIR, sprintf("depth_%s_predictions.tif", d))
  
  tryCatch({
    pred_stack <- rast(input_path)
    
    # Get transfer layer (usually band 3)
    layer_names <- names(pred_stack)
    transfer_idx <- grep("transfer|final", layer_names, ignore.case = TRUE)
    
    if (length(transfer_idx) > 0) {
      transfer_layer <- pred_stack[[transfer_idx[1]]]
    } else {
      transfer_layer <- pred_stack[[min(3, nlyr(pred_stack))]]
    }
    
    names(transfer_layer) <- sprintf("%.1f cm", d)
    transfer_rasters[[as.character(d)]] <- transfer_layer
    
  }, error = function(e) {
    cat(sprintf("  Warning: Could not load depth %s: %s\n", d, e$message))
  })
}

if (length(transfer_rasters) > 0) {
  
  # Compute common color scale across all depths
  all_vals <- unlist(lapply(transfer_rasters, function(r) {
    v <- values(r, na.rm = TRUE)
    if (length(v) > 5e5) sample(v, 5e5) else v
  }))
  
  common_limits <- quantile(all_vals, probs = c(QUANTILE_LOW, QUANTILE_HIGH), na.rm = TRUE)
  
  # Create overview plot
  overview_path <- file.path(OUTPUT_DIR, "all_depths_transfer_final.png")
  
  n_depths <- length(transfer_rasters)
  n_cols <- min(n_depths, 2)
  n_rows <- ceiling(n_depths / n_cols)
  
  tryCatch({
    png(overview_path, width = 700 * n_cols, height = 600 * n_rows, res = 150)
    
    par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 3, 5), oma = c(0, 0, 3, 0))
    
    for (i in seq_along(transfer_rasters)) {
      d <- names(transfer_rasters)[i]
      r <- transfer_rasters[[i]]
      
      plot(r,
           main = sprintf("Transfer Final: %s", names(r)),
           col = viridis(100),
           range = common_limits,
           axes = FALSE, box = FALSE,
           plg = list(title = "kg/m²"))
    }
    
    mtext("Transfer Learning Predictions by Depth", outer = TRUE, cex = 1.3, font = 2)
    
    dev.off()
    cat(sprintf("  Saved: %s\n", basename(overview_path)))
    
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    warning(sprintf("Failed to create overview: %s", e$message))
  })
}

# ============================================================================
# 6. EXPORT SUMMARY STATISTICS
# ============================================================================

cat("\n--- Exporting Summary Statistics ---\n")

if (length(summary_list) > 0) {
  summary_df <- bind_rows(summary_list)
  
  # Merge with validation statistics if available
  if (!is.null(stats_df)) {
    # Select relevant columns, handling different naming conventions
    stats_cols <- intersect(
      names(stats_df),
      c("depth_cm", "cv_r2", "cv_rmse", "transfer_r2", "transfer_rmse",
        "bias_correction", "bias_se", "residual_sd", "n_local", "n_global")
    )
    
    if (length(stats_cols) > 1) {
      summary_df <- summary_df %>%
        left_join(stats_df %>% select(all_of(stats_cols)), by = "depth_cm")
    }
  }
  
  out_csv <- file.path(OUTPUT_DIR, "prediction_statistics.csv")
  write_csv(summary_df, out_csv)
  
  cat("\n==================================================\n")
  cat("   PREDICTION SUMMARY\n")
  cat("==================================================\n")
  print(summary_df %>% mutate(across(where(is.numeric), ~ round(., 3))))
  
  cat(sprintf("\nExported: %s\n", out_csv))
}

# ============================================================================
# 7. FINAL SUMMARY
# ============================================================================

cat("\n==================================================\n")
cat("   VISUALIZATION COMPLETE\n")
cat("==================================================\n")

cat(sprintf("\nOutput directory: %s\n", normalizePath(OUTPUT_DIR, mustWork = FALSE)))

cat("\nGenerated files:\n")
for (d in available_depths) {
  cat(sprintf("  - comparison_plot_depth_%s.png\n", d))
}
if (length(transfer_rasters) > 0) {
  cat("  - all_depths_transfer_final.png\n")
}
cat("  - prediction_statistics.csv\n")

cat("\nInterpretation guide:\n")
cat("  Global Prior    - Predictions from Pacific-wide weighted model\n")
cat("  Local-Only      - Naive model using only local field data\n")
cat("  Transfer Final  - Bias-corrected transfer learning prediction\n")
cat("  Difference      - Spatial pattern of transfer learning adjustment\n")

cat("\nNext steps:\n")
cat("  1. Review difference maps for systematic spatial patterns\n
")
cat("  2. Compare transfer final uncertainty with local-only\n")
cat("  3. Export final predictions for VM0033 carbon accounting\n")

cat("\n==================================================\n")