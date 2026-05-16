# ============================================================================
# MODULE 06C: BAYESIAN POSTERIOR ESTIMATION (Part 4 - Optional)
# ============================================================================
# PURPOSE: Combine Bayesian priors with field data to generate posterior estimates
#          of carbon stocks (kg/m²) at VM0033 standard depths
#
# PREREQUISITES:
#   - Module 00C: Bayesian prior setup (carbon stocks kg/m²)
#   - Modules 01-03: Data collection and harmonization
#   - Modules 04-05: RF/Kriging predictions (likelihood - carbon stocks kg/m²)
#
# THEORY: Bayesian Update for Carbon Stocks
#   Prior × Likelihood → Posterior (all in kg/m² units)
#
#   Precision-weighted average:
#   τ_prior = 1 / σ²_prior
#   τ_field = 1 / σ²_field
#
#   μ_posterior = (τ_prior × μ_prior + τ_field × μ_field) / (τ_prior + τ_field)
#   σ²_posterior = 1 / (τ_prior + τ_field)
#
# IMPORTANT: All data are carbon stocks (kg/m²), NOT SOC concentrations (g/kg)
#
# INPUTS:
#   - data_prior/carbon_stock_prior_mean_*.tif (from Module 00C - carbon stocks kg/m²)
#   - data_prior/carbon_stock_prior_se_*.tif (from Module 00C - carbon stocks kg/m²)
#   - outputs/predictions/rf/carbon_stock_rf_*.tif (from Module 05 - carbon stocks kg/m²)
#   - outputs/predictions/rf/se_combined_*.tif (from Module 05 - uncertainty kg/m²)
#   OR
#   - outputs/predictions/kriging/carbon_stock_*_*.tif (from Module 04 - carbon stocks kg/m²)
#   - outputs/predictions/kriging/se_combined_*_*.tif (from Module 04 - uncertainty kg/m²)
#   - data_processed/cores_harmonized_bluecarbon.rds (sample locations)
#
# OUTPUTS:
#   - outputs/predictions/posterior/carbon_stock_posterior_mean_*.tif (kg/m²)
#   - outputs/predictions/posterior/carbon_stock_posterior_se_*.tif (kg/m²)
#   - outputs/predictions/posterior/carbon_stock_posterior_conservative_*.tif (kg/m²)
#   - diagnostics/bayesian/information_gain_*.tif
#   - diagnostics/bayesian/uncertainty_reduction.csv
#   - diagnostics/bayesian/prior_likelihood_posterior_comparison.png
#
# NOTE: All units are carbon stocks (kg/m²) for consistency with the updated workflow.
#       Prior and likelihood must both be in kg/m² for proper Bayesian updating.
#
# ============================================================================

# Clear workspace
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
rm(list = ls())

# Load required libraries
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(gridExtra)
})

# Load configuration
source("Analysis_Workflow/blue_carbon_config.R")

# ============================================================================
# SETUP LOGGING
# ============================================================================

log_message <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  cat(sprintf("%s %s: %s\n", timestamp, level, msg))
}

log_message("=== MODULE 06C: BAYESIAN POSTERIOR ESTIMATION ===")
log_message(sprintf("Project: %s", PROJECT_NAME))

# Check if Bayesian workflow is enabled
if (!USE_BAYESIAN) {
  stop("Bayesian workflow is disabled (USE_BAYESIAN = FALSE).\n",
       "Set USE_BAYESIAN <- TRUE in blue_carbon_config.R to enable Part 4.\n",
       "Or use standard Module 06 for non-Bayesian analysis.")
}

log_message("Bayesian posterior estimation enabled ✓")

# Create output directories
dir.create("outputs/predictions/posterior", recursive = TRUE, showWarnings = FALSE)
dir.create("diagnostics/bayesian", recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD FIELD DATA (for reference and logging)
# ============================================================================

log_message("\nLoading field sample locations...")

cores_file <- "data_processed/cores_harmonized_bluecarbon.rds"

if (!file.exists(cores_file)) {
  stop(sprintf("Harmonized cores file not found: %s\n", cores_file),
       "Please run Modules 01-03 first to harmonize field data.")
}

cores <- readRDS(cores_file)

# Get unique sample locations
sample_locations <- cores %>%
  select(core_id, longitude, latitude) %>%
  distinct()

log_message(sprintf("Loaded %d sample locations", nrow(sample_locations)))

# Convert to spatial
samples_sf <- st_as_sf(sample_locations,
                      coords = c("longitude", "latitude"),
                      crs = INPUT_CRS)

# Reproject to processing CRS
samples_sf <- st_transform(samples_sf, crs = PROCESSING_CRS)

# ============================================================================
# CHECK FOR LIKELIHOOD MAPS
# ============================================================================

log_message("\nChecking for likelihood maps (RF or Kriging)...")
log_message("Note: Looking for carbon stock predictions (kg/m²) from updated workflow")

# Check RF (updated to carbon stock files)
rf_dir <- "outputs/predictions/rf"
rf_files <- list.files(rf_dir, pattern = "^carbon_stock_rf_[0-9]+cm\\.tif$", full.names = TRUE)

# Check Kriging (updated to carbon stock files)
kriging_dir <- "outputs/predictions/kriging"
kriging_files <- list.files(kriging_dir, pattern = "^carbon_stock_.*_[0-9]+cm\\.tif$", full.names = TRUE)

use_rf <- length(rf_files) > 0
use_kriging <- length(kriging_files) > 0

if (!use_rf && !use_kriging) {
  stop("No likelihood maps found.\n",
       "Expected files:\n",
       "  RF: outputs/predictions/rf/carbon_stock_rf_*cm.tif\n",
       "  Kriging: outputs/predictions/kriging/carbon_stock_*_*cm.tif\n",
       "Please run Module 04 (Kriging) or Module 05 (RF) first to generate predictions.")
}

# Prefer RF if both available
if (use_rf && use_kriging) {
  log_message("Both RF and Kriging available - using RF for posterior", "INFO")
  likelihood_method <- "rf"
  likelihood_dir <- rf_dir
} else if (use_rf) {
  log_message(sprintf("Using RF for likelihood (%d carbon stock files found)", length(rf_files)))
  likelihood_method <- "rf"
  likelihood_dir <- rf_dir
} else {
  log_message(sprintf("Using Kriging for likelihood (%d carbon stock files found)", length(kriging_files)))
  likelihood_method <- "kriging"
  likelihood_dir <- kriging_dir
}

# ============================================================================
# LOAD PRIORS
# ============================================================================

log_message("\nLoading Bayesian priors (carbon stocks in kg/m²)...")

if (!dir.exists(BAYESIAN_PRIOR_DIR)) {
  stop(sprintf("Prior directory not found: %s\n", BAYESIAN_PRIOR_DIR),
       "Please run Module 00C first.")
}

# Get VM0033 depths
vm0033_depths <- c(7.5, 22.5, 40, 75)

priors <- list()

for (depth in vm0033_depths) {
  # Updated to carbon stock prior files (kg/m²)
  prior_mean_file <- file.path(BAYESIAN_PRIOR_DIR, sprintf("carbon_stock_prior_mean_%.1fcm.tif", depth))
  prior_se_file <- file.path(BAYESIAN_PRIOR_DIR, sprintf("carbon_stock_prior_se_%.1fcm.tif", depth))

  if (file.exists(prior_mean_file) && file.exists(prior_se_file)) {
    priors[[as.character(depth)]] <- list(
      mean = rast(prior_mean_file),
      se = rast(prior_se_file)
    )
    log_message(sprintf("  Loaded carbon stock prior for %.1f cm", depth))
  } else {
    log_message(sprintf("  Carbon stock prior not found for %.1f cm - skipping", depth), "WARNING")
    log_message(sprintf("    Expected: %s", basename(prior_mean_file)), "WARNING")
  }
}

if (length(priors) == 0) {
  stop("No carbon stock prior maps loaded.\n",
       "Please run Module 00C first to process carbon stock priors.\n",
       "Expected files: data_prior/carbon_stock_prior_mean_*.tif")
}

# ============================================================================
# SETUP COMPLETE - READY FOR BAYESIAN UPDATE
# ============================================================================

log_message("\nPreparing for Bayesian posterior estimation...")
log_message("  Prior source: GEE SoilGrids + Sothe et al. (if available)")
log_message("  Likelihood source: RF/Kriging spatial predictions")
log_message("  Method: Precision-weighted Bayesian update")
log_message(sprintf("  Sample locations: %d", nrow(sample_locations)))

# Note: We use RF/Kriging SE directly as likelihood uncertainty
# This already reflects model confidence based on cross-validation
# No additional sample density weighting needed - keeps it simple and robust

# ============================================================================
# BAYESIAN UPDATE FOR EACH DEPTH
# ============================================================================

log_message("\n=== BAYESIAN POSTERIOR ESTIMATION ===")

posterior_results <- list()
comparison_data <- data.frame()

for (depth_str in names(priors)) {

  depth <- as.numeric(depth_str)
  log_message(sprintf("\nProcessing depth: %.1f cm", depth))

  # === Load prior ===
  prior_mean <- priors[[depth_str]]$mean
  prior_se <- priors[[depth_str]]$se

  # Validate prior data
  log_message("  Validating prior data...")
  prior_mean_stats <- global(prior_mean, c("min", "max", "mean"), na.rm = TRUE)
  prior_se_stats <- global(prior_se, c("min", "max", "mean"), na.rm = TRUE)

  n_valid_prior <- sum(!is.na(values(prior_mean, mat = FALSE)))
  n_total_prior <- length(values(prior_mean, mat = FALSE))

  log_message(sprintf("    Prior mean: min=%.2f, max=%.2f, mean=%.2f (%d/%d valid cells)",
                     prior_mean_stats[1,1], prior_mean_stats[1,2], prior_mean_stats[1,3],
                     n_valid_prior, n_total_prior))
  log_message(sprintf("    Prior SE: min=%.2f, max=%.2f, mean=%.2f",
                     prior_se_stats[1,1], prior_se_stats[1,2], prior_se_stats[1,3]))

  if (n_valid_prior == 0) {
    log_message("  ERROR: Prior raster has no valid values - skipping depth", "ERROR")
    next
  }

  if (is.na(prior_se_stats[1,3]) || prior_se_stats[1,3] <= 0) {
    log_message("  ERROR: Prior SE has no valid values or all zeros - skipping depth", "ERROR")
    next
  }

  # === Load likelihood (field data - carbon stocks in kg/m²) ===
  # Note: Files use rounded depths (7.5→8, 22.5→22) in filenames
  depth_rounded <- round(depth)

  if (likelihood_method == "rf") {
    # RF carbon stock files: carbon_stock_rf_*cm.tif
    likelihood_mean_file <- file.path(likelihood_dir, sprintf("carbon_stock_rf_%dcm.tif", depth_rounded))
    # RF uncertainty files: se_combined_*cm.tif
    likelihood_se_file <- file.path(likelihood_dir, sprintf("se_combined_%dcm.tif", depth_rounded))
  } else {
    # Kriging carbon stock files: carbon_stock_*_*cm.tif (may have stratum name)
    # Need to find files matching the pattern for this depth
    kriging_pattern <- sprintf("carbon_stock_.*_%dcm\\.tif$", depth_rounded)
    kriging_matches <- list.files(likelihood_dir, pattern = kriging_pattern, full.names = TRUE)

    if (length(kriging_matches) > 0) {
      likelihood_mean_file <- kriging_matches[1]  # Take first match
    } else {
      likelihood_mean_file <- ""  # Will fail check below
    }

    # Kriging uncertainty files: se_combined_*_*cm.tif
    se_pattern <- sprintf("se_combined_.*_%dcm\\.tif$", depth_rounded)
    se_matches <- list.files(likelihood_dir, pattern = se_pattern, full.names = TRUE)

    if (length(se_matches) > 0) {
      likelihood_se_file <- se_matches[1]
    } else {
      likelihood_se_file <- ""
    }
  }

  if (!file.exists(likelihood_mean_file)) {
    log_message(sprintf("  ERROR: Likelihood mean not found: %s", basename(likelihood_mean_file)), "ERROR")
    next
  }

  likelihood_mean <- rast(likelihood_mean_file)

  if (file.exists(likelihood_se_file)) {
    likelihood_se <- rast(likelihood_se_file)
  } else {
    log_message("  WARNING: SE not found - using 15% of mean", "WARNING")
    likelihood_se <- likelihood_mean * 0.15
  }

  # Validate likelihood data
  log_message("  Validating likelihood data...")
  likelihood_mean_stats <- global(likelihood_mean, c("min", "max", "mean"), na.rm = TRUE)
  likelihood_se_stats <- global(likelihood_se, c("min", "max", "mean"), na.rm = TRUE)

  n_valid_likelihood <- sum(!is.na(values(likelihood_mean, mat = FALSE)))
  n_total_likelihood <- length(values(likelihood_mean, mat = FALSE))

  log_message(sprintf("    Likelihood mean: min=%.2f, max=%.2f, mean=%.2f (%d/%d valid cells)",
                     likelihood_mean_stats[1,1], likelihood_mean_stats[1,2], likelihood_mean_stats[1,3],
                     n_valid_likelihood, n_total_likelihood))
  log_message(sprintf("    Likelihood SE: min=%.2f, max=%.2f, mean=%.2f",
                     likelihood_se_stats[1,1], likelihood_se_stats[1,2], likelihood_se_stats[1,3]))

  if (n_valid_likelihood == 0) {
    log_message("  ERROR: Likelihood raster has no valid values - skipping depth", "ERROR")
    next
  }

  if (is.na(likelihood_se_stats[1,3]) || likelihood_se_stats[1,3] <= 0) {
    log_message("  ERROR: Likelihood SE has no valid values or all zeros - skipping depth", "ERROR")
    next
  }

  # Ensure spatial alignment
  if (!compareGeom(prior_mean, likelihood_mean, stopOnError = FALSE)) {
    log_message("  Resampling likelihood to match prior grid")
    likelihood_mean <- resample(likelihood_mean, prior_mean, method = "bilinear")
    likelihood_se <- resample(likelihood_se, prior_se, method = "bilinear")
  }

  # Check for spatial overlap after resampling
  log_message("  Checking spatial overlap...")
  # Create a mask of overlapping non-NA cells
  overlap_mask <- !is.na(prior_mean) & !is.na(likelihood_mean) &
                  !is.na(prior_se) & !is.na(likelihood_se)
  n_overlap <- sum(values(overlap_mask, mat = FALSE), na.rm = TRUE)

  log_message(sprintf("    Overlapping valid cells: %d", n_overlap))

  if (n_overlap == 0) {
    log_message("  ERROR: No spatial overlap between prior and likelihood - check study areas", "ERROR")
    log_message("    Prior extent:", "ERROR")
    log_message(sprintf("      %s", as.character(ext(prior_mean))), "ERROR")
    log_message("    Likelihood extent:", "ERROR")
    log_message(sprintf("      %s", as.character(ext(likelihood_mean))), "ERROR")
    next
  }

  # === Calculate precisions (inverse variance) ===
  # Add safeguards to prevent Inf/NaN
  log_message("  Computing Bayesian precisions...")

  prior_var <- prior_se^2
  prior_var <- clamp(prior_var, lower = 0.001, upper = 1000)  # Prevent extreme values
  tau_prior <- 1 / prior_var

  field_var <- likelihood_se^2
  field_var <- clamp(field_var, lower = 0.001, upper = 1000)
  tau_field <- 1 / field_var

  # Check for invalid precisions
  if (any(is.na(values(tau_prior, mat = FALSE)), na.rm = FALSE)) {
    log_message("  WARNING: Prior precision contains NA values", "WARNING")
  }
  if (any(is.na(values(tau_field, mat = FALSE)), na.rm = FALSE)) {
    log_message("  WARNING: Field precision contains NA values", "WARNING")
  }
  if (any(is.infinite(values(tau_prior, mat = FALSE)), na.rm = TRUE)) {
    log_message("  WARNING: Prior precision contains Inf values", "WARNING")
  }
  if (any(is.infinite(values(tau_field, mat = FALSE)), na.rm = TRUE)) {
    log_message("  WARNING: Field precision contains Inf values", "WARNING")
  }

  # === Diagnostic output ===
  log_message("  Diagnostics:")
  log_message(sprintf("    Prior SE range: %.2f - %.2f (mean: %.2f)",
                     global(prior_se, "min", na.rm = TRUE)[1,1],
                     global(prior_se, "max", na.rm = TRUE)[1,1],
                     global(prior_se, "mean", na.rm = TRUE)[1,1]))
  log_message(sprintf("    Likelihood SE range: %.2f - %.2f (mean: %.2f)",
                     global(likelihood_se, "min", na.rm = TRUE)[1,1],
                     global(likelihood_se, "max", na.rm = TRUE)[1,1],
                     global(likelihood_se, "mean", na.rm = TRUE)[1,1]))

  # === Bayesian posterior ===
  log_message("  Computing posterior...")

  # Posterior mean (precision-weighted average)
  posterior_mean <- (tau_prior * prior_mean + tau_field * likelihood_mean) /
                    (tau_prior + tau_field)

  # Posterior variance
  tau_total <- tau_prior + tau_field
  # Safeguard against division by zero
  tau_total <- clamp(tau_total, lower = 0.001, upper = Inf)
  posterior_var <- 1 / tau_total

  # Ensure variance is reasonable
  posterior_var <- clamp(posterior_var, lower = 0.001, upper = 1000)
  posterior_se <- sqrt(posterior_var)

  # Check for NaN/Inf in posterior
  if (any(is.na(values(posterior_se, mat = FALSE)), na.rm = FALSE)) {
    log_message("  ERROR: Posterior SE contains NA values - check input data", "ERROR")
  }
  if (any(is.infinite(values(posterior_se, mat = FALSE)), na.rm = TRUE)) {
    log_message("  ERROR: Posterior SE contains Inf values - check input data", "ERROR")
  }

  # Conservative estimate (95% CI lower bound for VM0033)
  posterior_conservative <- posterior_mean - qnorm((1 + CONFIDENCE_LEVEL) / 2) * posterior_se

  # Ensure non-negative carbon stocks
  posterior_conservative <- clamp(posterior_conservative, lower = 0, upper = Inf)

  # === Information gain ===
  # How much did field data reduce uncertainty?
  information_gain <- (1 / posterior_se^2) - (1 / prior_se^2)
  uncertainty_reduction_pct <- (1 - posterior_se / prior_se) * 100

  # === Save outputs ===
  log_message("  Saving posterior rasters...")

  # Format depth for filename (7.5 → 7_5) to match prior file naming
  depth_str_file <- gsub("\\.", "_", sprintf("%.1f", depth))

  # Output files: carbon stock posterior estimates (kg/m²)
  out_mean <- file.path("outputs/predictions/posterior",
                        sprintf("carbon_stock_posterior_mean_%scm.tif", depth_str_file))
  out_se <- file.path("outputs/predictions/posterior",
                      sprintf("carbon_stock_posterior_se_%scm.tif", depth_str_file))
  out_conservative <- file.path("outputs/predictions/posterior",
                               sprintf("carbon_stock_posterior_conservative_%scm.tif", depth_str_file))
  out_info_gain <- file.path("diagnostics/bayesian",
                            sprintf("information_gain_%scm.tif", depth_str_file))

  writeRaster(posterior_mean, out_mean, overwrite = TRUE)
  writeRaster(posterior_se, out_se, overwrite = TRUE)
  writeRaster(posterior_conservative, out_conservative, overwrite = TRUE)
  writeRaster(information_gain, out_info_gain, overwrite = TRUE)

  log_message(sprintf("  Saved: %s", basename(out_mean)))

  # === Extract statistics ===
  prior_vals <- values(prior_mean, mat = FALSE)
  prior_vals <- prior_vals[!is.na(prior_vals)]

  field_vals <- values(likelihood_mean, mat = FALSE)
  field_vals <- field_vals[!is.na(field_vals)]

  post_vals <- values(posterior_mean, mat = FALSE)
  post_vals <- post_vals[!is.na(post_vals)]

  prior_se_vals <- values(prior_se, mat = FALSE)
  prior_se_vals <- prior_se_vals[!is.na(prior_se_vals)]

  post_se_vals <- values(posterior_se, mat = FALSE)
  post_se_vals <- post_se_vals[!is.na(post_se_vals)]

  reduction_vals <- values(uncertainty_reduction_pct, mat = FALSE)
  reduction_vals <- reduction_vals[!is.na(reduction_vals)]

  # Store results
  posterior_results[[depth_str]] <- list(
    depth = depth,
    prior_mean = mean(prior_vals, na.rm = TRUE),
    prior_se = mean(prior_se_vals, na.rm = TRUE),
    field_mean = mean(field_vals, na.rm = TRUE),
    posterior_mean = mean(post_vals, na.rm = TRUE),
    posterior_se = mean(post_se_vals, na.rm = TRUE),
    uncertainty_reduction_pct = mean(reduction_vals, na.rm = TRUE)
  )

  # Comparison data for visualization
  comp_df <- data.frame(
    depth = depth,
    estimate = c(rep("Prior", length(prior_vals)),
                rep("Field", length(field_vals)),
                rep("Posterior", length(post_vals))),
    value = c(prior_vals, field_vals, post_vals)
  )

  comparison_data <- rbind(comparison_data, comp_df)

  log_message(sprintf("  Prior: %.2f ± %.2f kg/m²", mean(prior_vals), mean(prior_se_vals)))
  log_message(sprintf("  Field: %.2f kg/m²", mean(field_vals)))
  log_message(sprintf("  Posterior: %.2f ± %.2f kg/m²", mean(post_vals), mean(post_se_vals)))
  log_message(sprintf("  Uncertainty reduction: %.1f%%", mean(reduction_vals)))
}

# ============================================================================
# SAVE SUMMARY STATISTICS
# ============================================================================

log_message("\nSaving summary statistics...")

summary_df <- do.call(rbind, lapply(posterior_results, as.data.frame))

write_csv(summary_df, "diagnostics/bayesian/uncertainty_reduction.csv")
log_message("Saved: diagnostics/bayesian/uncertainty_reduction.csv")

# ============================================================================
# CREATE COMPARISON VISUALIZATIONS
# ============================================================================

log_message("\nCreating comparison visualizations...")

# Density plot comparing prior, field, and posterior
p1 <- ggplot(comparison_data, aes(x = value, fill = estimate)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~depth, scales = "free", ncol = 2,
            labeller = labeller(depth = function(x) sprintf("%.1f cm", as.numeric(x)))) +
  scale_fill_manual(values = c("Prior" = "#3498db", "Field" = "#e67e22", "Posterior" = "#2ecc71")) +
  labs(
    title = "Bayesian Update: Prior × Likelihood → Posterior",
    subtitle = sprintf("%s - %s (Carbon Stocks)", PROJECT_NAME, likelihood_method),
    x = "Carbon Stock (kg/m²)",
    y = "Density",
    fill = "Estimate"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Uncertainty reduction plot
p2 <- ggplot(summary_df, aes(x = factor(depth), y = uncertainty_reduction_pct)) +
  geom_col(fill = "#2ecc71") +
  geom_text(aes(label = sprintf("%.1f%%", uncertainty_reduction_pct)),
           vjust = -0.5, size = 4) +
  labs(
    title = "Uncertainty Reduction by Depth",
    x = "Depth (cm)",
    y = "Uncertainty Reduction (%)"
  ) +
  theme_minimal() +
  ylim(0, max(summary_df$uncertainty_reduction_pct) * 1.2)

# Combine plots
p_combined <- gridExtra::grid.arrange(p1, p2, nrow = 2, heights = c(2, 1))

ggsave("diagnostics/bayesian/prior_likelihood_posterior_comparison.png",
      p_combined, width = 12, height = 10, dpi = 300)

log_message("Saved: diagnostics/bayesian/prior_likelihood_posterior_comparison.png")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("BAYESIAN POSTERIOR ESTIMATION COMPLETE\n")
cat("========================================\n\n")

cat(sprintf("Likelihood method: %s\n", toupper(likelihood_method)))
cat(sprintf("Depths processed: %s\n", paste(summary_df$depth, collapse = ", ")))
cat(sprintf("Sample locations: %d\n\n", nrow(sample_locations)))

cat("Uncertainty Reduction:\n")
for (i in 1:nrow(summary_df)) {
  cat(sprintf("  %.1f cm: %.1f%% (%.2f → %.2f kg/m² SE)\n",
              summary_df$depth[i],
              summary_df$uncertainty_reduction_pct[i],
              summary_df$prior_se[i],
              summary_df$posterior_se[i]))
}

overall_reduction <- mean(summary_df$uncertainty_reduction_pct, na.rm = TRUE)
cat(sprintf("\nOverall mean reduction: %.1f%%\n\n", overall_reduction))

# Handle NaN/NA values in comparison
if (is.na(overall_reduction) || is.nan(overall_reduction)) {
  cat("⚠ WARNING: Could not calculate uncertainty reduction\n")
  cat("  Check input data quality (priors and likelihood SE)\n")
  cat("  Possible issues:\n")
  cat("    - Prior SE files may have incorrect units or values\n")
  cat("    - Likelihood SE files may be missing or invalid\n\n")
} else if (overall_reduction >= MIN_INFORMATION_GAIN_PCT) {
  cat(sprintf("✓ Information gain exceeds threshold (>%.0f%%)\n", MIN_INFORMATION_GAIN_PCT))
  cat("  Prior was informative - Bayesian update successful\n\n")
} else {
  cat(sprintf("⚠ Information gain below threshold (<%.0f%%)\n", MIN_INFORMATION_GAIN_PCT))
  cat("  Prior had limited information - field data dominated\n\n")
}

cat("Output files:\n")
cat(sprintf("  - %d posterior mean rasters\n", nrow(summary_df)))
cat(sprintf("  - %d posterior SE rasters\n", nrow(summary_df)))
cat(sprintf("  - %d posterior conservative rasters\n", nrow(summary_df)))
cat(sprintf("  - %d information gain rasters\n", nrow(summary_df)))
cat("  - 1 uncertainty reduction CSV\n")
cat("  - 1 comparison visualization\n\n")

cat("NEXT STEPS:\n")
cat("1. Review uncertainty_reduction.csv and comparison plots\n")
cat("2. Use posterior rasters for carbon stock calculation:\n")
cat("   - Replace RF/Kriging inputs in Module 06 with posterior/\n")
cat("   - Or continue with standard workflow and compare results\n")
cat("3. For temporal analysis: Run Module 08A-10 with posterior estimates\n")
cat("\n")

# ============================================================================
# PRESENTATION COPY: BAYESIAN ANALYSIS OUTPUTS
# ============================================================================

log_message("Copying key Bayesian outputs to outputs/Bayesian_analysis...")

dir.create("outputs/Bayesian_analysis", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Bayesian_analysis/Posterior_Distributions", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Bayesian_analysis/Prior_vs_Posterior_Maps", recursive = TRUE, showWarnings = FALSE)

get_depth_interval_label <- function(depth_cm) {
  depth_lookup <- data.frame(
    depth_midpoint = c(7.5, 22.5, 40, 75),
    depth_top = c(0, 15, 30, 50),
    depth_bottom = c(15, 30, 50, 100)
  )

  idx <- which.min(abs(depth_lookup$depth_midpoint - depth_cm))
  sprintf("%d_to_%dcm", depth_lookup$depth_top[idx], depth_lookup$depth_bottom[idx])
}

build_bayesian_map_name <- function(src_file) {
  src_name <- tools::file_path_sans_ext(basename(src_file))
  src_ext <- tools::file_ext(src_file)

  if (grepl("^carbon_stock_prior_mean_[0-9]+\\.[0-9]cm$", src_name)) {
    depth_cm <- as.numeric(gsub("^carbon_stock_prior_mean_([0-9]+\\.[0-9])cm$", "\\1", src_name))
    return(sprintf("Prior_Mean_Carbon_Stock_Map_%s.%s", get_depth_interval_label(depth_cm), src_ext))
  }
  if (grepl("^carbon_stock_rf_[0-9]+cm$", src_name)) {
    depth_cm <- as.numeric(gsub("^carbon_stock_rf_([0-9]+)cm$", "\\1", src_name))
    return(sprintf("Likelihood_RF_Carbon_Stock_Map_%s.%s", get_depth_interval_label(depth_cm), src_ext))
  }
  if (grepl("^carbon_stock_posterior_mean_[0-9]+_[0-9]cm$", src_name)) {
    depth_cm <- as.numeric(gsub("^carbon_stock_posterior_mean_([0-9]+_[0-9])cm$", "\\1", src_name))
    depth_cm <- as.numeric(gsub("_", ".", depth_cm))
    return(sprintf("Posterior_Mean_Carbon_Stock_Map_%s.%s", get_depth_interval_label(depth_cm), src_ext))
  }
  if (grepl("^carbon_stock_posterior_se_[0-9]+_[0-9]cm$", src_name)) {
    depth_cm <- as.numeric(gsub("^carbon_stock_posterior_se_([0-9]+_[0-9])cm$", "\\1", src_name))
    depth_cm <- as.numeric(gsub("_", ".", depth_cm))
    return(sprintf("Posterior_Standard_Error_Map_%s.%s", get_depth_interval_label(depth_cm), src_ext))
  }
  if (grepl("^carbon_stock_posterior_conservative_[0-9]+_[0-9]cm$", src_name)) {
    depth_cm <- as.numeric(gsub("^carbon_stock_posterior_conservative_([0-9]+_[0-9])cm$", "\\1", src_name))
    depth_cm <- as.numeric(gsub("_", ".", depth_cm))
    return(sprintf("Posterior_Conservative_Carbon_Stock_Map_%s.%s", get_depth_interval_label(depth_cm), src_ext))
  }
  if (grepl("^information_gain_[0-9]+_[0-9]cm$", src_name)) {
    depth_cm <- as.numeric(gsub("^information_gain_([0-9]+_[0-9])cm$", "\\1", src_name))
    depth_cm <- as.numeric(gsub("_", ".", depth_cm))
    return(sprintf("Bayesian_Information_Gain_Map_%s.%s", get_depth_interval_label(depth_cm), src_ext))
  }

  sprintf("Bayesian_Output_%s.%s", src_name, src_ext)
}

hero_candidates <- c(
  "outputs/predictions/stocks/stock_total_0-100cm_bayesian.tif",
  "outputs/predictions/stocks/carbon_stock_total_0-100cm_bayesian.tif",
  "outputs/predictions/posterior/carbon_stock_posterior_mean_75_0cm.tif"
)
hero_src <- hero_candidates[file.exists(hero_candidates)][1]
if (!is.na(hero_src)) {
  hero_dst <- "outputs/Bayesian_analysis/Bayesian_Final_Updated_Carbon_Stock_Map_0_to_100cm.tif"
  file.copy(hero_src, hero_dst, overwrite = TRUE)
  log_message(sprintf("Copied Bayesian hero output: %s", basename(hero_dst)))
} else {
  log_message("No Bayesian hero map found for root Bayesian_analysis copy step", "WARNING")
}

distribution_sources <- c(
  "diagnostics/bayesian/uncertainty_reduction.csv",
  "diagnostics/bayesian/prior_likelihood_posterior_comparison.png"
)
for (src in distribution_sources[file.exists(distribution_sources)]) {
  dst_name <- if (grepl("uncertainty_reduction\\.csv$", src)) {
    "Bayesian_Uncertainty_Reduction_Summary.csv"
  } else {
    "Bayesian_Prior_Likelihood_Posterior_Distribution_Plot.png"
  }
  dst <- file.path("outputs/Bayesian_analysis/Posterior_Distributions", dst_name)
  file.copy(src, dst, overwrite = TRUE)
  log_message(sprintf("Copied Bayesian distribution output: %s", basename(dst)))
}

bayes_map_sources <- c(
  list.files("data_prior", pattern = "^carbon_stock_prior_mean_[0-9]+\\.[0-9]cm\\.tif$", full.names = TRUE),
  list.files("outputs/predictions/rf", pattern = "^carbon_stock_rf_[0-9]+cm\\.tif$", full.names = TRUE),
  list.files("outputs/predictions/posterior", pattern = "^carbon_stock_posterior_(mean|se|conservative)_[0-9]+_[0-9]cm\\.tif$", full.names = TRUE),
  list.files("diagnostics/bayesian", pattern = "^information_gain_[0-9]+_[0-9]cm\\.tif$", full.names = TRUE)
)

if (length(bayes_map_sources) == 0) {
  log_message("No Bayesian map layers found for Prior_vs_Posterior_Maps copy step", "WARNING")
} else {
  for (src in bayes_map_sources) {
    dst_name <- build_bayesian_map_name(src)
    dst <- file.path("outputs/Bayesian_analysis/Prior_vs_Posterior_Maps", dst_name)
    file.copy(src, dst, overwrite = TRUE)
    log_message(sprintf("Copied Bayesian map layer: %s", basename(dst)))
  }
}

log_message("=== MODULE 06C COMPLETE ===")
