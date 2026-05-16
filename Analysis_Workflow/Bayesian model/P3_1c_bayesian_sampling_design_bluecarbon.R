# ============================================================================
# MODULE 01C: BAYESIAN SAMPLING DESIGN - NEYMAN ALLOCATION (Part 4 - Optional)
# ============================================================================
# PURPOSE: Design optimal sampling strategy using prior uncertainty
#
# PREREQUISITES:
#   - Run Module 00C first to process Bayesian priors
#
# INPUTS:
#   - data_prior/soc_prior_mean_*.tif (from Module 00C)
#   - data_prior/soc_prior_se_*.tif (from Module 00C)
#   - data_prior/uncertainty_strata.tif (from Module 00C, or will create)
#   - blue_carbon_config.R (configuration)
#
# OUTPUTS:
#   - sampling_locations_neyman.csv (optimized sample locations)
#   - sampling_allocation_neyman.csv (samples per stratum)
#   - sampling_map_neyman.png (visualization)
#   - data_processed/neyman_strata.tif (uncertainty strata raster)
#
# THEORY: Neyman Allocation
#   Allocate samples proportional to: n_h ∝ N_h × σ_h
#   Where N_h = area of stratum h, σ_h = standard deviation in stratum h
#   This minimizes total variance for fixed sample size
#
# ============================================================================

# Clear workspace
rm(list = ls())

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
# Load required libraries
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(spatstat)  # For spatial point pattern analysis
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

log_message("=== MODULE 01C: BAYESIAN SAMPLING DESIGN (NEYMAN ALLOCATION) ===")
log_message(sprintf("Project: %s", PROJECT_NAME))

# Check if Bayesian workflow is enabled
if (!USE_BAYESIAN) {
  stop("Bayesian workflow is disabled (USE_BAYESIAN = FALSE).\n",
       "Set USE_BAYESIAN <- TRUE in blue_carbon_config.R to enable Part 4.\n",
       "Or use standard Module 01 for non-Bayesian sampling.")
}

if (!USE_NEYMAN_SAMPLING) {
  stop("Neyman sampling is disabled (USE_NEYMAN_SAMPLING = FALSE).\n",
       "Set USE_NEYMAN_SAMPLING <- TRUE or use standard Module 01.")
}

log_message("Neyman allocation enabled ✓")

# ============================================================================
# LOAD BAYESIAN PRIORS
# ============================================================================

log_message("\nLoading Bayesian priors...")

if (!dir.exists(BAYESIAN_PRIOR_DIR)) {
  stop(sprintf("Prior directory not found: %s\n", BAYESIAN_PRIOR_DIR),
       "Please run Module 00C first to process Bayesian priors.")
}

# Check for prior files
prior_mean_files <- list.files(BAYESIAN_PRIOR_DIR,
                               pattern = "soc_prior_mean.*\\.tif$",
                               full.names = TRUE)

prior_se_files <- list.files(BAYESIAN_PRIOR_DIR,
                             pattern = "soc_prior_se.*\\.tif$",
                             full.names = TRUE)

if (length(prior_mean_files) == 0) {
  stop("No prior mean files found. Please run Module 00C first.")
}

log_message(sprintf("Found %d prior depth layers", length(prior_mean_files)))

# Load surface layer (7.5cm) for sampling design
surface_mean <- rast(file.path(BAYESIAN_PRIOR_DIR, "soc_prior_mean_7.5cm.tif"))
surface_se <- rast(file.path(BAYESIAN_PRIOR_DIR, "soc_prior_se_7.5cm.tif"))

log_message("Loaded surface layer (7.5cm) for sampling design")

# Calculate coefficient of variation (CV)
cv_raster <- (surface_se / surface_mean) * 100
names(cv_raster) <- "cv_percent"

log_message("Calculated coefficient of variation (CV)")

# ============================================================================
# CREATE OR LOAD UNCERTAINTY STRATA
# ============================================================================

log_message("\nCreating uncertainty strata for Neyman allocation...")

strata_file <- file.path(BAYESIAN_PRIOR_DIR, "uncertainty_strata.tif")

if (file.exists(strata_file)) {
  log_message("Loading existing uncertainty strata from Module 00C")
  uncertainty_strata <- rast(strata_file)
} else {
  log_message("Creating uncertainty strata from CV thresholds")

  # Create strata based on CV thresholds
  uncertainty_strata <- cv_raster
  uncertainty_strata[cv_raster < UNCERTAINTY_LOW_THRESHOLD] <- 1   # Low uncertainty
  uncertainty_strata[cv_raster >= UNCERTAINTY_LOW_THRESHOLD &
                    cv_raster < UNCERTAINTY_HIGH_THRESHOLD] <- 2  # Medium uncertainty
  uncertainty_strata[cv_raster >= UNCERTAINTY_HIGH_THRESHOLD] <- 3 # High uncertainty

  names(uncertainty_strata) <- "stratum"

  # Save for future use
  writeRaster(uncertainty_strata, strata_file, overwrite = TRUE)
  log_message("Saved uncertainty strata")
}

# ============================================================================
# CALCULATE STRATUM STATISTICS
# ============================================================================

log_message("\nCalculating statistics per uncertainty stratum...")

# Extract values
mean_vals <- values(surface_mean, mat = FALSE)
se_vals <- values(surface_se, mat = FALSE)
cv_vals <- values(cv_raster, mat = FALSE)
strata_vals <- values(uncertainty_strata, mat = FALSE)

# Create data frame
prior_data <- data.frame(
  mean = mean_vals,
  se = se_vals,
  cv = cv_vals,
  stratum = strata_vals
) %>%
  filter(!is.na(stratum))

# Calculate per-stratum statistics
stratum_stats <- prior_data %>%
  group_by(stratum) %>%
  summarise(
    n_pixels = n(),
    area_ha = n() * (res(surface_mean)[1]^2) / 10000,
    mean_soc = mean(mean, na.rm = TRUE),
    sd_soc = sd(mean, na.rm = TRUE),
    mean_cv = mean(cv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    stratum_label = case_when(
      stratum == 1 ~ "Low Uncertainty",
      stratum == 2 ~ "Medium Uncertainty",
      stratum == 3 ~ "High Uncertainty",
      TRUE ~ "Unknown"
    ),
    pct_area = 100 * area_ha / sum(area_ha)
  )

log_message("\nStratum Summary:")
for (i in 1:nrow(stratum_stats)) {
  log_message(sprintf("  %s: %.1f ha (%.1f%%), CV=%.1f%%, SD=%.1f g/kg",
                     stratum_stats$stratum_label[i],
                     stratum_stats$area_ha[i],
                     stratum_stats$pct_area[i],
                     stratum_stats$mean_cv[i],
                     stratum_stats$sd_soc[i]))
}

# ============================================================================
# NEYMAN ALLOCATION
# ============================================================================

log_message("\n=== NEYMAN ALLOCATION ===")

# Total samples needed (from config)
total_samples_base <- VM0033_MIN_CORES * length(VALID_STRATA)
total_samples_target <- ceiling(total_samples_base * NEYMAN_BUFFER_SAMPLES)

log_message(sprintf("Target total samples: %d (base=%d, buffer=%.1fx)",
                   total_samples_target, total_samples_base, NEYMAN_BUFFER_SAMPLES))

# Neyman allocation formula: n_h ∝ N_h × σ_h
stratum_allocation <- stratum_stats %>%
  mutate(
    # Allocation weight: area × standard deviation
    weight = area_ha * sd_soc,

    # Proportional allocation
    samples_neyman = (weight / sum(weight)) * total_samples_target,

    # Round up to ensure minimum samples per stratum
    samples_allocated = pmax(VM0033_MIN_CORES, ceiling(samples_neyman)),

    # Samples per hectare
    density_samples_per_ha = samples_allocated / area_ha
  )

# Adjust if total exceeds target due to rounding
actual_total <- sum(stratum_allocation$samples_allocated)
if (actual_total > total_samples_target) {
  excess <- actual_total - total_samples_target
  log_message(sprintf("Adjusting for %d excess samples from rounding", excess), "INFO")

  # Remove excess from largest strata first
  stratum_allocation <- stratum_allocation %>%
    arrange(desc(samples_allocated)) %>%
    mutate(
      reduction = pmin(samples_allocated - VM0033_MIN_CORES, excess),
      samples_allocated = samples_allocated - reduction
    ) %>%
    select(-reduction)
}

log_message("\nNeyman Allocation Results:")
for (i in 1:nrow(stratum_allocation)) {
  log_message(sprintf("  %s: %d samples (%.1f samples/ha)",
                     stratum_allocation$stratum_label[i],
                     stratum_allocation$samples_allocated[i],
                     stratum_allocation$density_samples_per_ha[i]))
}

log_message(sprintf("Total allocated: %d samples", sum(stratum_allocation$samples_allocated)))

# Save allocation table
write_csv(stratum_allocation, "sampling_allocation_neyman.csv")
log_message("\nSaved: sampling_allocation_neyman.csv")

# ============================================================================
# GENERATE SPATIALLY-BALANCED SAMPLE POINTS
# ============================================================================

log_message("\n=== GENERATING SAMPLE LOCATIONS ===")

all_sample_points <- list()

for (i in 1:nrow(stratum_allocation)) {

  stratum_id <- stratum_allocation$stratum[i]
  n_samples <- stratum_allocation$samples_allocated[i]
  stratum_name <- stratum_allocation$stratum_label[i]

  log_message(sprintf("\nGenerating %d points for %s...", n_samples, stratum_name))

  # Create mask for this stratum
  stratum_mask <- uncertainty_strata == stratum_id

  # Convert to polygon
  stratum_poly <- as.polygons(stratum_mask, values = TRUE, na.rm = TRUE)
  stratum_poly <- stratum_poly[stratum_poly[[1]] == 1, ]

  if (nrow(stratum_poly) == 0) {
    log_message(sprintf("  No valid area for %s", stratum_name), "WARNING")
    next
  }

  # Convert to sf
  stratum_sf <- st_as_sf(stratum_poly)
  stratum_sf <- st_union(stratum_sf)  # Merge polygons

  # Generate spatially-balanced points using systematic sampling
  # This ensures good spatial coverage
  sample_points <- st_sample(stratum_sf, size = n_samples, type = "regular")

  # If regular sampling fails (odd shapes), use random
  if (length(sample_points) < n_samples) {
    log_message("  Using random sampling (regular failed)", "INFO")
    sample_points <- st_sample(stratum_sf, size = n_samples, type = "random")
  }

  # Convert to data frame
  sample_coords <- st_coordinates(sample_points)

  sample_df <- data.frame(
    sample_id = sprintf("%s_%03d", gsub(" ", "_", stratum_name), 1:nrow(sample_coords)),
    longitude = sample_coords[, 1],
    latitude = sample_coords[, 2],
    uncertainty_stratum = stratum_id,
    stratum_label = stratum_name,
    allocation_method = "neyman"
  )

  all_sample_points[[i]] <- sample_df

  log_message(sprintf("  Generated %d spatially-balanced points", nrow(sample_df)))
}

# Combine all sample points
sample_locations <- bind_rows(all_sample_points)

log_message(sprintf("\nTotal sample locations generated: %d", nrow(sample_locations)))

# ============================================================================
# EXTRACT PRIOR VALUES AT SAMPLE LOCATIONS
# ============================================================================

log_message("\nExtracting prior values at sample locations...")

# Convert to spatial points
sample_sf <- st_as_sf(sample_locations,
                     coords = c("longitude", "latitude"),
                     crs = crs(surface_mean, proj = TRUE))

# Extract values
sample_locations$prior_soc_mean_gkg <- extract(surface_mean, vect(sample_sf))[, 2]
sample_locations$prior_soc_se_gkg <- extract(surface_se, vect(sample_sf))[, 2]
sample_locations$prior_cv_pct <- extract(cv_raster, vect(sample_sf))[, 2]

# Add metadata
sample_locations <- sample_locations %>%
  mutate(
    sampling_priority = case_when(
      uncertainty_stratum == 3 ~ "High",
      uncertainty_stratum == 2 ~ "Medium",
      uncertainty_stratum == 1 ~ "Low",
      TRUE ~ "Unknown"
    ),
    target_depth_cm = 100,  # VM0033 standard
    notes = sprintf("Neyman allocation - %s uncertainty zone", stratum_label)
  )

# ============================================================================
# SAVE SAMPLE LOCATIONS
# ============================================================================

log_message("\nSaving sample locations...")

# Save as CSV
write_csv(sample_locations, "sampling_locations_neyman.csv")
log_message("Saved: sampling_locations_neyman.csv")

# Save as shapefile for GIS
sample_sf_export <- st_as_sf(sample_locations,
                             coords = c("longitude", "latitude"),
                             crs = crs(surface_mean, proj = TRUE))

st_write(sample_sf_export, "sampling_locations_neyman.gpkg",
        delete_dsn = TRUE, quiet = TRUE)
log_message("Saved: sampling_locations_neyman.gpkg (GIS format)")

# ============================================================================
# CREATE VISUALIZATION
# ============================================================================

log_message("\nCreating sampling design map...")

# Convert strata to data frame for plotting
strata_df <- as.data.frame(uncertainty_strata, xy = TRUE)
colnames(strata_df) <- c("x", "y", "stratum")
strata_df$stratum_label <- factor(strata_df$stratum,
                                  levels = 1:3,
                                  labels = c("Low", "Medium", "High"))

# Create plot
p <- ggplot() +
  geom_raster(data = strata_df, aes(x = x, y = y, fill = stratum_label)) +
  scale_fill_manual(values = c("Low" = "#2ecc71", "Medium" = "#f39c12", "High" = "#e74c3c"),
                   name = "Prior Uncertainty") +
  geom_point(data = sample_locations, aes(x = longitude, y = latitude),
            color = "black", size = 2, shape = 21, fill = "white", stroke = 1.5) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = sprintf("Bayesian Sampling Design - Neyman Allocation (%d samples)", nrow(sample_locations)),
    subtitle = sprintf("%s - %s", PROJECT_NAME, PROJECT_LOCATION),
    x = "Easting (m)",
    y = "Northing (m)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

ggsave("sampling_map_neyman.png", p, width = 10, height = 8, dpi = 300)
log_message("Saved: sampling_map_neyman.png")

# ============================================================================
# COMPARISON WITH UNIFORM SAMPLING
# ============================================================================

log_message("\n=== EFFICIENCY COMPARISON ===")

# Expected variance reduction compared to uniform sampling
# Neyman optimal variance: V_neyman = (Σ N_h * σ_h)^2 / n
# Uniform sampling variance: V_uniform = Σ (N_h/N)^2 * σ_h^2 * n

total_area <- sum(stratum_allocation$area_ha)
n_total <- sum(stratum_allocation$samples_allocated)

# Neyman variance
v_neyman <- (sum(stratum_allocation$area_ha * stratum_allocation$sd_soc) / total_area)^2 / n_total

# Uniform variance
stratum_allocation <- stratum_allocation %>%
  mutate(
    uniform_samples = (area_ha / total_area) * n_total,
    var_contrib = (area_ha / total_area)^2 * (sd_soc^2 / uniform_samples)
  )

v_uniform <- sum(stratum_allocation$var_contrib, na.rm = TRUE)

efficiency_gain <- ((v_uniform - v_neyman) / v_uniform) * 100

log_message(sprintf("Neyman allocation is %.1f%% more efficient than uniform sampling",
                   efficiency_gain))
log_message(sprintf("Equivalent to %.1f%% reduction in sampling effort for same precision",
                   efficiency_gain))

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("BAYESIAN SAMPLING DESIGN COMPLETE\n")
cat("========================================\n\n")

cat(sprintf("Total samples allocated: %d\n", nrow(sample_locations)))
cat(sprintf("Study area: %.1f ha\n", total_area))
cat(sprintf("Overall sampling density: %.2f samples/ha\n\n", nrow(sample_locations) / total_area))

cat("Allocation by uncertainty stratum:\n")
for (i in 1:nrow(stratum_allocation)) {
  cat(sprintf("  %s: %d samples (%.1f%%)\n",
              stratum_allocation$stratum_label[i],
              stratum_allocation$samples_allocated[i],
              100 * stratum_allocation$samples_allocated[i] / nrow(sample_locations)))
}

cat(sprintf("\nEfficiency gain vs uniform: %.1f%%\n\n", efficiency_gain))

cat("Output files:\n")
cat("  - sampling_locations_neyman.csv (sample coordinates)\n")
cat("  - sampling_locations_neyman.gpkg (GIS format)\n")
cat("  - sampling_allocation_neyman.csv (allocation table)\n")
cat("  - sampling_map_neyman.png (visualization)\n\n")

cat("NEXT STEPS:\n")
cat("1. Review sampling_map_neyman.png to verify spatial coverage\n")
cat("2. Export sampling_locations_neyman.gpkg to field GPS units\n")
cat("3. Conduct field sampling at allocated locations\n")
cat("4. Run Modules 01-05 with collected data (standard workflow)\n")
cat("5. Run Module 06C for Bayesian posterior estimation\n")
cat("\n")

log_message("=== MODULE 01C COMPLETE ===")
