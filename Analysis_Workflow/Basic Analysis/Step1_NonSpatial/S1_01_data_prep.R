# ============================================================================
# MODULE 01: BLUE CARBON DATA PREPARATION
# ============================================================================
# PURPOSE: Load, clean, validate, and structure core data for blue carbon monitoring
#
# INPUTS:
#   - Pre-Analysis Data Preparation/data_raw/core_locations.csv (GPS + stratum assignments)
#   - Pre-Analysis Data Preparation/data_raw/core_samples.csv (depth profiles + SOC)
#   - Analysis_Workflow/blue_carbon_config.R (configuration)
#
# OUTPUTS:
#   Core Data (data_processed/):
#     - cores_clean_bluecarbon.rds/.csv (cleaned sample data)
#     - core_totals.rds/.csv (integrated carbon stocks by core)
#     - cores_summary_by_stratum.csv (stratum statistics)
#     - carbon_by_stratum_summary.csv (carbon summaries)
#
#   Diagnostics (diagnostics/data_prep/):
#     - sampling_adequacy_report.csv / vm0033_compliance_report.csv (sampling sufficiency check)
#     - core_depth_completeness.csv (depth coverage by core)
#     - depth_completeness_summary.csv (overall depth coverage)
#     - core_type_summary.csv (HR vs Composite comparison)
#     - core_type_statistical_tests.csv (statistical test results)
#
#   QA/QC (diagnostics/qaqc/):
#     - bd_transparency_report.csv (bulk density measured vs estimated)
#     - qa_report.rds (comprehensive QA summary)
#
#   Spatial Readiness Check (console only):
#     - Reports whether AOI_FILE and STRATA_DIR are configured and present.
#     - Provides instructions for enabling zonal statistics in Module 04.
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
    setwd(dirname(dirname(dirname(dirname(this_file)))))
    return()
  }
  # Method 2: running interactively in RStudio — rstudioapi gives active file
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NULL)
    if (!is.null(active) && nchar(active) > 0) {
      setwd(dirname(dirname(dirname(dirname(normalizePath(active))))))
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
  stop("Configuration file not found. Run 00_setup_bluecarbon.R first.")
}

# Verify required config variables are loaded
required_vars <- c("VM0033_MIN_CORES", "CONFIDENCE_LEVEL", "VALID_STRATA",
                   "ACADEMIC_CONFIDENCE_LEVEL", "ACADEMIC_MARGIN_OF_ERROR", "ACADEMIC_DEFAULT_CV",
                   "INPUT_CRS", "PROCESSING_CRS", "BD_DEFAULTS")
missing_vars <- required_vars[!sapply(required_vars, exists)]
if (length(missing_vars) > 0) {
  stop(sprintf("Configuration error: Missing required variables: %s\nPlease check blue_carbon_config.R",
               paste(missing_vars, collapse=", ")))
}

# Create required output directories
required_dirs <- c(
  "logs",
  "data_processed",
  "diagnostics/data_prep",
  "diagnostics/qaqc"
)

for (dir in required_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
}

# Initialize logging
log_file <- file.path("logs", paste0("data_prep_", Sys.Date(), ".log"))

log_message <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- sprintf("[%s] %s: %s", timestamp, level, msg)
  cat(log_entry, "\n")
  cat(log_entry, "\n", file = log_file, append = TRUE)
}

log_message("=== MODULE 01: BLUE CARBON DATA PREPARATION ===")

# Load packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(sf)
})

log_message("Packages loaded successfully")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Validate stratum names against valid list
validate_strata <- function(strata_vector, valid_strata = VALID_STRATA) {
  invalid <- setdiff(unique(strata_vector), valid_strata)
  if (length(invalid) > 0) {
    warning(sprintf("Invalid strata detected: %s", paste(invalid, collapse = ", ")))
    cat("\nValid strata options:\n")
    for (s in valid_strata) {
      cat("  -", s, "\n")
    }
    return(FALSE)
  }
  return(TRUE)
}

#' Validate core coordinates for spatial quality control
#'
#' Performs comprehensive validation of GPS coordinates including:
#' - Range validation (valid lat/lon)
#' - NA detection
#' - Duplicate location detection
#' - Spatial clustering analysis
#'
#' @param locations Data frame with longitude, latitude, core_id columns
#' @return Cleaned locations data frame with validation flags
#'
#' @details
#' Removes cores with:
#' - Invalid coordinates (outside valid range or NA)
#' - Coordinates outside QC thresholds
#'
#' Warns about:
#' - Multiple cores at same location (GPS precision issues)
#' - Spatially clustered cores (<10m apart)
#'
#' @examples
#' locations <- validate_coordinates(locations)
validate_coordinates <- function(locations) {
  n_initial <- nrow(locations)

  # Check for impossible coordinates
  invalid_coords <- locations %>%
    filter(
      is.na(longitude) | is.na(latitude) |
      longitude < QC_LON_MIN | longitude > QC_LON_MAX |
      latitude < QC_LAT_MIN | latitude > QC_LAT_MAX
    )

  if (nrow(invalid_coords) > 0) {
    log_message(sprintf("Found %d cores with invalid coordinates - removing",
                       nrow(invalid_coords)), level = "WARNING")
    log_message(sprintf("  Invalid core IDs: %s",
                       paste(head(invalid_coords$core_id, 5), collapse = ", ")),
               level = "WARNING")
    if (nrow(invalid_coords) > 5) {
      log_message(sprintf("  ... and %d more", nrow(invalid_coords) - 5), level = "WARNING")
    }

    # Remove invalid
    locations <- locations %>%
      filter(!(core_id %in% invalid_coords$core_id))
  }

  # Check for duplicate/clustered cores (same location within GPS precision)
  # Round to 5 decimal places (~1m precision)
  locations <- locations %>%
    mutate(
      lon_rounded = round(longitude, 5),
      lat_rounded = round(latitude, 5)
    ) %>%
    group_by(lon_rounded, lat_rounded) %>%
    mutate(n_at_location = n()) %>%
    ungroup()

  if (any(locations$n_at_location > 1)) {
    n_clustered <- sum(locations$n_at_location > 1)
    n_unique_clusters <- locations %>%
      filter(n_at_location > 1) %>%
      distinct(lon_rounded, lat_rounded) %>%
      nrow()

    log_message(sprintf("%d cores at %d duplicate locations - check GPS accuracy",
                       n_clustered, n_unique_clusters), level = "WARNING")

    # Save duplicate locations report
    duplicate_report <- locations %>%
      filter(n_at_location > 1) %>%
      mutate(
        lon_rounded = round(longitude, 5),
        lat_rounded = round(latitude, 5)
      ) %>%
      select(core_id, longitude, latitude, stratum, n_at_location, lon_rounded, lat_rounded) %>%
      arrange(lon_rounded, lat_rounded)

    dup_file <- file.path("diagnostics/qaqc",
                          sprintf("duplicate_locations_%s_%s.csv",
                                 PROJECT_SCENARIO, MONITORING_YEAR))
    write.csv(duplicate_report, dup_file, row.names = FALSE)
    log_message(sprintf("  Duplicate locations saved to: %s", dup_file), level = "INFO")
  }

  # Clean up temporary columns
  locations <- locations %>%
    select(-lon_rounded, -lat_rounded, -n_at_location)

  n_removed <- n_initial - nrow(locations)
  if (n_removed > 0) {
    log_message(sprintf("Coordinate validation removed %d cores (%d remaining)",
                       n_removed, nrow(locations)), level = "INFO")
  } else {
    log_message("All coordinates passed validation", level = "INFO")
  }

  return(locations)
}

#' Calculate SOC stock for a depth increment
#'
#' Standard formula: Stock (kg/m²) = SOC (g/kg) × BD (g/cm³) × depth (cm) / 1000
#'
#' @param soc_g_kg Soil organic carbon content (g/kg)
#' @param bd_g_cm3 Bulk density (g/cm³)
#' @param depth_top_cm Top depth (cm)
#' @param depth_bottom_cm Bottom depth (cm)
#' @return Carbon stock in kg/m²
#'
#' @details
#' Dimensional analysis:
#' SOC (g C / kg soil) × BD (g soil / cm³) × depth (cm) / 1000
#' = g C × g / (kg × cm³) × cm / 1000
#' = g² C / (kg × cm²) / 1000
#' For 1 m² = 10,000 cm²:
#' = SOC × BD × depth × 10,000 cm² / (1000 kg/g) / 1000
#' = SOC × BD × depth / 1000 kg/m²
#'
#' Note: To convert to Mg/ha for reporting, multiply by 10
#' (since 1 kg/m² = 10 Mg/ha)
calculate_soc_stock <- function(soc_g_kg, bd_g_cm3, depth_top_cm, depth_bottom_cm) {
  depth_increment <- depth_bottom_cm - depth_top_cm
  soc_stock_kg_m2 <- soc_g_kg * bd_g_cm3 * depth_increment / 100
  return(soc_stock_kg_m2)
}

#' Calculate SOC stock with full uncertainty propagation
#'
#' Propagates uncertainty from both SOC and bulk density measurements
#' using first-order Taylor approximation for error propagation.
#'
#' @param soc_g_kg Soil organic carbon content (g/kg)
#' @param soc_se Standard error of SOC (g/kg)
#' @param bd_g_cm3 Bulk density (g/cm³)
#' @param bd_se Standard error of bulk density (g/cm³)
#' @param depth_top_cm Top depth (cm)
#' @param depth_bottom_cm Bottom depth (cm)
#'
#' @return List with mean and se in kg/m²
#'
#' @details
#' Error propagation formula for multiplication f = x * y:
#' Var(f) = f² * [(Var(x)/x²) + (Var(y)/y²)]
#' where relative variance = (SE/mean)²
#'
#' For SOC stock = SOC * BD * depth:
#' rel_var_stock = rel_var_soc + rel_var_bd
#'
#' @examples
#' stock <- calculate_soc_stock_with_uncertainty(50, 5, 1.2, 0.1, 0, 15)
#' # Returns: list(mean = 0.9, se = ~0.12)
calculate_soc_stock_with_uncertainty <- function(soc_g_kg, soc_se, bd_g_cm3, bd_se,
                                                 depth_top_cm, depth_bottom_cm) {
  depth_increment <- depth_bottom_cm - depth_top_cm

  # Mean stock using correct formula: SOC × BD × depth / 100
  stock_mean <- soc_g_kg * bd_g_cm3 * depth_increment / 100

  # Error propagation (first-order Taylor approximation)
  # If SE not provided, use conservative defaults: 10% for SOC, 15% for BD
  rel_var_soc <- if (!is.na(soc_se) && soc_g_kg > 0) {
    (soc_se / soc_g_kg)^2
  } else {
    0.1^2  # Conservative default: 10% CV
  }

  rel_var_bd <- if (!is.na(bd_se) && bd_g_cm3 > 0) {
    (bd_se / bd_g_cm3)^2
  } else {
    0.15^2  # Conservative default: 15% CV
  }

  # Combined relative variance (assuming independence)
  stock_se <- stock_mean * sqrt(rel_var_soc + rel_var_bd)

  return(list(mean = stock_mean, se = stock_se))
}

#' Assign bulk density defaults by stratum if missing
assign_bd_defaults <- function(df, bd_col = "bulk_density_g_cm3",
                               stratum_col = "stratum") {
  df[[bd_col]] <- ifelse(
    is.na(df[[bd_col]]),
    sapply(df[[stratum_col]], function(s) {
      if (s %in% names(BD_DEFAULTS)) {
        BD_DEFAULTS[[s]]
      } else {
        1.0  # Generic default
      }
    }),
    df[[bd_col]]
  )
  return(df)
}

#' Calculate required sample size for statistical adequacy
#' Based on: n = (z * CV / target_precision)^2
calculate_required_n <- function(cv, target_precision = VM0033_TARGET_PRECISION,
                                confidence = CONFIDENCE_LEVEL) {
  z <- qnorm(1 - (1 - confidence) / 2)  # 1.96 for 95% CI
  n <- ceiling((z * cv / target_precision)^2)
  return(max(n, VM0033_MIN_CORES))  # Ensure at least minimum
}

#' Calculate achieved precision from sample size and CV
calculate_achieved_precision <- function(n, cv, confidence = CONFIDENCE_LEVEL) {
  if (n < 2) return(NA)
  z <- qnorm(1 - (1 - confidence) / 2)
  precision <- (z * cv) / sqrt(n)
  return(precision)
}

#' Calculate Cochran sample size recommendation for continuous data
#' n = (Z^2 * sigma^2) / e^2
calculate_cochran_n_continuous <- function(cv_percent,
                                           confidence = ACADEMIC_CONFIDENCE_LEVEL,
                                           margin_of_error = ACADEMIC_MARGIN_OF_ERROR,
                                           default_cv = ACADEMIC_DEFAULT_CV) {
  z <- qnorm(1 - (1 - confidence) / 2)
  sigma <- ifelse(is.na(cv_percent) || cv_percent <= 0, default_cv, cv_percent / 100)
  ceiling((z^2 * sigma^2) / (margin_of_error^2))
}

#' Calculate depth profile completeness (0-100%)
calculate_profile_completeness <- function(depth_top, depth_bottom, max_depth = MAX_CORE_DEPTH) {
  # Calculate total depth sampled
  total_sampled <- sum(depth_bottom - depth_top, na.rm = TRUE)
  completeness_pct <- (total_sampled / max_depth) * 100
  return(min(completeness_pct, 100))  # Cap at 100%
}

# ============================================================================
# LOAD CORE LOCATIONS
# ============================================================================

log_message("Loading core locations...")

locations_file <- CORE_LOCATIONS_FILE   # set in blue_carbon_config.R

if (!file.exists(locations_file)) {
  stop(sprintf("Core locations file not found: %s\n  → Place core_locations.csv in: %s",
               locations_file, DATA_RAW_DIR))
}

# Load with column name standardization
locations <- read_csv(locations_file, show_col_types = FALSE) %>%
  rename_with(tolower)

log_message(sprintf("Loaded %d core locations", nrow(locations)))

# Check required columns
required_cols_locations <- c("core_id", "longitude", "latitude", "stratum")
missing_cols <- setdiff(required_cols_locations, names(locations))

if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns in locations: %s", 
               paste(missing_cols, collapse = ", ")))
}

# Validate coordinates with comprehensive QC
locations <- validate_coordinates(locations)

# Validate strata
if (!validate_strata(locations$stratum)) {
  stop("Invalid stratum names detected. Please fix in source data.")
}

# Add VM0033 metadata if not present
if (!"scenario_type" %in% names(locations)) {
  locations$scenario_type <- PROJECT_SCENARIO
  log_message("Added scenario_type from config")
}

if (!"monitoring_year" %in% names(locations)) {
  locations$monitoring_year <- MONITORING_YEAR
  log_message("Added monitoring_year from config")
}

# Validate scenario_type (optional validation if VALID_SCENARIOS exists)
if ("scenario_type" %in% names(locations) && exists("VALID_SCENARIOS")) {
  invalid_scenarios <- setdiff(unique(locations$scenario_type), VALID_SCENARIOS)

  if (length(invalid_scenarios) > 0) {
    log_message(sprintf("WARNING: Invalid scenario types found: %s",
                       paste(invalid_scenarios, collapse = ", ")), "WARNING")
    log_message(sprintf("  Valid options: %s", paste(VALID_SCENARIOS, collapse = ", ")), "WARNING")
    log_message("  These will cause errors in temporal analysis modules (Module 08/09)", "WARNING")
  }
}

# Validate monitoring_year
if ("monitoring_year" %in% names(locations)) {
  current_year <- as.integer(format(Sys.Date(), "%Y"))
  invalid_years <- locations$monitoring_year[locations$monitoring_year > current_year |
                                               locations$monitoring_year < 1900]

  if (length(invalid_years) > 0) {
    log_message("WARNING: Suspicious monitoring years detected", "WARNING")
    log_message(sprintf("  Years outside reasonable range (1900-%d): %s",
                       current_year,
                       paste(unique(invalid_years), collapse = ", ")), "WARNING")
  }
}

# Ensure core_type exists
if (!"core_type" %in% names(locations)) {
  locations$core_type <- "unknown"
  log_message("core_type not specified, set to 'unknown'", "WARNING")
}

# Convert to spatial object
locations_sf <- st_as_sf(locations, 
                         coords = c("longitude", "latitude"),
                         crs = INPUT_CRS,
                         remove = FALSE)

# Transform to processing CRS
locations_sf <- st_transform(locations_sf, PROCESSING_CRS)

log_message(sprintf("Created spatial object with CRS %d", PROCESSING_CRS))

# ============================================================================
# LOAD CORE SAMPLES
# ============================================================================

log_message("Loading core samples...")

samples_file <- CORE_SAMPLES_FILE   # set in blue_carbon_config.R

if (!file.exists(samples_file)) {
  stop(sprintf("Core samples file not found: %s\n  → Place core_samples.csv in: %s",
               samples_file, DATA_RAW_DIR))
}

# Load with column name standardization
samples <- read_csv(samples_file, show_col_types = FALSE) %>%
  rename_with(tolower)

log_message(sprintf("Loaded %d samples", nrow(samples)))

# Check required columns
required_cols_samples <- c("core_id", "depth_top_cm", "depth_bottom_cm", "soc_g_kg")
missing_cols <- setdiff(required_cols_samples, names(samples))

if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns in samples: %s", 
               paste(missing_cols, collapse = ", ")))
}

# Calculate depth midpoint
samples <- samples %>%
  mutate(
    depth_cm = (depth_top_cm + depth_bottom_cm) / 2,
    interval_thickness_cm = depth_bottom_cm - depth_top_cm
  )

log_message("Calculated depth midpoints and interval thickness")

# ============================================================================
# DATA QUALITY CHECKS
# ============================================================================

log_message("Running quality checks...")

# Initialize QA flags
samples <- samples %>%
  mutate(
    qa_depth_valid = depth_top_cm >= 0 & 
                     depth_top_cm < depth_bottom_cm &
                     depth_bottom_cm <= MAX_CORE_DEPTH,
    
    qa_soc_valid = !is.na(soc_g_kg) & 
                   soc_g_kg >= QC_SOC_MIN & 
                   soc_g_kg <= QC_SOC_MAX
  )

# Check bulk density if present
if ("bulk_density_g_cm3" %in% names(samples)) {
  samples <- samples %>%
    mutate(
      qa_bd_valid = is.na(bulk_density_g_cm3) | 
                    (bulk_density_g_cm3 >= QC_BD_MIN & 
                     bulk_density_g_cm3 <= QC_BD_MAX),
      bd_measured = !is.na(bulk_density_g_cm3)
    )
} else {
  log_message("No bulk_density_g_cm3 column - will use defaults", "WARNING")
  samples$bulk_density_g_cm3 <- NA
  samples$qa_bd_valid <- TRUE
  samples$bd_measured <- FALSE
}

# Report QA results
n_depth_invalid <- sum(!samples$qa_depth_valid)
n_soc_invalid <- sum(!samples$qa_soc_valid)
n_bd_invalid <- sum(!samples$qa_bd_valid)

if (n_depth_invalid > 0) {
  log_message(sprintf("Invalid depths: %d samples", n_depth_invalid), "WARNING")
}
if (n_soc_invalid > 0) {
  log_message(sprintf("Invalid SOC values: %d samples", n_soc_invalid), "WARNING")
}
if (n_bd_invalid > 0) {
  log_message(sprintf("Invalid BD values: %d samples", n_bd_invalid), "WARNING")
}

# Filter to valid samples only
samples_clean <- samples %>%
  filter(qa_depth_valid & qa_soc_valid & qa_bd_valid)

log_message(sprintf("After QA: %d samples retained from %d cores",
                    nrow(samples_clean),
                    n_distinct(samples_clean$core_id)))

# ============================================================================
# MERGE LOCATIONS WITH SAMPLES
# ============================================================================

log_message("Merging locations with samples...")

# Drop geometry for merging (we'll add it back)
locations_df <- locations_sf %>%
  st_drop_geometry()

# Merge
cores_merged <- samples_clean %>%
  left_join(locations_df, by = "core_id")

# Check for cores without locations
cores_no_location <- samples_clean %>%
  anti_join(locations_df, by = "core_id") %>%
  pull(core_id) %>%
  unique()

if (length(cores_no_location) > 0) {
  log_message(sprintf("Warning: %d cores have samples but no location", 
                      length(cores_no_location)), "WARNING")
  log_message(sprintf("Missing cores: %s", 
                      paste(head(cores_no_location, 5), collapse = ", ")))
}

# Check for locations without samples
cores_no_samples <- locations_df %>%
  anti_join(samples_clean, by = "core_id") %>%
  pull(core_id) %>%
  unique()

if (length(cores_no_samples) > 0) {
  log_message(sprintf("Warning: %d cores have location but no samples",
                      length(cores_no_samples)), "WARNING")
}

# Filter to complete cases only
cores_complete <- cores_merged %>%
  filter(!is.na(longitude) & !is.na(latitude) & !is.na(stratum))

log_message(sprintf("Complete dataset: %d samples from %d cores",
                    nrow(cores_complete),
                    n_distinct(cores_complete$core_id)))

# ============================================================================
# STRATUM VALIDATION AND STATISTICS
# ============================================================================

log_message("Validating stratum assignments...")

# Validate all strata
if (!validate_strata(cores_complete$stratum)) {
  stop("Invalid stratum assignments in merged data")
}

# Calculate stratum statistics with uncertainty metrics (CHANGE #3)
stratum_stats <- cores_complete %>%
  group_by(stratum) %>%
  summarise(
    n_cores = n_distinct(core_id),
    n_samples = n(),
    mean_depth = mean(depth_cm),
    max_depth = max(depth_cm),
    mean_soc = mean(soc_g_kg, na.rm = TRUE),
    sd_soc = sd(soc_g_kg, na.rm = TRUE),
    cv_soc = (sd(soc_g_kg, na.rm = TRUE) / mean(soc_g_kg, na.rm = TRUE)) * 100,
    se_soc = sd(soc_g_kg, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cores))

log_message("Stratum summary with uncertainty metrics:")
print(stratum_stats)

# Save stratum summary
write_csv(stratum_stats, "data_processed/cores_summary_by_stratum.csv")
log_message("Saved stratum summary")

# ============================================================================
# STATISTICAL POWER & VM0033 ASSESSMENT
# ============================================================================

log_message(vm_label("Running Statistical Power & VM0033 Assessment...",
                     "Running Statistical Adequacy Assessment..."))

# Calculate power analysis for each stratum
vm0033_compliance <- stratum_stats %>%
  mutate(
    # Baseline requirements
    vm0033_min_requirement = VM0033_MIN_CORES,
    meets_vm0033_min_n = n_cores >= vm0033_min_requirement,

    # Calculate required N for target precision
    required_n_20pct = mapply(calculate_required_n, cv_soc, 20, CONFIDENCE_LEVEL),
    required_n_15pct = mapply(calculate_required_n, cv_soc, 15, CONFIDENCE_LEVEL),
    required_n_10pct = mapply(calculate_required_n, cv_soc, 10, CONFIDENCE_LEVEL),

    # Calculate achieved precision with current n
    achieved_precision_pct = mapply(calculate_achieved_precision, n_cores, cv_soc, CONFIDENCE_LEVEL),

    # Additional cores needed for VM0033 precision targets
    additional_for_20pct = pmax(0, required_n_20pct - n_cores),
    additional_for_15pct = pmax(0, required_n_15pct - n_cores),
    additional_for_10pct = pmax(0, required_n_10pct - n_cores),

    # Statistically robust recommendation (Cochran's formula)
    cochran_sigma_cv = ifelse(is.na(cv_soc) | cv_soc <= 0, ACADEMIC_DEFAULT_CV, cv_soc / 100),
    required_n_cochran = mapply(
      calculate_cochran_n_continuous,
      cv_soc,
      ACADEMIC_CONFIDENCE_LEVEL,
      ACADEMIC_MARGIN_OF_ERROR,
      ACADEMIC_DEFAULT_CV
    ),
    additional_for_cochran = pmax(0, required_n_cochran - n_cores),

    # VM0033 compliance flags
    meets_20pct_precision = achieved_precision_pct <= 20,
    meets_15pct_precision = achieved_precision_pct <= 15,
    meets_10pct_precision = achieved_precision_pct <= 10,

    # Overall VM0033 baseline assessment (min 3 cores AND ≤20% precision)
    vm0033_compliant = meets_vm0033_min_n & meets_20pct_precision,
    meets_cochran_recommendation = n_cores >= required_n_cochran,

    # Status assessment
    status = case_when(
      !meets_vm0033_min_n ~ "INSUFFICIENT (< 3 cores)",
      achieved_precision_pct <= 10 ~ "EXCELLENT (≤10%)",
      achieved_precision_pct <= 15 ~ "GOOD (≤15%)",
      achieved_precision_pct <= 20 ~ "ACCEPTABLE (≤20%)",
      achieved_precision_pct <= 30 ~ "MARGINAL (>20%, <30%)",
      TRUE ~ "POOR (≥30%)"
    ),
    assessment_summary = case_when(
      vm0033_compliant & meets_cochran_recommendation ~ "Meets VM0033 baseline and robust statistical recommendation",
      vm0033_compliant & !meets_cochran_recommendation ~ "Meets VM0033 baseline but under Cochran robust recommendation",
      !vm0033_compliant & meets_cochran_recommendation ~ "Meets Cochran recommendation but misses VM0033 precision baseline",
      TRUE ~ "Needs additional sampling under both frameworks"
    )
  )

cat("\n========================================\n")
cat(vm_label("STATISTICAL POWER & VM0033 ASSESSMENT\n",
             "STATISTICAL ADEQUACY ASSESSMENT\n"))
cat("========================================\n\n")

# Print adequacy by stratum
for (i in 1:nrow(vm0033_compliance)) {
  cat(sprintf("Stratum: %s\n", vm0033_compliance$stratum[i]))
  cat(sprintf("  Current samples: %d cores\n", vm0033_compliance$n_cores[i]))
  cat(sprintf("  CV: %.1f%%\n", vm0033_compliance$cv_soc[i]))
  cat(sprintf("  Achieved precision: %.1f%% (at %.0f%% CI)\n",
              vm0033_compliance$achieved_precision_pct[i],
              CONFIDENCE_LEVEL * 100))
  cat(sprintf("  %s: %d cores\n",
              vm_label("VM0033 Minimum Requirement",
                       "Statistical best-practice minimum"),
              vm0033_compliance$vm0033_min_requirement[i]))
  cat(sprintf("  Statistically Robust Recommendation (Cochran's): %d cores (sigma/CV=%.2f, e=%.2f)\n",
              vm0033_compliance$required_n_cochran[i],
              vm0033_compliance$cochran_sigma_cv[i],
              ACADEMIC_MARGIN_OF_ERROR))
  cat(sprintf("  Status: %s\n", vm0033_compliance$status[i]))
  cat(sprintf("  %s: %s\n",
              vm_label("VM0033 Baseline", "Adequacy threshold (n≥3, ≤20% precision)"),
              ifelse(vm0033_compliance$vm0033_compliant[i], "✓ YES", "✗ NO")))
  cat(sprintf("  Cochran Robust Recommendation Met: %s\n",
              ifelse(vm0033_compliance$meets_cochran_recommendation[i], "✓ YES", "✗ NO")))
  cat(sprintf("  Combined Assessment: %s\n", vm0033_compliance$assessment_summary[i]))

  # Recommendations
  if (!vm0033_compliance$vm0033_compliant[i] || !vm0033_compliance$meets_cochran_recommendation[i]) {
    cat("\n  Recommendations:\n")
    if (!vm0033_compliance$meets_vm0033_min_n[i]) {
      cat(sprintf("    • Add %d cores to satisfy the statistical minimum requirement (n=%d per stratum)\n",
                  vm0033_compliance$vm0033_min_requirement[i] - vm0033_compliance$n_cores[i],
                  vm0033_compliance$vm0033_min_requirement[i]))
    }
    if (vm0033_compliance$additional_for_20pct[i] > 0) {
      cat(sprintf("    • Add %d cores to achieve 20%% precision\n",
                  vm0033_compliance$additional_for_20pct[i]))
    }
    if (vm0033_compliance$additional_for_15pct[i] > 0) {
      cat(sprintf("    • Add %d cores to achieve 15%% precision\n",
                  vm0033_compliance$additional_for_15pct[i]))
    }
    if (vm0033_compliance$additional_for_cochran[i] > 0) {
      cat(sprintf("    • Add %d cores to reach Cochran robust recommendation\n",
                  vm0033_compliance$additional_for_cochran[i]))
    }
  }
  cat("\n")
}

# Overall project status
n_compliant <- sum(vm0033_compliance$vm0033_compliant)
n_total <- nrow(vm0033_compliance)

cat(sprintf("%s: %d/%d strata meet requirements\n",
            vm_label("Overall VM0033 Baseline", "Overall statistical adequacy"),
            n_compliant, n_total))
cat(sprintf("Overall Cochran Robust Recommendation: %d/%d strata meet recommendation\n\n",
            sum(vm0033_compliance$meets_cochran_recommendation), n_total))

if (n_compliant < n_total) {
  log_message(sprintf("WARNING: %d strata do not meet the statistical adequacy threshold",
                      n_total - n_compliant), "WARNING")
}

if (sum(vm0033_compliance$meets_cochran_recommendation) < n_total) {
  log_message(sprintf("WARNING: %d strata are below Cochran robust recommendation",
                      n_total - sum(vm0033_compliance$meets_cochran_recommendation)), "WARNING")
}

# Save sampling adequacy report
adequacy_report_file <- file.path(
  "diagnostics/data_prep",
  vm_label("vm0033_compliance_report.csv", "sampling_adequacy_report.csv")
)
write_csv(vm0033_compliance, adequacy_report_file)
log_message(sprintf("Saved %s",
                    vm_label("Statistical Power & VM0033 assessment report",
                             "Statistical adequacy assessment report")))

# ============================================================================
# BULK DENSITY HANDLING
# ============================================================================

log_message("Handling bulk density...")

n_bd_missing <- sum(is.na(cores_complete$bulk_density_g_cm3))
n_bd_measured <- sum(!is.na(cores_complete$bulk_density_g_cm3))

log_message(sprintf("BD measured: %d samples", n_bd_measured))
log_message(sprintf("BD missing: %d samples", n_bd_missing))

if (n_bd_missing > 0) {
  log_message("Applying stratum-specific BD defaults to missing values")

  # Show defaults being applied
  cat("\nBulk density defaults by stratum:\n")
  for (s in names(BD_DEFAULTS)) {
    cat(sprintf("  %s: %.2f g/cm³\n", s, BD_DEFAULTS[[s]]))
  }

  cores_complete <- assign_bd_defaults(cores_complete)

  # Flag which samples have estimated BD
  cores_complete <- cores_complete %>%
    mutate(bd_estimated = !bd_measured)
}

# ============================================================================
# BULK DENSITY TRANSPARENCY REPORT (CHANGE #4)
# ============================================================================

log_message("Generating bulk density transparency report...")

# Calculate BD statistics by stratum
bd_transparency <- cores_complete %>%
  group_by(stratum) %>%
  summarise(
    n_samples = n(),
    n_measured = sum(bd_measured),
    n_estimated = sum(!bd_measured),
    pct_measured = (sum(bd_measured) / n()) * 100,
    pct_estimated = (sum(!bd_measured) / n()) * 100,

    # Measured BD stats (where available)
    mean_bd_measured = ifelse(sum(bd_measured) > 0,
                               mean(bulk_density_g_cm3[bd_measured], na.rm = TRUE),
                               NA),
    sd_bd_measured = ifelse(sum(bd_measured) > 1,
                            sd(bulk_density_g_cm3[bd_measured], na.rm = TRUE),
                            NA),

    # Estimated BD (from defaults)
    mean_bd_estimated = ifelse(sum(!bd_measured) > 0,
                                mean(bulk_density_g_cm3[!bd_measured], na.rm = TRUE),
                                NA),

    # Overall BD
    mean_bd_all = mean(bulk_density_g_cm3, na.rm = TRUE),

    .groups = "drop"
  )

cat("\n========================================\n")
cat("BULK DENSITY TRANSPARENCY REPORT\n")
cat("========================================\n\n")

cat(sprintf("Total samples: %d\n", nrow(cores_complete)))
cat(sprintf("  Measured BD: %d (%.1f%%)\n",
            sum(cores_complete$bd_measured),
            100 * sum(cores_complete$bd_measured) / nrow(cores_complete)))
cat(sprintf("  Estimated BD: %d (%.1f%%)\n\n",
            sum(!cores_complete$bd_measured),
            100 * sum(!cores_complete$bd_measured) / nrow(cores_complete)))

cat("By stratum:\n")
for (i in 1:nrow(bd_transparency)) {
  cat(sprintf("\n%s:\n", bd_transparency$stratum[i]))
  cat(sprintf("  Measured: %d/%d (%.1f%%)\n",
              bd_transparency$n_measured[i],
              bd_transparency$n_samples[i],
              bd_transparency$pct_measured[i]))

  if (!is.na(bd_transparency$mean_bd_measured[i])) {
    cat(sprintf("  Mean measured BD: %.2f ± %.2f g/cm³\n",
                bd_transparency$mean_bd_measured[i],
                ifelse(is.na(bd_transparency$sd_bd_measured[i]), 0,
                       bd_transparency$sd_bd_measured[i])))
  }

  if (bd_transparency$n_estimated[i] > 0) {
    cat(sprintf("  Estimated BD (default): %.2f g/cm³\n",
                bd_transparency$mean_bd_estimated[i]))
  }

  cat(sprintf("  Overall mean BD: %.2f g/cm³\n",
              bd_transparency$mean_bd_all[i]))
}

cat("\n📝 Note: Estimated BD values are based on literature defaults.\n")
cat("   Carbon stock uncertainty will be higher for samples with estimated BD.\n")
cat("   Measuring BD for all cores is strongly recommended when possible.\n\n")

# Save BD transparency report
write_csv(bd_transparency, "diagnostics/qaqc/bd_transparency_report.csv")
log_message("Saved BD transparency report")

# ============================================================================
# CALCULATE CARBON STOCKS
# ============================================================================

log_message("Calculating carbon stocks...")

cores_complete <- cores_complete %>%
  mutate(
    # Carbon stock per sample (kg C/m²)
    # Standard unit for carbon stock reporting and aligns with prior data
    carbon_stock_kg_m2 = calculate_soc_stock(
      soc_g_kg,
      bulk_density_g_cm3,
      depth_top_cm,
      depth_bottom_cm
    )
  )

# Calculate total stocks per core
core_totals <- cores_complete %>%
  group_by(core_id, stratum) %>%
  summarise(
    total_carbon_stock = sum(carbon_stock_kg_m2, na.rm = TRUE),
    max_depth_sampled = max(depth_bottom_cm),
    n_samples = n(),
    .groups = "drop"
  )

log_message(sprintf("Calculated carbon stocks for %d cores", nrow(core_totals)))

# Summary by stratum
carbon_by_stratum <- core_totals %>%
  group_by(stratum) %>%
  summarise(
    n_cores = n(),
    mean_stock = mean(total_carbon_stock),
    sd_stock = sd(total_carbon_stock),
    min_stock = min(total_carbon_stock),
    max_stock = max(total_carbon_stock),
    .groups = "drop"
  )

log_message("Carbon stock summary by stratum:")
print(carbon_by_stratum)

# ============================================================================
# DEPTH PROFILE COMPLETENESS (CHANGE #7)
# ============================================================================

log_message("Calculating depth profile completeness...")

# Calculate completeness per core
core_depth_completeness <- cores_complete %>%
  group_by(core_id, stratum, core_type) %>%
  summarise(
    n_samples = n(),
    min_depth = min(depth_top_cm),
    max_depth = max(depth_bottom_cm),
    depth_range = max(depth_bottom_cm) - min(depth_top_cm),
    total_sampled = sum(depth_bottom_cm - depth_top_cm),
    completeness_pct = calculate_profile_completeness(depth_top_cm, depth_bottom_cm, MAX_CORE_DEPTH),

    # Check for depth gaps
    has_gaps = any(diff(sort(c(depth_top_cm, depth_bottom_cm))) > 5),

    # Classification
    profile_quality = case_when(
      completeness_pct >= 90 ~ "Complete (≥90%)",
      completeness_pct >= 70 ~ "Good (70-89%)",
      completeness_pct >= 50 ~ "Moderate (50-69%)",
      TRUE ~ "Incomplete (<50%)"
    ),

    .groups = "drop"
  )

# Summary by stratum
depth_completeness_summary <- core_depth_completeness %>%
  group_by(stratum) %>%
  summarise(
    n_cores = n(),
    mean_completeness = mean(completeness_pct),
    sd_completeness = sd(completeness_pct),
    min_completeness = min(completeness_pct),
    max_completeness = max(completeness_pct),
    n_complete = sum(completeness_pct >= 90),
    n_good = sum(completeness_pct >= 70 & completeness_pct < 90),
    n_moderate = sum(completeness_pct >= 50 & completeness_pct < 70),
    n_incomplete = sum(completeness_pct < 50),
    .groups = "drop"
  )

cat("\n========================================\n")
cat("DEPTH PROFILE COMPLETENESS\n")
cat("========================================\n\n")

for (i in 1:nrow(depth_completeness_summary)) {
  cat(sprintf("%s:\n", depth_completeness_summary$stratum[i]))
  cat(sprintf("  Mean completeness: %.1f%% ± %.1f%%\n",
              depth_completeness_summary$mean_completeness[i],
              depth_completeness_summary$sd_completeness[i]))
  cat(sprintf("  Range: %.1f%% - %.1f%%\n",
              depth_completeness_summary$min_completeness[i],
              depth_completeness_summary$max_completeness[i]))
  cat(sprintf("  Complete profiles (≥90%%): %d/%d\n",
              depth_completeness_summary$n_complete[i],
              depth_completeness_summary$n_cores[i]))
  cat(sprintf("  Good profiles (70-89%%): %d/%d\n",
              depth_completeness_summary$n_good[i],
              depth_completeness_summary$n_cores[i]))
  if (depth_completeness_summary$n_incomplete[i] > 0) {
    cat(sprintf("  ⚠ Incomplete profiles (<50%%): %d\n",
                depth_completeness_summary$n_incomplete[i]))
  }
  cat("\n")
}

# Save depth completeness report
write_csv(core_depth_completeness, "diagnostics/data_prep/core_depth_completeness.csv")
write_csv(depth_completeness_summary, "diagnostics/data_prep/depth_completeness_summary.csv")
log_message("Saved depth completeness reports")

# ============================================================================
# CORE TYPE COMPARISON: HR vs COMPOSITE (CHANGE #6)
# ============================================================================

log_message("Analyzing HR vs Composite core differences...")

# Check if core_type column exists and has values
if ("core_type" %in% names(cores_complete)) {

  # Standardize core type names
  cores_complete <- cores_complete %>%
    mutate(
      core_type_clean = case_when(
        tolower(core_type) %in% c("hr", "high-res", "high resolution", "high res") ~ "HR",
        tolower(core_type) %in% c("paired composite", "paired comp", "paired") ~ "Paired Composite",
        tolower(core_type) %in% c("unpaired composite", "unpaired comp", "unpaired", "composite", "comp") ~ "Unpaired Composite",
        TRUE ~ "Other"
      )
    )

  # Count by core type
  core_type_counts <- cores_complete %>%
    group_by(core_type_clean) %>%
    summarise(
      n_cores = n_distinct(core_id),
      n_samples = n(),
      .groups = "drop"
    )

  cat("\n========================================\n")
  cat("CORE TYPE DISTRIBUTION\n")
  cat("========================================\n\n")
  print(core_type_counts)

  # Only proceed with comparison if we have both HR and Paired Composite cores
  has_hr <- any(cores_complete$core_type_clean == "HR")
  has_paired_comp <- any(cores_complete$core_type_clean == "Paired Composite")

  if (has_hr && has_paired_comp) {

    # Summary by core type and stratum (only HR and Paired Composite)
    core_type_summary <- cores_complete %>%
      filter(core_type_clean %in% c("HR", "Paired Composite")) %>%
      group_by(stratum, core_type_clean) %>%
      summarise(
        n_cores = n_distinct(core_id),
        n_samples = n(),
        mean_soc = mean(soc_g_kg, na.rm = TRUE),
        sd_soc = sd(soc_g_kg, na.rm = TRUE),
        se_soc = sd(soc_g_kg, na.rm = TRUE) / sqrt(n()),
        mean_bd = mean(bulk_density_g_cm3, na.rm = TRUE),
        .groups = "drop"
      )

    # Statistical tests (t-test) by stratum - HR vs Paired Composite only
    statistical_tests <- list()

    for (s in unique(cores_complete$stratum)) {
      hr_data <- cores_complete %>%
        filter(stratum == s, core_type_clean == "HR") %>%
        pull(soc_g_kg)

      comp_data <- cores_complete %>%
        filter(stratum == s, core_type_clean == "Paired Composite") %>%
        pull(soc_g_kg)

      if (length(hr_data) >= 3 && length(comp_data) >= 3) {
        test_result <- t.test(hr_data, comp_data)

        statistical_tests[[s]] <- data.frame(
          stratum = s,
          n_hr = length(hr_data),
          n_paired_comp = length(comp_data),
          hr_mean = mean(hr_data, na.rm = TRUE),
          comp_mean = mean(comp_data, na.rm = TRUE),
          diff = mean(hr_data, na.rm = TRUE) - mean(comp_data, na.rm = TRUE),
          t_statistic = test_result$statistic,
          p_value = test_result$p.value,
          significant = test_result$p.value < 0.05,
          ci_lower = test_result$conf.int[1],
          ci_upper = test_result$conf.int[2]
        )
      }
    }

    if (length(statistical_tests) > 0) {
      statistical_tests_df <- bind_rows(statistical_tests)

      cat("\n========================================\n")
      cat("HR vs PAIRED COMPOSITE CORE COMPARISON\n")
      cat("========================================\n\n")

      cat("Statistical Tests (Two-sample t-tests):\n")
      cat("Purpose: Validate if paired composite cores can serve as proxy for HR cores\n\n")

      for (i in 1:nrow(statistical_tests_df)) {
        cat(sprintf("%s:\n", statistical_tests_df$stratum[i]))
        cat(sprintf("  Sample sizes: %d HR, %d Paired Composite\n",
                    statistical_tests_df$n_hr[i],
                    statistical_tests_df$n_paired_comp[i]))
        cat(sprintf("  HR mean SOC: %.1f g/kg\n", statistical_tests_df$hr_mean[i]))
        cat(sprintf("  Paired Composite mean SOC: %.1f g/kg\n", statistical_tests_df$comp_mean[i]))
        cat(sprintf("  Difference: %.1f g/kg\n", statistical_tests_df$diff[i]))
        cat(sprintf("  t-statistic: %.2f\n", statistical_tests_df$t_statistic[i]))
        cat(sprintf("  p-value: %.4f %s\n",
                    statistical_tests_df$p_value[i],
                    ifelse(statistical_tests_df$significant[i], "**", "")))
        cat(sprintf("  95%% CI: [%.1f, %.1f]\n",
                    statistical_tests_df$ci_lower[i],
                    statistical_tests_df$ci_upper[i]))

        if (statistical_tests_df$significant[i]) {
          cat("  ✗ Significant difference detected (p < 0.05)\n")
          cat("  → HR and Paired Composite cores differ significantly\n")
          cat("  → Recommendation: Analyze HR and Composite separately OR collect more paired samples\n")
        } else {
          cat("  ✓ No significant difference (p ≥ 0.05)\n")
          cat("  → Paired sampling assumption SUPPORTED\n")
          cat("  → Recommendation: Composite cores can be used as proxy for HR cores in this stratum\n")
        }
        cat("\n")
      }

      # Overall recommendation
      n_supported <- sum(!statistical_tests_df$significant)
      n_total <- nrow(statistical_tests_df)

      cat(sprintf("Overall: %d/%d strata support paired sampling approach\n\n", n_supported, n_total))

      if (n_supported == n_total) {
        cat("✓ EXCELLENT: All strata show paired composites can proxy for HR cores\n")
        cat("  → Cost-effective approach validated for this site\n")
      } else if (n_supported > 0) {
        cat("⚠ MIXED: Some strata support paired approach, others don't\n")
        cat("  → Consider stratum-specific strategies\n")
      } else {
        cat("✗ WARNING: No strata support paired approach\n")
        cat("  → HR and composite cores should be analyzed separately\n")
      }

      # Save core type comparison
      write_csv(core_type_summary, "diagnostics/data_prep/core_type_summary.csv")
      write_csv(statistical_tests_df, "diagnostics/data_prep/core_type_statistical_tests.csv")
      log_message("Saved core type comparison reports")

    } else {
      log_message("Insufficient data for HR vs Paired Composite statistical tests (need ≥3 each)", "WARNING")
      cat("\n⚠ Statistical tests require ≥3 samples per core type per stratum\n")
      cat("  Current data doesn't meet this threshold\n\n")
    }

  } else {
    # Report what we have
    if (!has_hr) {
      log_message("No HR cores found in dataset", "WARNING")
      cat("\n⚠ No HR cores detected in dataset\n")
    }
    if (!has_paired_comp) {
      log_message("No Paired Composite cores found in dataset", "WARNING")
      cat("\n⚠ No Paired Composite cores detected in dataset\n")
    }
    cat("  HR vs Paired Composite comparison requires both core types\n")
    cat("  See Pre-Analysis Data Preparation/data_raw/README_DATA_STRUCTURE.md for data requirements\n\n")
  }

} else {
  log_message("core_type column not found in data", "WARNING")
  cat("\n⚠ Core type comparison skipped: core_type column not found\n")
  cat("  See Pre-Analysis Data Preparation/data_raw/README_DATA_STRUCTURE.md for required data structure\n\n")
}

# ============================================================================
# ADD FINAL QA FLAGS
# ============================================================================

log_message("Adding final QA flags...")

cores_complete <- cores_complete %>%
  mutate(
    # Spatial validity
    qa_spatial_valid = !is.na(longitude) & !is.na(latitude),
    
    # Stratum validity
    qa_stratum_valid = stratum %in% VALID_STRATA,
    
    # Overall QA pass
    qa_pass = qa_spatial_valid & qa_depth_valid & qa_soc_valid & 
              qa_bd_valid & qa_stratum_valid,
    
    # Sample ID
    sample_id = paste0(core_id, "_", sprintf("%03d", row_number()))
  )

n_pass <- sum(cores_complete$qa_pass)
n_fail <- sum(!cores_complete$qa_pass)

log_message(sprintf("Final QA: %d samples passed, %d failed", n_pass, n_fail))

# ============================================================================
# CHECK FOR SPATIAL BOUNDARY FILES (AOI + STRATA)
# ============================================================================
# These files are OPTIONAL but unlock zonal statistics in Module 04 (kriging)
# and Module 05 (RF model), enabling per-AOI carbon stock summaries.
# ============================================================================

cat("\n========================================\n")
cat("CHECKING FOR SPATIAL BOUNDARY FILES\n")
cat("========================================\n\n")

# ── AOI boundary ────────────────────────────────────────────────────────────
aoi_configured <- exists("AOI_FILE") && !is.null(AOI_FILE)
aoi_found      <- aoi_configured && file.exists(AOI_FILE)

if (aoi_found) {
  cat(sprintf("  ✓ AOI boundary file found:\n      %s\n", AOI_FILE))
  # Try to read it and report geometry type
  tryCatch({
    aoi_sf <- sf::st_read(AOI_FILE, quiet = TRUE)
    cat(sprintf("      Features: %d  |  CRS: %s\n",
                nrow(aoi_sf),
                sf::st_crs(aoi_sf)$input))
    cat("    → Zonal statistics will be available in Module 04 (kriging)\n\n")
  }, error = function(e) {
    cat(sprintf("    ⚠ File found but could not be read: %s\n\n", e$message))
  })
} else if (aoi_configured && !aoi_found) {
  cat(sprintf("  ✗ AOI_FILE is set in config but file NOT found:\n      %s\n", AOI_FILE))
  cat("    → Check the path in blue_carbon_config.R\n\n")
} else {
  cat("  ⚠ No AOI boundary file configured (AOI_FILE = NULL).\n")
  cat("    To enable zonal statistics in Module 04 and 05:\n")
  cat(sprintf("    1. Export your study-area boundary as a shapefile, GeoJSON, or GPKG\n"))
  cat(sprintf("    2. Place it in: %s\n", DATA_RAW_DIR))
  cat( "    3. Update AOI_FILE in Analysis_Workflow/blue_carbon_config.R:\n")
  cat(sprintf("       AOI_FILE <- file.path(DATA_RAW_DIR, 'aoi_boundary.shp')\n\n"))
}

# ── GEE stratum masks ────────────────────────────────────────────────────────
strata_dir_used <- if (exists("STRATA_DIR")) STRATA_DIR else file.path(DATA_RAW_DIR, "gee_strata")
strata_dir_exists <- dir.exists(strata_dir_used)

strata_tifs <- if (strata_dir_exists) {
  list.files(strata_dir_used, pattern = "\\.tif$", full.names = FALSE)
} else {
  character(0)
}

if (length(strata_tifs) > 0) {
  cat(sprintf("  ✓ Found %d stratum mask TIF(s) in:\n      %s\n",
              length(strata_tifs), strata_dir_used))
  for (f in strata_tifs) cat(sprintf("      - %s\n", f))
  cat("    → Per-stratum spatial analysis enabled in Module 04\n\n")
} else {
  cat(sprintf("  ⚠ No GEE stratum mask TIFs found in:\n      %s\n", strata_dir_used))
  cat("    To enable per-stratum spatial mapping in Module 04:\n")
  cat("    1. Run the GEE covariate script to export individual stratum masks\n")
  cat(sprintf("    2. Place each stratum TIF in: %s\n", strata_dir_used))
  cat("    3. File naming (lowercase, underscores):\n")
  for (s in VALID_STRATA) {
    cat(sprintf("       '%s'  →  %s.tif\n", s,
                tolower(gsub("[^A-Za-z0-9]", "_", gsub("\\s+", "_", s)))))
  }
  cat("    Module 04 will proceed using point-based kriging only.\n\n")
}

# ── Summarise spatial readiness ─────────────────────────────────────────────
cat("  Spatial analysis readiness:\n")
cat(sprintf("    AOI boundary:   %s\n",
            ifelse(aoi_found, sprintf("✓ Ready (%s)", basename(AOI_FILE)), "⚠ Not configured")))
cat(sprintf("    Stratum masks:  %s\n\n",
            ifelse(length(strata_tifs) > 0,
                   sprintf("✓ Ready (%d mask(s))", length(strata_tifs)),
                   "⚠ Not found")))

log_message(sprintf("Spatial files check — AOI: %s | Strata masks: %d",
                    ifelse(aoi_found, "found", "not configured"),
                    length(strata_tifs)))

# ============================================================================
# EXPORT CLEANED DATA
# ============================================================================

log_message("Exporting cleaned data...")

# Save as RDS (preserves data types)
saveRDS(cores_complete, "data_processed/cores_clean_bluecarbon.rds")
log_message("Saved: cores_clean_bluecarbon.rds")

# Save as CSV (portable)
write_csv(cores_complete, "data_processed/cores_clean_bluecarbon.csv")
log_message("Saved: cores_clean_bluecarbon.csv")

# Save core totals
saveRDS(core_totals, "data_processed/core_totals.rds")
write_csv(core_totals, "data_processed/core_totals.csv")
log_message("Saved: core_totals")

# Save carbon by stratum summary
write_csv(carbon_by_stratum, "data_processed/carbon_by_stratum_summary.csv")
log_message("Saved: carbon_by_stratum_summary.csv")

# ============================================================================
# GENERATE QA REPORT
# ============================================================================

log_message("Generating QA report...")

qa_report <- list(
  # Overall statistics
  total_cores = n_distinct(cores_complete$core_id),
  total_samples = nrow(cores_complete),
  samples_passed_qa = n_pass,
  samples_failed_qa = n_fail,

  # By stratum
  cores_by_stratum = stratum_stats,
  carbon_by_stratum = carbon_by_stratum,

  # VM0033 compliance (NEW)
  vm0033_compliance = vm0033_compliance,
  n_compliant_strata = sum(vm0033_compliance$vm0033_compliant),
  n_total_strata = nrow(vm0033_compliance),

  # Bulk density
  bd_measured = n_bd_measured,
  bd_estimated = n_bd_missing,
  bd_transparency = bd_transparency,

  # Depth profile completeness (NEW)
  depth_completeness_summary = depth_completeness_summary,

  # QA flags
  qa_issues = list(
    invalid_depths = n_depth_invalid,
    invalid_soc = n_soc_invalid,
    invalid_bd = n_bd_invalid,
    cores_no_location = length(cores_no_location),
    cores_no_samples = length(cores_no_samples)
  ),

  # Metadata
  processing_date = Sys.Date(),
  project_name = PROJECT_NAME,
  scenario_type = PROJECT_SCENARIO,
  monitoring_year = MONITORING_YEAR
)

saveRDS(qa_report, "diagnostics/qaqc/qa_report.rds")
log_message("Saved: qa_report.rds")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("MODULE 01 COMPLETE\n")
cat("========================================\n\n")

cat("Data Summary:\n")
cat("----------------------------------------\n")
cat(sprintf("Cores processed: %d\n", n_distinct(cores_complete$core_id)))
cat(sprintf("Samples processed: %d\n", nrow(cores_complete)))
cat(sprintf("QA pass rate: %.1f%%\n", 100 * n_pass / nrow(cores_complete)))
cat(sprintf("\nStrata represented: %d\n", n_distinct(cores_complete$stratum)))

cat("\nSamples by stratum:\n")
for (i in 1:nrow(stratum_stats)) {
  cat(sprintf("  %s: %d cores, %d samples\n", 
              stratum_stats$stratum[i],
              stratum_stats$n_cores[i],
              stratum_stats$n_samples[i]))
}

cat("\nBulk density:\n")
cat(sprintf("  Measured: %d samples (%.1f%%)\n", n_bd_measured,
            100 * n_bd_measured / nrow(cores_complete)))
cat(sprintf("  Estimated: %d samples (%.1f%%)\n", n_bd_missing,
            100 * n_bd_missing / nrow(cores_complete)))

cat(sprintf("\n%s:\n",
            vm_label("Statistical Power & VM0033 Assessment",
                     "Statistical Adequacy Assessment")))
cat(sprintf("  Adequate strata: %d/%d\n",
            sum(vm0033_compliance$vm0033_compliant),
            nrow(vm0033_compliance)))
if (sum(vm0033_compliance$vm0033_compliant) < nrow(vm0033_compliance)) {
  cat(sprintf("  ⚠ Review %s for details\n",
              vm_label("vm0033_compliance_report.csv", "sampling_adequacy_report.csv")))
}

cat("\nOutputs saved to data_processed/:\n")
cat("  Core Data:\n")
cat("    - cores_clean_bluecarbon.rds\n")
cat("    - cores_clean_bluecarbon.csv\n")
cat("    - core_totals.csv\n")
cat("  Summaries:\n")
cat("    - cores_summary_by_stratum.csv\n")
cat("    - carbon_by_stratum_summary.csv\n")
cat(sprintf("  %s:\n",
            vm_label("Statistical Power & VM0033 Assessment",
                     "Statistical Adequacy Assessment")))
cat(sprintf("    - %s\n",
            vm_label("vm0033_compliance_report.csv", "sampling_adequacy_report.csv")))
cat("  Bulk Density:\n")
cat("    - bd_transparency_report.csv (NEW)\n")
cat("  Depth Profiles:\n")
cat("    - core_depth_completeness.csv (NEW)\n")
cat("    - depth_completeness_summary.csv (NEW)\n")
cat("  Core Type Comparison:\n")
cat("    - core_type_summary.csv (if applicable)\n")
cat("    - core_type_statistical_tests.csv (if applicable)\n")
cat("  QA Report:\n")
cat("    - qa_report.rds\n")

cat("\nNext steps:\n")
cat(sprintf("  1. Review the %s\n",
            vm_label("statistical power and VM0033 assessment report",
                     "statistical adequacy assessment report")))
cat("  2. Check BD transparency and depth completeness\n")
cat("  3. If needed, collect additional samples for strata below the adequacy threshold\n")
cat("  4. Run: source('02_exploratory_analysis_bluecarbon.R')\n\n")

# ============================================================================
# QUICK-ACCESS COPY: BASIC ANALYSIS OUTPUTS
# ============================================================================

log_message("Copying key Module 01 outputs to outputs/Basic_analysis/Step1_NonSpatial/Tables...")
step1_tables_dir <- "outputs/Basic_analysis/Step1_NonSpatial/Tables"
dir.create(step1_tables_dir, recursive = TRUE, showWarnings = FALSE)

basic_copy_map <- c(
  "data_processed/cores_clean_bluecarbon.csv" = "Cleaned_Core_Locations_and_Samples.csv",
  "data_processed/core_totals.csv" = "Aggregated_Core_Totals_by_Core.csv",
  "data_processed/cores_summary_by_stratum.csv" = "Core_Summary_Table_by_Stratum.csv",
  "data_processed/carbon_by_stratum_summary.csv" = "Carbon_Stock_Summary_by_Stratum.csv"
)
# Add adequacy/compliance report under the correct source filename
adequacy_src <- file.path("diagnostics/data_prep",
                          vm_label("vm0033_compliance_report.csv", "sampling_adequacy_report.csv"))
adequacy_dst <- vm_label("Sampling_Power_and_VM0033_Compliance_Assessment.csv",
                         "Sampling_Adequacy_Assessment.csv")
basic_copy_map[adequacy_src] <- adequacy_dst

for (src in names(basic_copy_map)) {
  dst <- file.path(step1_tables_dir, basic_copy_map[[src]])
  if (file.exists(src)) {
    file.copy(src, dst, overwrite = TRUE)
    log_message(sprintf("Copied to Step1_NonSpatial/Tables: %s", basename(dst)))
  } else {
    log_message(sprintf("Skipped missing file: %s", src), "WARNING")
  }
}

log_message("=== MODULE 01 COMPLETE ===")
