# ============================================================================
# BLUE CARBON PROJECT CONFIGURATION
# ============================================================================
# Edit these parameters for your specific project
# This file is sourced by analysis modules
#
# IMPORTANT: All path variables below are relative to BlueCarbon_Workflow_V1.0/
# (the project root). Each analysis script sets setwd() to that root via its
# PATH RESOLVER block before sourcing this config.
# ============================================================================

# ============================================================================
# PROJECT DATA PATHS
# ============================================================================

PRE_ANALYSIS_DIR <- "Pre-Analysis Data Preparation"

# Input data
DATA_RAW_DIR    <- file.path(PRE_ANALYSIS_DIR, "data_raw")     # field core CSVs
DATA_GLOBAL_DIR <- file.path(PRE_ANALYSIS_DIR, "data_global")  # global DB + GEE covariate CSVs
COVARIATES_DIR  <- file.path(PRE_ANALYSIS_DIR, "covariates")   # GEE 27-band multi-band TIF

# Core input filenames (inside DATA_RAW_DIR)
CORE_LOCATIONS_FILE <- file.path(DATA_RAW_DIR, "core_locations.csv")
CORE_SAMPLES_FILE   <- file.path(DATA_RAW_DIR, "core_samples.csv")

# GEE covariate raster — single multi-band TIF (27 canonical bands)
# Set COVARIATE_RASTER to the exact filename exported from GEE.
COVARIATE_RASTER <- file.path(COVARIATES_DIR, "BlueCarbon_Covariate_Snapshot_25m_2020_2023.tif")

# AOI boundary file (shapefile, GeoJSON, or GPKG) — used for zonal statistics.
# Set to NULL if you do not have a boundary file yet.
# Module 04 (kriging) and Module 05 (RF) will warn but proceed without it.
# Example: AOI_FILE <- file.path(DATA_RAW_DIR, "aoi_boundary.shp")
AOI_FILE <- file.path(DATA_RAW_DIR, "aoi_boundary.geojson")

# GEE stratum mask rasters directory — individual TIF per stratum (optional).
# Enables per-stratum zonal statistics and spatial validation in Module 04.
# File naming convention: "Upper Marsh" → upper_marsh.tif
# Export from the GEE covariate script and place TIFs in this folder.
STRATA_DIR <- file.path(DATA_RAW_DIR, "gee_strata")

# Column name in AOI_FILE that identifies ecosystem strata.
# The values in this column must match entries in VALID_STRATA (e.g. "IM", "NM", "MF").
# Set to NULL if your AOI has only one zone (whole-site analysis without stratum breakdown).
# Example: if your GeoJSON has a property called "ecosystem_type", set to "ecosystem_type".
AOI_STRATUM_FIELD <- NULL #"stratum"

# Confidence level for uncertainty intervals in Stage 2a Simple Extrapolation.
# Default 0.95 = 95% confidence interval (±1.96 × SE). Applies to the
# lower/upper bound columns in simple_extrapolation_by_stratum.csv outputs.
EXTRAP_CONFIDENCE <- 0.95

# ============================================================================
# PROJECT METADATA (VM0033 Required)
# ============================================================================

PROJECT_NAME <- "BC_Coastal_BlueCarbon_2026_Example"
PROJECT_SCENARIO <- "PROJECT"  # Options: BASELINE, PROJECT, CONTROL, DEGRADED
MONITORING_YEAR <- 2026

# Project location (for documentation)
PROJECT_LOCATION <- "Chemainus Test 1"
PROJECT_DESCRIPTION <- "Blue carbon calculation for coastal salt marsh and eelgrass and spatial mapping"

# ============================================================================
# ECOSYSTEM STRATIFICATION
# ============================================================================

# Valid ecosystem strata (must match GEE stratification tool)
#
# FILE NAMING CONVENTION:
#   Module 05 auto-detects GEE stratum masks using this pattern:
#   "Stratum Name" → stratum_name.tif in data_raw/gee_strata/
#
# Examples:
#   "Upper Marsh"           → upper_marsh.tif
#   "Underwater Vegetation" → underwater_vegetation.tif
#   "Emerging Marsh"        → emerging_marsh.tif
#
# CUSTOMIZATION OPTIONS:
#   1. Simple: Edit VALID_STRATA below and export GEE masks with matching names
#   2. Advanced: Create stratum_definitions.csv in project root for custom file names
#      and optional metadata (see stratum_definitions_EXAMPLE.csv template)
#
# See README section "Customizing Ecosystem Strata" for full details.
#
VALID_STRATA <- c(
  "IM",         
  "NM",             
  "MF"
)

# Stratum colors for plotting (match GEE tool)
STRATUM_COLORS <- c(
  "IM"              = "#FFFF99",
  "NM"            = "#99FF99",
  "MF"            = "#33CC33"
)

# ============================================================================
# DEPTH CONFIGURATION
# ============================================================================

# VM0033 standard depth intervals (cm) - depth midpoints for harmonization
# These correspond to VM0033 depth layers: 0-15, 15-30, 30-50, 50-100 cm
VM0033_DEPTH_MIDPOINTS <- c(7.5, 22.5, 40, 75)

# VM0033 depth intervals (cm) - for mass-weighted aggregation
VM0033_DEPTH_INTERVALS <- data.frame(
  depth_top = c(0, 15, 30, 50),
  depth_bottom = c(15, 30, 50, 100),
  depth_midpoint = c(7.5, 22.5, 40, 75),
  thickness_cm = c(15, 15, 20, 50)
)

# Standard depths for harmonization (VM0033 midpoints are default)
STANDARD_DEPTHS <- VM0033_DEPTH_MIDPOINTS

# Fine-scale depth intervals (optional, for detailed analysis)
FINE_SCALE_DEPTHS <- c(0, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100)

# Maximum core depth (cm)
MAX_CORE_DEPTH <- 100

# Key depth intervals for reporting (cm)
REPORTING_DEPTHS <- list(
  surface = c(0, 30),      # Top 30 cm (most active layer)
  subsurface = c(30, 100)  # 30-100 cm (long-term storage)
)

# ============================================================================
# COORDINATE SYSTEMS
# ============================================================================

# Input CRS (usually WGS84 for GPS data)
INPUT_CRS <- 4326  # EPSG:4326 (WGS84)

# Processing CRS (projected, equal-area for accurate calculations)
# Change this for your region:
PROCESSING_CRS <- 3347 # Canada Albers Equal Area (good for all Canada)
# Other options:
#   - 3005:  EPSG:3005 (NAD83 / BC Albers) - OPTIMIZED FOR BC
#   - 3347: Canada Albers Equal Area (good for all Canada)
#   - 32610: WGS 84 / UTM zone 10N (BC coast - Chemainus area)
#   - 32611: WGS 84 / UTM zone 11N (BC interior)

# ============================================================================
# BULK DENSITY DEFAULTS BY STRATUM
# ============================================================================
# Use these when bulk density is not measured
# Values in g/cm³ based on literature for BC coastal ecosystems
BD_DEFAULTS <- list(
  "IM" = 0.8,              # Lower density, more organic matter
  "NM" = 0.8,                # Moderate density
  "MF" = 0.8              # Higher 
)

# ============================================================================
# QUALITY CONTROL THRESHOLDS
# ============================================================================

# Soil Organic Carbon (SOC) thresholds (g/kg)
QC_SOC_MIN <- 0      # Minimum valid SOC
QC_SOC_MAX <- 500    # Maximum valid SOC (adjust for your ecosystem)

# Bulk Density thresholds (g/cm³)
QC_BD_MIN <- 0.1     # Minimum valid bulk density
QC_BD_MAX <- 3.0     # Maximum valid bulk density

# Depth thresholds (cm)
QC_DEPTH_MIN <- 0
QC_DEPTH_MAX <- MAX_CORE_DEPTH

# Coordinate validity (decimal degrees for WGS84)
QC_LON_MIN <- -180
QC_LON_MAX <- 180
QC_LAT_MIN <- -90
QC_LAT_MAX <- 90

# ============================================================================
# SAMPLING REQUIREMENTS
# (Thresholds follow both ecological best practice and VM0033 when applicable)
# ============================================================================

# Minimum cores per stratum — ecological best practice (also VM0033 minimum)
VM0033_MIN_CORES <- 3

# Academic/biometric baseline targets for statistically robust sampling
ACADEMIC_CONFIDENCE_LEVEL <- 0.90
ACADEMIC_MARGIN_OF_ERROR <- 0.20
ACADEMIC_DEFAULT_CV <- 0.5

# Target precision: 20% relative error at 95% CI
# (standard ecological target; also falls within VM0033 acceptable range of 10-20%)
VM0033_TARGET_PRECISION <- 20  # percent

# Target CV threshold (higher CV = higher uncertainty)
VM0033_CV_THRESHOLD <- 30  # percent

# Assumed CV for sample size calculation (conservative estimate)
VM0033_ASSUMED_CV <- 30  # percent

# ============================================================================
# TEMPORAL MONITORING & ADDITIONALITY PARAMETERS
# ============================================================================

# Valid scenario types for VM0033 (expanded for restoration chronosequence)
# Core scenarios:
# - BASELINE: Pre-restoration or current degraded condition (t0)
# - DEGRADED: Heavily degraded/lost ecosystem (lower bound)
# - DISTURBED: Moderately impacted ecosystem
# - REFERENCE: Natural healthy ecosystem (upper bound target)
# - CONTROL: No-intervention control site (tracks natural variation)
# Restoration trajectory scenarios:
# - PROJECT_Y0: Immediately post-restoration
# - PROJECT_Y1: 1 year post-restoration
# - PROJECT_Y5: 5 years post-restoration (VM0033 first verification)
# - PROJECT_Y10: 10 years post-restoration (VM0033 second verification)
# - PROJECT_Y15: 15+ years post-restoration
# - PROJECT: Generic project scenario (when year not specified)
VALID_SCENARIOS <- c("BASELINE", "DEGRADED", "DISTURBED", "REFERENCE", "CONTROL",
                     "PROJECT", "PROJECT_Y0", "PROJECT_Y1", "PROJECT_Y5",
                     "PROJECT_Y10", "PROJECT_Y15", "CUSTOM")

# Scenario hierarchy for modeling (relative carbon stock levels)
# Used by Module 08A to model missing scenarios from available data
SCENARIO_CARBON_LEVELS <- c(
  DEGRADED = 1.0,
  DISTURBED = 2.0,
  BASELINE = 3.0,
  PROJECT_Y0 = 3.0,
  PROJECT_Y1 = 4.0,
  PROJECT_Y5 = 6.0,
  PROJECT_Y10 = 8.0,
  PROJECT_Y15 = 9.5,
  REFERENCE = 10.0
)

# Minimum monitoring frequency (years) - VM0033 typically requires verification every 5 years
VM0033_MONITORING_FREQUENCY <- 5

# Minimum years for temporal change analysis
MIN_YEARS_FOR_CHANGE <- 3  # At least 3 years to establish trend

# Additionality test confidence level
ADDITIONALITY_CONFIDENCE <- 0.95  # 95% CI for statistical tests

# Conservative approach for additionality calculations
ADDITIONALITY_METHOD <- "lower_bound"  # Options: "mean", "lower_bound", "conservative"
# - "mean": Use mean difference between project and baseline
# - "lower_bound": Use 95% CI lower bound of difference (most conservative approach)
# - "conservative": Use mean - 1SD (moderately conservative)

# ============================================================================
# SCENARIO MODELING PARAMETERS (Module 08A)
# ============================================================================

# Enable scenario modeling (generate synthetic scenarios from reference trajectories)
SCENARIO_MODELING_ENABLED <- TRUE

# Canadian literature database for BC Coast ecosystems
CANADIAN_LITERATURE_DB <- "canadian_bluecarbon_parameters.csv"

# Scenario modeling configuration file
SCENARIO_CONFIG_FILE <- "scenario_modeling_config.csv"

# Recovery model types for reference trajectory method
# - "exponential": Fast initial recovery, slowing over time (most common)
# - "linear": Constant accumulation rate
# - "logistic": S-shaped curve with inflection point
# - "asymptotic": Approaches target asymptotically
RECOVERY_MODEL_TYPE <- "exponential"

# Uncertainty inflation for modeled scenarios (%)
# Adds additional uncertainty to account for modeling assumptions
MODELING_UNCERTAINTY_BUFFER <- 10  # percent

# Spatial resolution for modeled scenario rasters (if generating spatial outputs)
MODELED_RASTER_RESOLUTION <- 30  # meters

# ============================================================================
# WORKFLOW TOGGLES
# ============================================================================

# Enable Bayesian workflow (requires GEE prior maps)
USE_BAYESIAN <- FALSE  # Set to TRUE to enable Part 4

# Enable Transfer Learning workflow
# Prerequisites:
#   1. Run P4_3b_Depth_Harmonization_Global.R  →  data_processed/global_cores_harmonized_VM0033.csv
#   2. Run S1_01_data_prep.R                   →  data_processed/cores_harmonized_bluecarbon.csv
#   3. GEE covariate raster at COVARIATE_RASTER (for Model 1+3 spatial predictions)
#   4. Embedding raster at EMBEDDING_RASTER    (for Model 2+3 spatial predictions)
USE_TRANSFER_LEARNING <- TRUE  # Set to TRUE to enable Transfer Learning

# Embedding similarity GeoTIFF (5-band raster; band names = embedding_sim_to_cluster_0..4)
# Generated externally and placed in the covariates folder.
# Set to NULL if not yet available — script will fall back to nearest-neighbour assignment.
EMBEDDING_RASTER <- file.path(COVARIATES_DIR, "BlueCarbon_Embedding_Similarity_AOI.tif")

# ── Stage execution gates ─────────────────────────────────────────────────
# Control which workflow stages run from MasterScript_GitHub.R.
# Stage 1 (non-spatial) always runs — these gates only affect spatial stages.
#
# RUN_STAGE2 = TRUE:  Run Stage 2 (Simple Extrapolation + Kriging).
#   Requires: AOI_FILE set to a valid GeoJSON/shapefile path.
# RUN_STAGE3 = TRUE:  Run Stage 3 (Random Forest + Carbon Totals + Reporting).
#   Requires: AOI_FILE AND at least one .tif in COVARIATES_DIR.
#
# Set to FALSE to skip a stage even when its inputs are available.
RUN_STAGE2 <- TRUE
RUN_STAGE3 <- TRUE

# ── VM0033 Toggle ─────────────────────────────────
# VM0033_COMPLIANCE = FALSE (default): Uses neutral ecological/scientific
#   language throughout. Suitable when working with Indigenous Nations,
#   community partners, or any context where carbon market terminology
#   may be distrustful or premature.
#
# VM0033_COMPLIANCE = TRUE: Uses full VM0033/MMRV/carbon-credit framing.
#   Enables verification package outputs, ORRAA compliance checks, and
#   Verra-specific report titles. Use when preparing a formal carbon
#   credit submission under VM0033 (Tidal Wetland & Seagrass Restoration).
#
# NOTE: The underlying science, depth intervals, and statistical thresholds
#   are IDENTICAL in both modes — only the language and output file names
#   change.
VM0033_COMPLIANCE <- FALSE  # Set to TRUE for VM0033/MMRV carbon credit mode

# ── Language helpers (used by all analysis scripts) ────────────────────────
# Returns vm0033_text if VM0033_COMPLIANCE is TRUE, otherwise ecological_text.
vm_label <- function(vm0033_text, ecological_text) {
  if (isTRUE(VM0033_COMPLIANCE)) vm0033_text else ecological_text
}

# Returns the correct report output directory based on the compliance toggle.
#   VM0033 mode:    "outputs/mmrv_reports"
#   Ecological mode: "outputs/reports"
report_dir <- function() {
  if (isTRUE(VM0033_COMPLIANCE)) "outputs/mmrv_reports" else "outputs/reports"
}

# Returns the Stage 3 presentation report directory (always the same path,
# used for placing report copies in the stage-structured output folder).
stage3_report_dir <- function() "outputs/Advanced_analysis/Step3_Spatial_Advanced/Reports"
# ───────────────────────────────────────────────────────────────────────────

# Prior data directory
BAYESIAN_PRIOR_DIR <- "data_prior"

# GEE Data Sources (BC Coast)
# SoilGrids 250m - Global soil organic carbon maps
GEE_SOILGRIDS_ASSET <- "projects/soilgrids-isric/soc_mean"
GEE_SOILGRIDS_UNCERTAINTY <- "projects/soilgrids-isric/soc_uncertainty"

# Sothe et al. 2022 - BC Forest biomass and soil carbon
# Users should input their specific GEE asset paths here:
GEE_SOTHE_FOREST_BIOMASS <- ""  # User to provide
GEE_SOTHE_SOIL_CARBON <- ""     # User to provide
GEE_SOTHE_OTHER_BIOMASS <- ""   # User to provide

# Prior resolution (will be resampled to PREDICTION_RESOLUTION)
PRIOR_RESOLUTION <- 250  # meters (SoilGrids native resolution)
PREDICTION_RESOLUTION <- 30
# Bayesian sampling design (Neyman allocation)
USE_NEYMAN_SAMPLING <- TRUE  # Enable optimal allocation based on prior uncertainty
NEYMAN_STRATA <- 3           # Number of uncertainty strata (low/med/high)
NEYMAN_BUFFER_SAMPLES <- 1.2 # Oversample by 20% to account for inaccessible locations

# Uncertainty strata thresholds (coefficient of variation %)
UNCERTAINTY_LOW_THRESHOLD <- 10    # CV < 10% = low uncertainty
UNCERTAINTY_HIGH_THRESHOLD <- 30   # CV > 30% = high uncertainty
# Medium uncertainty = between thresholds

# Bayesian posterior weighting
# How to weight prior vs field data based on sample size
BAYESIAN_WEIGHT_METHOD <- "sqrt_samples"  # Options: "sqrt_samples", "linear", "fixed"
BAYESIAN_FIXED_WEIGHT <- 0.5              # Only used if method = "fixed"
BAYESIAN_TARGET_SAMPLES <- 30             # Target sample size for full field weight

# Precision adjustment
# Inflate prior uncertainty to account for potential bias/mismatch
PRIOR_UNCERTAINTY_INFLATION <- 1.2  # Multiply prior SE by this factor (conservative)

# Information gain threshold
# Minimum posterior uncertainty reduction to declare "informative prior"
MIN_INFORMATION_GAIN_PCT <- 20  # At least 20% uncertainty reduction

# ============================================================================
# DEPTH HARMONIZATION PARAMETERS
# ============================================================================

# Interpolation method: "equal_area_spline", "smoothing_spline", "linear", "all"
INTERPOLATION_METHOD <- "equal_area_spline"  # recommended for tidal wetland depth harmonization

# Spline smoothing parameters by core type
SPLINE_SPAR_HR <- 0.3           # Less smoothing for high-resolution cores
SPLINE_SPAR_COMPOSITE <- 0.5    # More smoothing for composite cores
SPLINE_SPAR_AUTO <- NULL        # NULL = automatic cross-validation

# Monotonicity parameters
ALLOW_DEPTH_INCREASES <- FALSE   # Allow slight SOC increases with depth (common in some ecosystems)
MAX_INCREASE_THRESHOLD <- 20    # Maximum % increase allowed between adjacent depths

# ============================================================================
# UNCERTAINTY PARAMETERS
# ============================================================================

# Confidence level for uncertainty estimation (95% — standard ecological practice)
CONFIDENCE_LEVEL <- 0.95

# Bootstrap parameters for spline uncertainty
BOOTSTRAP_ITERATIONS <- 100
BOOTSTRAP_SEED <- 42

# Cross-validation parameters
CV_FOLDS <- 3           # Number of folds for spatial CV (reduced for small datasets)
CV_SEED <- 42           # Random seed for reproducibility

# ============================================================================
# SPATIAL MODELING PARAMETERS
# ============================================================================

# Prediction resolution (meters)
KRIGING_CELL_SIZE <- 25
RF_CELL_SIZE <- 25

# Kriging parameters
KRIGING_MAX_DISTANCE <- 5000  # Maximum distance for variogram (meters)
KRIGING_CUTOFF <- NULL        # NULL = automatic
KRIGING_WIDTH <- 100          # Lag width for variogram (meters)

# Random Forest parameters
RF_NTREE <- 500              # Number of trees
RF_MTRY <- NULL              # NULL = automatic (sqrt of predictors)
RF_MIN_NODE_SIZE <- 5        # Minimum node size
RF_IMPORTANCE <- "permutation"  # Variable importance method

# ============================================================================
# AREA OF APPLICABILITY (AOA) PARAMETERS
# ============================================================================

# Enable AOA analysis (requires CAST package)
ENABLE_AOA <- FALSE

# AOA threshold (dissimilarity index)
AOA_THRESHOLD <- "default"  # "default" or numeric value

# ============================================================================
# REPORT GENERATION PARAMETERS
# ============================================================================

# Figure dimensions for saving (inches)
FIGURE_WIDTH <- 10
FIGURE_HEIGHT <- 6
FIGURE_DPI <- 300

# Table formatting
TABLE_DIGITS <- 2  # Decimal places for tables

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

#' Convert between carbon stock units
#'
#' @param value Numeric value to convert
#' @param from Source unit (kg_m2, Mg_ha, g_kg, pct)
#' @param to Target unit
#' @return Converted value
#' @examples
#' convert_units(1, "kg_m2", "Mg_ha")  # Returns 10
#' convert_units(10, "Mg_ha", "kg_m2") # Returns 1
convert_units <- function(value, from, to) {
  conversions <- list(
    "kg_m2_to_Mg_ha" = 10,
    "Mg_ha_to_kg_m2" = 0.1,
    "g_kg_to_pct" = 0.1,
    "pct_to_g_kg" = 10
  )

  key <- paste(from, "to", to, sep = "_")
  if (key %in% names(conversions)) {
    return(value * conversions[[key]])
  } else {
    stop(sprintf("Unknown conversion: %s to %s", from, to))
  }
}

# ============================================================================
# SESSION TRACKING
# ============================================================================

# Session tracking for reproducibility and unique output naming
SESSION_START <- Sys.time()
SESSION_ID <- format(SESSION_START, "%Y%m%d_%H%M%S")

# ============================================================================
# END OF CONFIGURATION
# ============================================================================

# Print confirmation when loaded
if (interactive()) {
  cat("Blue Carbon configuration loaded ✓
")
  cat(sprintf("  Project: %s
", PROJECT_NAME))
  cat(sprintf("  Location: %s
", PROJECT_LOCATION))
  cat(sprintf("  Scenario: %s
", PROJECT_SCENARIO))
  cat(sprintf("  Monitoring year: %d
", MONITORING_YEAR))
  cat(sprintf("  Session ID: %s
", SESSION_ID))
}
