## ============================================================================
## BLUE CARBON MMRV WORKFLOW: REMOTE IMPORT & STEP-BY-STEP EXECUTION
## ============================================================================
## LOCATION: Analysis_Workflow/MasterScript_GitHub.R
##
## HOW TO RUN:
##   Set your working directory to BlueCarbon_Workflow_V1.0/ before sourcing,
##   or use RStudio's "Source" button with this file open.
##
##   In R: setwd("/path/to/BlueCarbon_Workflow_V1.0")
##         source("Analysis_Workflow/MasterScript_GitHub.R")
##
## WORKFLOW STRUCTURE:
## ─────────────────────────────────────────────────────────────────────────
##  STEP 1 — Non-Spatial Analysis       (ALWAYS runs)
##    Data prep → EDA → Depth harmonization
##    Outputs: outputs/Basic_analysis/Step1_NonSpatial/
##    Requires: core_locations.csv + core_samples.csv
##
##  STEP 2 — Basic Spatial Analysis     (requires AOI_FILE in config)
##    2a. Simple Extrapolation (mean × area)
##    2b. Kriging (geostatistical prediction)
##    Outputs: outputs/Basic_analysis/Step2_Spatial_Basic/SimpleExtrapolation/ + outputs/Basic_analysis/Step2_Spatial_Basic/Kriging/
##    Requires: AOI GeoJSON/shapefile with stratum column (set AOI_FILE in config)
##
##  STEP 3 — Advanced Spatial Analysis  (requires AOI_FILE + covariates)
##    Random Forest → Carbon stock totals → HTML report
##    Outputs: outputs/Advanced_analysis/Step3_Spatial_Advanced/
##    Requires: AOI_FILE + GEE covariate TIF in covariates/
##
##  BAYESIAN (optional)  — set USE_BAYESIAN <- TRUE in blue_carbon_config.R
##  TRANSFER LEARNING (optional) — set USE_TRANSFER_LEARNING <- TRUE in config
## ─────────────────────────────────────────────────────────────────────────
##
## All scripts use paths relative to BlueCarbon_Workflow_V1.0/ (project root).
## ============================================================================

# ── Set working directory to project root (BlueCarbon_Workflow_V1.0/) ───────
rm(list = ls())

this_script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      dirname(normalizePath(rstudioapi::getActiveDocumentContext()$path))
    } else {
      stop("Could not detect script directory. Run setwd('/path/to/BlueCarbon_Workflow_V1.0') first.")
    }
  }
)
project_root <- dirname(this_script_dir)  # Analysis_Workflow/ → BlueCarbon_Workflow_V1.0/
if (dir.exists(project_root)) setwd(project_root)
cat(sprintf("Working directory: %s\n\n", getwd()))

# ============================================================================
# ── SETUP: REMOTE IMPORT FROM GITHUB ─────────────────────────────────────
# ============================================================================

repo_base <- "https://raw.githubusercontent.com/cathald/northstarproject_coastalbluecarbonmmrv/main/NorthStarProject_CoastalBlueCarbonMMRV/"

fetch_file <- function(file_path, overwrite = FALSE) {
  # SAFETY: Never overwrite an existing local file unless explicitly requested.
  # This prevents a missing/invalid GitHub URL from wiping locally-developed scripts.
  if (file.exists(file_path) && !overwrite) {
    # cat(sprintf("  [skip] Already exists locally: %s\n", file_path))
    return(invisible(NULL))
  }
  dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)
  tryCatch(
    download.file(paste0(repo_base, file_path), destfile = file_path, mode = "wb", quiet = TRUE),
    error = function(e) warning(sprintf("Could not download: %s\n  %s", file_path, e$message))
  )
}

# Core infrastructure
infra <- c(
  "Analysis_Workflow/blue_carbon_config.R",
  "Analysis_Workflow/Basic Analysis/Step1_NonSpatial/S1_00_setup_directories.R"
)
lapply(infra, fetch_file)

# Raw data templates
raw_files <- c(
  "Pre-Analysis Data Preparation/data_raw/core_locations.csv",
  "Pre-Analysis Data Preparation/data_raw/core_samples.csv",
  "Pre-Analysis Data Preparation/data_global/Global_Core_Locations.csv",
  "Pre-Analysis Data Preparation/data_global/Global_Core_Samples.csv",
  "Pre-Analysis Data Preparation/data_global/global_cores_with_gee_covariates.csv"
)
lapply(raw_files, fetch_file)

# Step 1 scripts
step1_scripts <- c(
  "Analysis_Workflow/Basic Analysis/Step1_NonSpatial/S1_01_data_prep.R",
  "Analysis_Workflow/Basic Analysis/Step1_NonSpatial/S1_02_exploratory_analysis.R",
  "Analysis_Workflow/Basic Analysis/Step1_NonSpatial/S1_03_depth_harmonization.R"
)

# Step 2 scripts
step2_scripts <- c(
  "Analysis_Workflow/Basic Analysis/Step2_Spatial_Basic/S2a_01_simple_extrapolation.R",
  "Analysis_Workflow/Basic Analysis/Step2_Spatial_Basic/S2b_02_kriging.R"
)

# Step 3 scripts
step3_scripts <- c(
  "Analysis_Workflow/Advanced Analysis/Step3_Spatial_Advanced/S3_01_random_forest.R",
  "Analysis_Workflow/Advanced Analysis/Step3_Spatial_Advanced/S3_02_carbon_stock_totals.R",
  "Analysis_Workflow/Advanced Analysis/Step3_Spatial_Advanced/S3_03_reporting.R"
)

lapply(c(step1_scripts, step2_scripts, step3_scripts), fetch_file)

cat("✓ SETUP COMPLETE: Scripts and raw data imported.\n\n")

# ============================================================================
# ── LOAD CONFIGURATION ───────────────────────────────────────────────────
# ============================================================================

source("Analysis_Workflow/blue_carbon_config.R")
source("Analysis_Workflow/Basic Analysis/Step1_NonSpatial/S1_00_setup_directories.R")

# ============================================================================
# ── STEP 1: NON-SPATIAL ANALYSIS (always runs) ───────────────────────────
# ── Requires: core_locations.csv + core_samples.csv
# ── Outputs:  outputs/Basic_analysis/Step1_NonSpatial/Tables/ and Plots/
# ============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  STEP 1 — NON-SPATIAL ANALYSIS                          ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

# S1_01: Data preparation, QA/QC, carbon stock calculations
# Inputs:  Pre-Analysis Data Preparation/data_raw/core_locations.csv + core_samples.csv
# Outputs: data_processed/cores_clean_bluecarbon.rds, carbon_by_stratum_summary.csv
source("Analysis_Workflow/Basic Analysis/Step1_NonSpatial/S1_01_data_prep.R")

# S1_02: Exploratory data analysis — depth profiles, SOC distributions, spatial map
# Inputs:  data_processed/cores_clean_bluecarbon.rds
# Outputs: outputs/plots/exploratory/*.png + outputs/Basic_analysis/Step1_NonSpatial/Plots/
source("Analysis_Workflow/Basic Analysis/Step1_NonSpatial/S1_02_exploratory_analysis.R")

# S1_03: Depth harmonization — standardize to VM0033 depths using equal-area splines
# Inputs:  data_processed/cores_clean_bluecarbon.rds
# Outputs: data_processed/cores_harmonized_bluecarbon.rds
source("Analysis_Workflow/Basic Analysis/Step1_NonSpatial/S1_03_depth_harmonization.R")

cat("\n✓ STEP 1 COMPLETE\n")
cat("  Non-spatial tables: outputs/Basic_analysis/Step1_NonSpatial/Tables/\n")
cat("  EDA plots:          outputs/Basic_analysis/Step1_NonSpatial/Plots/\n\n")

# ============================================================================
# ── STEP 2: BASIC SPATIAL ANALYSIS ───────────────────────────────────────
# ── Requires: AOI_FILE set in blue_carbon_config.R
# ──           AOI must have a stratum attribute column (AOI_STRATUM_FIELD)
# ── Outputs:  outputs/Basic_analysis/Step2_Spatial_Basic/SimpleExtrapolation/ + outputs/Basic_analysis/Step2_Spatial_Basic/Kriging/
# ============================================================================

aoi_available <- exists("AOI_FILE") && !is.null(AOI_FILE) && file.exists(AOI_FILE)

if (isTRUE(RUN_STAGE2) && aoi_available) {

  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  STEP 2 — BASIC SPATIAL ANALYSIS                        ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n\n")

  # S2a_01: Simple extrapolation — mean stock per stratum × stratum area
  # Inputs:  data_processed/carbon_by_stratum_summary.csv + AOI_FILE
  # Outputs: outputs/Basic_analysis/Step2_Spatial_Basic/SimpleExtrapolation/ (tables, maps, plots)
  cat("-- Step 2a: Simple Extrapolation --\n")
  source("Analysis_Workflow/Basic Analysis/Step2_Spatial_Basic/S2a_01_simple_extrapolation.R")

  # S2b_02: Kriging — variogram fitting and geostatistical spatial prediction
  # Inputs:  data_processed/cores_harmonized_bluecarbon.rds + AOI_FILE
  # Outputs: outputs/predictions/kriging/*.tif + outputs/Basic_analysis/Step2_Spatial_Basic/Kriging/
  cat("-- Step 2b: Kriging --\n")
  source("Analysis_Workflow/Basic Analysis/Step2_Spatial_Basic/S2b_02_kriging.R")

  cat("\n✓ STEP 2 COMPLETE\n")
  cat("  Simple extrapolation: outputs/Basic_analysis/Step2_Spatial_Basic/SimpleExtrapolation/\n")
  cat("  Kriging maps:         outputs/Basic_analysis/Step2_Spatial_Basic/Kriging/\n\n")

} else if (isTRUE(RUN_STAGE2)) {
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("  STEP 2 SKIPPED — AOI_FILE not set or file not found.\n\n")
  cat("  To run Step 2, set AOI_FILE in blue_carbon_config.R:\n")
  cat("    AOI_FILE <- 'Pre-Analysis Data Preparation/data_raw/aoi_boundary.geojson'\n")
  cat("  The AOI must contain a stratum attribute column.\n")
  cat("  Column name is specified by AOI_STRATUM_FIELD in config.\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")
} else {
  cat("  STEP 2 SKIPPED — RUN_STAGE2 = FALSE in config.\n\n")
}

# ============================================================================
# ── STEP 3: ADVANCED SPATIAL ANALYSIS ────────────────────────────────────
# ── Requires: AOI_FILE + at least one .tif in covariates/
# ── Outputs:  outputs/Advanced_analysis/Step3_Spatial_Advanced/
# ============================================================================

covariate_found <- length(list.files(COVARIATES_DIR, pattern = "\\.tif$")) > 0

if (isTRUE(RUN_STAGE3) && aoi_available && covariate_found) {

  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  STEP 3 — ADVANCED SPATIAL ANALYSIS (RANDOM FOREST)     ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n\n")

  # S3_01: Random Forest — ML prediction using GEE covariate TIFs
  # Inputs:  data_processed/cores_harmonized_bluecarbon.rds + covariates/*.tif
  # Outputs: outputs/predictions/rf/*.tif + outputs/Advanced_analysis/Step3_Spatial_Advanced/RF_Maps/*.png
  cat("-- Step 3a: Random Forest Modelling --\n")
  source("Analysis_Workflow/Advanced Analysis/Step3_Spatial_Advanced/S3_01_random_forest.R")

  # S3_02: Carbon stock totals — aggregate depth layers, calculate total kg per AOI
  # Inputs:  kriging + RF prediction TIFs + AOI_FILE
  # Outputs: outputs/carbon_stocks/*.csv (read by S3_03) + Stage3 totals table
  cat("-- Step 3b: Carbon Stock Totals --\n")
  source("Analysis_Workflow/Advanced Analysis/Step3_Spatial_Advanced/S3_02_carbon_stock_totals.R")

  # S3_03: Reporting — generate comprehensive HTML assessment/verification report
  # Inputs:  All outputs/carbon_stocks/*.csv + diagnostics + spatial outputs
  # Outputs: HTML report in report_dir() + outputs/Advanced_analysis/Step3_Spatial_Advanced/Reports/
  cat("-- Step 3c: Reporting --\n")
  source("Analysis_Workflow/Advanced Analysis/Step3_Spatial_Advanced/S3_03_reporting.R")

  cat("\n✓ STEP 3 COMPLETE\n")
  cat("  RF maps:         outputs/Advanced_analysis/Step3_Spatial_Advanced/RF_Maps/\n")
  cat("  Carbon totals:   outputs/Advanced_analysis/Step3_Spatial_Advanced/Carbon_Stock_Totals/\n")
  cat("  HTML report:     outputs/Advanced_analysis/Step3_Spatial_Advanced/Reports/\n\n")

} else if (isTRUE(RUN_STAGE3)) {
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("  STEP 3 SKIPPED — Missing inputs:\n\n")
  if (!aoi_available) {
    cat("  ✗ AOI_FILE not set → Set AOI_FILE in blue_carbon_config.R\n")
  }
  if (!covariate_found) {
    cat(sprintf("  ✗ No covariate TIFs found in: %s\n", COVARIATES_DIR))
    cat("    Export the 27-band covariate snapshot from Google Earth Engine\n")
    cat("    and place in the covariates/ directory.\n")
  }
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")
} else {
  cat("  STEP 3 SKIPPED — RUN_STAGE3 = FALSE in config.\n\n")
}

# ============================================================================
# ── BAYESIAN MODEL (optional — set USE_BAYESIAN = TRUE in config) ─────────
# ── Requires: GEE prior exports in data_prior/gee_exports/
# ============================================================================

if (exists("USE_BAYESIAN") && isTRUE(USE_BAYESIAN)) {
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  BAYESIAN MODEL                                          ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n\n")

  # Download Bayesian scripts from GitHub
  lapply(c(
    "Analysis_Workflow/Bayesian model/P3_0c_bayesian_prior_setup_bluecarbon.R",
    "Analysis_Workflow/Bayesian model/P3_1c_bayesian_sampling_design_bluecarbon.R",
    "Analysis_Workflow/Bayesian model/P3_06c_bayesian_posterior_estimation_bluecarbon.R",
    "Analysis_Workflow/Bayesian model/P3_07b_comprehensive_standards_report.R"
  ), fetch_file)

  # P3_0c: Extract global carbon stock priors from GEE rasters
  source("Analysis_Workflow/Bayesian model/P3_0c_bayesian_prior_setup_bluecarbon.R")

  # P3_1c: Neyman-allocated sampling design using prior uncertainty
  source("Analysis_Workflow/Bayesian model/P3_1c_bayesian_sampling_design_bluecarbon.R")

  # P3_06c: Bayesian posterior estimation (Normal-Normal conjugate)
  source("Analysis_Workflow/Bayesian model/P3_06c_bayesian_posterior_estimation_bluecarbon.R")

  # P3_07b: Comprehensive standards report combining all methods
  source("Analysis_Workflow/Bayesian model/P3_07b_comprehensive_standards_report.R")

} else {
  cat("  Bayesian model SKIPPED (USE_BAYESIAN = FALSE)\n")
  cat("  Set USE_BAYESIAN <- TRUE in blue_carbon_config.R to enable.\n\n")
}

# ============================================================================
# ── TRANSFER LEARNING (optional — set USE_TRANSFER_LEARNING = TRUE) ───────
# ── Requires: data_global/global_cores_with_gee_covariates.csv
# ============================================================================

if (exists("USE_TRANSFER_LEARNING") && isTRUE(USE_TRANSFER_LEARNING)) {
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  TRANSFER LEARNING                                       ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n\n")

  # Download Transfer Learning scripts from GitHub
  lapply(c(
    "Analysis_Workflow/Transfer Learning/P4_3b_Depth_Harmonization_Global.R",
    "Analysis_Workflow/Transfer Learning/P4_05_Transfer_Learning_Modelling.R",
    "Analysis_Workflow/Transfer Learning/P4_Transfer_Learning_Visualizations.R"
  ), fetch_file)

  # P4_3b: Harmonize global database profiles to VM0033 standard depths
  source("Analysis_Workflow/Transfer Learning/P4_3b_Depth_Harmonization_Global.R")

  # P4_05: Transfer learning modelling with Wadoux instance weighting
  source("Analysis_Workflow/Transfer Learning/P4_05_Transfer_Learning_Modelling.R")

  # P4: Transfer learning visualizations — compare local vs global predictions
  source("Analysis_Workflow/Transfer Learning/P4_Transfer_Learning_Visualizations.R")

} else {
  cat("  Transfer Learning SKIPPED (USE_TRANSFER_LEARNING = FALSE)\n")
  cat("  Set USE_TRANSFER_LEARNING <- TRUE in blue_carbon_config.R to enable.\n\n")
}

# ============================================================================
# ── WORKFLOW COMPLETE ──────────────────────────────────────────────────────
# ============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  WORKFLOW COMPLETE                                       ║\n")
cat("╠══════════════════════════════════════════════════════════╣\n")
cat("║  Step 1 outputs: outputs/Basic_analysis/Step1_NonSpatial/             ║\n")
if (aoi_available && isTRUE(RUN_STAGE2)) {
  cat("║  Step 2 outputs: outputs/Basic_analysis/Step2_Spatial_Basic/SimpleExtrapolation/  ║\n")
  cat("║                  outputs/Basic_analysis/Step2_Spatial_Basic/Kriging/               ║\n")
}
if (aoi_available && covariate_found && isTRUE(RUN_STAGE3)) {
  cat("║  Step 3 outputs: outputs/Advanced_analysis/Step3_Spatial_Advanced/        ║\n")
}
cat("╚══════════════════════════════════════════════════════════╝\n\n")
