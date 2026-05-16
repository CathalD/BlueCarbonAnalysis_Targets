## ============================================================================
## BLUE CARBON MMRV WORKFLOW: LOCAL EXECUTION
## ============================================================================
## LOCATION: Analysis_Workflow/MasterScript_Local.R
##
## USE THIS SCRIPT when running the workflow locally (desktop/laptop).
## It sources all scripts directly from the local filesystem — no GitHub
## connection required.
##
## HOW TO RUN:
##   Option A — RStudio:
##     Open this file and click "Source" (top-right of editor pane)
##
##   Option B — R console:
##     setwd("/path/to/BlueCarbon_Workflow_V1.0")
##     source("Analysis_Workflow/MasterScript_Local.R")
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
##    Requires: AOI GeoJSON/SHP with stratum column (set AOI_FILE in config)
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
## For GitHub-based deployment (downloads scripts from repo), use:
##   Analysis_Workflow/MasterScript_GitHub.R
## ============================================================================

# ── Set working directory to project root (BlueCarbon_Workflow_V1.0/) ────────
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
# ── LOAD CONFIGURATION ───────────────────────────────────────────────────
# ============================================================================

source("Analysis_Workflow/blue_carbon_config.R")

# ============================================================================
# ── STEP 1: NON-SPATIAL ANALYSIS (always runs) ───────────────────────────
# ── Requires: core_locations.csv + core_samples.csv
# ── Outputs:  outputs/Basic_analysis/Step1_NonSpatial/Tables/ and Plots/
# ============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  STEP 1 — NON-SPATIAL ANALYSIS                          ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

# S1_00: Directory setup — creates all output folders
# Run this once before the first time you run the full workflow
source("Analysis_Workflow/Basic Analysis/Step1_NonSpatial/S1_00_setup_directories.R")

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

  # S3_02: Carbon stock totals — aggregate depth layers, calculate total Mg C per AOI
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
# ============================================================================

if (exists("USE_BAYESIAN") && isTRUE(USE_BAYESIAN)) {
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  BAYESIAN MODEL                                          ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n\n")

  source("Analysis_Workflow/Bayesian model/P3_0c_bayesian_prior_setup_bluecarbon.R")
  source("Analysis_Workflow/Bayesian model/P3_1c_bayesian_sampling_design_bluecarbon.R")
  source("Analysis_Workflow/Bayesian model/P3_06c_bayesian_posterior_estimation_bluecarbon.R")
  source("Analysis_Workflow/Bayesian model/P3_07b_comprehensive_standards_report.R")

} else {
  cat("  Bayesian model SKIPPED (USE_BAYESIAN = FALSE)\n")
  cat("  Set USE_BAYESIAN <- TRUE in blue_carbon_config.R to enable.\n\n")
}

# ============================================================================
# ── TRANSFER LEARNING (optional — set USE_TRANSFER_LEARNING = TRUE) ───────
# ============================================================================

if (exists("USE_TRANSFER_LEARNING") && isTRUE(USE_TRANSFER_LEARNING)) {
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  TRANSFER LEARNING                                       ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n\n")

  # Step 4a: Harmonize global Janousek data to VM0033 depths
  # Inputs:  data_raw/janousek_layers.csv + data_raw/CorePoints_Covariates_BC_Canada.csv
  # Outputs: data_processed/global_cores_harmonized_VM0033.csv
  cat("-- Step 4a: Global Depth Harmonization --\n")
  source("Analysis_Workflow/Transfer Learning/P4_3b_Depth_Harmonization_Global.R")

  # Step 4b: 3-Model Transfer Learning Comparison (GEE / Embeddings / Combined)
  # Inputs:  data_processed/global_cores_harmonized_VM0033.csv
  #          data_raw/CorePoints_Covariates_BC_Canada.csv  (Models 1+3: GEE covariates)
  #          data_raw/CorePoints_EmbeddingSimilarity_BC_Canada.csv  (Models 2+3)
  #          data_processed/cores_harmonized_bluecarbon.csv  (local field data)
  #          COVARIATE_RASTER  (optional — spatial predictions for Models 1+3)
  #          EMBEDDING_RASTER  (optional — spatial predictions for Models 2+3)
  # Outputs: outputs/Transfer_learning/Model1_GEE/
  #          outputs/Transfer_learning/Model2_Embeddings/
  #          outputs/Transfer_learning/Model3_Combined/
  #          outputs/Transfer_learning/Comparison/model_comparison_all_depths.csv
  cat("-- Step 4b: Transfer Learning (3-Model Comparison) --\n")
  source("Analysis_Workflow/Transfer Learning/S4_01_transfer_learning.R")

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
