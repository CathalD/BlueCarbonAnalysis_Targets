# ============================================================================
# MODULE 00b: BLUE CARBON - DIRECTORY SETUP
# ============================================================================
# PURPOSE: Create only the directories actually used by workflow modules
# USAGE:   source("Analysis_Workflow/Basic Analysis/P1_0b_setup_directories.R")
#          from BlueCarbon_Workflow_V1.0/ working directory
# LOCATION: Analysis_Workflow/Basic Analysis/
# ============================================================================

# ── PATH RESOLVER ──────────────────────────────────────────────────────────
# Sets working directory to BlueCarbon_Workflow_V1.0/ (project root) so all
# relative data paths resolve correctly regardless of how the script is invoked.
local({
  # Method 1: called via source() — detect this script's location
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) {
    # Script is at Analysis_Workflow/Basic Analysis/ → go up 2 levels
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
  message("⚠ Could not auto-detect project root. Please run:\n",
          "  setwd('/path/to/BlueCarbon_Workflow_V1.0')\n",
          "before sourcing this script.")
})
cat(sprintf("Working directory: %s\n\n", getwd()))
# ──────────────────────────────────────────────────────────────────────────

cat("\n========================================\n")
cat("BLUE CARBON - DIRECTORY SETUP\n")
cat("========================================\n\n")

# ============================================================================
# PROJECT PATH DEFINITIONS
# ============================================================================
# All paths are relative to BlueCarbon_Workflow_V1.0/ (project root).
# Pre-Analysis data lives in a dedicated subfolder; analysis outputs live
# directly under the project root for easy access from all workflow modules.

PRE_ANALYSIS_DIR <- "Pre-Analysis Data Preparation"   # Pre-Analysis subfolder
DATA_RAW_DIR     <- file.path(PRE_ANALYSIS_DIR, "data_raw")      # field core CSVs
DATA_GLOBAL_DIR  <- file.path(PRE_ANALYSIS_DIR, "data_global")   # global DB + GEE outputs
COVARIATES_DIR   <- file.path(PRE_ANALYSIS_DIR, "covariates")    # GEE multi-band TIF

# ============================================================================
# CREATE DIRECTORY STRUCTURE
# ============================================================================

cat("Creating workflow directories...\n\n")

# Core processing directories (written by analysis modules)
core_dirs <- c(
  "data_processed",                   # Processed R objects (all modules write here)
  "logs"                              # Log files (all modules write here)
)

# Output directories (written by Modules 02–07)
output_dirs <- c(
  "outputs",
  "outputs/plots",
  "outputs/plots/by_stratum",                          # Module 03: harmonization plots
  "outputs/plots/exploratory",                         # Module 02: EDA figures
  "outputs/models",
  "outputs/models/kriging",                            # Module 04: saved kriging models
  "outputs/models/rf",                                 # Module 05: saved RF models
  "outputs/models/large_scale_bluecarbon",             # Transfer learning: global models
  "outputs/predictions",
  "outputs/predictions/kriging",                                  # Module 04: kriging (parent)
  "outputs/predictions/kriging/carbon_stock_maps",                # Module 04: kg/m² predictions
  "outputs/predictions/kriging/standard_error_maps",              # Module 04: SE maps
  "outputs/predictions/kriging/variance_maps",                    # Module 04: variance maps
  "outputs/predictions/kriging/aggregated_stocks",                # Module 04: Mg C/ha layer stocks
  "outputs/predictions/rf",                                       # Module 05: RF predictions
  "outputs/carbon_stocks",                             # Module 06: stock calculations
  "outputs/carbon_stocks/maps",                        # Module 06: stock raster maps
  report_dir(),                                         # Module 07: assessment/verification report package
  "outputs/Basic_analysis",
  "outputs/Basic_analysis/Carbon_Stock_Calculations",  # Tables & stock calculations
  "outputs/Basic_analysis/Spatial_Maps_by_Depth",      # Depth-specific kriging maps
  "outputs/Basic_analysis/Exploratory_Data_Plots",     # EDA plots and QA visuals
  "outputs/Advanced_spatial_analysis",
  "outputs/Advanced_spatial_analysis/RF_Maps_by_Depth",            # RF depth maps + AOA
  "outputs/Advanced_spatial_analysis/RF_Diagnostics_and_Importance", # RF diagnostics
  "outputs/Bayesian_analysis",
  "outputs/Bayesian_analysis/Posterior_Distributions", # Posterior distributions
  "outputs/Bayesian_analysis/Prior_vs_Posterior_Maps", # Prior/likelihood/posterior layers
  "outputs/Transfer_learning",
  "outputs/Transfer_learning/Global_vs_Local_Comparisons" # Transfer: comparison rasters
)

# Diagnostic directories (Modules 01, 03–05)
diagnostic_dirs <- c(
  "diagnostics",
  "diagnostics/data_prep",                    # Module 01: data prep diagnostics
  "diagnostics/qaqc",                         # Module 01: QA/QC reports
  "diagnostics/variograms",                   # Module 04: variogram plots
  "diagnostics/crossvalidation",              # Modules 04 & 05: CV metrics
  "diagnostics/large_scale_bluecarbon"        # Transfer learning diagnostics
)

# Combine all required directories
required_dirs <- c(core_dirs, output_dirs, diagnostic_dirs)

# Create directories (all relative to project root)
created <- 0
existed <- 0

for (dir_name in required_dirs) {
  if (!dir.exists(dir_name)) {
    dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)
    created <- created + 1
    cat(sprintf("  ✓ Created: %s\n", dir_name))
  } else {
    existed <- existed + 1
  }
}

cat(sprintf("\nDirectories created: %d\n", created))
cat(sprintf("Already existed:     %d\n", existed))
cat(sprintf("Total:               %d\n\n", created + existed))

# ============================================================================
# CHECK FOR CONFIGURATION FILE
# ============================================================================

cat("========================================\n")
cat("CHECKING CONFIGURATION\n")
cat("========================================\n\n")

config_path <- "Analysis_Workflow/blue_carbon_config.R"

if (file.exists(config_path)) {
  cat(sprintf("✓ Configuration file found: %s\n", config_path))
  tryCatch({
    source(config_path, local = TRUE)
    cat("✓ Configuration file validated successfully\n\n")
  }, error = function(e) {
    cat("⚠️  Configuration file has errors:\n")
    cat(sprintf("   %s\n\n", e$message))
  })
} else {
  cat(sprintf("✗ Configuration file NOT FOUND: %s\n", config_path))
  cat("  → Expected at: BlueCarbon_Workflow_V1.0/Analysis_Workflow/blue_carbon_config.R\n\n")
}

# ============================================================================
# CHECK FOR FIELD DATA FILES
# ============================================================================

cat("========================================\n")
cat("CHECKING FOR FIELD DATA FILES\n")
cat("========================================\n\n")

cat(sprintf("Looking in: %s\n\n", DATA_RAW_DIR))

# Expected field data files (in Pre-Analysis Data Preparation/data_raw/)
data_files_to_check <- list(
  core_locations = file.path(DATA_RAW_DIR, "core_locations.csv"),
  core_samples   = file.path(DATA_RAW_DIR, "core_samples.csv")
)

data_status <- list()

for (file_type in names(data_files_to_check)) {
  path  <- data_files_to_check[[file_type]]
  found <- file.exists(path)
  data_status[[file_type]] <- list(found = found, path = path)

  if (found) {
    n_rows <- tryCatch(nrow(read.csv(path, nrows = 1)), error = function(e) NA)
    cat(sprintf("  ✓ Found: %s\n", path))
  } else {
    cat(sprintf("  ✗ Missing: %s\n", path))
    cat(sprintf("    → Place your field data CSV here: %s\n", path))
  }
}

cat("\n")

# ============================================================================
# CHECK FOR GEE COVARIATE RASTER
# ============================================================================

cat("========================================\n")
cat("CHECKING FOR GEE COVARIATE RASTER\n")
cat("========================================\n\n")

cat(sprintf("Looking in: %s\n\n", COVARIATES_DIR))

# Single multi-band GeoTIFF exported from GEE JS script
# (27-band canonical stack: elevation, SAR, Sentinel-2, climate)
covariate_files <- list.files(COVARIATES_DIR,
                               pattern  = "\\.tif$",
                               recursive = FALSE,
                               full.names = TRUE)

if (length(covariate_files) > 0) {
  cat(sprintf("✓ Found %d covariate TIF file(s):\n", length(covariate_files)))
  for (f in covariate_files) {
    size_mb <- round(file.info(f)$size / 1e6, 1)
    cat(sprintf("   %s  (%.1f MB)\n", basename(f), size_mb))
  }

  # Try to read band count if terra is available
  if (requireNamespace("terra", quietly = TRUE)) {
    tryCatch({
      r <- terra::rast(covariate_files[1])
      cat(sprintf("\n  Bands detected: %d\n", terra::nlyr(r)))
      cat(sprintf("  Band names: %s\n",
                  paste(names(r)[seq_len(min(5, terra::nlyr(r)))],
                        collapse = ", "),
                  if (terra::nlyr(r) > 5) "..." else ""))
    }, error = function(e) NULL)
  }
  cat("\n")

} else {
  cat("⚠️  No GEE covariate TIF found.\n")
  cat("   Module 05 (Random Forest) requires the covariate raster.\n\n")
  cat("   To generate it:\n")
  cat("   1. Run GoogleEarthEngineAOICovariateAnalysis.js in the GEE Code Editor\n")
  cat("   2. Export the 27-band covariate snapshot from Step ④\n")
  cat(sprintf("   3. Place the downloaded TIF in: %s\n\n", COVARIATES_DIR))
  cat("   Module 04 (Kriging) can run without covariates.\n\n")
}

# ============================================================================
# ADDITIONAL DIRECTORIES CREATED BY MODULES (INFO ONLY)
# ============================================================================

cat("========================================\n")
cat("NOTE: ADDITIONAL DIRECTORIES\n")
cat("========================================\n\n")

cat("Some directories are created automatically by modules when needed:\n")
cat("  - diagnostics/variable_importance/  (Module 05)\n")
cat("  - outputs/predictions/aoa/          (Module 05, if AOA enabled)\n")
cat(sprintf("  - %s/spatial_exports/  (Module 07)\n\n", report_dir()))

# ============================================================================
# SAVE SETUP SUMMARY
# ============================================================================

setup_summary <- list(
  date                = Sys.Date(),
  r_version           = paste(R.version$major, R.version$minor, sep = "."),
  working_directory   = getwd(),
  directories_created = created,
  directories_existed = existed,
  config_path         = config_path,
  config_exists       = file.exists(config_path),
  data_raw_dir        = DATA_RAW_DIR,
  core_locations_found = data_status$core_locations$found,
  core_locations_path  = data_files_to_check$core_locations,
  core_samples_found  = data_status$core_samples$found,
  core_samples_path   = data_files_to_check$core_samples,
  covariates_dir      = COVARIATES_DIR,
  covariates_found    = length(covariate_files),
  covariate_files     = covariate_files
)

if (!dir.exists("data_processed")) dir.create("data_processed", recursive = TRUE)
saveRDS(setup_summary, "data_processed/setup_summary.rds")

# Write log
if (!dir.exists("logs")) dir.create("logs", showWarnings = FALSE)
log_file <- file.path("logs", paste0("setup_", Sys.Date(), ".txt"))
sink(log_file)
cat("Blue Carbon Workflow Setup Summary\n")
cat("===================================\n\n")
cat(sprintf("Date:                %s\n",   Sys.Date()))
cat(sprintf("R version:           %s\n",   setup_summary$r_version))
cat(sprintf("Working directory:   %s\n\n", setup_summary$working_directory))
cat(sprintf("Directories created: %d\n",   created))
cat(sprintf("Directories existed: %d\n",   existed))
cat(sprintf("Config file:         %s\n",   ifelse(setup_summary$config_exists, "Found", "Missing")))
cat(sprintf("Core locations:      %s\n",   ifelse(setup_summary$core_locations_found, setup_summary$core_locations_path, "Missing")))
cat(sprintf("Core samples:        %s\n",   ifelse(setup_summary$core_samples_found, setup_summary$core_samples_path, "Missing")))
cat(sprintf("Covariate TIF:       %d file(s) in %s\n", setup_summary$covariates_found, setup_summary$covariates_dir))
sink()

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("========================================\n")
cat("SETUP SUMMARY\n")
cat("========================================\n\n")

cat(sprintf("R version:         %s\n", setup_summary$r_version))
cat(sprintf("Working directory: %s\n\n", getwd()))

cat("Status:\n")
cat(sprintf("  %s Directory structure:  %d folders\n",
            "✓", created + existed))
cat(sprintf("  %s Configuration file:  %s\n",
            ifelse(setup_summary$config_exists, "✓", "✗"),
            ifelse(setup_summary$config_exists, "Found", "MISSING")))
cat(sprintf("  %s Field data (raw):    %s\n",
            ifelse(setup_summary$core_locations_found && setup_summary$core_samples_found, "✓", "⚠️"),
            ifelse(setup_summary$core_locations_found && setup_summary$core_samples_found,
                   sprintf("Found in %s", DATA_RAW_DIR),
                   sprintf("Incomplete — check %s", DATA_RAW_DIR))))
cat(sprintf("  %s GEE covariate TIF:  %s\n",
            ifelse(setup_summary$covariates_found > 0, "✓", "⚠️"),
            ifelse(setup_summary$covariates_found > 0,
                   sprintf("%d file(s) in %s", setup_summary$covariates_found, COVARIATES_DIR),
                   sprintf("Not found — export from GEE and place in %s", COVARIATES_DIR))))
cat("\n")

# ============================================================================
# NEXT STEPS
# ============================================================================

cat("========================================\n")
cat("NEXT STEPS\n")
cat("========================================\n\n")

ready_for_basic    <- setup_summary$config_exists &&
                      setup_summary$core_locations_found &&
                      setup_summary$core_samples_found
ready_for_advanced <- ready_for_basic && setup_summary$covariates_found > 0

if (ready_for_basic) {
  cat("✓ READY for Basic Analysis (Modules 01–04: data prep → kriging)\n\n")
  cat("  source('Analysis_Workflow/Basic Analysis/P1_01_data_prep_bluecarbon.R')\n")
  cat("  source('Analysis_Workflow/Basic Analysis/P2_02_exploratory_analysis_bluecarbon.R')\n")
  cat("  source('Analysis_Workflow/Basic Analysis/P2_3a_Depth_Harmonization_Local.R')\n")
  cat("  source('Analysis_Workflow/Basic Analysis/P2_04_raster_predictions_kriging_bluecarbon.R')\n\n")
} else {
  cat("⚠️  SETUP INCOMPLETE for Basic Analysis:\n\n")
  if (!setup_summary$config_exists)
    cat(sprintf("  ✗ Config missing → expected at: %s\n\n", config_path))
  if (!setup_summary$core_locations_found)
    cat(sprintf("  ✗ core_locations.csv missing → place at: %s\n\n", setup_summary$core_locations_path))
  if (!setup_summary$core_samples_found)
    cat(sprintf("  ✗ core_samples.csv missing → place at: %s\n\n", setup_summary$core_samples_path))
}

if (ready_for_advanced) {
  cat("✓ READY for Advanced Analysis (Module 05: Random Forest + covariates)\n\n")
} else if (ready_for_basic) {
  cat(sprintf("⚠️  Advanced Analysis (RF model) requires covariate TIF:\n"))
  cat(sprintf("   → Export from GEE JS script and place in: %s\n\n", COVARIATES_DIR))
}

cat(sprintf("Setup log saved to: %s\n\n", log_file))
cat("Done! 🌊\n\n")
