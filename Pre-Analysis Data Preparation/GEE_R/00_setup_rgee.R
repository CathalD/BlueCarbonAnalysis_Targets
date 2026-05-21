# =============================================================================
# 00_setup_rgee.R
# One-time setup script for rgee + Google Earth Engine authentication.
#
# Run this script ONCE before using the preanalysis pipeline:
#   source("Pre-Analysis Data Preparation/GEE_R/00_setup_rgee.R")
#
# After running this script successfully, you can run the pipeline:
#   targets::tar_make(script = "_targets_preanalysis.R",
#                     store  = "_targets_preanalysis")
#
# Prerequisites (install once if not already done):
#   install.packages(c("rgee", "reticulate"))
#   reticulate::install_miniconda()     # only if conda not already available
# =============================================================================

cat("=================================================================\n")
cat("  rgee + Google Earth Engine setup\n")
cat("=================================================================\n\n")

# ── Step 1: Check packages ────────────────────────────────────────────────────
cat("Step 1: Checking package availability...\n")

for (pkg in c("rgee", "reticulate", "sf", "dplyr", "readr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    install.packages(pkg)
  }
}
cat("  All required packages available.\n\n")


# ── Step 2: Python environment ────────────────────────────────────────────────
# rgee uses the Python earthengine-api package via reticulate.
# ee_install() creates an isolated conda environment named 'rgee_py'.

cat("Step 2: Setting up Python environment for rgee...\n")
library(rgee)

tryCatch({
  # Check if the rgee Python env already exists
  ee_check_python()
  cat("  Python environment OK (already configured).\n\n")
}, error = function(e) {
  cat("  Installing rgee Python dependencies (earthengine-api)...\n")
  cat("  This may take 2–5 minutes the first time.\n\n")
  rgee::ee_install(py_env = "rgee_py")
  cat("  Python environment created.\n\n")
})


# ── Step 3: GEE authentication ────────────────────────────────────────────────
# Opens a browser window to authenticate with your Google account.
# You must have access to the GEE project below.

GEE_PROJECT <- "north-star-project-470316"  # update if needed

cat("Step 3: Authenticating with Google Earth Engine...\n")
cat(sprintf("  Project: %s\n", GEE_PROJECT))
cat("  A browser window will open. Sign in with your Google account\n")
cat("  and grant access to Google Earth Engine.\n\n")

tryCatch({
  rgee::ee_Initialize(user    = "cathalpdoherty@gmail.com",
                      project = GEE_PROJECT,
                      drive   = FALSE,
                      gcs     = FALSE)
  cat("\n  Authentication successful!\n\n")
}, error = function(e) {
  cat(sprintf("\n  Authentication failed: %s\n\n", conditionMessage(e)))
  cat("  Troubleshooting:\n")
  cat("  1. Make sure your GEE account is activated at https://earthengine.google.com\n")
  cat("  2. Verify the project ID is correct: ", GEE_PROJECT, "\n")
  cat("  3. Try manually: library(rgee); ee_Initialize(user='your@email.com')\n")
  stop("Setup failed — see messages above.", call. = FALSE)
})


# ── Step 4: Verify connection ─────────────────────────────────────────────────
cat("Step 4: Verifying GEE connection...\n")

tryCatch({
  test_img <- ee$Image("NASA/NASADEM_HGT/001")
  band_names <- test_img$bandNames()$getInfo()
  cat(sprintf("  NASA NASADEM DEM accessible — bands: %s\n",
              paste(band_names, collapse = ", ")))

  test_s2 <- ee$ImageCollection("COPERNICUS/S2_SR_HARMONIZED")$
    filterDate("2023-01-01", "2023-01-07")$
    size()$getInfo()
  cat(sprintf("  Sentinel-2 collection accessible (%d scenes in test window)\n", test_s2))

  cat("\n  ✓ GEE connection verified. Ready to run the pipeline!\n\n")
}, error = function(e) {
  cat(sprintf("\n  Connection test failed: %s\n", conditionMessage(e)))
  cat("  Check your internet connection and GEE project access.\n")
  stop("Verification failed.", call. = FALSE)
})


# ── Step 5: Quick reference ───────────────────────────────────────────────────
cat("=================================================================\n")
cat("  Setup complete. To run the pre-analysis pipeline:\n\n")
cat("    targets::tar_make(\n")
cat("      script = '_targets_preanalysis.R',\n")
cat("      store  = '_targets_preanalysis'\n")
cat("    )\n\n")
cat("  To check pipeline status:\n")
cat("    targets::tar_visnetwork(\n")
cat("      script = '_targets_preanalysis.R',\n")
cat("      store  = '_targets_preanalysis'\n")
cat("    )\n\n")
cat("  To re-run GEE extraction for a single group (e.g. after a failure):\n")
cat("    targets::tar_make(\n")
cat("      names  = 'gee_s2_raw',\n")
cat("      script = '_targets_preanalysis.R',\n")
cat("      store  = '_targets_preanalysis'\n")
cat("    )\n")
cat("=================================================================\n")
