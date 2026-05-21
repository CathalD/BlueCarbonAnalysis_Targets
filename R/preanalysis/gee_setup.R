# =============================================================================
# R/preanalysis/gee_setup.R
# rgee initialization helpers and GEE session management
# =============================================================================
#
# Before running the preanalysis pipeline, run:
#   source("Pre-Analysis Data Preparation/GEE_R/00_setup_rgee.R")
# or interactively:
#   library(rgee); ee_Initialize(user = "your@email.com")
#
# The extraction functions call initialize_gee() internally so a stale
# session is reconnected without requiring the user to re-run ee_Initialize.
# =============================================================================

# ---------------------------------------------------------------------------
# initialize_gee()
# ---------------------------------------------------------------------------
# Ensures rgee is loaded and GEE is authenticated.  Safe to call multiple
# times — reuses an existing session rather than re-authenticating.
#
# project : GEE cloud project ID (NULL uses the project stored in the token)
# ---------------------------------------------------------------------------
initialize_gee <- function(project = NULL) {
  suppressPackageStartupMessages(library(rgee))

  tryCatch({
    if (!is.null(project)) {
      rgee::ee_Initialize(project = project, quiet = TRUE)
    } else {
      rgee::ee_Initialize(quiet = TRUE)
    }
    invisible(TRUE)
  }, error = function(e) {
    stop(
      "[GEE] Failed to initialize rgee.\n",
      "  Run the setup script first:\n",
      "    source('Pre-Analysis Data Preparation/GEE_R/00_setup_rgee.R')\n",
      "  Or manually:\n",
      "    library(rgee); ee_Initialize(user = 'your@email.com')\n",
      "  Error: ", conditionMessage(e),
      call. = FALSE
    )
  })
}


# Note: conversion of profiles to ee.FeatureCollection is handled directly
# in gee_covariates.R via .df_to_ee_fc() — no sf/geojsonio dependency needed.
