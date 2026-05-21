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


# ---------------------------------------------------------------------------
# profiles_to_sf()
# ---------------------------------------------------------------------------
# Converts a data.frame with latitude/longitude columns to an sf POINT
# object in WGS84, keeping only profile_id and dataset as properties
# (to keep GEE feature collections lean).
# ---------------------------------------------------------------------------
profiles_to_sf <- function(profiles_df) {
  suppressPackageStartupMessages(library(sf))

  sf::st_as_sf(
    profiles_df[, c("profile_id", "dataset", "longitude", "latitude")],
    coords = c("longitude", "latitude"),
    crs    = 4326L
  )
}
