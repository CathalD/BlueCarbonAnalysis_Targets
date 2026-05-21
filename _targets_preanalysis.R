# =============================================================================
# _targets_preanalysis.R
# Pre-analysis pipeline: global data preparation + GEE covariate extraction
#
# Run with:
#   targets::tar_make(script = "_targets_preanalysis.R",
#                     store  = "_targets_preanalysis")
#
# Prerequisites:
#   1. Install and authenticate rgee (run once):
#        source("Pre-Analysis Data Preparation/GEE_R/00_setup_rgee.R")
#   2. Run the main pipeline first to ensure packages are available.
#
# Pipeline structure:
#   Phase 1 — Global data prep (fast, no GEE needed):
#     janousek_harmonized   : Janousek layers + profiles (full harmonization)
#     profiles_for_gee      : EM + SG profiles with valid lat/lon
#
#   Phase 2 — GEE covariate extraction (requires rgee + authentication):
#     gee_climate           : TerraClimate MAT_C + MAP_mm (2000–2022)
#     gee_topo              : 7 topography/channel bands (30 m)
#     gee_sar               : Sentinel-1 VV/VH composite (2020–2023, 30 m)
#     gee_ndvi_stddev       : NDVI seasonal variability (2020–2023, 30 m)
#     gee_s2                : S2 raw (9 bands) + derived (5 bands) in one pass (2020–2023, 30 m)
#
#   Phase 3 — Combine + output:
#     global_covariates     : Merged 27-band canonical data.frame
#     covariates_file       : CorePoints_Covariates_BC_Canada.csv
#
# Extensibility:
#   To add a new data source (e.g. WoSIS):
#     1. Add ingest_wosis() to R/preanalysis/global_data.R
#     2. Add tar_target(wosis_harmonized, ...) below
#     3. Update filter_for_gee() to accept combined input, or add a separate
#        filter target and union the profile lists before GEE extraction.
# =============================================================================

library(targets)
library(tarchetypes)

# Source all preanalysis R modules
lapply(
  list.files("R/preanalysis", pattern = "\\.R$", full.names = TRUE),
  source
)

# ── Test mode ─────────────────────────────────────────────────────────────────
# Set TEST_MODE <- TRUE to run on 10 random points from the real dataset before
# committing to the full extraction. Flip to FALSE for the production run.
TEST_MODE <- TRUE
TEST_N    <- 50L   # 50 pts → 5 geographic batches of 10 after sort; keeps each bbox local
# Found in the Python notebook: ee.Initialize(project='...')
GEE_PROJECT <- "north-star-project-470316"

# ── Pipeline ──────────────────────────────────────────────────────────────────
list(

  # ── Phase 1: Raw data file inputs ──────────────────────────────────────────
  # Tracked as files — pipeline reruns if the source CSVs change.

  tar_target(
    jan_locations_file,
    file.path("Pre-Analysis Data Preparation", "data_global",
              "JANOUSEK_DATA", "Global_Core_Locations.csv"),
    format = "file"
  ),

  tar_target(
    jan_samples_file,
    file.path("Pre-Analysis Data Preparation", "data_global",
              "JANOUSEK_DATA", "Global_Core_Samples.csv"),
    format = "file"
  ),

  # ── Phase 1: Janousek harmonization ────────────────────────────────────────
  tar_target(
    janousek_harmonized,
    ingest_janousek(jan_locations_file, jan_samples_file)
  ),

  # ── Phase 1: Filter profiles for GEE (EM + SG ecosystems, valid lat/lon) ───
  tar_target(
    profiles_for_gee,
    filter_for_gee(janousek_harmonized, ecosystems = c("EM", "SG"))
  ),

  # ── Phase 1: Subsample for test mode ───────────────────────────────────────
  tar_target(
    profiles_for_extraction,
    if (TEST_MODE) {
      set.seed(42L)
      profiles_for_gee[sample(nrow(profiles_for_gee), min(TEST_N, nrow(profiles_for_gee))), ]
    } else {
      profiles_for_gee
    }
  ),

  # ── Phase 2: GEE covariate extraction ──────────────────────────────────────
  # Each target is independently cached — a re-run after a partial failure
  # only re-extracts the failed groups.
  #
  # NOTE: TerraClimate (cheap, 4 km) runs first so the user can inspect
  # climate values before committing to the expensive S2 extraction.

  tar_target(
    gee_climate,
    extract_climate(profiles_for_extraction, gee_project = GEE_PROJECT)
  ),

  tar_target(
    gee_topo,
    extract_topo(profiles_for_extraction, gee_project = GEE_PROJECT)
  ),

  tar_target(
    gee_sar,
    extract_sar(profiles_for_extraction, gee_project = GEE_PROJECT)
  ),

  tar_target(
    gee_ndvi_stddev,
    extract_ndvi_stddev(profiles_for_extraction, gee_project = GEE_PROJECT)
  ),

  # S2 raw reflectance (9 bands) + derived indices (5 bands) in a single target.
  # One S2 median composite is built per batch — previously two separate targets
  # doubled the number of GEE compute calls for the same imagery.
  tar_target(
    gee_s2,
    extract_s2_all(profiles_for_extraction, gee_project = GEE_PROJECT)
  ),

  # ── Phase 3: Combine all extractions ───────────────────────────────────────
  # Enforces canonical 27-band column order; fills failed bands with NA + warning.

  tar_target(
    global_covariates,
    combine_covariates(
      profiles_for_extraction,
      topo    = gee_topo,
      sar     = gee_sar,
      ndvi_sd = gee_ndvi_stddev,
      s2      = gee_s2,
      climate = gee_climate
    )
  ),

  # ── Phase 3: Write output CSV ───────────────────────────────────────────────
  # Overwrites the existing CorePoints_Covariates_BC_Canada.csv used by the
  # transfer learning pipeline (Step 4).

  tar_target(
    covariates_file,
    write_covariates_csv(
      global_covariates,
      path = file.path("Pre-Analysis Data Preparation", "data_raw",
                       "CorePoints_Covariates_BC_Canada.csv")
    ),
    format = "file"
  )

)
