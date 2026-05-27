# Blue Carbon Analysis — Step-by-step workflow
#
# Run each section in order in the RStudio Console.
#
# IMPORTANT: shiny::runApp() blocks R until you close the browser tab.
# Complete setup in the app, close it, then run the pipeline commands below.
# ──────────────────────────────────────────────────────────────────────────────


# ── 0. FIRST-TIME SETUP ───────────────────────────────────────────────────────

# Pull the latest version of the project from GitHub:
system("git pull")

# Open the setup wizard — complete all 4 steps, click Save Setup Files, close:
shiny::runApp("shiny")
# Shortcut: app()

# After closing the app, your config and CSV files are in place.
# Open the Run tab next time you launch the app to see pipeline status.


# ── 1. NON-SPATIAL ANALYSIS (always run first, ~5–15 min) ────────────────────
# Produces: per-stratum carbon stocks, depth profiles, VM0033 estimates
# Report:   reports/step1_nonspatial.html

targets::tar_make()
# Shortcut: tm()

# Check what will run before you start:
targets::tar_visnetwork()

# If something fails, check what went wrong:
targets::tar_meta() |> dplyr::filter(!is.na(error)) |> dplyr::select(name, error)


# ── 2. RF SPATIAL MAPS (optional — requires covariate raster, ~5–10 min) ─────
# Produces: 25-m carbon stock map, variable importance
# Report:   reports/step3_random_forest.html

targets::tar_make(script = "_targets_rf.R", store = "_targets_rf")
# Shortcut: tmrf()


# ── 3. GEE GLOBAL COVARIATE EXTRACTION (needed for TL methods, ~60 min) ──────
# Run once — safe to re-run, completed batches are skipped.

# Authenticate GEE first (a browser window will open — one-time only):
library(rgee)
ee_Initialize(user = "your.email@gmail.com", drive = TRUE)

targets::tar_make(script = "_targets_preanalysis.R", store = "_targets_preanalysis")


# ── 4. WADOUX TRANSFER LEARNING (~15 min) ─────────────────────────────────────
# Requires: Step 3 complete
# Produces: bias-corrected prediction maps with bootstrap uncertainty
# Report:   reports/step4_transfer_learning.html

targets::tar_make(script = "_targets_transfer.R", store = "_targets_transfer")


# ── 5. EMBEDDING TRANSFER LEARNING (~30 min, optional) ───────────────────────
# Requires: Step 3 complete
# Produces: foundation-model weighted prediction maps + comparison to Step 4
# Report:   reports/step5_embedding_tl.html

targets::tar_make(script = "_targets_embedding.R", store = "_targets_embedding")


# ── BLUECARBON PRE-PROCESSING (optional — requires core_compaction.csv) ───────
# Corrects sample depths for percussion core compaction, estimates OC stocks
# using the BlueCarbon method, and validates short-core extrapolation quality.
# Produces core_samples_bc_processed.csv for use in the spatial pipeline.
# Requires: install.packages("BlueCarbon")

targets::tar_make(script = "_targets_bluecarbon.R", store = "_targets_bluecarbon")
# Shortcut: tmbc()

# To use compaction-corrected data in the main spatial pipeline:
# Copy outputs/bluecarbon/core_samples_bc_processed.csv → data_raw/core_samples.csv
# then re-run targets::tar_make() to invalidate and rebuild all downstream targets.


# ── SEQUESTRATION RATES (optional — requires core_chronology.csv) ─────────────
# Estimates mean OC sequestration rates (g C/m²/yr) over 100 years for cores
# with radiometric chronology data (^210Pb or ^14C age-depth control points).
# See Pre-Analysis Data Preparation/data_raw/core_chronology.csv for format.
# Requires: install.packages("BlueCarbon")

targets::tar_make(script = "_targets_seqrates.R", store = "_targets_seqrates")
# Shortcut: tmsr()


# ── VIEW RESULTS ──────────────────────────────────────────────────────────────
# Open the app and go to the Outputs tab to download reports:
shiny::runApp("shiny")

# Or open HTML reports directly in your browser — they are in reports/
# Spatial rasters (.tif) are in outputs/rf/, outputs/transfer/, outputs/embedding/
