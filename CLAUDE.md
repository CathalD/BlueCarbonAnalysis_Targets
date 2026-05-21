# Blue Carbon Analysis — Project Guide

## What this project does

Spatial prediction of soil organic carbon stocks in coastal blue carbon
ecosystems (tidal marshes, seagrass meadows) following the **VM0033 standard**
for carbon market reporting.

The pipeline takes field cores → harmonizes depths → extracts remote sensing
covariates globally and locally → predicts stocks across a full site raster
→ estimates uncertainty at 90% confidence.

---

## Tech stack

| Tool | Role |
|------|------|
| `targets` + `tarchetypes` + `geotargets` | Reproducible pipeline orchestration |
| `terra` + `sf` | Raster and vector spatial operations |
| `rgee` | R interface to Google Earth Engine |
| `ranger` | Random forest (transfer learning — supports case weights) |
| `randomForest` | Random forest (main RF pipeline) |
| `dplyr` / `tidyr` / `readr` | Data wrangling |
| `ggplot2` | All plots |
| Quarto (HTML) | Reports |

---

## Repository layout

```
BlueCarbonAnalysis_Targets/
├── _targets.R                    # Main pipeline (Steps 1–3)
├── _targets_transfer.R           # Transfer learning pipeline (Step 4 — Wadoux)
├── _targets_preanalysis.R        # Pre-analysis: global GEE covariate extraction
├── _targets_embedding.R          # Embedding TL pipeline (Step 5 — Model 2)
├── _targets.yaml                 # Named configs: main / transfer / preanalysis / embedding
├── blue_carbon_config.R          # Site-specific settings (edit per project)
├── R/                            # One .R file per analysis step
│   ├── config.R                  # load_config() — wraps blue_carbon_config.R
│   ├── data_prep.R               # load_raw_data()
│   ├── depth_harmonization.R     # harmonize_depths(), fit_hybrid_profile()
│   ├── exploratory_analysis.R    # run_eda()
│   ├── random_forest.R           # prepare_rf_data(), train_rf(), predict_rf_rasters()
│   ├── simple_extrapolation.R    # simple_extrapolation()
│   ├── summarise.R               # summarise_strata()
│   ├── transfer_learning.R       # harmonize_global_layers(), prepare_tl_data(),
│   │                             # train_tl(), predict_tl_rasters(), plot_tl_maps()
│   ├── preanalysis/              # GEE extraction modules
│   │   ├── global_data.R         # ingest_janousek(), filter_for_gee()
│   │   ├── gee_covariates.R      # extract_*(). combine_covariates(), write_covariates_csv()
│   │   └── gee_setup.R           # initialize_gee()
│   └── embedding_tl/             # Foundation model similarity modules
│       ├── gee_embeddings.R      # extract_global_embeddings(), extract_aoi_embedding_raster()
│       ├── embedding_similarity.R# compute_embedding_weights(), compute_pixel_similarity()
│       └── embedding_tl_model.R  # prepare_emb_tl_data(), train_emb_tl(),
│                                 # predict_emb_tl_rasters(), plot_emb_tl_maps()
├── reports/
│   ├── step1_nonspatial.qmd
│   ├── step3_random_forest.qmd
│   ├── step4_transfer_learning.qmd
│   └── step5_embedding_tl.qmd
├── Pre-Analysis Data Preparation/
│   ├── data_raw/                 # Local field data + GEE exports
│   │   ├── core_locations.csv    # core_id, latitude, longitude, stratum
│   │   ├── core_samples.csv      # core_id, depth_cm, soc_g_kg, bulk_density_g_cm3
│   │   └── CorePoints_Covariates_BC_Canada.csv  # GEE covariates at global cores
│   ├── data_global/              # Global coastal wetland database
│   │   ├── combined_layers_filtered.csv   # Janousek EM+SG layers (~124K rows)
│   │   └── JANOUSEK_DATA/                 # Raw Janousek files
│   └── covariates/
│       └── BlueCarbon_Covariate_Snapshot_25m_2020_2023.tif  # Local GEE raster
└── outputs/
    ├── rf/                       # RF prediction rasters
    ├── transfer/                 # Wadoux TL prediction rasters
    └── embedding/                # AOI embedding raster + embedding TL rasters
```

---

## Pipeline steps

### Run order

```
1. tar_make()                                          # main pipeline
2. tar_make(script="_targets_preanalysis.R",
            store="_targets_preanalysis")              # global GEE extraction
3. tar_make(script="_targets_transfer.R",
            store="_targets_transfer")                 # Wadoux TL (Model 1)
4. tar_make(script="_targets_embedding.R",
            store="_targets_embedding")               # Embedding TL (Model 2)
```

Steps 3 and 4 read `cores_harmonized` from the main store and
`CorePoints_Covariates_BC_Canada.csv` from the preanalysis output.

---

### Main pipeline (`_targets.R`)

| Step | Key targets | What it produces |
|------|-------------|------------------|
| 1 — Data prep | `cores_raw`, `eda_plots`, `cores_harmonized` | VM0033-depth harmonized field cores |
| 2 — Simple extrapolation | `step2_extrapolation` | Per-stratum carbon densities (no raster) |
| 3 — Random forest | `rf_data`, `rf_models`, `rf_rasters`, `rf_maps` | Spatial carbon stock maps, variable importance |
| Reports | `report_nonspatial`, `report_rf` | HTML reports in `reports/` |

---

### Pre-analysis pipeline (`_targets_preanalysis.R`)

Extracts a canonical 26-band covariate stack at all Janousek EM+SG core
locations via Google Earth Engine. Output feeds into both TL pipelines.

| Phase | Key targets | What it produces |
|-------|-------------|------------------|
| 1 — Janousek ingest | `janousek_harmonized`, `profiles_for_gee` | Filtered EM+SG profiles |
| 2 — GEE extraction | `gee_climate`, `gee_topo`, `gee_sar`, `gee_s2` | Per-group covariate data.frames |
| 3 — Combine | `global_covariates`, `covariates_file` | 26-band CSV at 952 profiles |

**Key settings in `_targets_preanalysis.R`:**
- `TEST_MODE <- FALSE` — production run (952 profiles, ~191 S2 batches)
- `TEST_MODE <- TRUE` / `TEST_N <- 50L` — fast validation run
- `GEE_PROJECT <- "north-star-project-470316"`
- S2 batch size: 5 points, 5-minute timeout, 1s inter-batch sleep

**Canonical 26 bands:**
Topo (7): `elevation_m`, `slope`, `elevationRelMHW`, `twi`, `dist_to_channel_m`, `tidal_flat_prob`, `coastal_dist_m`
SAR (3): `VV_mean`, `VH_mean`, `VVVH_ratio`
S2 optical (9): `B`, `G`, `R`, `B5`, `B6`, `B7`, `NIR`, `SWIR1`, `SWIR2`
S2 derived (5): `NDVI_median`, `LSWI_median`, `mNDWI_median`, `SAVI_median`, `tidal_wetness`
Climate (2): `MAT_C`, `MAP_mm`

---

### Transfer learning pipeline — Model 1 (`_targets_transfer.R`)

Wadoux instance weighting: a domain classifier RF estimates how similar each
global core is to the local site, then a weighted global RF is trained and
bias-corrected using local residuals.

| Step | Key targets | What it produces |
|------|-------------|------------------|
| Global harmonization | `global_harmonized` | VM0033-harmonized Janousek cores |
| TL data prep | `tl_data` | Global + local joined with bridge variables |
| TL models | `tl_models` | Per-depth Wadoux RF + bias + bootstrap CI |
| TL rasters | `tl_rasters` | 4-band GeoTIFFs per depth |
| TL maps | `tl_maps` | Maps + weight distribution + similarity heatmap |
| Report | `report_tl` | HTML report at `reports/step4_transfer_learning.html` |

**Four-stage training per depth:**
1. Wadoux domain classifier (RF probability forest on 6 bridge variables) → instance weights
2. Weighted global RF (1000 trees, `ranger`)
3. Mean bias correction from local residuals
4. 500-replicate bootstrap → bias SE

---

### Embedding transfer learning pipeline — Model 2 (`_targets_embedding.R`)

Replaces the Wadoux RF domain classifier with cosine similarity in the 64-d
space of `GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL` (2023–2025 average).
Stages B–D (weighted global RF, bias correction, bootstrap) are identical.

| Step | Key targets | What it produces |
|------|-------------|------------------|
| GEE extraction | `global_embeddings` | 64-d embedding at 952 global cores |
| AOI raster | `aoi_embedding_raster` | 64-band raster at `outputs/embedding/` |
| Weights | `emb_weights` | cosine_sim + weight per global core |
| Models | `emb_tl_models` | Per-depth embedding-weighted RF + bias |
| Rasters | `emb_tl_rasters` | 4-band GeoTIFFs per depth |
| Maps | `emb_tl_maps` | Maps + weight distribution + similarity heatmap |
| Report | `report_emb_tl` | HTML with Step 4 vs Step 5 comparison table |

**Embedding weight computation:**
- AOI mean 64-d vector computed from `aoi_embedding_raster`
- Cosine similarity per global core against AOI mean
- Sharpening: `weight = sim^5`, normalised so `mean(weight) = 1`
- AOI raster downloaded in 4 chunks of 16 bands to stay under GEE's 48 MB limit

**Output raster bands** (both TL pipelines):
`dX_Global_Prior`, `dX_Transfer_Final`, `dX_Local_Only`, `dX_Difference`
where X = depth midpoint with `.` replaced by `_` (e.g. `d7_5`, `d22_5`)

---

## VM0033 standards — must never change without discussion

```r
VM0033_DEPTH_MIDPOINTS <- c(7.5, 22.5, 40, 75)   # cm
VM0033_DEPTH_INTERVALS:
  0–15 cm   (midpoint 7.5,  thickness 15 cm)
  15–30 cm  (midpoint 22.5, thickness 15 cm)
  30–50 cm  (midpoint 40,   thickness 20 cm)
  50–100 cm (midpoint 75,   thickness 50 cm)

Carbon stock formula:
  carbon_stock_kg_m2 = SOC(g/kg) × BD(g/cm³) × thickness(cm) / 100

Uncertainty requirement: 90% prediction intervals (not 95%)
CV strategy: leave-one-CORE-out (not leave-one-observation-out)
```

Always read these from `cfg$VM0033_DEPTH_MIDPOINTS` — never hardcode.

---

## Code conventions

### Adding a new pipeline step

1. Create `R/<step_name>.R` with pure functions (no side effects beyond `message()`)
2. Add targets to `_targets.R` following the existing pattern:
   - `tar_target(file, "path", format = "file")` for inputs
   - `tar_terra_rast(name, fn(...))` for SpatRaster outputs
   - `tar_quarto(report_X, path = "reports/stepN_X.qmd")` for reports
3. Create `reports/stepN_<name>.qmd` — see existing reports for template

### R module style

- No side effects in functions (no `setwd()`, no `write_csv()` outside designated output functions)
- Exception: output functions (`predict_*_rasters()`, `extract_aoi_embedding_raster()`) may write to `outputs/`
- Use `message("[prefix] ...")` for progress, not `cat()` or `print()`
- `suppressPackageStartupMessages({ library(...) })` inside each function
- The `%||%` null-coalescing operator is available from `depth_harmonization.R`

### Quarto report template

```yaml
engine: knitr
format:
  html:
    toc: true
    toc-depth: 3
    embed-resources: true
    theme: cosmo
execute:
  echo: false
  warning: false
  message: false
```

Setup chunk must include:
```r
if (!file.exists("_targets")) setwd("..")
```

Read from the correct store:
```r
tar_read(target_name)                               # main store
tar_read(target_name, store = "_targets_transfer")  # Wadoux TL store
tar_read(target_name, store = "_targets_embedding") # Embedding TL store
```

---

## Transfer learning specifics

**Bridge variables (reduced — used when N_local < 15):**
`NDVI_median`, `LSWI_median`, `mNDWI_median`, `VV_mean`, `elevation_m`, `elevationRelMHW`

**Bridge variables (full — used when N_local ≥ 15):**
Above plus `slope`, `twi`, `dist_to_channel_m`, `tidal_flat_prob`, `coastal_dist_m`,
`VH_mean`, `VVVH_ratio`, `B`, `G`, `R`, `NIR`, `SWIR1`, `SWIR2`, `SAVI_median`

**Global data compound key:** `paste(dataset, profile_id, sep = "_")`
e.g. `"Janousek_1"`, `"WOSIS 2023_1144105"`

**Similarity plots stored in `tl_maps` / `emb_tl_maps`:**
- `$weights` — Wadoux/cosine weight distribution per depth (log scale jitter)
- `$heatmap` — covariate × core heatmap, rows sorted by weight, local cores pinned at top
- `$maps` — Global Prior vs Transfer Final rasters per depth
- `$validation` — LOCO CV R² and RMSE bar chart

---

## Config values to know

```r
cfg$PROJECT_NAME         # "BC_Coastal_BlueCarbon_2026_Example"
cfg$VM0033_DEPTH_MIDPOINTS
cfg$VM0033_DEPTH_INTERVALS   # data.frame with depth_top/bottom/midpoint/thickness_cm
cfg$BD_DEFAULTS          # list(IM=0.8, NM=0.8, MF=0.8) g/cm³
cfg$VALID_STRATA         # c("IM", "NM", "MF")
cfg$COVARIATE_RASTER     # path to local GEE raster
cfg$DATA_GLOBAL_DIR      # "Pre-Analysis Data Preparation/data_global"
cfg$BAND_LABELS          # named vector: raster band name → human-readable label
```

---

## Current build status

| Step | Status | Notes |
|------|--------|-------|
| Step 1 — Data prep + harmonization | Complete | |
| Step 2 — Simple extrapolation | Complete | |
| Step 3 — Random forest | Complete | |
| Pre-analysis — Global GEE extraction | Complete | 952 profiles, 97.7% complete covariates |
| Step 4 — Transfer learning (Wadoux) | Complete | Validated on full Janousek dataset |
| Step 5 — Transfer learning (Embeddings) | In progress | GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL 2023–2025 |

---

## What to build next

- Confirm Step 5 embedding pipeline runs end-to-end; compare LOCO CV RMSE vs Step 4
- Option 2 (cluster-based pixel-wise embedding weights): implement `compute_pixel_similarity()` → k-means → per-cluster weighted RF
- Consider `tar_make_future()` for parallel depth processing if runtime becomes limiting
