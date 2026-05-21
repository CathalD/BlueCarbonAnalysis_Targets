# Blue Carbon Analysis — Project Guide

## What this project does

Spatial prediction of soil organic carbon stocks in coastal blue carbon
ecosystems (tidal marshes, mangroves) following the **VM0033 standard**
for carbon market reporting.

The pipeline takes field cores → harmonizes depths → predicts stocks across a
full site raster → estimates uncertainty at 90% confidence.

---

## Tech stack

| Tool | Role |
|------|------|
| `targets` + `tarchetypes` + `geotargets` | Reproducible pipeline orchestration |
| `terra` + `sf` | Raster and vector spatial operations |
| `ranger` | Random forest (transfer learning — supports case weights) |
| `randomForest` | Random forest (main RF pipeline) |
| `dplyr` / `tidyr` / `readr` | Data wrangling |
| `ggplot2` | All plots |
| Quarto (HTML) | Reports |

---

## Repository layout

```
BlueCarbonAnalysis_Targets/
├── _targets.R                  # Main pipeline
├── _targets_transfer.R         # Transfer learning pipeline
├── _targets.yaml               # Named configs: main / transfer
├── blue_carbon_config.R        # Site-specific settings (edit per project)
├── R/                          # One .R file per analysis step
│   ├── config.R                # load_config() — wraps blue_carbon_config.R
│   ├── data_prep.R             # load_raw_data()
│   ├── depth_harmonization.R   # harmonize_depths(), fit_hybrid_profile()
│   ├── exploratory_analysis.R  # run_eda()
│   ├── random_forest.R         # prepare_rf_data(), train_rf(), predict_rf_rasters()
│   ├── simple_extrapolation.R  # simple_extrapolation()
│   ├── summarise.R             # summarise_strata()
│   └── transfer_learning.R     # harmonize_global_layers(), prepare_tl_data(),
│                               # train_tl(), predict_tl_rasters(), plot_tl_maps()
├── reports/
│   ├── step1_nonspatial.qmd
│   ├── step3_random_forest.qmd
│   └── step4_transfer_learning.qmd
├── Pre-Analysis Data Preparation/
│   ├── data_raw/               # Local field data + GEE exports
│   │   ├── core_locations.csv  # core_id, latitude, longitude, stratum
│   │   ├── core_samples.csv    # core_id, depth_cm, soc_g_kg, bulk_density_g_cm3
│   │   └── CorePoints_Covariates_BC_Canada.csv  # GEE covariates at global cores
│   ├── data_global/            # Global coastal wetland database
│   │   ├── combined_layers_filtered.csv   # WoSIS + Janousek layers (~124K rows)
│   │   └── JANOUSEK_DATA/                 # Raw Janousek files
│   └── covariates/
│       └── BlueCarbon_Covariate_Snapshot_25m_2020_2023.tif  # Local GEE raster
└── outputs/
    ├── rf/                     # RF prediction rasters
    └── transfer/               # TL prediction rasters
```

---

## Pipeline steps

### Main pipeline (`tar_make()`)

| Step | Key targets | What it produces |
|------|-------------|------------------|
| 1 — Data prep | `cores_raw`, `eda_plots`, `cores_harmonized` | VM0033-depth harmonized field cores |
| 2 — Simple extrapolation | `step2_extrapolation` | Per-stratum carbon densities (no raster needed) |
| 3 — Random forest | `rf_data`, `rf_models`, `rf_rasters`, `rf_maps` | Spatial carbon stock maps, variable importance |
| Reports | `report_nonspatial`, `report_rf` | HTML reports in `reports/` |

### Transfer learning pipeline (`tar_make(script="_targets_transfer.R", store="_targets_transfer")`)

| Step | Key targets | What it produces |
|------|-------------|------------------|
| Global harmonization | `global_harmonized` | VM0033-depth harmonized global cores |
| TL data prep | `tl_data` | Global + local data joined with GEE covariates |
| TL models | `tl_models` | Per-depth Wadoux RF + bias correction + bootstrap CI |
| TL rasters | `tl_rasters` | 4-band GeoTIFFs per depth |
| TL maps | `tl_maps` | Comparison maps + LOCO CV validation plots |
| Report | `report_tl` | HTML report |

**Important:** run the main pipeline before the transfer pipeline. The transfer
pipeline reads `cores_harmonized` from the main store via
`tar_read(cores_harmonized, store = "_targets")`.

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
- Exception: output functions (`predict_*_rasters()`) may write to `outputs/` and return the object
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
tar_read(target_name)                              # main store
tar_read(target_name, store = "_targets_transfer") # transfer store
```

---

## Transfer learning specifics

**Bridge variables (reduced — used when N_local < 15):**
`NDVI_median`, `LSWI_median`, `mNDWI_median`, `VV_mean`, `elevation_m`, `elevationRelMHW`

**Bridge variables (full — used when N_local ≥ 15):**
Above plus `slope`, `twi`, `dist_to_channel_m`, `tidal_flat_prob`, `coastal_dist_m`,
`VH_mean`, `VVVH_ratio`, `B`, `G`, `R`, `NIR`, `SWIR1`, `SWIR2`, `SAVI_median`

**Four-stage training per depth:**
1. Wadoux domain classifier (RF, probability forest) → instance weights
2. Weighted global RF (1000 trees, `ranger`)
3. Mean bias correction from local residuals
4. 500-replicate bootstrap → bias SE

**Output raster bands** (`tl_carbon_stocks_kg_m2.tif`):
`dX_Global_Prior`, `dX_Transfer_Final`, `dX_Local_Only`, `dX_Difference`
where X = depth midpoint with `.` replaced by `_` (e.g. `d7_5`, `d22_5`)

**Global data compound key:** `paste(dataset, profile_id, sep = "_")`
e.g. `"Janousek_1"`, `"WOSIS 2023_1144105"`

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
| Step 4 — Transfer learning | Built, not yet tested | Awaiting run with real data |
| Step 4b — TL embeddings (Model 2) | Not started | Add after Step 4 is validated |

---

## What to build next

- Test Step 4 with real covariate data; verify raster band names match bridge vars
- Add embedding-similarity model (Model 2) to transfer pipeline once Step 4 passes
- Consider a combined `tar_make_future()` setup for parallel depth processing if runtime is slow
