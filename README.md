# Blue Carbon Analysis — Spatial Mapping and VM0033 Reporting

**Language:** R · **Orchestration:** `targets` · **Updated:** May 2026

---

## Table of Contents

1. [What Is Blue Carbon?](#what-is-blue-carbon)
2. [What This Workflow Does](#what-this-workflow-does)
3. [Scientific Methods at a Glance](#scientific-methods-at-a-glance)
4. [Getting Started](#getting-started)
   - [Install packages](#1-install-packages)
   - [Clone into RStudio](#2-clone-into-rstudio)
   - [Place your field data](#3-place-your-field-data)
   - [Configure your site](#4-configure-your-site)
   - [Authenticate Google Earth Engine](#5-authenticate-google-earth-engine)
5. [Pipeline 1 — Core Data Preparation](#pipeline-1--core-data-preparation-tar_make)
6. [Pipeline 2 — Global Covariate Extraction (Pre-Analysis)](#pipeline-2--global-covariate-extraction-pre-analysis)
7. [Pipeline 3 — Wadoux Transfer Learning](#pipeline-3--wadoux-transfer-learning)
8. [Pipeline 4 — Embedding Transfer Learning](#pipeline-4--embedding-transfer-learning-model-2)
9. [Understanding the Outputs](#understanding-the-outputs)
10. [How to Interpret Your Results](#how-to-interpret-your-results)
11. [VM0033 Compliance Notes](#vm0033-compliance-notes)
12. [Tips and Troubleshooting](#tips-and-troubleshooting)
13. [Scientific References](#scientific-references)

---

## What Is Blue Carbon?

**Blue carbon** is the carbon stored in the soils and biomass of coastal ecosystems —
tidal marshes, mangroves, and seagrasses. These habitats are exceptional carbon sinks:
they bury organic matter at rates 3–5× higher than terrestrial forests and can store
carbon for thousands of years.

When coastal wetlands are degraded or destroyed, this stored carbon is rapidly
oxidised and released as CO₂. Conversely, protecting or restoring these habitats
can generate verified, durable carbon credits under standards such as **Verra VM0033**
("Methodology for Tidal Wetland and Seagrass Restoration").

Quantifying exactly *how much* carbon is in the soil requires:

1. **Field sampling** — sediment cores extracted at representative locations across
   the site, with soil organic carbon (SOC) and bulk density (BD) measured in the lab
2. **Depth harmonization** — converting irregular sampling intervals to standard
   VM0033 depth layers (0–15, 15–30, 30–50, 50–100 cm)
3. **Spatial modeling** — using remote sensing covariates to interpolate between
   cores and produce a continuous stock map with uncertainty bounds

This workflow automates steps 2 and 3, and optionally leverages global databases
of thousands of coastal wetland cores to improve predictions at sites with limited
local sampling (transfer learning).

---

## What This Workflow Does

The analysis is organized as four reproducible pipelines, each with its own
[`targets`](https://docs.ropensci.org/targets/) store so they can be run and
updated independently.

| Pipeline | Script | Store | What it produces |
|----------|--------|-------|-----------------|
| **1. Core prep** | `_targets.R` | `_targets/` | QA-filtered cores → VM0033-depth harmonized stocks |
| **2. Pre-analysis** | `_targets_preanalysis.R` | `_targets_preanalysis/` | GEE covariates at 952 global Janousek cores |
| **3. Wadoux TL** | `_targets_transfer.R` | `_targets_transfer/` | RF-weighted transfer learning prediction maps |
| **4. Embedding TL** | `_targets_embedding.R` | `_targets_embedding/` | Cosine-similarity-weighted prediction maps (Model 2) |

**Run order:** Pipeline 1 → Pipeline 2 → Pipeline 3 and/or 4 (3 and 4 are independent
once 1 and 2 are complete).

All reports are rendered as self-contained HTML files in `reports/`.

---

## Scientific Methods at a Glance

### Carbon stock calculation

Each core sample contributes a carbon stock at each depth interval:

```
carbon_stock (kg C / m²) = SOC (g/kg) × BD (g/cm³) × thickness (cm) / 100
```

Where `thickness` is the VM0033 interval thickness (15, 15, 20, or 50 cm).
When bulk density was not measured in the lab, ecosystem-specific defaults are used:
Intertidal Marsh (IM) = 0.8 g/cm³, Near-Marsh (NM) = 0.8 g/cm³, Mudflat (MF) = 0.8 g/cm³.

### Depth harmonization

Field cores are sampled at irregular intervals (e.g., 0–10, 10–20 cm). VM0033 requires
standardized intervals (0–15, 15–30, 30–50, 50–100 cm). We use:

- **Equal-area splines** (Bishop et al. 1999) for cores with ≥ 3 samples — fits a
  continuous SOC–depth curve and re-integrates over the target interval, preserving
  total mass exactly
- **Linear interpolation** for 2-sample profiles
- **Exponential decay** to extrapolate below the deepest measured sample

### Random forest spatial prediction

For the local site, a random forest model (1000 trees, `ranger`) is trained on
field cores with 26 remote sensing covariates extracted from Google Earth Engine:
Sentinel-2 spectral bands, SAR backscatter (Sentinel-1), topographic indices
(elevation, slope, TWI, channel distance), and water indices (NDVI, LSWI, mNDWI,
SAVI). The model predicts carbon stock at each 25-metre pixel.

### Transfer learning (Wadoux method — Pipeline 3)

With only 6 local cores, a purely local model is severely data-limited. Transfer
learning supplements the local data with ~952 globally distributed wetland cores
from the Janousek et al. coastal carbon database, each described by the same 26
covariates.

The key question is: *which global cores are most similar to your local site?*
The Wadoux method trains a random forest probability classifier to distinguish
local from global samples, and converts the predicted probability into an
**instance weight** for each global core. Cores that look like your site get high
weights; dissimilar ones are down-weighted. The weighted global RF is then corrected
for local bias using the local residuals.

### Transfer learning (Embedding method — Pipeline 4)

Pipeline 4 uses Google's `GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL` foundation model —
a 64-dimensional per-pixel embedding derived from Sentinel-2 time series, encoding
spectral, textural, and phenological information learned from billions of pixels.

Instead of a hand-crafted similarity classifier, we compute the **cosine similarity**
between each global core's embedding vector and the mean embedding of the local AOI.
This single number captures the holistic "look" of each site in 64-dimensional
learned feature space. Similarity is raised to the power α = 5 (sharpening exponent)
so that only the closest analogues contribute substantially to the weighted model.

> **Key advantage:** No local labels are needed for Stage A — similarity is a
> direct geometric comparison, not a trained classifier.

---

## Getting Started

### 1. Install packages

Run this once in a fresh R session (R ≥ 4.3 recommended):

```r
install.packages(c(
  "targets", "tarchetypes", "geotargets",
  "dplyr", "readr", "tidyr", "ggplot2", "sf",
  "terra", "ranger", "randomForest",
  "quarto", "visNetwork"
))
```

For the transfer learning pipelines you also need `rgee`:

```r
install.packages("rgee")
# Python dependency (once per machine):
rgee::ee_install()
```

### 2. Clone into RStudio

**Via RStudio GUI:**
`File` → `New Project` → `Version Control` → `Git` → paste repository URL.

**Via R console:**
```r
usethis::create_from_github(
  "cathald/bluecarbonanalysis_targets",
  destdir = "~/projects"
)
```

The `.Rprofile` at the repo root loads automatically and gives you four console
shortcuts (see [Tips and Troubleshooting](#tips-and-troubleshooting)).

### 3. Place your field data

Put two CSVs in `Pre-Analysis Data Preparation/data_raw/`:

**`core_locations.csv`** — one row per core:

| Column | Type | Description |
|--------|------|-------------|
| `core_id` | character | Unique core identifier |
| `latitude` | numeric | Decimal degrees (WGS84) |
| `longitude` | numeric | Decimal degrees (WGS84) |
| `stratum` | character | Ecosystem stratum (must match `VALID_STRATA` in config) |

**`core_samples.csv`** — one row per depth sample:

| Column | Type | Description |
|--------|------|-------------|
| `core_id` | character | Links to core_locations.csv |
| `depth_top_cm` | numeric | Top of sampled interval (cm) |
| `depth_bottom_cm` | numeric | Bottom of sampled interval (cm) |
| `soc_g_kg` | numeric | Soil organic carbon (g C per kg dry soil) |
| `bulk_density_g_cm3` | numeric | Bulk density (g/cm³); can be `NA` if not measured |

You also need a **local covariate raster** — a multi-band GeoTIFF covering the
full site extent, exported from Google Earth Engine. Its path is set in
`blue_carbon_config.R` under `COVARIATE_RASTER`. This file should contain
the 26 canonical bands listed in `R/preanalysis/gee_covariates.R`.

### 4. Configure your site

Open `blue_carbon_config.R` at the project root and set:

```r
PROJECT_NAME    <- "YourSite_2026"
PROJECT_LOCATION <- "e.g. Chemainus Estuary, British Columbia, Canada"
VALID_STRATA    <- c("IM", "NM", "MF")   # must match stratum column in CSV
COVARIATE_RASTER <- "Pre-Analysis Data Preparation/covariates/your_raster.tif"
```

The VM0033 depth intervals, BD defaults, and reporting thresholds are pre-configured
and should only be changed after consulting the methodology documentation.

### 5. Authenticate Google Earth Engine

Required for Pipelines 2, 3, and 4. Run once per machine:

```r
library(rgee)
ee_Initialize(user = "your.email@gmail.com", drive = TRUE)
```

This opens a browser for OAuth. After authorising, credentials are cached and
you will not need to repeat this step. For the embedding pipeline (Pipeline 4),
Drive access is **not** required — embeddings are downloaded directly as GeoTIFFs
over HTTPS.

---

## Pipeline 1 — Core Data Preparation (`tar_make()`)

This is the foundation of the entire workflow. It must run to completion before
any other pipeline.

**What it does:**

1. Reads and merges `core_locations.csv` and `core_samples.csv`
2. Applies QA filters: flags missing bulk density, implausible SOC, implausible BD,
   and strata not in `VALID_STRATA`; fills missing BD with ecosystem defaults
3. Calculates carbon stocks at each measured depth
4. Harmonizes all cores to the four VM0033 depth intervals using equal-area splines
5. Runs exploratory data analysis (depth profiles, stratum comparisons, outlier plots)
6. Renders an HTML report (`reports/step1_nonspatial.html`)

**Run it:**

```r
targets::tar_make()
# or use the console shortcut:
tm()
```

**Check progress:**

```r
tv()               # visual dependency graph (green = done, orange = stale)
tar_meta() |> dplyr::select(name, seconds, error)   # per-target timing and errors
```

**Key targets you can inspect:**

```r
tar_read("cores_raw")         # raw field data after merging
tar_read("eda_plots")         # list of 6 exploratory plots
tar_read("cores_harmonized")  # harmonized carbon stocks at VM0033 depths
```

After this pipeline completes, `_targets/` contains a cached, validated version
of your field data. Any subsequent change to the input CSVs or config will
automatically mark downstream targets as stale.

---

## Pipeline 2 — Global Covariate Extraction (Pre-Analysis)

This pipeline extracts the same 26-band remote sensing covariates at each of the
~952 global Janousek wetland core locations. It runs once and its outputs are used
by both transfer learning pipelines.

**What it does:**

1. Reads global core lat/lon coordinates from the Janousek database
2. Sends batched requests to Google Earth Engine to extract:
   - Sentinel-2 spectral bands (B, G, R, NIR, SWIR1, SWIR2, NDVI, LSWI, mNDWI, SAVI)
   - Sentinel-1 SAR (VV, VH, VV/VH ratio)
   - Topography (elevation, slope, TWI, distance to channel, tidal flat probability,
     coastal distance, elevation relative to MHW)
   - Climate (MAP, MAT)
3. Combines all bands, validates against the canonical 26-band schema, and writes
   `Pre-Analysis Data Preparation/data_raw/CorePoints_Covariates_BC_Canada.csv`

**Run it:**

```r
targets::tar_make(
  script = "_targets_preanalysis.R",
  store  = "_targets_preanalysis"
)
```

> **Note on runtime:** GEE extraction for ~952 cores takes 30–90 minutes depending
> on network conditions. The pipeline uses batches of 5 cores with a 5-minute
> per-batch timeout and automatic retry. If it stops partway through, re-running
> `tar_make()` skips any batch that already succeeded.

**You only need to run this once.** The output CSV is committed to the repository;
if you have not modified the Janousek coordinate data or the band list, this
pipeline will report that all targets are up to date.

---

## Pipeline 3 — Wadoux Transfer Learning

This pipeline uses the Janousek global database, weighted by spectral similarity
to your local site, to improve carbon stock predictions at each VM0033 depth.

**What it does:**

1. Harmonizes the Janousek database to VM0033 depth intervals (same spline method
   as Pipeline 1)
2. Joins global and local cores with their GEE covariates (bridge variables)
3. **Stage A:** Trains a random forest probability classifier to separate local from
   global samples → converts predicted probabilities to instance weights
4. **Stage B:** Trains a weighted global RF (1000 trees) at each depth midpoint
5. **Stage C:** Computes mean bias from local residuals; estimates uncertainty via
   500-replicate bootstrap
6. **Stage D:** Applies bias correction to global RF predictions and predicts across
   the local site raster
7. Runs leave-one-core-out (LOCO) cross-validation to assess transferability
8. Produces comparison maps and a similarity heatmap
9. Renders `reports/step4_transfer_learning.html`

**Run it:**

```r
targets::tar_make(
  script = "_targets_transfer.R",
  store  = "_targets_transfer"
)
```

**Key targets:**

```r
targets::tar_read("tl_models",  store = "_targets_transfer")  # per-depth model objects
targets::tar_read("tl_maps",    store = "_targets_transfer")  # list of ggplot objects
targets::tar_read("tl_rasters", store = "_targets_transfer")  # SpatRaster (4 bands × 4 depths)
```

---

## Pipeline 4 — Embedding Transfer Learning (Model 2)

This is an alternative transfer learning approach that replaces the Wadoux
RF classifier with **embedding cosine similarity** from Google's foundation
model satellite embeddings. It runs independently from Pipeline 3.

**What it does:**

1. **Phase 1 — GEE extraction:**
   - Extracts a 64-dimensional embedding vector at each of the ~952 global core
     locations (`GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL`, averaged 2023–2025)
   - Downloads the same 64-band embedding raster over the local AOI extent
     (split into 4 chunks of 16 bands to stay under GEE's 48 MB download limit)
   - Writes `outputs/embedding/aoi_embedding_raster.tif`

2. **Phase 2 — Cosine similarity weights:**
   - Computes the mean embedding vector across all pixels in the local AOI
   - Computes cosine similarity between each global core's embedding and the AOI mean
   - Applies sharpening exponent α = 5: `weight = sim^5`, normalized so mean = 1

3. **Phase 3 — Weighted transfer learning:**
   - Same four-stage model as Pipeline 3, but instance weights come from step 2
     instead of a Wadoux RF classifier
   - Produces identical output format: 4-band GeoTIFFs per depth interval
   - Renders `reports/step5_embedding_tl.html`, including a side-by-side comparison
     of Pipeline 3 (Wadoux) vs Pipeline 4 (Embedding) LOCO CV performance

**Run it:**

```r
targets::tar_make(
  script = "_targets_embedding.R",
  store  = "_targets_embedding"
)
```

> **Prerequisites:** Pipeline 1 must be complete. The GEE project ID is set at the
> top of `_targets_embedding.R` — update `GEE_PROJECT` to your own project if needed.

---

## Understanding the Outputs

### Reports (HTML)

| File | From | Contents |
|------|------|----------|
| `reports/step1_nonspatial.html` | Pipeline 1 | QA flags, depth profiles, stratum summaries, harmonized stocks |
| `reports/step3_random_forest.html` | Pipeline 1 | Local RF variable importance, CV, prediction map |
| `reports/step4_transfer_learning.html` | Pipeline 3 | Wadoux weights, LOCO CV, spatial maps, similarity heatmap |
| `reports/step5_embedding_tl.html` | Pipeline 4 | Cosine similarity, LOCO CV, spatial maps, Step 4 vs 5 comparison |

### Prediction rasters (GeoTIFF)

Both transfer learning pipelines write one GeoTIFF per depth to `outputs/transfer/`.
Each file has **four bands**:

| Band name | Description |
|-----------|-------------|
| `dX_Global_Prior` | Weighted global RF prediction before bias correction |
| `dX_Transfer_Final` | Global prior + local bias correction → **use this for VM0033** |
| `dX_Local_Only` | Mean local carbon stock (spatial constant — no spatial model) |
| `dX_Difference` | Transfer Final minus Global Prior (magnitude of local correction) |

Where `X` is the depth midpoint with `.` replaced by `_` (e.g. `d7_5`, `d22_5`, `d40`, `d75`).

### Embedding raster

`outputs/embedding/aoi_embedding_raster.tif` — 64-band GeoTIFF at the local AOI
extent, with bands named `emb_1` … `emb_64`. Used internally for cosine similarity
computation; can also be used for unsupervised habitat clustering.

### Model objects (in `_targets` stores)

You can load any target directly from R:

```r
# Local RF models (Pipeline 1):
rf <- targets::tar_read("rf_models")           # list of per-depth ranger objects
rf$`7.5`$model                                 # ranger object for 0–15 cm

# Transfer learning models (Pipeline 3):
m  <- targets::tar_read("tl_models", store = "_targets_transfer")
m$models[["7.5"]]$bias_correction              # bias correction term in kg/m²
m$models[["7.5"]]$bias_se                      # bootstrap SE of bias
m$models[["7.5"]]$cv_r2_tl                     # LOCO CV R² (transfer-corrected)
m$models[["7.5"]]$cv_r2_global                 # LOCO CV R² (global prior only)
```

---

## How to Interpret Your Results

### LOCO cross-validation

The **leave-one-core-out (LOCO) CV** is the most important quality metric for
the transfer learning models. At each fold, one local core is held out and the
bias correction is estimated from the remaining local cores. This simulates
"would this model have performed well at a core we hadn't yet sampled?"

- **R² > 0.5** indicates good transferability at that depth
- **RMSE** is in kg C / m² — compare against the local mean carbon stock to
  interpret as a fraction
- **Transfer Final vs Global Prior:** if LOCO CV RMSE is lower for Transfer Final
  than Global Prior, bias correction is helping; if not, the global prior may already
  be well-calibrated for your site

### Cosine similarity (Pipeline 4)

Values range 0–1. Cores with **similarity > 0.9** are very close analogues to
your site in satellite appearance; cores below 0.7 contribute very little weight
after the α = 5 sharpening. The weight distribution plot in Step 5 shows how
concentrated or diffuse the weighting is.

### Choosing between Pipeline 3 and Pipeline 4

Compare the `cv_rmse_tl` and `cv_r2_tl` columns in the Step 5 HTML report's
comparison table:

- If embedding TL has **lower RMSE** → use Pipeline 4 outputs for VM0033 reporting
- If Wadoux TL has **lower RMSE** → use Pipeline 3 outputs
- If very similar → either is defensible; embeddings have the advantage of not
  requiring a local classifier, which is more robust with N_local < 10

### Prediction intervals for VM0033

VM0033 requires **90% prediction intervals** (not 95%). Compute them from any
depth model object:

```r
m   <- targets::tar_read("tl_models", store = "_targets_transfer")$models[["7.5"]]
se  <- sqrt(m$bias_se^2 + m$residual_sd^2)
z90 <- qnorm(0.95)   # one-sided 95% = two-sided 90% interval
# For a pixel with transfer_final prediction p:
lower_90 <- p - z90 * se
upper_90 <- p + z90 * se
```

Use the **lower bound** (`lower_90`) as the conservative VM0033 reporting value.

---

## VM0033 Compliance Notes

| Requirement | How this workflow meets it |
|-------------|--------------------------|
| Standard depth intervals | Hardcoded in `blue_carbon_config.R`: 0–15, 15–30, 30–50, 50–100 cm |
| Mass-preserving harmonization | Equal-area splines (Bishop et al. 1999) in `R/depth_harmonization.R` |
| Bulk density defaults | Ecosystem-specific defaults from `cfg$BD_DEFAULTS` applied only when measured BD is missing |
| Carbon stock formula | `SOC × BD × thickness / 100` in `R/qc.R` |
| 90% prediction intervals | Bias SE + residual SD propagated; `qnorm(0.95)` multiplier |
| Leave-one-core-out CV | Implemented in `train_tl()` and `train_emb_tl()` — not leave-one-observation-out |
| Conservative reporting | Lower 90% PI is the reporting value, not the mean prediction |

> **Do not change** `VM0033_DEPTH_MIDPOINTS` or the carbon stock formula without
> consulting your VM0033 methodology documentation and updating the verification package.

---

## Tips and Troubleshooting

### Console shortcuts (loaded from `.Rprofile`)

```r
tm()    # targets::tar_make()                   — run all outdated targets
tv()    # targets::tar_visnetwork()              — visual dependency graph
tl(x)  # targets::tar_load("x")               — load target x into session
tm1()  # run only Step 1 targets by name       — faster iteration during development
```

### Useful `targets` commands

```r
tar_outdated()                        # what will re-run and why
tar_make(names = "cores_harmonized")  # run one target + its dependencies
tar_meta() |> dplyr::select(name, seconds, error)  # timing and error status
tar_invalidate("global_embeddings")   # force one target to re-run next time
tar_destroy(); tar_make()             # clear entire cache and start fresh
```

For the non-main stores, add `script` and `store` arguments:

```r
targets::tar_make(
  names  = "tl_models",
  script = "_targets_transfer.R",
  store  = "_targets_transfer"
)
```

### GEE extraction is slow or hangs

- Each S2 batch uses a 5-minute timeout; a hang causes a catchable error and
  the batch is retried automatically on the next `tar_make()` call
- If a particular batch consistently fails, check `tar_meta()` for the error
  message and inspect the GEE JavaScript console for quota or collection issues
- The S2 extraction runs in batches of 5 cores — reducing this further is possible
  by editing `.BATCH_SIZE` in `R/preanalysis/gee_covariates.R`

### Raster band names look wrong after loading from store

`geotargets` serialises `SpatRaster` objects as GeoTIFF; band names can be
mangled during the round-trip. The `plot_tl_maps()` function reconstructs the
expected band names from the model object structure before indexing, so this
is handled automatically in the reports. If you load a raster manually and the
band names look wrong, run:

```r
r <- targets::tar_read("tl_rasters", store = "_targets_transfer")
# Check: names(r) should contain "d7_5_Global_Prior" etc.
# If not, reconstruct:
band_sfx <- c("Global_Prior", "Transfer_Final", "Local_Only", "Difference")
depths   <- c("7_5", "22_5", "40", "75")
names(r) <- paste0("d", rep(depths, each = 4), "_", band_sfx)
```

### "No embedding values returned" error in Pipeline 4

This usually means the `GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL` collection does
not have imagery for all requested years over some core locations. Check that
`EMB_YEARS` in `_targets_embedding.R` covers years where the collection has
global coverage (2017–present as of 2026).

### Changing the GEE project

All three GEE pipelines default to the project ID in their respective scripts:
- Pipeline 2: `GEE_PROJECT` at top of `_targets_preanalysis.R`
- Pipeline 3: same variable in `_targets_transfer.R`
- Pipeline 4: same variable in `_targets_embedding.R`

Update all three if moving to a different GCP project.

---

## Scientific References

**Depth harmonization**
- Bishop, T.F.A., McBratney, A.B., & Laslett, G.M. (1999). Modelling soil attribute depth functions with equal-area quadratic smoothing splines. *Geoderma*, 91, 27–45. https://doi.org/10.1016/S0016-7061(99)00003-8
- Malone, B.P. et al. (2009). Mapping continuous depth functions of soil carbon storage and available water capacity. *Geoderma*, 154, 138–152. https://doi.org/10.1016/j.geoderma.2009.10.007

**Blue carbon methodology and VM0033**
- Howard, J. et al. (2014). *Coastal Blue Carbon: Methods for Assessing Carbon Stocks and Emissions Factors in Mangroves, Tidal Salt Marshes, and Seagrass Meadows*. Conservation International / IOC-UNESCO / IUCN.
- IPCC (2014). *2013 Supplement to the 2006 IPCC Guidelines: Wetlands*. IPCC, Switzerland.

**Global coastal carbon database**
- Janousek, C.N. et al. (2025). Coastal wetland carbon stocks database. [Janousek et al. global synthesis]
- Holmquist, J.R. et al. (2018). Accuracy and Precision of Tidal Wetland Soil Carbon Mapping in the Conterminous United States. *Scientific Reports*, 8, 9478. https://doi.org/10.1038/s41598-018-26948-7

**Transfer learning and instance weighting**
- Wadoux, A.M.J.-C., Samuel-Rosa, A., Poggio, L., & Mulder, V.L. (2021). A note on knowledge discovery and machine learning in digital soil mapping. *European Journal of Soil Science*, 71, 133–136. https://doi.org/10.1111/ejss.12909

**Random forest**
- Breiman, L. (2001). Random Forests. *Machine Learning*, 45, 5–32. https://doi.org/10.1023/A:1010933404324
- Wright, M.N. & Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. *Journal of Statistical Software*, 77(1). https://doi.org/10.18637/jss.v077.i01

**Satellite embeddings**
- Google (2025). GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL — 64-band foundation model embeddings derived from Sentinel-2 annual composites. Google Earth Engine Data Catalog.

**Targets pipeline framework**
- Landau, W.M. (2021). The targets R package: a dynamic Make-like function-oriented pipeline toolkit for reproducibility and high-performance computing. *Journal of Open Source Software*, 6(57), 2959. https://doi.org/10.21105/joss.02959
