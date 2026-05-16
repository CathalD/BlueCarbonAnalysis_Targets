# Blue Carbon analysis, mapping and reporting workflow
 
**Date:** January, 2026  
**Language:** R  
**Platform:** RStudio  

-----

## Table of Contents

1. [Overview](#overview)
2. [Quick Guide](#quick-guide)
3. [Required Inputs](#required-inputs)
4. [Setting Up Your Project](#setting-up-your-project)
5. [Expected Outputs](#expected-outputs)
6. [Implementation Guide](#implementation-guide)
7. [References](#references)
8. [Appendix A - Scientific Foundation](#appendix-a---scientific-foundation)

---

## Overview

This workflow automates the harmonization, analysis, mapping, and reporting of carbon stocks from Coastal Blue carbon sediments. The system processes field sediment core data through spatial modeling to produce carbon stock estimates with quantified uncertainty.

### Analysis Types

There are 4 types of spatial analysis this workflow can perform:

**A) Basic Reporting** - Ingests individual core data and spatial boundaries, automatically calculates sample carbon data, harmonizes to defined depth intervals, averages carbon stocks across spatial boundaries, and runs spatial interpolation (kriging) to produce a carbon map.

**B) Spatial Extrapolation with Remote Sensing** - Advanced function requiring remote sensing covariates. Runs random-forest machine-learning algorithm to build spatial models reflecting smaller-scale variations across the project boundary.

**C) Bayesian-Based Analysis** - Can be run with A) or B). Requires prior carbon stock data within your project boundary. Updates prior probabilities based on new core data to produce refined carbon maps.

**D) Transfer Learning Analysis** - Can be used with A, B, or C. Uses previous core data to build predictive models by weighting similarity between your project area and existing datasets.

---

## Quick Guide

**Step 1** - Open R/RStudio

**Step 2 (Option 1)** - Download "MasterScript" and import into RStudio, run this script to pull documents directly from GitHub via URL

**Step 2 (Option 2)** - Download entire GitHub repo - Use "MasterScript_Offline" to run workflow

**Step 3** - Add your carbon data to folder "data_raw" in the same format as `core_locations_TEMPLATE.csv` and `core_samples_TEMPLATE.csv`

**Step 4** - Run "Modules P1 and P2" for Basic Analysis

### For Advanced Analysis with Remote Sensing:

**Step 5** - Add covariate files to `BlueCarbon_Workflow_V1.0/covariates` - these should be clipped to the project boundary

**Step 6** - Run Modules P1 and P2

### For Bayesian Analysis:

**Step 7** - Add prior carbon maps to folder `BlueCarbon_Workflow_V1.0/data_prior`

**Step 8** - In `blue_carbon_config.R` document, scroll to section "Bayesian" and change `Bayesian <- FALSE` to `Bayesian <- TRUE`

**Step 9** - Run modules P1, P2 and P3

### For Transfer Learning Analysis:

**Step 10** - Add covariate files to `BlueCarbon_Workflow_V1.0/covariates` - these should be clipped to the project boundary

**Step 11** - Ensure the covariates listed in "covariates" match those used in `global_cores_with_gee_covariates.csv` file

**Step 12** - Run P2 and P4 Modules

---

## Required Inputs

### A) Basic Reporting

1. **Carbon stock core locational data** as a spreadsheet in this format:

<img width="864" height="219" alt="Screenshot 2026-02-04 at 12 44 08 PM" src="https://github.com/user-attachments/assets/7e469415-67b8-41da-ab91-b4ce0592c757" />

2. **Carbon stock sample lab data** as a spreadsheet in this format:

<img width="751" height="317" alt="Screenshot 2026-02-04 at 12 45 17 PM" src="https://github.com/user-attachments/assets/caa4a89e-86cb-4530-a9f5-988d41c42425" />

3. **Project boundary and (optional) strata** in one of these formats: GeoJson, CSV, .shp, etc.

### B) Spatial Extrapolation with Remote Sensing Models

*All of above from "Basic reporting" +*

1. **Remote sensing covariates**
   - Can use any custom set of RS variables as desired
   - Option to follow Google Earth Engine workflow that automatically extracts covariates

### C) Bayesian-Based Analysis

*All from A) Basic Reporting (optional: can also handle RS input)*

1. **Prior carbon map** in a geoTiff (.tiff) format

### D) Transfer Learning Analysis

*All from B) Spatial extrapolation with remote sensing models*

1. Uses the global dataset from Janousek et al. database
   - **Option 1:** Use pre-loaded covariates from this script (Link: GEE Python API script)
   - **Option 2:** Add personal covariates (ensure local and global covariates align)
   - **Option 3:** Uses Google DeepMind's embedded learning to build a "similarity" matrix for local cores and global dataset

---

## Setting Up Your Project

The file `blue_carbon_config.R` is where you configure the workflow settings.

### 1. Change the Project Metadata

<details>
<summary><b>Click to view configuration code</b></summary>

```r
# =================
# PROJECT METADATA
# =================
PROJECT_NAME <- "BC_Coastal_BlueCarbon_2024"
PROJECT_SCENARIO <- "PROJECT"  # Options: BASELINE, PROJECT, CONTROL, DEGRADED
MONITORING_YEAR <- 2025

# Project location (for documentation)
PROJECT_LOCATION <- "Chemainus Estuary, British Columbia, Canada"
PROJECT_DESCRIPTION <- "Blue carbon monitoring to report to funder"
```

</details>

### 2. Define Your Stratified Boundaries

<details>
<summary><b>Click to view stratification configuration code</b></summary>

```r
# Valid ecosystem strata (must match GEE stratification tool)
# FILE NAMING CONVENTION:
#   Module 05 auto-detects GEE stratum masks using this pattern:
#   "Stratum Name" → stratum_name.tif in data_raw/gee_strata/
# Examples:
#   "Upper Marsh"           → upper_marsh.tif
#   "Underwater Vegetation" → underwater_vegetation.tif
#   "Emerging Marsh"        → emerging_marsh.tif

# CUSTOMIZATION OPTIONS:
#   1. Simple: Edit VALID_STRATA below and export GEE masks with matching names
#   2. Advanced: Create stratum_definitions.csv in project root for custom file names
#      and optional metadata (see stratum_definitions_EXAMPLE.csv template)

VALID_STRATA <- c(
  "Mid Marsh",         
  "Upper Marsh",             
  "Lower Marsh",
  "Underwater vegetation"
)

# Stratum colors for plotting (match GEE tool)
STRATUM_COLORS <- c(
  "Mid Marsh" = "#FFFF99",
  "Upper Marsh" = "#99FF99",
  "Lower Marsh" = "#33CC33",
  "Underwater vegetation" = "#FF7F50"
)
```

</details>

### 3. Configure Depth Harmonization

For data harmonization, define the depth interval midpoint values (default set to Verra tidal wetland standard depths).

<details>
<summary><b>Click to view depth configuration code</b></summary>

```r
# ============================================================================
# DEPTH CONFIGURATION
# ============================================================================
# VM0033 standard depth intervals (cm) - depth midpoints for harmonization
# These correspond to VM0033 depth layers: 0-15, 15-30, 30-50, 50-100 cm
VM0033_DEPTH_MIDPOINTS <- c(7.5, 22.5, 40, 75)

# VM0033 depth intervals (cm) - for mass-weighted aggregation
VM0033_DEPTH_INTERVALS <- data.frame(
  depth_top = c(0, 15, 30, 50),
  depth_bottom = c(15, 30, 50, 100),
  depth_midpoint = c(7.5, 22.5, 40, 75),
  thickness_cm = c(15, 15, 20, 50)
)

# Standard depths for harmonization (VM0033 midpoints are default)
STANDARD_DEPTHS <- VM0033_DEPTH_MIDPOINTS

# Fine-scale depth intervals (optional, for detailed analysis)
FINE_SCALE_DEPTHS <- c(0, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100)

# Maximum core depth (cm)
MAX_CORE_DEPTH <- 100

# Key depth intervals for reporting (cm)
REPORTING_DEPTHS <- list(
  surface = c(0, 30),      # Top 30 cm (most active layer)
  subsurface = c(30, 100)  # 30-100 cm (long-term storage)
)
```

</details>

---

## Expected Outputs

### A) Basic Reporting

**Module P2_02: Exploratory Data Analysis**
- Generates summary statistics by stratum
- Creates depth profile visualizations
- Identifies outliers and data quality issues
- Produces diagnostic plots

**Module P2_03: Depth Harmonization**
- Implements equal-area spline functions
- Standardizes to depth intervals (Default aligns with VM0033 Methods: 7.5, 22.5, 40, 75 cm)
- Extrapolates to maximum depth (100 cm)
- Validates mass balance
- Exports harmonized profiles for spatial modeling

**Module P2_04: Spatial Predictions - Kriging**
- Automated variogram modeling
- Universal kriging predictions
- Cross-validation assessment
- Uncertainty quantification via kriging variance
- Sample size power analysis

### B) Remote Sensing Informed Reporting

**Module P2_05: Spatial Predictions - Random Forest**
- Integration with Google Earth Engine covariates
- Random Forest model training
- Area of Applicability assessment
- Spatial cross-validation
- Variable importance analysis

**Module P2_06: Carbon Stock Calculation**
- Aggregates depth-specific predictions
- Calculates total carbon stocks (0-100 cm)
- Compares Kriging vs. Random Forest estimates
- Generates carbon stock maps
- Exports summary tables by stratum

**Module P2_07: MMRV Reporting**
- Generates VM0033 verification package
- Creates comprehensive HTML reports
- Exports spatial data for GIS verification
- Provides sampling recommendations
- Produces quality assurance documentation

---

## Implementation Guide

*(This section can be expanded with step-by-step instructions for running each module)*

---

## References

*(Add your references here as needed)*

---

## Appendix A - Scientific Foundation

### Depth Harmonization

**Equal-Area Spline Functions** (Bishop et al. 1999; Malone et al. 2009):
- Fits continuous depth functions to irregular sampling intervals
- Preserves mass balance (∫SOC·BD·dz constant)
- Prevents artifacts from interval averaging
- Enables standardized depth reporting for VM0033 compliance

**Hybrid Approach** (implemented in Module P2_03):
- Equal-area splines for profiles with ≥3 samples
- Mass-preserving linear interpolation for 2-sample profiles
- Exponential decay extrapolation beyond deepest sample (Holmquist et al. 2018)

<details>
<summary><b>Click to view references</b></summary>

- Bishop, T.F.A., McBratney, A.B., & Laslett, G.M. (1999). "Modelling soil attribute depth functions with equal-area quadratic smoothing splines." *Geoderma*, 91(1-2), 27-45. https://doi.org/10.1016/S0016-7061(99)00003-8
- Malone, B.P., McBratney, A.B., Minasny, B., & Laslett, G.M. (2009). "Mapping continuous depth functions of soil carbon storage and available water capacity." *Geoderma*, 154(1-2), 138-152. https://doi.org/10.1016/j.geoderma.2009.10.007
- Holmquist, J.R., Windham-Myers, L., Blaauw, M., et al. (2018). "Accuracy and Precision of Tidal Wetland Soil Carbon Mapping in the Conterminous United States." *Scientific Reports*, 8(1), 9478. https://doi.org/10.1038/s41598-018-26948-7

</details>

---

### Spatial Interpolation - Kriging

**Universal Kriging** (Goovaerts 1997; Webster & Oliver 2007):
- Geostatistical method accounting for spatial autocorrelation
- Variogram modeling captures spatial dependency structure
- Best Linear Unbiased Prediction (BLUP) with uncertainty estimates
- Suitable when covariates unavailable or sample size limited

**Implementation** (Module P2_04):
- Automated variogram fitting (exponential, spherical, Gaussian models)
- Cross-validation for model selection
- Prediction standard errors for uncertainty quantification

<details>
<summary><b>Click to view references</b></summary>

- Goovaerts, P. (1997). *Geostatistics for Natural Resources Evaluation*. Oxford University Press.
- Webster, R., & Oliver, M.A. (2007). *Geostatistics for Environmental Scientists*, 2nd Edition. Wiley. https://doi.org/10.1002/9780470517277

</details>

---

### Spatial Interpolation - Random Forest

**Machine Learning Approach** (Breiman 2001; Hengl et al. 2018):
- Ensemble decision tree method using remote sensing covariates
- Handles non-linear relationships and variable interactions
- Feature importance quantification
- Area of Applicability assessment (Meyer & Pebesma 2021)

**Implementation** (Module P2_05):
- Integrates Landsat, Sentinel-1 SAR, topographic data
- Ranger package for efficient random forest
- CAST package for spatial cross-validation and AOA
- Uncertainty via quantile regression forests

<details>
<summary><b>Click to view references</b></summary>

- Breiman, L. (2001). "Random Forests." *Machine Learning*, 45(1), 5-32. https://doi.org/10.1023/A:1010933404324
- Hengl, T., Nussbaum, M., Wright, M.N., Heuvelink, G.B.M., & Gräler, B. (2018). "Random forest as a generic framework for predictive modeling of spatial and spatio-temporal variables." *PeerJ*, 6, e5518. https://doi.org/10.7717/peerj.5518
- Meyer, H., & Pebesma, E. (2021). "Predicting into unknown space? Estimating the area of applicability of spatial prediction models." *Methods in Ecology and Evolution*, 12(9), 1620-1633. https://doi.org/10.1111/2041-210X.13650
- https://soilmapper.org/

</details>

---

### Bayesian Integration (Optional)

**Prior Knowledge Transfer** (Wadoux et al. 2021):
- Combines site-specific likelihood with global prior distribution
- Reduces uncertainty when local samples limited
- Weight calibration via empirical Bayes or cross-validation

**Implementation** (Part 3 modules):
- Global database as informative prior (Holmquist et al. 2018; Macreadie et al. 2019)
- Bayesian updating via conjugate Normal-Normal framework
- Posterior mean and credible intervals for conservative estimates

<details>
<summary><b>Click to view references</b></summary>

- Wadoux, A.M.J.-C., Brus, D.J., & Heuvelink, G.B.M. (2021). "Accounting for non-stationary variance in geostatistical mapping of soil properties." *Geoderma*, 324, 114138. https://doi.org/10.1016/j.geoderma.2018.03.010
- Holmquist, J.R., et al. (2018). "Coastal and Marine Ecological Classification Standard." NOAA Technical Memorandum NOS NCCOS 258.

</details>

---

### Transfer Learning (Optional)

**Global-to-Local Knowledge Transfer** (Module P4_05):
- Instance weighting (Wadoux et al. 2021) to identify relevant global samples
- Hierarchical modeling across depth intervals
- Bias correction for local conditions

**Global Database:**
- Coastal Carbon Network synthesis (Holmquist et al. 2018)
- Smithsonian MarineGEO global cores
- Harmonized to VM0033 depths

<details>
<summary><b>Click to view references</b></summary>

- Wadoux, A.M.J.-C., Samuel-Rosa, A., Poggio, L., & Mulder, V.L. (2021). "A note on knowledge discovery and machine learning in digital soil mapping." *European Journal of Soil Science*, 71(2), 133-136. https://doi.org/10.1111/ejss.12909

</details>

---

### Carbon Stock Calculation

**Methodology** (Howard et al. 2014; IPCC 2014):

<details>
<summary><b>Click to view calculation formula</b></summary>

```
Carbon Stock (kg C/m²) = ∫[SOC (g/kg) × BD (g/cm³) × 10] dz

Where:
  SOC = Soil organic carbon concentration (g C / kg dry soil)
  BD  = Bulk density (g dry soil / cm³)
  dz  = Depth increment (cm)
  10  = Unit conversion factor
```

</details>

**Bulk Density Handling:**
- Measured BD used when available
- Ecosystem-specific defaults (Morris et al. 2016) when missing:
  - Saltmarsh (EM): 0.52 g/cm³
  - Seagrass (SG): 0.89 g/cm³
  - Mangrove (FL): 0.38 g/cm³

<details>
<summary><b>Click to view references</b></summary>

- Howard, J., Hoyt, S., Isensee, K., Pidgeon, E., & Telszewski, M. (2014). *Coastal Blue Carbon: Methods for Assessing Carbon Stocks and Emissions Factors in Mangroves, Tidal Salt Marshes, and Seagrass Meadows*. Conservation International, Intergovernmental Oceanographic Commission of UNESCO, International Union for Conservation of Nature. Arlington, Virginia, USA.
- IPCC (2014). *2013 Supplement to the 2006 IPCC Guidelines for National Greenhouse Gas Inventories: Wetlands*. Hiraishi, T., Krug, T., Tanabe, K., et al. (eds). IPCC, Switzerland.
- Morris, J.T., Barber, D.C., Callaway, J.C., et al. (2016). "Contributions of organic and inorganic matter to sediment volume and accretion in tidal wetlands at steady state." *Earth's Future*, 4(4), 110-121. https://doi.org/10.1002/2015EF000334

</details>

---

### Uncertainty Quantification

**Conservative Approach** (IPCC 2006, VM0033):
- Propagation of measurement, spatial, and model uncertainties
- Lower 95% confidence bound for crediting
- Monte Carlo simulation when analytical propagation infeasible

**Sources of Uncertainty:**

<details>
<summary><b>Click to view uncertainty sources</b></summary>

1. **Measurement:** Lab analytical precision (typically 2-5% for SOC)
2. **Spatial:** Kriging variance or RF quantile spread
3. **Model:** Cross-validation RMSE, variogram uncertainty
4. **Temporal:** Inter-annual variability (if multi-year data)

</details>






