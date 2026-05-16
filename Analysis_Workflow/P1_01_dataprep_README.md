# Module 01: Blue Carbon Data Preparation

## Overview

This module is the foundation of the blue carbon analysis workflow. It ingests raw field core data, performs quality control, validates measurements against statistcal baselines, and prepares harmonized datasets for spatial analysis.

**Primary Goal:** Transform raw sediment core measurements into a clean, validated dataset ready for depth harmonization and spatial modeling. As well, inform the user how there sampling strategy compares to a theoretical optimal.

---

## Table of Contents

1. [Input Data Requirements](#input-data-requirements)
2. [Data Processing Workflow](#data-processing-workflow)
3. [Quality Control & Validation](#quality-control--validation)
4. [VM0033 Compliance Assessment](#vm0033-compliance-assessment)
5. [Bulk Density Handling](#bulk-density-handling)
6. [Carbon Stock Calculation](#carbon-stock-calculation)
7. [Core Type Analysis](#core-type-analysis)
8. [Outputs](#outputs)

---

## Input Data Requirements

### Required Files

**1. Core Locations (`data_raw/core_locations.csv`)**
- GPS coordinates (longitude, latitude)
- Stratum assignments (ecosystem classification)
- Core metadata (core_id, scenario_type, monitoring_year)

**2. Core Samples (`data_raw/core_samples.csv`)**
- Depth intervals (depth_top_cm, depth_bottom_cm)
- Soil organic carbon content (soc_g_kg)
- Bulk density measurements (bulk_density_g_cm3) - optional but recommended

**3. Configuration (`blue_carbon_config.R`)**
- Valid stratum names and colors
- *Change from VM0033 to statistcal baselines or "theoretical optimum sampling"
- Quality control thresholds
- Bulk density defaults by ecosystem (fill in the blanks)

---

## Data Processing Workflow

### 1. Load and Standardize Core Locations

**What it does:**
- Reads GPS coordinates for each core location
- Standardizes column names to lowercase for consistency
- Validates that all required columns are present

**Scientific rationale:** Accurate spatial information is critical for geostatistical modeling and stratified sampling design.

### 2. Coordinate Validation and Quality Control

**What it does:**
- Checks coordinates are within valid lat/lon ranges
- Identifies duplicate locations (cores within 1m GPS precision)
- Removes cores with invalid or missing coordinates
- Generates duplicate location reports for field verification

**Why this matters:** GPS precision issues can affect spatial autocorrelation analysis and kriging predictions. Cores at identical locations may indicate equipment errors or intentional paired sampling designs.

### 3. Stratum Validation

**What it does:**
- Verifies all stratum names match the valid list from configuration
- Flags invalid stratum assignments
- Ensures consistency with Google Earth Engine stratification

**Scientific rationale:** Ecosystem stratification is fundamental to blue carbon accounting. Invalid strata would compromise area-weighted carbon stock estimates.

### 4. Load and Merge Sample Data

**What it does:**
- Reads lab measurements (SOC, bulk density) for each depth interval
- Calculates depth midpoints and interval thickness
- Merges spatial coordinates with lab data
- Identifies orphaned records (samples without locations, locations without samples)

**Key calculations:**
```
depth_midpoint = (depth_top_cm + depth_bottom_cm) / 2
interval_thickness = depth_bottom_cm - depth_top_cm
```

---

## Quality Control & Validation

### Range Validation

**What it does:**
Flags samples that fall outside scientifically plausible ranges:

| Parameter | Valid Range | Flag Name |
|-----------|-------------|-----------|
| SOC | 0.1 - 600 g/kg | `qa_soc_valid` |
| Bulk Density | 0.05 - 2.5 g/cm³ | `qa_bd_valid` |
| Depth | 0 - 100 cm | `qa_depth_valid` |

**Scientific rationale:** 
- SOC >600 g/kg suggests pure organic matter or lab error
- BD <0.05 g/cm³ indicates measurement failure
- Depth inversions (top > bottom) are physically impossible

### Comprehensive QA Flagging

**What it does:**
Creates boolean flags for each validation check:
- `qa_spatial_valid`: Valid GPS coordinates
- `qa_depth_valid`: Depth intervals are logical and complete
- `qa_soc_valid`: SOC within plausible range
- `qa_bd_valid`: Bulk density within plausible range
- `qa_stratum_valid`: Stratum matches valid list
- `qa_pass`: All checks passed

**Only samples with `qa_pass = TRUE` proceed to analysis.**

---

## VM0033 Compliance Assessment ****** Change to Sampling design comparison ******

### Sample Size Requirements

**What it does:**
Evaluates whether the current sampling distribution to a theoretical optimum:

**Minimum Requirements:**
- Achieved precision ≤20% at 90% confidence interval or ≤30% at 95% confidence interval

**Calculations:**

<details>
<summary><b>Click to view statistical formulas</b></summary>

**Required sample size for target precision:**
```
n = (z × CV / target_precision)²

Where:
  z = 1.96 for 95% confidence
  CV = coefficient of variation (%)
  target_precision = 20% (VM0033 standard)
```

**Achieved precision with current sample size:**
```
precision = (z × CV) / √n
```

</details>

***** Also show lines of code where this should be *******


**Output metrics by stratum:**
- Current sample size
- Coefficient of variation (CV)
- Achieved precision at 95% CI
- Required n for 20%, 15%, 10% precision
- Additional cores needed
- Compliance status (EXCELLENT, GOOD, ACCEPTABLE, MARGINAL, POOR, INSUFFICIENT)

**Scientific rationale:** VM0033 requires demonstrating statistical confidence in carbon stock estimates. Insufficient sampling leads to high uncertainty and conservative (lower) creditable carbon estimates.

### Compliance Status Categories

| Status | Criteria | Interpretation |
|--------|----------|----------------|
| **EXCELLENT** | ≤10% precision | High confidence, minimal additional sampling needed |
| **GOOD** | ≤15% precision | Sufficient for most applications |
| **ACCEPTABLE** | ≤20% precision | Meets VM0033 minimum standard |
| **MARGINAL** | 20-30% precision | May require additional sampling |
| **POOR** | ≥30% precision | High uncertainty, more sampling strongly recommended |
| **INSUFFICIENT** | <3 cores | Does not meet VM0033 minimum |

---

## Bulk Density Handling

### Measurement vs. Default Values

**What it does:**
1. Identifies samples with measured bulk density
2. For samples without measurements, assigns ecosystem-specific defaults
3. Tracks which values are measured vs. estimated
4. Generates transparency reports

**Default bulk density values (from Morris et al. 2016):**

| Ecosystem | Default BD (g/cm³) | Source |
|-----------|-------------------|--------|
| Saltmarsh | 0.52 | Empirical average for emergent marsh |
| Seagrass | 0.89 | Global seagrass synthesis |
| Mangrove | 0.38 | Forest floor sediments |

**Transparency tracking:**
- `bd_measured`: Boolean flag (TRUE = field measured, FALSE = default applied)
- `bd_estimated`: Complementary flag for clarity

**Scientific implications:**
- Measured BD reduces uncertainty in carbon stock estimates
- Default BD introduces systematic uncertainty (typically 15-20% CV)
- VM0033 recommends measuring BD for all cores when feasible
- Mixed measurement approaches require uncertainty propagation in final estimates

**Output:** `bd_transparency_report.csv` shows percentage of measured vs. estimated BD by stratum.

---

## Carbon Stock Calculation

### Standard Formula

Carbon stock for each depth interval is calculated using:

```
Carbon Stock (kg C/m²) = SOC (g/kg) × BD (g/cm³) × depth_increment (cm) / 100

Where:
  SOC = Soil organic carbon content (g C / kg dry soil)
  BD = Bulk density (g dry soil / cm³)
  depth_increment = depth_bottom_cm - depth_top_cm
  100 = Conversion factor for units
```

**Unit conversion logic:**
```
g C/kg × g/cm³ × cm / 100 = kg C/m²

For VM0033 reporting (Mg C/ha):
kg C/m² × 10 = Mg C/ha
```

### Uncertainty Propagation

**What it does:**
Calculates standard error for carbon stock estimates using first-order Taylor approximation:

```
Relative Variance (C stock) = Relative Variance (SOC) + Relative Variance (BD)

Where:
  Relative Variance = (SE / mean)²
```

**Conservative defaults when SE not provided:**
- SOC: 10% coefficient of variation
- BD (measured): 15% CV
- BD (default): 20% CV (higher uncertainty)

**Scientific rationale:** Error propagation is required for VM0033 compliance. The 95% lower confidence bound is used for conservative carbon crediting.

### Core-Level Totals

**What it does:**
- Sums carbon stocks across all depth intervals for each core
- Tracks maximum depth sampled per core
- Calculates stratum-level statistics (mean, SD, min, max)

**Output:** `core_totals.csv` contains total carbon stock (0-100cm) for each core.

---

## Depth Profile Completeness

### Completeness Metrics

**What it does:**
Assesses how completely each core samples the 0-100cm target depth:

```
Completeness (%) = (Total depth sampled / 100 cm) × 100
```

**Classification:**
- **Complete:** ≥90% of target depth sampled
- **Good:** 70-89% sampled
- **Moderate:** 50-69% sampled
- **Incomplete:** <50% sampled

**Gap detection:**
Identifies depth intervals >5cm with no samples (potential core loss or sampling errors).

**Scientific rationale:** 
- VM0033 requires sampling to 100cm depth or refusal
- Incomplete profiles underestimate total carbon stocks
- Gaps in depth profiles may indicate:
  - Core loss during extraction
  - Sampling equipment limitations
  - Natural physical barriers (gravel layers, consolidated sediment)

**Output:** `core_depth_completeness.csv` shows completeness metrics and gap locations for each core.

---

## Core Type Analysis

### High-Resolution vs. Composite Cores

**What it does:**
Compares two common field sampling strategies:

**Core Types:**
1. **High-Resolution (HR):** Continuous 2cm intervals providing detailed depth profiles
2. **Paired Composite:** Composite samples at standard depths (0-15, 15-30, 30-50, 50-100cm), collected adjacent to HR core
3. **Unpaired Composite:** Composite samples at standalone locations

**Statistical validation:**
- Two-sample t-tests compare SOC between HR and Paired Composite cores within each stratum
- Tests null hypothesis: "Paired composite cores accurately represent HR measurements"

**Decision criteria:**
- **p ≥ 0.05 (not significant):** Paired composites can serve as proxy for HR cores → cost-effective sampling validated
- **p < 0.05 (significant):** HR and composite cores differ → analyze separately or collect more paired samples

**Scientific rationale:** 
High-resolution cores are expensive and time-intensive. If paired composite sampling provides statistically equivalent estimates, it offers a cost-effective monitoring approach. This analysis validates (or invalidates) that assumption for your specific site.

**Output:** 
- `core_type_summary.csv`: Descriptive statistics by type and stratum
- `core_type_statistical_tests.csv`: t-test results with recommendations

---

## Outputs

### Core Data Files

**Location:** `data_processed/`

| File | Description | Format |
|------|-------------|--------|
| `cores_clean_bluecarbon.rds` | Complete cleaned dataset (all QA flags) | R binary |
| `cores_clean_bluecarbon.csv` | Same as above | CSV (portable) |
| `core_totals.rds` | Total carbon stock per core (0-100cm) | R binary |
| `core_totals.csv` | Same as above | CSV |
| `cores_summary_by_stratum.csv` | Summary statistics by stratum | CSV |
| `carbon_by_stratum_summary.csv` | Carbon stock summaries | CSV |

### Diagnostic Reports

**Location:** `diagnostics/data_prep/`

| Report | Purpose |
|--------|---------|
| `vm0033_compliance_report.csv` | Sample size validation and power analysis |
| `core_depth_completeness.csv` | Depth coverage for each core |
| `depth_completeness_summary.csv` | Overall depth coverage statistics |
| `core_type_summary.csv` | HR vs Composite comparison (if applicable) |
| `core_type_statistical_tests.csv` | Statistical test results (if applicable) |

### Quality Assurance

**Location:** `diagnostics/qaqc/`

| Report | Purpose |
|--------|---------|
| `bd_transparency_report.csv` | Measured vs estimated bulk density |
| `duplicate_locations_[scenario]_[year].csv` | Cores at identical GPS locations |
| `qa_report.rds` | Comprehensive QA summary object |

### Log Files

**Location:** `logs/`

- `data_prep_[date].log`: Timestamped processing log with warnings and errors

---

## Key QA/QC Checks Performed

<details>
<summary><b>Spatial Quality Control</b></summary>

- ✓ Coordinates within valid lat/lon ranges
- ✓ No missing GPS data
- ✓ Duplicate location detection (±1m precision)
- ✓ Spatial clustering analysis

</details>

<details>
<summary><b>Measurement Validity</b></summary>

- ✓ SOC within 0.1-600 g/kg
- ✓ Bulk density within 0.05-2.5 g/cm³
- ✓ Depth intervals are logical (top < bottom)
- ✓ Depth values ≤100cm maximum

</details>

<details>
<summary><b>Data Completeness</b></summary>

- ✓ All cores have location data
- ✓ All samples have core assignments
- ✓ Required columns present
- ✓ Stratum assignments valid

</details>

<details>
<summary><b>VM0033 Compliance</b></summary>

- ✓ Minimum 3 cores per stratum
- ✓ Statistical power analysis
- ✓ Precision target assessment (20% at 95% CI)
- ✓ Sample size recommendations

</details>

<details>
<summary><b>Measurement Transparency</b></summary>

- ✓ Bulk density source tracking (measured vs default)
- ✓ Depth profile completeness
- ✓ Gap detection in core sequences
- ✓ Core type validation (HR vs Composite)

</details>

---

## Interpretation Guide

### Critical Warnings to Address

**🔴 High Priority:**
- **"INSUFFICIENT cores"** in VM0033 compliance → Collect more samples before proceeding
- **"Invalid coordinates"** → Fix GPS data in source files
- **"Poor precision (≥30%)"** → High uncertainty, additional sampling strongly recommended

**🟡 Medium Priority:**
- **"Duplicate locations detected"** → Verify if intentional paired sampling or GPS error
- **"Non-monotonic profiles"** → Review for measurement errors or natural stratification
- **"Marginal precision (20-30%)"** → Consider additional sampling for tighter confidence intervals

**🟢 Low Priority:**
- **"BD estimated for X samples"** → Acceptable, but measured BD preferred for lower uncertainty
- **"Depth gaps detected"** → Document in field notes, may affect interpolation

### When to Collect More Samples

**Absolute requirements:**
- Any stratum with <3 cores (VM0033 minimum)
- Any stratum with achieved precision >20% (fails VM0033 standard)

**Recommended additional sampling:**
- Strata with precision 20-30% (marginal compliance)
- High variability strata (CV >50%)
- Strata with <90% depth completeness
- Areas with significant spatial clustering

---

## Next Steps

After successful completion of Module 01:

1. **Review QA Reports:**
   - Check `vm0033_compliance_report.csv` for sampling adequacy
   - Review `bd_transparency_report.csv` for measurement coverage
   - Examine `qa_report.rds` for overall data quality

2. **Address Critical Issues:**
   - Collect additional samples for non-compliant strata
   - Verify and correct invalid data
   - Re-run Module 01 if source data was modified

3. **Proceed to Module 02:**
   ```r
   source('02_exploratory_analysis_bluecarbon.R')
   ```
   This will generate visualizations and deeper statistical summaries.

---

## Scientific References

**Bulk Density Defaults:**
- Morris, J.T., Barber, D.C., Callaway, J.C., et al. (2016). "Contributions of organic and inorganic matter to sediment volume and accretion in tidal wetlands at steady state." *Earth's Future*, 4(4), 110-121.

**VM0033 Standards:**
- Verra (2015). VM0033 Methodology for Tidal Wetland and Seagrass Restoration, v2.0.

**Carbon Stock Calculation:**
- Howard, J., Hoyt, S., Isensee, K., Pidgeon, E., & Telszewski, M. (2014). *Coastal Blue Carbon: Methods for Assessing Carbon Stocks and Emissions Factors*. Conservation International, UNESCO, IUCN.

**Statistical Power Analysis:**
- IPCC (2006). *Guidelines for National Greenhouse Gas Inventories*, Volume 4: Agriculture, Forestry and Other Land Use.

---

## Troubleshooting

**Problem:** "Missing required columns in locations/samples"  
**Solution:** Check CSV file headers match exactly: `core_id`, `longitude`, `latitude`, `stratum`, `depth_top_cm`, `depth_bottom_cm`, `soc_g_kg`

**Problem:** "Invalid stratum names detected"  
**Solution:** Review `blue_carbon_config.R` → `VALID_STRATA` and ensure source data uses exact names (case-sensitive)

**Problem:** "No samples passed QA"  
**Solution:** Check QC thresholds in config file. May need to adjust for your ecosystem (e.g., high-organic peat systems may have SOC >300 g/kg)

**Problem:** Module runs but produces no outputs  
**Solution:** Check file paths - ensure `data_raw/` folder exists with required CSVs. Review log file in `logs/` for specific errors.

---

## Summary

Module 01 transforms raw field measurements into a scientifically validated dataset by:

1. ✅ Loading and merging core location and sample data
2. ✅ Performing comprehensive quality control (spatial, measurement, completeness)
3. ✅ Validating VM0033 compliance (sample size, precision, power)
4. ✅ Handling bulk density (measured vs defaults with transparency)
5. ✅ Calculating carbon stocks with uncertainty propagation
6. ✅ Assessing depth profile completeness and gaps
7. ✅ Comparing core type strategies (HR vs Composite)
8. ✅ Generating comprehensive QA/QC reports

**Output:** Clean, validated dataset (`cores_clean_bluecarbon.rds/csv`) ready for depth harmonization and spatial modeling.
