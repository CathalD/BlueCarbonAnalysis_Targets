# Module 02: Blue Carbon Exploratory Data Analysis

## Overview

This module performs comprehensive exploratory data analysis (EDA) on cleaned blue carbon datasets. It generates diagnostic visualizations, statistical summaries, and quality assessments to reveal patterns in carbon distribution, validate data quality, and inform subsequent spatial modeling decisions.

**Primary Goal:** Understand spatial and vertical patterns in carbon stocks through visualization and statistical analysis, identifying potential data quality issues before proceeding to spatial prediction.

---

## Table of Contents

1. [Purpose and Scope](#purpose-and-scope)
2. [Input Requirements](#input-requirements)
3. [Statistical Quality Checks](#statistical-quality-checks)
4. [Visualization Outputs](#visualization-outputs)
5. [Interpretation Guide](#interpretation-guide)
6. [Decision Points](#decision-points)
7. [Outputs](#outputs)

---

## Purpose and Scope

### Why Exploratory Analysis Matters

Before building spatial models, it's critical to:

1. **Understand data structure** - Identify patterns, trends, and distributions
2. **Detect anomalies** - Find outliers and measurement errors that survived QA
3. **Assess variability** - Quantify spatial and vertical heterogeneity
4. **Validate assumptions** - Check if data meets modeling requirements
5. **Guide decisions** - Inform stratification, model selection, and sampling design

### What This Module Does

**Statistical Analysis:**
- Outlier detection using Isolation Forest algorithm
- Depth profile monotonicity checks
- Summary statistics by stratum
- Correlation analysis between variables

**Visual Diagnostics:**
- Spatial distribution of cores with compliance status
- SOC and bulk density distributions
- Depth profiles with confidence intervals
- Carbon stock comparisons
- Core type validation (HR vs Composite)
- VM0033 compliance dashboards
- Scenario comparisons (PROJECT vs BASELINE)

---

## Input Requirements

### Required Files

**From Module 01:**
- `data_processed/cores_clean_bluecarbon.rds` - Cleaned core dataset (primary input)

**Optional (Enhanced Visualizations):**
- `data_processed/vm0033_compliance.rds` - Compliance metrics
- `data_processed/bd_transparency.rds` - Bulk density source tracking
- `data_processed/depth_completeness.rds` - Profile completeness data
- `data_processed/hr_composite_comparison.rds` - Core type statistics

**Configuration:**
- `blue_carbon_config.R` - Stratum colors, figure dimensions, project metadata

---

## Statistical Quality Checks

### 1. Outlier Detection (Isolation Forest)

**What it does:**
Uses ensemble machine learning to identify statistical outliers in:
- Soil organic carbon (SOC)
- Bulk density (BD)
- Depth measurements

**Method:** Isolation Forest (Liu et al. 2008)
- Identifies samples that are "easy to isolate" from the main distribution
- Works with multivariate data (considers SOC, BD, and depth simultaneously)
- No assumptions about data distribution (robust to skewed distributions)

**Contamination rate:** Default 5% (expects ~5% of samples to be outliers)

**Output:**
- `outlier_flag`: Boolean flag for each sample
- `outlier_score`: Anomaly score (higher = more anomalous)
- Saves report: `diagnostics/qaqc/outliers_[scenario]_[year].csv`

**Scientific rationale:** 
Traditional outlier methods (Z-score, IQR) assume normal distributions and univariate analysis. Blue carbon data often shows:
- Skewed distributions (high SOC in organic-rich sediments)
- Multivariate relationships (SOC-BD inverse correlation)
- Stratum-specific patterns

Isolation Forest handles these characteristics better than classical methods.

**When to investigate outliers:**
- Outlier score >95th percentile
- SOC values that seem implausible for ecosystem type
- Samples flagged across multiple variables simultaneously

**Common causes of outliers:**
- Lab measurement errors (decimal point errors, unit conversions)
- Field sampling errors (contamination, mixed samples)
- Natural variability (buried organic layers, root masses)
- Data entry errors

### 2. Depth Profile Monotonicity

**What it does:**
Checks if SOC decreases with depth (typical pattern in most ecosystems).

**Typical pattern:**
- Surface layers (0-15cm): Highest SOC from active organic matter input
- Deeper layers (50-100cm): Lower SOC from older, more decomposed material

**Non-monotonic profiles flagged:**
Samples where SOC *increases* significantly with depth.

**Scientific interpretation:**

| Pattern | Likely Cause |
|---------|-------------|
| SOC increases 15-30cm | Possible buried organic layer or peat deposit |
| SOC increases 50-100cm | Historical deposition event or measurement error |
| Multiple increases | Stratified deposits or systematic measurement issue |

**Output:** `diagnostics/qaqc/non_monotonic_profiles_[scenario]_[year].csv`

**When to investigate:**
- Single cores with multiple reversals → Check lab measurements
- All cores in a stratum show reversals → May be real ecological pattern
- Large magnitude reversals (>50 g/kg) → Verify field notes and photos

**Important:** Non-monotonic profiles aren't always errors! Some ecosystems naturally show:
- Buried peat or organic horizons
- Historical high-productivity periods
- Mangrove root zones at depth
- Seagrass rhizome layers

---

## Visualization Outputs

### Plot 01: Spatial Distribution of Cores

**File:** `01_spatial_distribution.png`

**What it shows:**
- GPS locations of all cores
- **Shape** = Core type (HR, Paired Composite, Unpaired Composite, Unknown)
- **Fill color** = Stratum (from configuration)
- **Border color** = VM0033 compliance status
  - Green border = Excellent/Good compliance
  - Yellow border = Acceptable compliance
  - Red border = Poor/Insufficient compliance

**Scientific insights:**
- **Spatial clustering:** Are cores evenly distributed or clustered?
- **Stratification coverage:** Does each stratum have adequate spatial representation?
- **Compliance patterns:** Which areas need additional sampling?

**Decision points:**
- Clustered cores may violate spatial independence assumptions for kriging
- Gaps in spatial coverage reduce interpolation accuracy
- Non-compliant strata (red borders) require additional field sampling

---

### Plot 01B: VM0033 Compliance Dashboard

**File:** `01b_vm0033_compliance_dashboard.png`

**Two-panel visualization:**

**Panel 1: Sample Size vs Required**
- Gray bars = Required cores for 20% precision target
- Colored bars = Current cores (color indicates status)
- Red dashed line = VM0033 minimum (3 cores)

**Panel 2: Achieved Precision**
- Horizontal bars showing precision (%) for each stratum
- Green line = 20% target (acceptable)
- Blue line = 10% target (excellent)

**Scientific interpretation:**

| Achieved Precision | Interpretation | Action |
|-------------------|----------------|--------|
| ≤10% | Excellent - high confidence estimates | No action needed |
| 10-15% | Good - suitable for most applications | Consider a few more cores for excellence |
| 15-20% | Acceptable - meets VM0033 minimum | Adequate for crediting, but consider improvement |
| 20-30% | Marginal - high uncertainty | Additional sampling recommended |
| >30% | Poor - very high uncertainty | More sampling required before modeling |

**Power analysis interpretation:**
If "Current cores" bar is shorter than gray "Required" bar → collect more samples to achieve 20% precision target.

---

### Plot 02: SOC Distribution by Stratum

**File:** `02_soc_distribution_by_stratum.png`

**Two-panel visualization:**

**Panel 1: Box plots with 95% confidence intervals**
- Box shows quartiles (25th, median, 75th percentile)
- Whiskers show data range (excluding outliers)
- Black diamond = Mean SOC
- Black error bars = 95% confidence interval around mean

**Panel 2: Violin plots**
- Width shows probability density (wider = more samples at that value)
- Embedded box plot for reference

**Scientific insights:**

**Distribution shape:**
- **Symmetric distribution:** Normal variation expected
- **Right-skewed (long tail high):** A few very high SOC samples (organic-rich pockets)
- **Left-skewed (long tail low):** A few very low SOC samples (mineral sediment)
- **Bimodal (two peaks):** May indicate mixed ecosystem types or measurement issues

**Variability patterns:**
- Wide boxes = high within-stratum variability
- Narrow boxes = homogeneous conditions
- Large CI = high uncertainty (may need more samples)
- Small CI = precise estimates

**Stratum comparisons:**
- Which strata have highest/lowest carbon?
- Do strata show similar or different distributions?
- Are there unexpected patterns?

---

### Plot 03: Depth Profiles by Stratum

**File:** `03_depth_profiles_by_stratum.png`

**What it shows:**
- Mean SOC at 2cm depth intervals (0-2, 2-4, 4-6, etc.)
- Line shows trend through depth profile
- Point size = Number of samples at that depth
- Faceted by stratum for comparison

**Typical patterns:**

**Expected decreasing profile:**
```
Depth (cm)    SOC (g/kg)
0-15          High (100-300)
15-30         Moderate (50-150)
30-50         Low-Moderate (30-100)
50-100        Low (20-80)
```

**Scientific interpretation:**

| Pattern | Ecological Meaning |
|---------|-------------------|
| Steep decline 0-30cm | Active surface layer, rapid decomposition |
| Gradual decline | Slow decomposition, good preservation |
| Flat profile | Homogeneous deposits or well-mixed sediment |
| Increase at depth | Buried organic layer, historical deposition |

**Data quality indicators:**
- **Smooth curves** = good data quality
- **Erratic curves** = measurement noise or natural variability
- **Small points (low n)** = sparse sampling at that depth
- **Large points (high n)** = well-sampled depth interval

**Stratum comparison:**
- Do all strata show similar depth trends?
- Which strata have highest surface carbon?
- Which strata have best carbon preservation at depth?

---

### Plot 04: Bulk Density Patterns

**Files:** 
- `04_bulk_density_patterns.png` (BD by stratum, BD vs SOC)
- `04b_bulk_density_vs_depth.png` (BD depth profiles by stratum)

**What it shows:**

**Panel 1: BD by Stratum**
- Distribution of bulk density values within each stratum
- Identifies which ecosystems have highest/lowest density

**Panel 2: BD vs SOC Relationship**
- Scatter plot with linear regression
- Shows inverse relationship (high SOC → low BD)

**Scientific insights:**

**Typical BD ranges:**
- **Mangrove (0.2-0.6 g/cm³):** Low density, high organic content
- **Saltmarsh (0.4-0.8 g/cm³):** Moderate density
- **Seagrass (0.6-1.2 g/cm³):** Higher density, more mineral sediment

**BD-SOC relationship:**
- **Strong negative correlation** = expected pattern (organic matter is less dense than mineral matter)
- **Weak/no correlation** = may indicate measurement issues or unusual sediment composition
- **Positive correlation** = data quality issue (check measurements)

**BD depth trends:**
- **Increasing with depth** = compaction of older sediments
- **Decreasing with depth** = buried organic layers
- **Constant** = well-mixed or recent deposition

**Implications for carbon stocks:**
- Low BD + high SOC = very high carbon density (per unit volume)
- High BD + low SOC = low carbon density
- BD variability affects total stock uncertainty

---

### Plot 04C: Bulk Density Transparency

**File:** `04c_bd_transparency.png`

**What it shows:**
Stacked bar chart showing percentage of cores with:
- **Green bars:** Field-measured bulk density
- **Yellow bars:** Default bulk density (from literature)

**Target:** 100% measured BD for lowest uncertainty

**Scientific rationale:**

**Measured BD advantages:**
- Site-specific values
- Accounts for local sediment composition
- Reduces uncertainty in carbon stock estimates
- Typically 10-15% coefficient of variation

**Default BD limitations:**
- Generic ecosystem averages
- May not reflect local conditions
- Introduces systematic uncertainty
- Typically 15-20% coefficient of variation

**Interpretation guide:**

| % Measured | Assessment | Recommendation |
|-----------|------------|----------------|
| 80-100% | Excellent | Proceed with confidence |
| 60-79% | Good | Consider measuring BD for remaining cores |
| 40-59% | Moderate | Uncertainty will be higher; measure BD if possible |
| <40% | Poor | High reliance on defaults; prioritize BD measurement |

**Impact on uncertainty:**
- Each core with default BD adds ~5% to overall uncertainty
- Strata with <50% measured BD may not meet precision targets
- VM0033 verification may require conservative discounting

---

### Plot 04D: Depth Profile Completeness

**File:** `04d_depth_completeness.png`

**What it shows:**
Stacked bar chart by stratum showing cores categorized as:
- **Dark green:** Complete (≥90% of 0-100cm sampled)
- **Light green:** Good (70-89%)
- **Yellow:** Moderate (50-69%)
- **Orange:** Incomplete (<50%)

**Target:** All cores "Complete" or "Good"

**Scientific implications:**

**Incomplete profiles affect:**
1. **Total stock estimates:** Underestimate carbon if deep layers missed
2. **Depth harmonization:** Difficult to extrapolate beyond sampled depth
3. **Temporal comparisons:** Inconsistent depth coverage between years
4. **Model training:** Incomplete profiles reduce data quality for RF models

**Common causes:**
- **Refusal at depth:** Dense clay, gravel, or consolidated sediment
- **Core loss:** Material fell out during extraction
- **Equipment limitations:** Core barrel too short
- **Sampling protocol:** Intentional shallow sampling (0-50cm)

**Gap analysis:**
Text annotation shows cores with depth gaps >5cm:
- Small gaps (5-10cm) = acceptable, common at section breaks
- Large gaps (>20cm) = problematic, may indicate core loss or protocol issues

**Decision criteria:**
- If >20% of cores are "Incomplete" → consider re-sampling
- If gaps are systematic (e.g., all cores missing 40-60cm) → check field protocol
- If isolated gaps → document in metadata, may exclude from depth harmonization

---

### Plot 05: Carbon Stock by Stratum

**File:** `05_carbon_stock_by_stratum.png`

**Two-panel visualization:**

**Panel 1: Carbon Stock per Sample**
- Shows carbon density (kg C/m²) for each depth interval
- Black diamonds = Mean with 95% CI
- Compares sample-level variability across strata

**Panel 2: Total Carbon Stock per Core**
- Summed carbon stock (0-100cm) for each core
- Shows total ecosystem carbon storage capacity

**Scientific interpretation:**

**Sample-level patterns:**
- High variability = heterogeneous sediment composition
- Low variability = homogeneous conditions
- Outliers = potential measurement errors or natural hotspots

**Core-level totals:**
Typical ranges for 0-100cm:

| Ecosystem | Carbon Stock (kg C/m²) | Equivalent (Mg C/ha) |
|-----------|------------------------|----------------------|
| Mangrove | 20-80 | 200-800 |
| Saltmarsh | 15-60 | 150-600 |
| Seagrass | 10-50 | 100-500 |

**Stratum comparison:**
- Which strata are carbon "hotspots"?
- Is there overlap between strata (similar stocks)?
- Are differences statistically significant (check CI overlap)?

**Decision points:**
- Strata with very different stocks should remain separate in spatial models
- Strata with overlapping CIs might be combined to increase sample size
- High variability strata may need stratification refinement

---

### Plot 06: Core Type Comparison (HR vs Composite)

**File:** `06_core_type_comparison.png`

**Two-panel visualization:**

**Panel 1: Core Count by Type and Stratum**
- Stacked bars showing distribution of HR, Paired Composite, Unpaired Composite cores

**Panel 2: SOC Distribution by Core Type**
- Box plots comparing SOC values between core types
- Faceted by stratum
- P-values shown if statistical tests available

**Scientific rationale:**

**Research question:** 
Can composite cores (cheaper, faster) serve as statistically valid proxies for high-resolution cores (expensive, time-intensive)?

**Statistical validation:**
- **p ≥ 0.05 (not significant):** Core types are statistically equivalent → paired sampling approach validated
- **p < 0.05 (significant):** Core types differ significantly → analyze separately or collect more paired samples

**Interpretation by outcome:**

**Case 1: No significant difference (p ≥ 0.05)**
- **Meaning:** Paired composites accurately represent HR measurements
- **Implication:** Cost-effective monitoring strategy validated
- **Recommendation:** Continue paired composite approach for future sampling
- **Benefit:** ~60% cost reduction per core while maintaining accuracy

**Case 2: Significant difference (p < 0.05)**
- **Meaning:** Composite cores systematically differ from HR cores
- **Possible causes:**
  - Spatial heterogeneity (cores not truly paired)
  - Composite sampling bias (preferential sample selection)
  - Different depth integration (HR shows finer-scale variation)
- **Recommendation:** 
  - Analyze HR and Composite separately
  - Increase paired sample size to reassess
  - Use only HR cores for spatial modeling

**Economic implications:**
- HR cores: $150-300 per core (lab costs for ~50 depth intervals)
- Composite cores: $60-120 per core (lab costs for 4-5 intervals)
- If validated, paired approach saves ~$100 per core while maintaining statistical validity

---

### Plot 07: QA/QC Summary

**File:** `07_qa_summary.png`

**What it shows:**
Horizontal bar chart showing pass/fail rates for each QA check:
- Spatial validity (GPS coordinates)
- Depth validity (logical intervals)
- SOC validity (plausible range)
- BD validity (plausible range)
- Stratum validity (correct assignments)
- Overall pass (all checks passed)

**Target:** >95% pass rate for all checks

**Scientific interpretation:**

| Pass Rate | Assessment | Action |
|-----------|------------|--------|
| >95% | Excellent data quality | Proceed confidently |
| 90-95% | Good quality | Review failed samples |
| 80-90% | Moderate quality | Investigate systematic issues |
| <80% | Poor quality | Review field/lab protocols |

**Common failure patterns:**

**Systematic failures (many samples fail):**
- **All SOC invalid:** Check lab calibration or unit conversions
- **All BD invalid:** Verify measurement protocol
- **Many depth invalid:** Review field data entry

**Isolated failures (few samples fail):**
- **Random failures:** Likely data entry errors or one-off measurement issues
- **Single core failures:** Core-specific problem (contamination, poor preservation)

---

### Plot 07B: Scenario Comparison (PROJECT vs BASELINE)

**File:** `07b_scenario_comparison.png`

**Only created if:** Multiple scenarios present (e.g., PROJECT, BASELINE, CONTROL)

**Two-panel visualization:**

**Panel 1: SOC by Scenario and Stratum**
- Bars show mean SOC for each scenario
- Error bars = 95% confidence intervals

**Panel 2: Carbon Stock by Scenario and Stratum**
- Bars show mean carbon stock
- Error bars = 95% confidence intervals

**Scientific rationale:**

**Temporal monitoring applications:**
- **BASELINE:** Pre-restoration/pre-project condition
- **PROJECT:** Post-restoration/current condition
- **CONTROL:** Reference site (no intervention)

**Key comparisons:**

**PROJECT vs BASELINE:**
- Are PROJECT SOC/stocks higher than BASELINE? (restoration success)
- Is difference statistically significant? (CI bars don't overlap)
- Which strata show greatest improvement?

**PROJECT vs CONTROL:**
- Are PROJECT and CONTROL similar? (natural recovery)
- If PROJECT > CONTROL → restoration effect (additionality)
- If PROJECT = CONTROL → natural recovery trajectory

**Statistical significance:**
- **Overlapping error bars:** Not significantly different
- **Non-overlapping error bars:** Significantly different (conservative test)
- **For rigorous testing:** Use proper statistical tests (t-test, ANOVA) in Module 08

**Carbon credit implications (VM0033):**
- Credits = (PROJECT - BASELINE) per year
- Must demonstrate PROJECT > BASELINE with statistical confidence
- Control site helps separate project effects from natural variability

---

### Plot 07C: Correlation Matrix

**File:** `07c_correlation_matrix.png`

**What it shows:**
Heatmap showing Pearson correlations between:
- SOC (g/kg)
- Bulk density (g/cm³)
- Depth (cm)
- Carbon stock (kg/m²)

**Color scale:**
- **Green:** Positive correlation
- **White:** No correlation (r ≈ 0)
- **Red:** Negative correlation

**Expected patterns:**

| Variable Pair | Expected Correlation | Interpretation if Unexpected |
|--------------|---------------------|------------------------------|
| SOC ↔ BD | **Strong negative** (r = -0.6 to -0.8) | Organic matter reduces density |
| SOC ↔ Depth | **Moderate negative** (r = -0.3 to -0.5) | SOC decreases with depth |
| SOC ↔ C Stock | **Strong positive** (r = 0.7 to 0.9) | SOC drives carbon storage |
| BD ↔ Depth | **Weak positive** (r = 0.1 to 0.3) | Compaction with depth |

**Interpretation:**

**Strong correlations (|r| > 0.7):**
- Expected: SOC ↔ Carbon Stock (by definition)
- Unexpected: May indicate multicollinearity issues for modeling

**Weak correlations (|r| < 0.3):**
- Expected for some pairs
- If SOC ↔ BD is weak → check data quality (should be strong negative)

**Unexpected patterns:**
- **SOC ↔ BD positive:** Data quality issue (impossible physically)
- **SOC ↔ Depth positive:** Check for inverted depth or non-monotonic profiles
- **BD ↔ C Stock positive without SOC involvement:** Calculation error

---

### Plot 08: Summary Statistics Table

**File:** `08_summary_table.png`

**What it shows:**
Formatted table with key statistics by stratum:
- Number of cores
- Number of samples
- Mean core depth (cm)
- SOC (mean ± SD in g/kg)
- Bulk density (mean ± SD in g/cm³)
- Total carbon stock (kg/m²)

**Use cases:**
- Quick reference for manuscript tables
- Summary for reports and presentations
- Comparison across strata
- Input validation for spatial models

---

## Interpretation Guide

### Critical Patterns to Investigate

**🔴 High Priority (Stop and Investigate):**

1. **>10% outliers detected**
   - May indicate systematic measurement issues
   - Review outlier report and validate measurements
   - Consider excluding clear errors before modeling

2. **SOC-BD positive correlation**
   - Physically impossible (organic matter is less dense)
   - Check for data entry errors, unit conversions, or swapped columns

3. **>30% non-monotonic profiles**
   - If widespread, may be real ecological pattern (buried layers)
   - If isolated, likely measurement errors

4. **QA pass rate <90%**
   - Significant data quality issues
   - Review failed samples and source data
   - May need to return to Module 01 after corrections

5. **Stratum overlap in carbon stocks (CIs completely overlapping)**
   - May indicate poor stratification
   - Consider combining strata or refining boundaries

**🟡 Medium Priority (Document and Monitor):**

1. **5-10% outliers**
   - Typical for natural systems
   - Document outliers and justification for inclusion/exclusion

2. **High CV in key strata (>50%)**
   - Natural heterogeneity or measurement noise?
   - May need additional samples or stratification refinement

3. **Uneven spatial distribution**
   - Clustered samples may affect kriging
   - Gaps in coverage reduce prediction accuracy

4. **BD transparency <60% measured**
   - Higher uncertainty in carbon stocks
   - Prioritize BD measurement in future sampling

**🟢 Low Priority (Acceptable, Note for Discussion):**

1. **<5% outliers**
   - Expected for large datasets
   - Document reasoning for inclusion

2. **Some non-monotonic profiles**
   - May reflect real ecology (buried layers)
   - Document in site description

3. **BD transparency 60-80% measured**
   - Acceptable balance
   - Note uncertainty is slightly elevated

---

## Decision Points

### Should I Collect More Samples?

**Yes, if:**
- ✗ Any stratum shows "INSUFFICIENT" or "POOR" compliance
- ✗ VM0033 precision >20% in strata intended for crediting
- ✗ High spatial clustering (>30% of cores within 10m of another)
- ✗ Major gaps in spatial coverage
- ✗ QA pass rate <90%

**Consider, if:**
- ⚠ Precision 20-30% (marginal compliance)
- ⚠ High CV (>50%) suggests high natural variability
- ⚠ Core type comparison shows significant difference (need more paired samples)
- ⚠ Depth completeness <70% in many cores

**No additional sampling needed, if:**
- ✓ All strata show "ACCEPTABLE" or better compliance
- ✓ Precision ≤20% for all strata
- ✓ Even spatial distribution
- ✓ QA pass rate >95%
- ✓ Core types validated as equivalent

### Should I Refine My Stratification?

**Yes, if:**
- Strata show very similar carbon stocks (overlapping CIs)
- High within-stratum variability (CV >60%)
- Spatial patterns don't match stratum boundaries
- Depth profiles differ dramatically within a stratum

**No, if:**
- Strata show distinct carbon stocks (non-overlapping CIs)
- Within-stratum variability is moderate (CV 30-50%)
- Spatial patterns align with stratum boundaries

### Can I Use Composite Cores as Proxies for HR Cores?

**Yes (validated approach), if:**
- Statistical tests show p ≥ 0.05 (no significant difference)
- Visual distributions overlap substantially
- Sample sizes adequate (n ≥ 3 per type per stratum)

**No (analyze separately), if:**
- Statistical tests show p < 0.05 (significant difference)
- Systematic bias visible in plots
- Insufficient paired samples to validate

### Should I Proceed to Spatial Modeling?

**Yes, proceed to Module 03 (Depth Harmonization), if:**
- ✓ Data quality is acceptable (>90% pass rate)
- ✓ Outliers investigated and addressed
- ✓ VM0033 compliance adequate for intended use
- ✓ Patterns make scientific sense

**No, return to field/lab, if:**
- ✗ Sample sizes insufficient for any stratum
- ✗ Major data quality issues (QA <80%)
- ✗ Unexplainable patterns (e.g., SOC increases with depth everywhere)
- ✗ Precision targets not met and more samples feasible

---

## Outputs

### Visualization Files

**Location:** `outputs/plots/exploratory/`

| Plot | Description | Key Use |
|------|-------------|---------|
| `01_spatial_distribution.png` | Core locations with compliance status | Spatial sampling assessment |
| `01b_vm0033_compliance_dashboard.png` | Sample size and precision analysis | Compliance verification |
| `02_soc_distribution_by_stratum.png` | SOC box and violin plots | Variability assessment |
| `03_depth_profiles_by_stratum.png` | Mean SOC by depth | Vertical pattern analysis |
| `04_bulk_density_patterns.png` | BD distributions and relationships | BD quality check |
| `04b_bulk_density_vs_depth.png` | BD depth profiles | Compaction patterns |
| `04c_bd_transparency.png` | Measured vs default BD | Uncertainty assessment |
| `04d_depth_completeness.png` | Profile completeness by core | Sampling adequacy |
| `05_carbon_stock_by_stratum.png` | Carbon stock distributions | Stock estimation |
| `06_core_type_comparison.png` | HR vs Composite validation | Sampling strategy validation |
| `07_qa_summary.png` | QA pass/fail rates | Data quality overview |
| `07b_scenario_comparison.png` | PROJECT vs BASELINE | Temporal analysis |
| `07c_correlation_matrix.png` | Variable correlations | Relationship analysis |
| `08_summary_table.png` | Summary statistics table | Quick reference |

### Data Files

**Location:** `data_processed/`

| File | Description |
|------|-------------|
| `eda_summary.rds` | Comprehensive EDA results object |
| `eda_qc_summary.rds` | Quality check results |

**Location:** `diagnostics/qaqc/`

| File | Description |
|------|-------------|
| `outliers_[scenario]_[year].csv` | Detected outliers with scores |
| `non_monotonic_profiles_[scenario]_[year].csv` | Depth profile inversions |

### Log Files

**Location:** `logs/`

- `exploratory_analysis_[date].log` - Processing log with warnings

---

## Common Patterns and Their Meanings

### SOC Distribution Shapes

**Right-skewed (common):**
- A few very high SOC samples
- Typical in heterogeneous systems
- Consider log-transformation for modeling

**Bimodal:**
- Two distinct SOC populations
- May indicate:
  - Mixed ecosystem types within stratum
  - Measurement protocol inconsistency
  - Real ecological gradients (e.g., high/low marsh)

**Uniform:**
- Even distribution across range
- Unusual for natural systems
- May indicate data quality issues or highly homogeneous conditions

### Depth Profile Shapes

**Exponential decay (most common):**
```
Depth    SOC
0        ████████████ 200
20       ████████ 120
40       █████ 70
60       ███ 40
80       ██ 25
100      █ 15
```
- Typical pattern
- Surface carbon from recent deposition
- Decomposition/preservation gradient

**Linear decline:**
```
Depth    SOC
0        ████████████ 180
20       ██████████ 150
40       ████████ 120
60       ██████ 90
80       ████ 60
100      ██ 30
```
- Steady decomposition rate
- Common in stable, low-energy environments

**Subsurface maximum:**
```
Depth    SOC
0        ████████ 120
20       ██████████ 150  ← Peak
40       ████████ 120
60       █████ 80
80       ███ 50
100      ██ 30
```
- Buried organic layer
- Historical high-productivity event
- Investigate field notes and photos

---

## Next Steps

### After Successful EDA Completion:

1. **Review All Plots:**
   - Check `outputs/plots/exploratory/` directory
   - Identify any unexpected patterns
   - Document outliers and anomalies

2. **Compile Key Findings:**
   - Summarize stratum characteristics
   - Note compliance status
   - Document data quality issues

3. **Make Sampling Decisions:**
   - Determine if additional samples needed
   - Plan stratification refinements if needed
   - Validate core type approach

4. **Proceed to Module 03:**
   ```r
   source('03_depth_harmonization_bluecarbon.R')
   ```
   This will standardize depth profiles for spatial modeling.

### If Issues Detected:

**Data Quality Problems:**
- Return to source data
- Verify measurements
- Re-run Module 01 after corrections

**Insufficient Samples:**
- Plan additional field campaign
- Target non-compliant strata
- Fill spatial gaps

**Stratification Issues:**
- Refine stratum boundaries in GEE
- Re-export stratum masks
- Update configuration

---

## Scientific References

**Outlier Detection:**
- Liu, F.T., Ting, K.M., Zhou, Z.H. (2008). "Isolation Forest." *ICDM '08: Proceedings of the 2008 Eighth IEEE International Conference on Data Mining*, 413-422.

**Exploratory Data Analysis:**
- Tukey, J.W. (1977). *Exploratory Data Analysis*. Addison-Wesley.

**Carbon Profile Patterns:**
- Breithaupt, J.L., et al. (2012). "Organic carbon burial rates in mangrove sediments: Strengthening the global budget." *Global Biogeochemical Cycles*, 26(3).

**Statistical Methods:**
- Gotelli, N.J., & Ellison, A.M. (2004). *A Primer of Ecological Statistics*. Sinauer Associates.

---

## Summary

Module 02 provides comprehensive exploratory analysis by:

1. ✅ Detecting statistical outliers using machine learning
2. ✅ Checking depth profile quality and monotonicity
3. ✅ Visualizing spatial distribution and compliance status
4. ✅ Analyzing SOC and BD distributions across strata
5. ✅ Examining depth profiles and vertical patterns
6. ✅ Assessing carbon stock variability
7. ✅ Validating core type sampling strategies
8. ✅ Comparing scenarios (PROJECT vs BASELINE)
9. ✅ Evaluating measurement transparency and completeness
10. ✅ Generating comprehensive quality reports

**Output:** 14 diagnostic plots and statistical summaries that inform decisions about sampling adequacy, data quality, and readiness for spatial modeling.

**Key Decision:** Based on EDA results, either proceed to depth harmonization (Module 03) or return to field/lab to address identified issues.
