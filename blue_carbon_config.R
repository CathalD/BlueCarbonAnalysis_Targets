# ============================================================================
# BLUE CARBON PROJECT CONFIGURATION
# Edit the values below for your specific site and dataset.
# ============================================================================

# ── Project metadata ──────────────────────────────────────────────────────────
PROJECT_NAME     <- "BC_Coastal_BlueCarbon_2026_Example"
PROJECT_SCENARIO <- "PROJECT"   # BASELINE | PROJECT | CONTROL | DEGRADED
MONITORING_YEAR  <- 2026
PROJECT_LOCATION <- "Chemainus Test 1"

# ── Input file paths ──────────────────────────────────────────────────────────
PRE_ANALYSIS_DIR <- "Pre-Analysis Data Preparation"
DATA_RAW_DIR     <- file.path(PRE_ANALYSIS_DIR, "data_raw")
DATA_GLOBAL_DIR  <- file.path(PRE_ANALYSIS_DIR, "data_global")
COVARIATES_DIR   <- file.path(PRE_ANALYSIS_DIR, "covariates")

# Remote sensing covariate raster (multi-band GEE export)
COVARIATE_RASTER <- file.path(COVARIATES_DIR, "BlueCarbon_Covariate_Snapshot_25m_2020_2023.tif")

# Area of Interest boundary (GeoJSON, shapefile, or GPKG)
# Set to NULL if not yet available — Step 2 will return per-stratum densities only.
AOI_FILE <- file.path(DATA_RAW_DIR, "aoi_boundary.geojson")

# Column in AOI_FILE that identifies ecosystem strata.
# Set to NULL for a whole-site total (no per-stratum area breakdown).
AOI_STRATUM_FIELD <- NULL  # e.g. "ecosystem_type"

# ── Ecosystem strata ──────────────────────────────────────────────────────────
VALID_STRATA <- c("IM", "NM", "MF")

STRATUM_COLORS <- c(
  "IM" = "#FFFF99",   # Intertidal marsh
  "NM" = "#99FF99",   # Near-margin
  "MF" = "#33CC33"    # Marine fringe
)

# ── VM0033 depth intervals ────────────────────────────────────────────────────
# Standard depth midpoints (cm) for harmonization
VM0033_DEPTH_MIDPOINTS <- c(7.5, 22.5, 40, 75)

# Full interval table (used for thickness in carbon stock calculation)
VM0033_DEPTH_INTERVALS <- data.frame(
  depth_top      = c( 0, 15, 30,  50),
  depth_bottom   = c(15, 30, 50, 100),
  depth_midpoint = c(7.5, 22.5, 40, 75),
  thickness_cm   = c(15, 15, 20, 50)
)

# ── Bulk density defaults ─────────────────────────────────────────────────────
# Applied where bulk_density_g_cm3 is missing in the sample data (g/cm³).
# Based on literature values for BC coastal wetland ecosystems.
BD_DEFAULTS <- list(
  "IM" = 0.8,   # Intertidal marsh
  "NM" = 0.8,   # Near-margin
  "MF" = 0.8    # Marine fringe
)

# ── QC thresholds ─────────────────────────────────────────────────────────────
QC_SOC_MIN <- 0      # g/kg — minimum valid SOC
QC_SOC_MAX <- 500    # g/kg — maximum valid SOC (adjust for your ecosystem)
QC_BD_MIN  <- 0.1    # g/cm³ — minimum valid bulk density
QC_BD_MAX  <- 3.0    # g/cm³ — maximum valid bulk density

# ── Remote sensing band labels ────────────────────────────────────────────────
# Maps GEE raster band names to human-readable labels used in the RF importance
# plot. Bands not listed here will display their raw name.
#
# Literature context:
#   NIR/SWIR bands and NDVI are the strongest predictors of above-ground biomass
#   and SOC in tidal wetlands (Byrd et al. 2018; Macreadie et al. 2021).
#   SAR backscatter (VV/VH) captures canopy structure and inundation state
#   independent of cloud cover (Kasischke & Bourgeau-Chavez 1997).
#   Tidal inundation frequency is a primary control on carbon burial rate
#   in coastal wetlands (Chmura et al. 2003).
BAND_LABELS <- c(
  # Landsat 8/9 surface reflectance bands
  "SR_B1" = "Coastal aerosol (443 nm) — atmospheric/water",
  "SR_B2" = "Blue (482 nm) — water penetration",
  "SR_B3" = "Green (562 nm) — chlorophyll, turbidity",
  "SR_B4" = "Red (655 nm) — chlorophyll absorption",
  "SR_B5" = "NIR (865 nm) — vegetation biomass ★",
  "SR_B6" = "SWIR-1 (1609 nm) — soil & canopy moisture ★",
  "SR_B7" = "SWIR-2 (2201 nm) — soil organic matter ★",

  # Spectral indices
  "NDVI"  = "NDVI — live green vegetation density ★",
  "EVI"   = "EVI — canopy structure (atmosphere-corrected)",
  "SAVI"  = "SAVI — vegetation with soil-brightness correction",
  "NDWI"  = "NDWI — open water & canopy water content",
  "MNDWI" = "MNDWI — surface water extent",
  "NDMI"  = "NDMI — plant water stress & moisture ★",
  "BSI"   = "BSI — bare soil / unvegetated mudflat",
  "LSWI"  = "LSWI — land surface water index",

  # SAR Sentinel-1 (cloud-independent)
  "VV"         = "SAR VV — surface roughness & flooding ★",
  "VH"         = "SAR VH — vegetation structure & biomass ★",
  "VV_VH"      = "SAR VV/VH ratio — vegetation vs. open water",
  "HH"         = "SAR HH — surface scattering",
  "HV"         = "SAR HV — volume scattering / biomass",

  # Topographic & hydrological
  "elevation"  = "Elevation (m) — inundation risk & carbon depth ★",
  "slope"      = "Slope (°) — drainage & sediment transport",
  "aspect"     = "Aspect — solar exposure",
  "TPI"        = "Topographic position index — landform context",
  "TWI"        = "Topographic wetness index — soil moisture",

  # Tidal & coastal
  "tidal_freq"     = "Tidal inundation frequency — carbon burial rate ★",
  "dist_to_water"  = "Distance to water channel (m)",
  "dist_to_coast"  = "Distance to coastline (m)",

  # Climate
  "precip"    = "Mean annual precipitation (mm)",
  "temp_mean" = "Mean annual temperature (°C)"
)
