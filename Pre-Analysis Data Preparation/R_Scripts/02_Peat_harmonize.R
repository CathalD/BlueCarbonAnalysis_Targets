# =============================================================================
# SCRIPT 2: Canadian Peat Database — Harmonize Cores + Profiles → GEE Export
# =============================================================================
# Inputs:
#   - peat_cores.csv    : whole-core metadata (one row per CORE_ID)
#   - peat_profiles.csv : discrete depth intervals (one or more rows per CORE_ID)
#
# lower_depth = UPPER_SAMP_DEPTH + SAMP_THICK
#
# Outputs:
#   - peat_layers.csv   : one row per CORE_ID × SAMPLE_NO (layer-level)
#   - peat_profiles.csv : one row per CORE_ID (profile-level summary for GEE)
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# 0. USER SETTINGS
# -----------------------------------------------------------------------------

PATH_CORES    <- "CANPEAT_DATA/CANPEAT_CORES_DATA"
PATH_PROFILES <- "CANPEAT_DATA/CANPEAT_PROFILES_DATA"
PATH_OUT_LAYERS   <- "output/peat_layers.csv"
PATH_OUT_PROFILES <- "output/peat_profiles_gee.csv"

# -----------------------------------------------------------------------------
# 1. INGEST
# -----------------------------------------------------------------------------
message("--- Reading peat CSVs ---")

cores <- read_csv(PATH_CORES, show_col_types = FALSE, locale = locale(encoding = "UTF-8")) %>%
  # Remove junk BOM index column if present
  select(-any_of("...1")) %>%
  # Clean BOM artifact from CORE_ID column name (shows as "\xef..CORE_ID")
  rename_with(~ "CORE_ID", matches("CORE_ID"))

profs <- read_csv(PATH_PROFILES, show_col_types = FALSE, locale = locale(encoding = "UTF-8")) %>%
  select(-any_of("...1")) %>%
  rename_with(~ "CORE_ID", matches("CORE_ID"))

message(sprintf("Cores loaded:    %d rows | %d unique CORE_ID",
                nrow(cores), n_distinct(cores$CORE_ID)))
message(sprintf("Profiles loaded: %d rows | %d unique CORE_ID",
                nrow(profs), n_distinct(profs$CORE_ID)))

# -----------------------------------------------------------------------------
# 2. DATA AUDIT
# -----------------------------------------------------------------------------
message("\n--- Data audit ---")

# Check for CORE_IDs in profiles not present in cores (orphaned samples)
orphan_ids <- setdiff(profs$CORE_ID, cores$CORE_ID)
if (length(orphan_ids) > 0) {
  warning(sprintf("%d CORE_IDs in profiles have no matching core record: %s",
                  length(orphan_ids), paste(head(orphan_ids, 10), collapse = ", ")))
} else {
  message("OK: All profile CORE_IDs have a matching core record.")
}

# Check for CORE_IDs in cores with no profile data
no_profile <- setdiff(cores$CORE_ID, profs$CORE_ID)
message(sprintf("Cores with no profile data: %d", length(no_profile)))

# Duplicate SAMPLE_NO within a CORE_ID
dup_samples <- profs %>%
  group_by(CORE_ID, SAMPLE_NO) %>%
  filter(n() > 1)
if (nrow(dup_samples) > 0) {
  warning(sprintf("Duplicate CORE_ID × SAMPLE_NO combinations: %d rows", nrow(dup_samples)))
} else {
  message("OK: No duplicate CORE_ID × SAMPLE_NO.")
}

# NA summary for key variables
cat("\nKey variable NA counts in profiles:\n")
key_vars <- c("UPPER_SAMP_DEPTH", "SAMP_THICK", "BULK_DENSITY", "C_TOT_PCT",
              "C_ORG_PCT", "SAMP_CARB_MGHA")
print(colSums(is.na(profs[key_vars])))

# -----------------------------------------------------------------------------
# 3. RENAME & CALCULATE DERIVED COLUMNS IN PROFILES
# -----------------------------------------------------------------------------
message("\n--- Renaming and computing derived columns ---")

profs_clean <- profs %>%
  rename(
    upper_depth  = UPPER_SAMP_DEPTH,   # cm
    layer_thick  = SAMP_THICK,         # cm
    BDOD         = BULK_DENSITY,       # g/cm³  (BD_MEAS_EST kept as flag)
    TOTC         = C_TOT_PCT,          # %
    OrgC         = C_ORG_PCT,          # %
    TOTC_Stock   = SAMP_CARB_MGHA      # Mg C/ha (pre-calculated)
  ) %>%
  mutate(
    # Derive lower depth
    lower_depth       = upper_depth + layer_thick,          # cm
    layer_thickness_cm = layer_thick,                        # alias for clarity
    dataset           = "Peat Database"
  )

# -----------------------------------------------------------------------------
# 4. RENAME CORE METADATA
# -----------------------------------------------------------------------------
message("--- Renaming core metadata ---")

cores_clean <- cores %>%
  rename(
    Latitude   = LATITUDE,
    Longitude  = LONGITUDE,
    Year       = SAMPLING_YR,
    Org_C_MGHA = ORG_C_MGHA      # whole-core organic C (Mg C/ha)
  ) %>%
  select(CORE_ID, Latitude, Longitude, Year, Org_C_MGHA,
         CWCS_CLASS, PROV_TERR, ORG_DEPTH, PERMAFROST,
         SOURCE_SITEID, SOURCE_SUB)

# -----------------------------------------------------------------------------
# 5. MERGE PROFILES + CORES ON CORE_ID
# -----------------------------------------------------------------------------
message("--- Merging profiles with core metadata ---")

layers <- profs_clean %>%
  left_join(cores_clean, by = "CORE_ID") %>%
  # Reorder output columns cleanly
  select(
    dataset,
    CORE_ID,
    SAMPLE_NO,
    Latitude,
    Longitude,
    Year,
    PROV_TERR,
    upper_depth,
    lower_depth,
    layer_thickness_cm,
    BDOD,              # g/cm³
    BD_MEAS_EST,       # measured vs estimated flag
    TOTC,              # C_TOT_PCT (%)
    OrgC,              # C_ORG_PCT (%)
    TOTC_Stock,        # SAMP_CARB_MGHA (Mg C/ha) — pre-calculated
    ASH,
    VON_POST,
    MATERIAL_1, MATERIAL_2, MATERIAL_3,
    CSSC_HORIZON,
    SAMP_OM_CSSC,
    SAMP_PH,
    Org_C_MGHA,        # whole-core value from cores table
    CWCS_CLASS,
    ORG_DEPTH,
    PERMAFROST,
    SOURCE_SITEID
  )

message(sprintf("Layer table: %d rows | %d unique CORE_IDs",
                nrow(layers), n_distinct(layers$CORE_ID)))

write_csv(layers, PATH_OUT_LAYERS)
message(sprintf("Saved: %s", PATH_OUT_LAYERS))

# -----------------------------------------------------------------------------
# 6. PROFILE-LEVEL SUMMARY (one row per CORE_ID — for GEE import)
# -----------------------------------------------------------------------------
message("\n--- Building GEE profile summary ---")

profiles_gee <- layers %>%
  group_by(CORE_ID, dataset) %>%
  summarise(
    Latitude        = first(Latitude),
    Longitude       = first(Longitude),
    Year            = first(Year),
    PROV_TERR       = first(PROV_TERR),
    CWCS_CLASS      = first(CWCS_CLASS),
    ORG_DEPTH       = first(ORG_DEPTH),
    PERMAFROST      = first(PERMAFROST),
    Org_C_MGHA      = first(Org_C_MGHA),   # whole-core value
    n_samples       = n(),
    
    # Total depth = sum of all layer thicknesses (cm)
    Total_depth_cm  = sum(layer_thickness_cm, na.rm = TRUE),
    
    # Total_C = sum of per-layer pre-calculated TOTC_Stock (Mg C/ha)
    # NA if ALL layers are NA; partial sums use na.rm = TRUE
    Total_C_MgCha   = if (all(is.na(TOTC_Stock))) NA_real_
    else sum(TOTC_Stock, na.rm = TRUE),
    
    # Mean BDOD across layers (g/cm³)
    mean_BDOD       = mean(BDOD, na.rm = TRUE),
    
    # Mean TOTC and OrgC (%)
    mean_TOTC_pct   = mean(TOTC, na.rm = TRUE),
    mean_OrgC_pct   = mean(OrgC, na.rm = TRUE),
    
    .groups = "drop"
  )

cat("\nGEE profile summary preview:\n")
print(head(profiles_gee, 10))
cat(sprintf("\nTotal profiles for GEE: %d\n", nrow(profiles_gee)))

write_csv(profiles_gee, PATH_OUT_PROFILES)
message(sprintf("Saved: %s", PATH_OUT_PROFILES))

message("\n=== Peat harmonization complete ===")