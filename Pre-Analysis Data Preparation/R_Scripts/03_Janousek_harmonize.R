# =============================================================================
# SCRIPT 3: Janousek Coastal Wetland Database — Harmonize Layers + Profiles
# =============================================================================
# Inputs:
#   - Global_Core_Samples.csv   : one row per subsample (layer) per core
#   - Global_Core_Locations.csv : one row per core (location metadata)
#
# OrgC stock formula:
#   OrgC_Stock_kgm2 = OrgC_pct * bulk_density(g/cm³) * sample_length(cm) / 10
#
# Note: soc_percent is the average of perc_C_C and perc_C_OM.
#       Some rows contain "#VALUE!" in soc_percent due to Excel formula errors.
#       Fallback: use perc_C_C where soc_percent is invalid. Rows with neither
#       valid will have OrgC_pct = NA.
#
# Note: Global_Core_Locations.csv contains one header-artifact row
#       (core_id == "SampID") which is filtered before joining.
#
# Outputs:
#   - janousek_layers.csv   : one row per subsample (layer)
#   - janousek_profiles.csv : one row per core_id
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# 0. USER SETTINGS — update file paths as needed
# -----------------------------------------------------------------------------

PATH_SAMPLES   <- "JANOUSEK_DATA/Global_Core_Samples.csv"
PATH_LOCATIONS <- "JANOUSEK_DATA/Global_Core_Locations.csv"
PATH_OUT_LAYERS   <- "output/janousek_layers.csv"
PATH_OUT_PROFILES <- "output/janousek_profiles.csv"

# -----------------------------------------------------------------------------
# 1. INGEST
# -----------------------------------------------------------------------------
message("--- Reading Janousek CSVs ---")

# Read carbon columns as character to handle "#VALUE!" strings from Excel
samples <- read_csv(
  PATH_SAMPLES,
  show_col_types = FALSE,
  na = c("NA", ""),
  col_types = cols(
    soc_percent = col_character(),
    perc_C_OM   = col_character(),
    perc_C_C    = col_character(),
    .default    = col_guess()
  )
)

locations <- read_csv(PATH_LOCATIONS, show_col_types = FALSE, na = c("NA", ""))

message(sprintf("Samples loaded:   %d rows | %d unique core_id",
                nrow(samples), n_distinct(samples$core_id)))
message(sprintf("Locations loaded: %d rows | %d unique core_id",
                nrow(locations), n_distinct(locations$core_id)))

# -----------------------------------------------------------------------------
# 2. DATA AUDIT
# -----------------------------------------------------------------------------
message("\n--- Data audit ---")

# Remove header-artifact row from locations (core_id == "SampID")
artifact_rows <- locations %>% filter(core_id == "SampID")
if (nrow(artifact_rows) > 0) {
  message(sprintf("Removing %d header artifact row(s) from locations (core_id = 'SampID')",
                  nrow(artifact_rows)))
  locations <- locations %>% filter(core_id != "SampID")
}

message(sprintf("Locations after artifact removal: %d rows | %d unique core_id",
                nrow(locations), n_distinct(locations$core_id)))

# Check for orphaned sample core_ids (in samples but not in locations)
orphan_ids <- setdiff(samples$core_id, locations$core_id)
if (length(orphan_ids) > 0) {
  warning(sprintf("%d core_ids in samples have no matching location record: %s",
                  length(orphan_ids),
                  paste(head(orphan_ids, 10), collapse = ", ")))
} else {
  message("OK: All sample core_ids have a matching location record.")
}

# Count #VALUE! entries in soc_percent
n_value_error <- sum(samples$soc_percent == "#VALUE!", na.rm = TRUE)
n_soc_na      <- sum(is.na(samples$soc_percent))
message(sprintf("soc_percent '#VALUE!' entries: %d | NA: %d (of %d rows)",
                n_value_error, n_soc_na, nrow(samples)))

# NA summary for key measurement columns (bulk_density, sample_length already numeric)
cat("\nKey variable NA counts in samples:\n")
cat(sprintf("  bulk_density   NAs: %d\n", sum(is.na(samples$bulk_density))))
cat(sprintf("  sample_length  NAs: %d\n", sum(is.na(samples$sample_length))))
cat(sprintf("  depth_min      NAs: %d\n", sum(is.na(samples$depth_min))))
cat(sprintf("  depth_max      NAs: %d\n", sum(is.na(samples$depth_max))))

# Flag implausible bulk_density values (> 2.5 g/cm³ is unusual for wetland soils)
n_high_bd <- sum(!is.na(samples$bulk_density) & samples$bulk_density > 2.5)
if (n_high_bd > 0) {
  message(sprintf("WARNING: %d layers with bulk_density > 2.5 g/cm³ — may indicate mineral sediment or unit error",
                  n_high_bd))
}

# Flag implausible sample_length (> 100 cm is very thick for a single interval)
n_long_layer <- sum(!is.na(samples$sample_length) & samples$sample_length > 100)
if (n_long_layer > 0) {
  message(sprintf("WARNING: %d layers with sample_length > 100 cm — review these records",
                  n_long_layer))
}

# -----------------------------------------------------------------------------
# 3. CLEAN CARBON CONCENTRATION VALUES
# -----------------------------------------------------------------------------
message("\n--- Cleaning carbon concentration values ---")

# soc_percent = average(perc_C_OM, perc_C_C) — verified numerically
# Where soc_percent has "#VALUE!" (Excel error), fall back to perc_C_C
samples_clean <- samples %>%
  mutate(
    soc_numeric  = suppressWarnings(as.numeric(soc_percent)),
    perc_C_C_num = suppressWarnings(as.numeric(perc_C_C)),

    # OrgC_pct: prefer soc_percent; fall back to perc_C_C if soc is invalid
    OrgC_pct = case_when(
      !is.na(soc_numeric)    ~ soc_numeric,
      !is.na(perc_C_C_num)   ~ perc_C_C_num,
      TRUE                   ~ NA_real_
    )
  )

n_fallback    <- sum(!is.na(samples_clean$perc_C_C_num) & is.na(samples_clean$soc_numeric))
n_still_na    <- sum(is.na(samples_clean$OrgC_pct))
message(sprintf("Rows using perc_C_C fallback: %d", n_fallback))
message(sprintf("Rows where OrgC_pct is ultimately NA: %d", n_still_na))

cat("\nOrgC_pct summary (after cleaning):\n")
print(summary(samples_clean$OrgC_pct))

# -----------------------------------------------------------------------------
# 4. RENAME & DERIVE COLUMNS
# -----------------------------------------------------------------------------
message("\n--- Renaming and computing derived columns ---")

samples_renamed <- samples_clean %>%
  mutate(
    # Construct unique layer_id: core_id + SubSampID
    # SubSampID alone has duplicate core_id+SubSampID combos; compound key is unique
    layer_id          = paste(core_id, SubSampID, sep = "_"),

    # Depth columns
    upper_depth       = depth_min,        # cm
    lower_depth       = depth_max,        # cm
    layer_thickness_cm = sample_length,   # cm

    # Bulk density (already g/cm³)
    BDOD              = bulk_density,

    # Total carbon: not available in Janousek data
    TOTC_pct          = NA_real_,

    # Dataset label
    dataset           = "Janousek"
  )

# -----------------------------------------------------------------------------
# 5. CALCULATE CARBON STOCKS
# -----------------------------------------------------------------------------
message("--- Calculating carbon stocks ---")

# OrgC_Stock_kgm2 = OrgC(%) × BD(g/cm³) × thickness(cm) / 10
# Requires all three inputs to be non-NA
samples_renamed <- samples_renamed %>%
  mutate(
    OrgC_Stock_kgm2 = ifelse(
      !is.na(OrgC_pct) & !is.na(BDOD) & !is.na(layer_thickness_cm),
      OrgC_pct * BDOD * layer_thickness_cm / 10,
      NA_real_
    )
  )

cat("\nOrgC_Stock_kgm2 summary:\n")
print(summary(samples_renamed$OrgC_Stock_kgm2))
message(sprintf("Layers with calculable OrgC_Stock: %d of %d",
                sum(!is.na(samples_renamed$OrgC_Stock_kgm2)), nrow(samples_renamed)))

# -----------------------------------------------------------------------------
# 6. MERGE SAMPLES WITH LOCATION METADATA
# -----------------------------------------------------------------------------
message("\n--- Merging samples with location metadata ---")

locations_slim <- locations %>%
  mutate(core_id = as.character(core_id)) %>%
  select(core_id, latitude, longitude, study_id, ecosystem, state)

layers <- samples_renamed %>%
  mutate(core_id = as.character(core_id)) %>%
  left_join(locations_slim, by = "core_id")

message(sprintf("Merged layer table: %d rows | %d unique core_ids",
                nrow(layers), n_distinct(layers$core_id)))

# Verify no row multiplication from the join
if (nrow(layers) != nrow(samples_renamed)) {
  warning(sprintf("Row count changed after location join: %d → %d (check for duplicate core_ids in locations)",
                  nrow(samples_renamed), nrow(layers)))
} else {
  message("OK: Row count unchanged after location join.")
}

# -----------------------------------------------------------------------------
# 7. EXPORT LAYER-LEVEL CSV
# -----------------------------------------------------------------------------
message("\n--- Exporting janousek_layers.csv ---")

layers_out <- layers %>%
  select(
    dataset,
    profile_id        = core_id,
    layer_id,
    latitude,
    longitude,
    upper_depth,
    lower_depth,
    layer_thickness_cm,
    BDOD,
    OrgC_pct,
    TOTC_pct,
    OrgC_Stock_kgm2,
    country_name      = state,     # state/province code — closest available geographic field
    year              = StudyID    # StudyID is study number, not year; set year NA below
  ) %>%
  mutate(
    country_name = NA_character_,  # no country_name available
    year         = NA_integer_     # no sampling year available
  )

write_csv(layers_out, PATH_OUT_LAYERS)
message(sprintf("Saved: %s  (%d rows)", PATH_OUT_LAYERS, nrow(layers_out)))

# -----------------------------------------------------------------------------
# 8. PROFILE-LEVEL AGGREGATION
# -----------------------------------------------------------------------------
message("\n--- Building janousek_profiles.csv ---")

profiles <- layers %>%
  group_by(core_id, dataset) %>%
  summarise(
    latitude            = first(latitude),
    longitude           = first(longitude),
    total_depth_cm      = sum(layer_thickness_cm, na.rm = TRUE),
    n_layers            = n(),
    sum_OrgC_Stock_kgm2 = if (all(is.na(OrgC_Stock_kgm2))) NA_real_
                          else sum(OrgC_Stock_kgm2, na.rm = TRUE),
    mean_BDOD           = mean(BDOD, na.rm = TRUE),
    mean_OrgC_pct       = mean(OrgC_pct, na.rm = TRUE),
    country_name        = NA_character_,
    year                = NA_integer_,
    .groups = "drop"
  ) %>%
  rename(profile_id = core_id)

cat("\nProfile summary preview:\n")
print(head(profiles, 5))
cat(sprintf("\nTotal profiles: %d\n", nrow(profiles)))

cat("\nsum_OrgC_Stock_kgm2 summary:\n")
print(summary(profiles$sum_OrgC_Stock_kgm2))

# -----------------------------------------------------------------------------
# 9. EXPORT PROFILE-LEVEL CSV
# -----------------------------------------------------------------------------
write_csv(profiles, PATH_OUT_PROFILES)
message(sprintf("Saved: %s  (%d profiles)", PATH_OUT_PROFILES, nrow(profiles)))

message("\n=== Janousek harmonization complete ===")
