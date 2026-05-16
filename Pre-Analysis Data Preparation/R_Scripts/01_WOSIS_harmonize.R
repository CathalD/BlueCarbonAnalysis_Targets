# =============================================================================
# SCRIPT 1: WoSIS 2023 — Harmonize BDOD, OrgC, TOTC → Carbon Stocks
# =============================================================================
# Inputs (separate CSV files per variable):
#   - wosis_bdod.csv   : BDOD (kg/dm³)
#   - wosis_orgc.csv   : OrgC (g/kg)
#   - wosis_totc.csv   : TOTC (g/kg)
#
# Stock formula: Stock (kg/m²) = C (g/kg) × BDOD (kg/dm³) × thickness (cm) / 10
#
# Outputs:
#   - wosis_layers.csv   : one row per layer_id with all variables + stocks
#   - wosis_profiles.csv : one row per profile_id with summed stocks + total depth
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# 0. USER SETTINGS — update file paths as needed
# -----------------------------------------------------------------------------

PATH_BDOD  <- "WOSIS_DATA/WOSIS_BDOD_DATA"
PATH_ORGC  <- "WOSIS_DATA/WOSIS_ORGC_DATA"
PATH_TOTC  <- "WOSIS_DATA/WOSIS_TOTC_DATA"
PATH_OUT_LAYERS   <- "output/wosis_layers.csv"
PATH_OUT_PROFILES <- "output/wosis_profiles.csv"

# Key columns expected in every input file — updated to match actual WoSIS headers
KEY_COLS <- c("profile_id", "layer_id", "profile_code", "layer_name",
              "upper_depth", "lower_depth", "organic_surface",
              "longitude", "latitude", "country_name", "region",
              "continent", "date", "dataset_id")

# -----------------------------------------------------------------------------
# 1. INGEST
# -----------------------------------------------------------------------------
message("--- Reading raw CSVs ---")

read_wosis <- function(path, value_col) {
  # WoSIS raw files are comma-separated despite the .tsv extension
  df <- read_csv(path, show_col_types = FALSE)
  
  # Rename value_avg to the target variable name
  if (!"value_avg" %in% names(df)) stop(paste("'value_avg' column missing in:", path))
  df <- df %>% rename(!!value_col := value_avg)
  
  # Keep only the columns we need
  keep <- c(KEY_COLS, value_col)
  missing <- setdiff(keep, names(df))
  if (length(missing) > 0) stop(paste("Missing columns in", path, ":", paste(missing, collapse = ", ")))
  
  df %>% select(all_of(keep))
}

bdod <- read_wosis(PATH_BDOD, "BDOD")   # kg/dm³
orgc <- read_wosis(PATH_ORGC, "OrgC")   # g/kg
totc <- read_wosis(PATH_TOTC, "TOTC")   # g/kg

# -----------------------------------------------------------------------------
# 2. DATA AUDIT
# -----------------------------------------------------------------------------
message("--- Data audit ---")

audit_df <- function(df, name) {
  cat(sprintf("\n[%s] rows: %d | unique layer_id: %d | unique profile_id: %d\n",
              name, nrow(df), n_distinct(df$layer_id), n_distinct(df$profile_id)))
  cat("  NA summary:\n")
  print(colSums(is.na(df)))
  
  # Flag duplicate layer_ids
  dups <- df %>% filter(duplicated(layer_id) | duplicated(layer_id, fromLast = TRUE))
  if (nrow(dups) > 0) {
    cat(sprintf("  WARNING: %d rows with duplicate layer_id detected!\n", nrow(dups)))
  } else {
    cat("  OK: No duplicate layer_id.\n")
  }
}

audit_df(bdod, "BDOD")
audit_df(orgc, "OrgC")
audit_df(totc, "TOTC")

# -----------------------------------------------------------------------------
# 3. MERGE BY layer_id (full outer join to preserve all layers)
# -----------------------------------------------------------------------------
message("\n--- Merging by layer_id ---")

# Start from BDOD as the spine; outer join OrgC then TOTC
# KEY_COLS from BDOD are used as the reference geometry
layers <- bdod %>%
  full_join(orgc %>% select(layer_id, OrgC), by = "layer_id") %>%
  full_join(totc %>% select(layer_id, TOTC), by = "layer_id")

# Ensure profile_id is present for every layer
missing_profile <- layers %>% filter(is.na(profile_id))
if (nrow(missing_profile) > 0) {
  warning(sprintf("%d layer_id rows are missing a profile_id — inspect 'missing_profile_layers.csv'",
                  nrow(missing_profile)))
  write_csv(missing_profile, "output/missing_profile_layers.csv")
}

# Verify no duplicate layer_ids in merged result
if (any(duplicated(layers$layer_id))) {
  warning("Duplicate layer_ids found after merge — check input files for overlapping rows.")
}

message(sprintf("Merged layer table: %d rows, %d unique profiles",
                nrow(layers), n_distinct(layers$profile_id, na.rm = TRUE)))

# -----------------------------------------------------------------------------
# 4. CALCULATE LAYER THICKNESS & CARBON STOCKS
# -----------------------------------------------------------------------------
message("--- Calculating carbon stocks ---")

# layer_thickness in cm
layers <- layers %>%
  mutate(
    layer_thickness_cm = lower_depth - upper_depth,
    
    # Stock (kg/m²) = C (g/kg) × BDOD (kg/dm³) × thickness (cm) / 10
    # If either C or BDOD is NA, stock is NA
    OrgC_Stock  = ifelse(!is.na(OrgC)  & !is.na(BDOD) & !is.na(layer_thickness_cm),
                         OrgC  * BDOD * layer_thickness_cm / 10,
                         NA_real_),
    
    TOTC_Stock  = ifelse(!is.na(TOTC)  & !is.na(BDOD) & !is.na(layer_thickness_cm),
                         TOTC  * BDOD * layer_thickness_cm / 10,
                         NA_real_),
    
    # Dataset label
    dataset = "WOSIS 2023"
  )

# Stock summary
cat("\nOrgC_Stock (kg/m²) summary:\n"); print(summary(layers$OrgC_Stock))
cat("\nTOTC_Stock (kg/m²) summary:\n"); print(summary(layers$TOTC_Stock))

# -----------------------------------------------------------------------------
# 5. EXPORT LAYER-LEVEL CSV
# -----------------------------------------------------------------------------
message("--- Exporting layer CSV ---")

layers_out <- layers %>%
  select(dataset, profile_id, layer_id, profile_code, layer_name,
         longitude, latitude, country_name, region, continent, date, dataset_id,
         upper_depth, lower_depth, layer_thickness_cm, organic_surface,
         BDOD, OrgC, TOTC, OrgC_Stock, TOTC_Stock)

write_csv(layers_out, PATH_OUT_LAYERS)
message(sprintf("Saved: %s  (%d rows)", PATH_OUT_LAYERS, nrow(layers_out)))

# -----------------------------------------------------------------------------
# 6. PROFILE-LEVEL SUMMARY (one row per profile_id)
# -----------------------------------------------------------------------------
message("--- Building profile summary ---")

profiles <- layers %>%
  filter(!is.na(profile_id)) %>%
  group_by(profile_id, dataset) %>%
  summarise(
    longitude         = first(longitude),
    latitude          = first(latitude),
    country_name      = first(country_name),
    region            = first(region),
    continent         = first(continent),
    date              = first(date),
    total_depth_cm    = sum(layer_thickness_cm, na.rm = TRUE),
    n_layers          = n(),
    n_layers_orgc     = sum(!is.na(OrgC_Stock)),
    n_layers_totc     = sum(!is.na(TOTC_Stock)),
    sum_OrgC_Stock    = if (all(is.na(OrgC_Stock))) NA_real_ else sum(OrgC_Stock, na.rm = TRUE),
    sum_TOTC_Stock    = if (all(is.na(TOTC_Stock))) NA_real_ else sum(TOTC_Stock, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nProfile summary preview:\n")
print(head(profiles, 10))

write_csv(profiles, PATH_OUT_PROFILES)
message(sprintf("Saved: %s  (%d profiles)", PATH_OUT_PROFILES, nrow(profiles)))

# -----------------------------------------------------------------------------
# 7. CANADA-FILTERED EXPORTS
# -----------------------------------------------------------------------------
message("--- Exporting Canada-filtered CSVs ---")

layers_canada <- layers_out %>%
  filter(country_name == "Canada")

profiles_canada <- profiles %>%
  filter(country_name == "Canada")

message(sprintf("Canada layers:   %d rows", nrow(layers_canada)))
message(sprintf("Canada profiles: %d rows", nrow(profiles_canada)))

write_csv(layers_canada,   "output/wosis_layers_canada.csv")
write_csv(profiles_canada, "output/wosis_profiles_canada.csv")

message("Saved: output/wosis_layers_canada.csv")
message("Saved: output/wosis_profiles_canada.csv")

message("\n=== WoSIS harmonization complete ===")

