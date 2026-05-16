# =============================================================================
# SCRIPT 5: Geographic + Depth Filtering of Harmonized Soil Profiles
# =============================================================================
# Inputs:
#   - combined_profiles.csv : 17,097 profiles (WOSIS 2023 + Peat DB + Janousek)
#   - combined_layers.csv   : 685,173 layers (all three datasets)
#
# Filters applied:
#   1. GEOGRAPHIC — keep profiles within temperate/boreal Northern Hemisphere:
#        latitude  >= 30°N
#        longitude >= -170°E  (western Alaska limit)
#        longitude <=  +60°E  (excludes Central Asia, most of Middle East/Africa)
#      Rationale: retains North America (incl. US & Canada), Europe, and
#      western Russia as climate analogs for Canadian carbon systems.
#      Removes: Philippines, Zimbabwe, Pakistan, India, Indonesia, Puerto Rico
#               (~18°N), tropical Americas, SE Asia.
#
#   2. DEPTH — keep profiles with total_depth_cm >= 30
#      Rationale: profiles shallower than 30 cm lack sufficient depth to
#      represent meaningful carbon stocks for transfer learning.
#
#   Both filters are applied to profiles first; the filtered profile_id set is
#   then used to remove corresponding layers from combined_layers.csv.
#
# Outputs:
#   - combined_profiles_filtered.csv : profiles passing both filters
#   - combined_layers_filtered.csv   : layers whose profile_id passed both filters
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# 0. USER SETTINGS — update paths as needed
# -----------------------------------------------------------------------------

PATH_PROFILES     <- "output/combined_profiles.csv"
PATH_LAYERS       <- "output/combined_layers.csv"
PATH_OUT_PROFILES <- "output/combined_profiles_filtered.csv"
PATH_OUT_LAYERS   <- "output/combined_layers_filtered.csv"

# ── Geographic filter bounds ─────────────────────────────────────
# Keeps temperate/boreal Northern Hemisphere: N. America, Europe, W. Russia
GEO_LAT_MIN  <-  30    # degrees N
GEO_LON_MIN  <- -170   # degrees E (western Alaska)
GEO_LON_MAX  <-  60    # degrees E (excludes Central Asia / Middle East east of Urals)

# ── Depth filter ──────────────────────────────────────────────────
DEPTH_MIN_CM <- 30     # minimum total profile depth (cm)

# -----------------------------------------------------------------------------
# 1. LOAD INPUTS
# -----------------------------------------------------------------------------
message("--- Loading combined outputs ---")

profiles <- read_csv(PATH_PROFILES, show_col_types = FALSE)
layers   <- read_csv(PATH_LAYERS,   show_col_types = FALSE)

message(sprintf("Profiles loaded: %d rows | %d unique profile_ids",
                nrow(profiles), n_distinct(profiles$profile_id)))
message(sprintf("Layers   loaded: %d rows | %d unique profile_ids",
                nrow(layers), n_distinct(layers$profile_id)))

# -----------------------------------------------------------------------------
# 2. PRE-FILTER AUDIT
# -----------------------------------------------------------------------------
message("\n--- Pre-filter audit ---")

cat("\nProfile counts by dataset:\n")
print(count(profiles, dataset))

cat("\nDepth distribution (total_depth_cm):\n")
print(summary(profiles$total_depth_cm))
cat(sprintf("  Profiles with total_depth_cm < 30 cm: %d\n",
            sum(profiles$total_depth_cm < DEPTH_MIN_CM, na.rm = TRUE)))
cat(sprintf("  Profiles with total_depth_cm = NA:    %d\n",
            sum(is.na(profiles$total_depth_cm))))

cat("\nLatitude range by dataset:\n")
profiles %>%
  group_by(dataset) %>%
  summarise(lat_min = min(latitude, na.rm = TRUE),
            lat_max = max(latitude, na.rm = TRUE),
            lon_min = min(longitude, na.rm = TRUE),
            lon_max = max(longitude, na.rm = TRUE),
            .groups = "drop") %>%
  print()

# Profiles outside the geographic bounds (pre-filter preview)
geo_out <- profiles %>%
  filter(latitude < GEO_LAT_MIN |
           longitude < GEO_LON_MIN |
           longitude > GEO_LON_MAX)

cat(sprintf("\nProfiles outside geographic bounds (preview — will be removed): %d\n",
            nrow(geo_out)))
cat("By dataset:\n")
print(count(geo_out, dataset))
cat("WOSIS: top removed countries:\n")
geo_out %>%
  filter(dataset == "WOSIS 2023") %>%
  count(country_name, sort = TRUE) %>%
  head(15) %>%
  print()

# -----------------------------------------------------------------------------
# 3. GEOGRAPHIC FILTER
# -----------------------------------------------------------------------------
message("\n--- Applying geographic filter ---")
message(sprintf("  Bounds: lat >= %g°N | lon >= %g°E | lon <= %g°E",
                GEO_LAT_MIN, GEO_LON_MIN, GEO_LON_MAX))

profiles_geo <- profiles %>%
  filter(
    latitude  >= GEO_LAT_MIN,
    longitude >= GEO_LON_MIN,
    longitude <= GEO_LON_MAX
  )

n_removed_geo <- nrow(profiles) - nrow(profiles_geo)
message(sprintf("  Removed by geographic filter: %d profiles", n_removed_geo))
message(sprintf("  Remaining after geographic filter: %d profiles", nrow(profiles_geo)))

cat("\nRemaining by dataset after geographic filter:\n")
print(count(profiles_geo, dataset))

# -----------------------------------------------------------------------------
# 4. DEPTH FILTER
# -----------------------------------------------------------------------------
message("\n--- Applying depth filter (total_depth_cm >= 30 cm) ---")

profiles_filt <- profiles_geo %>%
  filter(total_depth_cm >= DEPTH_MIN_CM | is.na(total_depth_cm))
# Note: profiles with NA total_depth_cm are retained (conservative — don't remove
# if depth is unknown; they can be excluded later during modelling if needed)

n_removed_depth <- nrow(profiles_geo) - nrow(profiles_filt)
message(sprintf("  Removed by depth filter: %d profiles", n_removed_depth))
message(sprintf("  Remaining after depth filter: %d profiles", nrow(profiles_filt)))

cat("\nFinal profile counts by dataset:\n")
print(count(profiles_filt, dataset))

# Summarise what was removed at each stage
cat("\n── Filter summary ──────────────────────────────────────────────\n")
cat(sprintf("  Started with:              %6d profiles\n", nrow(profiles)))
cat(sprintf("  Removed (geographic):    - %6d profiles\n", n_removed_geo))
cat(sprintf("  Removed (depth < 30 cm): - %6d profiles\n", n_removed_depth))
cat(sprintf("  Final retained:            %6d profiles\n", nrow(profiles_filt)))
cat(sprintf("  Retention rate:            %.1f%%\n",
            100 * nrow(profiles_filt) / nrow(profiles)))

# Geographic range of retained profiles
cat("\nGeographic range of retained profiles:\n")
profiles_filt %>%
  summarise(lat_min = min(latitude, na.rm = TRUE),
            lat_max = max(latitude, na.rm = TRUE),
            lon_min = min(longitude, na.rm = TRUE),
            lon_max = max(longitude, na.rm = TRUE)) %>%
  print()

# -----------------------------------------------------------------------------
# 5. FILTER LAYERS TO MATCH FILTERED PROFILES
# -----------------------------------------------------------------------------
message("\n--- Filtering layers to match retained profiles ---")

# Get the set of profile_ids that passed both filters
# profile_id type must match layers$profile_id (both are character from script 04)
keep_ids <- profiles_filt %>%
  mutate(profile_id = as.character(profile_id)) %>%
  pull(profile_id)

layers_filt <- layers %>%
  mutate(profile_id = as.character(profile_id)) %>%
  filter(profile_id %in% keep_ids)

message(sprintf("  Layers before filter: %d", nrow(layers)))
message(sprintf("  Layers after filter:  %d", nrow(layers_filt)))
message(sprintf("  Layers removed:       %d", nrow(layers) - nrow(layers_filt)))

# Verify all filtered profiles have at least one layer
profiles_no_layers <- profiles_filt %>%
  mutate(profile_id = as.character(profile_id)) %>%
  filter(!profile_id %in% layers_filt$profile_id)
if (nrow(profiles_no_layers) > 0) {
  warning(sprintf("%d filtered profiles have no corresponding layers",
                  nrow(profiles_no_layers)))
} else {
  message("OK: All retained profiles have at least one layer.")
}

# Layer counts by dataset
cat("\nLayer counts by dataset after filtering:\n")
print(count(layers_filt, dataset))

# -----------------------------------------------------------------------------
# 6. POST-FILTER AUDIT
# -----------------------------------------------------------------------------
message("\n--- Post-filter QC ---")

# Depth range of retained profiles
cat("Depth range (total_depth_cm) of retained profiles:\n")
print(summary(profiles_filt$total_depth_cm))

# Check no profile_ids in layers that aren't in profiles
orphan_layer_ids <- setdiff(layers_filt$profile_id, as.character(profiles_filt$profile_id))
if (length(orphan_layer_ids) > 0) {
  warning(sprintf("WARNING: %d profile_ids in filtered layers have no matching filtered profile",
                  length(orphan_layer_ids)))
} else {
  message("OK: All filtered layer profile_ids match filtered profiles.")
}

# NA check on key columns
cat("\nNA counts in filtered profiles (key columns):\n")
profiles_filt %>%
  summarise(across(c(latitude, longitude, total_depth_cm, sum_OrgC_Stock_kgm2,
                     mean_BDOD, mean_OrgC_pct), ~ sum(is.na(.)))) %>%
  print()

# -----------------------------------------------------------------------------
# 7. EXPORT
# -----------------------------------------------------------------------------
message("\n--- Exporting filtered outputs ---")

write_csv(profiles_filt, PATH_OUT_PROFILES)
message(sprintf("Saved: %s  (%d profiles)", PATH_OUT_PROFILES, nrow(profiles_filt)))

write_csv(layers_filt,   PATH_OUT_LAYERS)
message(sprintf("Saved: %s  (%d layers)",   PATH_OUT_LAYERS,   nrow(layers_filt)))

message("\n=== Geographic + depth filtering complete ===")
message(sprintf("  combined_profiles_filtered.csv: %d profiles", nrow(profiles_filt)))
message(sprintf("  combined_layers_filtered.csv:   %d layers",   nrow(layers_filt)))
message("")
message("Next step: upload combined_profiles_filtered.csv to the Colab notebook")
message("  CoastalBlueCarbon_GlobalCoreCovariate_Extraction.ipynb")
message("  Climate (MAT/MAP) filtering will be applied there using TerraClimate.")
