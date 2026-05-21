# =============================================================================
# R/preanalysis/global_data.R
# Janousek coastal wetland database: ingest, harmonize, and filter profiles
# =============================================================================

# ---------------------------------------------------------------------------
# ingest_janousek()
# ---------------------------------------------------------------------------
# Reads Global_Core_Locations.csv + Global_Core_Samples.csv, harmonizes
# layer-level data, aggregates to profile level.
#
# Returns list(layers = <data.frame>, profiles = <data.frame>)
# Profiles include `ecosystem` column for downstream filtering.
# ---------------------------------------------------------------------------
ingest_janousek <- function(locations_file, samples_file) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
  })

  message("[Janousek] Reading raw CSVs...")

  samples <- read_csv(
    samples_file,
    show_col_types = FALSE,
    na = c("NA", ""),
    col_types = cols(
      soc_percent = col_character(),
      perc_C_OM   = col_character(),
      perc_C_C    = col_character(),
      .default    = col_guess()
    )
  )

  locations <- read_csv(locations_file, show_col_types = FALSE, na = c("NA", ""))

  message(sprintf("[Janousek] Samples: %d rows | %d unique core_ids",
                  nrow(samples), n_distinct(samples$core_id)))
  message(sprintf("[Janousek] Locations: %d rows | %d unique core_ids",
                  nrow(locations), n_distinct(locations$core_id)))

  # Remove header-artifact row (core_id == "SampID")
  locations <- locations |> filter(core_id != "SampID")

  # Check for orphaned core_ids
  orphans <- setdiff(samples$core_id, locations$core_id)
  if (length(orphans) > 0)
    warning(sprintf("[Janousek] %d core_ids in samples have no location record",
                    length(orphans)))

  # ── Clean carbon concentration ───────────────────────────────────────────
  samples_clean <- samples |>
    mutate(
      soc_numeric  = suppressWarnings(as.numeric(soc_percent)),
      perc_C_C_num = suppressWarnings(as.numeric(perc_C_C)),
      OrgC_pct = case_when(
        !is.na(soc_numeric)  ~ soc_numeric,
        !is.na(perc_C_C_num) ~ perc_C_C_num,
        TRUE                 ~ NA_real_
      )
    )

  # ── Rename + derive columns ──────────────────────────────────────────────
  samples_renamed <- samples_clean |>
    mutate(
      layer_id           = paste(core_id, SubSampID, sep = "_"),
      upper_depth        = depth_min,
      lower_depth        = depth_max,
      layer_thickness_cm = sample_length,
      BDOD               = bulk_density,
      TOTC_pct           = NA_real_,
      dataset            = "Janousek"
    )

  # ── Carbon stock per layer ───────────────────────────────────────────────
  # OrgC_Stock_kgm2 = OrgC(%) × BD(g/cm³) × thickness(cm) / 10
  samples_renamed <- samples_renamed |>
    mutate(
      OrgC_Stock_kgm2 = ifelse(
        !is.na(OrgC_pct) & !is.na(BDOD) & !is.na(layer_thickness_cm),
        OrgC_pct * BDOD * layer_thickness_cm / 10,
        NA_real_
      )
    )

  # ── Join location metadata (keep ecosystem for filtering) ────────────────
  locations_slim <- locations |>
    mutate(core_id = as.character(core_id)) |>
    select(core_id, latitude, longitude, study_id, ecosystem, state)

  layers_joined <- samples_renamed |>
    mutate(core_id = as.character(core_id)) |>
    left_join(locations_slim, by = "core_id")

  # ── Layer-level output ───────────────────────────────────────────────────
  layers_out <- layers_joined |>
    select(
      dataset,
      profile_id         = core_id,
      layer_id,
      ecosystem,
      latitude,
      longitude,
      upper_depth,
      lower_depth,
      layer_thickness_cm,
      BDOD,
      OrgC_pct,
      TOTC_pct,
      OrgC_Stock_kgm2
    ) |>
    mutate(
      country_name = NA_character_,
      year         = NA_integer_
    )

  # ── Profile-level aggregation ────────────────────────────────────────────
  profiles_out <- layers_joined |>
    group_by(core_id, dataset) |>
    summarise(
      ecosystem           = first(ecosystem),
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
    ) |>
    rename(profile_id = core_id)

  n_missing_latlon <- sum(is.na(profiles_out$latitude) | is.na(profiles_out$longitude))
  message(sprintf("[Janousek] Harmonization complete: %d profiles | %d layers",
                  nrow(profiles_out), nrow(layers_out)))
  message(sprintf("[Janousek] Ecosystem breakdown: %s",
                  paste(names(table(profiles_out$ecosystem)),
                        table(profiles_out$ecosystem), sep = "=", collapse = ", ")))
  if (n_missing_latlon > 0)
    message(sprintf("[Janousek] %d profiles missing lat/lon (will be dropped before GEE)", n_missing_latlon))

  list(layers = layers_out, profiles = profiles_out)
}


# ---------------------------------------------------------------------------
# filter_for_gee()
# ---------------------------------------------------------------------------
# Filters Janousek profiles to specified coastal wetland ecosystem types,
# drops profiles with missing coordinates, and returns a clean data.frame
# ready for GEE covariate extraction.
#
# ecosystems : character vector of Janousek ecosystem codes to retain.
#   "EM" = Estuarine Marsh (saltmarsh / tidal marsh)
#   "SG" = Seagrass meadow
#   "MG" = Mangrove (add here if needed in future)
# ---------------------------------------------------------------------------
filter_for_gee <- function(janousek_harmonized, ecosystems = c("EM", "SG")) {
  suppressPackageStartupMessages(library(dplyr))

  profiles <- janousek_harmonized$profiles

  n_start <- nrow(profiles)

  filtered <- profiles |>
    filter(ecosystem %in% ecosystems) |>
    filter(!is.na(latitude), !is.na(longitude)) |>
    mutate(profile_id = as.character(profile_id))

  n_eco_removed  <- n_start - nrow(profiles |> filter(ecosystem %in% ecosystems))
  n_coord_removed <- nrow(profiles |> filter(ecosystem %in% ecosystems)) - nrow(filtered)

  message(sprintf("[filter_for_gee] Input: %d profiles", n_start))
  message(sprintf("[filter_for_gee] Removed %d non-target ecosystems (kept: %s)",
                  n_eco_removed, paste(ecosystems, collapse = ", ")))
  message(sprintf("[filter_for_gee] Removed %d profiles with missing lat/lon", n_coord_removed))
  message(sprintf("[filter_for_gee] Output: %d profiles ready for GEE extraction", nrow(filtered)))

  filtered
}
