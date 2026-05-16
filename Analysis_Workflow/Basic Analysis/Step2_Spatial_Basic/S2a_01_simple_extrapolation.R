# ============================================================================
# STEP 2a: SIMPLE AREA-WEIGHTED EXTRAPOLATION
# ============================================================================
# PURPOSE: Estimate aerial and total carbon stocks by multiplying mean carbon
#          stock per stratum (kg/m²) by stratum area from the AOI boundary.
#          This is the simplest spatial upscaling method and requires only
#          the AOI GeoJSON/shapefile with a stratum attribute column.
#
# INPUTS:
#   - data_processed/carbon_by_stratum_summary.csv  (from S1_01_data_prep.R)
#   - AOI_FILE (set in blue_carbon_config.R) — GeoJSON/SHP/GPKG with stratum column
#   - AOI_STRATUM_FIELD (set in config) — column name identifying strata in AOI
#
# OUTPUTS:
#   - outputs/Stage2a_SimpleExtrapolation/tables/simple_extrapolation_by_stratum.csv
#   - outputs/Stage2a_SimpleExtrapolation/tables/simple_extrapolation_aoi_total.csv
#   - outputs/Stage2a_SimpleExtrapolation/maps/simple_extrapolation_carbon_stock_mean_kg_m2.tif
#   - outputs/Stage2a_SimpleExtrapolation/plots/simple_extrapolation_carbon_stock_map.png
#   - outputs/Stage2a_SimpleExtrapolation/plots/simple_extrapolation_stratum_summary_chart.png
#
# UNITS:
#   - Aerial carbon stock:  kg/m²  (= mean_stock from S1_01, directly)
#   - Aerial carbon stock:  Mg C/ha (×10 conversion from kg/m²)
#   - Total carbon stock:   kg  (aerial kg/m² × stratum area m²)
# ============================================================================

# ── PATH RESOLVER ────────────────────────────────────────────────────────────
# Script is 4 levels deep: BlueCarbon_Workflow_V1.0/Analysis_Workflow/
#   Basic Analysis/Step2_Spatial_Basic/S2a_01_simple_extrapolation.R
local({
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) {
    setwd(dirname(dirname(dirname(dirname(this_file)))))
    return()
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NULL)
    if (!is.null(active) && nchar(active) > 0) {
      setwd(dirname(dirname(dirname(dirname(normalizePath(active))))))
      return()
    }
  }
  message("⚠ Could not auto-detect project root. Please run:\n",
          "  setwd('/path/to/BlueCarbon_Workflow_V1.0')\nbefore sourcing this script.")
})
# ─────────────────────────────────────────────────────────────────────────────

if (file.exists("Analysis_Workflow/blue_carbon_config.R")) {
  source("Analysis_Workflow/blue_carbon_config.R")
} else {
  stop("Configuration file not found. Run S1_00_setup_directories.R first.")
}

# ── GATE: requires AOI_FILE ──────────────────────────────────────────────────
if (!exists("AOI_FILE") || is.null(AOI_FILE) || !file.exists(AOI_FILE)) {
  cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("  STEP 2a SKIPPED — AOI_FILE not set or file not found.\n\n")
  cat("  To run simple extrapolation, set AOI_FILE in:\n")
  cat("    Analysis_Workflow/blue_carbon_config.R\n")
  cat("  The AOI must have a stratum attribute column.\n")
  cat("  Column name is set by AOI_STRATUM_FIELD in config (default: 'stratum')\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")
  invisible(NULL)
} else {

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  STEP 2a — SIMPLE AREA-WEIGHTED EXTRAPOLATION           ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

# ── Logging ──────────────────────────────────────────────────────────────────
log_file <- file.path("logs", paste0("simple_extrapolation_", Sys.Date(), ".log"))
if (!dir.exists("logs")) dir.create("logs", recursive = TRUE)

log_msg <- function(msg, level = "INFO") {
  entry <- sprintf("[%s] %s: %s", format(Sys.time(), "%H:%M:%S"), level, msg)
  cat(entry, "\n")
  cat(entry, "\n", file = log_file, append = TRUE)
}

log_msg("=== STEP 2a: SIMPLE EXTRAPOLATION STARTED ===")

# ── Packages ─────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# ── Output directories ────────────────────────────────────────────────────────
out_tables <- "outputs/Basic_analysis/Step2_Spatial_Basic/SimpleExtrapolation/tables"
out_maps   <- "outputs/Basic_analysis/Step2_Spatial_Basic/SimpleExtrapolation/maps"
out_plots  <- "outputs/Basic_analysis/Step2_Spatial_Basic/SimpleExtrapolation/plots"
for (d in c(out_tables, out_maps, out_plots)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================================
# 1. LOAD STRATUM CARBON SUMMARY (from S1_01)
# ============================================================================

stratum_csv <- "data_processed/carbon_by_stratum_summary.csv"
if (!file.exists(stratum_csv)) {
  stop("carbon_by_stratum_summary.csv not found. Run S1_01_data_prep.R first.")
}

stratum_summary <- read_csv(stratum_csv, show_col_types = FALSE)
log_msg(sprintf("Loaded stratum summary: %d strata", nrow(stratum_summary)))

# Expected columns: stratum, n_cores, mean_stock, sd_stock, min_stock, max_stock
required_cols <- c("stratum", "n_cores", "mean_stock", "sd_stock")
missing <- setdiff(required_cols, names(stratum_summary))
if (length(missing) > 0) {
  stop(sprintf("carbon_by_stratum_summary.csv missing columns: %s", paste(missing, collapse = ", ")))
}

# ============================================================================
# 2. LOAD AOI BOUNDARY
# ============================================================================

aoi_stratum_field <- if (exists("AOI_STRATUM_FIELD") && !is.null(AOI_STRATUM_FIELD) &&
                         nchar(trimws(AOI_STRATUM_FIELD)) > 0) AOI_STRATUM_FIELD else NULL

# Determine run mode: stratified (with stratum column) or unstratified (whole AOI)
run_stratified <- !is.null(aoi_stratum_field)

log_msg(sprintf("AOI file: %s", AOI_FILE))
if (run_stratified) {
  log_msg(sprintf("Mode: STRATIFIED — using stratum column '%s'", aoi_stratum_field))
} else {
  log_msg("Mode: UNSTRATIFIED — AOI_STRATUM_FIELD is NULL; treating entire AOI as one zone")
}

aoi <- tryCatch(
  sf::st_read(AOI_FILE, quiet = TRUE),
  error = function(e) stop(sprintf("Failed to read AOI file: %s", e$message))
)

if (run_stratified) {
  # Check stratum column exists in AOI
  if (!aoi_stratum_field %in% names(aoi)) {
    stop(sprintf(
      "Stratum column '%s' not found in AOI. Available columns: %s\n  Update AOI_STRATUM_FIELD in blue_carbon_config.R\n  Set AOI_STRATUM_FIELD <- NULL to run without strata.",
      aoi_stratum_field,
      paste(names(aoi), collapse = ", ")
    ))
  }
  # Standardize stratum column name
  aoi <- aoi %>% rename(stratum = !!aoi_stratum_field)
} else {
  # Unstratified: assign a single dummy stratum so downstream code is uniform
  aoi <- aoi %>% mutate(stratum = "AOI")
}

log_msg(sprintf("AOI loaded: %d features, CRS: %s", nrow(aoi), sf::st_crs(aoi)$input))

# ============================================================================
# 3. PROJECT TO EQUAL-AREA CRS AND CALCULATE STRATUM AREAS
# ============================================================================

processing_crs <- if (exists("PROCESSING_CRS")) PROCESSING_CRS else 3347

aoi_proj <- sf::st_transform(aoi, processing_crs)
aoi_proj$area_m2 <- as.numeric(sf::st_area(aoi_proj))

# Aggregate area by stratum (in case AOI has multiple polygons per stratum)
stratum_areas <- aoi_proj %>%
  sf::st_drop_geometry() %>%
  group_by(stratum) %>%
  summarise(area_m2 = sum(area_m2, na.rm = TRUE), .groups = "drop") %>%
  mutate(area_ha = area_m2 / 10000)

log_msg("Stratum areas (projected):")
for (i in seq_len(nrow(stratum_areas))) {
  log_msg(sprintf("  %s: %.1f ha (%.0f m²)",
                  stratum_areas$stratum[i],
                  stratum_areas$area_ha[i],
                  stratum_areas$area_m2[i]))
}

# ============================================================================
# 4. JOIN STRATUM CARBON DATA WITH AREAS
# ============================================================================

if (!run_stratified) {
  # Unstratified mode: pool all cores into one row, apply to total AOI area
  log_msg("Unstratified mode: pooling all cores across strata for a single AOI-wide estimate")

  # Drop any rows with NA stock or core count before pooling
  ss_clean <- stratum_summary %>%
    filter(!is.na(mean_stock), !is.na(n_cores), !is.na(sd_stock))

  # Extract vectors explicitly to avoid dplyr column-name shadowing in summarise()
  .w  <- ss_clean$n_cores
  .x  <- ss_clean$mean_stock
  .sd <- ss_clean$sd_stock

  pooled_summary <- data.frame(
    stratum    = "AOI",
    n_cores    = sum(.w),
    mean_stock = weighted.mean(.x, w = .w),
    sd_stock   = sqrt(sum((.w - 1) * .sd^2) / (sum(.w) - 1)),
    min_stock  = min(ss_clean$min_stock),
    max_stock  = max(ss_clean$max_stock)
  )
  joined <- pooled_summary %>%
    mutate(
      area_m2 = sum(stratum_areas$area_m2),
      area_ha = sum(stratum_areas$area_ha)
    )

} else {
  # Stratified mode: case-insensitive join on stratum name
  stratum_areas_lower  <- stratum_areas %>%
    mutate(stratum_key = tolower(trimws(stratum)))
  stratum_summary_lower <- stratum_summary %>%
    mutate(stratum_key = tolower(trimws(stratum)))

  joined <- stratum_summary_lower %>%
    left_join(stratum_areas_lower %>% select(stratum_key, area_m2, area_ha),
              by = "stratum_key") %>%
    select(-stratum_key)

  # Warn about unmatched strata
  unmatched_data <- joined %>% filter(is.na(area_m2))
  if (nrow(unmatched_data) > 0) {
    warning(sprintf(
      "These strata have carbon data but no matching AOI polygon:\n  %s\n  Check stratum names match between data and AOI.",
      paste(unmatched_data$stratum, collapse = ", ")
    ))
  }

  unmatched_aoi <- stratum_areas %>%
    filter(!tolower(trimws(stratum)) %in% tolower(trimws(stratum_summary$stratum)))
  if (nrow(unmatched_aoi) > 0) {
    warning(sprintf(
      "These AOI strata have no matching core data:\n  %s",
      paste(unmatched_aoi$stratum, collapse = ", ")
    ))
  }
}

# ============================================================================
# 5. CALCULATE CARBON STOCKS
# ============================================================================

# Confidence interval z-score
extrap_confidence <- if (exists("EXTRAP_CONFIDENCE")) EXTRAP_CONFIDENCE else 0.95
z_score <- qnorm(1 - (1 - extrap_confidence) / 2)

results_by_stratum <- joined %>%
  filter(!is.na(area_m2)) %>%
  mutate(
    # Aerial carbon stock (per unit area)
    aerial_carbon_kg_m2   = mean_stock,             # already in kg/m²
    aerial_carbon_MgC_ha  = mean_stock * 10,        # kg/m² → Mg C/ha

    # Total carbon stock for the stratum
    total_carbon_kg       = mean_stock * area_m2,
    total_carbon_Mg       = total_carbon_kg / 1000,
    total_carbon_MgC_ha   = aerial_carbon_MgC_ha,   # same per-area value

    # Standard error of the mean (spatial)
    se_mean_kg_m2         = sd_stock / sqrt(n_cores),

    # Confidence interval on total carbon stock
    ci_half_kg            = z_score * se_mean_kg_m2 * area_m2,
    total_carbon_kg_lower = total_carbon_kg - ci_half_kg,
    total_carbon_kg_upper = total_carbon_kg + ci_half_kg,

    # VM0033 conservative estimate (90% one-sided lower bound)
    conservative_kg       = total_carbon_kg - 1.645 * se_mean_kg_m2 * area_m2,

    # Metadata
    confidence_level      = extrap_confidence,
    method                = "Simple Extrapolation (mean × area)"
  )

# AOI-wide totals
aoi_total <- results_by_stratum %>%
  summarise(
    n_strata              = n(),
    total_area_m2         = sum(area_m2),
    total_area_ha         = sum(area_ha),
    total_carbon_kg       = sum(total_carbon_kg),
    total_carbon_Mg       = sum(total_carbon_Mg),
    se_total_kg           = sqrt(sum((se_mean_kg_m2 * area_m2)^2)),  # propagated SE
    ci_half_kg            = z_score * sqrt(sum((se_mean_kg_m2 * area_m2)^2)),
    total_carbon_kg_lower = total_carbon_kg - ci_half_kg,
    total_carbon_kg_upper = total_carbon_kg + ci_half_kg,
    conservative_kg       = total_carbon_kg - 1.645 * sqrt(sum((se_mean_kg_m2 * area_m2)^2)),
    weighted_aerial_kg_m2 = sum(aerial_carbon_kg_m2 * area_m2) / sum(area_m2),
    weighted_aerial_MgC_ha = weighted_aerial_kg_m2 * 10,
    confidence_level      = extrap_confidence,
    method                = "Simple Extrapolation (mean × area)"
  )

# ============================================================================
# 6. SAVE TABLES
# ============================================================================

write_csv(results_by_stratum, file.path(out_tables, "simple_extrapolation_by_stratum.csv"))
write_csv(aoi_total,          file.path(out_tables, "simple_extrapolation_aoi_total.csv"))

log_msg(sprintf("Saved by-stratum table: %d strata", nrow(results_by_stratum)))
log_msg(sprintf("AOI total carbon stock: %.1f Mg C (%.0f kg)", aoi_total$total_carbon_Mg, aoi_total$total_carbon_kg))

# ============================================================================
# 7. RASTERIZE AOI: CARBON STOCK MAP (GeoTIFF)
# ============================================================================

log_msg("Rasterizing AOI stratum polygons...")

# Add carbon stock values back to AOI polygons
aoi_proj_with_stock <- aoi_proj %>%
  mutate(stratum_key = tolower(trimws(stratum))) %>%
  left_join(
    results_by_stratum %>%
      mutate(stratum_key = tolower(trimws(stratum))) %>%
      select(stratum_key, aerial_carbon_kg_m2, aerial_carbon_MgC_ha, total_carbon_kg),
    by = "stratum_key"
  ) %>%
  select(-stratum_key)

tryCatch({
  # Create reference raster from AOI extent
  aoi_vect   <- terra::vect(aoi_proj_with_stock)
  ref_raster <- terra::rast(aoi_vect, resolution = 100)  # 100m cells

  # Rasterize aerial carbon stock
  r_aerial <- terra::rasterize(aoi_vect, ref_raster,
                                field = "aerial_carbon_kg_m2",
                                fun = "mean")
  names(r_aerial) <- "aerial_carbon_kg_m2"

  tif_path <- file.path(out_maps, "simple_extrapolation_carbon_stock_mean_kg_m2.tif")
  terra::writeRaster(r_aerial, tif_path, overwrite = TRUE,
                     gdal = c("COMPRESS=LZW", "TILED=YES"))
  log_msg(sprintf("Saved raster: %s", basename(tif_path)))
}, error = function(e) {
  log_msg(sprintf("Rasterization skipped: %s", e$message), "WARNING")
})

# ============================================================================
# 8. PLOTS
# ============================================================================

fig_dpi   <- if (exists("FIGURE_DPI"))   FIGURE_DPI   else 300
fig_w     <- if (exists("FIGURE_WIDTH")) FIGURE_WIDTH else 10
fig_h     <- if (exists("FIGURE_HEIGHT")) FIGURE_HEIGHT else 7
str_colors <- if (exists("STRATUM_COLORS")) STRATUM_COLORS else NULL

# ── Plot A: Choropleth map ────────────────────────────────────────────────────
tryCatch({
  aoi_wgs84 <- sf::st_transform(aoi_proj_with_stock, 4326)

  p_map <- ggplot(aoi_wgs84) +
    geom_sf(aes(fill = aerial_carbon_kg_m2), color = "white", linewidth = 0.3) +
    scale_fill_viridis_c(option = "plasma", na.value = "grey80",
                         name = "Carbon Stock\n(kg C/m²)") +
    labs(
      title    = "Step 2a — Simple Extrapolation: Aerial Carbon Stock",
      subtitle = sprintf("Method: mean carbon stock per stratum × stratum area | n = %d strata",
                         nrow(results_by_stratum)),
      x = "Longitude", y = "Latitude",
      caption  = sprintf("Total AOI stock: %.1f Mg C (95%% CI: %.1f – %.1f Mg C)",
                         aoi_total$total_carbon_Mg,
                         aoi_total$total_carbon_kg_lower / 1000,
                         aoi_total$total_carbon_kg_upper / 1000)
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(out_plots, "simple_extrapolation_carbon_stock_map.png"),
         p_map, width = fig_w, height = fig_h, dpi = fig_dpi)
  log_msg("Saved: simple_extrapolation_carbon_stock_map.png")
}, error = function(e) log_msg(sprintf("Map plot skipped: %s", e$message), "WARNING"))

# ── Plot B: Stratum bar chart with CI ────────────────────────────────────────
tryCatch({
  p_bar <- ggplot(results_by_stratum,
                  aes(x = reorder(stratum, -aerial_carbon_kg_m2),
                      y = aerial_carbon_kg_m2)) +
    geom_col(aes(fill = stratum), alpha = 0.85, width = 0.7) +
    geom_errorbar(
      aes(ymin = aerial_carbon_kg_m2 - z_score * se_mean_kg_m2,
          ymax = aerial_carbon_kg_m2 + z_score * se_mean_kg_m2),
      width = 0.25, linewidth = 0.8, color = "black"
    ) +
    geom_text(aes(label = sprintf("n=%d\n%.0f ha", n_cores, area_ha)),
              vjust = -0.4, size = 3.2, color = "gray30") +
    labs(
      title    = "Step 2a — Mean Aerial Carbon Stock by Stratum",
      subtitle = sprintf("Error bars = %d%% CI | Stratum areas from AOI",
                         round(extrap_confidence * 100)),
      x        = "Stratum",
      y        = "Aerial Carbon Stock (kg C/m²)",
      fill     = "Stratum"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x    = element_text(angle = 35, hjust = 1),
      legend.position = "none",
      plot.title     = element_text(face = "bold")
    )

  if (!is.null(str_colors)) {
    p_bar <- p_bar + scale_fill_manual(values = str_colors)
  }

  ggsave(file.path(out_plots, "simple_extrapolation_stratum_summary_chart.png"),
         p_bar, width = fig_w, height = fig_h, dpi = fig_dpi)
  log_msg("Saved: simple_extrapolation_stratum_summary_chart.png")
}, error = function(e) log_msg(sprintf("Bar chart skipped: %s", e$message), "WARNING"))

# ============================================================================
# 9. CONSOLE SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("STEP 2a COMPLETE — SIMPLE EXTRAPOLATION\n")
cat("========================================\n\n")

cat(sprintf("Strata processed: %d\n", nrow(results_by_stratum)))
cat(sprintf("AOI total area:   %.1f ha\n\n", aoi_total$total_area_ha))

cat("Carbon stocks by stratum:\n")
for (i in seq_len(nrow(results_by_stratum))) {
  r <- results_by_stratum[i, ]
  cat(sprintf("  %-20s  %.2f kg/m²  |  %.1f Mg C/ha  |  %.1f Mg C total\n",
              r$stratum, r$aerial_carbon_kg_m2, r$aerial_carbon_MgC_ha, r$total_carbon_Mg))
}

cat(sprintf("\nAOI TOTAL: %.1f Mg C\n", aoi_total$total_carbon_Mg))
cat(sprintf("  %d%% CI:   [%.1f – %.1f] Mg C\n",
            round(extrap_confidence * 100),
            aoi_total$total_carbon_kg_lower / 1000,
            aoi_total$total_carbon_kg_upper / 1000))
cat(sprintf("  VM0033 conservative (90%% lower): %.1f Mg C\n\n",
            aoi_total$conservative_kg / 1000))

cat("Outputs:\n")
cat(sprintf("  Tables: %s\n", out_tables))
cat(sprintf("  Maps:   %s\n", out_maps))
cat(sprintf("  Plots:  %s\n\n", out_plots))

log_msg("=== STEP 2a COMPLETE ===")

} # end AOI gate
