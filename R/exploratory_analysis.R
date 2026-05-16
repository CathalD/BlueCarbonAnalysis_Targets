# R/exploratory_analysis.R
# ============================================================================
# PURPOSE: Generate Step 1 EDA plots as a named list of ggplot objects.
#
# INPUT:
#   cores_clean — QA-flagged data frame from run_qc()
#   cfg         — named list from load_config()
#
# OUTPUT — named list of ggplot objects:
#   eda_plots$depth_profiles    — SOC vs depth, one line per core, by stratum
#   eda_plots$soc_by_stratum    — SOC distribution (violin + jitter)
#   eda_plots$bd_distribution   — BD density, measured vs estimated flagged
#   eda_plots$carbon_stocks     — total stock per core, by stratum
#   eda_plots$spatial_map       — core locations coloured by stratum
#   eda_plots$qa_summary        — bar chart of QA failures by type
# ============================================================================
run_eda <- function(cores_clean, cfg) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(tidyr)
  })
  message("[eda] Building EDA plots...")
  cores_qa <- cores_clean |> filter(qa_pass)

  stratum_colours <- cfg$STRATUM_COLORS
  colour_scale <- if (!is.null(stratum_colours)) {
    scale_color_manual(values = stratum_colours)
  } else {
    scale_color_viridis_d(option = "D")
  }
  fill_scale <- if (!is.null(stratum_colours)) {
    scale_fill_manual(values = stratum_colours)
  } else {
    scale_fill_viridis_d(option = "D")
  }

  # ── 1. SOC depth profiles ─────────────────────────────────────────────────
  p_depth_profiles <- ggplot(
    cores_qa |> arrange(core_id, depth_cm),
    aes(x = soc_g_kg, y = depth_cm, group = core_id, colour = stratum)
  ) +
    geom_path(alpha = 0.5, linewidth = 0.4) +
    geom_point(alpha = 0.3, size = 0.8) +
    scale_y_reverse(name = "Depth (cm)") +
    scale_x_continuous(name = "SOC (g/kg)") +
    facet_wrap(~stratum, scales = "free_x") +
    colour_scale +
    theme_bw(base_size = 11) +
    theme(legend.position = "none") +
    labs(title = "SOC depth profiles by stratum",
         caption = "Each line = one core. QA-passed samples only.")

  # ── 2. SOC distribution by stratum ───────────────────────────────────────
  p_soc_by_stratum <- ggplot(
    cores_qa,
    aes(x = stratum, y = soc_g_kg, fill = stratum)
  ) +
    geom_violin(alpha = 0.5, draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    fill_scale +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(title = "SOC distribution by stratum", x = NULL, y = "SOC (g/kg)")

  # ── 3. BD distribution, measured vs estimated ─────────────────────────────
  p_bd_distribution <- ggplot(
    cores_qa |> filter(!is.na(bulk_density_g_cm3)),
    aes(x = bulk_density_g_cm3, fill = stratum, linetype = bd_estimated)
  ) +
    geom_density(alpha = 0.4) +
    fill_scale +
    scale_linetype_manual(
      values = c("FALSE" = "solid", "TRUE" = "dashed"),
      labels = c("FALSE" = "Measured", "TRUE" = "Default (estimated)"),
      name = "BD source"
    ) +
    theme_bw(base_size = 11) +
    labs(title = "Bulk density distribution",
         x = "Bulk density (g/cm³)", y = "Density",
         caption = "Dashed = literature defaults applied (bd_estimated = TRUE).")

  # ── 4. Total carbon stock per core ───────────────────────────────────────
  core_totals <- cores_qa |>
    group_by(core_id, stratum) |>
    summarise(total_stock = sum(carbon_stock_kg_m2, na.rm = TRUE), .groups = "drop")

  p_carbon_stocks <- ggplot(core_totals, aes(x = stratum, y = total_stock, fill = stratum)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21) +
    fill_scale +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(title = "Total carbon stock by stratum (raw, pre-harmonization)",
         x = NULL, y = expression("Carbon stock (kg C m"^{-2}*")"))

  # ── 5. Spatial map ────────────────────────────────────────────────────────
  p_spatial_map <- ggplot(
    cores_qa |> distinct(core_id, longitude, latitude, stratum),
    aes(x = longitude, y = latitude, colour = stratum, shape = stratum)
  ) +
    geom_point(size = 2.5, alpha = 0.8) +
    colour_scale +
    coord_sf() +
    theme_bw(base_size = 11) +
    labs(title = "Core locations", x = "Longitude", y = "Latitude",
         colour = "Stratum", shape = "Stratum")

  # ── 6. QA failures by type ────────────────────────────────────────────────
  qa_summary_df <- cores_clean |>
    summarise(
      "Depth invalid"      = sum(!qa_depth_valid,   na.rm = TRUE),
      "SOC out of range"   = sum(!qa_soc_valid,     na.rm = TRUE),
      "BD out of range"    = sum(!qa_bd_valid,       na.rm = TRUE),
      "Missing location"   = sum(!qa_has_location,  na.rm = TRUE)
    ) |>
    pivot_longer(everything(), names_to = "issue", values_to = "n_samples")

  p_qa_summary <- ggplot(qa_summary_df, aes(x = reorder(issue, -n_samples), y = n_samples)) +
    geom_col(fill = "#d95f02", alpha = 0.8) +
    theme_bw(base_size = 11) +
    labs(title = "QA/QC: samples flagged by issue type",
         x = NULL, y = "Number of samples flagged",
         caption = "A sample can be flagged for multiple issues.")

  message("[eda] Done. Returning named list of 6 plots.")
  list(
    depth_profiles  = p_depth_profiles,
    soc_by_stratum  = p_soc_by_stratum,
    bd_distribution = p_bd_distribution,
    carbon_stocks   = p_carbon_stocks,
    spatial_map     = p_spatial_map,
    qa_summary      = p_qa_summary
  )
}
