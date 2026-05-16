# R/exploratory_analysis.R
# PURPOSE: Generate Step 1 EDA plots as a named list of ggplot objects.
#
# OUTPUT — named list of 3 ggplot objects:
#   eda_plots$depth_profiles — SOC vs depth, one line per core, by stratum
#   eda_plots$spatial_map   — core locations coloured by stratum
#   eda_plots$carbon_stocks — total raw stock per core, by stratum
run_eda <- function(cores_raw, cfg) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
  })
  message("[eda] Building EDA plots...")

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
    cores_raw |> arrange(core_id, depth_cm),
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
         caption = "Each line = one core.")

  # ── 2. Spatial map ────────────────────────────────────────────────────────
  p_spatial_map <- ggplot(
    cores_raw |> distinct(core_id, longitude, latitude, stratum),
    aes(x = longitude, y = latitude, colour = stratum, shape = stratum)
  ) +
    geom_point(size = 2.5, alpha = 0.8) +
    colour_scale +
    coord_sf() +
    theme_bw(base_size = 11) +
    labs(title = "Core locations", x = "Longitude", y = "Latitude",
         colour = "Stratum", shape = "Stratum")

  # ── 3. Total carbon stock per core ───────────────────────────────────────
  core_totals <- cores_raw |>
    filter(!is.na(carbon_stock_kg_m2)) |>
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

  message("[eda] Done. Returning 3 plots.")
  list(
    depth_profiles = p_depth_profiles,
    spatial_map    = p_spatial_map,
    carbon_stocks  = p_carbon_stocks
  )
}
