# ── BlueCarbon pre-processing pipeline ───────────────────────────────────────
# Wraps the BlueCarbon package (EcologyR/BlueCarbon) to:
#   1. Estimate and correct for percussion core compaction
#   2. Estimate OC stocks with the BlueCarbon method (linear extrapolation)
#   3. Validate short-core extrapolation quality
#   4. Produce decompacted, gap-filled cores compatible with harmonize_depths()
#
# Requires: install.packages("BlueCarbon")
#
# Unit notes:
#   BlueCarbon expects OC as % (g OC / 100 g soil)  → soc_pct = soc_g_kg / 10
#   BlueCarbon stock output is g C / cm²             → × 10 = kg C / m²  (VM0033)

# ── 1. Compaction correction ─────────────────────────────────────────────────

run_compaction_correction <- function(samples_path, compaction_path) {
  samples   <- readr::read_csv(samples_path,    show_col_types = FALSE)
  comp_meas <- readr::read_csv(compaction_path, show_col_types = FALSE)

  # Estimate % compaction per core from field coring measurements
  comp_pct <- BlueCarbon::estimate_compaction(
    comp_meas,
    core              = "core_id",
    sampler_length    = "sampler_length",
    internal_distance = "internal_distance",
    external_distance = "external_distance"
  )

  message("[bluecarbon] Compaction estimates (%):")
  print(comp_pct[, c("core_id", "compaction")])

  # Join compaction onto samples; cores without a measurement get 0% (no correction)
  samples_prep <- samples |>
    dplyr::left_join(
      dplyr::select(comp_pct, core_id, compaction),
      by = "core_id"
    ) |>
    dplyr::mutate(
      compaction = dplyr::coalesce(compaction, 0),
      soc_pct    = soc_g_kg / 10    # g/kg → % for BlueCarbon
    )

  uncorrected <- samples_prep$core_id[samples_prep$compaction == 0] |> unique()
  if (length(uncorrected) > 0)
    message("[bluecarbon] No compaction measurement for: ",
            paste(uncorrected, collapse = ", "), " — using 0% (no correction).")

  # Correct depths and bulk density for compaction
  BlueCarbon::decompact(
    samples_prep,
    core       = "core_id",
    compaction = "compaction",
    mind       = "depth_top_cm",
    maxd       = "depth_bottom_cm",
    dbd        = "bulk_density_g_cm3"
  )
}

# ── 2. BlueCarbon OC stocks (linear extrapolation to 100 cm) ─────────────────

estimate_bc_stocks <- function(locations_path, decompacted_samples, depth = 100) {
  # estimate_oc_stock() calls estimate_h() internally for gap filling
  raw_stocks <- BlueCarbon::estimate_oc_stock(
    decompacted_samples,
    depth = depth,
    core  = "core_id",
    mind  = "mind_corrected",
    maxd  = "maxd_corrected",
    dbd   = "dbd_corrected",
    oc    = "soc_pct"
  )

  # BlueCarbon always names the ID column "core" — rename to project convention
  locs <- readr::read_csv(locations_path, show_col_types = FALSE)

  raw_stocks |>
    dplyr::rename(core_id = core) |>
    dplyr::left_join(
      dplyr::select(locs, core_id, stratum, latitude, longitude),
      by = "core_id"
    ) |>
    dplyr::mutate(
      stock_kg_m2    = stock    * 10,   # g/cm² → kg/m²
      stockwc_kg_m2  = stockwc  * 10,
      stock_se_kg_m2 = stock_se * 10
    )
}

# ── 3. Extrapolation validation ───────────────────────────────────────────────

run_extrapolation_test <- function(decompacted_samples, depth = 100) {
  # Only cores that reach `depth` can be used as reference — short-core dataset
  # will produce a message if no cores qualify.
  BlueCarbon::test_extrapolation(
    decompacted_samples,
    depth = depth,
    core  = "core_id",
    mind  = "mind_corrected",
    maxd  = "maxd_corrected",
    dbd   = "dbd_corrected",
    oc    = "soc_pct"
  )
}

# ── 4. Prepare decompacted cores for harmonize_depths() ──────────────────────
#
# Converts the BlueCarbon output back to the format expected by load_raw_data()
# output — i.e. the `cores_raw` structure used by harmonize_depths().
# This allows the decompacted/gap-filled cores to feed directly into the
# existing VM0033 spatial pipeline.
#
# Key changes from raw data:
#   depth_top/bottom  → corrected by decompaction factor
#   layer_thickness   → gap-filled by estimate_h() (midpoint split for non-contiguous samples)
#   bulk_density      → corrected for compaction (density increases when core expands)

prepare_bc_cores_for_harmonization <- function(locations_path, decompacted_samples) {
  # estimate_h() fills depth gaps between non-contiguous samples (midpoint split).
  # Returns emin, emax, h columns — h is the effective sample thickness.
  with_h <- BlueCarbon::estimate_h(
    decompacted_samples,
    core = "core_id",
    mind = "mind_corrected",
    maxd = "maxd_corrected"
  )

  locs <- readr::read_csv(locations_path, show_col_types = FALSE)

  with_h |>
    dplyr::mutate(
      depth_cm           = (mind_corrected + maxd_corrected) / 2,
      layer_thickness_cm = h,
      carbon_stock_kg_m2 = soc_g_kg * dbd_corrected * h / 100,
      bd_estimated       = FALSE
    ) |>
    dplyr::select(-depth_top_cm, -depth_bottom_cm, -bulk_density_g_cm3) |>
    dplyr::rename(
      depth_top_cm       = mind_corrected,
      depth_bottom_cm    = maxd_corrected,
      bulk_density_g_cm3 = dbd_corrected
    ) |>
    dplyr::select(
      core_id, depth_top_cm, depth_bottom_cm, depth_cm, layer_thickness_cm,
      soc_g_kg, bulk_density_g_cm3, carbon_stock_kg_m2, bd_estimated,
      compaction, soc_pct, h, emin, emax
    ) |>
    dplyr::left_join(
      dplyr::select(locs, core_id, stratum, latitude, longitude),
      by = "core_id"
    )
}

# Writes the harmonization-ready cores to a CSV file (tar format = "file")
write_bc_cores_csv <- function(bc_cores, output_dir = "outputs/bluecarbon") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(output_dir, "core_samples_bc_processed.csv")
  readr::write_csv(bc_cores, out_path)
  message("[bluecarbon] Processed cores written to: ", out_path)
  message("[bluecarbon] To use in main pipeline, set this file as your core_samples input.")
  out_path
}

# ── 5. Comparison with VM0033 main pipeline ───────────────────────────────────

compare_bc_vs_vm0033 <- function(bc_stocks, vm_stratum_summary) {
  if (is.null(vm_stratum_summary)) {
    message("[bluecarbon] VM0033 stratum summary not available — run main pipeline first.")
    return(NULL)
  }

  bc_by_stratum <- bc_stocks |>
    dplyr::group_by(stratum) |>
    dplyr::summarise(
      n_cores           = dplyr::n(),
      mean_bc_kg_m2     = mean(stock_kg_m2, na.rm = TRUE),
      sd_bc_kg_m2       = sd(stock_kg_m2,   na.rm = TRUE),
      .groups = "drop"
    )

  vm_totals <- vm_stratum_summary |>
    dplyr::group_by(stratum) |>
    dplyr::summarise(
      mean_vm0033_kg_m2 = sum(mean_stock, na.rm = TRUE),
      .groups = "drop"
    )

  dplyr::left_join(bc_by_stratum, vm_totals, by = "stratum") |>
    dplyr::mutate(
      diff_kg_m2 = mean_bc_kg_m2 - mean_vm0033_kg_m2,
      diff_pct   = round(100 * diff_kg_m2 / mean_vm0033_kg_m2, 1)
    )
}

plot_bc_comparison <- function(comparison_df) {
  if (is.null(comparison_df)) return(ggplot2::ggplot() + ggplot2::labs(
    title = "Comparison not available — run main pipeline first"))

  long <- comparison_df |>
    tidyr::pivot_longer(
      cols      = c(mean_bc_kg_m2, mean_vm0033_kg_m2),
      names_to  = "method",
      values_to = "stock_kg_m2"
    ) |>
    dplyr::mutate(method = dplyr::recode(method,
      mean_bc_kg_m2     = "BlueCarbon (decompacted)",
      mean_vm0033_kg_m2 = "VM0033 (raw depths)"
    ))

  ggplot2::ggplot(long,
    ggplot2::aes(x = stratum, y = stock_kg_m2, fill = method)) +
    ggplot2::geom_col(position = "dodge", width = 0.6) +
    ggplot2::geom_text(
      data = comparison_df,
      ggplot2::aes(
        x     = stratum,
        y     = pmax(mean_bc_kg_m2, mean_vm0033_kg_m2) + 0.3,
        label = paste0(ifelse(diff_pct > 0, "+", ""), diff_pct, "%")
      ),
      inherit.aes = FALSE, size = 3.5, colour = "grey30"
    ) +
    ggplot2::scale_fill_manual(values = c(
      "BlueCarbon (decompacted)" = "#2c7a4b",
      "VM0033 (raw depths)"      = "#6baed6"
    )) +
    ggplot2::labs(
      title    = "Carbon stock: BlueCarbon (decompacted) vs VM0033 (raw depths)",
      subtitle = "% label = effect of compaction correction on each stratum",
      x = "Stratum", y = "Mean total stock (kg C/m²)", fill = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}
