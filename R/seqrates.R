# ── Sequestration rate estimation ────────────────────────────────────────────
# Wraps BlueCarbon::estimate_seq_rate() for cores that have chronology data.
#
# Data requirements:
#   core_chronology.csv  — age-depth anchor points from radiometric dating
#     columns: core_id, depth_cm, age_ybp, dating_method (optional)
#     source:  ^210Pb (surface ~0–30 cm) or ^14C (deeper horizons)
#     at least 2 anchor points per core (surface + one date at depth)
#
# Ages are linearly interpolated between anchor points to assign an age
# to each sample midpoint. This is not a full age-depth model — for a
# proper Bayesian age model use Bacon/Bchron externally and import the
# modelled ages via this same CSV format.
#
# Unit notes:
#   BlueCarbon seq_rate output is g C / cm² / yr
#   × 10000 cm²/m² = g C / m² / yr  (typical reporting unit for blue carbon)
#   × 10             = kg C / m² / yr

# ── 1. Load chronology and assign ages to samples ────────────────────────────

load_and_assign_ages <- function(samples_path, compaction_path, chronology_path) {
  samples  <- readr::read_csv(samples_path,     show_col_types = FALSE)
  chrono   <- readr::read_csv(chronology_path,  show_col_types = FALSE)
  comp_meas <- readr::read_csv(compaction_path, show_col_types = FALSE)

  cores_with_chrono <- unique(chrono$core_id)
  message("[seqrates] Cores with chronology data: ",
          paste(cores_with_chrono, collapse = ", "))

  # Estimate compaction for chronology cores
  comp_pct <- BlueCarbon::estimate_compaction(
    comp_meas,
    core              = "core_id",
    sampler_length    = "sampler_length",
    internal_distance = "internal_distance",
    external_distance = "external_distance"
  )

  samples_sub <- samples |>
    dplyr::filter(core_id %in% cores_with_chrono) |>
    dplyr::left_join(dplyr::select(comp_pct, core_id, compaction), by = "core_id") |>
    dplyr::mutate(
      compaction = dplyr::coalesce(compaction, 0),
      soc_pct    = soc_g_kg / 10
    )

  decompacted <- BlueCarbon::decompact(
    samples_sub,
    core       = "core_id",
    compaction = "compaction",
    mind       = "depth_top_cm",
    maxd       = "depth_bottom_cm",
    dbd        = "bulk_density_g_cm3"
  )

  # Assign an age to each sample via linear interpolation of chronology anchors.
  # Each core's anchor points (depth_cm, age_ybp) define the age-depth relationship.
  # Samples outside the anchor range are extrapolated using the nearest slope.
  assign_ages_by_core(decompacted, chrono)
}

assign_ages_by_core <- function(decompacted, chrono) {
  cores <- split(decompacted, decompacted$core_id)

  result <- lapply(names(cores), function(cid) {
    df <- cores[[cid]]
    cc <- chrono[chrono$core_id == cid, ]

    df$depth_mid <- (df$mind_corrected + df$maxd_corrected) / 2

    if (nrow(cc) < 2) {
      message("[seqrates] Core ", cid, ": fewer than 2 chronology points — age set to NA.")
      df$age <- NA_real_
    } else {
      df$age <- stats::approx(
        x    = cc$depth_cm,
        y    = cc$age_ybp,
        xout = df$depth_mid,
        rule = 2    # extrapolate beyond range using the nearest-end slope
      )$y
    }
    df
  })

  do.call(rbind, result)
}

# ── 2. Estimate sequestration rates ──────────────────────────────────────────

estimate_sequestration_rates <- function(cores_with_ages, timeframe = 100) {
  eligible <- cores_with_ages |>
    dplyr::group_by(core_id) |>
    dplyr::filter(!is.na(age), max(age, na.rm = TRUE) >= timeframe) |>
    dplyr::ungroup()

  excluded <- setdiff(unique(cores_with_ages$core_id), unique(eligible$core_id))
  if (length(excluded) > 0)
    message("[seqrates] Cores excluded (max age < ", timeframe, " yr): ",
            paste(excluded, collapse = ", "))

  if (nrow(eligible) == 0) {
    message("[seqrates] No cores span ", timeframe, " years — cannot estimate rates.")
    return(NULL)
  }

  rates <- BlueCarbon::estimate_seq_rate(
    eligible,
    timeframe = timeframe,
    core      = "core_id",
    mind      = "mind_corrected",
    maxd      = "maxd_corrected",
    dbd       = "dbd_corrected",
    oc        = "soc_pct",
    age       = "age"
  )

  # BlueCarbon always names the ID column "core" in output — rename
  rates |>
    dplyr::rename(core_id = core) |>
    dplyr::mutate(
      seq_rate_g_m2_yr    = seq_rate    * 10000,  # g/cm²/yr → g C/m²/yr
      seq_rate_wc_g_m2_yr = seq_rate_wc * 10000,
      timeframe_yr        = timeframe
    )
}

# ── 3. Visualize sequestration rates ─────────────────────────────────────────

plot_seq_rates <- function(seq_rates, locations_path) {
  if (is.null(seq_rates) || nrow(seq_rates) == 0)
    return(ggplot2::ggplot() + ggplot2::labs(
      title = "No sequestration rates available",
      subtitle = "Ensure chronology cores have max age ≥ timeframe"))

  locs <- readr::read_csv(locations_path, show_col_types = FALSE)

  rates <- seq_rates |>
    dplyr::left_join(
      dplyr::select(locs, core_id, stratum, latitude, longitude),
      by = "core_id"
    )

  tf <- rates$timeframe_yr[1]

  ggplot2::ggplot(rates,
    ggplot2::aes(x = core_id, y = seq_rate_g_m2_yr, fill = stratum)) +
    ggplot2::geom_col(width = 0.6, alpha = 0.85) +
    ggplot2::geom_point(
      ggplot2::aes(y = seq_rate_wc_g_m2_yr),
      shape = 21, fill = "white", colour = "grey30", size = 3
    ) +
    ggplot2::scale_fill_manual(
      values = c(IM = "#2c7a4b", NM = "#6baed6", MF = "#fd8d3c"),
      na.value = "#aaaaaa"
    ) +
    ggplot2::labs(
      title    = paste0("Carbon sequestration rates — ", tf, "-year average"),
      subtitle = "Bar = rate to target depth  ·  Circle = whole-core rate",
      x        = NULL,
      y        = "Sequestration rate (g C/m²/yr)",
      fill     = "Stratum"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "right")
}

# ── 4. Summary table ─────────────────────────────────────────────────────────

summarise_seq_rates <- function(seq_rates, locations_path) {
  if (is.null(seq_rates) || nrow(seq_rates) == 0) return(NULL)

  locs <- readr::read_csv(locations_path, show_col_types = FALSE)

  seq_rates |>
    dplyr::left_join(dplyr::select(locs, core_id, stratum), by = "core_id") |>
    dplyr::group_by(stratum) |>
    dplyr::summarise(
      n_cores                  = dplyr::n(),
      mean_seq_rate_g_m2_yr    = round(mean(seq_rate_g_m2_yr,    na.rm = TRUE), 1),
      sd_seq_rate_g_m2_yr      = round(sd(seq_rate_g_m2_yr,      na.rm = TRUE), 1),
      mean_maxage_yr           = round(mean(maxage,               na.rm = TRUE), 0),
      .groups = "drop"
    )
}
