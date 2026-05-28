# ── Lead-210 age model comparison ─────────────────────────────────────────────
# Implements CRS (Constant Rate of Supply) and CIC (Constant Initial
# Concentration) age-depth models using the pb210 R package
# (paleolimbot/pb210), and compares them against the linear interpolation
# method already in seqrates.R.
#
# Install:
#   remotes::install_github("paleolimbot/pb210")
#
# Required data: core_pb210.csv
#   core_id                    — must match core_locations.csv
#   depth_top_cm / depth_bottom_cm — slice boundaries (cm)
#   pb210_total_Bq_kg          — total 210Pb activity (Bq/kg)
#   pb210_total_sd_Bq_kg       — 1-sigma measurement uncertainty (Bq/kg)
#   pb210_supported_Bq_kg      — background 210Pb (Ra-226 or deep-section mean)
#   pb210_supported_sd_Bq_kg   — uncertainty on supported activity
#   dry_mass_g                 — oven-dry mass of the slice (g)
#   core_area_cm2              — cross-sectional area of the corer tube (cm²)
#                                (π × r² — e.g. 19.63 for a 5 cm diameter corer)
#
# Model descriptions
# ──────────────────
# CRS (Constant Rate of Supply / Appleby & Oldfield 1978):
#   Assumes a constant flux of unsupported 210Pb to the sediment surface,
#   while allowing the sedimentation rate to vary over time.  Best for sites
#   with variable accumulation.  Ages derived from the cumulative inventory
#   of excess 210Pb remaining below each horizon.
#
# CIC (Constant Initial Concentration / Krishnaswami et al. 1971):
#   Assumes a constant initial 210Pb activity at the time of deposition,
#   implying a constant sedimentation rate.  Simpler model; less reliable
#   where accumulation has varied but provides a useful comparison.
#
# Linear interpolation (current default in seqrates.R):
#   Uses pre-assigned ages from core_chronology.csv, interpolated linearly
#   between radiometric anchor points.  No uncertainty propagation.

# ── 1. Load raw 210Pb measurements ───────────────────────────────────────────

load_pb210_data <- function(pb210_path) {
  df <- readr::read_csv(pb210_path, show_col_types = FALSE)

  required <- c(
    "core_id", "depth_top_cm", "depth_bottom_cm",
    "pb210_total_Bq_kg", "pb210_total_sd_Bq_kg",
    "pb210_supported_Bq_kg", "pb210_supported_sd_Bq_kg",
    "dry_mass_g", "core_area_cm2"
  )
  missing <- setdiff(required, names(df))
  if (length(missing) > 0)
    stop("[pb210] Missing required columns in core_pb210.csv: ",
         paste(missing, collapse = ", "))

  df |>
    dplyr::mutate(
      depth_mid_cm    = (depth_top_cm + depth_bottom_cm) / 2,
      pb210_excess    = pb210_total_Bq_kg    - pb210_supported_Bq_kg,
      pb210_excess_sd = sqrt(pb210_total_sd_Bq_kg^2 + pb210_supported_sd_Bq_kg^2)
    ) |>
    dplyr::filter(pb210_excess > 0)  # drop samples at or below detection limit
}

# ── 2. Fit CRS and CIC age models per core ───────────────────────────────────

fit_pb210_age_models <- function(pb210_data) {
  if (!requireNamespace("pb210", quietly = TRUE))
    stop("[pb210] Package 'pb210' is not installed.\n",
         "  Install with: remotes::install_github('paleolimbot/pb210')")

  cores <- split(pb210_data, pb210_data$core_id)

  lapply(cores, function(df) {
    df       <- dplyr::arrange(df, depth_top_cm)
    area_m2  <- df$core_area_cm2[1] / 10000           # cm² → m²
    cum_mass <- cumsum(df$dry_mass_g / 1000) / area_m2 # g → kg, per m²

    crs <- tryCatch(
      pb210::pb210_crs(
        cumulative_dry_mass   = cum_mass,
        excess_pb210_activity = df$pb210_excess,
        excess_pb210_sd       = df$pb210_excess_sd
      ),
      error = function(e) {
        message("[pb210] CRS model failed for ", df$core_id[1], ": ", e$message)
        NULL
      }
    )

    cic <- tryCatch(
      pb210::pb210_cic(
        cumulative_dry_mass   = cum_mass,
        excess_pb210_activity = df$pb210_excess,
        excess_pb210_sd       = df$pb210_excess_sd
      ),
      error = function(e) {
        message("[pb210] CIC model failed for ", df$core_id[1], ": ", e$message)
        NULL
      }
    )

    list(crs = crs, cic = cic, cum_mass = cum_mass, depth_mid_cm = df$depth_mid_cm)
  })
}

# Helper: extract numeric value and sd from a possibly-errors-package vector
.extract_val <- function(x) {
  if (is.null(x)) return(rep(NA_real_, 0))
  if (requireNamespace("errors", quietly = TRUE) && inherits(x, "errors"))
    return(errors::drop_errors(x))
  suppressWarnings(as.numeric(x))
}
.extract_sd <- function(x) {
  if (is.null(x)) return(rep(NA_real_, 0))
  if (requireNamespace("errors", quietly = TRUE) && inherits(x, "errors"))
    return(errors::errors(x))
  rep(NA_real_, length(x))
}

# ── 3. Extract age predictions from fitted models ─────────────────────────────

extract_pb210_ages <- function(pb210_models) {
  results <- lapply(names(pb210_models), function(cid) {
    m <- pb210_models[[cid]]
    if (is.null(m)) return(NULL)

    nd <- data.frame(cumulative_dry_mass = m$cum_mass)

    # CRS
    crs_row <- tryCatch({
      p <- predict(m$crs, newdata = nd)
      data.frame(
        age_crs    = .extract_val(p$age),
        age_crs_sd = if ("age_sd" %in% names(p)) p$age_sd else .extract_sd(p$age)
      )
    }, error = function(e) {
      message("[pb210] CRS predict failed for ", cid, ": ", e$message)
      data.frame(age_crs = NA_real_, age_crs_sd = NA_real_)
    })

    # CIC
    cic_row <- tryCatch({
      p <- predict(m$cic, newdata = nd)
      data.frame(
        age_cic    = .extract_val(p$age),
        age_cic_sd = if ("age_sd" %in% names(p)) p$age_sd else .extract_sd(p$age)
      )
    }, error = function(e) {
      message("[pb210] CIC predict failed for ", cid, ": ", e$message)
      data.frame(age_cic = NA_real_, age_cic_sd = NA_real_)
    })

    n <- length(m$depth_mid_cm)
    data.frame(
      core_id      = cid,
      depth_mid_cm = m$depth_mid_cm,
      cum_mass     = m$cum_mass,
      age_crs      = if (nrow(crs_row) == n) crs_row$age_crs    else rep(NA_real_, n),
      age_crs_sd   = if (nrow(crs_row) == n) crs_row$age_crs_sd else rep(NA_real_, n),
      age_cic      = if (nrow(cic_row) == n) cic_row$age_cic    else rep(NA_real_, n),
      age_cic_sd   = if (nrow(cic_row) == n) cic_row$age_cic_sd else rep(NA_real_, n)
    )
  })

  dplyr::bind_rows(results[!sapply(results, is.null)])
}

# ── 4. Assign pb210 model ages to decompacted SOC samples ────────────────────
# Takes sr_cores_with_ages (which has depth_mid and a linear-interp age column)
# and replaces the age column with CRS or CIC ages interpolated from pb210_ages.

assign_pb210_ages_to_cores <- function(cores_with_ages, pb210_ages, method = "crs") {
  age_col <- paste0("age_", method)
  if (!age_col %in% names(pb210_ages))
    stop("[pb210] Column '", age_col, "' not found in pb210_ages.")

  cores   <- split(cores_with_ages, cores_with_ages$core_id)
  results <- lapply(names(cores), function(cid) {
    df  <- cores[[cid]]
    pb  <- pb210_ages[pb210_ages$core_id == cid, ]

    if (nrow(pb) == 0 || all(is.na(pb[[age_col]]))) {
      message("[pb210] No ", toupper(method), " ages for core ", cid,
              " — keeping linear interpolation ages.")
      return(df)
    }

    df$age <- stats::approx(
      x    = pb$depth_mid_cm,
      y    = pb[[age_col]],
      xout = df$depth_mid,
      rule = 2
    )$y
    df
  })

  dplyr::bind_rows(results)
}

# ── 5. Compare sequestration rates across all three age models ────────────────

compare_seq_rate_methods <- function(linear_rates, crs_rates, cic_rates) {
  bind_method <- function(rates, label) {
    if (is.null(rates) || nrow(rates) == 0) return(NULL)
    rates |>
      dplyr::select(core_id, seq_rate_g_m2_yr, seq_rate_wc_g_m2_yr) |>
      dplyr::mutate(method = label)
  }

  dplyr::bind_rows(
    bind_method(linear_rates, "Linear interpolation"),
    bind_method(crs_rates,    "CRS (Appleby & Oldfield 1978)"),
    bind_method(cic_rates,    "CIC (Krishnaswami et al. 1971)")
  )
}

# ── 6. Plots ──────────────────────────────────────────────────────────────────

METHOD_COLOURS <- c(
  "Linear interpolation"         = "#6baed6",
  "CRS (Appleby & Oldfield 1978)"    = "#2c7a4b",
  "CIC (Krishnaswami et al. 1971)"   = "#e07b00"
)

plot_pb210_age_models <- function(pb210_ages, linear_cores) {
  # Build long-format age table for all three methods
  pb_long <- pb210_ages |>
    tidyr::pivot_longer(
      cols      = c(age_crs, age_cic),
      names_to  = "method",
      values_to = "age"
    ) |>
    dplyr::mutate(method = dplyr::recode(method,
      age_crs = "CRS (Appleby & Oldfield 1978)",
      age_cic = "CIC (Krishnaswami et al. 1971)"
    ))

  lin_long <- linear_cores |>
    dplyr::filter(!is.na(age)) |>
    dplyr::select(core_id, depth_mid_cm = depth_mid, age) |>
    dplyr::mutate(method = "Linear interpolation")

  all_ages <- dplyr::bind_rows(pb_long, lin_long)

  # CRS uncertainty ribbon (±1 SD)
  crs_ribbon <- pb210_ages |>
    dplyr::filter(!is.na(age_crs), !is.na(age_crs_sd)) |>
    dplyr::mutate(
      age_lo = age_crs - age_crs_sd,
      age_hi = age_crs + age_crs_sd
    )

  p <- ggplot2::ggplot(all_ages, ggplot2::aes(x = age, y = depth_mid_cm,
                                              colour = method, group = method)) +
    ggplot2::geom_ribbon(
      data        = crs_ribbon,
      ggplot2::aes(x = age_crs, xmin = age_lo, xmax = age_hi, y = depth_mid_cm),
      fill        = "#2c7a4b", alpha = 0.15, inherit.aes = FALSE
    ) +
    ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
    ggplot2::geom_point(size = 2, alpha = 0.8, na.rm = TRUE) +
    ggplot2::scale_y_reverse(name = "Decompacted depth (cm)") +
    ggplot2::scale_x_continuous(name = "Age (years BP)") +
    ggplot2::scale_colour_manual(values = METHOD_COLOURS, name = NULL) +
    ggplot2::facet_wrap(~core_id, scales = "free_x") +
    ggplot2::labs(
      title    = "Age-depth model comparison",
      subtitle = "Shaded band = CRS ±1 SD"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom",
                   legend.text = ggplot2::element_text(size = 9))

  p
}

plot_seq_rate_comparison <- function(comparison_df) {
  if (is.null(comparison_df) || nrow(comparison_df) == 0) {
    return(ggplot2::ggplot() +
      ggplot2::labs(title = "No sequestration rates to compare"))
  }

  ggplot2::ggplot(comparison_df,
    ggplot2::aes(x = core_id, y = seq_rate_g_m2_yr, fill = method)) +
    ggplot2::geom_col(position = "dodge", width = 0.7, alpha = 0.9) +
    ggplot2::geom_point(
      ggplot2::aes(y = seq_rate_wc_g_m2_yr, group = method),
      position    = ggplot2::position_dodge(width = 0.7),
      shape       = 21,
      fill        = "white",
      colour      = "grey30",
      size        = 2.5,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = METHOD_COLOURS, name = NULL) +
    ggplot2::labs(
      title    = "Sequestration rate comparison across age models",
      subtitle = "Bar = timeframe rate  ·  Circle = whole-core rate",
      x        = NULL,
      y        = "Sequestration rate (g C/m²/yr)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   legend.text = ggplot2::element_text(size = 9))
}
