library(targets)
library(tarchetypes)

# Sequestration rates pipeline
# ─────────────────────────────────────────────────────────────────────────────
# Estimates mean OC sequestration rates (g C/m²/yr) for cores that have
# radiometric chronology data (^210Pb or ^14C age-depth control points).
#
# Run:
#   targets::tar_make(script = "_targets_seqrates.R", store = "_targets_seqrates")
#   Shortcut: tmsr() (defined in .Rprofile)
#
# Requires:
#   install.packages("BlueCarbon")
#   remotes::install_github("paleolimbot/pb210")   # for CRS/CIC age models
#   Pre-Analysis Data Preparation/data_raw/core_chronology.csv
#   Pre-Analysis Data Preparation/data_raw/core_pb210.csv   (raw 210Pb activities)
#
# DATA FORMAT — core_chronology.csv:
#   core_id       — must match IDs in core_locations.csv
#   depth_cm      — depth of the dated horizon (cm), surface = 0
#   age_ybp       — age in years before present (0 = modern surface)
#   dating_method — informational label (Pb210, C14, etc.) — not used in calcs
#
# DATA FORMAT — core_pb210.csv (for CRS/CIC models):
#   core_id, depth_top_cm, depth_bottom_cm
#   pb210_total_Bq_kg, pb210_total_sd_Bq_kg
#   pb210_supported_Bq_kg, pb210_supported_sd_Bq_kg
#   dry_mass_g, core_area_cm2
#   See R/pb210_methods.R for full documentation.
#
# WHICH CORES TO DATE:
#   Only cores with max age ≥ timeframe (default 100 yr) produce a rate.
#   Cores without chronology data or with max age < timeframe are skipped.

tar_option_set(
  packages = c("dplyr", "readr", "tidyr", "ggplot2", "BlueCarbon", "pb210")
)

tar_source("R/")

data_raw <- "Pre-Analysis Data Preparation/data_raw"

list(

  # ── Input file tracking ─────────────────────────────────────────────────────
  tar_target(sr_locations_file,   file.path(data_raw, "core_locations.csv"),  format = "file"),
  tar_target(sr_samples_file,     file.path(data_raw, "core_samples.csv"),    format = "file"),
  tar_target(sr_compaction_file,  file.path(data_raw, "core_compaction.csv"), format = "file"),
  tar_target(sr_chronology_file,  file.path(data_raw, "core_chronology.csv"), format = "file"),
  tar_target(sr_pb210_file,       file.path(data_raw, "core_pb210.csv"),      format = "file"),

  # ── METHOD A: Linear interpolation (default) ─────────────────────────────
  # Applies compaction correction, assigns sample ages by linear interpolation
  # of chronology anchor points from core_chronology.csv.
  tar_target(sr_cores_with_ages,
    load_and_assign_ages(
      sr_samples_file,
      sr_compaction_file,
      sr_chronology_file
    )),

  tar_target(sr_rates,
    estimate_sequestration_rates(sr_cores_with_ages, timeframe = 100)),

  tar_target(sr_summary,
    summarise_seq_rates(sr_rates, sr_locations_file)),

  tar_target(sr_plot,
    plot_seq_rates(sr_rates, sr_locations_file)),

  # ── METHODS B & C: CRS and CIC from raw 210Pb activity measurements ────────
  # Fits Constant Rate of Supply (CRS) and Constant Initial Concentration (CIC)
  # age-depth models directly from core_pb210.csv, then re-estimates
  # sequestration rates using the same BlueCarbon method.
  # Requires: remotes::install_github("paleolimbot/pb210")
  tar_target(sr_pb210_data,
    load_pb210_data(sr_pb210_file)),

  tar_target(sr_pb210_models,
    fit_pb210_age_models(sr_pb210_data)),

  tar_target(sr_pb210_ages,
    extract_pb210_ages(sr_pb210_models)),

  tar_target(sr_pb210_age_plot,
    plot_pb210_age_models(sr_pb210_ages, sr_cores_with_ages)),

  # Assign CRS/CIC ages to the decompacted SOC samples (same structure as
  # sr_cores_with_ages, just with model-derived ages replacing linear interp)
  tar_target(sr_crs_cores,
    assign_pb210_ages_to_cores(sr_cores_with_ages, sr_pb210_ages, method = "crs")),

  tar_target(sr_cic_cores,
    assign_pb210_ages_to_cores(sr_cores_with_ages, sr_pb210_ages, method = "cic")),

  tar_target(sr_crs_rates,
    estimate_sequestration_rates(sr_crs_cores, timeframe = 100)),

  tar_target(sr_cic_rates,
    estimate_sequestration_rates(sr_cic_cores, timeframe = 100)),

  # ── Comparison ──────────────────────────────────────────────────────────────
  tar_target(sr_method_comparison,
    compare_seq_rate_methods(sr_rates, sr_crs_rates, sr_cic_rates)),

  tar_target(sr_method_comparison_plot,
    plot_seq_rate_comparison(sr_method_comparison)),

  # ── Report ───────────────────────────────────────────────────────────────────
  tar_quarto(report_seqrates, path = "reports/seqrates.qmd")

)
