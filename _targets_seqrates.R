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
#   Pre-Analysis Data Preparation/data_raw/core_chronology.csv
#     (see DATA FORMAT note below)
#
# DATA FORMAT — core_chronology.csv:
#   core_id       — must match IDs in core_locations.csv
#   depth_cm      — depth of the dated horizon (cm), surface = 0
#   age_ybp       — age in years before present (0 = modern surface)
#   dating_method — informational label (Pb210, C14, etc.) — not used in calcs
#
#   Provide at least 2 anchor points per core (surface 0,0 + at least one dated
#   horizon). Ages between anchor points are linearly interpolated; ages beyond
#   the deepest anchor are extrapolated using the nearest slope.
#
#   For a proper age model (Bayesian, with uncertainty), use Bacon or Bchron
#   externally and import the modelled ages in this same CSV format.
#
# WHICH CORES TO DATE:
#   Only cores with max age ≥ timeframe (default 100 yr) produce a rate.
#   Cores without chronology data or with max age < timeframe are skipped.

tar_option_set(
  packages = c("dplyr", "readr", "tidyr", "ggplot2", "BlueCarbon")
)

tar_source("R/")

data_raw <- "Pre-Analysis Data Preparation/data_raw"

list(

  # ── Input file tracking ─────────────────────────────────────────────────────
  tar_target(sr_locations_file,   file.path(data_raw, "core_locations.csv"),  format = "file"),
  tar_target(sr_samples_file,     file.path(data_raw, "core_samples.csv"),    format = "file"),
  tar_target(sr_compaction_file,  file.path(data_raw, "core_compaction.csv"), format = "file"),
  tar_target(sr_chronology_file,  file.path(data_raw, "core_chronology.csv"), format = "file"),

  # ── Step 1: Load samples + assign ages via chronology interpolation ─────────
  # Applies compaction correction (same as BlueCarbon pipeline) then assigns
  # an age to each sample depth by linear interpolation of the chronology
  # anchor points. Cores not in core_chronology.csv are excluded.
  tar_target(sr_cores_with_ages,
    load_and_assign_ages(
      sr_samples_file,
      sr_compaction_file,
      sr_chronology_file
    )),

  # ── Step 2: Estimate sequestration rates ─────────────────────────────────
  # Computes mean OC accumulation rate over the specified timeframe.
  # seq_rate    = accumulated OC in top `timeframe` years / timeframe
  # seq_rate_wc = whole-core rate (accumulated OC / max age in core)
  # Both are converted from g/cm²/yr to g C/m²/yr in the output.
  tar_target(sr_rates,
    estimate_sequestration_rates(sr_cores_with_ages, timeframe = 100)),

  # ── Step 3: Per-stratum summary ─────────────────────────────────────────────
  tar_target(sr_summary,
    summarise_seq_rates(sr_rates, sr_locations_file)),

  # ── Step 4: Plot ─────────────────────────────────────────────────────────────
  tar_target(sr_plot,
    plot_seq_rates(sr_rates, sr_locations_file)),

  # ── Report ───────────────────────────────────────────────────────────────────
  tar_quarto(report_seqrates, path = "reports/seqrates.qmd")

)
