library(targets)
library(tarchetypes)

# BlueCarbon pre-processing pipeline
# ─────────────────────────────────────────────────────────────────────────────
# Runs independently of the main pipeline. Requires:
#   install.packages("BlueCarbon")
#
# Run:
#   targets::tar_make(script = "_targets_bluecarbon.R", store = "_targets_bluecarbon")
#   Shortcut: tmbc() (defined in .Rprofile after adding below)
#
# What this produces:
#   1. Compaction-corrected, gap-filled cores (BlueCarbon method)
#   2. OC stocks estimated via linear extrapolation to 100 cm (BlueCarbon)
#   3. Extrapolation validation plot (how much error do short cores introduce?)
#   4. core_samples_bc_processed.csv — decompacted cores in cores_raw format,
#      ready to feed into harmonize_depths() and the full spatial pipeline
#   5. Side-by-side comparison with the VM0033 pipeline (run main pipeline first)
#
# Integration with spatial pipeline:
#   The file outputs/bluecarbon/core_samples_bc_processed.csv contains the
#   decompacted, gap-filled samples in the same format as core_samples.csv.
#   To run the full spatial pipeline on compaction-corrected data, copy or
#   symlink this file as your core_samples input, then re-run targets::tar_make().

tar_option_set(
  packages = c("dplyr", "readr", "tidyr", "ggplot2", "BlueCarbon")
)

tar_source("R/")

data_raw <- "Pre-Analysis Data Preparation/data_raw"

list(

  # ── Input file tracking ─────────────────────────────────────────────────────
  tar_target(bc_config_file,     "blue_carbon_config.R",                           format = "file"),
  tar_target(bc_locations_file,  file.path(data_raw, "core_locations.csv"),        format = "file"),
  tar_target(bc_samples_file,    file.path(data_raw, "core_samples.csv"),          format = "file"),
  tar_target(bc_compaction_file, file.path(data_raw, "core_compaction.csv"),       format = "file"),

  # ── Configuration ───────────────────────────────────────────────────────────
  tar_target(bc_cfg, load_config(bc_config_file)),

  # ── Step 1: Compaction correction ───────────────────────────────────────────
  # Estimates % compaction per core from field measurements, then corrects
  # sample depths and bulk density. Returns a data frame with:
  #   mind_corrected, maxd_corrected  — decompacted depths
  #   dbd_corrected                   — corrected bulk density
  #   compaction                      — % compaction per core
  #   soc_pct                         — SOC converted to % for BlueCarbon
  tar_target(bc_decompacted,
    run_compaction_correction(bc_samples_file, bc_compaction_file)),

  # ── Step 2: BlueCarbon OC stocks ─────────────────────────────────────────
  # Estimates total OC stock (g/cm² → kg/m²) to 100 cm per core.
  # Short cores are extrapolated via lm(cumulative_OC ~ depth).
  # Output includes stock_kg_m2, stockwc_kg_m2 (whole core), stock_se_kg_m2.
  tar_target(bc_stocks,
    estimate_bc_stocks(bc_locations_file, bc_decompacted, depth = 100)),

  # ── Step 3: Extrapolation validation ────────────────────────────────────────
  # Quantifies error introduced by the linear extrapolation in cores that
  # don't reach 100 cm. Uses the deepest cores as truth, truncates to
  # 90/75/50/25% depth, re-extrapolates, and measures % error.
  # Returns a list: $data (data frame) and $plot (ggplot).
  # NOTE: requires at least one core ≥ 100 cm to work. With current test data
  # (max depth 69.5 cm) this will message "no cores reach target depth" —
  # update depth argument or add deeper cores to test this validation.
  tar_target(bc_extrapolation_test,
    run_extrapolation_test(bc_decompacted, depth = 100)),

  # ── Step 4: Prepare cores for harmonize_depths() ────────────────────────────
  # Converts decompacted samples to cores_raw format:
  #   estimate_h() fills any depth gaps (midpoint split between non-contiguous samples)
  #   layer_thickness_cm = h (gap-filled, not simply maxd - mind)
  #   This output can replace core_samples.csv in the main pipeline.
  tar_target(bc_cores_harmonization_ready,
    prepare_bc_cores_for_harmonization(bc_locations_file, bc_decompacted)),

  # Write processed cores to CSV for use in main / spatial pipeline
  tar_target(bc_cores_csv,
    write_bc_cores_csv(bc_cores_harmonization_ready, "outputs/bluecarbon"),
    format = "file"),

  # ── Step 5: Compare with VM0033 main pipeline ─────────────────────────────
  # Reads stratum_summary from the main _targets store (run targets::tar_make()
  # first). Produces a comparison of BlueCarbon decompacted stocks vs VM0033
  # stocks computed from raw (uncorrected) depths.
  tar_target(bc_vm_comparison,
    compare_bc_vs_vm0033(
      bc_stocks,
      tryCatch(
        targets::tar_read(stratum_summary, store = "_targets"),
        error = function(e) {
          message("[bluecarbon] Main pipeline store not found — skipping comparison.")
          NULL
        }
      )
    )),

  tar_target(bc_comparison_plot,
    plot_bc_comparison(bc_vm_comparison))

)
