# ============================================================================
# MODULE 02: BLUE CARBON EXPLORATORY DATA ANALYSIS
# ============================================================================
# PURPOSE: Visualize and explore patterns in blue carbon data by stratum
# INPUTS: 
#   - data_processed/cores_clean_bluecarbon.rds
# OUTPUTS: 
#   - outputs/plots/exploratory/ (multiple figures)
#   - data_processed/eda_summary.rds
# ============================================================================

# ============================================================================
# SETUP
# ============================================================================

# Load configuration
# ── PATH RESOLVER ───────────────────────────────────────────────────────────────────
# Ensures working directory is BlueCarbon_Workflow_V1.0/ (project root)
# so all relative data paths (data_raw/, outputs/, etc.) resolve correctly.
local({
  # Method 1: called via source() — detect this script's location
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) {
    setwd(dirname(dirname(dirname(this_file))))
    return()
  }
  # Method 2: running interactively in RStudio — rstudioapi gives active file
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NULL)
    if (!is.null(active) && nchar(active) > 0) {
      setwd(dirname(dirname(dirname(normalizePath(active)))))
      return()
    }
  }
  # Method 3: fallback — warn and let the user set wd manually
  message("⚠ Could not auto-detect project root. Please run:
",
          "  setwd('/path/to/BlueCarbon_Workflow_V1.0')
",
          "before sourcing this script.")
})
# ──────────────────────────────────────────────────────────────────────────
if (file.exists("Analysis_Workflow/blue_carbon_config.R")) {
  source("Analysis_Workflow/blue_carbon_config.R")
} else {
  stop("Configuration file not found. Run 00b_setup_directories.R first.")
}

# Initialize logging
log_file <- file.path("logs", paste0("exploratory_analysis_", Sys.Date(), ".log"))
if (!dir.exists("logs")) dir.create("logs")

log_message <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- sprintf("[%s] %s: %s", timestamp, level, msg)
  cat(log_entry, "\n")
  cat(log_entry, "\n", file = log_file, append = TRUE)
}

log_message("=== MODULE 02: EXPLORATORY DATA ANALYSIS ===")

# Load packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
})

# Resolve namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

log_message("Packages loaded successfully")

# Set ggplot theme
theme_set(theme_minimal(base_size = 12))

# Create output directory
plot_dir <- "outputs/plots/exploratory"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Detect potential outliers using Isolation Forest
#'
#' Uses the isotree package to identify statistical outliers in SOC,
#' bulk density, and depth measurements. Outliers are flagged but not
#' removed - manual review is recommended.
#'
#' @param cores_data Data frame with soc_g_kg, bulk_density_g_cm3, depth_cm
#' @param contamination Expected proportion of outliers (default 0.05 = 5%)
#' @return Data with outlier_flag and outlier_score columns added
#'
#' @details
#'
#' Advantages:
#' - Works well with mixed distributions
#' - Handles multiple variables simultaneously
#' - No assumptions about data distribution
#' - Computationally efficient
#'
#' @examples
#' cores_with_flags <- detect_outliers(cores)
#' outliers <- cores_with_flags %>% filter(outlier_flag)
detect_outliers <- function(cores_data, contamination = 0.05) {

  # Check if isotree package is available
  if (!requireNamespace("isotree", quietly = TRUE)) {
    log_message("isotree package not available - skipping outlier detection",
               level = "WARNING")
    log_message("Install with: install.packages('isotree')", level = "INFO")

    # Return data with dummy flags
    cores_data$outlier_flag <- FALSE
    cores_data$outlier_score <- NA_real_
    return(cores_data)
  }

  # Select numeric columns for outlier detection
  numeric_cols <- c("soc_g_kg", "bulk_density_g_cm3", "depth_cm")

  # Check columns exist
  missing_cols <- setdiff(numeric_cols, names(cores_data))
  if (length(missing_cols) > 0) {
    log_message(sprintf("Missing columns for outlier detection: %s",
                       paste(missing_cols, collapse = ", ")),
               level = "WARNING")
    cores_data$outlier_flag <- FALSE
    cores_data$outlier_score <- NA_real_
    return(cores_data)
  }

  # Remove rows with NA in key columns
  complete_data <- cores_data[complete.cases(cores_data[, numeric_cols]), ]

  if (nrow(complete_data) < 10) {
    log_message("Too few complete cases for outlier detection (need >= 10)",
               level = "WARNING")
    cores_data$outlier_flag <- FALSE
    cores_data$outlier_score <- NA_real_
    return(cores_data)
  }

  # Fit Isolation Forest
  tryCatch({
    iso_model <- isotree::isolation.forest(
      complete_data[, numeric_cols],
      ntrees = 100,
      sample_size = min(256, nrow(complete_data)),
      ndim = length(numeric_cols),
      seed = 42
    )

    # Predict anomaly scores (higher = more anomalous)
    scores <- predict(iso_model, complete_data[, numeric_cols], type = "score")

    # Determine threshold based on contamination rate
    threshold <- quantile(scores, 1 - contamination, na.rm = TRUE)
    outliers <- scores > threshold

    # Add scores and flags to complete data
    complete_data$outlier_score <- scores
    complete_data$outlier_flag <- outliers

    # Merge back to original data
    cores_data$outlier_score <- NA_real_
    cores_data$outlier_flag <- FALSE

    match_idx <- match(
      paste(complete_data$core_id, complete_data$depth_cm),
      paste(cores_data$core_id, cores_data$depth_cm)
    )

    cores_data$outlier_score[match_idx] <- complete_data$outlier_score
    cores_data$outlier_flag[match_idx] <- complete_data$outlier_flag

    n_outliers <- sum(outliers)

    if (n_outliers > 0) {
      log_message(
        sprintf("Detected %d potential outliers (%.1f%% of complete cases)",
               n_outliers, 100 * n_outliers / nrow(complete_data)),
        level = "WARNING"
      )

      # Save outlier report
      outlier_report <- cores_data %>%
        filter(outlier_flag) %>%
        select(core_id, stratum, depth_cm, soc_g_kg, bulk_density_g_cm3,
               outlier_score) %>%
        arrange(desc(outlier_score))

      outlier_file <- file.path("diagnostics/qaqc",
                                sprintf("outliers_%s_%s.csv",
                                       PROJECT_SCENARIO, MONITORING_YEAR))

      if (!dir.exists("diagnostics/qaqc")) {
        dir.create("diagnostics/qaqc", recursive = TRUE)
      }

      write.csv(outlier_report, outlier_file, row.names = FALSE)
      log_message(sprintf("Outlier report saved to: %s", outlier_file),
                 level = "INFO")
    } else {
      log_message("No outliers detected", level = "INFO")
    }

    return(cores_data)

  }, error = function(e) {
    log_message(sprintf("Outlier detection failed: %s", e$message),
               level = "WARNING")
    cores_data$outlier_flag <- FALSE
    cores_data$outlier_score <- NA_real_
    return(cores_data)
  })
}

#' Check depth profile monotonicity
#'
#' Checks if SOC decreases monotonically with depth (typical pattern).
#' Non-monotonic profiles may indicate:
#' - Sampling/measurement errors
#' - Stratified soil layers
#' - Buried organic matter
#' - Natural variability
#'
#' @param core_data Data frame with core_id, depth_cm, soc_g_kg
#' @return Data frame of cores with non-monotonic profiles
#'
#' @details
#' Identifies cores where SOC increases with depth (non-monotonic).
#' This is flagged for review but may be legitimate in some ecosystems
#' (e.g., buried peat layers, stratified deposits).
#'
#' @examples
#' non_monotonic <- check_monotonicity(cores)
#' if (nrow(non_monotonic) > 0) {
#'   # Review these cores
#' }
check_monotonicity <- function(core_data) {

  non_monotonic <- core_data %>%
    group_by(core_id) %>%
    arrange(core_id, depth_cm) %>%
    mutate(
      soc_change = c(NA, diff(soc_g_kg)),
      depth_interval = c(NA, diff(depth_cm))
    ) %>%
    # Flag where SOC increases significantly with depth
    filter(soc_change > 0 & !is.na(soc_change)) %>%
    ungroup()

  if (nrow(non_monotonic) > 0) {
    n_cores_affected <- n_distinct(non_monotonic$core_id)
    log_message(
      sprintf("%d depth intervals in %d cores show non-monotonic SOC profiles",
             nrow(non_monotonic), n_cores_affected),
      level = "INFO"
    )

    # Save non-monotonic report
    monotonic_report <- non_monotonic %>%
      select(core_id, stratum, depth_cm, soc_g_kg, soc_change, depth_interval) %>%
      arrange(core_id, depth_cm)

    monotonic_file <- file.path("diagnostics/qaqc",
                                sprintf("non_monotonic_profiles_%s_%s.csv",
                                       PROJECT_SCENARIO, MONITORING_YEAR))

    if (!dir.exists("diagnostics/qaqc")) {
      dir.create("diagnostics/qaqc", recursive = TRUE)
    }

    write.csv(monotonic_report, monotonic_file, row.names = FALSE)
    log_message(sprintf("Non-monotonic profiles saved to: %s", monotonic_file),
               level = "INFO")
  } else {
    log_message("All profiles are monotonically decreasing", level = "INFO")
  }

  return(non_monotonic)
}

# ============================================================================
# LOAD DATA
# ============================================================================

log_message("Loading cleaned data...")

if (!file.exists("data_processed/cores_clean_bluecarbon.rds")) {
  stop("Cleaned data not found. Run 01_data_prep_bluecarbon.R first.")
}

cores <- readRDS("data_processed/cores_clean_bluecarbon.rds")

log_message(sprintf("Loaded: %d samples from %d cores",
                    nrow(cores),
                    n_distinct(cores$core_id)))

# Load VM0033 compliance data if available
vm0033_compliance <- NULL
if (file.exists("data_processed/vm0033_compliance.rds")) {
  vm0033_compliance <- readRDS("data_processed/vm0033_compliance.rds")
  log_message("Loaded VM0033 compliance data")
}

# Load BD transparency data if available
bd_transparency <- NULL
if (file.exists("data_processed/bd_transparency.rds")) {
  bd_transparency <- readRDS("data_processed/bd_transparency.rds")
  log_message("Loaded BD transparency data")
}

# Load depth completeness data if available
depth_completeness <- NULL
if (file.exists("data_processed/depth_completeness.rds")) {
  depth_completeness <- readRDS("data_processed/depth_completeness.rds")
  log_message("Loaded depth completeness data")
}

# Load HR vs Composite comparison if available
hr_comp_comparison <- NULL
if (file.exists("data_processed/hr_composite_comparison.rds")) {
  hr_comp_comparison <- readRDS("data_processed/hr_composite_comparison.rds")
  log_message("Loaded HR vs Composite comparison data")
}

# Filter to QA-passed samples only
cores_clean <- cores %>%
  filter(qa_pass)

log_message(sprintf("After QA filter: %d samples from %d cores",
                    nrow(cores_clean),
                    n_distinct(cores_clean$core_id)))

# ============================================================================
# DATA QUALITY CHECKS
# ============================================================================

log_message("Running data quality checks...")

# Ensure key numeric columns are actually numeric (guards against CSV-read character columns)
cores_clean <- cores_clean %>%
  mutate(
    depth_cm           = as.numeric(depth_cm),
    soc_g_kg           = as.numeric(soc_g_kg),
    bulk_density_g_cm3 = as.numeric(bulk_density_g_cm3)
  )

# Outlier detection
cores_clean <- detect_outliers(cores_clean, contamination = 0.05)

# Monotonicity check
non_monotonic_cores <- check_monotonicity(cores_clean)

# Save quality check results
qc_summary <- list(
  n_outliers = sum(cores_clean$outlier_flag, na.rm = TRUE),
  n_non_monotonic_intervals = nrow(non_monotonic_cores),
  n_non_monotonic_cores = n_distinct(non_monotonic_cores$core_id),
  outlier_threshold = quantile(cores_clean$outlier_score, 0.95, na.rm = TRUE),
  timestamp = Sys.time()
)

saveRDS(qc_summary, "data_processed/eda_qc_summary.rds")

log_message("Data quality checks complete")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

log_message("Calculating summary statistics...")

# Overall statistics
overall_stats <- cores_clean %>%
  summarise(
    n_cores = n_distinct(core_id),
    n_samples = n(),
    n_strata = n_distinct(stratum),
    mean_soc = mean(soc_g_kg, na.rm = TRUE),
    sd_soc = sd(soc_g_kg, na.rm = TRUE),
    min_soc = min(soc_g_kg, na.rm = TRUE),
    max_soc = max(soc_g_kg, na.rm = TRUE),
    mean_bd = mean(bulk_density_g_cm3, na.rm = TRUE),
    sd_bd = sd(bulk_density_g_cm3, na.rm = TRUE),
    mean_depth = mean(depth_cm, na.rm = TRUE),
    max_depth = max(depth_bottom_cm, na.rm = TRUE)
  )

cat("\n=== OVERALL STATISTICS ===\n")
print(overall_stats)

# Stratum-specific statistics
stratum_stats <- cores_clean %>%
  group_by(stratum) %>%
  summarise(
    n_cores = n_distinct(core_id),
    n_samples = n(),
    mean_soc = mean(soc_g_kg, na.rm = TRUE),
    sd_soc = sd(soc_g_kg, na.rm = TRUE),
    mean_bd = mean(bulk_density_g_cm3, na.rm = TRUE),
    sd_bd = sd(bulk_density_g_cm3, na.rm = TRUE),
    max_depth = sum(interval_thickness_cm/n_cores, na.rm = TRUE),
    mean_depth = mean(depth_cm, na.rm = TRUE),
    mean_carbon_stock = mean(carbon_stock_kg_m2, na.rm = TRUE),
    total_carbon_stoc = sum((carbon_stock_kg_m2/n_cores), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_soc))

cat("\n=== STATISTICS BY STRATUM ===\n")
print(stratum_stats)

# ============================================================================
# PLOT 1: SPATIAL DISTRIBUTION OF CORES
# ============================================================================

log_message("Creating spatial distribution plot...")

# Prepare spatial data with core type and compliance info
spatial_data <- cores_clean %>%
  distinct(core_id, .keep_all = TRUE) %>%
  mutate(
    core_type_clean = case_when(
      tolower(core_type) %in% c("hr", "high-res", "high resolution", "high res") ~ "HR",
      tolower(core_type) %in% c("paired composite", "paired comp", "paired") ~ "Paired Composite",
      tolower(core_type) %in% c("unpaired composite", "unpaired comp", "unpaired", "composite", "comp") ~ "Unpaired Composite",
      TRUE ~ ifelse(is.na(core_type), "Unknown", core_type)
    )
  )

# Add VM0033 compliance status per core if available
if (!is.null(vm0033_compliance) && "stratum_compliance" %in% names(vm0033_compliance)) {
  compliance_lookup <- vm0033_compliance$stratum_compliance %>%
    select(stratum, compliance_status) %>%
    mutate(
      border_color = case_when(
        compliance_status %in% c("EXCELLENT", "GOOD") ~ "#2E7D32",
        compliance_status == "ACCEPTABLE" ~ "#F9A825",
        TRUE ~ "#C62828"
      )
    )

  spatial_data <- spatial_data %>%
    left_join(compliance_lookup, by = "stratum") %>%
    mutate(border_color = ifelse(is.na(border_color), "black", border_color))
} else {
  spatial_data <- spatial_data %>%
    mutate(border_color = "black")
}

# Create enhanced spatial plot
p_spatial <- ggplot(spatial_data,
                    aes(x = longitude, y = latitude,
                        color = stratum,
                        shape = core_type_clean)) +
  geom_point(size = 4, alpha = 0.8, stroke = 1.5,
             aes(fill = stratum), color = spatial_data$border_color) +
  scale_color_manual(values = STRATUM_COLORS) +
  scale_fill_manual(values = STRATUM_COLORS) +
  scale_shape_manual(values = c("HR" = 21,
                                 "Paired Composite" = 22,
                                 "Unpaired Composite" = 23,
                                 "Unknown" = 24)) +
  labs(
    title = "Spatial Distribution of Core Locations",
    subtitle = sprintf("n = %d cores across %d strata\nShape = Core Type, Fill = Stratum, Border = VM0033 Compliance",
                       n_distinct(spatial_data$core_id),
                       n_distinct(spatial_data$stratum)),
    x = "Longitude",
    y = "Latitude",
    fill = "Stratum",
    shape = "Core Type"
  ) +
  guides(color = "none",  # Hide color legend (use fill instead)
         fill = guide_legend(override.aes = list(shape = 21, size = 4)),
         shape = guide_legend(override.aes = list(size = 4))) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14)
  )

ggsave(file.path(plot_dir, "01_spatial_distribution.png"),
       p_spatial, width = FIGURE_WIDTH * 1.2, height = FIGURE_HEIGHT, dpi = FIGURE_DPI)

log_message("Saved: 01_spatial_distribution.png")

# ============================================================================
# PLOT 1B: VM0033 COMPLIANCE DASHBOARD
# ============================================================================

if (!is.null(vm0033_compliance)) {
  log_message("Creating VM0033 compliance dashboard...")

  # Extract stratum-level data
  if ("stratum_compliance" %in% names(vm0033_compliance)) {
    compliance_data <- vm0033_compliance$stratum_compliance %>%
      mutate(
        status_color = case_when(
          compliance_status %in% c("EXCELLENT", "GOOD") ~ "#2E7D32",
          compliance_status == "ACCEPTABLE" ~ "#F9A825",
          compliance_status == "MARGINAL" ~ "#EF6C00",
          TRUE ~ "#C62828"
        )
      )

    # Plot 1: Sample size vs required
    p_sample_size <- ggplot(compliance_data,
                            aes(x = reorder(stratum, -current_n))) +
      geom_col(aes(y = required_n_20pct), fill = "gray70", alpha = 0.5, width = 0.7) +
      geom_col(aes(y = current_n, fill = status_color), alpha = 0.8, width = 0.7) +
      geom_hline(yintercept = VM0033_MIN_CORES, linetype = "dashed",
                 color = "red", size = 0.8) +
      scale_fill_identity() +
      labs(
        title = "VM0033 Sample Size Compliance",
        subtitle = sprintf("Minimum required: %d cores per stratum (red line)\nGray bars show target for 20%% precision",
                          VM0033_MIN_CORES),
        x = "Stratum",
        y = "Number of Cores"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold")
      ) +
      geom_text(aes(y = current_n, label = current_n),
                vjust = -0.5, size = 3.5, fontface = "bold")

    # Plot 2: Achieved precision vs target
    p_precision <- ggplot(compliance_data %>% filter(!is.na(achieved_precision_pct)),
                          aes(x = reorder(stratum, achieved_precision_pct))) +
      geom_col(aes(y = achieved_precision_pct, fill = status_color), alpha = 0.8) +
      geom_hline(yintercept = VM0033_TARGET_PRECISION, linetype = "dashed",
                 color = "#2E7D32", size = 0.8) +
      geom_hline(yintercept = 10, linetype = "dotted",
                 color = "#1565C0", size = 0.6) +
      scale_fill_identity() +
      coord_flip() +
      labs(
        title = "Achieved Precision at 95% Confidence",
        subtitle = "Green line = 20% target (acceptable), Blue line = 10% (excellent)",
        x = "Stratum",
        y = "Achieved Precision (%)"
      ) +
      theme(
        plot.title = element_text(face = "bold")
      ) +
      geom_text(aes(y = achieved_precision_pct,
                    label = sprintf("%.1f%%", achieved_precision_pct)),
                hjust = -0.2, size = 3.5)

    # Combine plots
    p_vm0033_dashboard <- grid.arrange(p_sample_size, p_precision, ncol = 2)

    ggsave(file.path(plot_dir, "01b_vm0033_compliance_dashboard.png"),
           p_vm0033_dashboard, width = FIGURE_WIDTH * 1.8, height = FIGURE_HEIGHT, dpi = FIGURE_DPI)

    log_message("Saved: 01b_vm0033_compliance_dashboard.png")
  }
} else {
  log_message("Skipping VM0033 compliance dashboard (data not available)", "WARNING")
}

# ============================================================================
# PLOT 2: SOC DISTRIBUTION BY STRATUM
# ============================================================================

log_message("Creating SOC distribution plots...")

# Calculate mean and 95% CI for SOC by stratum
soc_stats <- cores_clean %>%
  group_by(stratum) %>%
  summarise(
    mean_soc = mean(soc_g_kg, na.rm = TRUE),
    sd_soc = sd(soc_g_kg, na.rm = TRUE),
    n = n(),
    se = sd_soc / sqrt(n),
    ci_95 = qt(0.975, df = n - 1) * se,
    .groups = "drop"
  )

# Boxplot with 95% CI
p_soc_box <- ggplot(cores_clean, aes(x = reorder(stratum, -soc_g_kg),
                                      y = soc_g_kg,
                                      fill = stratum)) +
  geom_boxplot(alpha = 0.7) +
  geom_errorbar(data = soc_stats,
                aes(x = stratum, y = mean_soc,
                    ymin = mean_soc - ci_95, ymax = mean_soc + ci_95),
                width = 0.3, size = 1, color = "black", alpha = 0.8) +
  geom_point(data = soc_stats, aes(x = stratum, y = mean_soc),
             size = 3, color = "black", shape = 18) +
  scale_fill_manual(values = STRATUM_COLORS) +
  labs(
    title = "SOC Distribution by Stratum",
    subtitle = "Boxes show distribution, black diamonds show mean ± 95% CI",
    x = "Stratum",
    y = "SOC (g/kg)",
    fill = "Stratum"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

# Violin plot
p_soc_violin <- ggplot(cores_clean, aes(x = reorder(stratum, -soc_g_kg), 
                                         y = soc_g_kg, 
                                         fill = stratum)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.3) +
  scale_fill_manual(values = STRATUM_COLORS) +
  labs(
    title = "SOC Distribution by Stratum (Violin Plot)",
    x = "Stratum",
    y = "SOC (g/kg)",
    fill = "Stratum"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

# Combine
p_soc_combined <- grid.arrange(p_soc_box, p_soc_violin, ncol = 2)

ggsave(file.path(plot_dir, "02_soc_distribution_by_stratum.png"),
       p_soc_combined, width = FIGURE_WIDTH * 1.5, height = FIGURE_HEIGHT, dpi = FIGURE_DPI)

log_message("Saved: 02_soc_distribution_by_stratum.png")

# ============================================================================
# PLOT 3: DEPTH PROFILES BY STRATUM
# ============================================================================

log_message("Creating depth profile plots...")

# Calculate mean SOC by depth and stratum with 95% CI
depth_profiles <- cores_clean %>%
  mutate(
    # Coerce depth_cm to numeric in case it was read as character
    depth_cm = as.numeric(depth_cm),
    # Create 2cm depth bins (0-2, 2-4, 4-6, etc.)
    depth_bin = floor(depth_cm / 2) * 2,
    # Create interval labels for clarity
    depth_interval = paste0(depth_bin, "-", depth_bin + 2, " cm")
  ) %>%
  group_by(stratum, depth_bin) %>%
  summarise(
    mean_soc = mean(soc_g_kg, na.rm = TRUE),
    sd_soc = sd(soc_g_kg, na.rm = TRUE),
    n = n(),
    se_soc = sd_soc / sqrt(n),
    ci_95 = qt(0.975, df = n - 1) * se_soc,
    mean_bd = mean(bulk_density_g_cm3, na.rm =TRUE),
    .groups = "drop"
  ) %>%
  arrange(stratum, depth_bin)

p_depth_profiles <- ggplot(depth_profiles, aes(x = mean_soc, y = -depth_bin,
                                               color = stratum, group = stratum)) +
  geom_point(aes(size = n), alpha = 0.6) +
  geom_smooth(method = "lm") +
  #geom_errorbarh(aes(xmin = mean_soc - ci_95, xmax = mean_soc + ci_95),
                 #height = 1, alpha = 0.4) +  # Reduced height to 1cm for 2cm bins
  scale_color_manual(values = STRATUM_COLORS) +
  scale_size_continuous(range = c(2, 6)) +
  scale_y_continuous(
    breaks = seq(0, -max(depth_profiles$depth_bin), by = -10),  # Tick marks every 10cm
    labels = abs(seq(0, -max(depth_profiles$depth_bin), by = -10))  # Positive labels
  ) +
  scale_x_reverse()+
  labs(
    title = "SOC Depth Profiles by Stratum",
    subtitle = "Lines show mean SOC, error bars show 95% CI, point size shows sample size",
    x = "SOC (g/kg)",
    y = "Depth (cm)",
    color = "Stratum",
    size = "n samples"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.minor.y = element_blank()  # Cleaner y-axis
  ) +
  facet_wrap(~stratum, scales = "free_x")

ggsave(file.path(plot_dir, "03_depth_profiles_by_stratum.png"),
       p_depth_profiles, width = FIGURE_WIDTH * 1.5, height = FIGURE_HEIGHT * 1.2, dpi = FIGURE_DPI)
log_message("Saved: 03_depth_profiles_by_stratum.png")

# ============================================================================
# PLOT 4: BULK DENSITY PATTERNS
# ============================================================================

log_message("Creating bulk density plots...")
bulk_density_g_cm3 = depth_profiles$mean_bd
soc_g_kg = depth_profiles$mean_soc
# BD by stratum
p_bd_stratum <- ggplot(depth_profiles, aes(x = reorder(stratum, bulk_density_g_cm3), 
                                         y = bulk_density_g_cm3, 
                                         fill = stratum)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = STRATUM_COLORS) +
  labs(
    title = "Bulk Density by Stratum",
    x = "Stratum",
    y = "Bulk Density (g/cm³)",
    fill = "Stratum"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

# BD vs depth by stratum
p_bd_depth <- ggplot(depth_profiles, aes(x = bulk_density_g_cm3, y = -depth_bin, 
                                      color = stratum)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = STRATUM_COLORS) +
  labs(
    title = "Bulk Density vs Depth by Stratum",
    x = "Bulk Density (g/cm³)",
    y = "Depth (cm)",
    color = "Stratum"
  ) +
  theme(
    plot.title = element_text(face = "bold")
  ) +
  facet_wrap(~stratum, scales = "free_x")

# BD vs SOC
p_bd_soc <- ggplot(depth_profiles, aes(x = soc_g_kg, y = bulk_density_g_cm3, 
                                    color = stratum)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = STRATUM_COLORS) +
  labs(
    title = "Bulk Density vs SOC by Stratum",
    subtitle = "Lines show linear regression",
    x = "SOC (g/kg)",
    y = "Bulk Density (g/cm³)",
    color = "Stratum"
  ) +
  theme(
    plot.title = element_text(face = "bold")
  )

# Combine BD plots
p_bd_combined <- grid.arrange(
  p_bd_stratum, 
  p_bd_soc, 
  ncol = 2
)

ggsave(file.path(plot_dir, "04_bulk_density_patterns.png"),
       p_bd_combined, width = FIGURE_WIDTH * 1.5, height = FIGURE_HEIGHT, dpi = FIGURE_DPI)

ggsave(file.path(plot_dir, "04b_bulk_density_vs_depth.png"),
       p_bd_depth, width = FIGURE_WIDTH * 1.5, height = FIGURE_HEIGHT, dpi = FIGURE_DPI)

log_message("Saved: 04_bulk_density_patterns.png")

# ============================================================================
# PLOT 4B: BULK DENSITY TRANSPARENCY (MEASURED VS DEFAULT)
# ============================================================================

if (!is.null(bd_transparency)) {
  log_message("Creating bulk density transparency plot...")

  if ("bd_summary" %in% names(bd_transparency)) {
    bd_data <- bd_transparency$bd_summary %>%
      select(stratum, n_cores, n_measured, n_default) %>%
      mutate(
        pct_measured = (n_measured / n_cores) * 100,
        pct_default = (n_default / n_cores) * 100
      ) %>%
      pivot_longer(cols = c(pct_measured, pct_default),
                   names_to = "bd_source",
                   values_to = "percentage") %>%
      mutate(
        bd_source = factor(bd_source,
                          levels = c("pct_measured", "pct_default"),
                          labels = c("Measured", "Default"))
      )

    p_bd_transparency <- ggplot(bd_data, aes(x = reorder(stratum, -percentage),
                                              y = percentage,
                                              fill = bd_source)) +
      geom_col(position = "stack", alpha = 0.8) +
      scale_fill_manual(values = c("Measured" = "#2E7D32", "Default" = "#F9A825")) +
      labs(
        title = "Bulk Density Measurement Transparency",
        subtitle = "Green = Field-measured BD, Yellow = Stratum default BD\nMeasured BD reduces uncertainty in carbon stock estimates",
        x = "Stratum",
        y = "Percentage of Cores (%)",
        fill = "BD Source"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"),
        legend.position = "top"
      ) +
      geom_text(aes(label = sprintf("%.0f%%", percentage)),
                position = position_stack(vjust = 0.5),
                size = 3.5, fontface = "bold", color = "white")

    ggsave(file.path(plot_dir, "04c_bd_transparency.png"),
           p_bd_transparency, width = FIGURE_WIDTH, height = FIGURE_HEIGHT, dpi = FIGURE_DPI)

    log_message("Saved: 04c_bd_transparency.png")
  }
} else {
  log_message("Skipping BD transparency plot (data not available)", "WARNING")
}

# ============================================================================
# PLOT 4C: DEPTH PROFILE COMPLETENESS
# ============================================================================

if (!is.null(depth_completeness)) {
  log_message("Creating depth profile completeness plot...")

  if ("core_completeness" %in% names(depth_completeness)) {
    completeness_data <- depth_completeness$core_completeness %>%
      mutate(
        completeness_category = factor(completeness_category,
                                       levels = c("Complete", "Good", "Moderate", "Incomplete")),
        category_color = case_when(
          completeness_category == "Complete" ~ "#2E7D32",
          completeness_category == "Good" ~ "#66BB6A",
          completeness_category == "Moderate" ~ "#F9A825",
          completeness_category == "Incomplete" ~ "#EF6C00"
        )
      )

    # Count by category and stratum
    completeness_summary <- completeness_data %>%
      group_by(stratum, completeness_category, category_color) %>%
      summarise(n_cores = n(), .groups = "drop")

    p_completeness <- ggplot(completeness_summary,
                            aes(x = stratum, y = n_cores,
                                fill = category_color)) +
      geom_col(position = "stack", alpha = 0.8) +
      scale_fill_identity() +
      labs(
        title = "Depth Profile Completeness by Stratum",
        subtitle = "Complete (≥90%), Good (70-89%), Moderate (50-69%), Incomplete (<50%)\nTarget: All cores should be Complete or Good",
        x = "Stratum",
        y = "Number of Cores"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold")
      ) +
      geom_text(aes(label = completeness_category),
                position = position_stack(vjust = 0.5),
                size = 3, fontface = "bold", color = "white")

    # Add gap analysis text
    gap_summary <- completeness_data %>%
      filter(!is.na(depth_gaps)) %>%
      group_by(stratum) %>%
      summarise(
        cores_with_gaps = sum(depth_gaps != "", na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(cores_with_gaps > 0)

    if (nrow(gap_summary) > 0) {
      gap_text <- paste(
        "⚠ Depth Gaps Detected:",
        paste(sprintf("%s: %d cores with gaps",
                     gap_summary$stratum,
                     gap_summary$cores_with_gaps),
             collapse = "; ")
      )
    } else {
      gap_text <- "✓ No significant depth gaps detected"
    }

    # Add annotation
    p_completeness <- p_completeness +
      annotate("text", x = Inf, y = Inf,
               label = gap_text,
               hjust = 1.05, vjust = 1.5,
               size = 3, color = "gray30", fontface = "italic")

    ggsave(file.path(plot_dir, "04d_depth_completeness.png"),
           p_completeness, width = FIGURE_WIDTH * 1.2, height = FIGURE_HEIGHT, dpi = FIGURE_DPI)

    log_message("Saved: 04d_depth_completeness.png")
  }
} else {
  log_message("Skipping depth completeness plot (data not available)", "WARNING")
}

# ============================================================================
# PLOT 5: CARBON STOCK PATTERNS
# ============================================================================

log_message("Creating carbon stock plots...")

# Calculate mean and 95% CI for carbon stock by stratum
stock_stats <- cores_clean %>%
  group_by(stratum) %>%
  summarise(
    mean_stock = mean(carbon_stock_kg_m2, na.rm = TRUE),
    sd_stock = sd(carbon_stock_kg_m2, na.rm = TRUE),
    n = n(),
    se = sd_stock / sqrt(n),
    ci_95 = qt(0.975, df = n - 1) * se,
    .groups = "drop"
  )

# Carbon stock by stratum with 95% CI
p_stock_stratum <- ggplot(cores_clean, aes(x = reorder(stratum, -carbon_stock_kg_m2),
                                           y = carbon_stock_kg_m2,
                                           fill = stratum)) +
  geom_boxplot(alpha = 0.7) +
  geom_errorbar(data = stock_stats,
                aes(x = stratum, y = mean_stock,
                    ymin = mean_stock - ci_95, ymax = mean_stock + ci_95),
                width = 0.3, size = 1, color = "black", alpha = 0.8) +
  geom_point(data = stock_stats, aes(x = stratum, y = mean_stock),
             size = 3, color = "black", shape = 18) +
  scale_fill_manual(values = STRATUM_COLORS) +
  labs(
    title = "Carbon Stock per Sample by Stratum",
    subtitle = "Black diamonds show mean ± 95% CI",
    x = "Stratum",
    y = "Carbon Stock (kg C/m²)",
    fill = "Stratum"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

# Total carbon stock by core
core_totals <- cores_clean %>%
  group_by(core_id, stratum) %>%
  summarise(
    total_stock = sum(carbon_stock_kg_m2, na.rm = TRUE),
    max_depth = max(depth_bottom_cm),
    .groups = "drop"
  )

p_stock_total <- ggplot(core_totals, aes(x = reorder(stratum, -total_stock),
                                         y = total_stock,
                                         fill = stratum)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = STRATUM_COLORS) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  labs(
    title = "Total Carbon Stock by Core",
    subtitle = sprintf("Summed across depth (n = %d cores)", nrow(core_totals)),
    x = "Stratum",
    y = "Total Carbon Stock (kg C/m²)",
    fill = "Stratum"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

# Combine
p_stock_combined <- grid.arrange(p_stock_stratum, p_stock_total, ncol = 2)

ggsave(file.path(plot_dir, "05_carbon_stock_by_stratum.png"),
       p_stock_combined, width = FIGURE_WIDTH * 1.5, height = FIGURE_HEIGHT, dpi = FIGURE_DPI)

log_message("Saved: 05_carbon_stock_by_stratum.png")

# ============================================================================
# PLOT 6: CORE TYPE COMPARISON (HR vs PAIRED COMPOSITE)
# ============================================================================

log_message("Creating core type comparison plots...")

if ("core_type" %in% names(cores_clean) && n_distinct(cores_clean$core_type) > 1) {

  # Standardize core type names for clean plotting
  cores_plot <- cores_clean %>%
    mutate(
      core_type_clean = case_when(
        tolower(core_type) %in% c("hr", "high-res", "high resolution", "high res") ~ "HR",
        tolower(core_type) %in% c("paired composite", "paired comp", "paired") ~ "Paired Composite",
        tolower(core_type) %in% c("unpaired composite", "unpaired comp", "unpaired", "composite", "comp") ~ "Unpaired Composite",
        TRUE ~ core_type
      )
    )

  # Sample count by type and stratum
  type_counts <- cores_plot %>%
    group_by(stratum, core_type_clean) %>%
    summarise(n_cores = n_distinct(core_id), .groups = "drop")

  p_type_count <- ggplot(type_counts, aes(x = stratum, y = n_cores, fill = core_type_clean)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("HR" = "#1565C0",
                                  "Paired Composite" = "#43A047",
                                  "Unpaired Composite" = "#F9A825")) +
    labs(
      title = "Core Count by Type and Stratum",
      x = "Stratum",
      y = "Number of Cores",
      fill = "Core Type"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))

  # SOC comparison by core type
  p_type_soc <- ggplot(cores_plot, aes(x = core_type_clean, y = soc_g_kg, fill = core_type_clean)) +
    geom_boxplot(alpha = 0.8) +
    scale_fill_manual(values = c("HR" = "#1565C0",
                                  "Paired Composite" = "#43A047",
                                  "Unpaired Composite" = "#F9A825")) +
    labs(
      title = "SOC Distribution by Core Type",
      x = "Core Type",
      y = "SOC (g/kg)"
    ) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1)) +
    facet_wrap(~stratum, scales = "free_y")

  # Add statistical test results if available
  if (!is.null(hr_comp_comparison) && "statistical_tests" %in% names(hr_comp_comparison)) {
    test_results <- hr_comp_comparison$statistical_tests %>%
      mutate(
        significance = ifelse(p_value < 0.05, "*", "ns"),
        label = sprintf("p=%.3f %s", p_value, significance)
      )

    # Add p-values to the plot
    p_type_soc <- p_type_soc +
      geom_text(data = test_results,
                aes(x = 1.5, y = Inf, label = label),
                inherit.aes = FALSE,
                vjust = 1.5, hjust = 0.5,
                size = 3, fontface = "italic", color = "gray30")
  }

  p_type_combined <- grid.arrange(p_type_count, p_type_soc, ncol = 2)

  ggsave(file.path(plot_dir, "06_core_type_comparison.png"),
         p_type_combined, width = FIGURE_WIDTH * 1.8, height = FIGURE_HEIGHT * 1.2, dpi = FIGURE_DPI)

  log_message("Saved: 06_core_type_comparison.png")
} else {
  log_message("Skipping core type comparison (single type or missing data)", "WARNING")
}

# ============================================================================
# PLOT 7: DATA QUALITY FLAGS
# ============================================================================

log_message("Creating QA/QC summary plots...")

# Count QA issues
qa_summary <- cores %>%
  summarise(
    total_samples = n(),
    spatial_valid = sum(qa_spatial_valid, na.rm = TRUE),
    depth_valid = sum(qa_depth_valid, na.rm = TRUE),
    soc_valid = sum(qa_soc_valid, na.rm = TRUE),
    bd_valid = sum(qa_bd_valid, na.rm = TRUE),
    stratum_valid = sum(qa_stratum_valid, na.rm = TRUE),
    overall_pass = sum(qa_pass, na.rm = TRUE)
  )

qa_long <- qa_summary %>%
  select(-total_samples) %>%
  pivot_longer(everything(), names_to = "check", values_to = "passed") %>%
  mutate(
    failed = qa_summary$total_samples - passed,
    check = gsub("_", " ", check),
    check = tools::toTitleCase(check)
  )

p_qa <- ggplot(qa_long, aes(x = reorder(check, passed))) +
  geom_col(aes(y = passed), fill = "#2E7D32", alpha = 0.7) +
  geom_col(aes(y = failed), fill = "#C62828", alpha = 0.7) +
  coord_flip() +
  labs(
    title = "QA/QC Summary",
    subtitle = sprintf("Total samples: %d", qa_summary$total_samples),
    x = "",
    y = "Number of Samples"
  ) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(plot_dir, "07_qa_summary.png"),
       p_qa, width = FIGURE_WIDTH, height = FIGURE_HEIGHT * 0.8, dpi = FIGURE_DPI)

log_message("Saved: 07_qa_summary.png")

# ============================================================================
# PLOT 7B: SCENARIO COMPARISON (PROJECT vs BASELINE)
# ============================================================================

if ("scenario_type" %in% names(cores_clean) && n_distinct(cores_clean$scenario_type) > 1) {
  log_message("Creating scenario comparison plots...")

  # Calculate stats by scenario and stratum
  scenario_stats <- cores_clean %>%
    group_by(scenario_type, stratum) %>%
    summarise(
      mean_soc = mean(soc_g_kg, na.rm = TRUE),
      sd_soc = sd(soc_g_kg, na.rm = TRUE),
      n = n(),
      se = sd_soc / sqrt(n),
      ci_95 = qt(0.975, df = n - 1) * se,
      mean_stock = mean(carbon_stock_kg_m2, na.rm = TRUE),
      sd_stock = sd(carbon_stock_kg_m2, na.rm = TRUE),
      se_stock = sd_stock / sqrt(n),
      ci_95_stock = qt(0.975, df = n - 1) * se_stock,
      .groups = "drop"
    )

  # SOC comparison
  p_scenario_soc <- ggplot(scenario_stats,
                           aes(x = stratum, y = mean_soc,
                               fill = scenario_type)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_soc - ci_95, ymax = mean_soc + ci_95),
                  position = position_dodge(width = 0.9),
                  width = 0.3) +
    scale_fill_manual(values = c("PROJECT" = "#2E7D32",
                                  "BASELINE" = "#1565C0",
                                  "CONTROL" = "#F9A825",
                                  "DEGRADED" = "#C62828")) +
    labs(
      title = "SOC by Scenario and Stratum",
      subtitle = "Error bars show 95% CI",
      x = "Stratum",
      y = "Mean SOC (g/kg)",
      fill = "Scenario"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )

  # Carbon stock comparison
  p_scenario_stock <- ggplot(scenario_stats,
                             aes(x = stratum, y = mean_stock,
                                 fill = scenario_type)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_stock - ci_95_stock,
                      ymax = mean_stock + ci_95_stock),
                  position = position_dodge(width = 0.9),
                  width = 0.3) +
    scale_fill_manual(values = c("PROJECT" = "#2E7D32",
                                  "BASELINE" = "#1565C0",
                                  "CONTROL" = "#F9A825",
                                  "DEGRADED" = "#C62828")) +
    labs(
      title = "Carbon Stock by Scenario and Stratum",
      subtitle = "Error bars show 95% CI",
      x = "Stratum",
      y = "Mean Carbon Stock (Mg C/ha)",
      fill = "Scenario"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )

  p_scenario_combined <- grid.arrange(p_scenario_soc, p_scenario_stock, ncol = 2)

  ggsave(file.path(plot_dir, "07b_scenario_comparison.png"),
         p_scenario_combined, width = FIGURE_WIDTH * 1.8, height = FIGURE_HEIGHT, dpi = FIGURE_DPI)

  log_message("Saved: 07b_scenario_comparison.png")
} else {
  log_message("Skipping scenario comparison (single scenario or missing data)", "WARNING")
}

# ============================================================================
# PLOT 7C: CORRELATION MATRIX
# ============================================================================

log_message("Creating correlation matrix...")

# Select numeric variables for correlation
cor_vars <- cores_clean %>%
  select(soc_g_kg, bulk_density_g_cm3, depth_cm, carbon_stock_kg_m2) %>%
  na.omit()

if (nrow(cor_vars) > 0) {
  # Calculate correlation matrix
  cor_matrix <- cor(cor_vars, use = "complete.obs")

  # Convert to long format for ggplot
  cor_long <- as.data.frame(cor_matrix) %>%
    mutate(var1 = rownames(.)) %>%
    pivot_longer(cols = -var1, names_to = "var2", values_to = "correlation")

  # Create correlation heatmap
  p_correlation <- ggplot(cor_long, aes(x = var1, y = var2, fill = correlation)) +
    geom_tile(color = "white", size = 1) +
    geom_text(aes(label = sprintf("%.2f", correlation)),
              color = "black", size = 4, fontface = "bold") +
    scale_fill_gradient2(low = "#C62828", mid = "white", high = "#2E7D32",
                         midpoint = 0, limits = c(-1, 1),
                         name = "Correlation") +
    labs(
      title = "Correlation Matrix: Key Variables",
      subtitle = "SOC, Bulk Density, Depth, Carbon Stock",
      x = "",
      y = ""
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "right"
    ) +
    coord_fixed()

  ggsave(file.path(plot_dir, "07c_correlation_matrix.png"),
         p_correlation, width = FIGURE_WIDTH, height = FIGURE_HEIGHT * 0.9, dpi = FIGURE_DPI)

  log_message("Saved: 07c_correlation_matrix.png")
} else {
  log_message("Skipping correlation matrix (insufficient data)", "WARNING")
}

# ============================================================================
# PLOT 8: SUMMARY STATISTICS TABLE
# ============================================================================

log_message("Creating summary statistics table...")

# Format stratum stats for display
stratum_table <- stratum_stats %>%
  mutate(
    SOC = sprintf("%.1f ± %.1f", mean_soc, sd_soc),
    BD = sprintf("%.2f ± %.2f", mean_bd, sd_bd),
    `C Stock` = sprintf("%.2f", total_carbon_stoc),
    `Core depth (cm)` = sprintf("%.2f", max_depth)
  ) %>%
  select(Stratum = stratum,
         `N Cores` = n_cores,
         `N Samples` = n_samples,
         `Core depth (cm)` = `Core depth (cm)` ,
         `SOC (g/kg)` = SOC,
         `BD (g/cm³)` = BD,
         `C Stock (kg/m²)` = `C Stock`)

# Create table plot
table_grob <- gridExtra::tableGrob(stratum_table, rows = NULL)

p_table <- grid.arrange(
  table_grob,
  top = "Summary Statistics by Stratum\n(Mean ± SD)"
)

ggsave(file.path(plot_dir, "08_summary_table.png"),
       p_table, width = FIGURE_WIDTH * 1.2, height = FIGURE_HEIGHT * 0.8, dpi = FIGURE_DPI)

log_message("Saved: 08_summary_table.png")

# ============================================================================
# SAVE EDA SUMMARY
# ============================================================================

log_message("Saving EDA summary...")

eda_summary <- list(
  overall_stats = overall_stats,
  stratum_stats = stratum_stats,
  core_totals = core_totals,
  qa_summary = qa_summary,
  depth_profiles = depth_profiles,
  processing_date = Sys.Date(),
  n_plots_created = length(list.files(plot_dir, pattern = "\\.png$"))
)

saveRDS(eda_summary, "data_processed/eda_summary.rds")

log_message("Saved: eda_summary.rds")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("MODULE 02 COMPLETE\n")
cat("========================================\n\n")

cat("Exploratory Analysis Summary:\n")
cat("----------------------------------------\n")
cat(sprintf("Cores analyzed: %d\n", n_distinct(cores_clean$core_id)))
cat(sprintf("Samples analyzed: %d\n", nrow(cores_clean)))
cat(sprintf("Strata represented: %d\n", n_distinct(cores_clean$stratum)))
cat(sprintf("\nMean SOC: %.1f ± %.1f g/kg\n", 
            overall_stats$mean_soc, overall_stats$sd_soc))
cat(sprintf("Mean BD: %.2f ± %.2f g/cm³\n", 
            overall_stats$mean_bd, overall_stats$sd_bd))

cat("\nSOC by stratum (highest to lowest):\n")
for (i in 1:nrow(stratum_stats)) {
  cat(sprintf("  %s: %.1f g/kg (n=%d cores)\n",
              stratum_stats$stratum[i],
              stratum_stats$mean_soc[i],
              stratum_stats$n_cores[i]))
}

cat(sprintf("\nPlots created: %d\n", eda_summary$n_plots_created))
cat(sprintf("Output directory: %s\n", plot_dir))

cat("\nPlots created:\n")
cat("  Core plots:\n")
cat("    01_spatial_distribution.png - Enhanced with core type shapes & VM0033 borders\n")
cat("    01b_vm0033_compliance_dashboard.png - Sample size & precision analysis\n")
cat("  \n")
cat("  Carbon & SOC:\n")
cat("    02_soc_distribution_by_stratum.png - With 95% CI\n")
cat("    03_depth_profiles_by_stratum.png - With 95% CI\n")
cat("  \n")
cat("  Bulk Density:\n")
cat("    04_bulk_density_patterns.png\n")
cat("    04b_bulk_density_vs_depth.png\n")
cat("    04c_bd_transparency.png - Measured vs default BD\n")
cat("    04d_depth_completeness.png - With gap analysis\n")
cat("  \n")
cat("  Carbon Stocks:\n")
cat("    05_carbon_stock_by_stratum.png - With 95% CI\n")
cat("  \n")
cat("  Comparisons:\n")
cat("    06_core_type_comparison.png - HR vs Paired/Unpaired Composite\n")
cat("  \n")
cat("  Quality & Analysis:\n")
cat("    07_qa_summary.png\n")
cat("    07b_scenario_comparison.png - PROJECT vs BASELINE (if applicable)\n")
cat("    07c_correlation_matrix.png\n")
cat("    08_summary_table.png\n")

cat("\nKey Enhancements:\n")
cat("  ✓ Added VM0033 compliance dashboard\n")
cat("  ✓ All plots show 95% confidence intervals\n")
cat("  ✓ BD transparency visualization (measured vs default)\n")
cat("  ✓ Depth completeness with gap analysis\n")
cat("  ✓ Enhanced spatial plot (shape = core type, color = stratum, border = compliance)\n")
cat("  ✓ Scenario comparison (PROJECT vs BASELINE)\n")
cat("  ✓ Correlation matrix for key variables\n")

cat("\nNext steps:\n")
cat("  1. Review plots in outputs/plots/exploratory/\n")
cat("  2. Check VM0033 compliance and sample size recommendations\n")
cat("  3. Review depth gaps and BD transparency\n")
cat("  4. Run: source('03_depth_harmonization_bluecarbon.R')\n\n")

# ============================================================================
# QUICK-ACCESS COPY: BASIC ANALYSIS OUTPUTS
# ============================================================================

log_message("Copying key Module 02 outputs to outputs/Basic_analysis...")
dir.create("outputs/Basic_analysis", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Basic_analysis/Exploratory_Data_Plots", recursive = TRUE, showWarnings = FALSE)

basic_module2_files <- c(
  "outputs/plots/exploratory/01_spatial_distribution.png" = "Exploratory_Map_of_Sampling_Locations_by_Stratum.png",
  "outputs/plots/exploratory/02_soc_distribution_by_stratum.png" = "Exploratory_SOC_Distribution_by_Stratum.png",
  "outputs/plots/exploratory/03_depth_profiles_by_stratum.png" = "Exploratory_Depth_Profiles_by_Stratum.png",
  "outputs/plots/exploratory/05_carbon_stock_by_stratum.png" = "Exploratory_Carbon_Stock_Comparison_by_Stratum.png",
  "outputs/plots/exploratory/08_summary_table.png" = "Exploratory_Summary_Table_of_Key_Statistics.png"
)

for (src in names(basic_module2_files)) {
  dst <- file.path("outputs/Basic_analysis/Exploratory_Data_Plots", basic_module2_files[[src]])
  if (file.exists(src)) {
    file.copy(src, dst, overwrite = TRUE)
    log_message(sprintf("Copied to Basic_analysis: %s", basename(dst)))
  } else {
    log_message(sprintf("Skipped missing file for Basic_analysis: %s", src), "WARNING")
  }
}

log_message("=== MODULE 02 COMPLETE ===")
