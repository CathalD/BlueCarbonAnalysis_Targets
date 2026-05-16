# ============================================================================
# MODULE 03a: LOCAL DEPTH HARMONIZATION
# ============================================================================

# --- SETUP ---
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
  # Defaults if config missing
  QC_SOC_MAX <- 550   # g/kg (55%)
  QC_BD_MIN <- 0.05   # g/cm3
  QC_BD_MAX <- 2.0    # g/cm3
  BOOTSTRAP_SEED <- 123
}

# Create directories
dir.create("outputs/plots/by_stratum", recursive = TRUE, showWarnings = FALSE)
dir.create("diagnostics", recursive = TRUE, showWarnings = FALSE)
dir.create("data_processed", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

# Initialize Logging (FIXED FUNCTION DEFINITION)
log_file <- file.path("logs", paste0("depth_harmonization_", Sys.Date(), ".log"))

log_message <- function(msg, level = "INFO") {
  # Now accepts 'level' argument so it won't crash on warnings
  entry <- sprintf("[%s] %s: %s", format(Sys.time(), "%H:%M:%S"), level, msg)
  cat(entry, "\n")
  cat(entry, "\n", file = log_file, append = TRUE)
}

log_message("=== MODULE 03: LOCAL DEPTH HARMONIZATION STARTED ===")

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(splines); library(readr); library(tidyr)
})
# --- LOAD DATA ---
if (!file.exists("data_processed/cores_clean_bluecarbon.rds")) {
  stop("Cleaned data not found. Run Module 01 first.")
}
cores <- readRDS("data_processed/cores_clean_bluecarbon.rds")
cores_clean <- cores %>% filter(qa_pass)

log_message(sprintf("Loaded %d samples from %d cores.", nrow(cores_clean), n_distinct(cores_clean$core_id)))

# Configuration
STANDARD_DEPTHS <- c(7.5, 22.5, 40, 75) # Midpoints
VM0033_THICKNESS <- c(15, 15, 20, 50)   # Thickness

# ============================================================================
# 1. MODELING & DIAGNOSTIC FUNCTIONS
# ============================================================================

#' Calculate fit diagnostics (RMSE, R2) for the spline model
get_diagnostics <- function(core_df, fitted_vals, target_depths) {
  # We need to predict at the ORIGINAL depths to calculate error
  # Re-fit the spline just for diagnostic comparison
  tryCatch({
    spline_fun <- splinefun(x = core_df$depth_cm, y = core_df$soc_g_kg, method = "monoH.FC")
    preds_at_orig <- spline_fun(core_df$depth_cm)
    
    residuals <- core_df$soc_g_kg - preds_at_orig
    rmse <- sqrt(mean(residuals^2))
    r2 <- 1 - sum(residuals^2) / sum((core_df$soc_g_kg - mean(core_df$soc_g_kg))^2)
    
    return(list(rmse = rmse, r2 = r2, n_samples = nrow(core_df)))
  }, error = function(e) return(list(rmse=NA, r2=NA, n_samples=nrow(core_df))))
}

#' Hybrid Fit: Spline (Interp) + Decay (Extrap)
fit_hybrid_profile <- function(depths, values, target_depths) {
  if (length(depths) < 2) return(NULL)
  max_measured <- max(depths)
  
  # 1. Interpolation (Spline)
  spline_fun <- tryCatch(splinefun(x=depths, y=values, method="monoH.FC"), error=function(e) NULL)
  if(is.null(spline_fun)) return(NULL)
  
  # 2. Extrapolation (Decay)
  decay_model <- NULL
  is_decreasing <- cor(depths, values, method="spearman", use="complete.obs") < -0.3
  
  if (!is.na(is_decreasing) && is_decreasing && length(depths) >= 3) {
    try({ decay_model <- lm(log(values + 0.1) ~ depths) }, silent=TRUE)
  }
  
  preds <- numeric(length(target_depths))
  
  for (i in seq_along(target_depths)) {
    d <- target_depths[i]
    if (d <= max_measured) {
      preds[i] <- spline_fun(d) # Interpolate
    } else {
      # Extrapolate (Limit to 2.5x depth)
      if (d > max_measured * 2.5) { 
        preds[i] <- NA 
      } else if (!is.null(decay_model) && coef(decay_model)[2] < 0) {
        preds[i] <- exp(predict(decay_model, newdata=data.frame(depths=d)))
      } else {
        preds[i] <- values[which.max(depths)] # Constant fallback
      }
    }
  }
  return(pmax(0, preds))
}

process_core <- function(core_df) {
  # Fit SOC & BD
  soc_pred <- fit_hybrid_profile(core_df$depth_cm, core_df$soc_g_kg, STANDARD_DEPTHS)
  
  # BD extrapolation: Constant is safer
  bd_spline <- tryCatch(splinefun(core_df$depth_cm, core_df$bulk_density_g_cm3, method="monoH.FC"), error=function(e) NULL)
  if(is.null(bd_spline)) return(NULL)
  
  bd_last <- core_df$bulk_density_g_cm3[which.max(core_df$depth_cm)]
  bd_pred <- sapply(STANDARD_DEPTHS, function(d) {
    if(d <= max(core_df$depth_cm)) bd_spline(d) else bd_last
  })
  
  if (is.null(soc_pred)) return(NULL)
  
  # Calculate Diagnostics
  diag <- get_diagnostics(core_df, soc_pred, STANDARD_DEPTHS)
  
  # Build Result
  res <- data.frame(
    core_id = unique(core_df$core_id)[1],
    stratum = unique(core_df$stratum)[1],
    depth_cm_midpoint = STANDARD_DEPTHS,
    thickness_cm = VM0033_THICKNESS,
    soc_harmonized = soc_pred,
    bd_harmonized = bd_pred,
    latitude = unique(core_df$latitude)[1],
    longitude = unique(core_df$longitude)[1],
    is_extrapolated = STANDARD_DEPTHS > max(core_df$depth_cm),
    rmse = diag$rmse,
    r2 = diag$r2
  )
  
  # Calculate Stock (kg/m2)
  # Corrected Formula: SOC(g/kg) * BD(g/cm3) * Thick(cm) / 100
  res$carbon_stock_kg_m2 <- (res$soc_harmonized * res$bd_harmonized * res$thickness_cm) / 100
  
  return(res)
}

# ============================================================================
# 2. EXECUTION LOOP
# ============================================================================

log_message("Starting harmonization loop...")
results_list <- list()

for (cid in unique(cores_clean$core_id)) {
  sub <- cores_clean %>% filter(core_id == cid) %>% arrange(depth_cm)
  if (nrow(sub) < 2) next
  
  out <- process_core(sub)
  if (!is.null(out)) results_list[[cid]] <- out
}

harmonized_cores <- bind_rows(results_list) %>% 
  filter(!is.na(carbon_stock_kg_m2))

# ============================================================================
# 3. RESTORED QA/QC FLAGS & DIAGNOSTICS
# ============================================================================

log_message("Applying QA/QC flags...")

# 1. Realistic Range Checks
harmonized_cores <- harmonized_cores %>%
  mutate(
    qa_realistic_soc = soc_harmonized >= 0 & soc_harmonized <= QC_SOC_MAX,
    qa_realistic_bd = bd_harmonized >= QC_BD_MIN & bd_harmonized <= QC_BD_MAX,
    qa_realistic = qa_realistic_soc & qa_realistic_bd
  )

# 2. Monotonicity Check (Is carbon decreasing with depth?)
# We calculate correlation per core
monotonicity_flags <- harmonized_cores %>%
  group_by(core_id) %>%
  summarise(
    cor_with_depth = cor(depth_cm_midpoint, soc_harmonized, use="complete.obs"),
    is_monotonic = cor_with_depth < -0.3, # Threshold for "decreasing"
    max_soc = max(soc_harmonized),
    min_soc = min(soc_harmonized)
  ) %>%
  mutate(qa_monotonic = is_monotonic) # Rename for consistency

# Merge flags back
harmonized_cores <- harmonized_cores %>%
  left_join(monotonicity_flags %>% select(core_id, qa_monotonic), by = "core_id")

# 3. Extract Diagnostics Table (One row per core)
diagnostics_df <- harmonized_cores %>%
  distinct(core_id, stratum, rmse, r2, qa_realistic, qa_monotonic)

# Save Diagnostics
saveRDS(diagnostics_df, "diagnostics/harmonization_diagnostics.rds")
write_csv(diagnostics_df, "diagnostics/harmonization_diagnostics.csv")
saveRDS(monotonicity_flags, "diagnostics/monotonicity_summary.rds")

# Count issues
n_unrealistic <- sum(!harmonized_cores$qa_realistic, na.rm=TRUE)
n_non_monotonic <- sum(!monotonicity_flags$qa_monotonic, na.rm=TRUE)

if (n_unrealistic > 0) log_message(sprintf("WARNING: %d unrealistic predictions detected.", n_unrealistic), "WARNING")
if (n_non_monotonic > 0) log_message(sprintf("WARNING: %d cores have non-monotonic profiles.", n_non_monotonic), "WARNING")

# ============================================================================
# 4. VISUALIZATION (Grouped by Stratum)
# ============================================================================

log_message("Generating grouped diagnostic plots...")

strata_list <- unique(harmonized_cores$stratum)

for (strat in strata_list) {
  strat_safe <- gsub("[^a-zA-Z0-9]", "_", strat)
  
  # Data subset
  harm_strat <- harmonized_cores %>% filter(stratum == strat)
  raw_strat <- cores_clean %>% filter(stratum == strat, core_id %in% harm_strat$core_id)
  
  # Dynamic Height
  n_cores <- n_distinct(harm_strat$core_id)
  plot_height <- min(max(5, ceiling(n_cores / 3) * 3), 25)
  
  # A. SOC Profiles
  p_soc <- ggplot() +
    geom_path(data=harm_strat, aes(x=soc_harmonized, y=depth_cm_midpoint), color="blue", alpha=0.6) +
    geom_point(data=harm_strat, aes(x=soc_harmonized, y=depth_cm_midpoint, color=is_extrapolated), size=2) +
    geom_point(data=raw_strat, aes(x=soc_g_kg, y=depth_cm), color="black", shape=1, size=2) +
    scale_y_reverse(limits=c(100,0)) +
    scale_color_manual(values=c("FALSE"="blue", "TRUE"="red"), name="Extrapolated") +
    facet_wrap(~core_id, ncol=4) +
    labs(title=paste(strat, "- SOC Profiles"), x="SOC (g/kg)", y="Depth") + theme_bw()
  
  ggsave(sprintf("outputs/plots/by_stratum/%s_SOC_Profiles.png", strat_safe), p_soc, width=10, height=plot_height, limitsize=FALSE)
  
  # B. Bulk Density
  p_bd <- ggplot() +
    geom_path(data=harm_strat, aes(x=bd_harmonized, y=depth_cm_midpoint), color="darkgreen", alpha=0.6) +
    geom_point(data=harm_strat, aes(x=bd_harmonized, y=depth_cm_midpoint), color="darkgreen", size=2) +
    geom_point(data=raw_strat, aes(x=bulk_density_g_cm3, y=depth_cm), color="black", shape=1) +
    scale_y_reverse(limits=c(100,0)) +
    facet_wrap(~core_id, ncol=4) +
    labs(title=paste(strat, "- Bulk Density"), x="BD (g/cm3)", y="Depth") + theme_bw()
  
  ggsave(sprintf("outputs/plots/by_stratum/%s_BD_Profiles.png", strat_safe), p_bd, width=10, height=plot_height, limitsize=FALSE)
  
  # C. Carbon Stock
  p_stock <- ggplot(harm_strat, aes(y=depth_cm_midpoint, x=carbon_stock_kg_m2)) +
    geom_col(orientation="y", fill="orange", alpha=0.5) +
    scale_y_reverse(limits=c(100,0)) +
    facet_wrap(~core_id, ncol=4) +
    labs(title=paste(strat, "- Carbon Stock"), x="Stock (kg/m2)", y="Depth") + theme_bw()
  
  ggsave(sprintf("outputs/plots/by_stratum/%s_Stock_Profiles.png", strat_safe), p_stock, width=10, height=plot_height, limitsize=FALSE)
}

# ============================================================================
# 5. SAVE & SUMMARY REPORTS
# ============================================================================

# Save Method Metadata (Crucial for Verifiers)
harmonization_metadata <- list(
  method = "Hybrid Spline + Log-Linear Decay",
  standard_depths = STANDARD_DEPTHS,
  vm0033_intervals = VM0033_THICKNESS,
  processing_date = Sys.Date(),
  qc_thresholds = list(max_soc = QC_SOC_MAX, min_bd = QC_BD_MIN),
  n_cores_harmonized = n_distinct(harmonized_cores$core_id)
)
saveRDS(harmonization_metadata, "data_processed/harmonization_metadata.rds")

# Save Final Data
saveRDS(harmonized_cores, "data_processed/cores_harmonized_bluecarbon.rds")
write_csv(harmonized_cores, "data_processed/cores_harmonized_bluecarbon.csv")

# --- CONSOLE SUMMARY ---
cat("\n========================================\n")
cat("MODULE 03 COMPLETE: HARMONIZATION SUMMARY\n")
cat("========================================\n")
cat(sprintf("Cores Processed: %d\n", n_distinct(harmonized_cores$core_id)))
cat(sprintf("Unrealistic Values: %d\n", n_unrealistic))
cat(sprintf("Non-Monotonic Cores: %d\n", n_non_monotonic))

cat("\n--- Mean Carbon Stock (kg/m2) per Layer ---\n")
summary_tab <- harmonized_cores %>%
  filter(qa_realistic) %>%
  group_by(stratum, depth_cm_midpoint) %>%
  summarise(
    mean_stock = mean(carbon_stock_kg_m2, na.rm=TRUE),
    sd_stock = sd(carbon_stock_kg_m2, na.rm=TRUE),
    n = n(),
    .groups = 'drop'
  ) %>%
  mutate(se_stock = sd_stock / sqrt(n))

print(as.data.frame(summary_tab))

cat("\nOutputs saved to: data_processed/ and outputs/plots/by_stratum/\n")

# ============================================================================
# QUICK-ACCESS COPY: BASIC ANALYSIS OUTPUTS
# ============================================================================

log_message("Copying key Module 03a outputs to outputs/Basic_analysis...")
dir.create("outputs/Basic_analysis", recursive = TRUE, showWarnings = FALSE)

basic_module3_files <- c(
  "data_processed/cores_harmonized_bluecarbon.csv" = "basic_03a_harmonized_cores.csv",
  "diagnostics/harmonization_diagnostics.csv" = "basic_03a_harmonization_diagnostics.csv"
)

for (src in names(basic_module3_files)) {
  dst <- file.path("outputs/Basic_analysis", basic_module3_files[[src]])
  if (file.exists(src)) {
    file.copy(src, dst, overwrite = TRUE)
    log_message(sprintf("Copied to Basic_analysis: %s", basename(dst)))
  } else {
    log_message(sprintf("Skipped missing file for Basic_analysis: %s", src), "WARNING")
  }
}

log_message("=== MODULE 03 COMPLETE ===")
