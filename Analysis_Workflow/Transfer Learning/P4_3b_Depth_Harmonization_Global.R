# ============================================================================
# MODULE 03b: GLOBAL DATA HARMONIZATION (JANOUSEK / GENERIC)
# ============================================================================
# PURPOSE:   1. Merge Global Core Locations (Metadata) with Sample Data (Depth).
#            2. Harmonize profiles to VM0033 depths (7.5, 22.5, 40, 75cm).
#            3. Apply the SAME Hybrid Spline/Decay logic as Local Data.
#            4. Export for GEE Covariate Extraction.
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
if (file.exists("Analysis_Workflow/blue_carbon_config.R")) source("Analysis_Workflow/blue_carbon_config.R")

# Create directories
dir.create("data_processed", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

# Logging
log_file <- file.path("logs", paste0("global_harmonization_", Sys.Date(), ".log"))
log_msg <- function(msg) {
  entry <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg)
  cat(entry, "\n"); cat(entry, "\n", file=log_file, append=TRUE)
}

log_msg("=== STARTED GLOBAL DATA HARMONIZATION ===")

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(stringr)
})

# Configuration (MUST MATCH LOCAL SETTINGS)
STANDARD_DEPTHS <- c(7.5, 22.5, 40, 75)
VM0033_THICKNESS <- c(15, 15, 20, 50)

# ============================================================================
# 1. LOAD RAW FILES
# ============================================================================

# You provided filenames in the prompt, assuming they are in 'data_global/'
# DATA_GLOBAL_DIR set in blue_carbon_config.R → Pre-Analysis Data Preparation/data_global/
loc_file <- list.files(DATA_GLOBAL_DIR, pattern = "Global_Core.*Locations*\\.csv", full.names = TRUE)
dat_file <- list.files(DATA_GLOBAL_DIR, pattern = "Global_Core.*Samples*\\.csv", full.names = TRUE)

if (length(loc_file) == 0 || length(dat_file) == 0) {
  stop(sprintf("Could not find Global CSV files in '%s'.\n  Expected: Global_Core_Locations.csv and Global_Core_Samples.csv",
               DATA_GLOBAL_DIR))
}

log_msg("Loading Global Data...")
raw_locations <- read_csv(loc_file[1], show_col_types = FALSE)
raw_samples   <- read_csv(dat_file[1], show_col_types = FALSE)

# Standardize Column Names (Lower case)
names(raw_locations) <- tolower(names(raw_locations))
names(raw_samples)   <- tolower(names(raw_samples))

# ============================================================================
# 2. CLEAN & MERGE (TYPE-SAFE VERSION)
# ============================================================================
log_msg("Cleaning and Merging Data...")

# A. Prepare Locations (Force IDs to Character)
locations_clean <- raw_locations %>%
  mutate(
    study_id = as.character(study_id),
    core_id = as.character(core_id)
  ) %>%
  select(study_id, core_id, 
         latitude, longitude, 
         ecosystem, state, country = estuary_id) %>% 
  distinct(study_id, core_id, .keep_all = TRUE)

# B. Prepare Samples (With Safety Checks & Type Conversion)
# 1. Define potential column names
required_cols <- c("study_id", "studyid", 
                   "depth_min", "depth_top_cm", 
                   "depth_max", "depth_bottom_cm", 
                   "soc_g_kg", "soc_percent", 
                   "bd_g_cm3", "bulk_density")

# 2. Safety Loop: Create missing columns with NA
samples_prep <- raw_samples
for (col in required_cols) {
  if (!col %in% names(samples_prep)) {
    samples_prep[[col]] <- NA
  }
}

# 3. Apply Logic
samples_clean <- samples_prep %>%
  mutate(
    # ID Logic: Force to character immediately
    final_study_id = coalesce(as.character(study_id), as.character(studyid)),
    core_id = as.character(core_id),
    
    # Depth Logic
    depth_top = coalesce(as.numeric(depth_min), as.numeric(depth_top_cm)),
    depth_bottom = coalesce(as.numeric(depth_max), as.numeric(depth_bottom_cm)),
    depth_cm = (depth_top + depth_bottom) / 2,
    
    # Carbon Logic (suppressWarnings allows "N/A" text to become NA without clogging console)
    soc_g_kg_final = case_when(
      !is.na(soc_g_kg) ~ suppressWarnings(as.numeric(soc_g_kg)),
      !is.na(soc_percent) ~ suppressWarnings(as.numeric(soc_percent)) * 10,
      TRUE ~ NA_real_
    ),
    
    # Bulk Density Logic
    bd_final = case_when(
      !is.na(bd_g_cm3) ~ suppressWarnings(as.numeric(bd_g_cm3)),
      !is.na(bulk_density) ~ suppressWarnings(as.numeric(bulk_density)),
      TRUE ~ NA_real_
    )
  ) %>%
  select(study_id = final_study_id, 
         core_id, 
         depth_cm, 
         soc_g_kg = soc_g_kg_final, 
         bulk_density_g_cm3 = bd_final) %>%
  filter(!is.na(depth_cm)) 

# C. Merge
# Now both sides have character IDs, so this will succeed
global_merged <- left_join(samples_clean, locations_clean, by = c("study_id", "core_id")) %>%
  filter(!is.na(latitude), !is.na(longitude)) %>% 
  filter(!is.na(soc_g_kg) | !is.na(bulk_density_g_cm3))

# Create unique ID
global_merged$unique_id <- paste(global_merged$study_id, global_merged$core_id, sep="_")

log_msg(sprintf("Merged Dataset: %d samples from %d unique cores.", 
                nrow(global_merged), n_distinct(global_merged$unique_id)))
# ============================================================================
# 3. HARMONIZATION FUNCTION (MATCHING LOCAL LOGIC)
# ============================================================================

# Copying the exact logic from Module 03a
fit_hybrid_profile <- function(depths, values, target_depths) {
  if (length(depths) < 2) return(NULL)
  
  max_measured <- max(depths)
  
  # 1. Spline (Interpolation)
  tryCatch({
    spline_fun <- splinefun(x = depths, y = values, method = "monoH.FC")
  }, error = function(e) return(NULL))
  
  # 2. Decay (Extrapolation)
  decay_model <- NULL
  # Only fit decay if we have enough points and it's actually decreasing
  is_decreasing <- cor(depths, values, method="spearman", use="complete.obs") < -0.3
  if (!is.na(is_decreasing) && is_decreasing && length(depths) >= 3) {
    try({ decay_model <- lm(log(values + 0.1) ~ depths) }, silent=TRUE)
  }
  
  preds <- numeric(length(target_depths))
  
  for (i in seq_along(target_depths)) {
    d <- target_depths[i]
    if (d <= max_measured) {
      preds[i] <- spline_fun(d)
    } else {
      # Extrapolate (Limit 2.5x)
      if (d > max_measured * 2.5) {
        preds[i] <- NA
      } else if (!is.null(decay_model) && coef(decay_model)[2] < 0) {
        preds[i] <- exp(predict(decay_model, newdata = data.frame(depths = d)))
      } else {
        preds[i] <- values[which.max(depths)] # Constant
      }
    }
  }
  return(pmax(0, preds))
}

process_global_core <- function(core_df) {
  # SOC Fit
  if (sum(!is.na(core_df$soc_g_kg)) < 2) return(NULL)
  soc_pred <- fit_hybrid_profile(core_df$depth_cm, core_df$soc_g_kg, STANDARD_DEPTHS)
  
  # BD Fit (Use mean if missing, or spline if available)
  if (sum(!is.na(core_df$bulk_density_g_cm3)) >= 2) {
    bd_pred <- fit_hybrid_profile(core_df$depth_cm, core_df$bulk_density_g_cm3, STANDARD_DEPTHS)
  } else {
    # Global data often misses BD. Use a default or mean if available.
    # Safe fallback: 0.5 g/cm3 (approx average for marsh) or mean of core
    bd_val <- mean(core_df$bulk_density_g_cm3, na.rm=TRUE)
    if (is.nan(bd_val)) bd_val <- 0.5 
    bd_pred <- rep(bd_val, 4)
  }
  
  if (is.null(soc_pred)) return(NULL)
  if (is.null(bd_pred)) bd_pred <- rep(0.5, 4) # Final safety net
  
  # Result
  res <- data.frame(
    core_id = unique(core_df$unique_id),
    study_id = unique(core_df$study_id),
    latitude = unique(core_df$latitude),
    longitude = unique(core_df$longitude),
    ecosystem = unique(core_df$ecosystem),
    depth_cm_midpoint = STANDARD_DEPTHS,
    thickness_cm = VM0033_THICKNESS,
    soc_harmonized = soc_pred,
    bd_harmonized = bd_pred
  )
  
  # Calculate Stock
  res$carbon_stock_kg_m2 <- (res$soc_harmonized * res$bd_harmonized * res$thickness_cm) / 100
  
  return(res)
}

# ============================================================================
# 4. EXECUTION LOOP
# ============================================================================

log_msg("Starting harmonization loop (this may take a minute)...")

# Split data by core for processing
core_list <- split(global_merged, global_merged$unique_id)
results <- list()

counter <- 0
for (id in names(core_list)) {
  out <- tryCatch({
    process_global_core(core_list[[id]])
  }, error = function(e) NULL)
  
  if (!is.null(out)) results[[id]] <- out
  
  counter <- counter + 1
  if (counter %% 500 == 0) print(sprintf("Processed %d / %d cores...", counter, length(core_list)))
}

final_global <- bind_rows(results) %>%
  filter(!is.na(carbon_stock_kg_m2))

# ============================================================================
# 5. SAVE FOR GEE
# ============================================================================

# Add Source Tag
final_global$data_source <- "global"

# 1. Save Full Object (RDS)
saveRDS(final_global, "data_processed/global_cores_harmonized_VM0033.rds")

# 2. Save CSV for Google Earth Engine (Colab)
# GEE only needs Location + ID. It doesn't technically need the stock data yet,
# but we keep it to ensure row-matching later.
write_csv(final_global, "data_processed/global_cores_harmonized_VM0033.csv")

log_msg(sprintf("SUCCESS. Harmonized %d Global Cores.", n_distinct(final_global$core_id)))
log_msg(sprintf("Total Data Points: %d", nrow(final_global)))
log_msg("Output saved to: data_processed/global_cores_harmonized_VM0033.csv")


# ============================================================================
# 7. OPTIONAL: GLOBAL VS LOCAL COMPARISON
# ============================================================================
# Checks if local harmonization (Module 03a) has been run.
# If yes, compares distributions to check for domain shift.

local_file <- "data_processed/cores_harmonized_bluecarbon.csv"

if (file.exists(local_file)) {
  
  log_msg("\nLocal data found. Generating comparison metrics...")
  
  # 1. Load and Prep Local Data
  local_data <- read_csv(local_file, show_col_types = FALSE) %>%
    select(core_id, depth_cm_midpoint, 
           carbon_stock_kg_m2, soc_harmonized, bd_harmonized) %>%
    mutate(source = "Local (Site)")
  
  # 2. Prep Global Data (Subset to same columns)
  global_prep <- final_global %>%
    select(core_id, depth_cm_midpoint, 
           carbon_stock_kg_m2, soc_harmonized, bd_harmonized) %>%
    mutate(source = "Global (Database)")
  
  # 3. Combine
  comparison_df <- bind_rows(local_data, global_prep)
  
  # --- A. SUMMARY TABLE ---
  comp_stats <- comparison_df %>%
    group_by(source, depth_cm_midpoint) %>%
    summarise(
      n = n(),
      mean_stock = mean(carbon_stock_kg_m2, na.rm = TRUE),
      sd_stock = sd(carbon_stock_kg_m2, na.rm = TRUE),
      mean_soc = mean(soc_harmonized, na.rm = TRUE),
      mean_bd = mean(bd_harmonized, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(se_stock = sd_stock / sqrt(n)) %>%
    arrange(depth_cm_midpoint, source)
  
  # Save Table
  write_csv(comp_stats, "diagnostics/global_vs_local_summary.csv")
  print(as.data.frame(comp_stats))
  
  # --- B. VISUALIZATION 1: STOCK PROFILES ---
  p_comp_profile <- ggplot(comp_stats, 
                           aes(x = mean_stock, y = depth_cm_midpoint, 
                               color = source, shape = source)) +
    geom_path(linewidth = 1) +
    geom_point(size = 4) +
    geom_errorbarh(aes(xmin = mean_stock - se_stock, 
                       xmax = mean_stock + se_stock), 
                   height = 5) +
    scale_y_reverse(breaks = STANDARD_DEPTHS) +
    scale_color_manual(values = c("Global (Database)" = "gray60", 
                                  "Local (Site)" = "blue")) +
    labs(title = "Carbon Stock Profile Comparison",
         subtitle = "Mean ± Standard Error",
         x = "Carbon Stock (kg C / m²)", y = "Depth (cm)") +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank())
  
  ggsave("outputs/plots/global_summary/global_vs_local_profile.png", 
         p_comp_profile, width = 6, height = 8, bg = "white")
  
  # --- C. VISUALIZATION 2: DISTRIBUTION BOXPLOTS ---
  # This helps identify if your local data is an outlier compared to global
  p_boxplot <- ggplot(comparison_df, 
                      aes(x = factor(depth_cm_midpoint), 
                          y = carbon_stock_kg_m2, 
                          fill = source)) +
    geom_boxplot(outlier.alpha = 0.3, outlier.size = 1) +
    scale_fill_manual(values = c("Global (Database)" = "gray90", 
                                 "Local (Site)" = "lightblue")) +
    labs(title = "Distribution of Carbon Stocks",
         subtitle = "Global Variability vs Local Site Variability",
         x = "Depth (cm)", y = "Carbon Stock (kg C / m²)") +
    theme_minimal() +
    theme(legend.position = "bottom", legend.title = element_blank())
  
  ggsave("outputs/plots/global_summary/global_vs_local_boxplot.png", 
         p_boxplot, width = 8, height = 6, bg = "white")
  
  log_msg("Comparison plots saved to outputs/plots/global_summary/")
  
} else {
  log_msg("Local data not found. Skipping comparison step.")
}

