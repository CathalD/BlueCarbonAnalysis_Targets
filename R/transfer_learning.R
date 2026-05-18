# R/transfer_learning.R
# ============================================================================
# PURPOSE: Transfer learning for blue carbon spatial prediction.
#
# Leverages a global coastal wetland database while correcting for local bias
# using Wadoux instance weighting (domain adaptation).
#
# WORKFLOW:
#   harmonize_global_layers() — VM0033-depth harmonization of global data
#   prepare_tl_data()         — join GEE covariates for global and local cores
#   train_tl()                — per-depth: Wadoux → weighted RF → bias → bootstrap CI
#   predict_tl_rasters()      — 4-band GeoTIFF per VM0033 depth
#   plot_tl_maps()            — comparison maps and LOCO CV validation summary
#
# REFERENCES:
#   Wadoux et al. (2021) Sampling design optimisation for soil mapping with
#   machine learning. Geoderma 383, 114708.
# ============================================================================

# ── Module-level constants ────────────────────────────────────────────────────

.BRIDGE_VARS_FULL <- c(
  "elevation_m", "slope", "elevationRelMHW", "twi",
  "dist_to_channel_m", "tidal_flat_prob", "coastal_dist_m",
  "VV_mean", "VH_mean", "VVVH_ratio",
  "B", "G", "R", "NIR", "SWIR1", "SWIR2",
  "NDVI_median", "LSWI_median", "mNDWI_median", "SAVI_median"
)

# Scientifically motivated reduced set for N_local < .N_THRESHOLD_REDUCED:
# vegetation density (NDVI), wetness (LSWI/mNDWI), SAR structure (VV),
# tidal position (elevationRelMHW), and absolute elevation.
.BRIDGE_VARS_REDUCED <- c(
  "NDVI_median", "LSWI_median", "mNDWI_median",
  "VV_mean", "elevation_m", "elevationRelMHW"
)

.N_THRESHOLD_REDUCED <- 15L   # switch to reduced vars when N_local < this
.N_BOOTSTRAP         <- 500L  # bootstrap replicates for bias SE
.CI_LEVEL            <- 0.90  # prediction interval width (VM0033 target)
.TL_SEED             <- 42L


# ── 1. Harmonize global layer data to VM0033 depths ──────────────────────────
harmonize_global_layers <- function(global_layers_file, cfg) {
  suppressPackageStartupMessages({ library(dplyr); library(readr) })

  message("[tl] Reading global layer data...")
  layers <- read_csv(global_layers_file, show_col_types = FALSE)

  n_profiles_raw <- n_distinct(paste(layers$dataset, layers$profile_id))
  message(sprintf("[tl] %d layers from %d profiles (before filtering)",
                  nrow(layers), n_profiles_raw))

  # Reshape into the format expected by harmonize_depths():
  #   core_id, stratum, latitude, longitude, depth_cm, soc_g_kg, bulk_density_g_cm3
  #
  # compound core_id = "dataset_profileid" keeps both attributes in a single
  # key without requiring a separate join column downstream.
  cores_global <- layers |>
    filter(
      !is.na(upper_depth), !is.na(lower_depth),
      !is.na(BDOD),        !is.na(OrgC_pct),
      BDOD     > 0,
      OrgC_pct >= 0
    ) |>
    mutate(
      core_id            = paste(dataset, as.character(profile_id), sep = "_"),
      depth_cm           = (upper_depth + lower_depth) / 2,
      soc_g_kg           = OrgC_pct * 10,   # % → g/kg
      bulk_density_g_cm3 = BDOD,
      stratum            = "global"
    ) |>
    select(core_id, stratum, latitude, longitude, depth_cm, soc_g_kg, bulk_density_g_cm3)

  message(sprintf("[tl] Passing %d complete layers (%d profiles) to harmonize_depths()",
                  nrow(cores_global), n_distinct(cores_global$core_id)))

  harmonize_depths(cores_global, cfg)
}


# ── 2. Assemble combined TL training data ─────────────────────────────────────
prepare_tl_data <- function(cores_harmonized, global_harmonized,
                             global_covar_file, local_covar_file) {
  suppressPackageStartupMessages({
    library(dplyr); library(readr); library(terra); library(sf)
  })

  # ── Global: join harmonized depth data with GEE covariates ──────────────
  message("[tl] Loading global GEE covariates...")
  global_covars <- read_csv(global_covar_file, show_col_types = FALSE) |>
    mutate(
      profile_id = as.character(profile_id),
      core_id    = paste(dataset, profile_id, sep = "_")
    )

  covar_cols <- intersect(
    union(.BRIDGE_VARS_FULL, .BRIDGE_VARS_REDUCED),
    names(global_covars)
  )
  if (length(covar_cols) == 0)
    stop("[tl] No bridge variables found in global covariate file.")

  global_tl <- global_harmonized |>
    inner_join(
      global_covars |> select(core_id, all_of(covar_cols)),
      by = "core_id"
    ) |>
    mutate(data_source = "global")

  message(sprintf("[tl] Global TL dataset: %d rows from %d profiles",
                  nrow(global_tl), n_distinct(global_tl$core_id)))

  # Store median per covariate from global data for raster scale-checking later
  global_medians <- vapply(
    covar_cols,
    function(v) median(global_tl[[v]], na.rm = TRUE),
    numeric(1)
  )

  # ── Local: extract covariate values at core locations ────────────────────
  message("[tl] Extracting local covariate values from raster...")
  local_rast <- rast(local_covar_file)
  local_locs <- cores_harmonized |> distinct(core_id, latitude, longitude)
  local_vect <- vect(local_locs, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  local_vect <- project(local_vect, crs(local_rast))
  covar_vals <- extract(local_rast, local_vect, ID = FALSE)

  avail_local <- intersect(covar_cols, names(covar_vals))
  if (length(avail_local) == 0)
    stop("[tl] No bridge variables found in local covariate raster.")

  local_covar_df <- bind_cols(
    local_locs,
    covar_vals[, avail_local, drop = FALSE]
  )

  # Scale correction: GEE sometimes exports reflectance as integers (×10000).
  # Detect via median ratio and rescale if needed.
  scale_factors <- setNames(rep(1, length(avail_local)), avail_local)
  for (v in avail_local) {
    if (grepl("elevation", v, ignore.case = TRUE)) next
    g_med <- global_medians[[v]]
    l_med <- median(local_covar_df[[v]], na.rm = TRUE)
    if (is.na(g_med) || is.na(l_med)) next
    if (abs(g_med) < 0.001 || abs(l_med) < 0.001) next
    ratio <- abs(l_med / g_med)
    if (ratio > 50 && ratio < 20000) {
      message(sprintf("[tl] Rescaling '%s' by 1/10000 (local/global median ratio = %.0f)", v, ratio))
      local_covar_df[[v]] <- local_covar_df[[v]] / 10000
      scale_factors[[v]]  <- 1 / 10000
    }
  }

  local_tl <- cores_harmonized |>
    left_join(local_covar_df |> select(-latitude, -longitude), by = "core_id") |>
    mutate(data_source = "local")

  message(sprintf("[tl] Local TL dataset: %d rows from %d cores",
                  nrow(local_tl), n_distinct(local_tl$core_id)))

  # ── Available bridge variable sets ───────────────────────────────────────
  both_sides     <- intersect(covar_cols, avail_local)
  bridge_full    <- intersect(.BRIDGE_VARS_FULL,    both_sides)
  bridge_reduced <- intersect(.BRIDGE_VARS_REDUCED, both_sides)

  if (length(bridge_reduced) < 2)
    warning("[tl] Fewer than 2 reduced bridge variables available — TL may be unreliable.")

  message(sprintf("[tl] Bridge variables available: %d full, %d reduced",
                  length(bridge_full), length(bridge_reduced)))

  list(
    global         = global_tl,
    local          = local_tl,
    bridge_full    = bridge_full,
    bridge_reduced = bridge_reduced,
    global_medians = global_medians,
    scale_factors  = scale_factors
  )
}


# ── 3. Train per-depth transfer learning models ───────────────────────────────
train_tl <- function(tl_data, cfg) {
  suppressPackageStartupMessages({ library(dplyr); library(ranger) })
  set.seed(.TL_SEED)

  depths     <- cfg$VM0033_DEPTH_MIDPOINTS %||% c(7.5, 22.5, 40, 75)
  bv_full    <- tl_data$bridge_full
  bv_red     <- tl_data$bridge_reduced
  g_medians  <- tl_data$global_medians
  s_factors  <- tl_data$scale_factors

  models <- vector("list", length(depths))
  names(models) <- as.character(depths)

  for (d in depths) {
    dchar <- as.character(d)
    message(sprintf("\n[tl] === Depth %.1f cm ===", d))

    g_data <- tl_data$global |>
      filter(depth_cm_midpoint == d, !is.na(carbon_stock_kg_m2)) |>
      filter(if_all(all_of(bv_red), ~ !is.na(.)))

    l_data <- tl_data$local |>
      filter(depth_cm_midpoint == d, !is.na(carbon_stock_kg_m2)) |>
      filter(if_all(all_of(bv_red), ~ !is.na(.)))

    n_global <- nrow(g_data)
    n_local  <- nrow(l_data)
    n_cores  <- n_distinct(l_data$core_id)

    message(sprintf("[tl]   Global N = %d | Local N = %d (%d cores)",
                    n_global, n_local, n_cores))

    if (n_global < 5) {
      message("[tl]   SKIPPING: fewer than 5 global observations at this depth")
      next
    }
    if (n_local < 2) {
      message("[tl]   SKIPPING: fewer than 2 local observations at this depth")
      next
    }

    # ── Stage A: Wadoux instance weighting ──────────────────────────────────
    # Train a domain classifier to estimate P(local | covariates) for each
    # global sample. Up-weight global samples that look like the local site.
    message("[tl]   Stage A: Wadoux domain weighting...")

    domain_df <- bind_rows(
      g_data |> select(all_of(bv_red)) |> mutate(is_target = 0L),
      l_data |> select(all_of(bv_red)) |> mutate(is_target = 1L)
    ) |> drop_na()

    rf_domain <- ranger(
      is_target ~ .,
      data          = domain_df,
      num.trees     = 500,
      probability   = TRUE,
      min.node.size = 5,
      seed          = .TL_SEED
    )

    pred_domain <- predict(rf_domain,
                           data = g_data |> select(all_of(bv_red)) |> drop_na())$predictions
    p_target <- if (is.matrix(pred_domain)) {
      col_idx <- which(colnames(pred_domain) == "1")
      if (!length(col_idx)) col_idx <- 2L
      pred_domain[, col_idx]
    } else {
      pred_domain
    }

    p_target      <- pmin(pmax(p_target, 0.01), 0.99)
    wadoux_w      <- p_target / (1 - p_target)
    wadoux_w      <- wadoux_w / mean(wadoux_w)   # normalize so mean weight = 1
    g_data$weight <- wadoux_w

    eff_n <- sum(wadoux_w)^2 / sum(wadoux_w^2)
    message(sprintf("[tl]   Weights [%.2f, %.2f], ESS = %.0f / %d (%.0f%%)",
                    min(wadoux_w), max(wadoux_w), eff_n, n_global,
                    100 * eff_n / n_global))

    # ── Stage B: Weighted global RF ─────────────────────────────────────────
    # Use the reduced covariate set when N_local is small to avoid overfitting
    # the bias correction step.
    covars <- if (n_local < .N_THRESHOLD_REDUCED) bv_red else bv_full
    message(sprintf("[tl]   Stage B: global RF (%d covariates, %s set)...",
                    length(covars),
                    if (identical(covars, bv_red)) "reduced" else "full"))

    g_complete <- g_data |>
      filter(if_all(all_of(covars), ~ !is.na(.)))

    rf_global <- ranger(
      formula       = as.formula(
        paste("carbon_stock_kg_m2 ~", paste(covars, collapse = " + "))
      ),
      data          = g_complete,
      case.weights  = g_complete$weight,
      num.trees     = 1000,
      mtry          = max(2L, floor(length(covars) / 3L)),
      min.node.size = 5,
      importance    = "permutation",
      seed          = .TL_SEED
    )

    var_imp <- sort(rf_global$variable.importance, decreasing = TRUE)
    message(sprintf("[tl]   Top predictors: %s",
                    paste(names(head(var_imp, 3)), collapse = ", ")))

    # ── Stage C: Bias estimation and LOCO cross-validation ──────────────────
    message("[tl]   Stage C: bias estimation + LOCO CV...")

    l_complete <- l_data |>
      filter(if_all(all_of(covars), ~ !is.na(.)))

    l_complete$global_pred <- predict(rf_global, data = l_complete)$predictions
    l_complete$residual    <- l_complete$carbon_stock_kg_m2 - l_complete$global_pred

    bias_mean   <- mean(l_complete$residual)
    local_mean  <- mean(l_complete$carbon_stock_kg_m2, na.rm = TRUE)

    message(sprintf("[tl]   Bias: %+.3f kg/m²  (local mean: %.3f)", bias_mean, local_mean))

    # Leave-one-core-out CV
    unique_cores <- unique(l_complete$core_id)
    cv_rows      <- vector("list", length(unique_cores))
    for (j in seq_along(unique_cores)) {
      hld     <- unique_cores[j]
      train_l <- l_complete[l_complete$core_id != hld, ]
      test_l  <- l_complete[l_complete$core_id == hld, ]
      if (nrow(train_l) < 2) next
      cv_bias          <- mean(train_l$residual)
      test_l$cv_pred   <- test_l$global_pred + cv_bias
      cv_rows[[j]]     <- test_l[, c("core_id", "carbon_stock_kg_m2",
                                     "global_pred", "cv_pred")]
    }
    cv_df <- bind_rows(cv_rows)

    r2_tl     <- NA_real_; r2_global  <- NA_real_
    rmse_tl   <- NA_real_; rmse_global <- NA_real_
    if (nrow(cv_df) >= 2) {
      ss_tot     <- sum((cv_df$carbon_stock_kg_m2 - mean(cv_df$carbon_stock_kg_m2))^2)
      ss_tl      <- sum((cv_df$carbon_stock_kg_m2 - cv_df$cv_pred)^2)
      ss_glob    <- sum((cv_df$carbon_stock_kg_m2 - cv_df$global_pred)^2)
      r2_tl      <- 1 - ss_tl   / ss_tot
      r2_global  <- 1 - ss_glob / ss_tot
      rmse_tl    <- sqrt(mean((cv_df$carbon_stock_kg_m2 - cv_df$cv_pred)^2))
      rmse_global <- sqrt(mean((cv_df$carbon_stock_kg_m2 - cv_df$global_pred)^2))
      message(sprintf(
        "[tl]   LOCO CV  TL: R²=%.3f RMSE=%.4f  |  Global: R²=%.3f RMSE=%.4f",
        r2_tl, rmse_tl, r2_global, rmse_global
      ))
    } else {
      message("[tl]   LOCO CV: insufficient data for cross-validation")
    }

    # ── Stage D: Bootstrap uncertainty on bias correction ───────────────────
    message("[tl]   Stage D: bootstrap uncertainty...")

    boot_biases  <- replicate(.N_BOOTSTRAP, {
      idx <- sample(nrow(l_complete), replace = TRUE)
      mean(l_complete$residual[idx])
    })
    bias_se      <- sd(boot_biases)
    residual_var <- var(l_complete$residual)

    message(sprintf("[tl]   Bias SE = %.4f  Residual SD = %.4f",
                    bias_se, sqrt(residual_var)))

    # ── Assemble model object ────────────────────────────────────────────────
    models[[dchar]] <- list(
      depth_cm        = d,
      global_model    = rf_global,
      predictors      = covars,
      bias_correction = bias_mean,
      bias_se         = bias_se,
      residual_sd     = sqrt(residual_var),
      local_mean      = local_mean,
      n_global        = n_global,
      n_local         = n_local,
      n_covariates    = length(covars),
      cv_r2_tl        = r2_tl,
      cv_rmse_tl      = rmse_tl,
      cv_r2_global    = r2_global,
      cv_rmse_global  = rmse_global,
      var_importance  = var_imp,
      global_medians  = g_medians[names(g_medians) %in% covars],
      scale_factors   = s_factors[names(s_factors) %in% covars],
      method          = "Wadoux_weighted_global_RF_plus_bias_correction"
    )
  }

  valid <- Filter(Negate(is.null), models)
  if (length(valid) == 0)
    stop("[tl] No TL models were trained. Check data completeness.")

  summary_df <- bind_rows(lapply(valid, function(m) {
    data.frame(
      depth_cm        = m$depth_cm,
      n_global        = m$n_global,
      n_local         = m$n_local,
      n_covariates    = m$n_covariates,
      bias_correction = round(m$bias_correction, 4),
      bias_se         = round(m$bias_se,         4),
      residual_sd     = round(m$residual_sd,      4),
      cv_r2_tl        = round(m$cv_r2_tl,         3),
      cv_rmse_tl      = round(m$cv_rmse_tl,        4),
      cv_r2_global    = round(m$cv_r2_global,      3),
      cv_rmse_global  = round(m$cv_rmse_global,    4)
    )
  }))

  message("\n[tl] --- Validation summary ---")
  message(paste(capture.output(print(summary_df, row.names = FALSE)), collapse = "\n"))

  list(models = models, summary = summary_df)
}


# ── 4. Generate spatial prediction rasters ────────────────────────────────────
predict_tl_rasters <- function(tl_models, covar_file) {
  suppressPackageStartupMessages({ library(terra); library(ranger) })

  valid <- Filter(Negate(is.null), tl_models$models)
  if (length(valid) == 0)
    stop("[tl] No trained TL models available for raster prediction.")

  message("[tl] Loading covariate raster...")
  covar_rast <- rast(covar_file)

  depth_rasters <- list()

  for (dchar in names(valid)) {
    m <- valid[[dchar]]
    d <- m$depth_cm
    message(sprintf("[tl] Predicting raster for depth %.1f cm...", d))

    missing_bands <- setdiff(m$predictors, names(covar_rast))
    if (length(missing_bands) > 0) {
      message(sprintf("[tl]   SKIPPING: raster missing bands: %s",
                      paste(missing_bands, collapse = ", ")))
      next
    }

    pred_stack <- covar_rast[[m$predictors]]

    # Apply same scale correction used during point extraction in prepare_tl_data
    for (v in names(m$scale_factors)) {
      sf <- m$scale_factors[[v]]
      if (sf != 1 && v %in% names(pred_stack)) {
        pred_stack[[v]] <- pred_stack[[v]] * sf
      }
    }

    # Band 1: Global Prior — weighted RF prediction, no local correction
    global_prior <- terra::predict(
      pred_stack,
      m$global_model,
      fun   = function(model, newdata, ...) predict(model, data = newdata)$predictions,
      na.rm = TRUE
    )
    names(global_prior) <- "Global_Prior"

    # Band 2: Transfer Final — Global Prior + mean local bias correction
    transfer_final <- global_prior + m$bias_correction
    names(transfer_final) <- "Transfer_Final"

    # Band 3: Local Only — constant raster equal to mean of local observations
    local_only <- global_prior * 0 + m$local_mean
    names(local_only) <- "Local_Only"

    # Band 4: Difference — Transfer Final minus Local Only
    difference <- transfer_final - local_only
    names(difference) <- "Difference"

    depth_label <- gsub("\\.", "_", dchar)
    r <- c(global_prior, transfer_final, local_only, difference)
    names(r) <- paste0("d", depth_label, "_", names(r))
    depth_rasters[[dchar]] <- r

    message(sprintf(
      "[tl]   Global Prior mean: %.3f kg/m²  |  Transfer Final mean: %.3f kg/m² (bias: %+.3f)",
      global(global_prior,    "mean", na.rm = TRUE)[1, 1],
      global(transfer_final,  "mean", na.rm = TRUE)[1, 1],
      m$bias_correction
    ))
  }

  if (length(depth_rasters) == 0)
    stop("[tl] No rasters produced. Check that raster band names match model predictors.")

  combined <- do.call(c, depth_rasters)

  out_dir <- "outputs/transfer"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(out_dir, "tl_carbon_stocks_kg_m2.tif")
  writeRaster(combined, out_path, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  message(sprintf("[tl] Rasters written to %s", out_path))

  combined
}


# ── 5. Maps and validation plots ─────────────────────────────────────────────
plot_tl_maps <- function(tl_rasters, tl_models, cfg) {
  suppressPackageStartupMessages({
    library(terra); library(ggplot2); library(dplyr); library(tidyr)
  })

  depth_labels <- c(
    "7_5"  = "0–15 cm",
    "22_5" = "15–30 cm",
    "40"   = "30–50 cm",
    "75"   = "50–100 cm"
  )

  bands_to_show <- c("Global_Prior", "Transfer_Final")

  map_df <- bind_rows(lapply(names(depth_labels), function(dl) {
    lyr_names <- paste0("d", dl, "_", bands_to_show)
    avail     <- intersect(lyr_names, names(tl_rasters))
    if (!length(avail)) return(NULL)
    bind_rows(lapply(avail, function(ln) {
      df <- as.data.frame(tl_rasters[[ln]], xy = TRUE, na.rm = TRUE)
      names(df)[3] <- "carbon_stock_kg_m2"
      df$depth <- depth_labels[[dl]]
      df$model <- sub(paste0("^d", dl, "_"), "", ln)
      df
    }))
  }))

  p_maps <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = "No raster layers to map",
             size = 5, colour = "grey50")

  if (nrow(map_df) > 0) {
    map_df$model <- factor(map_df$model,
                           levels = c("Global_Prior", "Transfer_Final"),
                           labels = c("Global prior", "Transfer (bias-corrected)"))
    map_df$depth <- factor(map_df$depth, levels = unname(depth_labels))

    p_maps <- ggplot(map_df, aes(x = x, y = y, fill = carbon_stock_kg_m2)) +
      geom_raster() +
      facet_grid(model ~ depth) +
      scale_fill_distiller(name = "kg C/m²", palette = "YlOrRd",
                           direction = 1, na.value = "grey90") +
      coord_equal() +
      theme_bw(base_size = 10) +
      theme(
        axis.title = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 9)
      ) +
      labs(
        title    = "Transfer learning predictions by VM0033 depth interval",
        subtitle = "Rows: model type  |  Columns: depth interval"
      )
  }

  # Validation bar chart: LOCO CV R² and RMSE, TL vs global-only
  summary_df  <- tl_models$summary
  p_validation <- ggplot() + theme_void()

  if (nrow(summary_df) > 0) {
    val_long <- summary_df |>
      select(depth_cm, cv_r2_tl, cv_r2_global, cv_rmse_tl, cv_rmse_global) |>
      pivot_longer(-depth_cm, names_to = "key", values_to = "value") |>
      mutate(
        stat  = ifelse(grepl("r2",   key), "R² (LOCO CV)",    "RMSE kg/m² (LOCO CV)"),
        model = ifelse(grepl("global", key), "Global only",         "Transfer (bias-corrected)")
      ) |>
      filter(!is.na(value))

    if (nrow(val_long) > 0) {
      p_validation <- ggplot(val_long,
                             aes(x = factor(depth_cm), y = value, fill = model)) +
        geom_col(position = "dodge", width = 0.65) +
        facet_wrap(~stat, scales = "free_y") +
        scale_fill_manual(
          values = c("Global only"              = "#d9534f",
                     "Transfer (bias-corrected)" = "#5bc0de")
        ) +
        theme_bw(base_size = 11) +
        labs(
          title = "LOCO cross-validation: transfer learning vs global-only",
          x = "Depth (cm midpoint)", y = NULL, fill = NULL
        )
    }
  }

  list(maps = p_maps, validation = p_validation, summary = summary_df)
}
