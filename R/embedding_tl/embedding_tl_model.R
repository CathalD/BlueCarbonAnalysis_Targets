# =============================================================================
# R/embedding_tl/embedding_tl_model.R
# Embedding-weighted transfer learning pipeline.
#
# Mirrors the four-stage structure of transfer_learning.R but replaces
# Stage A (Wadoux RF domain classifier) with embedding cosine similarity
# weights from compute_embedding_weights().
#
# Stages:
#   A — Embedding cosine similarity weights  (this file; replaces Wadoux RF)
#   B — Weighted global RF                   (same as Wadoux pipeline)
#   C — Bias correction + LOCO CV            (same)
#   D — Bootstrap uncertainty                (same)
# =============================================================================

# ── 1. Assemble combined TL training data ─────────────────────────────────────
# Thin wrapper around the shared prepare_tl_data() — reuses bridge variable
# logic and scale correction from transfer_learning.R without duplication.
prepare_emb_tl_data <- function(cores_harmonized, global_harmonized,
                                 global_covar_file, local_covar_file) {
  message("[emb-tl] Preparing training data (delegates to prepare_tl_data)...")
  prepare_tl_data(cores_harmonized, global_harmonized,
                  global_covar_file, local_covar_file)
}


# ── 2. Train per-depth embedding-weighted TL models ───────────────────────────
train_emb_tl <- function(tl_data, emb_weights, cfg) {
  suppressPackageStartupMessages({ library(dplyr); library(ranger) })
  set.seed(.TL_SEED)

  depths    <- cfg$VM0033_DEPTH_MIDPOINTS %||% c(7.5, 22.5, 40, 75)
  bv_full   <- tl_data$bridge_full
  bv_red    <- tl_data$bridge_reduced
  g_medians <- tl_data$global_medians
  s_factors <- tl_data$scale_factors

  # emb_weights keyed by core_id for fast lookup
  w_lookup <- setNames(emb_weights$weight,     emb_weights$core_id)
  s_lookup <- setNames(emb_weights$cosine_sim, emb_weights$core_id)

  models <- vector("list", length(depths))
  names(models) <- as.character(depths)

  for (d in depths) {
    dchar <- as.character(d)
    message(sprintf("\n[emb-tl] === Depth %.1f cm ===", d))

    g_data <- tl_data$global |>
      filter(depth_cm_midpoint == d, !is.na(carbon_stock_kg_m2)) |>
      filter(if_all(all_of(bv_red), ~ !is.na(.)))

    l_data <- tl_data$local |>
      filter(depth_cm_midpoint == d, !is.na(carbon_stock_kg_m2)) |>
      filter(if_all(all_of(bv_red), ~ !is.na(.)))

    n_global <- nrow(g_data)
    n_local  <- nrow(l_data)
    n_cores  <- n_distinct(l_data$core_id)

    message(sprintf("[emb-tl]   Global N = %d | Local N = %d (%d cores)",
                    n_global, n_local, n_cores))

    if (n_global < 5) {
      message("[emb-tl]   SKIPPING: fewer than 5 global observations")
      next
    }
    if (n_local < 2) {
      message("[emb-tl]   SKIPPING: fewer than 2 local observations")
      next
    }

    # ── Stage A: Embedding cosine similarity weights ─────────────────────────
    # Look up pre-computed weights by core_id. Any core not in emb_weights
    # (e.g. missing embedding due to cloud) falls back to weight = 0.
    message("[emb-tl]   Stage A: embedding cosine similarity weights...")

    emb_w   <- w_lookup[g_data$core_id]
    emb_w[is.na(emb_w)] <- 0
    emb_sim <- s_lookup[g_data$core_id]
    emb_sim[is.na(emb_sim)] <- 0

    # Normalise so mean(weight) = 1 — consistent with Wadoux convention
    w_mean <- mean(emb_w[emb_w > 0])
    if (w_mean > 0) emb_w <- emb_w / w_mean

    g_data$weight <- emb_w

    eff_n <- sum(emb_w)^2 / sum(emb_w^2)
    message(sprintf("[emb-tl]   Weights [%.2f, %.2f], ESS = %.0f / %d (%.0f%%)",
                    min(emb_w), max(emb_w), eff_n, n_global,
                    100 * eff_n / n_global))

    # ── Stage B: Weighted global RF ──────────────────────────────────────────
    covars <- if (n_local < .N_THRESHOLD_REDUCED) bv_red else bv_full
    message(sprintf("[emb-tl]   Stage B: global RF (%d covariates, %s set)...",
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
    message(sprintf("[emb-tl]   Top predictors: %s",
                    paste(names(head(var_imp, 3)), collapse = ", ")))

    # ── Stage C: Bias estimation + LOCO CV ──────────────────────────────────
    message("[emb-tl]   Stage C: bias estimation + LOCO CV...")

    l_complete <- l_data |>
      filter(if_all(all_of(covars), ~ !is.na(.)))

    l_complete$global_pred <- predict(rf_global, data = l_complete)$predictions
    l_complete$residual    <- l_complete$carbon_stock_kg_m2 - l_complete$global_pred

    bias_mean  <- mean(l_complete$residual)
    local_mean <- mean(l_complete$carbon_stock_kg_m2, na.rm = TRUE)
    message(sprintf("[emb-tl]   Bias: %+.3f kg/m²  (local mean: %.3f)",
                    bias_mean, local_mean))

    unique_cores <- unique(l_complete$core_id)
    cv_rows      <- vector("list", length(unique_cores))
    for (j in seq_along(unique_cores)) {
      hld     <- unique_cores[j]
      train_l <- l_complete[l_complete$core_id != hld, ]
      test_l  <- l_complete[l_complete$core_id == hld, ]
      if (nrow(train_l) < 2) next
      cv_bias        <- mean(train_l$residual)
      test_l$cv_pred <- test_l$global_pred + cv_bias
      cv_rows[[j]]   <- test_l[, c("core_id", "carbon_stock_kg_m2",
                                   "global_pred", "cv_pred")]
    }
    cv_df <- bind_rows(cv_rows)

    r2_tl <- r2_global <- rmse_tl <- rmse_global <- NA_real_
    if (nrow(cv_df) >= 2) {
      ss_tot      <- sum((cv_df$carbon_stock_kg_m2 - mean(cv_df$carbon_stock_kg_m2))^2)
      r2_tl       <- 1 - sum((cv_df$carbon_stock_kg_m2 - cv_df$cv_pred)^2)    / ss_tot
      r2_global   <- 1 - sum((cv_df$carbon_stock_kg_m2 - cv_df$global_pred)^2) / ss_tot
      rmse_tl     <- sqrt(mean((cv_df$carbon_stock_kg_m2 - cv_df$cv_pred)^2))
      rmse_global <- sqrt(mean((cv_df$carbon_stock_kg_m2 - cv_df$global_pred)^2))
      message(sprintf(
        "[emb-tl]   LOCO CV  TL: R²=%.3f RMSE=%.4f  |  Global: R²=%.3f RMSE=%.4f",
        r2_tl, rmse_tl, r2_global, rmse_global
      ))
    } else {
      message("[emb-tl]   LOCO CV: insufficient data")
    }

    # ── Stage D: Bootstrap uncertainty ──────────────────────────────────────
    message("[emb-tl]   Stage D: bootstrap uncertainty...")
    boot_biases  <- replicate(.N_BOOTSTRAP, {
      idx <- sample(nrow(l_complete), replace = TRUE)
      mean(l_complete$residual[idx])
    })
    bias_se      <- sd(boot_biases)
    residual_var <- var(l_complete$residual)
    message(sprintf("[emb-tl]   Bias SE = %.4f  Residual SD = %.4f",
                    bias_se, sqrt(residual_var)))

    # Store domain data for similarity plots (cosine_sim per core)
    domain_data <- g_data |>
      select(core_id, all_of(bv_red)) |>
      mutate(
        weight      = emb_w,
        p_similarity = emb_sim,
        source      = "global"
      ) |>
      bind_rows(
        l_data |>
          select(core_id, all_of(bv_red)) |>
          mutate(weight = NA_real_, p_similarity = NA_real_, source = "local")
      )

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
      domain_data     = domain_data,
      global_medians  = g_medians[names(g_medians) %in% covars],
      scale_factors   = s_factors[names(s_factors) %in% covars],
      method          = "embedding_cosine_similarity_weighted_RF_plus_bias_correction"
    )
  }

  valid <- Filter(Negate(is.null), models)
  if (length(valid) == 0)
    stop("[emb-tl] No models trained — check data completeness.")

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

  message("\n[emb-tl] --- Validation summary ---")
  message(paste(capture.output(print(summary_df, row.names = FALSE)), collapse = "\n"))

  list(models = valid, summary = summary_df)
}


# ── 3. Spatial prediction rasters ─────────────────────────────────────────────
# Delegates to the shared predict_tl_rasters() — model object structure is
# identical so no duplication needed.
predict_emb_tl_rasters <- function(emb_tl_models, covar_file) {
  message("[emb-tl] Predicting rasters (delegates to predict_tl_rasters)...")
  predict_tl_rasters(emb_tl_models, covar_file)
}


# ── 4. Maps and validation plots ──────────────────────────────────────────────
# Delegates to the shared plot_tl_maps() — same raster structure and model
# summary format. The similarity heatmap uses domain_data which is stored
# in the model object with the same schema.
plot_emb_tl_maps <- function(emb_tl_rasters, emb_tl_models, cfg) {
  message("[emb-tl] Building maps (delegates to plot_tl_maps)...")
  plot_tl_maps(emb_tl_rasters, emb_tl_models, cfg)
}
