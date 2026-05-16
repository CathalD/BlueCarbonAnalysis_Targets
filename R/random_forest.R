# R/random_forest.R
# PURPOSE: Step 3 — Random Forest spatial prediction of carbon stocks.
#
# Four functions, one per targets step:
#   prepare_rf_data()     — extract raster covariates at core locations
#   train_rf()            — 80/20 split, train one RF per depth + total
#   predict_rf_rasters()  — predict across full covariate raster → SpatRaster
#   plot_rf_importance()  — ggplot variable importance (total stock model)
#   plot_rf_maps()        — 4-panel depth map + separate total stock map

# ── 1. Prepare training data ──────────────────────────────────────────────────
prepare_rf_data <- function(cores_harmonized, covar_file) {
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(terra); library(sf)
  })
  message("[rf] Extracting covariate values at core locations...")

  # One column per depth interval + one total column
  depth_stocks <- cores_harmonized |>
    mutate(depth_col = paste0("d", gsub(".", "_", as.character(depth_cm_midpoint),
                                        fixed = TRUE), "cm")) |>
    select(core_id, depth_col, carbon_stock_kg_m2) |>
    pivot_wider(names_from = depth_col, values_from = carbon_stock_kg_m2)

  total_stocks <- cores_harmonized |>
    group_by(core_id) |>
    summarise(total_stock = sum(carbon_stock_kg_m2, na.rm = TRUE), .groups = "drop")

  core_locs <- cores_harmonized |>
    distinct(core_id, stratum, latitude, longitude)

  training_wide <- core_locs |>
    left_join(total_stocks, by = "core_id") |>
    left_join(depth_stocks,  by = "core_id")

  # Extract covariate values at core locations
  covar_rast  <- rast(covar_file)
  cores_vect  <- vect(training_wide, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  cores_vect  <- project(cores_vect, crs(covar_rast))
  covar_vals  <- extract(covar_rast, cores_vect, ID = FALSE)

  result <- cbind(training_wide, covar_vals)

  message(sprintf("[rf] Training data: %d cores x %d predictor bands.",
                  nrow(result), ncol(covar_vals)))
  result
}

# ── 2. Train RF models ────────────────────────────────────────────────────────
train_rf <- function(rf_data) {
  suppressPackageStartupMessages({ library(dplyr) })

  meta_cols     <- c("core_id", "stratum", "latitude", "longitude")
  outcome_cols  <- c("total_stock", grep("^d[0-9]", names(rf_data), value = TRUE))
  predictor_cols <- setdiff(names(rf_data), c(meta_cols, outcome_cols))

  set.seed(42)
  n         <- nrow(rf_data)
  train_idx <- sample(seq_len(n), floor(0.8 * n))
  train_df  <- rf_data[train_idx, ]
  test_df   <- rf_data[-train_idx, ]

  message(sprintf("[rf] Training on %d cores, testing on %d.", length(train_idx),
                  n - length(train_idx)))

  models <- lapply(outcome_cols, function(y) {
    df <- na.omit(train_df[, c(y, predictor_cols)])
    if (nrow(df) < 5) {
      message(sprintf("[rf] Skipping %s — fewer than 5 complete rows.", y))
      return(NULL)
    }
    randomForest::randomForest(
      as.formula(paste(y, "~ .")),
      data      = df,
      ntree     = 500,
      importance = TRUE
    )
  })
  names(models) <- outcome_cols
  models <- Filter(Negate(is.null), models)

  # CV metrics on held-out 20%
  cv_metrics <- bind_rows(lapply(names(models), function(y) {
    test_data <- na.omit(test_df[, c(y, predictor_cols)])
    if (nrow(test_data) == 0)
      return(data.frame(outcome = y, n_test = 0L, rmse = NA_real_, r2 = NA_real_))
    preds  <- predict(models[[y]], test_data)
    actual <- test_data[[y]]
    rmse   <- sqrt(mean((preds - actual)^2))
    r2     <- 1 - sum((preds - actual)^2) / sum((actual - mean(actual))^2)
    data.frame(outcome = y, n_test = nrow(test_data),
               rmse = round(rmse, 3), r2 = round(r2, 3))
  }))

  message("[rf] CV metrics:")
  message(paste(capture.output(print(cv_metrics)), collapse = "\n"))

  list(
    models         = models,
    cv_metrics     = cv_metrics,
    predictor_cols = predictor_cols,
    outcome_cols   = names(models)
  )
}

# ── 3. Predict across full raster ────────────────────────────────────────────
predict_rf_rasters <- function(rf_models, covar_file) {
  suppressPackageStartupMessages({ library(terra) })
  message("[rf] Predicting carbon stocks across full AOI raster...")

  covar_rast <- rast(covar_file)

  pred_layers <- lapply(rf_models$outcome_cols, function(y) {
    message(sprintf("[rf]   Predicting %s...", y))
    terra::predict(covar_rast, rf_models$models[[y]], type = "response", na.rm = TRUE)
  })

  result <- rast(pred_layers)
  names(result) <- rf_models$outcome_cols

  # Write GeoTIFFs to outputs/rf/
  out_dir <- "outputs/rf"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeRaster(result, file.path(out_dir, "rf_carbon_stocks_kg_m2.tif"), overwrite = TRUE)
  message(sprintf("[rf] Rasters written to %s/rf_carbon_stocks_kg_m2.tif", out_dir))

  result
}

# ── 4. Variable importance plot ───────────────────────────────────────────────
plot_rf_importance <- function(rf_models, top_n = 20) {
  suppressPackageStartupMessages({ library(ggplot2) })

  m <- rf_models$models[["total_stock"]]
  if (is.null(m)) {
    m <- rf_models$models[[1]]
    message("[rf] total_stock model not found — using first available model for importance plot.")
  }

  imp <- randomForest::importance(m, type = 1)  # %IncMSE
  imp_df <- data.frame(variable = rownames(imp), importance = imp[, 1]) |>
    (\(d) d[order(-d$importance), ])() |>
    head(top_n)

  ggplot(imp_df, aes(x = reorder(variable, importance), y = importance)) +
    geom_col(fill = "#2c7bb6", alpha = 0.8) +
    coord_flip() +
    theme_bw(base_size = 11) +
    labs(title = "Random Forest variable importance",
         subtitle = "Total stock model — % increase in MSE when variable permuted",
         x = NULL, y = "% increase in MSE")
}

# ── 5. Prediction maps ────────────────────────────────────────────────────────
plot_rf_maps <- function(rf_rasters, cfg) {
  suppressPackageStartupMessages({ library(terra); library(ggplot2); library(dplyr) })

  rast_df <- function(r, lyr) {
    df <- as.data.frame(r[[lyr]], xy = TRUE, na.rm = TRUE)
    names(df)[3] <- "carbon_stock_kg_m2"
    df$layer <- lyr
    df
  }

  depth_cols <- grep("^d[0-9]", rf_rasters$outcome_cols %||% names(rf_rasters), value = TRUE)
  if (length(depth_cols) == 0)
    depth_cols <- grep("^d[0-9]", names(rf_rasters), value = TRUE)

  depth_df <- bind_rows(lapply(depth_cols, rast_df, r = rf_rasters))

  p_depths <- ggplot(depth_df, aes(x = x, y = y, fill = carbon_stock_kg_m2)) +
    geom_raster() +
    facet_wrap(~layer, ncol = 2) +
    scale_fill_viridis_c(name = "kg C/m²", option = "YlOrRd", na.value = "grey90") +
    coord_equal() +
    theme_bw(base_size = 11) +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank()) +
    labs(title = "Predicted carbon stock by VM0033 depth interval")

  total_df <- rast_df(rf_rasters, "total_stock")

  p_total <- ggplot(total_df, aes(x = x, y = y, fill = carbon_stock_kg_m2)) +
    geom_raster() +
    scale_fill_viridis_c(name = "kg C/m²", option = "YlOrRd", na.value = "grey90") +
    coord_equal() +
    theme_bw(base_size = 11) +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank()) +
    labs(title = "Total predicted carbon stock to 1 m depth")

  list(depth_panel = p_depths, total = p_total)
}
