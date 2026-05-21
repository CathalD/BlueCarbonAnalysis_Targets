# =============================================================================
# R/embedding_tl/embedding_similarity.R
# Cosine similarity between the local AOI and global core embeddings.
#
# Option 1a (implemented here): AOI-mean embedding → one weight per global
# core.  Straightforward to validate, directly comparable to Wadoux weights.
#
# Option 2 (cluster-based, future): cluster AOI pixels → k centroids → k
# weight vectors → k weighted models. Enabled via compute_pixel_similarity().
# =============================================================================


# -----------------------------------------------------------------------------
# compute_embedding_weights()
# -----------------------------------------------------------------------------
# Computes a Wadoux-compatible weight for each global core based on cosine
# similarity between that core's embedding and the mean AOI embedding.
#
# Steps:
#   1. Average the 64-band AOI raster → one 64-d AOI centroid vector
#   2. Cosine similarity: dot(core_vec_normed, aoi_vec_normed) ∈ [0, 1]
#   3. Sharpening: weight = sim^alpha  (alpha=5 strongly favours similar cores)
#   4. Normalise: mean(weight) = 1  (matches Wadoux normalisation convention)
#
# Returns a data.frame: core_id, profile_id, dataset, cosine_sim, weight.
# -----------------------------------------------------------------------------
compute_embedding_weights <- function(aoi_emb_rast, global_emb_df, alpha = 5) {
  suppressPackageStartupMessages({ library(terra); library(dplyr) })

  emb_cols <- grep("^emb_", names(global_emb_df), value = TRUE)
  if (length(emb_cols) == 0L)
    stop("[EMB] No 'emb_*' columns in global_emb_df — run extract_global_embeddings() first.")

  # ── Step 1: AOI mean embedding vector ─────────────────────────────────────
  message("[EMB] Computing AOI mean embedding vector...")
  aoi_means <- terra::global(aoi_emb_rast[[emb_cols]], fun = "mean", na.rm = TRUE)
  aoi_vec   <- as.numeric(aoi_means[, "mean"])
  aoi_norm  <- sqrt(sum(aoi_vec^2))
  if (aoi_norm == 0) stop("[EMB] AOI mean embedding is a zero vector.")
  aoi_unit  <- aoi_vec / aoi_norm

  # ── Step 2: Per-core cosine similarity ─────────────────────────────────────
  message(sprintf("[EMB] Computing cosine similarity for %d global cores...",
                  nrow(global_emb_df)))

  core_mat   <- as.matrix(global_emb_df[, emb_cols])
  na_rows    <- apply(core_mat, 1, anyNA)
  core_norms <- sqrt(rowSums(core_mat^2))
  core_norms[core_norms == 0 | na_rows] <- NA_real_

  # Normalise rows then dot with AOI unit vector
  core_unit <- sweep(core_mat, 1, core_norms, "/")
  sim        <- as.numeric(core_unit %*% aoi_unit)
  sim[is.na(sim) | na_rows] <- 0
  sim <- pmin(pmax(sim, 0), 1)   # clamp: embeddings should be non-negative, but be safe

  if (any(na_rows))
    message(sprintf("[EMB]   %d cores had NA embeddings → cosine_sim = 0", sum(na_rows)))

  # ── Step 3 & 4: Sharpen and normalise ──────────────────────────────────────
  weights <- sim^alpha
  pos_mean <- mean(weights[weights > 0])
  if (pos_mean > 0) weights <- weights / pos_mean

  result <- global_emb_df |>
    select(profile_id, dataset) |>
    mutate(
      core_id    = paste(dataset, as.character(profile_id), sep = "_"),
      cosine_sim = round(sim, 4),
      weight     = round(weights, 4)
    )

  eff_n <- sum(weights)^2 / sum(weights^2)
  message(sprintf(
    "[EMB] Similarity: min=%.3f  median=%.3f  max=%.3f  (alpha=%d)",
    min(sim), median(sim), max(sim), as.integer(alpha)
  ))
  message(sprintf(
    "[EMB] Weights: ESS = %.0f / %d (%.0f%%) | %d cores with weight > 1",
    eff_n, nrow(result), 100 * eff_n / nrow(result),
    sum(weights > 1)
  ))

  result
}


# -----------------------------------------------------------------------------
# compute_pixel_similarity()
# -----------------------------------------------------------------------------
# For every AOI pixel, compute cosine similarity to every global core.
# Returns a SpatRaster with nlyr = n_cores (memory-intensive for large AOIs).
# Intended for Option 2 (cluster-based weighting) — call with care.
# -----------------------------------------------------------------------------
compute_pixel_similarity <- function(aoi_emb_rast, global_emb_df) {
  suppressPackageStartupMessages({ library(terra); library(dplyr) })

  emb_cols   <- grep("^emb_", names(global_emb_df), value = TRUE)
  aoi_vals   <- terra::values(aoi_emb_rast[[emb_cols]])   # n_pixels × 64
  valid_px   <- complete.cases(aoi_vals)

  aoi_norms  <- sqrt(rowSums(aoi_vals^2))
  aoi_norms[aoi_norms == 0 | !valid_px] <- NA_real_
  aoi_unit   <- aoi_vals / aoi_norms

  core_mat   <- as.matrix(global_emb_df[, emb_cols])
  core_norms <- sqrt(rowSums(core_mat^2))
  core_norms[core_norms == 0] <- NA_real_
  core_unit  <- sweep(core_mat, 1, core_norms, "/")

  n_px    <- nrow(aoi_vals)
  n_cores <- nrow(core_mat)
  message(sprintf("[EMB] Computing %d × %d cosine similarity matrix (~%.0f MB)...",
                  n_px, n_cores, n_px * n_cores * 4 / 1e6))

  sim_mat <- aoi_unit %*% t(core_unit)     # n_pixels × n_cores
  sim_mat[is.na(sim_mat)] <- 0
  sim_mat <- pmin(pmax(sim_mat, 0), 1)

  template  <- aoi_emb_rast[[1]]
  sim_rast  <- rep(template, n_cores)
  terra::values(sim_rast) <- sim_mat
  names(sim_rast) <- paste0("sim_", global_emb_df$profile_id)

  sim_rast
}
