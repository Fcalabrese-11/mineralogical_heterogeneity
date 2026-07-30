# =============================================================
# 09_raoq.R
# Compute Rao's Q quadratic entropy per pixel for all samples.
#
# Rao Q measures the mean pairwise Euclidean distance between NMF
# component-weight vectors in a sliding spatial window, capturing
# local compositional diversity in physical space.
#
# For each sample and each FoV the script:
#   1. Loads the post-NMF .RData (from script 05b)
#   2. Reads the corrected FoV file to recover spatial dimensions (nrow × ncol)
#   3. Applies the same normalisation + 10%-threshold + dominant-removal
#      as script 07, but per FoV (needed because RaoQ is a spatial metric)
#   4. Builds 8 component-weight subsets:
#        all / nobio / nores / nobio_nores
#        × each with and without dominant-component removal (nodom variants)
#   5. Computes the edge-corrected Rao Q for each subset using
#        edge_corrected_raoq() — properly clips the window at FoV borders
#        instead of padding, eliminating the edge artefact present in the
#        original spectralrao package implementation.
#   6. Saves RaoQ matrices and plots pheatmap heatmaps per FoV in a PDF.
#
# NOTE: RaoQ is computed by a double pixel loop and is O(n²) per window.
#       Expect several minutes per FoV. For large maps, consider running
#       per sample and saving intermediate .RData files.
#
# Prerequisite: config.R and functions/raman_functions.R loaded.
#               Post-NMF .RData must exist (from script 05b).
# Output:
#   paths$diversity_output / raoq_all_samples_<date>.RData  (combined table)
#   paths$plots            / <sample>_09_raoq_heatmaps_<date>.pdf
# =============================================================

options(max.print = 1000)
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(data.table)
library(pheatmap)
library(viridis)   # magma() colour palette

# ====================================================================
# Helper 1: edge-corrected Rao Q
#
# Computes Rao's quadratic entropy within a sliding window. Border pixels
# use a clipped (smaller) window rather than padding, avoiding the edge
# artefact in the original spectralrao package.
#
# @param mat_list  named list of nrow × ncol matrices, one per NMF component.
#                  mat_list[[k]][i, j] = normalised weight of component k at
#                  pixel (i, j).
# @param window    odd integer window size (default 3 = 3×3 neighbourhood)
# @return          matrix of Rao Q values, same nrow × ncol as input matrices
# ====================================================================
edge_corrected_raoq <- function(mat_list, window = 3L) {
  n_row   <- nrow(mat_list[[1L]])
  n_col   <- ncol(mat_list[[1L]])
  n_bands <- length(mat_list)

  # Stack into a 3D array: row × col × band
  array3d <- array(unlist(mat_list), dim = c(n_row, n_col, n_bands))

  pad     <- floor(window / 2L)
  rao_mat <- matrix(NA_real_, n_row, n_col)

  for (i in seq_len(n_row)) {
    for (j in seq_len(n_col)) {

      row_min <- max(1L, i - pad)
      row_max <- min(n_row, i + pad)
      col_min <- max(1L, j - pad)
      col_max <- min(n_col, j + pad)

      local_win <- array3d[row_min:row_max, col_min:col_max, , drop = FALSE]

      # Reshape to (n_pixels_in_window) × n_bands, column-major within window
      n_px_win   <- prod(dim(local_win)[1L:2L])
      local_flat <- matrix(aperm(local_win, c(2L, 1L, 3L)),
                           nrow = n_px_win, ncol = n_bands)

      # Drop all-NA pixels (e.g. masked pixels)
      local_flat <- local_flat[!apply(is.na(local_flat), 1L, all), , drop = FALSE]
      if (nrow(local_flat) < 2L) next

      # Rao Q = mean pairwise Euclidean distance (method = "euclidean" is the default
      # but stated explicitly to confirm correctness per analysis design)
      d <- dist(local_flat, method = "euclidean")
      rao_mat[i, j] <- sum(d) / choose(nrow(local_flat), 2L)
    }
  }
  rao_mat
}

# ====================================================================
# Helper 2: threshold + dominant-removal (mirrors script 07 steps 5–6)
#
# @param comp_mat  components × pixels matrix for one FoV.
#                  rownames must be set to the component labels.
# @param alpha     relative threshold (default 0.10 = same as script 07)
# @return list with:
#   $threshold   pixels × components data.table  (normalised + thresholded)
#   $minus       pixels × components data.table  (dominant component zeroed)
#   $comp_names  character vector of column names
# ====================================================================
threshold_and_minus <- function(comp_mat, alpha = 0.10) {

  X          <- as.data.table(t(comp_mat))   # pixels × components
  comp_names <- names(X)                      # = rownames(comp_mat) = comp labels

  # Step 1: row-normalise (each pixel sums to 1)
  X[, (comp_names) := lapply(.SD, function(v) v / rowSums(.SD)),
    .SDcols = comp_names]

  # Step 2: 10%-of-max threshold (zero out everything below alpha × row maximum)
  X[, row_max := do.call(pmax, .SD), .SDcols = comp_names]
  X[, (comp_names) := lapply(.SD, function(v) fifelse(v < alpha * row_max, 0, v)),
    .SDcols = comp_names]
  X[, row_max := NULL]

  # Step 3: re-normalise after thresholding
  X[, (comp_names) := lapply(.SD, function(v)
    ifelse(rowSums(.SD) > 0, v / rowSums(.SD), 0)), .SDcols = comp_names]

  # Step 4: dominant-removal copy
  X_minus <- copy(X)
  for (k in seq_len(nrow(X_minus))) {
    row_vals <- as.numeric(X_minus[k, ..comp_names])
    if (sum(row_vals) > 0) {
      row_vals[which.max(row_vals)] <- 0
      X_minus[k, (comp_names) := as.list(row_vals)]
    }
  }
  X_minus[, (comp_names) := {
    rs <- rowSums(.SD)
    lapply(seq_along(.SD), function(j) ifelse(rs > 0, .SD[[j]] / rs, 0))
  }, .SDcols = comp_names]

  list(threshold = X, minus = X_minus, comp_names = comp_names)
}

# ====================================================================
# Main loop
# ====================================================================
raoq_all <- list()   # accumulates per-sample per-FoV results

for (sample_name in names(samples)) {

  message("\n=== RaoQ: ", sample_name, " ===\n")
  cfg <- samples[[sample_name]]

  if (!file.exists(cfg$rdata_file)) {
    warning("Post-NMF .RData not found: ", cfg$rdata_file, " — skipping ", sample_name)
    next
  }
  load(cfg$rdata_file)   # loads: fov_h_FoV, fov.meta, fov_w.identity.cos, ...

  # Kept component IDs and their labels (from config.R + 05b output)
  kept_ids    <- as.numeric(fov_w.identity.cos$id)
  comp_labels <- cfg$comp_names[as.character(kept_ids)]
  bio_label   <- cfg$comp_names[as.character(cfg$biomass_comp)]
  res_label   <- cfg$comp_names[as.character(cfg$resin_comp)]
  ker_label   <- cfg$comp_names[as.character(cfg$kerogen_comp)]

  # ----------------------------------------------------------------
  # Read corrected FoV files to recover per-FoV spatial dimensions.
  # (Pixel count = n_row × n_col; row/col read from V1/V2 coordinates.)
  # ----------------------------------------------------------------
  fov_files <- list.files(path = cfg$fov_dir, pattern = cfg$fov_pattern,
                          full.names = TRUE)
  if (length(fov_files) == 0) {
    warning("No FoV files found for ", sample_name, " — skipping.")
    next
  }

  fov_dims <- list()
  for (fp in fov_files) {
    fname <- make.names(gsub("_corrected_SG_BL.*.txt", "", basename(fp)))
    raw_fov <- fread(fp, select = c(1L, 2L))
    fov_dims[[fname]] <- c(
      n_row = length(unique(raw_fov$V1)) - 1L,  # -1 for the wavelength header row
      n_col = length(unique(raw_fov$V2)) - 1L
    )
  }

  # ----------------------------------------------------------------
  # Open PDF for heatmaps
  # ----------------------------------------------------------------
  pdf_path <- file.path(paths$plots,
                        paste0(sample_name, "_09_raoq_heatmaps_", Sys.Date(), ".pdf"))
  pdf(pdf_path, width = 10, height = 8)

  sample_raoq_list <- list()

  for (fov_name in names(fov_h_FoV)) {
    message("  FoV: ", fov_name)

    if (!fov_name %in% names(fov_dims)) {
      warning("Spatial dimensions not found for FoV '", fov_name, "' — skipping.")
      next
    }

    n_row <- fov_dims[[fov_name]]["n_row"]
    n_col <- fov_dims[[fov_name]]["n_col"]
    message("    nrow=", n_row, "  ncol=", n_col, "  pixels=", n_row * n_col)

    # Extract kept components × pixels for this FoV
    fov_sub <- fov_h_FoV[[fov_name]][kept_ids, , drop = FALSE]
    rownames(fov_sub) <- comp_labels

    # Apply threshold and dominant-removal
    tm         <- threshold_and_minus(fov_sub)
    X_thres    <- tm$threshold   # pixels × components, thresholded
    X_minus    <- tm$minus       # pixels × components, dominant zeroed
    cn         <- tm$comp_names  # component labels (column names)

    # Helper: convert one column of a data.table to an nrow × ncol spatial matrix
    col_to_mat <- function(dt, col_name) {
      matrix(as.numeric(dt[[col_name]]), nrow = n_row, ncol = n_col, byrow = TRUE)
    }

    # Build 8 component subsets: {all, nobio, nores, nobio_nores} × {with/without nodom}
    excl_nobio       <- c(res_label, bio_label)
    excl_nores       <- res_label
    excl_nobio_nores <- c(bio_label, res_label,ker_label)

    subsets <- list(
      all               = list(dt = X_thres, cols = cn),
      nobio             = list(dt = X_thres, cols = setdiff(cn, excl_nobio)),
      nores             = list(dt = X_thres, cols = setdiff(cn, excl_nores)),
      nobio_nores       = list(dt = X_thres, cols = setdiff(cn, excl_nobio_nores)),
      all_nodom         = list(dt = X_minus, cols = cn),
      nobio_nodom       = list(dt = X_minus, cols = setdiff(cn, excl_nobio)),
      nores_nodom       = list(dt = X_minus, cols = setdiff(cn, excl_nores)),
      nobio_nores_nodom = list(dt = X_minus, cols = setdiff(cn, excl_nobio_nores))
    )

    fov_raoq <- list()

    for (subset_name in names(subsets)) {
      s    <- subsets[[subset_name]]
      cols <- s$cols

      if (length(cols) < 2L) {
        message("    Subset '", subset_name, "': fewer than 2 components — returning NA matrix.")
        fov_raoq[[subset_name]] <- matrix(NA_real_, n_row, n_col)
        next
      }

      # Build list of nrow × ncol matrices (one per component)
      mat_list        <- lapply(cols, col_to_mat, dt = s$dt)
      names(mat_list) <- cols

      message("    RaoQ: ", subset_name, "  (", length(cols), " components)")
      fov_raoq[[subset_name]] <- edge_corrected_raoq(mat_list, window = 3L)
    }

    # ----------------------------------------------------------------
    # Heatmap plots for this FoV
    # ----------------------------------------------------------------
    for (subset_name in names(fov_raoq)) {
      rao_mat <- fov_raoq[[subset_name]]
      if (all(is.na(rao_mat))) next

      pheatmap::pheatmap(
        rao_mat,
        scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
        cellwidth  = 4, cellheight = 4,
        legend = TRUE, show_rownames = FALSE, show_colnames = FALSE,
        border_color = FALSE,
        col = magma (50), rev = TRUE,
        #breaks = seq(0,1,0.02), # uncomment to set scale 0-1
        #color  = hcl.colors(50, "magma", rev = TRUE),
        main   = paste(sample_name, fov_name, "-", gsub("_", " ", subset_name))
      )
    }

    # ----------------------------------------------------------------
    # Flatten each RaoQ matrix to a per-pixel data.table row
    # Row-major unrolling (byrow = TRUE was used to build matrices,
    # so c(t(mat)) recovers the original pixel order).
    # ----------------------------------------------------------------
    safe_vec <- function(m) if (is.null(m)) rep(NA_real_, n_row * n_col) else c(t(m))

    fov_raoq_dt <- data.table(
      sample                = sample_name,
      fov                   = fov_name,
      px_id                 = seq_len(n_row * n_col),
      rao_all               = safe_vec(fov_raoq$all),
      rao_nobio             = safe_vec(fov_raoq$nobio),
      rao_nores             = safe_vec(fov_raoq$nores),
      rao_nobio_nores       = safe_vec(fov_raoq$nobio_nores),
      rao_all_nodom         = safe_vec(fov_raoq$all_nodom),
      rao_nobio_nodom       = safe_vec(fov_raoq$nobio_nodom),
      rao_nores_nodom       = safe_vec(fov_raoq$nores_nodom),
      rao_nobio_nores_nodom = safe_vec(fov_raoq$nobio_nores_nodom)
    )

    sample_raoq_list[[fov_name]] <- fov_raoq_dt
  }

  dev.off()
  message("Heatmap PDF saved: ", pdf_path)

  raoq_all[[sample_name]] <- rbindlist(sample_raoq_list)

  # Save per-sample .RData (allows resuming if a later sample fails)
  out_sample_rdata <- file.path(
    paths$diversity_output,
    paste0(sample_name, "_raoq_", Sys.Date(), ".RData")
  )
  raoq_sample <- raoq_all[[sample_name]]
  save(raoq_sample, file = out_sample_rdata)
  message("Per-sample RaoQ saved: ", out_sample_rdata)
}

# ====================================================================
# Combine across all samples
# ====================================================================
raoq_combined <- rbindlist(raoq_all)

out_rdata <- file.path(paths$diversity_output,
                       paste0("raoq_all_samples_", Sys.Date(), ".RData"))
save(raoq_combined, file = out_rdata)

message("\nRaoQ combined table saved to: ", out_rdata)
message("Columns: ", paste(names(raoq_combined), collapse = ", "))
message("=== 09_raoq.R complete ===")
