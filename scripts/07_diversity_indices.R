# =============================================================
# 07_diversity_indices.R
# Compute diversity indices (Gini, SD, Shannon, etc.) for all samples.
#
# For each sample this script:
#   1. Loads the post-NMF .RData (from script 05)
#   2. Builds the components × pixels matrix and names it from config
#   3. Computes original (non-normalised) Gini
#   4. Normalises by row, applies 10%-of-max threshold, re-normalises
#   5. Removes the dominant component per pixel
#   6. Assembles a per-pixel table of the retained diversity metrics
#      (Shannon entropy and Gini; the heterogeneity score and Rao's Q
#       are produced by scripts 06 and 09)
#   7. Adds sample labels
#   8. Combines across all samples and saves combined .RData
#      → used by 08_plot_diversity.R
#
# Also computes the Gini drop-loop (effect of removing each component
# on overall Gini).
#
# Prerequisite: config.R and functions/raman_functions.R loaded.
# =============================================================

options(max.print = 1000)
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(data.table)
library(vegan)   # diversity()
library(ineq)    # ineq() for Gini
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# ====================================================================
# Helper: compute all diversity metrics for one row (proportions)
# ====================================================================
compute_metrics_all <- function(x) {
  if (sum(x) == 0)
    return(list(shannon = 0, simpson = 0, gini_simpson = 0, gini = 0,
                sd = 0, hill_q1 = 0, hill_q2 = 0, pielou = 0,
                richness = 0, hill_infinite = 0))

  # Only components that are actually present (non-zero after thresholding).
  # Gini and SD must not include absent components — zeros inflate both.
  x_pos <- x[x > 0]
  S     <- length(x_pos)

  # vegan::diversity() applies the 0*log(0)=0 convention, so zeros are
  # correctly excluded from Shannon and Simpson without any epsilon trick.
  H         <- vegan::diversity(x, index = "shannon")   # Shannon entropy (nats)
  gini_simp <- vegan::diversity(x, index = "simpson")   # 1 - sum(p^2)
  sum_p2    <- 1 - gini_simp                            # Simpson concentration

  list(
    shannon       = H,
    simpson       = sum_p2,
    gini_simpson  = gini_simp,
    # Gini and SD computed on present components only
    gini          = if (S > 1) ineq::ineq(x_pos, type = "Gini") else 0,
    sd            = if (S > 1) sd(x_pos) else 0,
    hill_q1       = exp(H),                             # Hill q=1
    hill_q2       = if (sum_p2 > 0) 1 / sum_p2 else 0, # Hill q=2
    pielou        = if (S > 1) H / log(S) else 0,
    richness      = S,
    hill_infinite = 1 / max(x)                          # reciprocal Berger-Parker
  )
}

extract_metric <- function(res, name) sapply(res, `[[`, name)

# ====================================================================
# Helper: classify each pixel by its dominant mineral identity
# (defined here so it is available both inside the loop and in the
#  combining section below)
# ====================================================================
classify_identity <- function(dominant, biomass_pattern = "Biomass",
                              kerogen_pattern = "Kerogen",
                              resin_pattern   = "Resin") {
  duo <- ifelse(grepl(biomass_pattern, dominant, ignore.case = TRUE), #grepl(kerogen_pattern, dominant, ignore.case = TRUE) |
                "biomass", "other")
  duo <- ifelse(grepl(resin_pattern, dominant, ignore.case = TRUE), "resin", duo)
  
  bioker <- ifelse(grepl(kerogen_pattern, dominant, ignore.case = TRUE) |
                     grepl(biomass_pattern, dominant, ignore.case = TRUE),
                   "biomass", "mineral")
  bioker <- ifelse(grepl(resin_pattern, dominant, ignore.case = TRUE), "resin", bioker)
  
  list(duo_identity = duo, bioker_identity = bioker)
}
# ====================================================================
# Helper: PCA on a pixels × components matrix, coloured by dominant
# ====================================================================
# mat       : numeric matrix, pixels × components (all-zero columns are dropped)
# dom_vec   : character vector length nrow(mat), dominant component per pixel
# title_str : plot title
# out_file  : path to output PDF
pca_plot <- function(mat, dom_vec, title_str, out_file) {
  # Drop constant columns (zero-variance after thresholding)
  keep <- apply(mat, 2, var) > 0
  if (sum(keep) < 2) {
    message("  PCA skipped (fewer than 2 variable columns): ", title_str)
    return(invisible(NULL))
  }
  pca_res <- prcomp(mat[, keep, drop = FALSE], center = TRUE, scale. = TRUE)
  scores  <- as.data.frame(pca_res$x[, 1:2])
  colnames(scores) <- c("PC1", "PC2")
  scores$dominant  <- dom_vec

  var_exp <- summary(pca_res)$importance[2, 1:2] * 100   # % variance explained

  p <- ggplot(scores, aes(x = PC1, y = PC2, colour = dominant)) +
    geom_point(alpha = 0.35, size = 0.7) +
    labs(
      title  = title_str,
      x      = sprintf("PC1 (%.1f %%)", var_exp[1]),
      y      = sprintf("PC2 (%.1f %%)", var_exp[2]),
      colour = "Dominant\ncomponent"
    ) +
    theme_bw(base_size = 11) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme(legend.text  = element_text(size = 8),
          legend.title = element_text(size = 9))

  ggsave(out_file, plot = p, width = 7, height = 5, dpi = 150)
  message("  PCA plot saved: ", out_file)
  invisible(pca_res)
}

# ====================================================================
# Containers for combined results across samples
# ====================================================================
gini_sd_hs_all    <- list()   # original (non-normalised) gini/hscore per sample
metric_all_list   <- list()   # post-threshold diversity metrics (Shannon, Gini) per sample
gini_drop_loop_all <- list()  # gini_without named list per sample (for 08 drop-loop plots)
gini_ref_all      <- list()   # gini_above_cosine_all_comps per sample

# ====================================================================
# Main loop
# ====================================================================
for (sample_name in names(samples)) {

  message("\n=== Diversity indices: ", sample_name, " ===\n")
  cfg <- samples[[sample_name]]

  if (!file.exists(cfg$rdata_file)) {
    warning("Post-NMF .RData not found: ", cfg$rdata_file, " — skipping ", sample_name)
    next
  }
  load(cfg$rdata_file)   # loads fov_h, comp_above_cosine, biomass_comp_ids, resin_ids, etc.

  # ----------------------------------------------------------------
  # 1. Build X (pixels × components) from comp_above_cosine
  #    comp_above_cosine is standardised as components × pixels (script 05)
  # ----------------------------------------------------------------
  kept_comp_ids <- setdiff(seq_len(cfg$n_components), cfg$exclude_comps)
  col_labels    <- cfg$comp_names[as.character(kept_comp_ids)]

  X <- as.data.table(t(comp_above_cosine))   # pixels × components
  setnames(X, old = names(X), new = col_labels)
  comp_cols <- col_labels

  message("X dim (pixels × components): ", paste(dim(X), collapse = " x "))

  # ----------------------------------------------------------------
  # 2. Original (non-normalised) Gini and SD
  #    all_X, noresin_X, nobio_noresin_X  (components × pixels)
  # ----------------------------------------------------------------
  all_X          <- comp_above_cosine
  noresin_X      <- fov_h[-c(cfg$exclude_comps, cfg$resin_comp), ]
  nobio_noresin_X <- fov_h[-c(cfg$exclude_comps, cfg$resin_comp,
                                cfg$biomass_comp, cfg$kerogen_comp), ]

  message("dim all_X: ", paste(dim(all_X), collapse = " x "))
  message("dim noresin_X: ", paste(dim(noresin_X), collapse = " x "))
  message("dim nobio_noresin_X: ", paste(dim(nobio_noresin_X), collapse = " x "))

  # Helper that restricts to non-zero component weights per pixel so that
  # absent components do not inflate Gini.
  gini_nz <- function(v) { v <- v[v > 0]; if (length(v) > 1) ineq::ineq(v, type = "Gini") else 0 }

  gini_above_cosine_all_comps     <- apply(all_X,           2, gini_nz)
  gini_above_cosine_noresin       <- apply(noresin_X,       2, gini_nz)
  gini_above_cosine_nobio_noresin <- apply(nobio_noresin_X, 2, gini_nz)

  # ----------------------------------------------------------------
  # 3. Gini drop-loop: effect of removing each component on Gini.
  #    Uses gini_nz (non-zero only) defined in section 2 above.
  # ----------------------------------------------------------------
  gini_without <- list()
  for (i in seq_len(nrow(all_X))) {
    gini_without[[i]] <- apply(all_X[-i, ], 2, gini_nz)
  }
  names(gini_without) <- paste0("without_", col_labels)

  # ----------------------------------------------------------------
  # 4. Load h-scores
  # ----------------------------------------------------------------
  het_score_all_px <- NULL
  if (file.exists(cfg$hscore_file)) {
    h_scores <- readRDS(cfg$hscore_file)
    het_score_all_px <- do.call(c, sapply(h_scores, function(x) values(x$h_score)))
    message("H-score pixels loaded: ", length(het_score_all_px))
  } else {
    warning("H-score file not found: ", cfg$hscore_file)
    het_score_all_px <- rep(NA_real_, ncol(all_X))
  }

  # ----------------------------------------------------------------
  # 5. Normalise X per pixel (row), threshold, re-normalise
  # ----------------------------------------------------------------
  X_norm <- copy(X)
  X_norm[, (comp_cols) := lapply(.SD, function(x) x / rowSums(.SD)), .SDcols = comp_cols]

  alpha <- 0.10   # keep components contributing > 10% of dominant weight
  X_norm_threshold <- as.data.table(X_norm)
  X_norm_threshold[, row_max := do.call(pmax, .SD), .SDcols = comp_cols]
  X_norm_threshold[, (comp_cols) := lapply(.SD, function(x)
    fifelse(x < alpha * row_max, 0, x)), .SDcols = comp_cols]
  X_norm_threshold[, row_max := NULL]

  X_norm_threshold[, (comp_cols) := lapply(.SD, function(x)
    ifelse(rowSums(.SD) > 0, x / rowSums(.SD), 0)), .SDcols = comp_cols]

  # ----------------------------------------------------------------
  # 6. Remove dominant component per pixel
  # ----------------------------------------------------------------
  X_norm_threshold[, dominant := comp_cols[max.col(.SD)], .SDcols = comp_cols]
  dominant_comp <- X_norm_threshold$dominant

  X_norm_threshold_minus <- copy(X_norm_threshold)
  for (i in seq_len(nrow(X_norm_threshold_minus))) {
    row_vals <- as.numeric(X_norm_threshold_minus[i, ..comp_cols])
    if (sum(row_vals) > 0) {
      row_vals[which.max(row_vals)] <- 0
      X_norm_threshold_minus[i, (comp_cols) := as.list(row_vals)]
    }
  }
  X_norm_threshold_minus[, (comp_cols) := {
    rs <- rowSums(.SD)
    lapply(seq_along(.SD), function(j) ifelse(rs > 0, .SD[[j]] / rs, 0))
  }, .SDcols = comp_cols]

  # ----------------------------------------------------------------
  # 6b. gini_original_minus: non-thresholded Gini with per-pixel
  #     dominant component zeroed (mirrors gini_thres_minus but on
  #     raw NMF H weights in all_X rather than the thresholded matrix)
  # ----------------------------------------------------------------
  dom_row_idx <- match(dominant_comp, col_labels)
  gini_original_minus <- sapply(seq_len(ncol(all_X)), function(j) {
    v       <- all_X[, j]
    dom_idx <- dom_row_idx[j]
    if (!is.na(dom_idx)) v[dom_idx] <- 0
    gini_nz(v)
  })

  # ----------------------------------------------------------------
  # 7. Diversity metrics on thresholded matrix
  # ----------------------------------------------------------------
  # Names of biomass and resin columns (for sub-setting)
  bio_col    <- cfg$comp_names[as.character(cfg$biomass_comp)]
  resin_col  <- cfg$comp_names[as.character(cfg$resin_comp)]
  kero_col   <- cfg$comp_names[as.character(cfg$kerogen_comp)]

  remove_for_noresin      <- c(resin_col, "dominant")
  remove_for_nobio_noresin <- c(bio_col, resin_col, kero_col, "dominant")

  X_norm_thres_noresin     <- X_norm_threshold[, !remove_for_noresin,  with = FALSE]
  X_norm_thres_nobio_noresin <- X_norm_threshold[, !remove_for_nobio_noresin, with = FALSE]

  results_all         <- apply(as.matrix(X_norm_threshold[, ..comp_cols]), 1, compute_metrics_all)
  results_minus       <- apply(as.matrix(X_norm_threshold_minus[, ..comp_cols]), 1, compute_metrics_all)
  results_noresin     <- apply(as.matrix(X_norm_thres_noresin),      1, compute_metrics_all)
  results_nobio_noresin <- apply(as.matrix(X_norm_thres_nobio_noresin), 1, compute_metrics_all)

  # ----------------------------------------------------------------
  # 8. Assemble per-pixel metrics data table for this sample
  # ----------------------------------------------------------------
  n <- nrow(X_norm_threshold)

  # Metrics retained in the analysis: Shannon entropy and Gini
  # (plus the heterogeneity score and Rao's Q, handled elsewhere).
  metric_all <- as.data.table(lapply(
    c("shannon","gini"),
    function(m) extract_metric(results_all, m)
  ))
  setnames(metric_all, c("shannon","gini"))
  metric_all[, px_id   := .I]
  metric_all[, sample  := sample_name]
  metric_all[, identity := dominant_comp]

  metric_all_list[[sample_name]] <- metric_all

  # ----------------------------------------------------------------
  # 9. Original Gini/SD/H-score combined table for this sample.
  #   Bug B fix: constructed here (after step 8) so dominant_comp and
  #              results_noresin/results_nobio_noresin are all available.
  #   Bug A fix: add threshold noresin/nobio_noresin gini+sd from those results.
  #   Bug B fix: add identity/duo_identity/bioker_identity columns.
  # ----------------------------------------------------------------
  idents_px <- classify_identity(dominant_comp)
  n_px <- ncol(all_X)
  gini_sd_hs <- data.table(
    sample                   = sample_name,
    px_id                    = seq_len(n_px),
    identity                 = dominant_comp,
    duo_identity             = idents_px$duo_identity,
    bioker_identity          = idents_px$bioker_identity,
    gini                     = gini_above_cosine_all_comps,
    gini_nores               = gini_above_cosine_noresin,
    gini_nobio_nores         = gini_above_cosine_nobio_noresin,
    gini_original_minus      = gini_original_minus,
    het_score                = het_score_all_px,
    gini_thres_noresin            = extract_metric(results_noresin,        "gini"),
    gini_thres_nobio_noresin      = extract_metric(results_nobio_noresin,  "gini"),
    shannon_thres_noresin         = extract_metric(results_noresin,        "shannon"),
    shannon_thres_nobio_noresin   = extract_metric(results_nobio_noresin,  "shannon"),
    shannon_thres_minus           = extract_metric(results_minus,          "shannon"),
    gini_thres_minus              = extract_metric(results_minus,          "gini")
  )
  gini_sd_hs_all[[sample_name]] <- gini_sd_hs

  # ----------------------------------------------------------------
  # 10. Save per-sample diversity tables as text files
  # ----------------------------------------------------------------

  out_dir <- file.path(paths$diversity_output, sample_name)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  write.table(metric_all,
              file  = file.path(out_dir, paste0(sample_name, "_diversity_all_threshold.txt")),
              quote = FALSE, sep = "\t", row.names = FALSE)

  write.table(gini_sd_hs,
              file  = file.path(out_dir, paste0(sample_name, "_gini_sd_hetscore_original.txt")),
              quote = FALSE, sep = "\t", row.names = FALSE)

  # ----------------------------------------------------------------
  # 11. PCA on normalised-thresholded matrices (pixels × components)
  #     Three versions: all components / no-resin / no-bio-no-resin.
  #     Pixels are coloured by their dominant component.
  # ----------------------------------------------------------------
  message("  Running PCA plots for ", sample_name, " ...")

  # (a) All components
  pca_plot(
    mat       = as.matrix(X_norm_threshold[, ..comp_cols]),
    dom_vec   = dominant_comp,
    title_str = paste0(sample_name, " — PCA all components (thresholded)"),
    out_file  = file.path(out_dir, paste0(sample_name, "_pca_all.pdf"))
  )

  # (b) No-resin: columns already selected in X_norm_thres_noresin
  pca_plot(
    mat       = as.matrix(X_norm_thres_noresin),
    dom_vec   = dominant_comp,
    title_str = paste0(sample_name, " — PCA no-resin (thresholded)"),
    out_file  = file.path(out_dir, paste0(sample_name, "_pca_noresin.pdf"))
  )

  # (c) No-bio, no-resin: columns already selected in X_norm_thres_nobio_noresin
  pca_plot(
    mat       = as.matrix(X_norm_thres_nobio_noresin),
    dom_vec   = dominant_comp,
    title_str = paste0(sample_name, " — PCA no-bio no-resin (thresholded)"),
    out_file  = file.path(out_dir, paste0(sample_name, "_pca_nobio_noresin.pdf"))
  )

  # (d) Dominant-removed: X_norm_threshold_minus (dominant set to 0 per pixel,
  #     remainder re-normalised). dom_vec is still the original dominant so the
  #     colour legend matches the other three plots.
  pca_plot(
    mat       = as.matrix(X_norm_threshold_minus[, ..comp_cols]),
    dom_vec   = dominant_comp,
    title_str = paste0(sample_name, " — PCA dominant-removed (thresholded minus)"),
    out_file  = file.path(out_dir, paste0(sample_name, "_pca_minus.pdf"))
  )

  # Store drop-loop results for 08_plot_diversity.R
  gini_drop_loop_all[[sample_name]] <- gini_without
  gini_ref_all[[sample_name]]       <- gini_above_cosine_all_comps

  message("Diversity tables saved to: ", out_dir)
  message("=== Done: ", sample_name, " ===\n")
}

# ====================================================================
# Combine across all samples for cross-sample plotting (08_plot_diversity.R)
# ====================================================================
gini_sd_hs_combined <- rbindlist(gini_sd_hs_all)
metric_combined     <- rbindlist(metric_all_list)

# ---- Build combined long-format tables for gini and SD ----
# (these are the objects expected by 08_plot_diversity.R)

# Original (non-normalised) Gini comparison: nores vs. nobio_nores
# Bug C fix: include identity/duo_identity/bioker_identity in id.vars (from gini_sd_hs_combined)
#            so that downstream facet plots can use these columns directly.
#            The four broken classify_identity(variable) calls are removed.
het_stats_melt_gini_original <- melt(
  gini_sd_hs_combined[, .(sample_original = sample, px_id, identity, duo_identity, bioker_identity,
                           gini_nores_original       = 1 - gini_nores,
                           gini_nobio_nores_original = 1 - gini_nobio_nores,
                           het_score)],
  id.vars      = c("sample_original", "px_id", "identity", "duo_identity", "bioker_identity",
                   "het_score"),
  measure.vars = c("gini_nores_original", "gini_nobio_nores_original"),
  variable.name = "variable",
  value.name    = "value"
)

# Load per-sample threshold metrics (Shannon used in 08_plot_diversity.R p10)
threshold_list <- lapply(names(samples), function(sn) {
  f <- file.path(paths$diversity_output, sn,
                 paste0(sn, "_diversity_all_threshold.txt"))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f)
  dt[, sample_original := sn]
  dt
})
threshold_combined <- rbindlist(Filter(Negate(is.null), threshold_list))
idents_tc <- classify_identity(threshold_combined$identity)
threshold_combined[, duo_identity    := idents_tc$duo_identity]
threshold_combined[, bioker_identity := idents_tc$bioker_identity]

#   after-threshold Gini melt tables use gini_sd_hs_combined columns which were
#   computed from results_noresin / results_nobio_noresin (not from metric_all for both)

het_stats_melt_gini_threshold <- melt(
  gini_sd_hs_combined[, .(sample_original = sample, px_id, identity, duo_identity, bioker_identity,
                           gini_thres_noresin       = 1 - gini_thres_noresin,
                           gini_thres_nobio_noresin = 1 - gini_thres_nobio_noresin)],
  id.vars      = c("sample_original", "px_id", "identity", "duo_identity", "bioker_identity"),
  measure.vars = c("gini_thres_noresin", "gini_thres_nobio_noresin"),
  variable.name = "variable",
  value.name    = "value"
)

# Gini minus (threshold): dominant-removed vs original, both as 1-gini
het_stats_melt_gini_minus <- melt(
  gini_sd_hs_combined[, .(sample_original = sample, px_id, identity, duo_identity, bioker_identity,
                           gini_thres_noresin = 1 - gini_thres_noresin,
                           gini_thres_minus   = 1 - gini_thres_minus)],
  id.vars       = c("sample_original", "px_id", "identity", "duo_identity", "bioker_identity"),
  measure.vars  = c("gini_thres_noresin", "gini_thres_minus"),
  variable.name = "variable",
  value.name    = "value"
)

# Gini minus (original, non-thresholded): dominant-removed vs original, both as 1-gini
het_stats_melt_gini_original_minus <- melt(
  gini_sd_hs_combined[, .(sample_original = sample, px_id, identity, duo_identity, bioker_identity,
                           gini_original      = 1 - gini,
                           gini_original_minus = 1 - gini_original_minus)],
  id.vars       = c("sample_original", "px_id", "identity", "duo_identity", "bioker_identity"),
  measure.vars  = c("gini_original", "gini_original_minus"),
  variable.name = "variable",
  value.name    = "value"
)

# ====================================================================
# Save combined objects for 08_plot_diversity.R
# ====================================================================
diversity_rdata <- file.path(paths$diversity_output,
                              paste0("diversity_all_samples_", Sys.Date(), ".RData"))

save(gini_sd_hs_combined,
     metric_combined,
     het_stats_melt_gini_original,
     het_stats_melt_gini_threshold,
     het_stats_melt_gini_minus,
     het_stats_melt_gini_original_minus,
     threshold_combined,
     gini_drop_loop_all,
     gini_ref_all,
     file = diversity_rdata)

message("\nCombined diversity objects saved to: ", diversity_rdata)
