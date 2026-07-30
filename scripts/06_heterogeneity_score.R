# =============================================================
# 06_heterogeneity_score.R
# Compute per-pixel spectral heterogeneity (H) scores for all samples.
#
# For each FoV in each sample, the script:
#   - Reads the corrected FoV file
#   - Interpolates to common_x
#   - Calls corNeigh() (cosine dissimilarity to 8 neighbours)
#   - Saves an .rds list per sample containing the h-score rasters
# Then plots the H-score heatmaps for each FoV (pheatmap).
#
# Prerequisite: config.R and functions/raman_functions.R loaded.
# Output: one .rds file per sample in paths$h_scores.
# =============================================================

options(max.print = 100)
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(data.table)
library(raster)
library(lsa)
library(pheatmap)

# ====================================================================
# Main loop over samples
# ====================================================================
for (sample_name in names(samples)) {

  message("\n=== H-score: ", sample_name, " ===\n")
  cfg   <- samples[[sample_name]]

  fovs <- list.files(path = cfg$fov_dir, pattern = cfg$fov_pattern, full.names = TRUE)
  if (length(fovs) == 0) {
    warning("No FoV files found for ", sample_name, " — skipping.")
    next
  }

  h_scores <- list()

  for (i in fovs) {
    message(i)

    # Bug fix: unified file-name pattern (was "SGL_BL" typo for VR and MN)
    tmp_name   <- make.names(gsub("_corrected_SG_BL.*.txt", "", basename(i)))
    raw_sample <- fread(i)

    my_rows   <- length(unique(raw_sample$V1)) - 1
    my_cols   <- length(unique(raw_sample$V2)) - 1
    my_pixels <- my_rows * my_cols
    message("rows: ", my_rows, "  cols: ", my_cols, "  pixels: ", my_pixels)

    # Interpolation (common_x defined globally in config.R)
    mysample <- raw_sample[, -c(1, 2)]
    mysample <- as.data.table(t(apply(mysample[-1, ], 1, interpolation,
                                       x       = as.vector(as.matrix(mysample[1, ])),
                                       common_x = common_x)))
    mysample[, id := 1:nrow(mysample)]
    setcolorder(mysample, c("id", setdiff(names(mysample), "id")))

    # Heterogeneity score
    h_score     <- corNeigh(r.sample = copy(mysample), n.neigh = 8,
                             n.row = my_rows, n.col = my_cols, by.row = TRUE)
    h_score_cum <- sum(values(h_score) / my_pixels)
    message("Cumulative H-score: ", round(h_score_cum, 4))

    h_scores[[tmp_name]] <- list(
      FoV         = tmp_name,
      mysample    = mysample,
      h_score     = h_score,
      h_score_cum = h_score_cum,
      my_rows     = my_rows,
      my_cols     = my_cols,
      my_pixels   = my_pixels
    )
  }

  # Save h-scores for this sample
  if (!dir.exists(paths$h_scores)) dir.create(paths$h_scores, recursive = TRUE)
  
  out_rds <- file.path(paths$h_scores,
                       paste0(sample_name, "_FoVs.hscore.100-3000range.corrected_SG_BL_",
                              Sys.Date(), ".rds"))
  saveRDS(h_scores, file = out_rds)
  message("H-scores saved: ", out_rds)

  # ------------------------------------------------------------------
  # Plot heatmaps of H-score per FoV
  # ------------------------------------------------------------------
  for (fov_name in names(h_scores)) {
    fov_data <- h_scores[[fov_name]]
    my_rows  <- fov_data$my_rows
    my_cols  <- fov_data$my_cols

    pheatmap::pheatmap(
      matrix(values(fov_data$h_score), nrow = my_rows, ncol = my_cols, byrow = TRUE),
      scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
      cellwidth = 6, cellheight = 6,
      legend = TRUE, show_rownames = FALSE, show_colnames = FALSE,
      border_color = FALSE,
      color  = hcl.colors(10, "Blue-Yellow", rev = TRUE),
      breaks = seq(0, 0.25, 0.025),
      main   = paste(sample_name, fov_name)
    )
  }

  # ------------------------------------------------------------------
  # H-score density plots (requires post-NMF biomass_comp_ids in environment)
  # Loaded from the post-NMF .RData if it exists
  # ------------------------------------------------------------------
  if (file.exists(cfg$rdata_file)) {

    tmp_env <- new.env()
    load(cfg$rdata_file, envir = tmp_env)

    if (exists("biomass_comp_ids", envir = tmp_env)) {

      biomass_ids         <- tmp_env$biomass_comp_ids
      hscore_all_px       <- do.call(c, sapply(h_scores, function(x) values(x$h_score)))
      hscore_bio_px       <- hscore_all_px[biomass_ids]

      sorted_HS     <- sort(hscore_all_px,  decreasing = TRUE)
      sorted_HS_bio <- sort(hscore_bio_px, decreasing = TRUE)

      par(mfrow = c(2, 2), mar = c(4, 4, 4, 2))
      plot(sorted_HS,     type = "l", ylim = c(0, 1), main = paste(sample_name, "— all px"),
           xlab = "pixels", ylab = "H score")
      plot(sorted_HS_bio, type = "l", ylim = c(0, 1), main = paste(sample_name, "— biomass px"),
           xlab = "pixels", ylab = "H score")
      plot(density(sorted_HS),     col = "cyan",  xlim = c(0, 1),
           main = paste(sample_name, "— density all px"))
      plot(density(sorted_HS_bio), col = "blue",  xlim = c(0, 1),
           main = paste(sample_name, "— density biomass px"))
    }
  }

  message("=== Done: ", sample_name, " ===\n")
}
