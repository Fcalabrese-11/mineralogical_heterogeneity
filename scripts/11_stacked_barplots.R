# =============================================================
# 11_stacked_barplots.R
#
# Stacked barplots and categorical identity maps for all samples.
#
# Produces three types of plots per sample, saved in one PDF:
#   A. Mean NMF-weight proportion per dominant component
#      — answers "when a pixel is dominated by mineral X,
#        what weights do the other minerals carry on average?"
#   B. Categorical spatial maps per FoV
#      — each pixel coloured by its top-1 dominant mineral
#   C. Co-occurrence stacked barplots (top-1 vs top-2)
#      — for each dominant mineral, what is the distribution
#        of second-ranked minerals across its pixels?
#   D. Colour legend page
#
# Prerequisite: scripts 05b and 06 must have run so that
#   cfg$rdata_file contains fov_h_long_top1, fov_h_long_top2,
#   fov_h_long_top_third, comp_above_cosine, fov.meta,
#   fov_w.identity.cos, fov_h_long_top3_occ.
#   cfg$hscore_file contains per-FoV my_rows / my_cols.
# =============================================================

options(stringsAsFactors = FALSE)
options(scipen = 999)

library(data.table)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# ---- Load pipeline configuration if not already in environment ----
if (!exists("paths")) source("config.R")

# ====================================================================
# Shared mineral colour palette — single source of truth in config.R
# (`color_vector_minerals`). Add colours for new minerals there.
# Unlisted minerals fall back to grey in the plotting code below.
# ====================================================================
if (!exists("color_vector_minerals")) source("config.R")

# ====================================================================
# Main loop over samples
# ====================================================================
for (sample_name in names(samples)) {

  message("\n=== 11_stacked_barplots: ", sample_name, " ===\n")
  cfg <- samples[[sample_name]]

  # Skip if config.R component fields are still NULL
  if (is.null(cfg$comp_names)) {
    message("  comp_names is NULL — skipping ", sample_name,
            ".\n  Fill config.R after running 05a.")
    next
  }

  # ------------------------------------------------------------------
  # 1. Load post-NMF objects (from 05b_post_nmf_assign.R)
  # ------------------------------------------------------------------
  if (!file.exists(cfg$rdata_file)) {
    warning("  RData not found: ", cfg$rdata_file, " — skipping ", sample_name)
    next
  }
  load(cfg$rdata_file)
  # Now available: fov_h_long_top1, fov_h_long_top2, fov_h_long_top_third,
  #   fov_h_long_top3, fov_h_long_top3_occ, comp_above_cosine,
  #   fov_w.identity.cos, fov.meta, fov_h_FoV
  message("  Loaded: ", cfg$rdata_file)

  # ------------------------------------------------------------------
  # 2. Load h-scores for per-FoV grid dimensions (my_rows, my_cols)
  # ------------------------------------------------------------------
  if (!file.exists(cfg$hscore_file)) {
    warning("  h-score file not found: ", cfg$hscore_file, " — skipping ", sample_name)
    next
  }
  h_scores <- readRDS(cfg$hscore_file)
  message("  Loaded h-scores: ", cfg$hscore_file)

  # ------------------------------------------------------------------
  # 3. Build per-sample colour vector
  #    Map each kept component → its mineral label → shared colour.
  #    Components not found in color_vector_minerals get grey (#AAAAAA).
  # ------------------------------------------------------------------
  all_comp_ids <- fov_w.identity.cos$id   # character, e.g. "1","2",...

  # mineral_for_comp: named by comp.id, value = human-readable mineral label
  mineral_for_comp <- sapply(all_comp_ids, function(id) {
    nm <- cfg$comp_names[id]
    if (is.na(nm) || is.null(nm)) paste0("comp", id) else nm
  })

  # col_mineral: named by comp.id
  col_mineral <- setNames(
    ifelse(mineral_for_comp %in% names(color_vector_minerals),
           color_vector_minerals[mineral_for_comp],
           "#AAAAAA"),
    all_comp_ids
  )

  # ------------------------------------------------------------------
  # 4. Open per-sample PDF
  # ------------------------------------------------------------------
  if (!dir.exists(paths$plots)) dir.create(paths$plots, recursive = TRUE)

  pdf_path <- file.path(
    paths$plots,
    paste0(sample_name, "_stacked_barplots_",
           format(Sys.Date(), "%Y-%m-%d"), ".pdf")
  )
  pdf(file = pdf_path, width = 10, height = 7)
  message("  Writing PDF: ", pdf_path)

  # ==================================================================
  # PLOT A — Mean NMF-weight proportion per dominant component
  # ==================================================================
  # Build df_components: pixels × kept-components
  # comp_above_cosine has orientation: components × pixels (from 05b)
  df_components <- as.data.table(t(comp_above_cosine))

  # Column names = mineral labels in the same order as rows of comp_above_cosine
  kept_ids_sorted <- as.character(sort(as.numeric(all_comp_ids)))
  col_labels       <- sapply(kept_ids_sorted, function(id) {
    nm <- cfg$comp_names[id]
    if (is.na(nm) || is.null(nm)) paste0("comp", id) else nm
  })
  setnames(df_components, col_labels)

  # Add per-pixel metadata from fov_h_long_top1
  df_components[, comp_id := fov_h_long_top1$comp.id]
  df_components[, mineral  := fov_h_long_top1$mineral]

  # Normalize each pixel's component weights to proportions (row sums → 1)
  weight_cols <- col_labels
  df_components[, total_weight := rowSums(.SD), .SDcols = weight_cols]
  df_components[, (weight_cols) := lapply(.SD, function(x) x / total_weight),
                .SDcols = weight_cols]

  # Average proportions across pixels grouped by dominant component
  df_avg <- df_components[, lapply(.SD, mean), by = comp_id,
                          .SDcols = weight_cols]

  # Add human-readable x-axis label (mineral name of the dominant component)
  df_avg[, dominant_mineral := sapply(comp_id, function(id) {
    nm <- cfg$comp_names[id]
    if (is.na(nm) || is.null(nm)) paste0("comp", id) else nm
  })]

  # Melt to long format for ggplot
  dt_long <- melt(df_avg,
                  id.vars      = c("comp_id", "dominant_mineral"),
                  variable.name = "component",
                  value.name   = "proportion")

  # Order x-axis by dominant mineral name
  dt_long[, dominant_mineral := factor(dominant_mineral,
                                       levels = unique(dominant_mineral))]

  # Colour lookup by component label (not comp.id)
  col_by_label <- setNames(
    ifelse(col_labels %in% names(color_vector_minerals),
           color_vector_minerals[col_labels],
           "#AAAAAA"),
    col_labels
  )

  p_prop <- ggplot(dt_long,
                   aes(x = dominant_mineral, y = proportion, fill = component)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = col_by_label) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste(sample_name, "— mean weight proportion by dominant component")) +
    ylab("Mean proportion of weight") +
    xlab("Dominant component (top 1 mineral)")

  print(p_prop)

  # ==================================================================
  # PLOT B — Categorical spatial maps per FoV
  #   Each pixel coloured by its dominant (top-1) mineral.
  # ==================================================================
  # Compute per-FoV pixel boundaries from fov.meta
  group_lengths <- sapply(fov.meta, function(x) length(x[["id"]]))
  end_idx        <- cumsum(group_lengths)
  start_idx      <- c(1L, head(end_idx, -1L) + 1L)

  for (fi in seq_along(names(fov_h_FoV))) {

    fov_name <- names(fov_h_FoV)[fi]

    if (!fov_name %in% names(h_scores)) {
      warning("    FoV '", fov_name, "' not in h_scores — skipping map.")
      next
    }
    rows <- h_scores[[fov_name]]$my_rows
    cols <- h_scores[[fov_name]]$my_cols

    # Slice the top-1 assignments for this FoV's pixels
    fov_top1 <- fov_h_long_top1[start_idx[fi]:end_idx[fi], ]

    # Encode comp.id as integer factor with ALL levels for consistent colours
    comps_factor <- factor(fov_top1$comp.id, levels = all_comp_ids)
    comp_matrix  <- matrix(as.integer(comps_factor),
                           nrow  = rows, ncol = cols, byrow = TRUE)

    # pheatmap colour vector must be length = number of factor levels
    pheatmap_cols <- col_mineral[levels(comps_factor)]

    pheatmap::pheatmap(
      comp_matrix,
      color         = pheatmap_cols,
      breaks        = seq(0.5, length(all_comp_ids) + 0.5, by = 1),
      cluster_rows  = FALSE,
      cluster_cols  = FALSE,
      scale         = "none",
      cellwidth     = 5,
      cellheight    = 5,
      show_rownames = FALSE,
      show_colnames = FALSE,
      border_color  = NA,
      legend        = FALSE,
      main          = paste0(sample_name, " — dominant mineral  |  ", fov_name)
    )
  }

  # Do the same for the top2 identity
  group_lengths <- sapply(fov.meta, function(x) length(x[["id"]]))
  end_idx        <- cumsum(group_lengths)
  start_idx      <- c(1L, head(end_idx, -1L) + 1L)
  
  for (fi in seq_along(names(fov_h_FoV))) {
    
    fov_name <- names(fov_h_FoV)[fi]
    
    if (!fov_name %in% names(h_scores)) {
      warning("    FoV '", fov_name, "' not in h_scores — skipping map.")
      next
    }
    rows <- h_scores[[fov_name]]$my_rows
    cols <- h_scores[[fov_name]]$my_cols
    
    # Slice the top-1 assignments for this FoV's pixels
    fov_top2 <- fov_h_long_top2[start_idx[fi]:end_idx[fi], ]
    
    # Encode comp.id as integer factor with ALL levels for consistent colours
    comps_factor <- factor(fov_top2$comp.id, levels = all_comp_ids)
    comp_matrix  <- matrix(as.integer(comps_factor),
                           nrow  = rows, ncol = cols, byrow = TRUE)
    
    # pheatmap colour vector must be length = number of factor levels
    pheatmap_cols <- col_mineral[levels(comps_factor)]
    
    pheatmap::pheatmap(
      comp_matrix,
      color         = pheatmap_cols,
      breaks        = seq(0.5, length(all_comp_ids) + 0.5, by = 1),
      cluster_rows  = FALSE,
      cluster_cols  = FALSE,
      scale         = "none",
      cellwidth     = 5,
      cellheight    = 5,
      show_rownames = FALSE,
      show_colnames = FALSE,
      border_color  = NA,
      legend        = FALSE,
      main          = paste0(sample_name, " — top 2nd mineral  |  ", fov_name)
    )
  }  
  
  # Do the same for the top3 identity
  group_lengths <- sapply(fov.meta, function(x) length(x[["id"]]))
  end_idx        <- cumsum(group_lengths)
  start_idx      <- c(1L, head(end_idx, -1L) + 1L)
  
  for (fi in seq_along(names(fov_h_FoV))) {
    
    fov_name <- names(fov_h_FoV)[fi]
    
    if (!fov_name %in% names(h_scores)) {
      warning("    FoV '", fov_name, "' not in h_scores — skipping map.")
      next
    }
    rows <- h_scores[[fov_name]]$my_rows
    cols <- h_scores[[fov_name]]$my_cols
    
    # Slice the top-1 assignments for this FoV's pixels
    fov_top3 <- fov_h_long_top_third[start_idx[fi]:end_idx[fi], ]
    
    # Encode comp.id as integer factor with ALL levels for consistent colours
    comps_factor <- factor(fov_top3$comp.id, levels = all_comp_ids)
    comp_matrix  <- matrix(as.integer(comps_factor),
                           nrow  = rows, ncol = cols, byrow = TRUE)
    
    # pheatmap colour vector must be length = number of factor levels
    pheatmap_cols <- col_mineral[levels(comps_factor)]
    
    pheatmap::pheatmap(
      comp_matrix,
      color         = pheatmap_cols,
      breaks        = seq(0.5, length(all_comp_ids) + 0.5, by = 1),
      cluster_rows  = FALSE,
      cluster_cols  = FALSE,
      scale         = "none",
      cellwidth     = 5,
      cellheight    = 5,
      show_rownames = FALSE,
      show_colnames = FALSE,
      border_color  = NA,
      legend        = FALSE,
      main          = paste0(sample_name, " — top 3rd mineral  |  ", fov_name)
    )
  }
  
  # ==================================================================
  # PLOT C — Co-occurrence barplots (top-1 vs top-2, then top-1 vs top-3)
  # ==================================================================
  # fov_h_long_top3_occ: pixels × top (1,2,3), saved by 05b
  fov_top3_occ_dt <- as.data.table(fov_h_long_top3_occ)
  # Rename columns to top1, top2, top3 for clarity
  setnames(fov_top3_occ_dt, c("top1", "top2", "top3"))

  # Count top2 occurrences for each top1 value
  occ_bar <- fov_top3_occ_dt[
    , .(occ = .N), by = .(top1, top2)
  ][order(top1, -occ)]

  occ_wide <- dcast(occ_bar, top1 ~ top2, value.var = "occ", fill = 0L)

  # Bar order: x-axis = top1 mineral ID; stacked fill = top2 mineral
  top2_ids  <- setdiff(names(occ_wide), "top1")
  mat_occ   <- as.matrix(occ_wide[, -1])
  rownames(mat_occ) <- occ_wide$top1

  # Map top1 x-labels to mineral names
  x_labels <- sapply(occ_wide$top1, function(id) {
    nm <- cfg$comp_names[id]; if (is.na(nm)) paste0("comp",id) else nm
  })
  # Map top2 fill colours
  fill_cols <- col_mineral[top2_ids]
  fill_cols[is.na(fill_cols)] <- "#AAAAAA"

  par(mar = c(8, 5, 4, 10), mfrow = c(2, 1))

  # Raw counts
  barplot(t(mat_occ),
          names.arg = x_labels,
          col       = fill_cols,
          las       = 2,
          ylab      = "Number of pixels",
          main      = paste(sample_name, "— top-2 mineral co-occurring with top-1"))
  legend("topright", inset = c(-0.25, 0),
         legend = sapply(top2_ids, function(id) {
           nm <- cfg$comp_names[id]; if (is.na(nm)) paste0("comp",id) else nm
         }),
         fill = fill_cols, xpd = TRUE, bty = "n", cex = 0.8)

  # Percentage: normalise row-wise so each top-1 bar sums to 100 %
  row_s   <- rowSums(mat_occ)
  mat_rel <- mat_occ / ifelse(row_s == 0, 1, row_s) * 100
  barplot(t(mat_rel),
          names.arg = x_labels,
          col       = fill_cols,
          las       = 2,
          ylab      = "% of co-occurring pixels",
          main      = paste(sample_name, "— top-2 mineral co-occurring with top-1 (%)"))

  par(mfrow = c(1, 1))

  # ==================================================================
  # PLOT D — Colour legend page
  # ==================================================================
  par(mar = c(2, 2, 3, 2))
  plot(NULL, xaxt = "n", yaxt = "n", bty = "n",
       ylab = "", xlab = "",
       xlim = 0:1, ylim = 0:1,
       main = paste(sample_name, "— mineral colour legend"))
  legend("topleft",
         legend  = sapply(all_comp_ids, function(id) {
           nm <- cfg$comp_names[id]; if (is.na(nm)) paste0("comp",id) else nm
         }),
         pch     = 19,
         col     = col_mineral[all_comp_ids],
         bty     = "n",
         ncol    = 3,
         y.intersp = 1,
         cex     = 1)

  dev.off()
  message("  PDF saved: ", pdf_path)
  message("=== Done: ", sample_name, " ===\n")
}

message("\n=== 11_stacked_barplots complete ===")
