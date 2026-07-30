# =============================================================
# 10_nmf_component_maps.R
#
# Generate spatial heatmaps of NMF component weights for every
# kept component × every FoV, for all samples.
#
# Prerequisite: run AFTER scripts 05b and 06 so that
#   - cfg$rdata_file exists (post-NMF .RData from 05b)
#   - cfg$hscore_file exists  (h-score .rds from 06, used for
#     FoV grid dimensions: my_rows, my_cols)
#
# One PDF per sample is written to paths$plots.
# Each page shows one component × one FoV spatial map.
#
# Run from main.R (add a source() call after Step 6) or
# interactively after sourcing config.R + raman_functions.R.
# =============================================================

options(stringsAsFactors = FALSE)
options(scipen = 999)

library(pheatmap)

# ---- Load pipeline configuration if not already in environment ----
if (!exists("paths")) source("config.R")

# Color ramp shared across all maps
# Uses YlOrRd (yellow → orange → red) to show low → high weights
PALETTE <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)

# ====================================================================
# Main loop over samples
# ====================================================================
for (sample_name in names(samples)) {

  message("\n=== nmf_component_maps: ", sample_name, " ===\n")
  cfg <- samples[[sample_name]]

  # Skip if config.R component fields are still NULL
  if (is.null(cfg$comp_names)) {
    message("  comp_names is NULL in config.R — skipping ", sample_name,
            ".\n  Run 05a, fill config.R, then re-run this script.")
    next
  }

  # ------------------------------------------------------------------
  # 1. Load post-NMF objects (saved by 05b_post_nmf_assign.R)
  #    Loads: fov_h_FoV, fov_w.identity.cos, fov_h, fov.meta, ...
  # ------------------------------------------------------------------
  if (!file.exists(cfg$rdata_file)) {
    warning("  RData not found: ", cfg$rdata_file, " — skipping ", sample_name)
    next
  }
  load(cfg$rdata_file)
  message("  Loaded: ", cfg$rdata_file)

  # ------------------------------------------------------------------
  # 2. Load h-score object (saved by 06_heterogeneity_score.R)
  #    Used only to recover per-FoV grid dimensions (my_rows, my_cols).
  # ------------------------------------------------------------------
  if (!file.exists(cfg$hscore_file)) {
    warning("  h-score file not found: ", cfg$hscore_file, " — skipping ", sample_name)
    next
  }
  h_scores <- readRDS(cfg$hscore_file)
  message("  Loaded h-scores: ", cfg$hscore_file)

  # ------------------------------------------------------------------
  # 3. Build a lookup: fov_name -> (rows, cols)
  #    fov_h_FoV names and h_scores names must match — both derive
  #    from the same corrected-file basenames (e.g. "CR_FoV.1_325nm").
  # ------------------------------------------------------------------
  fov_dims <- lapply(names(fov_h_FoV), function(fn) {
    if (!fn %in% names(h_scores)) {
      warning("    FoV '", fn, "' not found in h_scores — will skip.")
      return(NULL)
    }
    list(rows = h_scores[[fn]]$my_rows,
         cols = h_scores[[fn]]$my_cols)
  })
  names(fov_dims) <- names(fov_h_FoV)

  # ------------------------------------------------------------------
  # 4. Open per-sample PDF
  # ------------------------------------------------------------------
  if (!dir.exists(paths$plots)) dir.create(paths$plots, recursive = TRUE)

  pdf_path <- file.path(
    paths$plots,
    paste0(sample_name, "_nmf_component_maps_",
           format(Sys.Date(), "%Y-%m-%d"), ".pdf")
  )
  pdf(file = pdf_path, width = 8, height = 7)
  message("  Writing PDF: ", pdf_path)

  # ------------------------------------------------------------------
  # 5. Loop: kept components × FoVs
  #    fov_w.identity.cos$id lists the kept component IDs (from 05b).
  #    cfg$comp_names maps id -> human-readable mineral label.
  # ------------------------------------------------------------------
  kept_ids <- fov_w.identity.cos$id   # character vector, e.g. "1", "2", ...

  for (comp_id_chr in kept_ids) {

    comp_id  <- as.integer(comp_id_chr)
    # Human-readable label from config.R; fall back to "comp<id>" if missing
    label    <- cfg$comp_names[comp_id_chr]
    if (is.na(label) || is.null(label)) label <- paste0("comp", comp_id_chr)

    # Color breaks: use a tighter range for low-weight components
    # (all maps use 0 – 0.2 with 20 steps; easy to adjust per component here)
    map_breaks <- seq(0, 0.2, 0.01)    # 21 breakpoints → 20 color steps
    col_pal    <- PALETTE(length(map_breaks) - 1)

    for (fov_name in names(fov_h_FoV)) {

      dim_info <- fov_dims[[fov_name]]
      if (is.null(dim_info)) next   # h_scores entry was missing

      rows <- dim_info$rows
      cols <- dim_info$cols

      # Reshape the component's weight vector into the spatial grid.
      # Pixels are stored row-wise (left→right, top→bottom).
      weight_mat <- matrix(
        as.vector(fov_h_FoV[[fov_name]][comp_id, ]),
        nrow  = rows,
        ncol  = cols,
        byrow = TRUE
      )

      pheatmap::pheatmap(
        weight_mat,
        scale         = "none",
        cluster_rows  = FALSE,
        cluster_cols  = FALSE,
        cellwidth     = 5,
        cellheight    = 5,
        legend        = TRUE,
        show_rownames = FALSE,
        show_colnames = FALSE,
        border_color  = FALSE,
        color         = col_pal,
        breaks        = map_breaks,
        main          = paste0(sample_name,
                               "  |  comp ", comp_id_chr,
                               "  —  ", label,
                               "\n", fov_name)
      )
    }   # end FoV loop
  }   # end component loop

  dev.off()
  message("  PDF saved: ", pdf_path)
  message("=== Done: ", sample_name, " ===\n")
}

message("\n=== nmf_component_maps complete ===")
