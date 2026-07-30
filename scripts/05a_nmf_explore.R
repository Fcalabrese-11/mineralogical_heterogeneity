# =============================================================
# 05a_nmf_explore.R
# NMF component explorer — identify which component is which mineral.
#
# Run this BEFORE filling in the component fields in config.R.
# It produces one diagnostic PDF per sample in paths$plots.
#
# For each sample the script:
#   1. Loads and interpolates the corrected FoV files.
#   2. Loads the saved NMF result (.rds from script 04).
#   3. Computes cosine similarity between every NMF component and
#      every spectrum in the reference library.
#   4. Writes a diagnostic PDF with:
#        - All component spectra (stacked overview)
#        - For each component: its spectrum overlaid with its
#          top-3 best-matching reference spectra
#        - A cosine similarity bar chart per component
#   5. Prints a summary table to the console with:
#        Component ID | Best reference match | Cosine similarity
#
# After reviewing the PDF and the console table, open config.R and
# fill in the NULL fields (biomass_comp, kerogen_comp, resin_comp,
# exclude_comps, comp_names) for each sample.
#
# Prerequisite: config.R and functions/raman_functions.R loaded.
#               Reference library .RData must exist (from script 01).
#               NMF .rds files must exist (from script 04).
# =============================================================

options(stringsAsFactors = FALSE)
options(scipen = 999)

library(data.table)
library(NMF)
library(lsa)

# ---- Load reference library ----
ref_files <- list.files(paths$ref_library, pattern = "\\.RData$", full.names = TRUE)
if (length(ref_files) == 0) stop("No reference library .RData found in: ", paths$ref_library)
ref_file  <- ref_files[which.max(file.mtime(ref_files))]
message("Loading reference library: ", ref_file)
load(ref_file)   # loads: rP_all_dt, common_x (or bins)

# Scale reference to [0, 1]
rP_all_dt <- t(apply(rP_all_dt, 1, range01))

# Ensure output directory exists
if (!dir.exists(paths$plots)) dir.create(paths$plots, recursive = TRUE)

# ====================================================================
# Loop over samples
# ====================================================================
for (sample_name in names(samples)) {

  message("\n=== 05a Exploring NMF components: ", sample_name, " ===\n")
  cfg <- samples[[sample_name]]

  # ------------------------------------------------------------------
  # 1. Load and interpolate FoV files
  # ------------------------------------------------------------------
  fovs <- list.files(path = cfg$fov_dir, pattern = cfg$fov_pattern,
                     full.names = TRUE)
  if (length(fovs) == 0) {
    warning("No FoV files found for ", sample_name, " — skipping.")
    next
  }

  fov <- list()
  for (i in fovs) {
    tmp_name    <- make.names(gsub("_corrected_SG_BL.*.txt", "", basename(i)))
    tmp_fov     <- fread(i)
    wvl         <- as.vector(as.matrix(tmp_fov[1, -c(1, 2)]))
    tmp_fov_mtx <- as.matrix(tmp_fov[-1, -c(1, 2)])
    fov[[tmp_name]] <- as.data.table(
      t(apply(tmp_fov_mtx, 1, interpolation, x = wvl, common_x = common_x))
    )
  }

  fov_mtx <- as.matrix(rbindlist(fov))

  # ------------------------------------------------------------------
  # 2. Load NMF result
  # ------------------------------------------------------------------
  if (!file.exists(cfg$nmf_file)) {
    warning("NMF file not found: ", cfg$nmf_file, " — skipping ", sample_name)
    next
  }
  res.nmf <- readRDS(cfg$nmf_file)
  fov_w   <- basis(res.nmf)    # wavelengths × components
  fov_h   <- coef(res.nmf)     # components × pixels

  n_comp <- ncol(fov_w)
  message("Number of NMF components: ", n_comp)

  # ------------------------------------------------------------------
  # 3. Cosine similarity: every component vs every reference spectrum
  # ------------------------------------------------------------------
  # Build a data.table of component spectra (one row per component)
  fov_w.intrp     <- as.data.table(t(fov_w))
  fov_w.intrp[, id := as.character(seq_len(n_comp))]
  setcolorder(fov_w.intrp, c("id", setdiff(names(fov_w.intrp), "id")))

  cos_all <- cosinPix(r.sample = copy(fov_w.intrp), reference = copy(rP_all_dt))

  # Best match per component
  cos_max <- cos_all[, .SD[which.max(value)], by = id]
  # Top 3 matches per component (for plotting)
  cos_top3 <- cos_all[order(id, -value), .SD[1:3], by = id]

  # ------------------------------------------------------------------
  # 4. Diagnostic PDF
  # ------------------------------------------------------------------
  pdf_path <- file.path(paths$plots,
                        paste0(sample_name, "_05a_component_diagnostics_", Sys.Date(), ".pdf"))
  pdf(pdf_path, width = 12, height = 8)

  # --- 4a. Stacked overview of all component spectra ---
  par(mar = c(4, 4, 4, 2), mfrow = c(1, 1))
  y_max <- n_comp * 9
  plot(seq_along(common_x), type = "n",
       ylim = c(0, y_max), xlab = "Wavenumber index", ylab = "",
       main = paste(sample_name, "— all NMF component spectra (stacked)"),
       axes = TRUE)
  for (i in seq_len(n_comp)) {
    y <- as.vector(as.matrix(fov_w[, i])) + (i - 1) * 8.5
    lines(seq_along(common_x), y, col = i)
    text(x = 5, y = (i - 1) * 8.5 + 0.5, labels = paste("Comp", i),
         cex = 0.6, col = i)
  }

  # --- 4b. Individual component plots: component + top-3 reference matches ---
  for (comp_id in seq_len(n_comp)) {

    comp_id_chr   <- as.character(comp_id)
    top3_matches  <- cos_top3[id == comp_id_chr]
    best_match    <- cos_max[id == comp_id_chr]

    comp_spectrum <- as.vector(as.matrix(fov_w[, comp_id]))
    comp_scaled   <- as.vector(scale(comp_spectrum, center = FALSE, scale = TRUE))

    # 4-row layout: component on top, then one row per top-3 reference (no overlay)
    par(mar = c(3, 4, 3, 2), mfrow = c(4, 1))

    # Row 1 — component spectrum alone
    plot(common_x, comp_scaled, type = "l",
         xlab = "cm-1", ylab = "scaled intensity",
         main = paste0("Component ", comp_id,
                       "  (best: ", best_match$variable,
                       "  cos=", round(best_match$value, 3), ")"),
         col = "blue", xaxt="n")
    axis(side=1, at=seq(0, 3000, by=100))
    abline(v = seq(0, 3000, by = 10)) #add vertical lines

    # Rows 2-4 — top-3 reference spectra, each separately (no component overlay)
    for (k in seq_len(nrow(top3_matches))) {
      ref_name     <- as.vector(top3_matches$variable[k])
      ref_spectrum <- as.vector(as.matrix(rP_all_dt[ref_name, , drop = FALSE]))
      ref_scaled   <- as.vector(scale(ref_spectrum, center = FALSE, scale = TRUE))
      cos_val      <- round(top3_matches$value[k], 3)

      plot(common_x, ref_scaled, type = "l",
           xlab = "cm-1", ylab = "scaled intensity",
           main = paste0("Ref #", k, ": ", ref_name, "  cos=", cos_val),
           col = "darkred", xaxt="n")
      axis(side=1, at=seq(0, 3000, by=100))

    }
  }

  # --- 4c. Cosine similarity heatmap: components (rows) vs top-20 references ---
  top_refs <- cos_all[, .SD[which.max(value)], by = variable][order(-value)][1:min(20, .N)]$variable
  cos_wide <- dcast(cos_all[variable %in% top_refs], id ~ variable, value.var = "value")
  cos_mat  <- as.matrix(cos_wide[, -1])
  rownames(cos_mat) <- paste("Comp", cos_wide$id)

  par(mar = c(12, 4, 4, 2), mfrow = c(1, 1))
  heatmap(cos_mat, scale = "none",
          main = paste(sample_name, "— cosine similarity (components vs top-20 refs)"),
          col  = hcl.colors(50, "YlOrRd", rev = TRUE),
          cexCol = 0.6, cexRow = 0.8, Rowv = NA)

  dev.off()
  message("Diagnostic PDF saved to: ", pdf_path)

  # ------------------------------------------------------------------
  # 5. Console summary table — use this to fill in config.R
  # ------------------------------------------------------------------
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("  COMPONENT SUMMARY FOR:", sample_name, "\n")
  cat("  Use this table to fill in config.R for script 05b.\n")
  cat(strrep("=", 70), "\n")
  cat(sprintf("  %-6s  %-6s  %-55s  %s\n",
              "Comp", "Rank", "Best reference match", "Cosine"))
  cat(strrep("-", 70), "\n")

  for (comp_id in seq_len(n_comp)) {
    comp_id_chr <- as.character(comp_id)
    top3        <- cos_top3[id == comp_id_chr]
    for (k in seq_len(nrow(top3))) {
      cat(sprintf("  %-6s  %-6s  %-55s  %.3f\n",
                  ifelse(k == 1, comp_id_chr, ""),
                  paste0("#", k),
                  substr(top3$variable[k], 1, 55),
                  top3$value[k]))
    }
    cat(strrep("-", 70), "\n")
  }

  cat("\n  >> Now open config.R and fill in these fields for", sample_name, ":\n")
  cat("       biomass_comp  — component ID whose spectrum matches your biomass reference\n")
  cat("       kerogen_comp  — component ID matching kerogen / burnt organic material\n")
  cat("       resin_comp    — component ID matching embedding resin\n")
  cat("       exclude_comps — component IDs with cosine < threshold or ambiguous identity\n")
  cat("       comp_names    — named vector: component ID (as string) -> mineral name\n")
  cat(strrep("=", 70), "\n\n")
}

message("\n=== 05a complete — review the PDFs in paths$plots, then fill in config.R ===")
message("=== When config.R is filled, run scripts/05b_post_nmf_assign.R ===\n")
