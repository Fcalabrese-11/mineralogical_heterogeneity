# =============================================================
# 05b_post_nmf_assign.R
# Assign mineral identities to NMF components and extract
# per-pixel identity maps for all samples (CR, VR, MN).
#
# Run this AFTER:
#   1. Running 05a_nmf_explore.R
#   2. Reviewing the diagnostic PDFs
#   3. Filling in the component fields in config.R:
#        biomass_comp, kerogen_comp, resin_comp,
#        exclude_comps, comp_names
#
# For each sample in `samples` (config.R) this script:
#   1.  Validates that all component fields in config.R are filled.
#   2.  Loads and interpolates FoV files.
#   3.  Loads the NMF result (.rds from script 04).
#   4.  Matches each kept component to its reference mineral
#       using cosine similarity and the config.R threshold.
#   5.  Builds comp_above_cosine (components × pixels).
#   6.  Pearson & Spearman correlation heatmaps of kept components.
#   7.  Assigns pixel-level mineral identities (fov_h_long, top3).
#   8.  Extracts biomass / resin / kerogen pixel-ID vectors.
#   9.  Splits fov_h by FoV → fov_h_FoV.
#  10.  Saves .RData (fov_h_FoV is always defined before save).
#  11.  Exports per-component per-FoV CSVs.
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
library(raster)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(readxl)

# ---- Validate config.R component fields ----
# All component assignment fields must be non-NULL before this script runs.
required_fields <- c("biomass_comp", "kerogen_comp", "resin_comp",
                     "exclude_comps", "comp_names")

config_ok <- TRUE
for (s in names(samples)) {
  missing <- required_fields[sapply(required_fields, function(f) is.null(samples[[s]][[f]]))]
  if (length(missing) > 0) {
    message("config.R — sample '", s, "': these fields are still NULL:\n  ",
            paste(missing, collapse = ", "))
    config_ok <- FALSE
  }
}
if (!config_ok) {
  stop(
    "\nSome component fields in config.R are still NULL.\n",
    "Run scripts/05a_nmf_explore.R first, review the diagnostic PDFs,\n",
    "then fill in the NULL fields in config.R before re-running this script.\n"
  )
}

# ---- Load reference library ----
ref_files <- list.files(paths$ref_library, pattern = "\\.RData$", full.names = TRUE)
if (length(ref_files) == 0) stop("No reference library .RData found in: ", paths$ref_library)
ref_file  <- ref_files[which.max(file.mtime(ref_files))]
message("Loading reference library: ", ref_file)
load(ref_file)   # loads: rP_all_dt, common_x

rP_all_dt <- t(apply(rP_all_dt, 1, range01))

# ---- Load mineral annotation table ----
ref_anno <- readxl::read_xlsx(paths$mineral_id_xlsx, sheet = mineral_id_sheet, col_names = FALSE)
colnames(ref_anno) <- c("ref_name", "formula", "class", "subclass", "name")

# ====================================================================
# Main loop over samples
# ====================================================================
for (sample_name in names(samples)) {

  message("\n=== 05b Processing sample: ", sample_name, " ===\n")
  cfg <- samples[[sample_name]]

  # ------------------------------------------------------------------
  # 2. Load and interpolate FoV files
  # ------------------------------------------------------------------
  fovs <- list.files(path = cfg$fov_dir, pattern = cfg$fov_pattern,
                     full.names = TRUE)
  if (length(fovs) == 0) {
    warning("No FoV files found for ", sample_name, " in ", cfg$fov_dir, " — skipping.")
    next
  }
  message("FoV files found: ", length(fovs))

  fov      <- list()
  fov.meta <- list()

  for (i in fovs) {
    message("  Interpolating: ", basename(i))
    tmp_name    <- make.names(gsub("_corrected_SG_BL.*.txt", "", basename(i)))
    tmp_fov     <- fread(i)
    wvl         <- as.vector(as.matrix(tmp_fov[1, -c(1, 2)]))
    tmp_fov_mtx <- as.matrix(tmp_fov[-1, -c(1, 2)])

    fov[[tmp_name]] <- as.data.table(
      t(apply(tmp_fov_mtx, 1, interpolation, x = wvl, common_x = common_x))
    )
    fov.meta[[tmp_name]] <- list(FoV = tmp_name, id = seq_len(nrow(tmp_fov_mtx)))
  }

  fov_mtx   <- as.matrix(rbindlist(fov))
  fov_mtx01 <- t(apply(fov_mtx, 1, range01)) + 0.0001  # slight offset keeps NMF non-negative

  # ------------------------------------------------------------------
  # 3. Load NMF result
  # ------------------------------------------------------------------
  if (!file.exists(cfg$nmf_file)) {
    warning("NMF file not found: ", cfg$nmf_file, " — skipping ", sample_name)
    next
  }
  res.nmf <- readRDS(cfg$nmf_file)

  fov_w <- basis(res.nmf)   # wavelengths × components
  fov_h <- coef(res.nmf)    # components × pixels

  message("fov_w: ", paste(dim(fov_w), collapse = " x "),
          "   fov_h: ", paste(dim(fov_h), collapse = " x "))

  # ------------------------------------------------------------------
  # 4. Cosine similarity — confirm component identities
  # ------------------------------------------------------------------
  fov_w.intrp <- as.data.table(t(fov_w))
  fov_w.intrp[, id := as.character(seq_len(ncol(fov_w)))]
  setcolorder(fov_w.intrp, c("id", setdiff(names(fov_w.intrp), "id")))

  fov_w.intrp.cos     <- cosinPix(r.sample = copy(fov_w.intrp), reference = copy(rP_all_dt))
  fov_w.intrp.cos_max <- fov_w.intrp.cos[, .SD[value == max(value)], by = "id"]

  # Apply threshold and user-specified exclusions
  fov_w.identity.cos <- fov_w.intrp.cos_max[
    value > cfg$cosine_threshold][, .SD[1], by = "id"]
  fov_w.identity.cos <- fov_w.identity.cos[
    !id %in% as.character(cfg$exclude_comps)]

  # Ensure identity table is sorted by numeric component id for label alignment
  fov_w.identity.cos <- fov_w.identity.cos[order(as.numeric(id))]
  message("Kept components: ", paste(fov_w.identity.cos$id, collapse = ", "))

  # Top-10 cosine matches per component (saved in .RData for reference)
  top10_cosine <- fov_w.intrp.cos[order(id, -value), .SD[1:10], by = id]

  # ------------------------------------------------------------------
  # 5. Build comp_above_cosine (components × pixels)
  #    Consistent orientation: rows = components, columns = pixels.
  # ------------------------------------------------------------------
  kept_comp_ids <- as.numeric(fov_w.identity.cos$id)
  comp_above_cosine <- fov_h[kept_comp_ids, , drop = FALSE]

  # ------------------------------------------------------------------
  # 6. Correlation heatmaps of kept components
  # ------------------------------------------------------------------
  # t(comp_above_cosine) → pixels × components, so cor() gives component × component
  pearson_coeff_heatmap   <- cor(t(comp_above_cosine))
  spearmann_coeff_heatmap <- cor(t(comp_above_cosine), method = "spearman")

  heatmap(cor(t(fov_h)), scale = "none",
          main = paste(sample_name, "— all components (Pearson)"),
          zlim = c(-1, 1),
          col  = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))

  heatmap(pearson_coeff_heatmap, scale = "none",
          main = paste(sample_name, "— kept components (Pearson)"),
          zlim = c(-1, 1), cexRow = 0.6, cexCol = 0.6,
          col  = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
          labRow = fov_w.identity.cos$variable,
          labCol = fov_w.identity.cos$id)

  heatmap(spearmann_coeff_heatmap, scale = "none",
          main = paste(sample_name, "— kept components (Spearman)"),
          zlim = c(-1, 1), cexRow = 0.6, cexCol = 0.6,
          col  = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
          labRow = fov_w.identity.cos$variable,
          labCol = fov_w.identity.cos$id)

  # ------------------------------------------------------------------
  # 7. Assign pixel identities
  # ------------------------------------------------------------------
  fov_h.dt <- data.table(fov_h)
  fov_h.dt[, comp.id := fov_w.intrp.cos_max$id]

  fov_h_long <- melt(fov_h.dt, id.vars = "comp.id",
                     value.name = "coeff", variable.name = "px")
  fov_h_long <- fov_h_long[comp.id %in% fov_w.identity.cos$id]

  fov_h_long <- merge(fov_h_long,
                      fov_w.identity.cos[, c("id", "variable")],
                      by.x = "comp.id", by.y = "id")

  fov_h_long[, mineral := gsub("_.*", "", variable)]
  # Special case: the kerogen reference name contains underscores (which the
  # gsub above would truncate) — preserve it verbatim. Name set in config.R
  # (`kerogen_ref_name`); leave that NULL if your library has no such spectrum.
  if (!is.null(kerogen_ref_name)) {
    fov_h_long[, mineral := ifelse(
      grepl(kerogen_ref_name, variable, fixed = TRUE),
      kerogen_ref_name,
      mineral
    )]
  }
  fov_h_long$class <- tolower(ref_anno$class[match(fov_h_long$mineral, ref_anno$ref_name)])

  # Top-3 components by weight for each pixel
  fov_h_long_top3 <- fov_h_long[order(-coeff), .SD[1:3], by = px][order(px)]
  fov_h_long_top3[, top := rep(1:3, length(unique(fov_h_long_top3$px)))]

  fov_h_long_top3_occ <- dcast(fov_h_long_top3, px ~ top, value.var = "comp.id")[, -1]
  fov_h_long_top3_occ_freq <- as.data.table(fov_h_long_top3_occ)[
    , .(Freq = .N), by = c("1", "2", "3")]
  colnames(fov_h_long_top3_occ_freq) <- make.names(colnames(fov_h_long_top3_occ_freq))

  # Sub-tables by rank (convenient for downstream analyses)
  fov_h_long_top1     <- fov_h_long_top3[top == 1]
  fov_h_long_top2     <- fov_h_long_top3[top == 2]
  fov_h_long_top_third <- fov_h_long_top3[top == 3]

  # ------------------------------------------------------------------
  # 8. Extract pixel-ID vectors for key biological and technical components
  # ------------------------------------------------------------------
  px_to_ids <- function(px_vec) as.numeric(gsub("V", "", px_vec))

  biomass_comp     <- fov_h_long_top3[comp.id == as.character(cfg$biomass_comp) & top == 1]
  biomass_comp_ids <- px_to_ids(biomass_comp$px)
  message("Biomass pixels (top-1): ", length(biomass_comp_ids))

  resin_comp <- fov_h_long_top3[comp.id == as.character(cfg$resin_comp) & top == 1]
  resin_ids  <- px_to_ids(resin_comp$px)
  message("Resin pixels (top-1):   ", length(resin_ids))

  kero_comp <- fov_h_long_top3[comp.id == as.character(cfg$kerogen_comp) & top == 1]
  kero_ids  <- px_to_ids(kero_comp$px)
  message("Kerogen pixels (top-1): ", length(kero_ids))

  all_px_ids <- px_to_ids(unique(fov_h_long_top3$px))

  no_biomass_comp_ids <- setdiff(all_px_ids, biomass_comp_ids)

  # Bug fix: no_bio_no.res_comp_ids defined AFTER both biomass_comp_ids and
  # resin_ids exist (was used before definition in the original MN script).
  no_bio_no.res_comp_ids <- setdiff(all_px_ids, c(biomass_comp_ids, resin_ids))

  # Bar plots: which mineral classes co-occur most with biomass pixels
  fov_h_long_no_kerogen <- if (is.null(kerogen_ref_name))
    fov_h_long_top3 else fov_h_long_top3[variable != kerogen_ref_name]
  px_top1_bio <- fov_h_long_no_kerogen[
    top == 1 & comp.id == as.character(cfg$biomass_comp)]$px

  par(mar = c(2, 4, 4, 2), mfrow = c(2, 1))
  barplot(sort(table(fov_h_long_no_kerogen[px %in% px_top1_bio][top == 2]$class),
               decreasing = TRUE),
          las = 2,
          main = paste(sample_name, "top 2nd co-occurring mineral with biomass — no kerogen"),
          ylab = "number of pixels")
  barplot(sort(table(fov_h_long_no_kerogen[px %in% px_top1_bio][top == 3]$class),
               decreasing = TRUE),
          las = 2,
          main = paste(sample_name, "top 3rd co-occurring mineral with biomass — no kerogen"),
          ylab = "number of pixels")

  # ------------------------------------------------------------------
  # 9. Split fov_h by FoV
  #    Derived from fov.meta — no hard-coded pixel index ranges.
  # ------------------------------------------------------------------
  group_lengths <- sapply(fov.meta, function(x) length(x[["id"]]))
  end_indices   <- cumsum(group_lengths)
  start_indices <- c(1, head(end_indices, -1) + 1)

  fov_h_FoV <- mapply(
    function(s, e) fov_h[, s:e, drop = FALSE],
    start_indices, end_indices,
    SIMPLIFY = FALSE
  )
  names(fov_h_FoV) <- names(fov.meta)

  # ------------------------------------------------------------------
  # 10. Save .RData  — fov_h_FoV is guaranteed to exist at this point
  # ------------------------------------------------------------------
  if (!dir.exists(dirname(cfg$rdata_file)))
    dir.create(dirname(cfg$rdata_file), recursive = TRUE)

  save(
    fov_w.identity.cos, fov_h_FoV, fov_h_long_top3, fov_h_long,
    fov_h_long_top1, fov_h_long_top2, fov_h_long_top_third,
    spearmann_coeff_heatmap, pearson_coeff_heatmap,
    fov_h, comp_above_cosine, fov_w.intrp.cos_max, top10_cosine,
    biomass_comp, biomass_comp_ids,
    resin_comp,   resin_ids,
    kero_comp,    kero_ids,
    no_biomass_comp_ids, no_bio_no.res_comp_ids,
    fov_h_long_top3_occ, fov_h_long_top3_occ_freq,
    fov.meta,
    file = cfg$rdata_file
  )
  message("Saved .RData: ", cfg$rdata_file)

  # ------------------------------------------------------------------
  # 11. Export per-component per-FoV CSVs
  # ------------------------------------------------------------------
  out_csv_dir <- file.path(paths$csv_output, sample_name)
  if (!dir.exists(out_csv_dir)) dir.create(out_csv_dir, recursive = TRUE)

  for (comp_id_chr in fov_w.identity.cos$id) {
    comp_id    <- as.numeric(comp_id_chr)
    comp_label <- cfg$comp_names[comp_id_chr]
    if (is.na(comp_label)) comp_label <- paste0("comp", comp_id_chr)

    # All FoVs combined
    write.csv(
      fov_h[comp_id, ],
      file = file.path(out_csv_dir,
                       paste0(sample_name, "_allFoV_comp", comp_id_chr,
                              "_", comp_label, ".csv"))
    )

    # Per FoV
    for (fov_name in names(fov_h_FoV)) {
      write.csv(
        fov_h_FoV[[fov_name]][comp_id, ],
        file = file.path(out_csv_dir,
                         paste0(sample_name, "_", fov_name,
                                "_comp", comp_id_chr, "_", comp_label, ".csv"))
      )
    }
  }

  message("CSVs written to: ", out_csv_dir)
  message("=== Done: ", sample_name, " ===\n")
}

message("\n=== 05b complete — .RData objects and CSVs are ready ===")
