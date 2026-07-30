# =============================================================
# 13_merge_and_figures.R
# Merge the per-sample outputs, write a combined table, and produce
# the final summary figures with statistics.
#
# Merges the per-sample outputs written by earlier steps (so samples
# may be processed in separate sessions), saves a combined CSV, and
# draws the summary boxplots with Wilcoxon significance brackets.
#
# Inputs merged (one set per sample, under paths$diversity_output):
#   <sample>/<sample>_gini_sd_hetscore_original.txt   (step 07)
#   <sample>_raoq_<date>.RData                          (step 09, optional)
#
# Retained metrics: Gini, Shannon entropy, heterogeneity score, Rao's Q.
# Cross-sample comparisons are built from names(samples), so the figures
# generalise to any number of samples.
#
# Outputs (in paths$diversity_output and paths$plots):
#   all_samples_indices_<date>.csv     — merged per-pixel table
#   figures_<date>.pdf                 — final figures
#
# Prerequisite: config.R loaded; steps 07 (and ideally 09) run per sample.
# =============================================================

options(stringsAsFactors = FALSE)

library(data.table)
library(ggplot2)
library(ggpubr)

if (!exists("paths")) source("config.R")

# ====================================================================
# 1. Merge the per-sample diversity tables (step 07 output)
# ====================================================================
per_sample <- lapply(names(samples), function(sn) {
  f <- file.path(paths$diversity_output, sn,
                 paste0(sn, "_gini_sd_hetscore_original.txt"))
  if (!file.exists(f)) {
    warning("Per-sample diversity file not found, skipping: ", f)
    return(NULL)
  }
  fread(f)
})
dt.indices.new <- rbindlist(Filter(Negate(is.null), per_sample), fill = TRUE)
if (nrow(dt.indices.new) == 0)
  stop("No per-sample diversity tables found under ", paths$diversity_output,
       "\nRun 07_diversity_indices.R first.")
message("Merged diversity tables for: ",
        paste(unique(dt.indices.new$sample), collapse = ", "))

# ====================================================================
# 2. Merge the per-sample Rao's Q tables (step 09 output), if present
# ====================================================================
raoq_per_sample <- lapply(names(samples), function(sn) {
  fs <- list.files(paths$diversity_output,
                   pattern = paste0("^", sn, "_raoq_.*\\.RData$"),
                   full.names = TRUE)
  if (length(fs) == 0) return(NULL)
  e <- new.env()
  load(fs[which.max(file.mtime(fs))], envir = e)   # -> raoq_sample
  e$raoq_sample
})
raoq_all <- rbindlist(Filter(Negate(is.null), raoq_per_sample), fill = TRUE)
raoq_available <- nrow(raoq_all) > 0

if (raoq_available) {
  # within-sample sequential index matches the per-pixel ordering of step 07
  raoq_all[, global_px_id := seq_len(.N), by = sample]
  dt.indices.new <- merge(
    dt.indices.new,
    raoq_all[, .(sample, global_px_id, rao_nores, rao_nobio_nores, rao_nores_nodom)],
    by.x = c("sample", "px_id"), by.y = c("sample", "global_px_id"),
    all.x = TRUE
  )
  message("Merged Rao's Q for: ", paste(unique(raoq_all$sample), collapse = ", "))
} else {
  message("No per-sample RaoQ .RData found — Rao's Q figures skipped. Run 09_raoq.R first.")
}

# ====================================================================
# 3. Save the merged combined table
# ====================================================================
if (!dir.exists(paths$diversity_output)) dir.create(paths$diversity_output, recursive = TRUE)
combined_csv <- file.path(paths$diversity_output,
                          paste0("all_samples_indices_", Sys.Date(), ".csv"))
fwrite(dt.indices.new, combined_csv)
message("Combined per-pixel table saved: ", combined_csv)

# ====================================================================
# 4. Summary figures
# ====================================================================
# Sample comparisons built from config (no hardcoded sample names)
sample_levels      <- intersect(names(samples), unique(dt.indices.new$sample))
sample_comparisons <- if (length(sample_levels) >= 2)
  combn(sample_levels, 2, simplify = FALSE) else list()

fig_theme <- theme_light() +
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.title.y = element_text(size = 14))

no_resin <- dt.indices.new[bioker_identity != "resin"]

if (!dir.exists(paths$plots)) dir.create(paths$plots, recursive = TRUE)
pdf(file.path(paths$plots, paste0("figures_", Sys.Date(), ".pdf")),
    width = 8, height = 6)

# Fig 2 — spectral heterogeneity score across samples
print(
  ggplot(dt.indices.new, aes(sample, het_score)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 0.2, alpha = 0.4) +
    stat_compare_means(comparisons = sample_comparisons,
                       method = "wilcox.test", label = "p.signif") +
    labs(title = "Spectral Heterogeneity Score across samples",
         x = "", y = "Spectral Heterogeneity Score") +
    fig_theme
)

# Fig 4 — Gini (1 - gini) for biomass vs mineral pixels, with/without biomass
bio_mv <- intersect(c("gini_nores", "gini_nobio_nores"), names(dt.indices.new))
if (length(bio_mv) >= 2) {
  bio_gini <- melt(dt.indices.new[bioker_identity == "biomass"],
                   id.vars = "sample", measure.vars = bio_mv,
                   variable.name = "group", value.name = "gini")
  print(
    ggplot(bio_gini, aes(group, 1 - gini)) +
      geom_boxplot(outlier.shape = NA) +
      facet_wrap(~ sample) +
      stat_compare_means(method = "wilcox.test", label = "p.signif") +
      labs(title = "Gini (biomass pixels)", y = "1 - Gini") +
      fig_theme
  )
}

min_mv <- intersect(c("gini_nores", "gini_nobio_nores", "gini_original_minus"),
                    names(dt.indices.new))
if (length(min_mv) >= 2) {
  min_gini <- melt(dt.indices.new[bioker_identity == "mineral"],
                   id.vars = "sample", measure.vars = min_mv,
                   variable.name = "group", value.name = "gini")
  print(
    ggplot(min_gini, aes(group, 1 - gini)) +
      geom_boxplot(outlier.shape = NA) +
      facet_wrap(~ sample) +
      stat_compare_means(comparisons = combn(min_mv, 2, simplify = FALSE),
                         method = "wilcox.test", label = "p.signif") +
      labs(title = "Gini (mineral pixels)", y = "1 - Gini") +
      fig_theme
  )
}

# Gini panels: biomass vs mineral (several transforms)
for (yv in c("gini_original_minus", "gini_nores", "gini_thres_noresin", "gini_thres_minus")) {
  if (!yv %in% names(no_resin)) next
  print(
    ggplot(no_resin, aes(bioker_identity, 1 - get(yv))) +
      geom_boxplot(outlier.shape = NA) +
      facet_wrap(~ sample) +
      stat_compare_means(comparisons = list(c("biomass", "mineral")),
                         method = "wilcox.test", label = "p.signif") +
      labs(title = paste0("1 - ", yv), x = NULL, y = paste0("1 - ", yv)) +
      fig_theme
  )
}

# Fig 5 — Shannon (thresholded and dominant-omitted)
for (yv in c("shannon_thres_noresin", "shannon_thres_minus")) {
  if (!yv %in% names(no_resin)) next
  print(
    ggplot(no_resin, aes(bioker_identity, get(yv))) +
      geom_boxplot(outlier.shape = NA) +
      ylim(c(0, 2.5)) +
      facet_wrap(~ sample) +
      stat_compare_means(comparisons = list(c("biomass", "mineral")),
                         method = "wilcox.test", label = "p.signif") +
      labs(title = yv, x = NULL, y = "Shannon entropy") +
      fig_theme
  )
}

# Fig 6 — heterogeneity score, biomass vs mineral
print(
  ggplot(no_resin, aes(bioker_identity, het_score)) +
    geom_boxplot() +
    facet_wrap(~ sample) +
    stat_compare_means(comparisons = list(c("biomass", "mineral")),
                       method = "wilcox.test", label = "p.signif") +
    labs(title = "Heterogeneity score (biomass vs mineral)",
         x = NULL, y = "Spectral Heterogeneity Score") +
    fig_theme
)

# Fig 7 — Rao's Q panels (only if RaoQ available)
if (raoq_available) {
  for (yv in c("rao_nores", "rao_nobio_nores", "rao_nores_nodom")) {
    if (!yv %in% names(no_resin)) next
    print(
      ggplot(no_resin, aes(bioker_identity, get(yv))) +
        geom_boxplot(outlier.shape = NA) +
        facet_wrap(~ sample) +
        stat_compare_means(comparisons = list(c("biomass", "mineral")),
                           method = "wilcox.test", label = "p.signif") +
        labs(title = yv, x = NULL, y = "Rao's Q") +
        fig_theme
    )
  }
}

dev.off()
message("Figures written to ", paths$plots)
