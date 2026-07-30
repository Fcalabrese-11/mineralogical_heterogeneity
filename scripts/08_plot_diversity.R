# =============================================================
# 08_plot_diversity.R
# Generate all cross-sample diversity comparison plots.
#
# Loads the combined .RData produced by 07_diversity_indices.R
# and (optionally) the RaoQ .RData from 09_raoq.R, then generates
# publication-ready boxplots for the retained metrics:
#   - Gini index (original and after threshold)
#   - Shannon entropy
#   - Heterogeneity score (H-score)
#   - Rao's Q quadratic entropy (edge-corrected)
#
# Each section uses either:
#   het_stats_melt_*       — long format (variable × value), for paired
#                            with-biomass / without-biomass comparisons
#   all_het_stats_noresin  — wide format, resin pixels excluded, for
#                            identity-stratified "standard" plots
#
# Prerequisite: 07_diversity_indices.R has been run.
#               09_raoq.R has been run (optional; RaoQ plots skipped if absent).
# =============================================================

options(max.print = 1000)
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# ====================================================================
# 1. Load combined diversity objects from 07_diversity_indices.R
# ====================================================================
div_files <- list.files(paths$diversity_output,
                        pattern = "diversity_all_samples.*\\.RData$",
                        full.names = TRUE)
if (length(div_files) == 0)
  stop("No diversity .RData found in ", paths$diversity_output,
       "\nRun 07_diversity_indices.R first.")

load(div_files[which.max(file.mtime(div_files))])
# Loaded: gini_sd_hs_combined, metric_combined,
#         het_stats_melt_gini_original, het_stats_melt_gini_threshold,
#         het_stats_melt_gini_minus, het_stats_melt_gini_original_minus,
#         threshold_combined, gini_drop_loop_all, gini_ref_all

message("Diversity objects loaded.")
message("Samples: ",
        paste(unique(gini_sd_hs_combined$sample), collapse = ", "))

# ====================================================================
# 2. Load RaoQ objects from 09_raoq.R (optional)
# ====================================================================
raoq_available <- FALSE
raoq_files <- list.files(paths$diversity_output,
                         pattern = "raoq_all_samples.*\\.RData$",
                         full.names = TRUE)
if (length(raoq_files) > 0) {
  load(raoq_files[which.max(file.mtime(raoq_files))])
  # Loaded: raoq_combined  (columns: sample, fov, px_id,
  #         rao_all, rao_nobio, rao_nores, rao_nobio_nores,
  #         rao_all_nodom, rao_nobio_nodom, rao_nores_nodom, rao_nobio_nores_nodom)
  raoq_available <- TRUE
  message("RaoQ objects loaded.")
} else {
  message("No RaoQ .RData found — RaoQ plots will be skipped. Run 09_raoq.R first.")
}

# ====================================================================
# 3. Build all_het_stats_noresin
#
# Wide-format per-pixel table: resin pixels excluded, with transformed
# Gini columns and (if available) RaoQ columns.
#
# Column mapping:
#   gini_nores_original       = 1 - gini_nores       (original, with bio)
#   gini_nobio_nores_original = 1 - gini_nobio_nores  (original, without bio)
#   gini_thres_noresin        = 1 - gini_thres_noresin   (after threshold)
#   gini_thres_nobio_noresin  = 1 - gini_thres_nobio_noresin
#   het_score_original        = het_score
#
# RaoQ columns (added if raoq_available):
#   rao_nores         — edge-corrected, noresin subset (thresholded)
#   rao_nobio_nores   — edge-corrected, nobio_nores subset (thresholded)
# ====================================================================
all_het_stats_noresin <- copy(
  gini_sd_hs_combined[!grepl("Resin", identity, ignore.case = TRUE)]
)
all_het_stats_noresin[, sample_original := sample]
all_het_stats_noresin[, `:=`(
  gini_nores_original       = 1 - gini_nores,
  gini_nobio_nores_original = 1 - gini_nobio_nores,
  gini_thres_noresin        = 1 - gini_thres_noresin,
  gini_thres_nobio_noresin  = 1 - gini_thres_nobio_noresin,
  gini_thres_minus          = 1 - gini_thres_minus,
  gini_original_minus       = 1 - gini_original_minus,
  het_score_original        = het_score
)]

if (raoq_available) {
  # Add a within-sample sequential index matching gini_sd_hs px_id ordering.
  # Both 07 and 09 iterate fov_h_FoV in the same names() order, so
  # sequential row index within each sample group aligns correctly.
  raoq_combined[, global_px_id := seq_len(.N), by = sample]

  raoq_sub <- raoq_combined[, .(
    sample,
    global_px_id,
    rao_all,
    rao_nobio,
    rao_nores,
    rao_nobio_nores,
    rao_all_nodom,
    rao_nobio_nodom,
    rao_nores_nodom,
    rao_nobio_nores_nodom
  )]

  all_het_stats_noresin <- merge(
    all_het_stats_noresin,
    raoq_sub,
    by.x = c("sample", "px_id"),
    by.y = c("sample", "global_px_id"),
    all.x = TRUE
  )
  message("RaoQ columns joined to all_het_stats_noresin.")
}

# ====================================================================
# Shared theme
# ====================================================================
raman_theme <- function(legend = "none") {
  theme_light() +
    theme(
      text            = element_text(size = 14),
      axis.text       = element_text(size = 14),
      strip.text      = element_text(size = 14),
      axis.title.y    = element_text(size = 14),
      legend.position = legend
    )
}

raman_theme_small <- function(legend = "none") {
  theme_light() +
    theme(
      text            = element_text(size = 10),
      axis.text       = element_text(size = 12),
      strip.text      = element_text(size = 14),
      axis.title.y    = element_text(size = 10),
      legend.position = legend
    )
}

# Pairwise sample comparisons for Wilcoxon significance brackets, built from
# whichever samples are actually present in the loaded data (ordered by the
# config `samples` list). No hardcoded sample names — add a sample in config.R
# and it is compared automatically.
sample_levels <- intersect(names(samples), unique(all_het_stats_noresin$sample))
wilcox_comparisons <- if (length(sample_levels) >= 2)
  combn(sample_levels, 2, simplify = FALSE) else list()

# ====================================================================
# Open PDF output
# ====================================================================
if (!dir.exists(paths$plots)) dir.create(paths$plots, recursive = TRUE)
pdf_path <- file.path(paths$plots,
                      paste0("diversity_plots_", Sys.Date(), ".pdf"))
pdf(pdf_path, width = 10, height = 8)
message("Writing plots to: ", pdf_path)

# ====================================================================
# SECTION A — Gini index, ORIGINAL (no threshold)
#   Data: het_stats_melt_gini_original (long format)
#   variables: gini_nores_original, gini_nobio_nores_original
# ====================================================================
message("\n--- Section A: Gini original ---")

# A1: Per-sample paired comparison (with vs. without biomass)
p_a1 <- ggplot(het_stats_melt_gini_original,
               aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(staplewidth = 0.25, fill = "lightgreen", color = "darkgreen") +
  stat_compare_means(
    aes(group = variable),
    comparisons = list(c("gini_nores_original", "gini_nobio_nores_original")),
    method = "wilcox.test", label = "p.signif", bracket.size = 0.6
  ) +
  facet_wrap(~ sample_original, scales = "fixed") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "darkred") +
  labs(title = "Gini index with and without biomass (no threshold)",
       x = "with vs without biomass (no resin pixels)",
       y = "Gini index (1 - gini)") +
  raman_theme()
print(p_a1)

# A2: Facet grid: sample × duo_identity (bio vs. other)
p_a2 <- ggplot(het_stats_melt_gini_original,
               aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(staplewidth = 0.25, fill = "lightgreen", color = "darkgreen") +
  stat_summary(geom = "crossbar", fun = "median", fun.min = "median", fun.max = "median",
               linewidth = 0.25, width = 0.75) +
  scale_y_continuous(limits = c(0, NA)) +
  stat_compare_means(
    aes(group = variable),
    comparisons = list(c("gini_nores_original", "gini_nobio_nores_original")),
    method = "wilcox.test", label = "p.signif", bracket.size = 0.6
  ) +
  facet_grid(duo_identity ~ sample_original, scales = "fixed") +
  labs(title = "Gini index with and without biomass (others = min + ker)",
       y = "Gini index (1 - gini)") +
  raman_theme()
print(p_a2)

# A3: Facet grid: sample × bioker_identity (bio+ker vs. minerals)
p_a3 <- ggplot(het_stats_melt_gini_original,
               aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(staplewidth = 0.25, fill = "palegreen", color = "darkgreen") +
  stat_summary(geom = "crossbar", fun = "median", fun.min = "median", fun.max = "median",
               linewidth = 0.25, width = 0.75) +
  scale_y_continuous(limits = c(0, 1.2)) +
  stat_compare_means(
    aes(group = variable),
    comparisons = list(c("gini_nores_original", "gini_nobio_nores_original")),
    method = "wilcox.test", label = "p.signif", bracket.size = 0.6
  ) +
  facet_grid(bioker_identity ~ sample_original, scales = "fixed") +
  labs(title = "Gini index with and without biomass (bio + ker)",
       y = "Gini index (1 - gini)") +
  raman_theme()
print(p_a3)

# A4: Cross-sample: bio vs. other (ker + min) facet
p_a4 <- ggplot(het_stats_melt_gini_original,
               aes(x = variable, y = value, fill = sample_original)) +
  geom_boxplot(staplewidth = 0.25) +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~ duo_identity, scales = "fixed") +
  ggtitle("Gini index: biomass vs. other (ker+min) across samples — no threshold") +
  raman_theme()
print(p_a4)

# A5: Cross-sample: bio vs. other (ker + min) facet
p_a5 <- ggplot(het_stats_melt_gini_original,
               aes(x = variable, y = value, fill = sample_original)) +
  geom_boxplot(staplewidth = 0.25) +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~ bioker_identity, scales = "fixed") +
  ggtitle("Gini index: biomass (ker+min) vs. other across samples — no threshold") +
  raman_theme()
print(p_a5)

# ====================================================================
# SECTION C — Gini index, AFTER THRESHOLD
#   Data: het_stats_melt_gini_threshold (long format)
# ====================================================================
message("\n--- Section C: Gini threshold ---")

# C1: Per-sample paired comparison
p_c1 <- ggplot(het_stats_melt_gini_threshold,
               aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(staplewidth = 0.25, fill = "lightgreen", color = "darkgreen") +
  scale_y_continuous(limits = c(0, NA)) +
  stat_compare_means(
    aes(group = variable),
    comparisons = list(c("gini_thres_noresin", "gini_thres_nobio_noresin")),
    method = "wilcox.test", label = "p.signif", bracket.size = 0.6
  ) +
  facet_wrap(~ sample_original, scales = "fixed") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "darkred") +
  labs(title = "Gini index with and without biomass (after threshold)",
       x = "with vs without biomass (no resin pixels)",
       y = "Gini index (1 - gini)") +
  raman_theme()
print(p_c1)

# C2: Facet grid: sample × duo_identity
p_c2 <- ggplot(het_stats_melt_gini_threshold,
               aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(staplewidth = 0.25, fill = "lightgreen", color = "darkgreen") +
  stat_summary(geom = "crossbar", fun = "median", fun.min = "median", fun.max = "median",
               linewidth = 0.25, width = 0.75) +
  scale_y_continuous(limits = c(0, NA)) +
  stat_compare_means(
    aes(group = variable),
    comparisons = list(c("gini_thres_noresin", "gini_thres_nobio_noresin")),
    method = "wilcox.test", label = "p.signif", bracket.size = 0.6
  ) +
  facet_grid(sample_original ~ duo_identity, scales = "fixed") +
  labs(title = "Gini index with and without biomass (others = min+ker) after threshold",
       y = "Gini index (1 - gini)") +
  raman_theme()
print(p_c2)

# C3: Facet grid: sample × duo_identity
p_c3 <- ggplot(het_stats_melt_gini_threshold,
               aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(staplewidth = 0.25, fill = "lightgreen", color = "darkgreen") +
  stat_summary(geom = "crossbar", fun = "median", fun.min = "median", fun.max = "median",
               linewidth = 0.25, width = 0.75) +
  scale_y_continuous(limits = c(0, NA)) +
  stat_compare_means(
    aes(group = variable),
    comparisons = list(c("gini_thres_noresin", "gini_thres_nobio_noresin")),
    method = "wilcox.test", label = "p.signif", bracket.size = 0.6
  ) +
  facet_grid(sample_original ~ bioker_identity, scales = "fixed") +
  labs(title = "Gini index with and without biomass (bio+ker) after threshold",
       y = "Gini index (1 - gini)") +
  raman_theme()
print(p_c3)

# ====================================================================
# SECTION E — Shannon entropy (after threshold, all components)
#   Data: threshold_combined
# ====================================================================
message("\n--- Section E: Shannon entropy ---")

# E1: Shannon by dominant component — global Kruskal-Wallis per facet
p_e1 <- ggplot(threshold_combined,
               aes(x = identity, y = shannon, fill = sample_original)) +
  geom_boxplot(staplewidth = 0.25) +
  stat_compare_means(method = "kruskal.test", label.y.npc = "top",
                     label = "p.format") +
  facet_wrap(~ sample_original, scales = "fixed") +
  ggtitle("Shannon entropy by dominant component — global Kruskal-Wallis (after threshold)") +
  raman_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_e1)

# E2: Shannon by bioker_identity — pairwise Wilcoxon (biomass vs mineral)
p_e2 <- ggplot(threshold_combined,
               aes(x = bioker_identity, y = shannon, fill = bioker_identity)) +
  geom_boxplot(staplewidth = 0.25) +
  stat_compare_means(comparisons = list(c("biomass", "mineral")),
                     method = "wilcox.test", label = "p.signif") +
  facet_wrap(~ sample_original, scales = "fixed") +
  labs(title = "Shannon entropy: biomass vs mineral — pairwise Wilcoxon (after threshold)",
       x = NULL, y = "Shannon entropy") +
  raman_theme_small()
print(p_e2)

# ====================================================================
# SECTION F — STANDARD PLOTS: Gini original (all_het_stats_noresin)
# ====================================================================
message("\n--- Section F: Standard plots — Gini original ---")

# F1: duo_identity × facet sample, Gini with biomass
p_f1 <- ggplot(all_het_stats_noresin,
               aes(x = duo_identity, y = gini_nores_original, fill = duo_identity)) +
  geom_boxplot(staplewidth = 0.25) +
  stat_compare_means(comparisons = list(c("biomass", "other")),
                     method = "wilcox.test", label = "p.signif") +
  facet_wrap(~ sample_original, scales = "fixed") +
  labs(title = "Gini index - with bio (vs. ker + min) no threshold",
       x = "biomass component included in the calculation, resin pixels excluded from plot",
       y = "Gini coefficient (1-gini)") +
  raman_theme_small()
print(p_f1)

# F2: bioker_identity × facet sample, Gini with biomass (ker+bio)
p_f2 <- ggplot(all_het_stats_noresin,
               aes(x = bioker_identity, y = gini_nores_original, fill = bioker_identity)) +
  geom_boxplot(staplewidth = 0.25) +
  stat_compare_means(aes(group = bioker_identity),
                     comparisons = list(c("biomass", "mineral")),
                     method = "wilcox.test", label = "p.signif") +
  facet_wrap(~ sample_original, scales = "fixed") +
  labs(title = "Gini index - with bio (ker+bio) no threshold",
       x = "biomass component included in the calculation, resin pixels excluded from plot",
       y = "Gini coefficient (1-gini)") +
  raman_theme_small()
print(p_f2)

# F3: bioker_identity × facet sample, Gini without biomass
p_f3 <- ggplot(all_het_stats_noresin,
               aes(x = bioker_identity, y = gini_nobio_nores_original, fill = bioker_identity)) +
  geom_boxplot(staplewidth = 0.25) +
  stat_compare_means(aes(group = bioker_identity),
                     comparisons = list(c("biomass", "mineral")),
                     method = "wilcox.test", label = "p.signif") +
  facet_wrap(~ sample_original, scales = "fixed") +
  labs(title = "Gini index - without bio (ker+bio) no threshold",
       x = "biomass component excluded from calculation, resin pixels excluded from plot",
       y = "Gini coefficient (1-gini)") +
  raman_theme_small()
print(p_f3)

# F4: Cross-sample, facet by bioker_identity, Gini with biomass
p_f4 <- ggplot(all_het_stats_noresin,
               aes(x = sample_original, y = gini_nores_original, fill = bioker_identity)) +
  geom_boxplot(staplewidth = 0.25) +
  ylim(c(0.1, NA)) +
  stat_compare_means(aes(group = sample_original),
                     comparisons = wilcox_comparisons,
                     method = "wilcox.test", label = "p.signif") +
  facet_wrap(~ bioker_identity, scales = "fixed") +
  labs(title = "Gini coefficient across samples with biomass — biomass (bio+ker) vs. minerals",
       x = "",
       y = "Gini coefficient (1-gini)") +
  raman_theme_small(legend = "right")
print(p_f4)

# F5: Cross-sample, facet by bioker_identity, Gini without biomass
p_f5 <- ggplot(all_het_stats_noresin,
               aes(x = sample_original, y = gini_nobio_nores_original, fill = bioker_identity)) +
  geom_boxplot(staplewidth = 0.25) +
  ylim(c(0.1, NA)) +
  stat_compare_means(aes(group = sample_original),
                     comparisons = wilcox_comparisons,
                     method = "wilcox.test", label = "p.signif") +
  facet_wrap(~ bioker_identity, scales = "fixed") +
  labs(title = "Gini coefficient across samples without biomass — biomass (bio+ker) vs. minerals",
       x = "biomass component subtracted from calculation, resin pixels excluded from plot",
       y = "Gini coefficient (1-gini)") +
  raman_theme_small(legend = "right")
print(p_f5)

# F6: Across-sample, no pixel identity, Gini with biomass
p_f6 <- ggplot(all_het_stats_noresin,
               aes(x = sample_original, y = gini_nores_original)) +
  geom_boxplot(staplewidth = 0.25) +
  stat_compare_means(aes(group = sample_original),
                     comparisons = wilcox_comparisons,
                     method = "wilcox.test", label = "p.signif") +
  labs(title = "Gini index across samples with biomass",
       x = "",
       y = "Gini index (1-gini)") +
  raman_theme()
print(p_f6)

# F7: Across-sample, no pixel identity, Gini without biomass
p_f7 <- ggplot(all_het_stats_noresin,
               aes(x = sample_original, y = gini_nobio_nores_original)) +
  geom_boxplot(staplewidth = 0.25) +
  stat_compare_means(aes(group = sample_original),
                     comparisons = wilcox_comparisons,
                     method = "wilcox.test", label = "p.signif") +
  labs(title = "Gini index across samples without biomass",
       x = "biomass component subtracted from calculation, resin pixels excluded from plot",
       y = "Gini index (1-gini)") +
  raman_theme()
print(p_f7)

# ====================================================================
# SECTION G — STANDARD PLOTS: Gini threshold (all_het_stats_noresin)
# ====================================================================
message("\n--- Section G: Standard plots — Gini threshold ---")

# G1: duo_identity × facet sample (with biomass)
p_g1 <- ggplot(all_het_stats_noresin,
               aes(x = duo_identity, y = gini_thres_noresin, fill = duo_identity)) +
  geom_boxplot(staplewidth = 0.25) +
  stat_compare_means(comparisons = list(c("biomass", "other")),
                     method = "wilcox.test", label = "p.signif") +
  facet_wrap(~ sample_original, scales = "fixed") +
  labs(title = "Gini index - with bio (vs. ker + min) after threshold",
       x = "biomass component included in the calculation, resin pixels excluded from plot",
       y = "Gini coefficient (1-gini)") +
  raman_theme_small()
print(p_g1)

# G2: bioker_identity × facet sample (with biomass, ker+bio)
p_g2 <- ggplot(all_het_stats_noresin,
               aes(x = bioker_identity, y = gini_thres_noresin, fill = bioker_identity)) +
  geom_boxplot(staplewidth = 0.25) +
  stat_compare_means(aes(group = bioker_identity),
                     comparisons = list(c("biomass", "mineral")),
                     method = "wilcox.test", label = "p.signif") +
  facet_wrap(~ sample_original, scales = "fixed") +
  labs(title = "Gini index - with bio (ker+bio) after threshold",
       x = "biomass component included in the calculation, resin pixels excluded from plot",
       y = "Gini coefficient (1-gini)") +
  raman_theme_small()
print(p_g2)

# G3: bioker_identity × facet sample (without biomass)
p_g3 <- ggplot(all_het_stats_noresin,
               aes(x = bioker_identity, y = gini_thres_nobio_noresin, fill = bioker_identity)) +
  geom_boxplot(staplewidth = 0.25) +
  stat_compare_means(aes(group = bioker_identity),
                     comparisons = list(c("biomass", "mineral")),
                     method = "wilcox.test", label = "p.signif") +
  facet_wrap(~ sample_original, scales = "fixed") +
  labs(title = "Gini index - without bio (ker+bio) after threshold",
       x = "biomass component excluded from calculation, resin pixels excluded from plot",
       y = "Gini coefficient (1-gini)") +
  raman_theme_small()
print(p_g3)

# ====================================================================
# SECTION J — Heterogeneity score (H-score)
#   Data: all_het_stats_noresin$het_score_original
# ====================================================================
message("\n--- Section J: Heterogeneity score ---")

# J1: duo_identity × facet sample (bio vs. other)
p_j1 <- ggplot(all_het_stats_noresin,
               aes(x = duo_identity, y = het_score_original, fill = duo_identity)) +
  geom_boxplot(staplewidth = 0.25) +
  ylim(c(0, 0.8)) +
  stat_compare_means(comparisons = list(c("biomass", "other")),
                     method = "wilcox.test", label = "p.signif") +
  facet_wrap(~ sample_original, scales = "fixed") +
  labs(title = "Heterogeneity score biomass (min+ker = mineral)",
       x = "bio vs mineral pixels across samples - resin pixels excluded",
       y = "Heterogeneity score (1 - HS)") +
  raman_theme_small()
print(p_j1)

# J2: bioker_identity × facet sample (bio+ker vs. minerals)
p_j2 <- ggplot(all_het_stats_noresin,
               aes(x = bioker_identity, y = het_score_original, fill = bioker_identity)) +
  geom_boxplot(staplewidth = 0.25) +
  ylim(c(0, 0.8)) +
  stat_compare_means(comparisons = list(c("biomass", "mineral")),
                     method = "wilcox.test", label = "p.signif") +
  facet_wrap(~ sample_original, scales = "fixed") +
  labs(title = "Heterogeneity score biomass (ker+bio = biomass)",
       x = "bio vs mineral pixels across samples - resin pixels excluded",
       y = "Heterogeneity score (1 - HS)") +
  raman_theme_small()
print(p_j2)

# J3: Cross-sample, facet by duo_identity
p_j3 <- ggplot(all_het_stats_noresin,
               aes(x = sample_original, y = het_score_original, fill = duo_identity)) +
  geom_boxplot(staplewidth = 0.25) +
  scale_y_continuous(limits = c(0, 0.8)) +
  stat_compare_means(comparisons = wilcox_comparisons,
                     method = "wilcox.test", label = "p.signif") +
  facet_wrap(~ duo_identity, scales = "fixed") +
  labs(title = "Heterogeneity score (other = min + ker)",
       x = "bio vs mineral pixels across samples",
       y = "Heterogeneity score (1 - HS)") +
  raman_theme()
print(p_j3)

# J4: Cross-sample, facet by bioker_identity
p_j4 <- ggplot(all_het_stats_noresin,
               aes(x = sample_original, y = het_score_original, fill = bioker_identity)) +
  geom_boxplot(staplewidth = 0.25) +
  scale_y_continuous(limits = c(0, 0.8)) +
  stat_compare_means(comparisons = wilcox_comparisons,
                     method = "wilcox.test", label = "p.signif") +
  facet_wrap(~ bioker_identity, scales = "fixed") +
  labs(title = "Heterogeneity score (biomass = bio + ker)",
       x = "bio vs mineral pixels across samples",
       y = "Heterogeneity score (1 - HS)") +
  raman_theme()
print(p_j4)

# J5: Across-sample, no pixel identity
p_j5 <- ggplot(all_het_stats_noresin,
               aes(x = sample_original, y = het_score_original)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  #geom_boxplot(staplewidth = 0.05, outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 0.8)) +
  stat_compare_means(
    comparisons = wilcox_comparisons,
    method = "wilcox.test", label = "p.signif"
  ) +
  labs(title = "Heterogeneity score across samples",
       x = "",
       y = "Heterogeneity score (1 - HS)") +
  raman_theme()
print(p_j5)

# median.quartile <- function(x){
#   out <- quantile(x, probs = c(0.25,0.5,0.75))
#   names(out) <- c("ymin","y","ymax")
#   return(out) 
# }
# 
p_j6<- ggplot(all_het_stats_noresin,
       aes(x = sample_original, y = het_score_original)) +
  geom_violin()+
  #stat_summary(fun = median.quartile,fun.y = median, geom='point')+
scale_y_continuous(limits = c(0, 0.8)) +
  stat_compare_means(comparisons = wilcox_comparisons,
    method = "wilcox.test", label = "p.signif"  ) +
  labs(title = "Heterogeneity score across samples",
       x = "",
       y = "Heterogeneity score (1 - HS)") +
  raman_theme()
#print(p_j6)

# ====================================================================
# SECTION K — Rao's Q quadratic entropy (edge-corrected)
#   Data: all_het_stats_noresin with RaoQ columns from 09_raoq.R
#   Skipped if 09_raoq.R has not been run.
#
#   Column mapping from 09_raoq.R:
#     rao_nores         = edge-corrected, noresin subset (with bio, thresholded)
#     rao_nobio_nores   = edge-corrected, nobio_nores subset (thresholded)
# ====================================================================
if (raoq_available && "rao_nores" %in% names(all_het_stats_noresin)) {

  message("\n--- Section K: Rao's Q ---")

  # K1: duo_identity × facet sample (with biomass)
  p_k1 <- ggplot(all_het_stats_noresin,
                 aes(x = duo_identity, y = rao_nores, fill = duo_identity)) +
    geom_boxplot(staplewidth = 0.25) +
    ylim(c(0, NA)) +
    stat_compare_means(comparisons = list(c("biomass", "other")),
                       method = "wilcox.test", label = "p.signif") +
    facet_wrap(~ sample_original, scales = "fixed") +
    labs(title = "RaoQ (edge-corrected) with biomass — bio vs. ker+min",
         x = "bio vs mineral pixels across samples",
         y = "Rao's Q quadratic entropy index (edge corrected)") +
    raman_theme_small()
  print(p_k1)

  # K2: bioker_identity × facet sample (with biomass, ker+bio)
  p_k2 <- ggplot(all_het_stats_noresin,
                 aes(x = bioker_identity, y = rao_nores, fill = bioker_identity)) +
    geom_boxplot(staplewidth = 0.25) +
    ylim(c(0, NA)) +
    stat_compare_means(comparisons = list(c("biomass", "mineral")),
                       method = "wilcox.test", label = "p.signif") +
    facet_wrap(~ sample_original, scales = "fixed") +
    labs(title = "RaoQ (edge-corrected) with biomass — bio+ker vs. minerals",
         x = "bio vs mineral pixels across samples - resin pixels excluded",
         y = "Rao's Q quadratic entropy index (edge corrected)") +
    raman_theme_small()
  print(p_k2)

  # K3: bioker_identity × facet sample (without biomass)
  p_k3 <- ggplot(all_het_stats_noresin,
                 aes(x = bioker_identity, y = rao_nobio_nores, fill = bioker_identity)) +
    geom_boxplot(staplewidth = 0.25) +
    ylim(c(0, NA)) +
    stat_compare_means(comparisons = list(c("biomass", "mineral")),
                       method = "wilcox.test", label = "p.signif") +
    facet_wrap(~ sample_original, scales = "fixed") +
    labs(title = "RaoQ (edge-corrected) without biomass — bio+ker vs. minerals",
         x = "bio vs mineral pixels across samples - resin pixels excluded",
         y = "Rao's Q quadratic entropy index (edge corrected)") +
    raman_theme_small()
  print(p_k3)

  # K4: Cross-sample, facet by duo_identity (with biomass)
  p_k4 <- ggplot(all_het_stats_noresin,
                 aes(x = sample_original, y = rao_nores, fill = duo_identity)) +
    geom_boxplot(staplewidth = 0.25) +
    scale_y_continuous(limits = c(0, NA)) +
    stat_compare_means(comparisons = wilcox_comparisons,
                       method = "wilcox.test", label = "p.signif") +
    facet_wrap(~ duo_identity, scales = "fixed") +
    labs(title = "RaoQ (edge-corrected) with biomass — other = min+ker",
         x = "bio vs mineral pixels across samples",
         y = "Rao's Q quadratic entropy index (edge corrected)") +
    raman_theme()
  print(p_k4)

  # K5: Cross-sample, facet by bioker_identity (with biomass)
  p_k5 <- ggplot(all_het_stats_noresin,
                 aes(x = sample_original, y = rao_nores, fill = bioker_identity)) +
    geom_boxplot(staplewidth = 0.25) +
    scale_y_continuous(limits = c(0, NA)) +
    stat_compare_means(comparisons = wilcox_comparisons,
                       method = "wilcox.test", label = "p.signif") +
    facet_wrap(~ bioker_identity, scales = "fixed") +
    labs(title = "RaoQ (edge-corrected) with biomass — biomass = bio+ker",
         x = "bio vs mineral pixels across samples",
         y = "Rao's Q quadratic entropy index (edge corrected)") +
    raman_theme()
  print(p_k5)

  # K6: Across-sample, no pixel identity, with biomass
  p_k6 <- ggplot(all_het_stats_noresin,
                 aes(x = sample_original, y = rao_nores)) +
    geom_boxplot(staplewidth = 0.25, outlier.shape = NA) +
    scale_y_continuous(limits = c(0, NA)) +
    stat_compare_means(
      comparisons = wilcox_comparisons,
      method = "wilcox.test", label = "p.signif"
    ) +
    labs(title = "RaoQ (edge-corrected) with biomass",
         x = "",
         y = "Rao's Q quadratic entropy index(edge corr))") +
    raman_theme()
  print(p_k6)

  # K7: Across-sample, no pixel identity, without biomass
  p_k7 <- ggplot(all_het_stats_noresin,
                 aes(x = sample_original, y = rao_nobio_nores)) +
    geom_boxplot(staplewidth = 0.25, outlier.shape = NA) +
    scale_y_continuous(limits = c(0, NA)) +
    stat_compare_means(
      comparisons = wilcox_comparisons,
      method = "wilcox.test", label = "p.signif"
    ) +
    labs(title = "RaoQ (edge-corrected) without biomass",
         x = "",
         y = "Rao's Q quadratic entropy index (edge corrected)") +
    raman_theme()
  print(p_k7)
  
  # K8: bioker_identity × facet sample (without biomass)
  p_k8 <- ggplot(all_het_stats_noresin,
                 aes(x = bioker_identity, y = rao_all_nodom, fill = bioker_identity)) +
    geom_boxplot(staplewidth = 0.25) +
    ylim(c(0, NA)) +
    stat_compare_means(comparisons = list(c("biomass", "mineral")),
                       method = "wilcox.test", label = "p.signif") +
    facet_wrap(~ sample_original, scales = "fixed") +
    labs(title = "RaoQ (edge-corrected) without dominant — bio+ker vs. minerals",
         x = "bio vs mineral pixels across samples - resin pixels excluded",
         y = "Rao's Q quadratic entropy index (edge corrected)") +
    raman_theme_small()
  print(p_k8)
  
} else {
  message("Skipping Section K (RaoQ): data not available.")
}


# ====================================================================
# SECTION N — Gini drop-loop diagnostic plots
#   For each sample: shows how Gini changes when each component is
#   removed one at a time (leave-one-out component analysis).
#
#   Uses gini_drop_loop_all / gini_ref_all saved by 07_diversity_indices.R.
#   Skipped if those objects are absent.
#
#   Per sample, 5 plots:
#     1. Density distributions — reference + all leave-one-out vectors
#     2. Boxplot of raw values per removal
#     3. Boxplot of Δ (removal minus reference)
#     4. Boxplot of normalised Δ (relative to reference)
#     5. Summary: median |relative Δ| per component (most influential = tallest bar)
# ====================================================================
if (exists("gini_drop_loop_all") && length(gini_drop_loop_all) > 0) {
  message("\n--- Section N: Gini drop-loop plots ---")

  # Helper: build the five diagnostic plots for one index (SD or Gini)
  # mat_ref : named numeric vector (one value per pixel), the reference
  # drop_list : named list of such vectors, one per leave-one-out removal
  # index_name : "SD" or "Gini" (used in plot titles)
  # sample_name: passed in for titles
  plot_drop_loop <- function(mat_ref, drop_list, index_name, sn) {

    n_comp   <- length(drop_list)
    comp_nms <- names(drop_list)
    cols     <- rainbow(n_comp + 1)

    # Build full matrix: reference + all leave-one-out vectors (pixels × columns)
    full_mat  <- do.call(cbind, c(list(reference = mat_ref), drop_list))
    # Build comp-only matrix (without the reference column)
    comp_mat  <- do.call(cbind, drop_list)

    # 1. Density distributions -----------------------------------------------
    par(mar = c(4, 4, 4, 2))
    all_vals <- unlist(full_mat)
    x_range  <- range(all_vals[is.finite(all_vals)], na.rm = TRUE)
    dens_list <- lapply(seq_len(ncol(full_mat)), function(j) {
      v <- full_mat[, j]
      density(v[is.finite(v)], from = x_range[1], to = x_range[2], n = 256)
    })
    y_max <- max(sapply(dens_list, function(d) max(d$y)))
    plot(dens_list[[1]]$x, dens_list[[1]]$y,
         type = "l", lty = 2, lwd = 1.5, col = cols[1],
         xlim = x_range, ylim = c(0, y_max * 1.1),
         xlab = paste(index_name, "per pixel"), ylab = "Density",
         main = paste(sn, "—", index_name, "distributions (drop-loop)"))
    for (j in seq_along(dens_list)[-1])
      lines(dens_list[[j]]$x, dens_list[[j]]$y, lty = 2, lwd = 1.5, col = cols[j])
    legend("topright", c("reference", comp_nms), col = cols,
           lty = 2, cex = 0.6, ncol = 2)

    # 2. Boxplot of raw values -------------------------------------------------
    par(mar = c(9, 4, 4, 2))
    boxplot(as.data.frame(full_mat), outline = FALSE, col = cols,
            names = colnames(full_mat), las = 2,
            ylab = paste(index_name, "per pixel"),
            main = paste(sn, "— Effect of removing one component on", index_name))

    # 3. Delta boxplot ---------------------------------------------------------
    delta <- sweep(comp_mat, 1, mat_ref, "-")
    par(mar = c(9, 5, 3, 3))
    boxplot(as.data.frame(delta), outline = FALSE,
            names = colnames(delta), las = 2,
            ylab = bquote(Delta ~ .(index_name)),
            main = paste(sn, "— Influence of removing each component on", index_name))
    abline(h = 0, col = "red", lty = 2)

    # 4. Normalised delta boxplot ----------------------------------------------
    #NA_real_ ensures the result stays numeric (double type), avoiding type coercion issues.
    safe_ref <- ifelse(mat_ref == 0, NA_real_, mat_ref)
    delta_norm <- sweep(delta, 1, safe_ref, "/")
    par(mar = c(9, 5, 3, 3))
    boxplot(as.data.frame(delta_norm), outline = FALSE,
            names = colnames(delta_norm), las = 2,
            ylab = paste("Relative", bquote(Delta), index_name,
                         "(fraction of reference)"),
            main = paste(sn, "— Relative influence (normalised)"))
    abline(h = 0, col = "red", lty = 2)

    # 5. Summary: median |relative Δ| per component ---------------------------
    rel_change <- apply(delta_norm, 2, function(x) median(abs(x), na.rm = TRUE))
    par(mar = c(9, 5, 4, 3))
    barplot(sort(rel_change, decreasing = TRUE),
            las = 2, col = "steelblue",
            ylab = paste("Median |relative", paste0("\u0394", index_name), "|"),
            main = paste(sn, "— Relative contribution of each component to", index_name))
  }

  for (sn in names(gini_drop_loop_all)) {
    # Gini drop-loop: transform to 1-gini convention before plotting
    gini_ref_1mg  <- 1 - gini_ref_all[[sn]]
    gini_drop_1mg <- lapply(gini_drop_loop_all[[sn]], function(v) 1 - v)

    # Save raw delta (gini_without_i - gini_original) per component per pixel as CSV
    gini_delta_mat <- do.call(cbind, lapply(gini_drop_loop_all[[sn]],
                                            function(v) v - gini_ref_all[[sn]]))
    colnames(gini_delta_mat) <- names(gini_drop_loop_all[[sn]])
    out_delta_dir <- file.path(paths$diversity_output, sn)
    if (!dir.exists(out_delta_dir)) dir.create(out_delta_dir, recursive = TRUE)
    write.table(as.data.frame(gini_delta_mat),
                file  = file.path(out_delta_dir,
                                  paste0(sn, "_gini_minus_original_droploop.txt")),
                quote = FALSE, sep = "\t", row.names = FALSE)
    message("  Gini drop-loop delta saved: ", sn)

    plot_drop_loop(gini_ref_1mg, gini_drop_1mg, index_name = "1 - Gini", sn = sn)
  }

} else {
  message("Skipping Section N (Gini drop-loop plots): objects not found.",
          " Re-run 07_diversity_indices.R to generate them.")
}

# ====================================================================
# SECTION O — Gini minus (threshold): dominant-removed vs original
#   Data: het_stats_melt_gini_minus (long format, both as 1-gini)
#   variables: gini_thres_noresin, gini_thres_minus
# ====================================================================
if (exists("het_stats_melt_gini_minus")) {
  message("\n--- Section O: Gini minus (threshold) ---")

  # O1: Per-sample paired comparison
  p_o1 <- ggplot(het_stats_melt_gini_minus,
                 aes(x = variable, y = value, fill = variable)) +
    geom_boxplot(staplewidth = 0.25, fill = "lightgreen", color = "darkgreen") +
    stat_compare_means(
      aes(group = variable),
      comparisons = list(c("gini_thres_noresin", "gini_thres_minus")),
      method = "wilcox.test", label = "p.signif", bracket.size = 0.6
    ) +
    facet_wrap(~ sample_original, scales = "fixed") +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "darkred") +
    labs(title = "Gini (threshold): original vs dominant-removed",
         x = "original vs dominant-removed (no resin pixels)",
         y = "Gini index (1 - gini)") +
    raman_theme()
  print(p_o1)

  # O2: Facet grid: sample × bioker_identity
  p_o2 <- ggplot(het_stats_melt_gini_minus,
                 aes(x = variable, y = value, fill = variable)) +
    geom_boxplot(staplewidth = 0.25, fill = "lightgreen", color = "darkgreen") +
    stat_summary(geom = "crossbar", fun = "median",
                 fun.min = "median", fun.max = "median",
                 linewidth = 0.25, width = 0.75) +
    scale_y_continuous(limits = c(0, NA)) +
    stat_compare_means(
      aes(group = variable),
      comparisons = list(c("gini_thres_noresin", "gini_thres_minus")),
      method = "wilcox.test", label = "p.signif", bracket.size = 0.6
    ) +
    facet_grid(bioker_identity ~ sample_original, scales = "fixed") +
    labs(title = "Gini (threshold): original vs dominant-removed — by identity",
         y = "Gini index (1 - gini)") +
    raman_theme()
  print(p_o2)
}

# ====================================================================
# SECTION P — Gini minus (original, non-thresholded):
#             dominant-removed vs original
#   Data: het_stats_melt_gini_original_minus (long format, both as 1-gini)
#   variables: gini_original, gini_original_minus
# ====================================================================
if (exists("het_stats_melt_gini_original_minus")) {
  message("\n--- Section P: Gini minus (original, non-thresholded) ---")

  # P1: Per-sample paired comparison
  p_p1 <- ggplot(het_stats_melt_gini_original_minus,
                 aes(x = variable, y = value, fill = variable)) +
    geom_boxplot(staplewidth = 0.25, fill = "palegreen", color = "darkgreen") +
    stat_compare_means(
      aes(group = variable),
      comparisons = list(c("gini_original", "gini_original_minus")),
      method = "wilcox.test", label = "p.signif", bracket.size = 0.6
    ) +
    facet_wrap(~ sample_original, scales = "fixed") +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "darkred") +
    labs(title = "Gini (non-thresholded): original vs dominant-removed",
         x = "original vs dominant-removed (no resin pixels)",
         y = "Gini index (1 - gini)") +
    raman_theme()
  print(p_p1)

  # P2: Facet grid: sample × bioker_identity
  p_p2 <- ggplot(het_stats_melt_gini_original_minus,
                 aes(x = variable, y = value, fill = variable)) +
    geom_boxplot(staplewidth = 0.25, fill = "palegreen", color = "darkgreen") +
    stat_summary(geom = "crossbar", fun = "median",
                 fun.min = "median", fun.max = "median",
                 linewidth = 0.25, width = 0.75) +
    scale_y_continuous(limits = c(0, NA)) +
    stat_compare_means(
      aes(group = variable),
      comparisons = list(c("gini_original", "gini_original_minus")),
      method = "wilcox.test", label = "p.signif", bracket.size = 0.6
    ) +
    facet_grid(bioker_identity ~ sample_original, scales = "fixed") +
    labs(title = "Gini (non-thresholded): original vs dominant-removed — by identity",
         y = "Gini index (1 - gini)") +
    raman_theme()
  print(p_p2)
}

dev.off()
message("PDF saved: ", pdf_path)

message("\n=== 08_plot_diversity.R complete ===")
message("Plots generated:")
message("  A1-A5  : Gini original (long format)")
message("  C1-C2  : Gini threshold (long format)")
message("  E1-E2  : Shannon entropy")
message("  F1-F7  : Standard plots — Gini original")
message("  G1-G3  : Standard plots — Gini threshold")
message("  J1-J5  : Heterogeneity score")
if (raoq_available) message("  K1-K8  : Rao's Q (edge corrected)")
message("  N      : Gini drop-loop diagnostics (per sample, 5 plots each)")
message("  O-P    : Gini minus (dominant-removed vs original)")
