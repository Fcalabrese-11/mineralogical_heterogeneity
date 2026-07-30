# ============================================================
# RAMAN ANALYSIS PIPELINE — MAIN ENTRY POINT
# ============================================================
#
# HOW TO USE:
#   1. Open config.R and fill in all folder paths under `paths`.
#      Also review the per-sample NMF parameters (component IDs,
#      excluded components, mineral names) in the `samples` list.
#
#   2. Run this file from top to bottom inside RStudio, or from
#      the terminal with:
#        Rscript main.R
#
#   3. Each step is a separate source() call.  You can run any
#      individual step by highlighting it and pressing Ctrl+Enter,
#      without re-running the whole pipeline.
#
#   4. NMF (Step 4) is SLOW (hours for large maps).  Run it once,
#      then set RUN_NMF <- FALSE on subsequent runs.
#
# WHAT EACH STEP DOES:
#   01  Build the interpolated reference spectral library
#       (RRUFF + custom references → Ref*.RData)
#   02  Correct raw biomass reference spectra
#       (SG smoothing + ALS baseline → .txt files)
#   03  Correct raw sample FoV maps
#       (SG first, then ALS baseline → corrected_SG_BL_*.txt)
#   04  Run NMF decomposition on corrected FoV maps
#       (per-sample, per-rank → *.rds)
#   05a Explore NMF components interactively
#       (cosine similarity table + diagnostic PDF → you fill in config.R)
#   05b Assign pixel identities using your config.R settings
#       (cosine similarity vs. reference → *_NMF*.RData + CSVs)
#   06  Compute spectral heterogeneity scores per pixel
#       (cosine dissimilarity to neighbours → *.rds)
#   07  Compute diversity indices (Gini, Shannon)
#       (per pixel, per sample → combined .RData)
#   08  Generate cross-sample diversity comparison plots
#       (boxplots by sample and mineral identity)
#   09  Rao's Q quadratic entropy (spatial diversity, 3x3 window)
#   10  Spatial heatmaps of NMF component weights per FoV
#   11  Stacked barplots + categorical mineral identity maps
#   12  Merge per-sample outputs + summary figures (boxplots + stats)
#
# ============================================================

# ---- Anchor the working directory to this file's folder ----
# Lets the pipeline be launched from anywhere (`Rscript path/to/main.R`, or
# source()d from another directory): all the relative source() calls and the
# relative paths in config.R resolve against the pipeline folder regardless.
.this_file <- tryCatch({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 1 && nzchar(f)) normalizePath(f)
  else normalizePath(sys.frame(1)$ofile)
}, error = function(e) NA_character_)
if (!is.na(.this_file)) setwd(dirname(.this_file))

# ---- Load shared configuration and helper functions ----
source("config.R") # paths, common_x, samples, palette, annotation params
source("functions/raman_functions.R") # interpolation, range01, etc.

# ============================================================
# Step 2: Correct biomass reference spectra   [runs before step 1]
# ============================================================
# Applies Savitzky-Golay + ALS baseline to raw biomass spectra and
# writes them to paths$corrected_biomass. When custom biomass
# references are used they are folded into the reference library, so
# this correction must happen BEFORE the library is built. Script 02
# is therefore sourced ahead of script 01.
source("scripts/02_preprocess_biomass.R")

# ============================================================
# Step 1: Build the reference spectral library
# ============================================================
# Reads the RRUFF references (paths$rruff_refs) and the custom
# references (paths$custom_refs — the corrected biomass from step 2),
# interpolates to common_x, merges, and saves an .RData library.
# Re-run whenever references change.
source("scripts/01_build_reference_library.R")

# ============================================================
# Step 3: Correct sample spectra (SG → baseline)
# ============================================================
# Reads raw FoV .txt files from paths$raw_data and applies
# SG smoothing (p=2, n=15) followed by ALS (lambda=4, p=0.0001).
source("scripts/03_preprocess_samples.R")

# ============================================================
# Step 4: Run NMF decomposition
# ============================================================
# NOTE: This step is computationally expensive (hours per sample).
#       Set RUN_NMF <- FALSE to skip it once .rds files exist.
#
# The NMF script accepts command-line arguments:
#   Rscript scripts/04_run_nmf.R <sample_name> <rank>
# For example:
#   Rscript scripts/04_run_nmf.R CR 15
#   Rscript scripts/04_run_nmf.R VR 15
#   Rscript scripts/04_run_nmf.R MN 15
#
# When sourced from here, it reads the first sample in `samples`
# with the default rank from that sample's n_components.
# For a proper run, execute the Rscript commands above in a terminal.
RUN_NMF <- FALSE   # <--- change to TRUE to run NMF from this script
if (RUN_NMF) {
  for (sn in names(samples)) {
    r <- samples[[sn]]$n_components
    args <- c(sn, r)
    source("scripts/04_run_nmf.R")
  }
}

# ============================================================
# Step 5a: Explore NMF components — identify which is which mineral
# ============================================================
# For each sample this script:
#   - loads the NMF result and the reference library
#   - computes cosine similarity between every component and
#     every reference spectrum
#   - saves a diagnostic PDF in paths$plots showing each
#     component overlaid with its top-3 reference matches
#   - prints a summary table to the console
#
# ACTION REQUIRED after this step:
#   Open config.R, scroll to the per-sample sections, and fill
#   in the four NULL fields for each sample:
#     biomass_comp, kerogen_comp, resin_comp,
#     exclude_comps, comp_names
#   Use the console table and the PDF to guide your choices.
#   Then comment out or delete the stop() line below and continue.
source("scripts/05a_nmf_explore.R")

stop(paste(
  "\n",
  "ACTION REQUIRED — read before continuing:\n",
  "  1. Open the diagnostic PDF(s) saved in paths$plots.\n",
  "  2. Open config.R and fill in the NULL component fields\n",
  "     (biomass_comp, kerogen_comp, resin_comp, exclude_comps,\n",
  "      comp_names) for each sample.\n",
  "  3. Save config.R.\n",
  "  4. Comment out or delete this stop() call.\n",
  "  5. Re-run main.R from Step 5b onwards.\n"
))

# ============================================================
# Step 5b: Assign pixel identities using your config.R settings
# ============================================================
# Loops over all samples, validates that the config.R component
# fields are filled, then:
#   - assigns each pixel to its dominant mineral component
#   - saves .RData objects and per-component per-FoV CSVs
source("scripts/05b_post_nmf_assign.R")

# ============================================================
# Step 6: Compute heterogeneity scores per pixel
# ============================================================
# Measures spectral variability at each pixel relative to its
# 8 nearest neighbours (cosine dissimilarity).
# Saves one .rds per sample in paths$h_scores.
source("scripts/06_heterogeneity_score.R")

# ============================================================
# Step 7: Compute diversity indices
# ============================================================
# Computes the retained per-pixel diversity metrics (Gini and
# Shannon entropy) for all samples, plus the Gini drop-loop analysis
# (effect of removing each NMF component on total Gini).
# Saves a combined .RData in paths$diversity_output.
source("scripts/07_diversity_indices.R")

# ============================================================
# Step 8: Generate diversity comparison plots
# ============================================================
# Loads the combined .RData from Step 7 and produces all
# cross-sample boxplots comparing diversity metrics.
source("scripts/08_plot_diversity.R")

# ============================================================
# Step 9: Compute Rao's Q quadratic entropy (spatial diversity)
# ============================================================
# Measures local compositional diversity by computing the mean
# pairwise Euclidean distance between NMF weight vectors in a
# 3×3 sliding window across each FoV map.  Produces per-pixel
# RaoQ values for 8 subsets (all/nobio/nores/nobio_nores,
# each with and without dominant-component removal).
#
# NOTE: This step is slow (double pixel loop per FoV).
#       Set RUN_RAOQ <- FALSE to skip on re-runs once the
#       raoq_all_samples_<date>.RData is already saved.
RUN_RAOQ <- TRUE
if (RUN_RAOQ) source("scripts/09_raoq.R")

# ============================================================
# Step 10: Spatial heatmaps of NMF component weights
# ============================================================
# For every kept component × every FoV, produces a pheatmap
# showing the spatial distribution of NMF weights.
# One PDF per sample is saved in paths$plots.
# Requires: 05b (rdata_file) and 06 (hscore_file) must exist.
source("scripts/10_nmf_component_maps.R")

# ============================================================
# Step 11: Stacked barplots and categorical identity maps
# ============================================================
# For each sample produces three plot types in one PDF:
#   A. Mean NMF-weight proportion per dominant component
#   B. Categorical spatial maps (dominant mineral per pixel)
#      for each FoV
#   C. Co-occurrence barplots (top-1 vs top-2 mineral, raw
#      and relative)
# Requires: 05b (rdata_file) and 06 (hscore_file) must exist.
source("scripts/11_stacked_barplots.R")

# ============================================================
# Step 12: Merge per-sample outputs + summary figures
# ============================================================
# Merges the per-sample diversity tables (Step 7) and Rao's Q tables
# (Step 9), writes a combined CSV (all_samples_indices_<date>.csv),
# and draws the summary boxplots (heterogeneity score, Gini, Shannon,
# Rao's Q) with Wilcoxon significance brackets. One PDF in paths$plots.
# Works even if samples were processed in separate sessions.
source("scripts/12_merge_and_figures.R")

message("\n=== Pipeline complete ===")
