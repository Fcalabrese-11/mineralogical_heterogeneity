# =============================================================
# config.R — Central configuration for the Raman analysis pipeline
# =============================================================
# EDIT THIS FILE BEFORE RUNNING main.R.
# All file paths and per-sample parameters are set here — this is the
# ONLY file you need to edit to run the pipeline on your own data or to
# add a new sample.
# =============================================================

# ---- Shared spectral parameters ----
common_x <- seq(100, 3000, by = 3)   # common wavenumber axis (cm-1)
min.row  <- 100                       # minimum rows a reference spectrum must have

# ---- Reference-annotation parameters ----
# mineral_id_sheet : worksheet index in mineral_id_xlsx holding the
#                    ref_name / formula / class / subclass / name columns (script 05b).
# kerogen_ref_name : full reference-spectrum name for the kerogen/burnt-organic
#                    reference. Its underscores would otherwise be truncated when
#                    deriving mineral names in 05b, so it is preserved verbatim.
#                    Set to NULL if your reference library has no kerogen spectrum.
mineral_id_sheet <- 1
kerogen_ref_name <- "Shewanella.oneidensis_on_glass_spectrum_kerogen_05"

# ---- Shared mineral colour palette ----
# Single source of truth for scripts 08, 11, 12 (sourced automatically).
# Add an entry for any new mineral you identify; minerals not listed here
# fall back to grey in the plotting code.
color_vector_minerals <- c(
  Annite          = "#CC6677",   # rose
  Anorthite       = "#AA4499",   # purple
  PhillipsiteCa   = "#332288",   # indigo
  Ferrosilite     = "#CC79A7",   # reddish purple
  Vesuvianite     = "#D55E00",   # vermillion
  Kutnohorite     = "#009E73",   # bluish green
  Aragonite       = "#56B4E9",   # sky blue
  Opal            = "#E69F00",   # orange
  Stibnite        = "#DDCC77",   # sand
  Feitknechtite   = "#F0E442",   # yellow
  Jacobsite       = "#0072B2",   # blue
  Plattnerite     = "#88CCEE",   # light blue
  Lithiophorite   = "#44AA99",   # teal
  LithiophoriteII = "#0000FF",   # pure blue
  Biomass         = "#117733",   # dark green
  Kerogen         = "#000000",   # black
  Resin           = "#999999"    # grey
)

# ---- Directory paths ----
# All inputs live under data/ and all outputs under results/. Paths are
# relative to the repo root (the folder containing this file); main.R sets
# the working directory there automatically, so they resolve however the
# pipeline is launched (RStudio or `Rscript main.R`).
#
# INPUT folders (under data/, shipped with the MN example):
#   raw_data          — raw, uncorrected FoV .txt files, one sub-folder per sample
#                       (e.g. data/raw_data/MN/)
#   corrected_data    — where script 03 writes SG+BL-corrected FoV .txt files
#                       (one sub-folder per sample)
#   raw_biomass       — raw biomass reference spectra (.asc or .txt, two-column)
#   corrected_biomass — where script 02 writes baseline-corrected biomass files
#   rruff_refs        — RRUFF database .txt files (two-column: wavenumber, intensity)
#   custom_refs       — any additional custom reference spectra (.txt, two-column)
#
# OUTPUT folders (under results/, created by the scripts if missing):
#   ref_library, r_objects, nmf_output, h_scores, csv_output,
#   diversity_output, plots
#
# SINGLE FILE:
#   mineral_id_xlsx   — Excel file with mineral identity / cluster annotations.

paths <- list(
  raw_data          = "data/raw_data/",
  corrected_data    = "data/corrected_data/",
  raw_biomass       = "data/reference_biomass/",
  corrected_biomass = "data/corrected_biomass",
  rruff_refs        = "data/reference_rruff/",
  custom_refs       = "data/corrected_biomass",
  ref_library       = "data/reference_library/",
  r_objects         = "results/r_objects/",
  nmf_output        = "results/nmf_output",
  h_scores          = "results/heterogeneity_scores",
  csv_output        = "results/csv_output",
  diversity_output  = "results/diversity_output",
  plots             = "results/plots",
  mineral_id_xlsx   = "data/mineral.id.cluster.xlsx"
)

# ---- Per-sample NMF parameters ----
# Each entry in `samples` drives scripts 03, 04, 05a, 05b, 06, 07, 09-12.
# The pipeline loops over names(samples), so adding a sample here is all that
# is required — no script edits.
#
# Fields set BEFORE running anything:
#   fov_dir          path to the corrected FoV files for this sample
#   fov_pattern      file-name pattern to match corrected FoV files
#   nmf_file         path to the saved NMF .rds for this sample (script 04 output)
#   rdata_file       path to save/load the post-NMF .RData object (script 05b output)
#   hscore_file      path to the h-score .rds for this sample (script 06 output)
#   n_components     number of NMF components (rank r) used in script 04
#   cosine_threshold minimum cosine similarity to accept a component identity
#
# Fields filled in AFTER running 05a and reviewing the diagnostic PDF:
#   biomass_comp     integer: component ID for the biomass component
#   kerogen_comp     integer: component ID for the kerogen / burnt-organic component
#   resin_comp       integer: component ID for the embedding-resin component
#   exclude_comps    integer vector: component IDs to drop (below threshold / ambiguous)
#   comp_names       named character vector: "component_id" -> "mineral name"
#                    (only IDs NOT in exclude_comps need a name; leave the whole
#                     field NULL until you have run 05a)
#
# The MN entry below is a complete, runnable example (its data ships with the
# repo). To add your own sample, copy the MN block, rename the key, and update
# the paths and component fields. See the CR / VR template blocks at the bottom.

samples <- list(

  MN = list(
    fov_dir          = file.path(paths$corrected_data, "MN"),
    fov_pattern      = ".*corrected_SG_BL.*.txt",
    nmf_file         = file.path(paths$nmf_output, "MN_nmf_corrected_15_2026-07-21.rds"),
    rdata_file       = file.path(paths$r_objects, "MN_all-FoVs_interp_3step_NMF_15_strong_corr.RData"),
    hscore_file      = file.path(paths$h_scores, "MN_FoVs.hscore.100-3000range.corrected_SG_BL_2026-07-21.rds"),
    n_components     = 15,
    cosine_threshold = 0.75,

    # ---- filled in after running 05a ----
    biomass_comp  = 14,
    kerogen_comp  = 8,
    resin_comp    = 9,
    exclude_comps = c(1, 3, 4, 10, 12, 13),
    comp_names    = c("2"="Lithiophorite", "5"="Jacobsite",
                      "6"="LithiophoriteII", "7"="Annite",
                      "8"="Kerogen", "9"="Resin", "11"="Plattnerite",
                      "14"="Biomass", "15"="Opal")
  )

)

# =============================================================
# TEMPLATE — additional samples
# =============================================================
# Copy a block below into the `samples` list above and adjust the values to
# add a sample. CR and VR are additional sample blocks kept as templates;
# their data is not shipped, so they are left OUT of `samples` above to keep
# the default run to the example (MN).
#
#   CR = list(
#     fov_dir          = file.path(paths$corrected_data, "CR"),
#     fov_pattern      = ".*corrected_SG_BL.*.txt",
#     nmf_file         = file.path(paths$nmf_output, "CR_nmf_corrected_15_2026-01-19.rds"),
#     rdata_file       = file.path(paths$r_objects, "CR_all-FoVs_interp_3step_NMF_15_strong_corr.RData"),
#     hscore_file      = file.path(paths$h_scores, "CR_FoVs.hscore.100-3000range.corrected_SG_BL_2026-03-09.rds"),
#     n_components     = 15,
#     cosine_threshold = 0.56,
#     biomass_comp  = 15,
#     kerogen_comp  = 14,
#     resin_comp    = 5,
#     exclude_comps = c(3, 6, 7, 8, 12, 13),
#     comp_names    = c("1"="Stibnite", "2"="PhillipsiteCa",
#                       "4"="Vesuvianite", "5"="Resin",
#                       "9"="Aragonite", "10"="Feitknechtite",
#                       "11"="Kutnohorite", "14"="Kerogen", "15"="Biomass")
#   ),
#
#   VR = list(
#     fov_dir          = file.path(paths$corrected_data, "VR"),
#     fov_pattern      = ".*corrected_SG_BL.*.txt",
#     nmf_file         = file.path(paths$nmf_output, "VR_nmf_corrected_15_2026-01-19.rds"),
#     rdata_file       = file.path(paths$r_objects, "VR_all-FoVs_interp_3step_NMF_15_strong_corr.RData"),
#     hscore_file      = file.path(paths$h_scores, "VR_FoVs.hscore.100-3000range.corrected_SG_BL_2026-03-09.rds"),
#     n_components     = 15,
#     cosine_threshold = 0.75,
#     biomass_comp  = 3,
#     kerogen_comp  = 8,
#     resin_comp    = 12,
#     exclude_comps = c(4, 5, 6, 7, 11, 14),
#     comp_names    = c("1"="Ferrosilite", "2"="Annite", "3"="Biomass",
#                       "8"="Kerogen", "9"="Vesuvianite", "10"="Stibnite",
#                       "12"="Resin", "13"="Anorthite", "15"="Lithiophorite")
#   )
