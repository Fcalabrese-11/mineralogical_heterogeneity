# Raman Spectroscopy Analysis Pipeline

## Overview

This pipeline processes hyperspectral Raman maps and turns raw pixel spectra
into quantitative mineralogy and diversity metrics. It uses a ready-to-run
small example (sample `MN`). Running it on new data means editing a single
file, `config.R`.

Pipeline steps (matching the numbered scripts):

1. Build an interpolated reference spectral library (RRUFF + custom spectra).
2. Correct raw biomass reference spectra (Savitzky-Golay + ALS baseline).
   These feed the custom references in step 1, so step 2 runs **before** step 1
   whenever custom references are used; `main.R` sources them in that order.
3. Correct raw sample spectra (Savitzky-Golay + ALS baseline).
4. Decompose each sample map with Non-negative Matrix Factorization (NMF).
5. Match NMF components to reference spectral library via cosine similarity, then
   assign a mineral identity to each pixel.
6. Quantify spatial heterogeneity with Spectral heterogeneity score (SHS).
7. Compute diversity indices (Gini index, Shannon entropy).
8. Compare diversity metrics across samples.
9. Compute Rao's Q quadratic entropy in a sliding spatial window.
10. Generate per-component spatial weight maps for each FoV.
11. Generate stacked barplots and categorical mineral identity maps.
12. Merge the per-sample outputs and draw the summary figures.

The diversity metrics carried through the analysis are the **Gini index**,
**Shannon entropy**, the **spectral heterogeneity score (SHS)**, and
**Rao's Quadratic entropy**.

## Repository layout

```
.
├── main.R          # sources the steps in order
├── config.R        # the only file to edit (paths, samples, palette, params)
├── functions/      # shared helper functions
├── scripts/        # numbered pipeline steps 01–12
├── data/           # inputs 
│   ├── raw_data/ corrected_data/ reference_rruff/ reference_biomass/
│   ├── corrected_biomass/ reference_library/ mineral.id.cluster.xlsx
└── results/        # outputs (nmf_output/ heterogeneity_scores/ plots/ …)
```

All paths in `config.R` are relative to the repository root, and `main.R` sets
the working directory to its own folder, so the pipeline runs the same from
RStudio or with `Rscript main.R`.

## Requirements

- **R** ≥ 4.2
- Packages:

```r
install.packages(c(
  "data.table", "baseline", "signal", "lsa", "raster",
  "NMF", "ggplot2", "ggpubr", "pheatmap", "RColorBrewer",
  "readxl", "vegan", "ineq",
  "terra", "spdep"   
))
```

## Configuration

`config.R` holds everything sample- and machine-specific: `paths`, the shared
mineral colour palette, the annotation parameters, and the per-sample `samples`
list. Because the pipeline loops over `names(samples)` and builds all
cross-sample comparisons from that list, adding or removing a sample is a
`config.R` edit only — no script changes.

Adding a sample: copy the  `samples`, rename it, and set its paths. The component fields
(`biomass_comp`, `resin_comp`, `comp_names`, …) stay at their defaults / `NULL`
until step `05a` has been run and its diagnostic PDF reviewed; they are then
filled from that output before running `05b`.

## The example

Sample `MN` ships with the repository: its raw and corrected FoV, NMF result,
reference library, and SH-score are included, and `config.R` is pre-filled for
it. The downstream analysis and figures can be reproduced without the slow
steps by sourcing, from the repository root:

```r
source("config.R"); source("functions/raman_functions.R")
source("scripts/05b_post_nmf_assign.R")   # pixel identities
source("scripts/06_heterogeneity_score.R")
source("scripts/07_diversity_indices.R")  # Gini + Shannon
source("scripts/09_raoq.R")               # Rao's Q
source("scripts/12_merge_and_figures.R")  # merge + figures
```

Outputs land under `results/` (per-pixel tables and CSV in
`results/diversity_output/`, figures in `results/plots/`). A full run from raw
spectra steps through `main.R` from the top; steps `01`–`04` rebuild the
reference library, preprocess, and run NMF (NMF is slow).

Step `05a` is a deliberate pause: it saves a diagnostic PDF matching each NMF
component to its top reference spectra. The component fields in `config.R` are
filled from that output before step `05b` continues.

## Pipeline steps

| Script | Description |
|--------|-------------|
| `01_build_reference_library.R` | Reads RRUFF database files and the custom reference spectra, interpolates every spectrum onto `common_x`, and saves a merged `.RData` library. When custom references are used it consumes the corrected biomass from `02`, so `02` runs first. |
| `02_preprocess_biomass.R` | Applies Savitzky-Golay (p=3, n=13) and ALS baseline (λ=4, p=0.0001) to raw biomass spectra; writes corrected `.txt` files that feed the library in `01`. |
| `03_preprocess_samples.R` | Applies SG (p=2, n=15) then ALS (λ=4, p=0.0001) to all raw sample FoV maps. Output files are named `*_corrected_SG_BL_<date>.txt`. |
| `04_run_nmf.R` | Runs NMF decomposition on one sample at a specified rank. Invoked from the command line: `Rscript scripts/04_run_nmf.R MN 15`. |
| `05a_nmf_explore.R` | Per sample: computes cosine similarity between every NMF component and every reference spectrum, saves a diagnostic PDF (component spectra vs. top-3 reference matches) in `paths$plots`, and prints a summary table. Its output drives the component fields in `config.R`, filled before `05b`. |
| `05b_post_nmf_assign.R` | Validates that the component fields in `config.R` are filled, then assigns pixel-level mineral identities and saves `.RData` objects and per-component per-FoV CSVs. |
| `06_heterogeneity_score.R` | Computes the heterogeneity score (H) per pixel as the mean cosine dissimilarity to its 8 neighbours. Saves `.rds` per sample; plots pheatmap grids. |
| `07_diversity_indices.R` | Normalises NMF weights, applies a 10%-of-max threshold, removes the dominant component, then computes the retained per-pixel diversity metrics — **Shannon entropy** and the **Gini index** — plus a Gini drop-loop. Writes per-sample tables and a combined `.RData`. |
| `08_plot_diversity.R` | Loads the combined `.RData` from step 07 (and optional Rao's Q from step 09) and generates boxplots comparing the retained metrics (Gini, Shannon, heterogeneity score, Rao's Q) across all samples with Wilcoxon significance annotations. Comparisons are generated from `names(samples)`. |
| `09_raoq.R` | Computes Rao's Q quadratic entropy per pixel using a 3×3 sliding window. Applies the same normalisation + threshold + dominant-removal as step 07, then runs `edge_corrected_raoq()` for 8 component subsets. Saves spatial heatmap PDFs and a combined `.RData`. |
| `10_nmf_component_maps.R` | For every kept NMF component × every FoV, generates a spatial pheatmap of component weights (fixed YlOrRd scale 0–0.2). One PDF per sample in `paths$plots`. Requires 05b and 06. |
| `11_stacked_barplots.R` | Per sample: mean NMF-weight proportion per dominant component, categorical spatial identity maps per FoV, and top-1 vs top-2 co-occurrence barplots, plus a colour legend page. One PDF per sample in `paths$plots`. Requires 05b and 06. |
| `12_merge_and_figures.R` | Merges the per-sample outputs — `<sample>_gini_sd_hetscore_original.txt` (step 07) and `<sample>_raoq_<date>.RData` (step 09) — writes a combined `all_samples_indices_<date>.csv`, then draws the summary boxplots (heterogeneity score, Gini, Shannon, Rao's Q) with Wilcoxon brackets. Comparisons come from `names(samples)`. Works even when samples are processed in separate sessions. One PDF (`figures_<date>.pdf`) in `paths$plots`. Requires step 07 (and ideally 09). |

## Output files

| Location | Contents |
|----------|----------|
| `paths$ref_library` | `Ref.interp_step3.100.3000.*.RData` — interpolated reference matrix |
| `paths$corrected_data/<sample>/` | `*_corrected_SG_BL_<date>.txt` — corrected FoV files |
| `paths$nmf_output/` | `<sample>_nmf_corrected_<rank>_<date>.rds` — NMF result |
| `paths$r_objects/` | `<sample>_all-FoVs_interp_3step_NMF_*.RData` — post-NMF objects |
| `paths$h_scores/` | `<sample>_FoVs.hscore.*.rds` — per-FoV H-score rasters |
| `paths$csv_output/<sample>/` | Per-component per-FoV weight CSVs |
| `paths$diversity_output/<sample>/` | `*_diversity_all_threshold.txt`, `*_gini_sd_hetscore_original.txt` |
| `paths$diversity_output/` | `diversity_all_samples_<date>.RData` — combined diversity tables |
| `paths$diversity_output/` | `raoq_all_samples_<date>.RData` — combined RaoQ table; also `<sample>_raoq_<date>.RData` per sample |
| `paths$diversity_output/` | `all_samples_indices_<date>.csv` — merged per-pixel table across samples (step 12) |
| `paths$plots/` | PDFs: preprocessing diagnostics, 05a component identification, component weight maps (step 10), stacked barplots and identity maps (step 11), diversity result plots (step 08), and summary figures (`figures_<date>.pdf`, step 12) |

## Notes

- **No FoV files found** — `cfg$fov_dir` / `cfg$fov_pattern` in `config.R` must
  point at the corrected files and match their names (e.g. `.*corrected_SG_BL.*.txt`).
- **Spectral range** — all scripts use `common_x <- seq(100, 3000, by = 3)`.
  Changing it requires re-running scripts 01–07 from the beginning.
- **NMF convergence** — increasing `nrun` in `04_run_nmf.R` or adjusting the
  rank `r` improves component quality; components can be inspected after step 05.
- **Adding reference minerals** — custom spectra are two-column `.txt` files
  (wavenumber, intensity) placed in `paths$custom_refs`; re-run script 01.
- **NMF component ordering is not deterministic** — after each new NMF run,
  `05a_nmf_explore.R` is re-run to regenerate the diagnostic PDF, and the
  component fields in `config.R` are re-identified before `05b`.
