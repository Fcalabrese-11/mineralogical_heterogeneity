# =============================================================
# 04_run_nmf.R
# Run Non-negative Matrix Factorisation (NMF) on one sample.
#
# HOW TO RUN:
#   Rscript scripts/04_run_nmf.R <sample_name> <rank>
#   e.g.  Rscript scripts/04_run_nmf.R CR 15
#
# This script is intentionally NOT sourced from main.R by default
# because NMF is computationally expensive and is usually run once.
# Set RUN_NMF <- TRUE in main.R to enable it.
#
# Prerequisite: config.R must be sourced so that paths and common_x are available.
#   When run via Rscript, source config.R explicitly at the top.
#
# Output: <sample>_nmf_corrected_<rank>_<date>.rds in paths$nmf_output
# =============================================================

options(max.print = 1000)
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(data.table)
library(NMF)

# ---- Source config if running standalone (run from the pipeline root) ----
if (!exists("paths"))         source("config.R")
if (!exists("interpolation")) source("functions/raman_functions.R")

# ---- Command-line arguments ----
# When run via Rscript, read <sample_name> <rank> from the command line.
# When sourced from main.R's RUN_NMF loop, `args` is already set by the caller.
if (!exists("args")) args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: Rscript scripts/04_run_nmf.R <sample_name> <rank>")
message("Sample: ", args[1], "  Rank: ", args[2])
sample_name <- args[1]
r           <- as.numeric(args[2])

# ---- Load corrected FoV files for this sample ----

cfg  <- samples[[sample_name]]

# Bug fix #9: use common_x from config.R (seq(100, 3000, by=3)) not the old seq(100, 3300, by=3)
# Bug fix #10: both load path and save path now unified via paths from config.R

fovs <- list.files(
  path       = cfg$fov_dir,
  pattern    = cfg$fov_pattern,
  full.names = TRUE,
  recursive  = TRUE
)

if (length(fovs) == 0) stop("No FoV files found in: ", cfg$fov_dir)
message("FoV files found: ", length(fovs))

fov      <- list()
fov.meta <- list()

for (i in fovs) {
  message(i)

  tmp_name   <- gsub("\\.txt$", "", basename(i))
  tmp_fov    <- fread(i)
  wvl        <- as.vector(as.matrix(tmp_fov[1, -c(1, 2)]))
  tmp_fov_mtx <- as.matrix(tmp_fov[-1, -c(1, 2)])

  tmp.fov.intrp <- as.data.table(t(apply(tmp_fov_mtx, 1, interpolation, x = wvl, common_x = common_x)))
  fov[[tmp_name]] <- tmp.fov.intrp

  fov.meta[[tmp_name]] <- c(length(unique(tmp_fov[-1]$V1)), length(unique(tmp_fov[-1]$V2)))
}

fov     <- rbindlist(fov)
fov_mtx <- as.matrix(fov)

# ---- Scale to [0,1] ----
fov_mtx01 <- t(apply(fov_mtx, 1, range01)) + 0.0001

# ---- Run NMF ----
res.nmf <- nmf(t(fov_mtx01), r, nrun = 5, seed = 123, .opt = "vp4")

# ---- Save result ----
if (!dir.exists(paths$nmf_output)) dir.create(paths$nmf_output, recursive = TRUE)

out_file <- file.path(paths$nmf_output,
                      paste0(sample_name, "_nmf_corrected_", r, "_", Sys.Date(), ".rds"))
saveRDS(res.nmf, file = out_file)
message("NMF result saved to: ", out_file)
