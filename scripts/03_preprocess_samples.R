# =============================================================
# 03_preprocess_samples.R
# Preprocess raw Raman sample maps (FoV files).
#
# Preprocessing order (SG first, then BL — as decided):
#   1. Savitzky-Golay smoothing: p=2, n=15
#   2. ALS baseline correction: lambda=4, p=0.0001
#
# Writes corrected files to paths$corrected_data, preserving
# the sample-name subfolder structure.
#
# Prerequisite: config.R loaded (uses paths$raw_data and paths$corrected_data).
# =============================================================

options(max.print = 1000)
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(data.table)
library(baseline)
library(signal)

# ---- Find raw sample files ----

# Discover raw files for every sample defined in config.R. The match pattern is
# built from names(samples), so adding a sample there is all that is needed —
# no edit to this script. Raw file names are expected to start with the sample
# id (e.g. "MN_FoV1_..._raw.txt"), which is also how the output subfolder is
# derived below.
sample_pattern <- paste0("^(", paste(names(samples), collapse = "|"), ")")
raw_files <- list.files(paths$raw_data, pattern = sample_pattern,
                        full.names = TRUE, recursive = TRUE)
message("Sample files found: ", length(raw_files))

mynames    <- basename(raw_files)
raw_sample <- lapply(raw_files, fread)
names(raw_sample) <- mynames

# ---- Correction: SG smoothing first, then ALS baseline ----

pdf(file.path(paths$plots, paste0("Raw_samples_preprocessing_SGparam_", Sys.Date(), ".pdf")))

corrected_sample <- lapply(names(raw_sample), function(x) {

  # Remove pixel-coordinate columns (V1, V2) and the wavelength row (row 1)
  raw_tmp_mtx <- as.matrix(raw_sample[[x]][-1, -c(1, 2)])

  # Step 1: Savitzky-Golay smoothing (p=2 polynomial, n=15 window)
  bc.sg.tmp <- t(apply(raw_tmp_mtx, 1, sgolayfilt, p = 2, n = 15))

  par(mfrow = c(3, 1))
  px <- 1
  plot(as.vector(raw_tmp_mtx[px, ]),   type = "l", sub = "raw",              main = x)
  plot(as.vector(bc.sg.tmp[px, ]),     type = "l", sub = "sgolay p=2, n=15", main = x)

  # Step 2: ALS baseline correction
  bc.als.tmp <- baseline(bc.sg.tmp, lambda = 4, p = 0.0001, method = "als", maxit = 500)
  lines(as.vector(bc.als.tmp@baseline[px, ]), type = "l", col = "red")
  plot(as.vector(bc.als.tmp@corrected[px, ]), type = "l",
       sub = "als corrected lambda=4, p=0.0001")

  # Re-attach pixel coordinates and wavelength row
  bc.complete.tmp <- cbind(as.matrix(raw_sample[[x]][-1, c(1, 2)]), bc.als.tmp@corrected)
  bc.complete.tmp <- rbind(as.matrix(raw_sample[[x]][1, ]), bc.complete.tmp)

  return(bc.complete.tmp)
})

dev.off()
names(corrected_sample) <- mynames

# ---- Write corrected files ----

for (i in names(corrected_sample)) {
  tmp_dd  <- corrected_sample[[i]]
  smp_name <- strsplit(i, "_")[[1]][1]   # CR, CRD, VR, or MN

  out_dir <- file.path(paths$corrected_data, smp_name)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  out_file <- file.path(out_dir,
                        gsub("raw", paste0("corrected_SG_BL_", Sys.Date()), i))

  write.table(tmp_dd, file = out_file,
              sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}

message("Corrected sample files written to ", paths$corrected_data)
