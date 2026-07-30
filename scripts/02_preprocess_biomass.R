# =============================================================
# 02_preprocess_biomass.R
# Preprocess raw biomass reference spectra.
#
# Applies Savitzky-Golay smoothing followed by asymmetric
# least-squares (ALS) baseline correction to each biomass file.
# Writes the corrected spectra as .txt files to paths$corrected_biomass
# for use by 01_build_reference_library.R.
#
# Prerequisite: config.R loaded (uses paths$raw_biomass and paths$corrected_biomass).
# =============================================================

options(max.print = 1000)
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(data.table)
library(baseline)
library(signal)
library(lsa)

# ---- Load raw biomass files ----

# Kanno-format .asc files and our own .txt files may coexist
raw_files_kanno <- list.files(paths$raw_biomass, pattern = "\\.asc$", recursive = TRUE, full.names = TRUE)
raw_files_ours  <- list.files(paths$raw_biomass, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
raw_files       <- c(raw_files_kanno, raw_files_ours)
message("Total biomass files found: ", length(raw_files))

# Build display names
mynames_kanno <- gsub("/", "_", gsub(paste0(paths$raw_biomass, "/"), "", raw_files_kanno))
mynames_kanno <- sub("_", ".", mynames_kanno)
mynames_ours  <- basename(raw_files_ours)
mynames       <- c(mynames_kanno, mynames_ours)

# ---- Read into a combined long data.table ----

raw_biomass <- lapply(raw_files, fread, select = c(1, 2))
names(raw_biomass) <- mynames
raw_biomass <- rbindlist(raw_biomass, idcol = "sample")

# ---- Savitzky-Golay + ALS correction ----

if (!dir.exists(paths$plots)) dir.create(paths$plots, recursive = TRUE)

pdf(file.path(paths$plots, paste0("Biomass_preprocessing_all_", Sys.Date(), ".pdf")))

corrected_biomass <- lapply(unique(raw_biomass$sample), function(x) {
  print(x)
  raw_tmp     <- raw_biomass[sample == x]
  raw_tmp_mtx <- t(as.matrix(raw_tmp[, V2]))
  colnames(raw_tmp_mtx) <- raw_tmp$V1

  # Savitzky-Golay smoothing (stronger parameters for biomass)
  bc.sg.tmp <- sgolayfilt(raw_tmp_mtx, p = 3, n = 13)

  par(mfrow = c(3, 1))
  plot(as.vector(raw_tmp_mtx),  type = "l", sub = "raw",    main = x)
  plot(as.vector(bc.sg.tmp),    type = "l", sub = "sgolay", main = x)

  bc.sg.tmp <- as.matrix(bc.sg.tmp)

  # ALS baseline correction
  bc.als.tmp <- baseline(t(bc.sg.tmp), lambda = 4, p = 0.0001, method = "als", maxit = 500)
  lines(as.vector(bc.als.tmp@baseline), type = "l", col = "red")
  plot(as.vector(bc.als.tmp@corrected), type = "l", sub = "als corrected")

  return(bc.als.tmp)
})

dev.off()
names(corrected_biomass) <- mynames

# ---- Write corrected files ----

if (!dir.exists(paths$corrected_biomass)) dir.create(paths$corrected_biomass, recursive = TRUE)

for (i in unique(raw_biomass$sample)) {
  header <- paste(
    paste0("###", i),
    paste0("###corrected file ", Sys.Date()),
    sep = "\n"
  )

  tmp_dd <- data.frame(
    wvl = raw_biomass[sample == i]$V1,
    I   = as.vector(corrected_biomass[[i]]@corrected)
  )

  data_string <- capture.output(
    write.table(tmp_dd, file = "", sep = "\t", col.names = FALSE,
                row.names = FALSE, quote = FALSE)
  )

  out_file <- file.path(
    paths$corrected_biomass,
    paste0(gsub("\\.asc$|\\.txt$", "", i), ".txt")
  )

  writeLines(c(header, data_string), out_file)
}

message("Corrected biomass spectra written to ", paths$corrected_biomass)
