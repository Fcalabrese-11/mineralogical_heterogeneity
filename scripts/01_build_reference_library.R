# =============================================================
# 01_build_reference_library.R
# Build the interpolated reference spectral library.
#
# Reads RRUFF database files + custom reference files,
# interpolates every spectrum onto common_x, merges them,
# and saves an .RData file that downstream scripts load.
#
# Prerequisite: config.R and functions/raman_functions.R loaded.
# Output: Ref*.RData in paths$ref_library
# =============================================================

options(max.print = 1000)
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(data.table)

# ---- RRUFF reference files ----

rruff_files <- list.files(paths$rruff_refs, full.names = TRUE)

rP_all <- list()

for (i in rruff_files) {
  message(paste("processing .. ", basename(i)))

  # Read all lines (including metadata header)
  rP <- fread(i, fill = TRUE, blank.lines.skip = TRUE, header = FALSE)
  rP_metadata <- rP[grepl("^##", rP$V1)]

  # Skip spectra not confirmed by XRD + chemical analysis
  if (!(any(grepl("X-ray diffraction and chemical analysis", rP_metadata$V1)))) {
    message(paste("Not certain! skip .. ", basename(i)))
    next
  }

  # Keep only data rows, split on comma, coerce to numeric
  rP <- rP[!grepl("^##", rP$V1)][, tstrsplit(V1, ",", fixed = TRUE)][, lapply(.SD, as.numeric)]

  if (dim(rP)[2] < 2) stop(paste("Something is wrong, check", i))

  if (dim(rP)[1] <= min.row) {
    message(paste("less than", min.row, ". skip .. ", basename(i)))
    next
  }

  rP <- rP[complete.cases(rP)]
  # Bug fix #7: V1 is already numeric after tstrsplit + as.numeric — no gsub needed

  rP_y_interpolated <- interpolation(rP$V1, rP$V2, common_x)
  rP_all[[basename(i)]] <- rP_y_interpolated
}

rP_all_dt <- do.call(rbind, rP_all)
message("RRUFF spectra loaded: ", nrow(rP_all_dt))

# Strip acquisition-metadata suffix from row names
rownames(rP_all_dt) <- gsub("__Raman_Data.*", "", rownames(rP_all_dt))

# ---- Custom (biomass + other) reference files ----

myfiles <- list.files(paths$custom_refs, full.names = TRUE)

myref_all <- list()

for (i in myfiles) {
  message(paste("processing .. ", basename(i)))

  rP <- fread(i, blank.lines.skip = TRUE, header = FALSE)
  rP <- rP[!grepl("^#", rP$V1)][, lapply(.SD, as.numeric)]

  if (dim(rP)[2] < 2) stop(paste("Something is wrong, check", i))

  if (dim(rP)[1] <= min.row) {
    message(paste("less than", min.row, ". skip .. ", basename(i)))
    next
  }

  rP <- rP[complete.cases(rP)]
  # Bug fix #7: already numeric

  rP_y_interpolated <- interpolation(rP$V1, rP$V2, common_x)
  myref_all[[basename(i)]] <- rP_y_interpolated
}

myref_all_dt <- do.call(rbind, myref_all)
message("Custom reference spectra loaded: ", nrow(myref_all_dt))

# Bug fix #8: escape dot in regex so ".txt" is not treated as "any char + txt"
rownames(myref_all_dt) <- gsub("\\.txt$", "", rownames(myref_all_dt))

# ---- Merge and save ----

stopifnot(ncol(rP_all_dt) == ncol(myref_all_dt))  # must have same number of wavelength points
rP_all_dt <- rbind(rP_all_dt, myref_all_dt)
message("Total reference spectra: ", nrow(rP_all_dt))

# Quick sanity check: confirm biomass is present
message("Biomass entries found: ", sum(grepl("Therm|biomass|Kanno|LRW", rownames(rP_all_dt), ignore.case = TRUE)))

if (!dir.exists(paths$ref_library)) dir.create(paths$ref_library, recursive = TRUE)

save(rP_all_dt, common_x,
     file = file.path(paths$ref_library,
                      paste0("Ref.interp_step3.100.3000.XRDandChem.analysis-filtered.",
                             Sys.Date(), ".RData")))

message("Reference library saved to ", paths$ref_library)
