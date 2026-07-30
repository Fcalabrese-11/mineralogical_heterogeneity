## raman_functions.R
## Shared helper functions for the Raman spectroscopy analysis pipeline.
## Source this file at the top of every script, or let main.R load it once.

library(data.table)
library(lsa)       # cosine()
library(raster)    # raster(), adjacent(), values()

# ============================================================
# 1. Spectral utility functions
# ============================================================

#' Interpolate a single spectrum onto a common wavenumber axis.
#'
#' @param x       numeric vector — measured wavenumbers
#' @param y       numeric vector — measured intensities
#' @param common_x numeric vector — target wavenumber axis
#' @return numeric vector of interpolated intensities (length = length(common_x))
interpolation <- function(x, y, common_x) {
  interp <- approx(x, y, xout = common_x, na.rm = TRUE, yleft = 0.1, yright = 0.1)
  interp$y <- pmax(interp$y, 0)  # clip negative values to 0
  return(interp$y)
}

#' Scale a numeric vector to the [0, 1] range.
#'
#' @param x numeric vector
#' @return rescaled numeric vector
range01 <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# ============================================================
# 2. Reference library builders
# ============================================================

#' Read RRUFF-format reference files, bin to a common grid, and save.
#'
#' @param ref.dir   path to folder containing reference .txt files
#' @param bins.start,bins.end,bins.by  bin grid parameters (default 0–4000 by 5)
#' @param min.row   skip files with fewer rows than this (default 100)
#' @param ref.out   output folder (trailing slash or no slash both fine)
#' @param ref.prefix filename prefix for the saved .RData (default "Ref")
#' @param average   if TRUE (default), average replicates by mineral name;
#'                  if FALSE, keep every spectrum individually
createReference <- function(ref.dir,
                            bins.start = 0, bins.end = 4000, bins.by = 5,
                            min.row    = 100,
                            ref.out,
                            ref.prefix = "Ref",
                            average    = TRUE) {

  myfiles <- list.files(ref.dir, full.names = TRUE)
  bins    <- seq(bins.start, bins.end, bins.by)
  rP_all  <- list()

  for (i in myfiles) {
    message(paste("processing .. ", basename(i)))

    rP <- fread(i, skip = "#")

    if (dim(rP)[1] <= min.row) {
      message(paste("less than", min.row, ". skip .. ", basename(i)))
      next
    }

    rP <- rP[, c(1, 2)]
    colnames(rP) <- c("V1", "V2")
    rP <- rP[complete.cases(rP)]
    rP[, V1 := as.numeric(V1)]   # Bug fix: V1 is already numeric after fread

    rP_bin  <- rP[, bin := cut(V1, bins)][, .(V2 = mean(V2)), by = bin]
    dt      <- data.table(bin = cut(bins, bins)[-1], V2 = 0)
    rP_bin2 <- merge(dt, rP_bin, by = "bin", all.x = TRUE)[, V2 := V2.y][, .(bin, V2)][order(bin)]

    rP_all[[basename(i)]] <- rP_bin2$V2
  }

  rP_all_dt <- do.call(rbind, rP_all)

  if (average) {
    # Bug fix: use Group.1 (the actual aggregated key) as rownames BEFORE deleting it
    rP_all_dt_m <- aggregate(rP_all_dt, list(gsub("_.*", "", rownames(rP_all_dt))), mean)
    rownames(rP_all_dt_m) <- rP_all_dt_m$Group.1  # correct order from aggregate()
    rP_all_dt_m$Group.1  <- NULL
  } else {
    rP_all_dt_m <- rP_all_dt
  }

  save(rP_all_dt_m, bins,
       file = paste(paste0(ref.out, ref.prefix), bins.start, bins.end, bins.by, "RData", sep = "."))
}

# ============================================================
# 3. Sample binning
# ============================================================

#' Bin a single corrected sample file to a common grid.
#'
#' @param sample.dir path to a corrected .txt FoV file
#' @param bins       numeric vector defining bin edges (must be provided explicitly)
#' @return data.table with columns: id (pixel) + one column per bin
binSample <- function(sample.dir, bins) {  # Bug fix: bins is now an explicit parameter

  message(paste0("Processing ", sample.dir))

  crf1    <- fread(sample.dir, colClasses = "numeric")
  crf1_dt <- crf1[-1, -c(1, 2)]
  names(crf1_dt) <- as.character(cut(as.numeric(crf1[1, -c(1, 2)]), breaks = bins))
  crf1_dt[, id := 1:nrow(crf1_dt)]

  crf1_dt_long <- melt(crf1_dt, id.vars = "id", variable.name = "bin", value.name = "V2")
  crf1_dt_long[, V2 := pmax(V2, 0)]
  crf1_dt_long <- crf1_dt_long[, .(V2 = mean(V2)), by = c("id", "bin")]

  crf1_dt_long <- merge(
    data.table(
      id  = rep(unique(crf1_dt_long$id), each  = length(cut(bins, bins)[-1])),
      bin = rep(cut(bins, bins)[-1], length(unique(crf1_dt_long$id))),
      V2  = 0
    ),
    crf1_dt_long, by = c("id", "bin"), all.x = TRUE
  )
  crf1_dt_long[, V2 := V2.y]

  dcast(crf1_dt_long, id ~ bin, value.var = "V2")
}

# ============================================================
# 4. Cosine similarity
# ============================================================

#' Compute a cosine-similarity matrix between two sets of spectra.
#'
#' @param A numeric matrix (n_a × p)
#' @param B numeric matrix (n_b × p)
#' @return numeric matrix (n_a × n_b)
#'
#' Bug fix: replaces a double for-loop with a vectorised computation.
cosine_similarity_matrix <- function(A, B) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  norms_A <- sqrt(rowSums(A^2))
  norms_B <- sqrt(rowSums(B^2))
  sim <- (A %*% t(B)) / outer(norms_A, norms_B)
  return(sim)
}

#' Correlate each pixel spectrum with each reference spectrum (Pearson).
#'
#' @param r.sample data.table with column `id` + spectral columns
#' @param reference numeric matrix (reference spectra × wavelengths), named rows
#' @return long data.table: id, variable (reference name), value (correlation)
corPix <- function(r.sample, reference) {
  crf1_dt_cor <- matrix(
    apply(r.sample[, -1], 1, cor, y = t(reference), use = "pairwise.complete.obs"),
    ncol = dim(reference)[1], byrow = TRUE
  )
  colnames(crf1_dt_cor) <- rownames(reference)
  crf1_dt_cor <- as.data.table(crf1_dt_cor)
  crf1_dt_cor[, id := r.sample[, 1]]
  data.table::melt(crf1_dt_cor, id.vars = "id")
}

#' Correlate each pixel spectrum with each reference spectrum (cosine).
#'
#' @param r.sample data.table with column `id` + spectral columns
#' @param reference numeric matrix, named rows
#' @return long data.table: id, variable (reference name), value (cosine similarity)
cosinPix <- function(r.sample, reference) {
  sample_ids <- r.sample[[1]]
  r.mat      <- as.matrix(r.sample[, -1])
  ref.mat    <- as.matrix(reference)

  sim_matrix <- cosine_similarity_matrix(r.mat, ref.mat)

  sim_dt <- as.data.table(sim_matrix)
  colnames(sim_dt) <- rownames(reference)
  sim_dt[, id := sample_ids]
  melt(sim_dt, id.vars = "id")
}

# ============================================================
# 5. Spatial heterogeneity functions
# ============================================================

#' Compute per-pixel cosine-dissimilarity to neighbours and return a raster.
#'
#' @param r.sample   data.table with column `id` + spectral columns
#' @param n.neigh    neighbourhood size passed to raster::adjacent()
#' @param n.row,n.col  raster dimensions
#' @param t          number of parallel threads (unused — kept for API compatibility)
#' @param by.row     whether pixels are filled row-by-row (default TRUE)
#' @return RasterLayer with mean cosine dissimilarity per pixel
corNeigh <- function(r.sample, n.neigh, n.row, n.col, t = 10, by.row = TRUE) {

  mt  <- matrix(r.sample$id, nrow = n.row, ncol = n.col, byrow = by.row)
  s   <- raster(mt)
  a   <- adjacent(s, 1:ncell(s), directions = n.neigh, sorted = TRUE, include = FALSE)

  out      <- data.frame(a)
  out$cors <- NA

  do.cor <- function(x) {
    1 - cosine(t(as.matrix(r.sample[r.sample$id %in% values(s)[x[c(1, 2)]], -1])))[1, 2]
  }
  out$cors <- apply(a, 1, do.cor)

  r_out_vals <- aggregate(unlist(out$cors), by = list(out$from), FUN = median)
  r_out      <- s[[1]]
  r_out[]    <- r_out_vals$x
  return(r_out)
}

#' Find neighbourhood co-occurrence of mineral identities at two spatial scales.
#'
#' @param r.sample   data.table with column `id` + spectral columns
#' @param neigh.d    integer vector of neighbourhood window sizes (e.g. c(3, 5))
#' @param cor.mtx    long data.table of cosine similarities (id, variable, value)
#' @param n.row,n.col  raster dimensions
#' @param t          unused (kept for API compatibility)
#' @param threshold  cosine threshold for accepting an identity (default 0.6)
#' @param by.row     whether pixels fill row-by-row (default TRUE)
#' @return data.table with co-occurrence fractions by scale and identity pair
whichNeigh <- function(r.sample, neigh.d = c(3, 5), cor.mtx,
                       n.row, n.col, t = 8, threshold = 0.6, by.row = TRUE) {

  neighs_out <- list()
  a_pre      <- matrix(ncol = 2, nrow = 0)

  for (n.neigh in neigh.d) {

    message(paste("Processing ", n.neigh))

    m.neigh <- matrix(1, ncol = n.neigh, nrow = n.neigh)
    m.neigh[ceiling(n.neigh / 2), ceiling(n.neigh / 2)] <- 0

    mt <- matrix(r.sample$id, nrow = n.row, ncol = n.col, byrow = by.row)
    s  <- raster(mt)
    a  <- adjacent(s, 1:ncell(s), directions = m.neigh, sorted = TRUE, include = FALSE)

    out          <- data.frame(a)
    out$from_real <- values(s)[a[, 1]]
    out$to_real   <- values(s)[a[, 2]]

    identity <- cor.mtx[, max := max(value), by = id][value == max][value > threshold]
    identity[, variable := as.character(variable)]

    neighs <- data.table(
      from = identity$variable[match(out$from_real, identity$id)],
      to   = identity$variable[match(out$to_real,   identity$id)]
    )

    neighs_out[[as.character(n.neigh)]] <- (neighs[, .(table(from, to))])
  }

  neighs_out <- rbindlist(neighs_out, idcol = TRUE)
  neighs_out <- neighs_out[order(from, to)][,
    N := N - data.table::shift(N, fill = 0), by = c("from", "to")]
  neighs_out <- neighs_out[, frac := N / sum(N), by = list(from, .id)]

  return(neighs_out)
}
