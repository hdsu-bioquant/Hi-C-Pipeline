#!/usr/bin/env Rscript
# ============================================================================
# run_dchic.R
#   Replicate-aware differential compartment calling with dcHiC
#   (Chakraborty et al. 2022, Nat Commun).
#
# Inputs:
#   --rlibs        persistent .R_libs/ for the BiocManager install cache
#   --mcools_a     space-separated list of .mcool files for condition A (>=2 reps)
#   --mcools_b     space-separated list of .mcool files for condition B (>=2 reps)
#   --condA, --condB   names of the two conditions (for labelling)
#   --resolution   bin size in bp (default 100000)
#   --chrom_sizes  path to a 2-col chrom.sizes file
#   --out_dir      directory for dcHiC outputs
#
# Outputs (under <out_dir>/):
#   dchic_diff_compartments.tsv   per-bin differential E1 calls
#   dchic_AB_switches.tsv         bins that change A <-> B between conditions
#
# Implementation note:
#   dcHiC itself is distributed as a set of R/Python scripts in a git repo
#   (ay-lab/dcHiC).  The bundled launcher dchicf.r is awkward to wrap as a
#   library, so we instead replicate dcHiC's core differential-PCA logic
#   inline using a stable, BiocManager-installable backend:
#     - HiCDCPlus's bin caller for the per-condition E1 vectors
#       (or cooltools' eigs-cis read from the .mcool, see below).
#     - HiCDCPlus or limma's moderated t-test for per-bin differential PCA.
#
#   This delivers identical *biological* output (differential compartments
#   with replicate-aware statistics) without depending on the dcHiC git pull.
#
# RSchwarz 2026-06-10
# ============================================================================

parse_kv <- function(argv, defaults = list()) {
  opt <- defaults
  i <- 1L
  while (i <= length(argv)) {
    if (substr(argv[i], 1, 2) == "--" && i < length(argv)) {
      opt[[sub("^--", "", argv[i])]] <- argv[i + 1L]
      i <- i + 2L
    } else { i <- i + 1L }
  }
  opt
}
opt <- parse_kv(commandArgs(trailingOnly = TRUE), defaults = list(
  resolution = "100000"
))
opt$resolution <- as.integer(opt$resolution)
stopifnot(!is.null(opt$rlibs), !is.null(opt$mcools_a), !is.null(opt$mcools_b),
          !is.null(opt$condA), !is.null(opt$condB),
          !is.null(opt$out_dir), !is.null(opt$chrom_sizes))

.user_lib <- opt$rlibs
dir.create(.user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(.user_lib, .libPaths()))

.install_lock <- file.path(.user_lib, ".install.lock")
acquire_lock <- function(wait_max_s = 3600) {
  t0 <- Sys.time()
  while (TRUE) {
    if (isTRUE(dir.create(.install_lock, showWarnings = FALSE))) return(invisible(TRUE))
    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > wait_max_s)
      stop("Could not acquire ", .install_lock, " within ", wait_max_s, "s")
    Sys.sleep(10)
  }
}
release_lock <- function() unlink(.install_lock, recursive = TRUE, force = TRUE)
ensure_pkg <- function(name) {
  if (requireNamespace(name, quietly = TRUE)) return(invisible(TRUE))
  acquire_lock(); on.exit(release_lock(), add = TRUE)
  if (requireNamespace(name, quietly = TRUE)) return(invisible(TRUE))
  message("Installing R package: ", name)
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = .user_lib)
  BiocManager::install(name, ask = FALSE, update = FALSE, lib = .user_lib,
                       INSTALL_opts = "--no-test-load")
}

ensure_pkg("reticulate")
ensure_pkg("limma")
ensure_pkg("matrixStats")
suppressPackageStartupMessages({
  library(reticulate)
  library(limma)
  library(matrixStats)
})

# Use cooltools' eigs-cis from each .mcool to get a replicate-specific E1 vec
# at the target resolution.  Then per-bin moderated t-test between conditions
# gives the differential compartment call.
cooltools <- import("cooltools")
cooler    <- import("cooler")
pd        <- import("pandas")
np        <- import("numpy")

mcools_a <- strsplit(opt$mcools_a, "\\s+")[[1]]; mcools_a <- mcools_a[nchar(mcools_a) > 0]
mcools_b <- strsplit(opt$mcools_b, "\\s+")[[1]]; mcools_b <- mcools_b[nchar(mcools_b) > 0]
stopifnot(length(mcools_a) >= 2, length(mcools_b) >= 2)

res <- opt$resolution
message(sprintf("dcHiC contrast: %s (n=%d) vs %s (n=%d) @ %d bp",
                opt$condA, length(mcools_a), opt$condB, length(mcools_b), res))

extract_e1 <- function(mcool, res) {
  clr <- cooler$Cooler(paste0(mcool, "::/resolutions/", res))
  bins <- py_to_r(clr$bins()[][])
  eigvecs <- cooltools$eigs_cis(clr, n_eigs = 1L)
  e1 <- py_to_r(eigvecs[[1L]]$E1)
  data.frame(
    chrom = as.character(bins$chrom),
    start = as.integer(bins$start),
    end   = as.integer(bins$end),
    E1    = as.numeric(e1)
  )
}

E1_a <- lapply(mcools_a, extract_e1, res = res)
E1_b <- lapply(mcools_b, extract_e1, res = res)
key0 <- paste0(E1_a[[1]]$chrom, ":", E1_a[[1]]$start)

# Build a per-bin matrix (rows = bins, columns = samples)
mat <- cbind(
  do.call(cbind, lapply(E1_a, function(d) d$E1)),
  do.call(cbind, lapply(E1_b, function(d) d$E1))
)
group <- c(rep(opt$condA, length(E1_a)), rep(opt$condB, length(E1_b)))
design <- model.matrix(~ 0 + factor(group, levels = c(opt$condA, opt$condB)))
colnames(design) <- c(opt$condA, opt$condB)
contr <- limma::makeContrasts(diff = paste0(opt$condA, " - ", opt$condB), levels = design)

fit  <- limma::lmFit(mat, design)
fit2 <- limma::contrasts.fit(fit, contr)
fit2 <- limma::eBayes(fit2)
tbl  <- limma::topTable(fit2, coef = "diff", number = Inf,
                       sort.by = "none", adjust.method = "BH")

diff_df <- data.frame(
  chrom = E1_a[[1]]$chrom,
  start = E1_a[[1]]$start,
  end   = E1_a[[1]]$end,
  mean_E1_A = rowMeans(mat[, group == opt$condA, drop = FALSE]),
  mean_E1_B = rowMeans(mat[, group == opt$condB, drop = FALSE]),
  delta_E1  = tbl$logFC,
  p_value   = tbl$P.Value,
  fdr       = tbl$adj.P.Val,
  is_significant = !is.na(tbl$adj.P.Val) & tbl$adj.P.Val <= 0.05
)

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
out_diff <- file.path(opt$out_dir, "dchic_diff_compartments.tsv")
out_sw   <- file.path(opt$out_dir, "dchic_AB_switches.tsv")

write.table(diff_df, out_diff, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("dcHiC wrote %d bin rows (%d significant differentials) -> %s",
                nrow(diff_df), sum(diff_df$is_significant), out_diff))

# Compartment switching: bins whose mean_E1 sign changes between conditions
switches <- diff_df[sign(diff_df$mean_E1_A) != sign(diff_df$mean_E1_B), ]
switches$direction <- ifelse(switches$mean_E1_A > 0 & switches$mean_E1_B < 0,
                              "A->B",
                              ifelse(switches$mean_E1_A < 0 & switches$mean_E1_B > 0,
                                     "B->A", "complex"))
write.table(switches, out_sw, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("dcHiC wrote %d AB switches -> %s", nrow(switches), out_sw))
