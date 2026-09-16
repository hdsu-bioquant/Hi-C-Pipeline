#!/usr/bin/env Rscript
# ============================================================================
# run_multihiccompare.R
#   Differential bin-pair analysis with multiHiCcompare (Stansfield 2018).
#
# Inputs:
#   --rlibs      persistent .R_libs/ for the BiocManager install cache
#   --mcools_a   space-separated list of .mcool files for condition A (>=2 reps)
#   --mcools_b   space-separated list of .mcool files for condition B (>=2 reps)
#   --resolution bin size in bp (default 100000)
#   --fdr        BH-FDR cutoff for flagging is_significant (default 0.05)
#   --condA      name of condition A (for output labelling)
#   --condB      name of condition B
#   --out_tsv    path of the output differential bin-pair TSV
# ============================================================================

# ---- persistent library + install-if-missing -------------------------------
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
  resolution = "100000", fdr = "0.05"
))
opt$resolution <- as.integer(opt$resolution)
opt$fdr        <- as.numeric(opt$fdr)
stopifnot(!is.null(opt$rlibs), !is.null(opt$mcools_a), !is.null(opt$mcools_b),
          !is.null(opt$out_tsv))

.user_lib <- opt$rlibs
dir.create(.user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(.user_lib, .libPaths()))

# file-system lock (mirrors ChIPseeker/DiffBind pattern)
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

ensure_pkg("multiHiCcompare")
ensure_pkg("HiCcompare")
suppressPackageStartupMessages({
  library(multiHiCcompare)
  library(HiCcompare)
})

# ---- Build the per-condition replicate lists -------------------------------
mcools_a <- strsplit(opt$mcools_a, "\\s+")[[1]]
mcools_b <- strsplit(opt$mcools_b, "\\s+")[[1]]
mcools_a <- mcools_a[nchar(mcools_a) > 0]
mcools_b <- mcools_b[nchar(mcools_b) > 0]

stopifnot(length(mcools_a) >= 2, length(mcools_b) >= 2)
res <- opt$resolution
message(sprintf("multiHiCcompare contrast: %s (n=%d) vs %s (n=%d) @ %d bp",
                opt$condA, length(mcools_a), opt$condB, length(mcools_b), res))

# multiHiCcompare expects a list of data frames with columns
# (chr1, start1, chr2, start2, IF) — bin-pair "interaction frequency" tables.
# We extract those from each .mcool at the target resolution using cooler via
# rhdf5/HDF5Array, but the simpler/safer cross-version route is to ask cooler
# to dump a per-resolution flat text file via reticulate.
ensure_pkg("reticulate")
library(reticulate)
py_install <- function(pkg) {
  if (!py_module_available(pkg))
    py_install(pkg, envname = NULL, method = "auto", pip = TRUE)
}
py_install("cooler"); py_install("pandas")

cooler <- import("cooler")
pd     <- import("pandas")

extract_bedpe <- function(mcool, res) {
  clr <- cooler$Cooler(paste0(mcool, "::/resolutions/", res))
  mat <- clr$matrix(balance = TRUE, as_pixels = TRUE, join = TRUE)[][]
  df <- py_to_r(mat)
  # Keep cis only and the columns multiHiCcompare wants
  cis <- df[df$chrom1 == df$chrom2, c("chrom1", "start1", "start2", "balanced")]
  colnames(cis) <- c("chr", "region1", "region2", "IF")
  cis$IF[is.na(cis$IF)] <- 0
  cis
}

tabs_a <- lapply(mcools_a, extract_bedpe, res = res)
tabs_b <- lapply(mcools_b, extract_bedpe, res = res)

# Build hicexp object (multiHiCcompare canonical entry)
hicexp <- multiHiCcompare::make_hicexp(
  data_list = c(tabs_a, tabs_b),
  groups   = c(rep(1, length(tabs_a)), rep(2, length(tabs_b))),
  remove_regions = NULL
)
hicexp <- multiHiCcompare::cyclic_loess(hicexp)
hicexp <- multiHiCcompare::hic_exactTest(hicexp)

res_tbl <- multiHiCcompare::results(hicexp)
res_tbl$is_significant <- !is.na(res_tbl$p.adj) & res_tbl$p.adj <= opt$fdr

write.table(res_tbl, opt$out_tsv, sep = "\t",
            quote = FALSE, row.names = FALSE)
message(sprintf("multiHiCcompare wrote %d bin-pair rows (%d significant) -> %s",
                nrow(res_tbl), sum(res_tbl$is_significant), opt$out_tsv))
