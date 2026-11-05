## Helpers to build small synthetic datasets ----------------------------------

## A custom marker profile with `n` markers in one set called "toy".
toy_markers <- function(n = 10, dir = tempfile("magscot_markers_")) {
  dir.create(dir, showWarnings = FALSE)
  path <- file.path(dir, "toy_markers.tsv")
  writeLines(c("marker\tset", sprintf("M%02d.hmm\ttoy", seq_len(n))), path)
  path
}

## Write contig-to-bin rows and marker hits to temporary files.
## `bins`: data.frame(bin, contig, set); `hits`: data.frame(contig, marker)
toy_files <- function(bins, hits, dir = tempfile("magscot_toy_")) {
  dir.create(dir, showWarnings = FALSE)
  ctb <- file.path(dir, "ctb.tsv")
  hmm <- file.path(dir, "hits.hmm")
  write.table(bins[, c("bin", "contig", "set")], ctb, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  genes <- paste0(hits$contig, "_", ave(seq_along(hits$contig), hits$contig, FUN = seq_along))
  write.table(data.frame(genes, hits$marker, 1e-50), hmm, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  list(ctb = ctb, hmm = hmm, dir = dir)
}

## One contig per marker: contigs <prefix>_1..n carrying M01..Mn
genome <- function(prefix, markers) {
  data.frame(contig = paste0(prefix, "_", seq_along(markers)),
             marker = sprintf("M%02d", markers))
}

quietly <- function(expr) {
  old <- options(MAGScoT.verbose = FALSE)
  on.exit(options(old))
  suppressMessages(expr)
}

example_dir <- function() {
  d <- testthat::test_path("..", "..", "example")
  if (!file.exists(file.path(d, "example.hmm"))) testthat::skip("example data not available")
  d
}
