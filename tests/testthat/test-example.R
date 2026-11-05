## Regression test on the bundled example data against MAGScoT v1.1 output.

test_that("example: scores match v1.1 and refined bins pass the cutoff", {
  ex <- example_dir()
  out <- file.path(tempfile("magscot_ex_"), "example")
  res <- quietly(magscot_refine(file.path(ex, "example.contigs_to_bin.tsv"),
                                file.path(ex, "example.hmm"), out = out, n_iterations = 1))
  v11 <- fread(test_path("fixtures", "example_v1.1_scores.tsv"))

  orig_new <- res$scores[!grepl("^MAGScoT_", bintool)]
  orig_old <- v11[!grepl("^MAGScoT_", bintool)]
  m <- merge(orig_old, orig_new, by = c("bin", "bintool"), suffixes = c(".old", ".new"))
  expect_equal(nrow(m), nrow(orig_old))
  expect_equal(m$score.new, m$score.old, tolerance = 1e-12)
  expect_equal(m$set.new, m$set.old)

  ## hybrid bins of the first merging round are the same as in v1.1
  expect_equal(sort(res$scores[bintool == "MAGScoT_1", score]),
               sort(v11[bintool == "MAGScoT_1", score]), tolerance = 1e-12)

  expect_equal(nrow(res$refined), 15L)
  expect_true(all(res$refined$score >= 0.5))
  expect_false(anyDuplicated(res$assignments$contig) > 0)
})

test_that("example: bins are written as FASTA files", {
  ex <- example_dir()
  dir <- tempfile("magscot_fa_")
  a <- data.table(binnew = c("b1", "b1", "b2"),
                  contig = c("example_contig_13", "example_contig_20282", "example_contig_30422"))
  quietly(write_bin_fastas(a, file.path(ex, "example.contigs.fasta.gz"), dir, chunk_lines = 50))
  hdr <- function(f) sub("^>([^ ]+).*", "\\1", grep("^>", readLines(file.path(dir, f)), value = TRUE))
  expect_setequal(hdr("b1.fasta"), c("example_contig_13", "example_contig_20282"))
  expect_equal(hdr("b2.fasta"), "example_contig_30422")
})
