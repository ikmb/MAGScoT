## These tests need hmmsearch on the PATH (skipped otherwise).

test_that("a protein FASTA gives the same refinement as precomputed marker hits", {
  skip_if(!nzchar(Sys.which("hmmsearch")), "hmmsearch not available")
  ex <- example_dir()
  dir <- tempfile("magscot_prot_")
  from_hmm <- quietly(magscot_refine(file.path(ex, "example.contigs_to_bin.tsv"),
                                     hmm = file.path(ex, "example.hmm"),
                                     out = file.path(dir, "hmm"), n_iterations = 1))
  from_faa <- quietly(magscot_refine(file.path(ex, "example.contigs_to_bin.tsv"),
                                     proteins = file.path(ex, "example.prodigal.faa.gz"),
                                     out = file.path(dir, "faa"), n_iterations = 1, threads = 2))
  expect_true(file.exists(file.path(dir, "faa.hmm")))
  expect_equal(from_faa$scores, from_hmm$scores)
  expect_equal(from_faa$refined[, !"refined_id"], from_hmm$refined[, !"refined_id"])
})

test_that("annotate accepts a protein FASTA and writes the marker file", {
  skip_if(!nzchar(Sys.which("hmmsearch")), "hmmsearch not available")
  ex <- example_dir()
  out <- file.path(tempfile("magscot_ann_"), "ex")
  quietly(magscot_annotate(proteins = file.path(ex, "example.prodigal.faa.gz"), out = out))
  hits <- read_hmm(paste0(out, ".hmm"))
  ref <- read_hmm(file.path(ex, "example.hmm"))
  expect_equal(nrow(hits), nrow(ref))
  expect_setequal(paste(hits$gene, hits$marker), paste(ref$gene, ref$marker))
  expect_false(file.exists(paste0(out, ".prodigal.faa")))
})

test_that("exactly one marker source must be given", {
  expect_error(magscot_refine("x", out = tempfile()), "either a marker file")
  expect_error(magscot_refine("x", hmm = "a", proteins = "b", out = tempfile()), "either a marker file")
  expect_error(magscot_annotate(out = tempfile()), "either an assembly")
})
