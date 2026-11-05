run_toy <- function(bins, hits, n_markers = 10, ...) {
  f <- toy_files(bins, hits)
  quietly(magscot_refine(f$ctb, f$hmm, out = file.path(f$dir, "out"),
                         profile = toy_markers(n_markers), ...))
}

test_that("bin names shared between binning tools are kept apart (issue #4)", {
  g1 <- genome("g1", 1:10)
  g2 <- genome("g2", 1:10)
  bins <- rbind(data.frame(bin = "bin_1", contig = g1$contig, set = "toolA"),
                data.frame(bin = "bin_1", contig = g2$contig, set = "toolB"))
  res <- run_toy(bins, rbind(g1, g2))
  expect_equal(nrow(res$refined), 2L)
  expect_true(all(res$refined$Completeness == 1 & res$refined$Contamination == 0))
  expect_setequal(res$refined$bintool, c("toolA", "toolB"))
})

test_that("refined bins never fall below the score cutoff", {
  ## toolA: two clean bins; plus a contaminated bin from toolA that would have
  ## been accepted by the same-tool shortcut in v1.1
  g1 <- genome("g1", 1:10)
  g2 <- genome("g2", 1:9)
  g3 <- rbind(genome("g3", 1:6), genome("g3x", 1:4))   # completeness .6, contamination .67
  bins <- rbind(data.frame(bin = "a1", contig = g1$contig, set = "toolA"),
                data.frame(bin = "a2", contig = g2$contig, set = "toolA"),
                data.frame(bin = "a3", contig = g3$contig, set = "toolA"),
                data.frame(bin = "b1", contig = g1$contig[1:5], set = "toolB"))
  res <- run_toy(bins, rbind(g1, g2, g3))
  expect_true(all(res$refined$score >= 0.5))
  expect_setequal(res$refined$bin, c("a1", "a2"))
})

test_that("a bin consisting of a single contig can be selected", {
  hits <- data.frame(contig = "circ", marker = sprintf("M%02d", 1:10))
  bins <- rbind(data.frame(bin = "x", contig = "circ", set = "toolA"),
                data.frame(bin = "y", contig = "circ", set = "toolB"))
  res <- run_toy(bins, hits, skip_merge_bins = TRUE)
  expect_equal(nrow(res$refined), 1L)
  expect_equal(res$assignments$contig, "circ")
})

test_that("a single binning tool is refined without merging", {
  g1 <- genome("g1", 1:10)
  g2 <- genome("g2", 1:3)
  bins <- rbind(data.frame(bin = "a1", contig = g1$contig, set = "toolA"),
                data.frame(bin = "a2", contig = g2$contig, set = "toolA"))
  res <- run_toy(bins, rbind(g1, g2))
  expect_equal(res$refined$bin, "a1")
  expect_equal(nrow(res$scores), 2L)
})

test_that("score_only writes scores and nothing else", {
  g1 <- genome("g1", 1:10)
  bins <- rbind(data.frame(bin = "a1", contig = g1$contig, set = "toolA"),
                data.frame(bin = "b1", contig = g1$contig, set = "toolB"))
  f <- toy_files(bins, g1)
  out <- file.path(f$dir, "out")
  quietly(magscot_refine(f$ctb, f$hmm, out = out, profile = toy_markers(), score_only = TRUE))
  expect_true(file.exists(paste0(out, ".scores.out")))
  expect_false(file.exists(paste0(out, ".refined.out")))
})

test_that("refined bin names use the basename of the output prefix", {
  g1 <- genome("g1", 1:10)
  bins <- rbind(data.frame(bin = "a1", contig = g1$contig, set = "toolA"),
                data.frame(bin = "b1", contig = g1$contig, set = "toolB"))
  f <- toy_files(bins, g1)
  res <- quietly(magscot_refine(f$ctb, f$hmm, out = file.path(f$dir, "sub", "S1"),
                                profile = toy_markers()))
  expect_equal(res$refined$refined_id, "S1_cleanbin_000001")
  expect_true(file.exists(file.path(f$dir, "sub", "S1.refined.contig_to_bin.out")))
})

test_that("markers that do not match the profile are reported", {
  g1 <- genome("g1", 1:10)
  g1$marker <- paste0("X", g1$marker)
  bins <- data.frame(bin = "a1", contig = g1$contig, set = "toolA")
  f <- toy_files(bins, g1)
  expect_error(quietly(magscot_refine(f$ctb, f$hmm, out = file.path(f$dir, "o"),
                                      profile = toy_markers())),
               "None of the markers")
})
