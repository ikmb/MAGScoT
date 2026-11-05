test_that("completeness and contamination follow the DAS Tool formulas", {
  markers <- data.table(marker = sprintf("M%02d", 1:10), set = "toy")
  ## bin 1: markers 1-8 once, marker 1 on two more contigs (3 copies),
  ## marker 2 on one more contig (2 copies)
  cm <- data.table(contig = c(paste0("c", 1:8), "c9", "c10", "c11"),
                   marker = sprintf("M%02d", c(1:8, 1, 1, 2)))
  membership <- data.table(bin_id = 1L, contig = c(paste0("c", 1:11), "no_marker"))
  s <- score_bins(membership, cm, markers)

  expect_equal(s$uniqueSCGs, 8)
  expect_equal(s$multipleSCGs, 2)
  expect_equal(s$sumSCGs, 11)
  expect_equal(s$n_contigs, 12)
  expect_equal(s$Completeness, 0.8)
  expect_equal(s$Contamination, 2 / 8)
  expect_equal(s$additionalSCGs, 1)
  expect_equal(s$MultiContam, 1 / 10)
  expect_equal(s$score, 0.8 - 0.5 * 0.25 - 0.5 * 0.1)
})

test_that("bins are scored against every marker set and the best set is kept", {
  markers <- data.table(marker = c("A1", "A2", "B1", "B2", "B3", "B4"),
                        set = c("arc", "arc", "bac", "bac", "bac", "bac"))
  cm <- data.table(contig = c("x1", "x2", "x3"), marker = c("A1", "A2", "B1"))
  membership <- data.table(bin_id = 1L, contig = c("x1", "x2", "x3"))
  s <- score_bins(membership, cm, markers)
  expect_setequal(s$set, c("arc", "bac"))
  bins <- data.table(bin_id = 1L, bin = "b", bintool = "t", is_merged = FALSE)
  r <- rank_bins(s, bins)
  expect_equal(nrow(r), 1L)
  expect_equal(r$set, "arc")
  expect_equal(r$Completeness, 1)
})

test_that("scoring presets and explicit overrides", {
  expect_equal(score_weights(), c(a = 1, b = 0.5, c = 0.5))
  expect_equal(score_weights("metawrap"), c(a = 1, b = 5, c = 0.5))
  expect_equal(score_weights("metawrap", b = 2), c(a = 1, b = 2, c = 0.5))
  expect_equal(score_weights("dastool", a = 2, c = 0), c(a = 2, b = 0.5, c = 0))
  expect_error(score_weights("checkm"), "Unknown scoring preset")
})

test_that("the metawrap preset is stricter on contamination", {
  ## a: complete but 20% contaminated; b: 80% complete and clean
  g <- genome("g", 1:10)
  extra <- genome("x", 1:2)
  bins <- rbind(data.frame(bin = "a", contig = c(g$contig, extra$contig), set = "toolA"),
                data.frame(bin = "b", contig = g$contig[1:8], set = "toolB"))
  f <- toy_files(bins, rbind(g, extra))
  run <- function(...) quietly(magscot_refine(f$ctb, f$hmm, out = file.path(f$dir, "o"),
                                              profile = toy_markers(), skip_merge_bins = TRUE, ...))
  expect_equal(run()$refined$bin, "a")
  expect_equal(run(preset = "metawrap")$refined$bin, "b")
  expect_equal(run(preset = "metawrap", score_b = 0.5)$refined$bin, "a")
})
