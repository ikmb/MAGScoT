## A case where a hybrid bin from iteration 1 is needed to absorb a third bin:
##   A = a1-a6 + b1-b24, B = b1-b100      -> A and B merge (A shares 80% of its markers)
##   C = a1-a6 + b1-b3 + b25-b100 + c1-c30 -> C overlaps neither A nor B enough,
##   but the hybrid A+B shares 85/106 of its markers with C -> A+B+C in iteration 2
make_overlap <- function() {
  g <- genome("g", 1:136)               # a = 1-6, b = 7-106, c = 107-136
  ct <- g$contig
  bins <- rbind(
    data.frame(bin = "A", contig = ct[1:30], set = "toolA"),
    data.frame(bin = "B", contig = ct[7:106], set = "toolB"),
    data.frame(bin = "C", contig = ct[c(1:9, 31:136)], set = "toolC"))
  list(bins = bins, hits = g)
}

prep <- function(bins, hits) {
  ctb <- as.data.table(bins)
  b <- unique(ctb[, .(bin, bintool = set)])
  b[, `:=`(bin_id = seq_len(.N), is_merged = FALSE)]
  setcolorder(b, c("bin_id", "bin", "bintool", "is_merged"))
  m <- b[ctb, on = c(bin = "bin", bintool = "set"), .(bin_id, contig)]
  list(membership = m, bins = b, cm = as.data.table(hits))
}

test_that("hybrid bins of later iterations include the seed's contigs", {
  d <- make_overlap()
  p <- prep(d$bins, d$hits)
  res <- quietly(merge_bins(p$membership, p$bins, p$cm, n_iter = 2, min_markers = 10))
  expect_equal(nrow(res$bins[bintool == "MAGScoT_1"]), 1L)
  iter2 <- res$bins[bintool == "MAGScoT_2"]
  expect_gt(nrow(iter2), 0L)
  ## the union of all three bins covers the whole genome
  sizes <- res$membership[bin_id %in% iter2$bin_id, .N, by = bin_id]
  expect_true(any(sizes$N == 136L))
  ## no two candidate bins have identical contig sets
  keys <- contig_set_keys(res$membership)
  expect_false(anyDuplicated(keys$key) > 0)
})

test_that("an iteration without overlaps stops cleanly (issue #9)", {
  g1 <- genome("g1", 1:30)
  g2 <- genome("g2", 1:30)
  bins <- rbind(data.frame(bin = "A", contig = g1$contig, set = "toolA"),
                data.frame(bin = "B", contig = g2$contig, set = "toolB"),
                data.frame(bin = "C", contig = g1$contig[1:2], set = "toolC"))
  p <- prep(bins, rbind(g1, g2))
  res <- quietly(merge_bins(p$membership, p$bins, p$cm, n_iter = 2, min_markers = 10))
  expect_equal(nrow(res$bins), 3L)

  f <- toy_files(bins, rbind(g1, g2))
  expect_no_error(quietly(magscot_refine(f$ctb, f$hmm, out = file.path(f$dir, "o"),
                                         profile = toy_markers(30), min_markers = 10)))
})
