test_that("tab-separated, CRLF and whitespace-aligned inputs are read", {
  dir <- tempfile(); dir.create(dir)
  tab <- file.path(dir, "tab.tsv"); crlf <- file.path(dir, "crlf.tsv"); ws <- file.path(dir, "ws.tsv")
  writeLines(c("b1\tc1\tt1", "b1\tc2\tt1", "b2\tc3\tt2"), tab)
  writeBin(charToRaw("b1\tc1\tt1\r\nb1\tc2\tt1\r\nb2\tc3\tt2\r\n"), crlf)
  writeLines(c("b1   c1  t1", "b1   c2  t1", "b2   c3  t2"), ws)

  ref <- quietly(read_contig_to_bin(tab))
  expect_equal(quietly(read_contig_to_bin(crlf)), ref)
  expect_warning(x <- quietly(read_contig_to_bin(ws)), "not tab-separated")
  expect_equal(x, ref)
})

test_that("malformed contig-to-bin files give clear errors", {
  dir <- tempfile(); dir.create(dir)
  two <- file.path(dir, "two.tsv"); writeLines(c("b1\tc1", "b1\tc2"), two)
  expect_error(read_contig_to_bin(two), "must have 3 tab-separated columns")

  dup <- file.path(dir, "dup.tsv"); writeLines(c("b1\tc1\tt1", "b2\tc1\tt1", "b2\tc3\tt1", "b3\tc2\tt1"), dup)
  expect_error(read_contig_to_bin(dup), "more than one bin by the same binning tool")

  swapped <- file.path(dir, "swapped.tsv"); writeLines(c("c1\tb1\tt1", "c2\tb1\tt1"), swapped)
  expect_error(read_contig_to_bin(swapped), "more than one bin|column order")

  expect_error(read_contig_to_bin(file.path(dir, "missing.tsv")), "File not found")
})

test_that("marker hits keep the best hit per protein and map proteins to contigs", {
  f <- tempfile()
  writeLines(c("ctg_1_1\tM01\t1e-10", "ctg_1_1\tM02\t1e-50",
               "ctg_1_2\tM02\t1e-20", "ctg_2_1\tM03\t1e-5"), f)
  h <- hits_to_contigs(read_hmm(f))
  expect_equal(h[order(contig, marker)],
               data.table(contig = c("ctg_1", "ctg_2"), marker = c("M02", "M03")))
})

test_that("built-in and custom marker profiles load", {
  d <- read_markers("default")
  expect_true(all(c("gtdb_rel207_bac120", "gtdb_rel207_ar53") %in% d$set))
  expect_false(any(grepl("[.]hmm$", d$marker, ignore.case = TRUE)))
  expect_equal(nrow(read_markers("bac120")), 120L)
  expect_equal(nrow(read_markers(toy_markers(5))), 5L)
  expect_error(read_markers("nonsense"), "Unknown profile")
})

test_that("hmmsearch tblout is parsed using accessions", {
  f <- tempfile()
  writeLines(c(
    "# target name  accession  query name  accession  E-value",
    "ctg_1_1 - Ribosomal_L2 PF00380.20 1.2e-40 100.0 0.1 ...",
    "ctg_2_3 - TIGR00006 TIGR00006 3e-10 50.0 0.1 ...",
    "ctg_3_1 - SomeHMM - 1e-3 10.0 0.1 ..."), f)
  x <- parse_tblout(f)
  expect_equal(x$gene, c("ctg_1_1", "ctg_2_3", "ctg_3_1"))
  expect_equal(x$marker, c("PF00380.20", "TIGR00006", "SomeHMM"))
  expect_equal(x$evalue, c(1.2e-40, 3e-10, 1e-3))
})
