## Main workflow -------------------------------------------------------------------

## Score, merge and refine metagenomic bins
##
## Scores all bins with GTDB r207 marker genes, optionally creates hybrid bins
## from overlapping bins of different binning tools, and iteratively selects
## the best non-redundant set of bins. Results are written to files starting
## with `out`.
##
## @param input Contig-to-bin file: bin, contig, binning tool (tab-separated, no header).
## @param hmm Marker file: protein id, marker, e-value (tab-separated, no header).
##   Protein ids must be `<contig>_<n>` (prodigal convention).
## @param proteins Protein FASTA from prodigal or pyrodigal (may be gzipped), as
##   an alternative to `hmm`: the marker genes are then searched with hmmsearch
##   and the hits written to `<out>.hmm`.
## @param threads Number of hmmsearch CPUs (with `proteins`).
## @param out Output prefix; writes `<out>.scores.out`, `<out>.refined.out` and
##   `<out>.refined.contig_to_bin.out`.
## @param profile Marker profile: `"default"` (bac120 + ar53), `"bac120"`,
##   `"ar53"`, or the path to a TSV with columns `marker` and `set`.
## @param bin_separator Infix of refined bin names (`<out>_<infix>_000001`).
## @param preset Scoring weight preset: `"dastool"` (a = 1, b = 0.5, c = 0.5;
##   default) or `"metawrap"` (a = 1, b = 5, c = 0.5).
## @param score_a,score_b,score_c Scoring weights for completeness,
##   contamination and surplus marker copies; override the preset if given.
## @param threshold Minimum score (and completeness) of refined bins.
## @param max_contamination Maximum contamination of bins considered for
##   refinement, as a fraction (0.1 = 10%) of MAGScoT's `Contamination` column
##   (markers present in more than one copy / unique markers). The default 1
##   applies no filter.
## @param score_only Only score the input bins; no merging or refinement.
## @param skip_merge_bins Skip bin merging.
## @param min_markers Minimum number of unique markers of a seed bin, and of
##   shared markers, for bin merging.
## @param min_sharing Minimum fraction of a seed bin's markers that must be
##   found in the partner bin for merging.
## @param n_iterations Number of bin-merging iterations.
## @param max_marker_hits Maximum total number of marker hits of a seed bin.
## @param contigs Optional assembly FASTA (may be gzipped); if given, refined
##   bins are written as FASTA files to `bins_dir`.
## @param bins_dir Directory for the refined bin FASTA files.
## @return Invisibly, a list with `scores` (all candidate bins), `refined`
##   (refined bin statistics) and `assignments` (contig to refined bin).

magscot_refine <- function(input, hmm = NULL, out = "MAGScoT", profile = "default",
                           proteins = NULL, threads = 1L,
                           bin_separator = "cleanbin",
                           preset = "dastool",
                           score_a = NULL, score_b = NULL, score_c = NULL,
                           threshold = 0.5, max_contamination = 1,
                           score_only = FALSE, skip_merge_bins = FALSE,
                           min_markers = 25, min_sharing = 0.8, n_iterations = 2,
                           max_marker_hits = 150,
                           contigs = NULL, bins_dir = paste0(out, "_bins")) {
  out_dir <- dirname(out)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  if (is.null(hmm) == is.null(proteins)) {
    stop("Give either a marker file (--hmm) or a protein FASTA (--proteins).", call. = FALSE)
  }
  w <- score_weights(preset, score_a, score_b, score_c)
  log_msg(sprintf("Scoring weights: a = %g, b = %g, c = %g (preset '%s'%s).",
                  w[["a"]], w[["b"]], w[["c"]], preset,
                  if (any(!vapply(list(score_a, score_b, score_c), is.null, logical(1))))
                    ", with overrides" else ""))
  ctb <- read_contig_to_bin(input)
  markers <- read_markers(profile)
  if (is.null(proteins)) {
    hits <- read_hmm(hmm)
  } else {
    tmp <- tempfile("magscot_", tmpdir = out_dir)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    hits <- search_markers(proteins, threads, tmp)
    hmm <- paste0(out, ".hmm")
    fwrite(hits, hmm, sep = "\t", col.names = FALSE, quote = FALSE)
    log_msg(sprintf("Wrote %d marker hits to: %s (reuse with --hmm)", nrow(hits), hmm))
  }
  contig_markers <- hits_to_contigs(hits)

  matched <- mean(unique(contig_markers$marker) %in% markers$marker)
  if (matched == 0) {
    stop("None of the markers in '", hmm, "' are part of the '", profile,
         "' marker profile. Were the GTDB r207 HMMs (and their accessions) used?",
         call. = FALSE)
  }
  in_bins <- mean(unique(contig_markers$contig) %in% ctb$contig)
  if (in_bins == 0) {
    stop("No marker-carrying contig from '", hmm, "' occurs in '", input,
         "'. Do the protein IDs follow the prodigal convention <contig>_<n>?",
         call. = FALSE)
  }

  ## internal integer ids; bins are identified by binning tool + name
  bins <- unique(ctb[, .(bin, bintool = set)])
  bins[, `:=`(bin_id = seq_len(.N), is_merged = FALSE)]
  setcolorder(bins, c("bin_id", "bin", "bintool", "is_merged"))
  membership <- bins[ctb, on = c(bin = "bin", bintool = "set"), .(bin_id, contig)]

  n_sets <- uniqueN(bins$bintool)
  log_msg(sprintf("Read %d bins from %d binning tool(s) covering %d contigs; %d marker hits.",
                  nrow(bins), n_sets, uniqueN(ctb$contig), nrow(contig_markers)))

  if (n_sets == 1L && !score_only) {
    log_msg("Only one binning tool provided: skipping bin merging.")
    skip_merge_bins <- TRUE
  }
  n_iterations <- as.integer(n_iterations)
  if (!skip_merge_bins && !score_only && n_iterations >= n_sets) {
    log_msg("Bin-merging iterations must be fewer than the number of binning tools; using ",
            n_sets - 1L, ".")
    n_iterations <- n_sets - 1L
  }

  if (!skip_merge_bins && !score_only && n_iterations > 0L) {
    merged <- merge_bins(membership, bins, contig_markers, n_iter = n_iterations,
                         min_markers = min_markers, min_sharing = min_sharing,
                         max_marker_hits = max_marker_hits)
    membership <- merged$membership
    bins <- merged$bins
  }

  res <- refine_bins(membership, bins, contig_markers, markers,
                     a = w[["a"]], b = w[["b"]], c = w[["c"]],
                     cutoff = threshold, max_contamination = max_contamination,
                     score_only = score_only)

  scores_file <- paste0(out, ".scores.out")
  scores <- format_scores(res$initial_scores)
  write_tsv(scores, scores_file)
  log_msg("Scores for all candidate bins written to: ", scores_file)
  if (score_only) return(invisible(list(scores = scores)))

  if (!nrow(res$selected)) {
    log_msg("No bins with score >= ", threshold, " were found.")
    return(invisible(list(scores = scores, refined = NULL, assignments = NULL)))
  }

  selected <- res$selected
  selected[, refined_id := sprintf("%s_%s_%06d", basename(out), bin_separator, seq_len(.N))]
  stats <- cbind(format_scores(selected), selected[, .(refined_id)])
  assignments <- selected[, .(bin_id, binnew = refined_id)][res$assignments, on = "bin_id",
                                                             .(binnew, contig)]

  stats_file <- paste0(out, ".refined.out")
  ctb_file <- paste0(out, ".refined.contig_to_bin.out")
  write_tsv(stats, stats_file)
  write_tsv(assignments, ctb_file)
  log_msg(sprintf("Refinement yielded %d bins with score >= %s.", nrow(stats), threshold))
  log_msg("Refined bin statistics written to: ", stats_file)
  log_msg("Contig-to-refined-bin mapping written to: ", ctb_file)

  if (!is.null(contigs)) {
    write_bin_fastas(assignments, contigs, bins_dir)
  }
  invisible(list(scores = scores, refined = stats, assignments = assignments))
}
