## Iterative refinement ------------------------------------------------------------
##
## Repeatedly pick the highest-scoring bin, remove its contigs from all other
## candidate bins, re-score the affected bins, and continue until no bin
## reaches the score cutoff.
##
## Speed-up: bins of the same binning tool never share contigs, so as long as
## the next-best bins come from the same tool as the current winner they can be
## accepted without re-scoring (this does not apply to hybrid bins, which can
## overlap each other).

## @param membership data.table(bin_id, contig) of all candidate bins
## @param bins data.table(bin_id, bin, bintool, is_merged)
## @return list(initial_scores, selected, assignments)
refine_bins <- function(membership, bins, contig_markers, markers,
                        a = 1, b = 0.5, c = 0.5, cutoff = 0.5, max_contamination = 1,
                        score_only = FALSE) {
  remain <- membership
  scores <- NULL
  recalc <- unique(remain$bin_id)
  initial_scores <- NULL
  selected <- list()
  assigned <- list()
  i <- 1L

  while (nrow(remain) > 0L) {
    if (length(recalc)) {
      fresh <- score_bins(remain[bin_id %in% recalc], contig_markers, markers, a, b, c)
      if (!is.null(scores)) {
        scores <- scores[!bin_id %in% recalc & bin_id %in% remain$bin_id]
      }
      scores <- rbind(scores, fresh)
    } else {
      scores <- scores[bin_id %in% remain$bin_id]
    }
    ranked <- rank_bins(scores, bins)

    if (i == 1L) {
      initial_scores <- ranked
      if (score_only) break
    }
    if (!nrow(ranked) || ranked$score[1] < cutoff) break

    ## bins below the completeness cutoff can never recover; high-contamination
    ## bins are dropped as requested by --max_contamination
    cand <- ranked[Completeness >= cutoff & Contamination <= max_contamination]

    take <- 0L
    if (nrow(cand)) {
      tool <- cand$bintool[1]
      while (take < nrow(cand) && cand$bintool[take + 1L] == tool &&
             cand$score[take + 1L] >= cutoff) {
        take <- take + 1L
        if (cand$is_merged[take]) break
      }
    }

    if (take > 0L) {
      winners <- cand[seq_len(take)]
      win_members <- remain[bin_id %in% winners$bin_id]
      selected[[length(selected) + 1L]] <- winners
      assigned[[length(assigned) + 1L]] <- win_members
      recalc <- unique(remain[contig %in% win_members$contig, bin_id])
      recalc <- setdiff(recalc, winners$bin_id)
      log_msg(sprintf("Refinement iteration %d: selected %d bin(s) from '%s' (best: %s, score %.3f)",
                      i, take, winners$bintool[1], winners$bin[1], winners$score[1]))
      cand <- cand[-seq_len(take)]
      remain <- remain[bin_id %in% cand$bin_id & !contig %in% win_members$contig]
    } else {
      recalc <- integer(0)
      remain <- remain[bin_id %in% cand$bin_id]
    }
    i <- i + 1L
  }

  list(initial_scores = initial_scores,
       selected = rbindlist(selected),
       assignments = rbindlist(assigned))
}
