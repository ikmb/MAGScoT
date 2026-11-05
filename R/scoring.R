## Bin scoring ------------------------------------------------------------------
##
## For every bin and every marker set (e.g. bac120, ar53):
##   Completeness        = unique markers / markers in set
##   Contamination       = markers present in >1 copy / unique markers
##   MultiContam         = surplus copies beyond the second / markers in set
##   Contamination_Score = b * Contamination + c * MultiContam
##   score               = a * Completeness - Contamination_Score
## The weights come from a preset (see SCORE_PRESETS) and can be overridden
## individually.

## "dastool": DAS Tool defaults (MAGScoT default since v1).
## "metawrap": contamination weighted 10x more strongly, following MetaWRAP's
##   completeness - 5 x contamination. MetaWRAP uses CheckM contamination, so
##   this approximates its stringency rather than reproducing it exactly.
SCORE_PRESETS <- list(
  dastool  = c(a = 1, b = 0.5, c = 0.5),
  metawrap = c(a = 1, b = 5,   c = 0.5)
)

## Weights from `preset`, with any non-NULL a/b/c taking precedence.
score_weights <- function(preset = "dastool", a = NULL, b = NULL, c = NULL) {
  if (!preset %in% names(SCORE_PRESETS)) {
    stop("Unknown scoring preset '", preset, "'. Use one of: ",
         paste(names(SCORE_PRESETS), collapse = ", "), ".", call. = FALSE)
  }
  w <- SCORE_PRESETS[[preset]]
  if (!is.null(a)) w[["a"]] <- a
  if (!is.null(b)) w[["b"]] <- b
  if (!is.null(c)) w[["c"]] <- c
  w
}

SCORE_COLUMNS <- c("bin_id", "set", "uniqueSCGs", "multipleSCGs", "sumSCGs",
                   "n_markers", "n_contigs", "Completeness", "additionalSCGs",
                   "Contamination", "MultiContam", "Contamination_Score", "score")

## Score bins against each marker set.
##
## @param membership data.table(bin_id, contig)
## @param contig_markers data.table(contig, marker), unique rows
## @param markers data.table(marker, set)
## @return one row per (bin_id, marker set) that has at least one marker hit
score_bins <- function(membership, contig_markers, markers, a = 1, b = 0.5, c = 0.5) {
  n_contigs <- membership[, .(n_contigs = .N), by = bin_id]
  hits <- contig_markers[membership, on = "contig", nomatch = NULL,
                         allow.cartesian = TRUE, .(bin_id, marker)]
  counts <- hits[, .(count = .N), by = .(bin_id, marker)]
  counts <- markers[counts, on = "marker", nomatch = NULL, allow.cartesian = TRUE]

  res <- counts[, .(uniqueSCGs = .N,
                    multipleSCGs = sum(count > 1L),
                    sumSCGs = sum(count)),
                by = .(bin_id, set)]
  res <- markers[, .(n_markers = .N), by = set][res, on = "set"]
  res <- n_contigs[res, on = "bin_id"]
  res[, `:=`(
    Completeness   = uniqueSCGs / n_markers,
    additionalSCGs = sumSCGs - uniqueSCGs - multipleSCGs,
    Contamination  = multipleSCGs / uniqueSCGs
  )]
  res[, MultiContam := additionalSCGs / n_markers]
  res[, Contamination_Score := b * Contamination + c * MultiContam]
  res[, score := a * Completeness - Contamination_Score]
  setcolorder(res, SCORE_COLUMNS)
  res[]
}

## Keep the best-scoring marker set per bin and rank bins. Ties are broken by
## fewer contigs, then original over merged bins, then name.
rank_bins <- function(scores, bins) {
  x <- bins[scores, on = "bin_id"]
  setorderv(x, c("score", "n_contigs", "is_merged", "bin", "set"),
            order = c(-1L, 1L, 1L, 1L, 1L))
  x[!duplicated(bin_id)]
}

## Format a ranked score table for output (columns as in MAGScoT v1).
format_scores <- function(ranked) {
  out <- ranked[, c("bin", setdiff(SCORE_COLUMNS, "bin_id"), "bintool"), with = FALSE]
  out[]
}
