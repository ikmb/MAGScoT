## Bin merging --------------------------------------------------------------------
##
## Bins from different binning tools that share most of their marker-carrying
## contigs are combined into hybrid bins. These are added to the candidate pool
## (as binning tool "MAGScoT_<iteration>") if their contig set is new.
##
## Iteration 1 uses the original bins as seeds. Each later iteration uses the
## hybrid bins created in the previous iteration as seeds and merges them with
## original bins again.

## Canonical key of a bin's contig set, used to drop duplicate bins.
contig_set_keys <- function(membership) {
  membership[order(contig), .(key = paste(contig, collapse = "\r")), by = bin_id]
}

## Create hybrid bins from overlapping bins of different binning tools.
##
## @param membership data.table(bin_id, contig) of the original bins
## @param bins data.table(bin_id, bin, bintool, is_merged)
## @param contig_markers data.table(contig, marker)
## @return list(membership, bins) including the original and new hybrid bins
merge_bins <- function(membership, bins, contig_markers, n_iter = 2L,
                       min_markers = 25, min_sharing = 0.8, max_marker_hits = 150) {
  all_members <- membership
  all_bins <- bins
  known_keys <- contig_set_keys(membership)$key
  seeds <- membership
  parents <- data.table(bin_id = integer(), parent = integer())
  next_id <- max(bins$bin_id) + 1L

  for (it in seq_len(n_iter)) {
    log_msg("Bin merging, iteration ", it)

    seed_hits <- contig_markers[seeds, on = "contig", nomatch = NULL,
                                allow.cartesian = TRUE, .(bin_id, contig, marker)]
    seed_stats <- seed_hits[, .(n_unique = uniqueN(marker), n_total = .N), by = bin_id]
    candidates <- seed_stats[n_unique >= min_markers & n_total <= max_marker_hits]
    if (!nrow(candidates)) {
      log_msg("  no candidate seed bins; stopping")
      break
    }

    ## markers of each seed that lie on contigs of an original bin (partner)
    sh <- seed_hits[bin_id %in% candidates$bin_id]
    setnames(sh, "bin_id", "seed")
    pm <- membership[sh, on = "contig", nomatch = NULL, allow.cartesian = TRUE,
                     .(seed, partner = bin_id, marker)]
    pm <- unique(pm)[seed != partner]
    if (nrow(parents)) {
      pm <- pm[!parents, on = c(seed = "bin_id", partner = "parent")]
    }
    shared <- pm[, .(shared_markers = .N), by = .(seed, partner)]
    shared <- candidates[, .(seed = bin_id, n_unique)][shared, on = "seed"]
    shared[, shared_markers_rel := shared_markers / n_unique]
    shared <- shared[shared_markers_rel >= min_sharing]

    ## one entry per unordered pair (first by seed name, then partner name);
    ## then require an absolute number of shared markers
    shared <- all_bins[, .(seed = bin_id, seed_name = bin)][shared, on = "seed"]
    shared <- all_bins[, .(partner = bin_id, partner_name = bin)][shared, on = "partner"]
    setorder(shared, seed_name, partner_name)
    shared[, pair := paste(pmin(seed, partner), pmax(seed, partner))]
    shared <- shared[!duplicated(pair) & shared_markers >= min_markers]
    if (!nrow(shared)) {
      log_msg("  no (more) overlapping bins found; stopping")
      break
    }

    shared[, `:=`(new_id = next_id + seq_len(.N) - 1L,
                  new_name = paste0("MAGScoT_", it, "_", seq_len(.N)))]
    new_members <- unique(rbind(
      seeds[shared[, .(bin_id = seed, new_id)], on = "bin_id", allow.cartesian = TRUE, .(bin_id = new_id, contig)],
      membership[shared[, .(bin_id = partner, new_id)], on = "bin_id", allow.cartesian = TRUE, .(bin_id = new_id, contig)]
    ))

    keys <- contig_set_keys(new_members)
    keys <- keys[!key %in% known_keys & !duplicated(key)]
    n_new <- nrow(keys)
    log_msg(sprintf("  %d overlapping pairs, %d new hybrid bins", nrow(shared), n_new))
    if (!n_new) break

    keep <- keys$bin_id
    known_keys <- c(known_keys, keys$key)
    new_members <- new_members[bin_id %in% keep]
    kept <- shared[new_id %in% keep]
    all_members <- rbind(all_members, new_members)
    all_bins <- rbind(all_bins, data.table(bin_id = kept$new_id, bin = kept$new_name,
                                           bintool = paste0("MAGScoT_", it),
                                           is_merged = TRUE))
    parents <- rbind(parents,
                     kept[, .(bin_id = new_id, parent = seed)],
                     kept[, .(bin_id = new_id, parent = partner)])
    seeds <- new_members
    next_id <- max(all_bins$bin_id) + 1L
  }
  list(membership = all_members, bins = all_bins)
}
