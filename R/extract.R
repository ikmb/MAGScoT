## Write refined bins as FASTA files -----------------------------------------------

## Write refined bins as FASTA files
##
## Writes one FASTA file per refined bin. The assembly is streamed in chunks, so
## memory use stays low for large assemblies. Sequence headers are matched on
## their first word.
##
## @param assignments A data frame with columns `binnew` and `contig`, or the
##   path to a `<out>.refined.contig_to_bin.out` file.
## @param fasta Assembly FASTA (may be gzipped).
## @param outdir Output directory (created if needed).
## @param extension File extension of the bin files.
## @param chunk_lines Number of FASTA lines to read at a time.
## @return Invisibly, the paths of the written files.

write_bin_fastas <- function(assignments, fasta, outdir, extension = ".fasta",
                             chunk_lines = 1e6) {
  if (is.character(assignments)) {
    assignments <- fread(assignments, sep = "\t", header = TRUE, colClasses = "character")
  }
  assignments <- as.data.table(assignments)
  if (!all(c("binnew", "contig") %in% names(assignments))) {
    stop("Bin assignments need the columns 'binnew' and 'contig'.", call. = FALSE)
  }
  if (anyDuplicated(assignments$contig)) {
    stop("A contig is assigned to more than one refined bin.", call. = FALSE)
  }
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  files <- file.path(outdir, paste0(unique(assignments$binnew), extension))
  existing <- files[file.exists(files)]
  if (length(existing)) unlink(existing)

  lookup <- setNames(assignments$binnew, assignments$contig)
  con <- open_text(fasta)
  on.exit(close(con))
  current <- NA_character_
  found <- 0L
  repeat {
    lines <- readLines(con, n = chunk_lines, warn = FALSE)
    if (!length(lines)) break
    is_header <- startsWith(lines, ">")
    ids <- sub("^>([^[:space:]]*).*$", "\\1", lines[is_header])
    hdr_bin <- unname(lookup[ids])
    found <- found + sum(!is.na(hdr_bin))
    ## carry the bin of the last header in the previous chunk into this one
    line_bin <- c(current, hdr_bin)[cumsum(is_header) + 1L]
    current <- line_bin[length(line_bin)]
    keep <- !is.na(line_bin)
    if (any(keep)) {
      by_bin <- split(lines[keep], line_bin[keep])
      for (bn in names(by_bin)) {
        cat(by_bin[[bn]], file = file.path(outdir, paste0(bn, extension)),
            sep = "\n", append = TRUE)
      }
    }
  }
  missing <- length(lookup) - found
  if (missing > 0L) {
    warning(missing, " assigned contig(s) were not found in '", fasta, "'.", call. = FALSE)
  }
  log_msg(sprintf("Wrote %d bin FASTA files to: %s", uniqueN(assignments$binnew), outdir))
  invisible(files)
}
