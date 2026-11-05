## Input / output helpers -----------------------------------------------------

## Location of the bundled marker profiles (profiles/) and HMMs (hmm/): the
## MAGScoT directory, which MAGScoT.R stores in the option `MAGScoT.data_dir`.
magscot_data_dir <- function() {
  d <- getOption("MAGScoT.data_dir", "")
  if (!nzchar(d) || !dir.exists(file.path(d, "profiles"))) {
    stop("Cannot locate MAGScoT data directory (marker profiles / HMMs).", call. = FALSE)
  }
  d
}

log_msg <- function(...) {
  if (isTRUE(getOption("MAGScoT.verbose", TRUE))) message(...)
}

open_text <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  if (grepl("[.]gz$", path)) gzfile(path, "rt") else file(path, "rt")
}

## Read a headerless table that should have `n` columns. Tab-separated input
## is expected, but whitespace-separated files (e.g. column-aligned output) are
## accepted as a fallback because that is a common source of user errors.
read_columns <- function(path, col_names, what) {
  con <- open_text(path)
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  n <- length(col_names)
  if (!length(lines)) stop("The ", what, " file '", path, "' is empty.", call. = FALSE)

  fields <- strsplit(lines, "\t", fixed = TRUE)
  nf <- lengths(fields)
  if (all(nf == 1L)) {
    fields <- strsplit(trimws(lines), "[[:space:]]+")
    nf <- lengths(fields)
    if (all(nf == n)) {
      warning("The ", what, " file '", path, "' is not tab-separated; ",
              "splitting columns on whitespace instead.", call. = FALSE)
    }
  }
  if (any(nf != n)) {
    bad <- which(nf != n)[1]
    stop(sprintf(paste0(
      "The %s file '%s' must have %d tab-separated columns without a header (%s).\n",
      "Line %d has %d column(s): %s"),
      what, path, n, paste(col_names, collapse = ", "), bad, nf[bad],
      substr(lines[bad], 1, 200)), call. = FALSE)
  }
  dt <- as.data.table(do.call(rbind, fields))
  setnames(dt, col_names)
  dt[]
}

## Read the contig-to-bin table (bin, contig, binning tool; no header).
read_contig_to_bin <- function(path) {
  ctb <- read_columns(path, c("bin", "contig", "set"), "contig-to-bin")
  ctb <- ctb[!is.na(set) & nzchar(set)]
  ctb <- unique(ctb)

  ## within one binning tool there cannot be more bins than contigs, unless
  ## the bin and contig columns are swapped
  per_set <- ctb[, .(n_bins = uniqueN(bin), n_contigs = uniqueN(contig)), by = set]
  if (any(per_set$n_bins > per_set$n_contigs)) {
    stop("The contig-to-bin file has more distinct bins than contigs for binning tool '",
         per_set[n_bins > n_contigs, set][1],
         "'. Expected column order: bin, contig, binning tool.", call. = FALSE)
  }

  dup <- ctb[, .N, by = .(set, contig)][N > 1L]
  if (nrow(dup)) {
    stop(sprintf(paste0(
      "%d contig(s) are assigned to more than one bin by the same binning tool ",
      "(e.g. contig '%s' in set '%s'). Each binning tool must assign a contig to at ",
      "most one bin. If these are different runs of one tool, give each run its ",
      "own name in the third column."),
      nrow(dup), dup$contig[1], dup$set[1]), call. = FALSE)
  }

  shared <- unique(ctb[, .(bin, set)])[, .N, by = bin][N > 1L]
  if (nrow(shared)) {
    log_msg(sprintf(
      "Note: %d bin name(s) occur in more than one binning tool; bins are identified by tool + name.",
      nrow(shared)))
  }
  ctb[]
}

## Read marker annotations (protein id, marker, e-value; no header).
read_hmm <- function(path) {
  hmm <- read_columns(path, c("gene", "marker", "evalue"), "HMM")
  hmm[, evalue := suppressWarnings(as.numeric(evalue))]
  if (anyNA(hmm$evalue)) {
    stop("The third column of the HMM file '", path, "' must contain numeric e-values.",
         call. = FALSE)
  }
  hmm[]
}

## Map marker hits (gene, marker, evalue) to contigs. Only the best hit per
## protein is kept, and a marker found on several genes of the same contig is
## counted once. Protein ids are <contig>_<n> (prodigal / pyrodigal convention).
hits_to_contigs <- function(hits, gene_suffix = "_[0-9]+$") {
  hits <- copy(hits)
  setorder(hits, evalue)
  hits <- hits[!duplicated(gene)]
  hits[, contig := sub(gene_suffix, "", gene)]
  unique(hits[, .(contig, marker)])
}

## Read a marker profile: "default" (bac120 + ar53), "bac120", "ar53", or the
## path to a custom two-column TSV with header `marker<TAB>set`.
read_markers <- function(profile = "default") {
  builtin <- c(default = "gtdb_rel207_default_markers.tsv",
               bac120  = "gtdb_rel207_bac120_markers.tsv",
               ar53    = "gtdb_rel207_ar53_markers.tsv")
  path <- if (profile %in% names(builtin)) {
    file.path(magscot_data_dir(), "profiles", builtin[[profile]])
  } else if (file.exists(profile)) {
    profile
  } else {
    stop("Unknown profile '", profile, "'. Use one of: ",
         paste(names(builtin), collapse = ", "), ", or a path to a marker TSV.",
         call. = FALSE)
  }
  markers <- fread(path, header = TRUE, sep = "\t", colClasses = "character")
  if (!all(c("marker", "set") %in% names(markers))) {
    stop("Marker profile '", path, "' needs the columns 'marker' and 'set'.", call. = FALSE)
  }
  markers <- unique(markers[, .(marker = sub("[.](hmm|HMM)$", "", marker), set)])
  markers[]
}

write_tsv <- function(x, path) {
  fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
}
