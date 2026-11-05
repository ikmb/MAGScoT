MAGSCOT_VERSION <- "2.0.0"

## Command-line interface -------------------------------------------------------------

CLI_USAGE <- "Usage: MAGScoT.R [refine|annotate|extract] [options]

Commands:
  refine    Score, merge and refine bins (default if no command is given)
  annotate  Find GTDB r207 marker genes (prodigal + hmmsearch) without refining
  extract   Write FASTA files for refined bins

Run 'MAGScoT.R <command> --help' for the options of each command."

refine_options <- function() list(
  optparse::make_option(c("-i", "--input"), type = "character",
    help = "Tab-separated contig-to-bin file with three columns and no header: bin, contig, binning tool"),
  optparse::make_option("--hmm", type = "character",
    help = "Tab-separated marker file with three columns and no header: protein id, marker, e-value. Protein ids must be <contig>_<n> (prodigal convention); see 'MAGScoT.R annotate'"),
  optparse::make_option("--proteins", type = "character",
    help = "Protein FASTA from prodigal or pyrodigal (may be gzipped), instead of --hmm; marker genes are searched with hmmsearch and the hits written to <out>.hmm"),
  optparse::make_option("--threads", type = "integer", default = 1L,
    help = "Number of hmmsearch CPUs (with --proteins) [default: %default]"),
  optparse::make_option(c("-p", "--profile"), type = "character", default = "default",
    help = "Marker profile (GTDB r207): default (bac120 + ar53), bac120, ar53, or a path to a TSV with columns marker, set [default: %default]"),
  optparse::make_option(c("-o", "--out"), type = "character", default = "MAGScoT",
    help = "Output prefix [default: %default]"),
  optparse::make_option("--bin_separator", type = "character", default = "cleanbin",
    help = "Infix for refined bin names: <prefix>_<infix>_000001 [default: %default]"),
  optparse::make_option("--preset", type = "character", default = "dastool",
    help = "Scoring weight preset: dastool (a=1, b=0.5, c=0.5) or metawrap (a=1, b=5, c=0.5) [default: %default]"),
  optparse::make_option(c("-a", "--score_a"), type = "double",
    help = "Scoring weight a (completeness); overrides the preset"),
  optparse::make_option(c("-b", "--score_b"), type = "double",
    help = "Scoring weight b (contamination); overrides the preset"),
  optparse::make_option(c("-c", "--score_c"), type = "double",
    help = "Scoring weight c (surplus marker copies); overrides the preset"),
  optparse::make_option("--max_cont", type = "double", default = 1,
    help = "Maximum contamination; bins above this are not considered [default: %default]"),
  optparse::make_option(c("-t", "--threshold"), type = "double", default = 0.5,
    help = "Minimum score (and completeness) of refined bins [default: %default]"),
  optparse::make_option("--score_only", action = "store_true", default = FALSE,
    help = "Only score the input bins; no merging or refinement"),
  optparse::make_option("--skip_merge_bins", action = "store_true", default = FALSE,
    help = "Skip bin merging"),
  optparse::make_option(c("-m", "--min_markers"), type = "integer", default = 25L,
    help = "Minimum number of unique markers of a seed bin, and of shared markers, for bin merging [default: %default]"),
  optparse::make_option("--max_marker_hits", type = "integer", default = 150L,
    help = "Maximum total marker hits of a seed bin for bin merging [default: %default]"),
  optparse::make_option(c("-s", "--min_sharing"), type = "double", default = 0.8,
    help = "Minimum fraction of a seed bin's markers shared with a partner bin for merging [default: %default]"),
  optparse::make_option(c("-n", "--n_iterations"), type = "integer", default = 2L,
    help = "Number of bin-merging iterations [default: %default]"),
  optparse::make_option("--contigs", type = "character",
    help = "Optional assembly FASTA (may be gzipped); refined bins are written as FASTA files"),
  optparse::make_option("--bins_dir", type = "character",
    help = "Directory for refined bin FASTA files [default: <out>_bins]"),
  optparse::make_option(c("-q", "--quiet"), action = "store_true", default = FALSE,
    help = "Suppress progress messages")
)

annotate_options <- function() list(
  optparse::make_option("--contigs", type = "character",
    help = "Assembly FASTA (may be gzipped); genes are predicted with prodigal"),
  optparse::make_option("--proteins", type = "character",
    help = "Protein FASTA from prodigal or pyrodigal (may be gzipped), instead of --contigs"),
  optparse::make_option(c("-o", "--out"), type = "character", default = "MAGScoT",
    help = "Output prefix; writes <out>.hmm (and <out>.prodigal.faa with --contigs) [default: %default]"),
  optparse::make_option("--threads", type = "integer", default = 1L,
    help = "Number of threads [default: %default]"),
  optparse::make_option("--keep_tmp", action = "store_true", default = FALSE,
    help = "Keep intermediate files"),
  optparse::make_option(c("-q", "--quiet"), action = "store_true", default = FALSE,
    help = "Suppress progress messages")
)

extract_options <- function() list(
  optparse::make_option(c("-i", "--input"), type = "character",
    help = "Refined contig-to-bin file (<out>.refined.contig_to_bin.out)"),
  optparse::make_option("--contigs", type = "character",
    help = "Assembly FASTA (may be gzipped)"),
  optparse::make_option(c("-o", "--out"), type = "character", default = "MAGScoT_bins",
    help = "Output directory [default: %default]"),
  optparse::make_option(c("-q", "--quiet"), action = "store_true", default = FALSE,
    help = "Suppress progress messages")
)

require_args <- function(opt, parser, names) {
  missing <- names[vapply(names, function(n) is.null(opt[[n]]), logical(1))]
  if (length(missing)) {
    optparse::print_help(parser)
    stop("Missing required option(s): ", paste0("--", missing, collapse = ", "), call. = FALSE)
  }
}

require_one_of <- function(opt, parser, names) {
  given <- names[!vapply(names, function(n) is.null(opt[[n]]), logical(1))]
  if (length(given) != 1L) {
    optparse::print_help(parser)
    stop("Give exactly one of: ", paste0("--", names, collapse = ", "), call. = FALSE)
  }
}

## MAGScoT command-line interface
##
## Entry point used by `MAGScoT.R`. Run with `"--help"` for usage.
##
## @param args Command-line arguments, e.g. `c("refine", "-i", "bins.tsv", "--hmm", "markers.hmm")`.
## @return Invisibly, the exit status (0 on success).

cli_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  commands <- c("refine", "annotate", "extract")
  cmd <- "refine"
  if (length(args) && args[1] %in% commands) {
    cmd <- args[1]
    args <- args[-1]
  } else if (length(args) && args[1] %in% c("help", "--version")) {
    if (args[1] == "--version") cat("MAGScoT", magscot_version(), "\n") else cat(CLI_USAGE, "\n")
    return(invisible(0L))
  }

  opts <- switch(cmd, refine = refine_options(), annotate = annotate_options(),
                 extract = extract_options())
  prog <- if (cmd == "refine") "MAGScoT.R [refine]" else paste("MAGScoT.R", cmd)
  parser <- optparse::OptionParser(usage = paste("%prog [options]"), option_list = opts,
                                   prog = prog,
                                   epilogue = if (cmd == "refine") CLI_USAGE else "")
  opt <- optparse::parse_args(parser, args = args)
  old <- options(MAGScoT.verbose = !isTRUE(opt$quiet))
  on.exit(options(old))

  status <- tryCatch({
    if (cmd == "refine") {
      require_args(opt, parser, "input")
      require_one_of(opt, parser, c("hmm", "proteins"))
      magscot_refine(
        input = opt$input, hmm = opt$hmm, proteins = opt$proteins, threads = opt$threads,
        out = opt$out, profile = opt$profile,
        bin_separator = opt$bin_separator,
        preset = opt$preset,
        score_a = opt$score_a, score_b = opt$score_b, score_c = opt$score_c,
        threshold = opt$threshold, max_cont = opt$max_cont,
        score_only = opt$score_only, skip_merge_bins = opt$skip_merge_bins,
        min_markers = opt$min_markers, min_sharing = opt$min_sharing,
        n_iterations = opt$n_iterations, max_marker_hits = opt$max_marker_hits,
        contigs = opt$contigs,
        bins_dir = if (is.null(opt$bins_dir)) paste0(opt$out, "_bins") else opt$bins_dir)
    } else if (cmd == "annotate") {
      require_one_of(opt, parser, c("contigs", "proteins"))
      magscot_annotate(contigs = opt$contigs, proteins = opt$proteins, out = opt$out,
                       threads = opt$threads, keep_tmp = opt$keep_tmp)
    } else {
      require_args(opt, parser, c("input", "contigs"))
      write_bin_fastas(opt$input, opt$contigs, opt$out)
    }
    0L
  }, error = function(e) {
    message("Error: ", conditionMessage(e))
    1L
  })
  if (status != 0L && !interactive()) quit(save = "no", status = status)
  invisible(status)
}

magscot_version <- function() MAGSCOT_VERSION
