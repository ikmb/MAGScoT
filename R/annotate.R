## Marker annotation: prodigal + hmmsearch ------------------------------------------

require_tool <- function(tool) {
  path <- Sys.which(tool)
  if (!nzchar(path)) {
    stop("'", tool, "' was not found on the PATH. Install it (e.g. `conda install -c bioconda ",
         tool, "`) and try again.", call. = FALSE)
  }
  path
}

run_tool <- function(cmd, args) {
  status <- system2(cmd, args, stdout = FALSE, stderr = FALSE)
  if (!identical(status, 0L)) {
    stop(basename(cmd), " failed (exit status ", status, "): ",
         paste(c(cmd, args), collapse = " "), call. = FALSE)
  }
}

## Split a (gzipped) FASTA into `n` plain-text chunks of roughly equal size,
## assigning whole records round-robin.
split_fasta <- function(fasta, n, dir, chunk_lines = 1e6) {
  files <- file.path(dir, sprintf("chunk_%03d.fna", seq_len(n)))
  cons <- lapply(files, file, open = "wt")
  on.exit(lapply(cons, close))
  con <- open_text(fasta)
  on.exit(close(con), add = TRUE)
  rec <- 0L
  repeat {
    lines <- readLines(con, n = chunk_lines, warn = FALSE)
    if (!length(lines)) break
    is_header <- startsWith(lines, ">")
    rec_no <- rec + cumsum(is_header)
    rec <- rec_no[length(rec_no)]
    target <- ((rec_no - 1L) %% n) + 1L
    target[rec_no == 0L] <- 1L
    for (k in unique(target)) writeLines(lines[target == k], cons[[k]])
  }
  files
}

## Parse hmmsearch --tblout output into (gene, marker, evalue). The marker is
## the HMM accession (e.g. PF00380.20, TIGR00006), or its name if it has none.
parse_tblout <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!startsWith(lines, "#") & nzchar(lines)]
  if (!length(lines)) {
    return(data.table(gene = character(), marker = character(), evalue = numeric()))
  }
  f <- strsplit(lines, "[[:space:]]+")
  acc <- vapply(f, `[`, "", 4L)
  name <- vapply(f, `[`, "", 3L)
  data.table(gene = vapply(f, `[`, "", 1L),
             marker = ifelse(acc == "-", name, acc),
             evalue = as.numeric(vapply(f, `[`, "", 5L)))
}

## Search the GTDB r207 marker HMMs (TIGRFAM + Pfam, trusted cutoffs) in a
## protein FASTA (prodigal / pyrodigal output, may be gzipped) with hmmsearch.
## Returns data.table(gene, marker, evalue).
search_markers <- function(proteins, threads = 1L, tmp = tempfile("magscot_hmm_")) {
  hmmsearch <- require_tool("hmmsearch")
  if (!file.exists(proteins)) stop("File not found: ", proteins, call. = FALSE)
  threads <- max(1L, as.integer(threads))
  hmms <- file.path(magscot_data_dir(), "hmm",
                    c("gtdbtk_rel207_tigrfam.hmm", "gtdbtk_rel207_Pfam-A.hmm"))
  if (!dir.exists(tmp)) dir.create(tmp, recursive = TRUE)

  faa <- proteins
  if (grepl("[.]gz$", proteins)) {
    faa <- file.path(tmp, "proteins.faa")
    split_fasta(proteins, 1L, tmp)
    file.rename(file.path(tmp, "chunk_001.fna"), faa)
  }
  first <- readLines(faa, n = 1L, warn = FALSE)
  if (!length(first) || !startsWith(first, ">")) {
    stop("'", proteins, "' does not look like a FASTA file.", call. = FALSE)
  }

  log_msg("Searching GTDB r207 marker HMMs with hmmsearch (", threads, " CPU(s))...")
  hits <- rbindlist(lapply(seq_along(hmms), function(k) {
    tbl <- file.path(tmp, paste0("hits_", k, ".tbl"))
    run_tool(hmmsearch, c("--tblout", shQuote(tbl), "-o", "/dev/null", "--noali", "--notextw",
                          "--cut_nc", "--cpu", threads, shQuote(hmms[k]), shQuote(faa)))
    parse_tblout(tbl)
  }))
  if (!nrow(hits)) warning("hmmsearch found no marker genes in '", proteins, "'.", call. = FALSE)
  hits
}

## Annotate GTDB r207 marker genes
##
## Predicts genes with prodigal (unless a protein FASTA is given) and searches
## the GTDB r207 marker HMMs with hmmsearch. Writes `<out>.prodigal.faa` (when
## predicting genes) and `<out>.hmm`, the marker file used by
## magscot_refine(). `prodigal` and `hmmsearch` must be on the PATH.
##
## @param contigs Assembly FASTA (may be gzipped).
## @param proteins Protein FASTA from prodigal or pyrodigal (may be gzipped);
##   use instead of `contigs` to skip gene prediction.
## @param out Output prefix.
## @param threads Number of parallel prodigal jobs and hmmsearch CPUs.
## @param keep_tmp Keep intermediate files.
## @return Invisibly, the path of the marker file.
magscot_annotate <- function(contigs = NULL, proteins = NULL, out = "MAGScoT",
                             threads = 1L, keep_tmp = FALSE) {
  if (is.null(contigs) == is.null(proteins)) {
    stop("Give either an assembly (--contigs) or a protein FASTA (--proteins).", call. = FALSE)
  }
  threads <- max(1L, as.integer(threads))
  out_dir <- dirname(out)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  tmp <- tempfile("magscot_", tmpdir = out_dir)
  dir.create(tmp)
  if (!keep_tmp) on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  if (is.null(proteins)) proteins <- predict_genes(contigs, out, threads, tmp)
  hits <- search_markers(proteins, threads, tmp)
  hmm_file <- paste0(out, ".hmm")
  fwrite(hits, hmm_file, sep = "\t", col.names = FALSE, quote = FALSE)
  log_msg(sprintf("Wrote %d marker hits to: %s", nrow(hits), hmm_file))
  invisible(hmm_file)
}

## Predict genes with prodigal in `threads` parallel chunks; returns the path
## of <out>.prodigal.faa.
predict_genes <- function(contigs, out, threads, tmp) {
  prodigal <- require_tool("prodigal")
  log_msg("Predicting genes with prodigal (", threads, " job(s))...")
  chunks <- split_fasta(contigs, threads, tmp)
  faas <- sub("[.]fna$", ".faa", chunks)
  run_one <- function(k) {
    run_tool(prodigal, c("-p", "meta", "-q", "-i", shQuote(chunks[k]),
                         "-a", shQuote(faas[k]), "-o", shQuote(sub("[.]fna$", ".gbk", chunks[k]))))
    TRUE
  }
  if (threads > 1L && .Platform$OS.type == "unix") {
    ok <- parallel::mclapply(seq_along(chunks), run_one, mc.cores = threads)
    errs <- vapply(ok, inherits, logical(1), "try-error")
    if (any(errs)) stop(attr(ok[[which(errs)[1]]], "condition"))
  } else {
    lapply(seq_along(chunks), run_one)
  }
  faa <- paste0(out, ".prodigal.faa")
  file.copy(faas[1], faa, overwrite = TRUE)
  if (length(faas) > 1L) file.append(faa, faas[-1])
  log_msg("Proteins written to: ", faa)
  faa
}
