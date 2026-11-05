# MAGScoT 2.0.0

Rewrite of the R implementation. The code is split into modules in `R/`,
loaded by `MAGScoT.R`. `Rscript MAGScoT.R -i ... --hmm ...` works as before
and produces the same output files and columns.

## Bug fixes

* Bin merging iterations after the first had no effect: hybrid bins were built
  from the partner bin only, without the contigs of the (hybrid) seed bin.
* A merging iteration without overlapping bins crashed with
  "Join columns in `x` must be present in the data" (#9).
* Bin names that occur in more than one binning tool were pooled into one bin
  (#4). Bins are now identified by binning tool and name.
* Refinement could select bins below the score cutoff when several bins of the
  same binning tool were accepted in a row.
* A last remaining bin consisting of a single contig was never evaluated.
* Hashes used to detect duplicate hybrid bins could collide (contig indices were
  concatenated without separator); contig sets are now compared directly.
* Refined bin names used the full `--out` path instead of its basename.

## New features

* `--proteins`: give a prodigal / pyrodigal protein FASTA instead of a marker
  file; MAGScoT runs hmmsearch itself and saves the hits as `<out>.hmm` (#1).
* `MAGScoT.R annotate`: gene prediction (prodigal, parallelised) and marker
  annotation (hmmsearch) without refinement; also accepts `--proteins`.
* `--contigs` / `MAGScoT.R extract`: write refined bins as FASTA files (#10).
* `--preset dastool|metawrap`: scoring weight presets. The default stays
  `dastool` (a = 1, b = 0.5, c = 0.5); `metawrap` weights contamination with
  b = 5. Explicit `-a`/`-b`/`-c` override the preset.
* `--max_cont` (added in v1.1) is renamed to `--max_contamination`, and its
  help explains that it is a fraction of the `Contamination` column. The old
  name still works and prints a deprecation note.
* `--profile` accepts a custom marker TSV; new option `--max_marker_hits`
  (previously fixed at 150) and `--quiet`.
* Input checks with clear error messages: wrong number of columns, swapped
  columns, contigs assigned to several bins of one tool, markers that do not
  match the profile, protein IDs that do not match contig IDs. Whitespace-
  separated and CRLF files are accepted.
* With only one binning tool, bins are now refined (filtered and renamed)
  instead of only scored.

## Other changes

* Dependencies reduced to `data.table` and `optparse` (no dplyr, readr, funr,
  digest). Bin merging is vectorised.
* Progress messages go to stderr.
* Test suite (testthat), including a regression test against v1.1 output.
