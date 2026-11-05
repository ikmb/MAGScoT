#!/usr/bin/env Rscript
### MAGScoT - MAG Scoring and Refinement Tool
### https://github.com/ikmb/MAGScoT
###
### Command-line entry point; the implementation lives in R/ next to this script.

local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_dir <- if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  } else {
    getwd()
  }

  missing <- c("data.table", "optparse")[
    !vapply(c("data.table", "optparse"), requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Error: MAGScoT needs the R package(s): ", paste(missing, collapse = ", "),
            ". Install them, e.g. with: conda install -c conda-forge ",
            paste0("r-", tolower(missing), collapse = " "))
    quit(save = "no", status = 1)
  }
  suppressPackageStartupMessages(library(data.table))

  env <- new.env()
  for (f in sort(list.files(file.path(script_dir, "R"), pattern = "[.]R$", full.names = TRUE))) {
    sys.source(f, envir = env)
  }
  options(MAGScoT.data_dir = script_dir)
  env$cli_main(commandArgs(trailingOnly = TRUE))
})
