## Load MAGScoT from the R/ directory of this repository.
library(data.table)
magscot_root <- normalizePath(testthat::test_path("..", ".."))
for (f in sort(list.files(file.path(magscot_root, "R"), pattern = "[.]R$", full.names = TRUE))) {
  source(f, local = TRUE)
}
options(MAGScoT.data_dir = magscot_root)
