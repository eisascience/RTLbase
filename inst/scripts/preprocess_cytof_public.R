# Example preprocessing workflow for public CyTOF exports
# ------------------------------------------------------
# This script downloads a tabular CyTOF export, applies arcsinh normalization,
# harmonizes metadata, validates the result with RTLbase helpers, and writes
# CSV outputs that can be bundled with the package or used directly.

# Required packages: data.table (for fast I/O), utils (base), RTLbase

library(utils)
library(data.table)
library(RTLbase)

set.seed(7)

# Configure locations -------------------------------------------------------
# Replace these URLs/paths with the dataset you want to process.
cytof_url <- "https://example.org/public-cytof-export.csv"
metadata_url <- "https://example.org/public-cytof-metadata.csv"
working_dir <- tempdir()
exprs_path <- file.path(working_dir, "cytof_exprs.csv")
metadata_path <- file.path(working_dir, "cytof_metadata.csv")
feature_metadata_path <- file.path(working_dir, "cytof_feature_metadata.csv")

# Download raw files -------------------------------------------------------
# Comment these lines if you already have local files.
try(download.file(cytof_url, exprs_path, quiet = TRUE), silent = TRUE)
try(download.file(metadata_url, metadata_path, quiet = TRUE), silent = TRUE)

# Load the expression matrix (cells x markers)
exprs_dt <- fread(exprs_path)

# Derive or load feature annotations
marker_columns <- setdiff(colnames(exprs_dt), c("File", "EventLength", "Time"))
feature_metadata <- data.frame(
  marker = marker_columns,
  channel = marker_columns,
  measurement_type = "protein",
  stringsAsFactors = FALSE
)
write.csv(feature_metadata, feature_metadata_path, row.names = FALSE)

# Load and tidy cell-level metadata
cell_metadata <- fread(metadata_path)
if (!"cell_id" %in% colnames(cell_metadata)) {
  cell_metadata[, cell_id := paste0("cell_", .I)]
}

# Minimal QC: drop empty events and arcsinh-transform markers
exprs_dt <- exprs_dt[rowSums(is.na(.SD)) < ncol(.SD) * 0.1, ..marker_columns]
exprs_dt <- as.data.frame(exprs_dt)
exprs_dt <- asinh(exprs_dt / 5)

# Validate and write outputs ----------------------------------------------
validated <- load_cytof_assay(
  exprs = exprs_dt,
  cell_metadata = as.data.frame(cell_metadata),
  marker_columns = marker_columns,
  class_column = "label",
  cell_id_column = "cell_id"
)

write.csv(validated$exprs, exprs_path, row.names = FALSE)
write.csv(validated$cell_metadata, metadata_path, row.names = FALSE)
write.csv(validated$feature_metadata, feature_metadata_path, row.names = FALSE)

message("CyTOF preprocessing complete: ", exprs_path)
