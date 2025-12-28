\name{CyTOF_updatedExample}
\alias{CyTOF_updatedExample}
\title{Updated CyTOF example dataset}
\format{Three CSV files:
\itemize{
\item \verb{CyTOF_updatedExample_exprs.csv}: marker intensities with one row per
  cell and one column per marker.
\item \verb{CyTOF_updatedExample_metadata.csv}: cell-level metadata including
  \code{cell_id}, \code{sample_id}, \code{condition}, and \code{label}.
\item \verb{CyTOF_updatedExample_feature_metadata.csv}: marker/channel
  annotations.
}}
\description{
A lightweight CyTOF table modeled after updated panel exports. The CSV trio
pairs cell-level metadata with marker intensities suitable for the RTL transfer
learning helpers and vignettes.
}
\usage{
cytof <- load_cytof_assay(
  exprs = system.file("data", "CyTOF_updatedExample_exprs.csv", package = "RTLbase"),
  cell_metadata = read.csv(system.file("data", "CyTOF_updatedExample_metadata.csv", package = "RTLbase")),
  class_column = "label"
)
}
\details{
Values are untransformed exports; use \code{load_cytof_assay()} to apply the
default arcsinh transform (cofactor = 5) and harmonize metadata. Marker names
reflect common lineage and activation channels used in immune profiling.
}
