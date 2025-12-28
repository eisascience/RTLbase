\name{scRNAseq_updatedExample}
\alias{scRNAseq_updatedExample}
\title{Updated scRNA-seq example dataset}
\format{Three CSV files:
\itemize{
\item \verb{scRNAseq_updatedExample_counts.csv}: gene-by-cell raw counts with
  a \code{gene_id} column.
\item \verb{scRNAseq_updatedExample_cell_metadata.csv}: cell-level annotations
  including \code{cell_id}, \code{sample_id}, \code{cluster}, and \code{label}.
\item \verb{scRNAseq_updatedExample_feature_metadata.csv}: gene-level metadata.
}}
\description{
A compact scRNA-seq count matrix and accompanying metadata, shaped to mirror
updated public data releases and ready for downstream validation with
\code{load_scrnaseq_counts()}.
}
\usage{
scrna <- load_scrnaseq_counts(
  counts = system.file("data", "scRNAseq_updatedExample_counts.csv", package = "RTLbase"),
  cell_metadata = read.csv(system.file("data", "scRNAseq_updatedExample_cell_metadata.csv", package = "RTLbase")),
  feature_metadata = read.csv(system.file("data", "scRNAseq_updatedExample_feature_metadata.csv", package = "RTLbase"))
)
}
\details{
Counts were generated from a balanced mix of three sample IDs and two response
labels. Low-variance and low-detection features were retained to keep the
example small while still exercising preprocessing steps in the vignettes.
}
