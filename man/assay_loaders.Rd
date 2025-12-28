\name{load_cytof_assay}
\alias{load_cytof_assay}
\alias{load_scrnaseq_counts}
\alias{validate_assay_inputs}
\title{Load and validate single-cell assay inputs}
\usage{
load_cytof_assay(
  exprs,
  cell_metadata = NULL,
  marker_columns = NULL,
  arcsinh_cofactor = 5,
  class_column = "label",
  cell_id_column = "cell_id",
  qc_max_na = 0.05
)

load_scrnaseq_counts(
  counts,
  cell_metadata = NULL,
  feature_metadata = NULL,
  normalization = c("log1p", "cpm", "none"),
  cell_id_column = "cell_id",
  gene_id_column = "gene_id",
  min_cells = 5,
  min_features = 200,
  class_column = "label"
)

validate_assay_inputs(
  feature_matrix,
  cell_metadata,
  cell_id_column = "cell_id",
  class_column = NULL,
  required_metadata = NULL
)
}
\arguments{
\item{exprs}{A data.frame or matrix of CyTOF expression values, or a path to a
delimited file containing them.}

\item{cell_metadata}{Optional data.frame of per-cell annotations. For CyTOF
inputs, non-marker columns in \code{exprs} are treated as metadata when this
argument is \code{NULL}.}

\item{marker_columns}{Optional character vector listing the columns that
represent markers. Defaults to all numeric columns after metadata is removed.}

\item{arcsinh_cofactor}{Cofactor used to arcsinh-transform CyTOF intensities.
Set to \code{NULL} to skip transformation.}

\item{class_column}{Metadata column containing class labels. Created and filled
with \code{NA} values when missing.}

\item{cell_id_column}{Metadata column containing unique cell identifiers.
Created automatically when absent.}

\item{qc_max_na}{Maximum allowed fraction of missing values per cell. Cells
exceeding this value are removed.}

\item{counts}{Matrix or data.frame of raw scRNA-seq counts (genes as rows,
cells as columns), or a path to a delimited file containing them.}

\item{feature_metadata}{Optional data.frame describing genes or features.}

\item{normalization}{Normalization to apply to counts. Supported options are
\code{"log1p"}, \code{"cpm"}, and \code{"none"}.}

\item{gene_id_column}{Column name storing gene identifiers when reading from
tabular counts files.}

\item{min_cells}{Minimum number of cells a gene must be expressed in to be
retained.}

\item{min_features}{Minimum number of detected genes required for a cell to be
retained.}

\item{feature_matrix}{Matrix or data.frame where rows correspond to cells and
columns correspond to features.}

\item{required_metadata}{Optional character vector of metadata columns that
must exist. Missing entries are created with \code{NA}.}
}
\value{
For \code{load_cytof_assay}, a list containing transformed expression values
(\code{exprs}), harmonized \code{cell_metadata}, generated \code{feature_metadata},
and an \code{assay_type} flag. For \code{load_scrnaseq_counts}, a list containing
raw \code{counts} (genes x cells), normalized counts (\code{normalized_counts}, cells x genes),
\code{cell_metadata}, \code{feature_metadata}, and an \code{assay_type} flag.
For \code{validate_assay_inputs}, a list with harmonized \code{features} and
\code{cell_metadata}.
}
\description{
Helpers for importing updated CyTOF and scRNA-seq assay formats with light
quality control, normalization, and metadata validation. A shared validator
ensures consistent identifiers and required columns across assay types.
}
\examples{
\dontrun{
cytof <- load_cytof_assay(
  exprs = system.file("data", "CyTOF_updatedExample_exprs.csv", package = "RTLbase"),
  cell_metadata = read.csv(system.file("data", "CyTOF_updatedExample_metadata.csv", package = "RTLbase"))
)

rna <- load_scrnaseq_counts(
  counts = system.file("data", "scRNAseq_updatedExample_counts.csv", package = "RTLbase"),
  cell_metadata = read.csv(system.file("data", "scRNAseq_updatedExample_cell_metadata.csv", package = "RTLbase")),
  feature_metadata = read.csv(system.file("data", "scRNAseq_updatedExample_feature_metadata.csv", package = "RTLbase")),
  normalization = "log1p"
)

validated <- validate_assay_inputs(cytof$exprs, cytof$cell_metadata)
}
}
