\name{alg1_baselineClass}
\alias{alg1_baselineClass}
\title{Train baseline SVM classifiers for multiple source datasets}
\usage{
alg1_baselineClass(
  TrainXls,
  TrainYls,
  TestXls,
  TestYls,
  K_forCrossV,
  svmGamma,
  svmCost,
  prnt2scr,
  X_cols2Keep,
  transX = FALSE,
  sampleRed = FALSE,
  doParalellSVM = FALSE,
  datatyp = "FC"
)
}
\arguments{
\item{TrainXls}{List of training feature matrices (one per source dataset).}

\item{TrainYls}{List of training labels aligned to \code{TrainXls}.}

\item{TestXls}{List of test feature matrices (one per target task).}

\item{TestYls}{List of test labels aligned to \code{TestXls}.}

\item{K_forCrossV}{Number of folds for internal SVM cross-validation.}

\item{svmGamma}{Radial basis gamma parameter passed to \code{e1071::svm}.}

\item{svmCost}{Cost parameter passed to \code{e1071::svm}.}

\item{prnt2scr}{Logical; print progress to the console.}

\item{X_cols2Keep}{Column indices to retain from the feature matrices.}

\item{transX}{Logical; apply additional transformations via \code{RTL::AllDataManipulations}.}

\item{sampleRed}{Optional integer sample size for down-sampling during training.}

\item{doParalellSVM}{Logical; use a parallelized SVM fit when \code{TRUE}.}

\item{datatyp}{Character flag indicating data type (e.g., \code{"FC"} for flow cytometry).}
}
\value{
A list containing the baseline hyperplanes, per-dataset results, and metadata
(number of datasets and dimensions).
}
\description{
Fits a linear SVM for each source dataset and returns the learned hyperplane
parameters along with cross-validation diagnostics. This is the entry point for
initializing the RTL workflow before transfer learning adjustments are applied.
}
