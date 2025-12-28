\name{alg3_shiftComp}
\alias{alg3_shiftComp}
\title{Apply shift compensation using maximum cross-correlation}
\usage{
alg3_shiftComp(
  source_list,
  task_list,
  alg1_result,
  alg2_result,
  print2screen,
  save2file,
  maximumLag,
  ImpFeats,
  ADM = FALSE,
  datatyp = "FC",
  useAbsCor = TRUE,
  medianMediansBL = FALSE,
  CoreClassifier
)
}
\arguments{
\item{source_list}{List of source feature matrices used to train the classifiers.}

\item{task_list}{List of target task feature matrices to be aligned.}

\item{alg1_result}{Output list from \code{\link{alg1_baselineClass}}.}

\item{alg2_result}{Output list from \code{\link{alg2_rob_meanNCov}}.}

\item{print2screen}{Logical; print progress updates.}

\item{save2file}{Logical; save diagnostic plots to disk.}

\item{maximumLag}{Maximum lag to consider in cross-correlation.}

\item{ImpFeats}{Optional vector of feature indices to retain.}

\item{ADM}{Logical; apply additional data manipulations.}

\item{datatyp}{Character flag indicating data type (e.g., \code{"FC"}).}

\item{useAbsCor}{Logical; use absolute correlations when selecting lags.}

\item{medianMediansBL}{Logical; use the median of medians strategy for shift estimation.}

\item{CoreClassifier}{Core classifier identifier (e.g., \code{"LinSVM"}).}
}
\value{
A numeric vector of updated intercepts per task.
}
\description{
Aligns target task projections to source datasets by estimating robust lag
shifts via maximum cross-correlation between kernel density estimates. The
updated intercepts are returned for downstream bias correction.
}
