\name{alg6_NormalVectorUpdate}
\alias{alg6_NormalVectorUpdate}
\title{Update normal vectors using Algorithm 6}
\usage{
alg6_NormalVectorUpdate(
  task_list,
  alg1_result,
  alg2_result,
  alg3_result,
  alg4_result,
  X_feat_cols,
  Marg,
  save2file,
  ADM,
  datatyp = "FC",
  RCSmodeBL = FALSE,
  RCSfreqSet = c(0, 0),
  CoreClassifier = "LinSVM"
)
}
\arguments{
\item{task_list}{List of target task feature matrices.}

\item{alg1_result}{Output list from \code{\link{alg1_baselineClass}}.}

\item{alg2_result}{Output list from \code{\link{alg2_rob_meanNCov}}.}

\item{alg3_result}{Vector of bias updates from \code{\link{alg3_shiftComp}}.}

\item{alg4_result}{Data frame of bias-adjusted intercepts from \code{\link{alg4_BiasUpdate}}.}

\item{X_feat_cols}{Optional feature indices to retain.}

\item{Marg}{Margin parameter controlling density evaluation width.}

\item{save2file}{Logical; save diagnostic figures.}

\item{ADM}{Logical; apply additional data manipulations.}

\item{datatyp}{Character flag indicating data type (e.g., \code{"FC"}).}

\item{RCSmodeBL}{Logical; enable rare-cell subset adjustments.}

\item{RCSfreqSet}{Numeric vector specifying frequency thresholds for RCS mode.}

\item{CoreClassifier}{Core classifier identifier (e.g., \code{"LinSVM"}).}
}
\value{
A list containing updated normal vectors, smoothing metadata, and key optima
explored during the search.
}
\description{
Rotates the decision boundary normal vectors toward directions that maximize
density support on target tasks. The routine evaluates candidate steps along
principal directions and returns updated hyperplanes.
}
