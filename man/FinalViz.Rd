\name{FinalViz}
\alias{FinalViz}
\title{Visualize RTL hyperplane performance across algorithms}
\usage{
FinalViz(
  TrainTestSet.ls,
  alg1_result,
  alg2_result,
  alg3_result,
  alg4_result,
  alg6_result,
  datatyp,
  ADM
)
}
\arguments{
\item{TrainTestSet.ls}{List of train/test splits with \code{X.test} and \code{Y.test} elements.}

\item{alg1_result}{Output list from \code{\link{alg1_baselineClass}}.}

\item{alg2_result}{Output list from \code{\link{alg2_rob_meanNCov}}.}

\item{alg3_result}{Vector of bias updates from \code{\link{alg3_shiftComp}}.}

\item{alg4_result}{Data frame of bias-adjusted intercepts from \code{\link{alg4_BiasUpdate}}.}

\item{alg6_result}{List from \code{\link{alg6_NormalVectorUpdate}} containing updated weights.}

\item{datatyp}{Character flag indicating data type (e.g., \code{"FC"}).}

\item{ADM}{Logical; apply additional data manipulations.}
}
\value{
A list of ggplot objects, combined statistics, and per-task predictions for each
algorithm stage.
}
\description{
Compares baseline, bias-updated, and normal-vector–updated hyperplanes on
held-out target tasks, generating summary statistics and ggplot visualizations
for downstream reporting.
}
