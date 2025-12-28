\name{alg4_BiasUpdate}
\alias{alg4_BiasUpdate}
\title{Update classification bias on target tasks}
\usage{
alg4_BiasUpdate(
  task_list,
  alg1_result,
  alg2_result,
  alg3_result,
  goodColumns,
  save2file,
  Marg,
  alg4MinFx,
  ADM = FALSE,
  useMedian = TRUE,
  ZnormMappingBL = FALSE,
  datatyp = "FC",
  RCSmodeBL = FALSE,
  RCSfreqSet = c(0, 0),
  CoreClassifier = "LinSVM",
  use_parallel = FALSE,
  parallel_cores = NULL,
  wide_data_threshold = 200
)
}
\arguments{
\item{task_list}{List of target task feature matrices.}

\item{alg1_result}{Output list from \code{\link{alg1_baselineClass}}.}

\item{alg2_result}{Output list from \code{\link{alg2_rob_meanNCov}}.}

\item{alg3_result}{Vector of bias updates from \code{\link{alg3_shiftComp}}.}

\item{goodColumns}{Optional column indices to retain.}

\item{save2file}{Logical; persist diagnostic plots to disk.}

\item{Marg}{Margin parameter controlling neighborhood width.}

\item{alg4MinFx}{Mode selection for minima detection (\code{"gd"}, \code{"mean"}, \code{"win"}).}

\item{ADM}{Logical; apply additional data manipulations.}

\item{useMedian}{Logical; smooth minima search using medians.}

\item{ZnormMappingBL}{Logical; use z-normalized mapping for bias updates.}

\item{datatyp}{Character flag indicating data type (e.g., \code{"FC"}).}

\item{RCSmodeBL}{Logical; enable rare-cell subset adjustments.}

\item{RCSfreqSet}{Numeric vector specifying frequency thresholds for RCS mode.}

\item{CoreClassifier}{Core classifier identifier (e.g., \code{"LinSVM"}).}

\item{use_parallel}{Logical; parallelize per-task updates when \code{TRUE}.}

\item{parallel_cores}{Optional integer overriding detected core count.}

\item{wide_data_threshold}{Integer; coerce wide feature sets via \code{data.table} to reduce copies.}
}
\value{
A list with updated predictions, density summaries, and adjusted intercepts for
each target task.
}
\description{
Implements Algorithm 4 to refine decision boundaries on target tasks by
estimating density-informed bias adjustments. Supports optional data
manipulations and RCS (rare-cell subset) focusing.
}
