\name{alg2_functions}
\alias{alg2_rob_meanNCov}
\alias{alg2_CalcMore}
\title{Estimate and refine robust hyperplane moments}
\usage{
alg2_rob_meanNCov(alg1_result_baselineSVM, print2screen = FALSE)

alg2_CalcMore(alg2_res)
}
\arguments{
\item{alg1_result_baselineSVM}{Matrix of hyperplane parameters from
\code{\link{alg1_baselineClass}} where each row represents a source model.}

\item{print2screen}{Logical; print progress updates.}

\item{alg2_res}{List output from \code{\link{alg2_rob_meanNCov}}.}
}
\value{
Both functions return augmented lists containing robust/simple means,
covariance matrices, normalized vectors, and helper quantities used by
downstream algorithms.
}
\description{
\code{alg2_rob_meanNCov()} computes robust mean and covariance estimates for
baseline SVM hyperplanes, falling back to empirical moments when robust
estimation fails. \code{alg2_CalcMore()} post-processes those values to derive
normalized hyperplane parameters, eigenvalues, and helper vectors for later RTL
steps.
}
