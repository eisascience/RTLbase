\name{rtl_benchmark}
\alias{rtl_benchmark}
\title{Benchmark RTLbase training and update routines}
\usage{
rtl_benchmark(
  sizes = c(small = 250, medium = 1500, large = 3200),
  features = c(6, 24, 64),
  repetitions = 2,
  parallel_cores = NULL,
  margin = 0.5,
  seed = 99
)
}
\arguments{
\item{sizes}{Named integer vector of dataset sizes to benchmark.}

\item{features}{Integer vector of feature counts aligned to \code{sizes} (recycled when length-one).}

\item{repetitions}{Number of repeated runs per dataset size and mode.}

\item{parallel_cores}{Optional integer overriding detected cores for parallel runs.}

\item{margin}{Margin value forwarded to \code{\link[=alg4_BiasUpdate]{alg4_BiasUpdate()}}.}

\item{seed}{Integer seed for reproducible synthetic data generation.}
}
\value{
A \code{data.table} summarizing elapsed seconds per dataset and mode.
}
\description{
Runs the baseline classifier and downstream updates on synthetic datasets of
varying sizes to quantify speedups from the new vectorized and parallel
execution paths.
}
