\name{secondary_tools}
\alias{ClosingWindowAproxFunOptimizedMinima}
\alias{SpecializedGD2D}
\alias{find_peaks}
\alias{SmartMinimaAk}
\alias{CriticalPointsX}
\alias{inflect}
\alias{argmax}
\title{Secondary utilities for minima/peak finding}
\usage{
ClosingWindowAproxFunOptimizedMinima(xyDF, NWindows = 50)

SpecializedGD2D(
  xyDF,
  quants = c(0.1, 0.9),
  MinMaxRange = FALSE,
  learningRate = 0.1,
  maximumIters4GD = 100,
  print2screen = FALSE,
  plot2screen = FALSE
)

find_peaks(x, m = 4)

SmartMinimaAk(
  SmoothedXY,
  print2screen = FALSE,
  datatyp = "FC",
  learnRate = 0.1,
  mode = "NZ",
  cleanQuantsX = TRUE,
  cleanQuantsY = FALSE
)

CriticalPointsX(InputNumVec, n = 3, plot2screen = TRUE)

inflect(x, threshold = 1)

argmax(x, y, w = 1, ...)
}
\arguments{
\item{xyDF}{Data frame with \code{x} and \code{y} columns used by
\code{ClosingWindowAproxFunOptimizedMinima()} and \code{SpecializedGD2D()}.}

\item{NWindows}{Integer number of window pairs to evaluate (closing-window minima).}

\item{quants}{Optional quantile bounds used to trim the search region (gradient descent).}

\item{MinMaxRange}{Optional numeric range restricting \code{x} values (gradient descent).}

\item{learningRate}{Learning rate for \code{gradDescent::SAGD} (gradient descent).}

\item{maximumIters4GD}{Maximum number of iterations for gradient descent.}

\item{print2screen}{Logical; print progress updates.}

\item{plot2screen}{Logical; produce diagnostic plots.}

\item{x}{Numeric vector for \code{find_peaks()} or the x-axis input to \code{argmax()}.}

\item{m}{Neighborhood size for peak prominence.}

\item{SmoothedXY}{Data frame of smoothed \code{x/y} coordinates for \code{SmartMinimaAk()}.}

\item{datatyp}{Character flag indicating data type (e.g., \code{"FC"}).}

\item{learnRate}{Learning rate passed to helper routines inside \code{SmartMinimaAk()}.}

\item{mode}{Minima selection strategy for \code{SmartMinimaAk()}.}

\item{cleanQuantsX}{Logical; trim extreme \code{x} quantiles before searching.}

\item{cleanQuantsY}{Logical; trim low \code{y} density values before searching.}

\item{InputNumVec}{Numeric vector to evaluate in \code{CriticalPointsX()}.}

\item{n}{Integer sensitivity threshold for inflection detection.}

\item{threshold}{Integer window size for detecting changes in slope (\code{inflect()}).}

\item{y}{Numeric vector of y-values supplied to \code{argmax()}.}

\item{w}{Window size supplied to \code{zoo::rollapply} in \code{argmax()}.}

\item{...}{Additional arguments passed to \code{stats::loess} in \code{argmax()}.}
}
\value{
Each function returns helper minima/peak locations or supporting statistics for
Algorithm 4 and Algorithm 6 workflows.
}
\description{
Utilities used throughout the RTL pipeline for finding minima and peaks,
smoothing density curves, and extracting critical points from noisy vectors.
}
