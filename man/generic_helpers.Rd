\name{generic_helpers}
\alias{ColorTheme}
\alias{Mode}
\alias{IndicatorFX}
\alias{DF2TrainTestls}
\title{General helper utilities for RTLbase}
\usage{
ColorTheme()

Mode(x)

IndicatorFX(x)

DF2TrainTestls(X, Y)
}
\arguments{
\item{x}{Vector of values for \code{Mode()}, or logical flag for \code{IndicatorFX()}.}

\item{X}{Data frame or matrix of features.}

\item{Y}{Vector or factor of class labels aligned to \code{X}.}
}
\value{
\itemize{
\item \code{ColorTheme()} returns a list containing \code{col_vector} and \code{scaleyellowred} palettes.
\item \code{Mode()} returns the most frequent value in \code{x}.
\item \code{IndicatorFX()} returns \code{1} for \code{TRUE} and \code{0} for \code{FALSE}.
\item \code{DF2TrainTestls()} returns a list with \code{TrainXY.ls} and \code{TestXY.ls} stratified splits.
}
}
\description{
Small helper functions for plotting palettes, indicator conversion, computing
modes, and generating stratified train/test splits compatible with RTL
algorithms.
}
