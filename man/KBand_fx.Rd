\name{KBand_fx}
\alias{KBand_fx}
\title{Kernel bandwidth estimation for Algorithm 5}
\usage{
KBand_fx(s_k, c_k)
}
\arguments{
\item{s_k}{Vector of grid locations.}

\item{c_k}{Vector of counts or weights at each grid location.}
}
\value{
Numeric bandwidth estimate.
}
\description{
Computes a Silverman-inspired bandwidth to smooth kernel density estimates of
projection counts used in Algorithm 5.
}
