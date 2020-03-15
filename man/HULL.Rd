% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HULL.R
\name{HULL}
\alias{HULL}
\title{Hull method for determining the number of factors to retain}
\source{
Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. (2011). The Hull method for selecting the number of common factors. Multivariate behavioral research, 46(2), 340-364.
}
\usage{
HULL(x, cors = TRUE, n_cases = NA, n_factors = NA, ...)
}
\arguments{
\item{x}{matrix or data.frame. Correlation matrix or raw data.}

\item{cors}{logical. Whether x is a correlation matrix. Default is set to TRUE.}

\item{n_cases}{numeric. Number of cases in the data. This is passed to \link{PARALLEL}.
Only has to be specified if x is a correlation matrix.}

\item{n_factors}{numeric. Theoretical number of factors to retain. The maximum
of this number and the number of factors suggested by \link{PARALLEL} olus
onewill be used in the Hull method.}

\item{...}{Further arguments passed to \link{PARALLEL}.}
}
\value{
A list of class HULL containing the following objects
\item{n_factors}{The number of factors to retain according to the Hull method.}
\item{solutions}{A matrix containing the goodness of fit indices (CAF), degrees of freedom, and for the factors lying on the hull, the st values of the hull solution (see Lorenzo-Seva, Timmereman, and Kiers 2011 for details).}
}
\description{
Implementation of the Hull method suggested by Lorenzo-Seva, Timmereman,
and Kiers (2011), but using principal axis factoring. See details for
parallelization.
}
\details{
The \link{PARALLEL} function and the principal axis factoring of the
  different number of factors can be parallelized using the future framework,
  by calling the \link{future}{plan} function. The examples provide example code
  on how to enable parallel processing.
}
\examples{
\dontrun{
HULL(IDS2_R, n_cases = 2000)

# using parallel processing (Note: plans can be adapted, see the future
# package for details)
future::plan(future::multisession)
HULL(IDS2_R, n_cases = 2000)
}
}