% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.EFA_AVERAGE.R
\name{plot.EFA_AVERAGE}
\alias{plot.EFA_AVERAGE}
\title{Plot EFA_AVERAGE object}
\usage{
\method{plot}{EFA_AVERAGE}(x, ...)
}
\arguments{
\item{x}{list. An output from the \link{EFA_AVERAGE} function.}

\item{...}{not used.}
}
\description{
Plot method showing a summarized output of the \link{EFA_AVERAGE} function
}
\examples{
\dontrun{
EFA_aver <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500)
EFA_aver
}

}
