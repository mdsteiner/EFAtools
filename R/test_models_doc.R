#' test_models
#'
#' Correlation matrices created from simulated data from four of the
#' \code{population_models} cases, each with strong factor intercorrealtions.
#' These are used to compare the psych and SPSS implementations in Grieder and
#' Steiner (2019) with the actual implementations. For details on the cases, see
#' \link{population_models}.
#'
#'
#' @format A list of 4 lists "baseline", "case_1a", "case_6b", and"case_11b", each with the following elements.
#' \describe{
#'   \item{cormat}{(matrix) - The correlation matrix of the simulated data.}
#'   \item{n_factors}{(numeric) - The true number of factors.}
#'   \item{N}{(numeric) - The sample size of the generated data.}
#'  }
#' @source Grieder, S., & Steiner, M.D.(2019). Algorithmic Jingle Jungle: Comparison of Implementations of an EFA Procedure in R psych Versus SPSS, MacOrtho, and Omega. Submitted Manuscript.
"test_models"
