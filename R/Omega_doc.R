#' omega_Watkins
#'
#' @description Output of the omega computation performed in Grieder and Steiner (2019),
#' using Watkins' omega program (Watkins, 2013) on the Schmid-Leiman
#' solution of the IDS-2 data.
#'
#'
#' @format A 6 x 3 data.frame
#' \describe{
#'   \item{omega_tot}{(numeric) - Omega total of the factors. For details see \code{\link{OMEGA}} and Grieder and Steiner (2019).}
#'   \item{omega_h}{(numeric) - Omega hierarchical of the general factor. For details see \code{\link{OMEGA}} Grieder and Steiner (2019).}
#'   \item{omega_sub}{(numeric) - Omega subscale of the general factor. For details see \code{\link{OMEGA}} Grieder and Steiner (2019).}
#'  }
#' @source Grieder, S., & Steiner, M.D.(2019). Algorithmic Jingle Jungle: Comparison of Implementations of an EFA Procedure in R psych Versus SPSS, MacOrtho, and Omega. Submitted Manuscript.
#' @source Watkins, M. W. (2013). Omega [Computer Software]. Phoenix, AZ: Ed & Psych Associates.
"omega_Watkins"
