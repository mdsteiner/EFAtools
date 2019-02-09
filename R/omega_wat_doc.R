#' omega_Watkins
#'
#' Output of the omega computation performed in Steiner & Grieder (2019), using Watkins' omega program (Watkins, 2013) on the Schmid-Leiman solution of the \code{\link[EFAdiff:IDS2_R]{IDS-2}}} data.
#'
#'
#' @format A 6 x 3 data.frame
#' \describe{
#'   \item{omega_tot}{(numeric) - Omega total of the factors. For details see \code{\link{OMEGA}} and Steiner and Grieder (2019).}
#'   \item{omega_h}{(numeric) - Omega hierarchical of the general factor. For details see \code{\link{OMEGA}} Steiner and Grieder (2019).}
#'   \item{omega_sub}{(numeric) - Omega subscale of the general factor. For details see \code{\link{OMEGA}} Steiner and Grieder (2019).}
#'  }
#' @source Steiner, M.D., & Grieder, S. (2019). Resolving differences in principal axis factor analysis and promax rotation in SPSS and R psych. Submitted Manuscript.
#' @source Watkins, M. W. (2013). Omega [Computer Software]. Phoenix, AZ: Ed & Psych Associates.
"omega_Watkins"
