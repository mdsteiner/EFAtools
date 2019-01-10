#' omega_wat
#'
#' Output of the omega computation performed in Steiner & Grieder (2019), using Watkins' omega program (Watkins, 2013) on the Schmid-Leiman solution of the \code{\link[EFAdiff:IDS2_R]{IDS-2}}} data.
#'
#'
#' @format A 8 x 4 data.frame
#' \describe{
#'   \item{Factor}{(character) - Abbreviation of the Factor names. g = general factor. VP = visual processing. PS = processing speed. ASM = auditory short term memory. VSSM = visual spatial short term memory. AR = abstract reasoning. VR = verbal reasoning. LTM = long term memory. For details see \code{\link{IDS2_R}} and Grieder and Grob (2019).}
#'   \item{omega_tot}{(numeric) - Omega total of the factors. For details see \code{\link{OMEGA}} and Steiner and Grieder (2019).}
#'   \item{omega_h}{(numeric) - Omega hierarchical of the general factor. For details see \code{\link{OMEGA}} Steiner and Grieder (2019).}
#'   \item{omega_sub}{(numeric) - Omega subscale of the general factor. For details see \code{\link{OMEGA}} Steiner and Grieder (2019).}
#'  }
#' @source Steiner, M.D., & Grieder, S. (2019). Resolving differences in principal axis factor analysis and promax rotation in SPSS and R psych. Submitted Manuscript.
"omega_wat"
