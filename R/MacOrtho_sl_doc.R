#' MacOrtho
#'
#' Output of the Schmid Leiman transformation performed in Grieder & Steiner (2019), using the MacOrtho (Watkins, 2004) program.
#'
#'
#' @format A list of 2
#' \describe{
#'   \item{sl}{(matrix) - g = General / second order factor of the Schmid-Leiman solution. F1 to F5  = First order factors of the Schmid-Leiman solution. This Schmid-Leiman solution was found using the MacOrtho program (Watkins, 2004).}
#'   \item{L2}{(matrix) - Second order loadings used for the Schmid-Leiman transformation. This Schmid-Leiman solution was found using the MacOrtho program (Watkins, 2004).}
#'  }
#' @source Grieder, S., & Steiner, M.D.(2019). Resolving differences in principal axis factor analysis and promax rotation in SPSS and R psych. Submitted Manuscript.
#' @source Watkins, M. W. (2004). Macortho [Computer Software]. Phoenix, AZ: EdPsych Associates, Inc.
"MacOrtho"
