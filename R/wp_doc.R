#' wp
#'
#' Loadings and Schmid Leiman solutions reported in Wolff and Preising (2005).
#'
#'
#' @format A list of 5
#' \describe{
#'   \item{L1}{(matrix) - First order loadings (pattern matrix), reported in the Appendix of Wolff and Preising (2005). Note that the loadings reported in Table 2 of Wolff and Preising (2005) contain typos and are less precise.}
#'   \item{Phi1}{(matrix) - First order factor intercorrelations reported in Wolff and Preising (2005)}
#'   \item{sls_L1}{(matrix) - First order Schmid Leiman solution reported in Wolff and Preising (2005). Note that the solution reported in the paper contains a typo in the third item of factor three, correct loadings can be found in the Appendix, table A2.}
#'   \item{sls_L2}{(matrix) - Second order Schmid Leiman solution reported in Wolff and Preising (2005). Note that the solution reported in the paper contains a typo in the third item of the general factor, correct loadings can be found in the Appendix, table A2.}
#'   \item{sls_MacOrtho}{(data.frame) - subtest = the nine subtest names. General = second order Schmid Leiman solution obtained with MacOrtho (Watkins, 2004). 1 to 3 = First order Schmid Leiman solution obtained with MacOrtho (Watkins, 2004).}
#'  }
#' @source Steiner, M.D., & Grieder, S. (2019). Resolving differences in principal axis factor analysis and promax rotation in SPSS and R psych. Submitted Manuscript.
#' Wolff, H.G & Preising, K. (2005). Exploring item and higher order factor structure with the Schmid–Leiman solution: Syntax codes for SPSS and SAS. Behavior Research Methods, 37 (1), 48-58.
"wp"
