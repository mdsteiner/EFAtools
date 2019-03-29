#' spss_6f
#'
#' Various outputs from SPSS FACTOR used in Grieder & Steiner (2019). See below for details
#'
#'
#' @format A list of 8
#' \describe{
#'   \item{final_communalities}{(vector) - The final communalities obtained with the FACTOR algorithm with PAF and no rotation. For details see Grieder and Grob (2019).}
#'   \item{final_eigenvalues}{(vector) - Final eigenvalues obtained with the FACTOR algorithm with PAF and no rotation.}
#'   \item{loadings_unrotated}{(matrix) - F1 to F6 = unrotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names. For details see Grieder and Grob (2019).}
#'   \item{paf_iter}{(numeric) - Number of iterations needed for the principal axis factoring to converge.}
#'   \item{loadings_varimax}{(matrix) - F1 to F6 = varimax rotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names. For details see Grieder and Grob (2019).}
#'   \item{loadings_promax}{(matrix) - F1 to F6 = promax rotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names. For details see Grieder and Grob (2019).}
#'   \item{structure_promax}{(matrix) - F1 to F6  = promax rotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names. For details see Grieder and Grob (2019).}
#'   \item{phi_promax}{(matrix) - F1 to F6 = intercorrelations of the promax rotated loadings.}
#'  }
#' @source Grieder, S., & Steiner, M.D. (2019). Resolving differences in principal axis factor analysis and promax rotation in SPSS and R psych. Submitted Manuscript.
"spss_6f"
