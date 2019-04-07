#' spss_7f
#'
#' Various outputs from SPSS FACTOR used in Grieder & Steiner (2019). See below for details
#'
#'
#' @format A list of 8
#' \describe{
#'   \item{final_communalities}{(vector) - The final communalities obtained with the FACTOR algorithm with PAF and no rotation. For details see Grieder and Grob (2019).}
#'   \item{final_eigenvalues}{(vector) - Final eigenvalues obtained with the FACTOR algorithm with PAF and no rotation.}
#'   \item{loadings_unrotated}{(matrix) - F1 to F7 = unrotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names. For details see Grieder and Grob (2019).}
#'   \item{paf_iter}{(numeric) - Number of iterations needed for the principal axis factoring to converge.}
#'   \item{loadings_varimax}{(matrix) - F1 to F7 = varimax rotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names. For details see Grieder and Grob (2019).}
#'   \item{loadings_promax}{(matrix) - F1 to F7 = promax rotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names. For details see Grieder and Grob (2019).}
#'   \item{structure_promax}{(matrix) - F1 to F7 = promax rotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names. For details see Grieder and Grob (2019).}
#'   \item{phi_promax}{(matrix) - F1 to F7 = intercorrelations of the promax rotated loadings.}
#'  }
#' @source Grieder, S., & Steiner, M.D.(2019). Algorithmic Jingle Jungle: Comparison of Implementations of an EFA Procedure in R psych Versus SPSS, MacOrtho, and Omega. Submitted Manuscript.
#' @source Grieder, S., & Grob, A. (in Press). Exploratory factor analysis of the intelligence and development scales–2 (IDS-2): Implications for theory and practice. Assessment.
"spss_7f"
