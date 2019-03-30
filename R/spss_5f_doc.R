#' spss_5f
#'
#' Various outputs from SPSS FACTOR used in Grieder and Steiner (2019). See below for details
#'
#'
#' @format A list of 10
#' \describe{
#'   \item{final_communalities}{(vector) - The final communalities obtained with the FACTOR algorithm with PAF and no rotation. For details see Grieder and Grob (2019).}
#'   \item{final_eigenvalues}{(vector) - Final eigenvalues obtained with the FACTOR algorithm with PAF and no rotation.}
#'   \item{loadings_unrotated}{(matrix) - F1 to F5 = unrotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names. For details see Grieder and Grob (2019).}
#'   \item{paf_iter}{(numeric) - Number of iterations needed for the principal axis factoring to converge.}
#'   \item{loadings_varimax}{(matrix) - F1 to F5 = varimax rotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names. For details see Grieder and Grob (2019).}
#'   \item{loadings_promax}{(matrix) - F1 to F5 = promax rotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names. For details see Grieder and Grob (2019).}
#'   \item{structure_promax}{(matrix) - F1 to F5  = promax rotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names. For details see Grieder and Grob (2019).}
#'   \item{phi_promax}{(matrix) - F1 to F5  = intercorrelations of the promax rotated loadings.}
#'   \item{sl}{(matrix) - g = General / second order factor of the Schmid-Leiman solution. F1 to F5  = First order factors of the Schmid-Leiman solution. h2 = Communalities of the Schmid-Leiman solution. This Schmid-Leiman solution was found using the SPSS Syntax provided by Wolff and Preising (2005).}
#'   \item{L2}{(matrix) - Second order loadings used for the Schmid-Leiman transformation. This Schmid-Leiman solution was found using the SPSS Syntax provided by Wolff and Preising (2005).}
#'  }
#' @source Grieder, S., & Steiner, M.D.(2019). Algorithmic Jingle Jungle: Comparison of Implementations of an EFA Procedure in R psych Versus SPSS, MacOrtho, and Omega. Submitted Manuscript.
#' @source Wolff, H.G., & Preising, K. (2005). Exploring item and higher order factor structure with the schmid-leiman solution: Syntax codes for spss and sas. Behavior Research Methods, 37, 48–58. doi: 10.3758/BF03206397
#' @source Grieder, S., & Grob, A. (in Press). Exploratory factor analysis of the intelligence and development scales–2 (IDS-2): Implications for theory and practice. Assessment.
"spss_5f"
