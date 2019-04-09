#' IDS2_7f
#'
#' Outputs from SPSS FACTOR for the IDS-2 seven-factor solution used in Grieder & Steiner (2019).
#'
#'
#' @format A list of 6
#' \describe{
#'   \item{paf_comm}{(vector) - The final communalities obtained with the FACTOR algorithm with PAF and no rotation. For details see Grieder and Grob (2019).}
#'   \item{paf_load}{(matrix) - F1 to F7 = unrotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names. For details see Grieder and Grob (2019).}
#'   \item{paf_iter}{(numeric) - Number of iterations needed for the principal axis factoring to converge.}
#'   \item{var_load}{(matrix) - F1 to F7 = varimax rotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names. For details see Grieder and Grob (2019).}
#'   \item{pro_load}{(matrix) - F1 to F7 = promax rotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names. For details see Grieder and Grob (2019).}
#'   \item{pro_phi}{(matrix) - F1 to F7 = intercorrelations of the promax rotated loadings.}
#'  }
#' @source Grieder, S., & Steiner, M.D.(2019). Algorithmic Jingle Jungle: Comparison of Implementations of an EFA Procedure in R psych Versus SPSS, MacOrtho, and Omega. Submitted Manuscript.
#' @source Grieder, S., & Grob, A. (in Press). Exploratory factor analysis of the intelligence and development scales–2 (IDS-2): Implications for theory and practice. Assessment.
"IDS2_7f"
