#' SPSS
#'
#' Various outputs from SPSS FACTOR for the IDS-2, the WJIV (3 to 5 and 20 to 39 years), the DOSPERT, and the NEO-PI-R used in Grieder and Steiner (2019).
#'
#'
#' @format A list of 5
#' \describe{
#'   \item{IDS_2}{(list) - A list of 8}
#'    \describe{
#'   \item{paf_comm}{(vector) - The final communalities obtained with the FACTOR algorithm with PAF and no rotation. For details see Grieder and Grob (2019).}
#'   \item{paf_load}{(matrix) - F1 to F7 = unrotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names. For details see Grieder and Grob (2019).}
#'   \item{paf_iter}{(numeric) - Number of iterations needed for the principal axis factoring to converge.}
#'   \item{var_load}{(matrix) - F1 to F7 = varimax rotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names. For details see Grieder and Grob (2019).}
#'   \item{pro_load}{(matrix) - F1 to F7 = promax rotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names. For details see Grieder and Grob (2019).}
#'   \item{pro_phi}{(matrix) - F1 to F7 = intercorrelations of the promax rotated loadings.}
#'   \item{sl}{(matrix) - g = General / second order factor of the Schmid-Leiman solution. F1 to F5  = First order factors of the Schmid-Leiman solution. h2 = Communalities of the Schmid-Leiman solution. This Schmid-Leiman solution was found using the SPSS Syntax provided by Wolff and Preising (2005).}
#'   \item{L2}{(matrix) - Second order loadings used for the Schmid-Leiman transformation. This Schmid-Leiman solution was found using the SPSS Syntax provided by Wolff and Preising (2005).}
#'  }
#'  \item{WJIV3_5}{(list) - A list of 8}
#'  \item{WJIV20_39}{(list) - A list of 8}
#'  \item{DOSPERT}{(list) - A list of 8}
#'  \item{NEO}{(list) - A list of 6}
#'  }
#'
#' @source Grieder, S., & Steiner, M.D.(2019). Algorithmic Jingle Jungle: Comparison of Implementations of an EFA Procedure in R psych Versus SPSS, MacOrtho, and Omega. Submitted Manuscript.
#' @source Wolff, H.G., & Preising, K. (2005). Exploring item and higher order factor structure with the schmid-leiman solution: Syntax codes for spss and sas. Behavior Research Methods, 37, 48–58. doi: 10.3758/BF03206397
#' @source Grieder, S., & Grob, A. (in Press). Exploratory factor analysis of the intelligence and development scales–2 (IDS-2): Implications for theory and practice. Assessment.
"SPSS"
