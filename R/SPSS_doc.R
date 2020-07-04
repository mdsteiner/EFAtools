#' Various outputs from SPSS FACTOR
#'
#' Various outputs from SPSS FACTOR for the IDS-2 (Grob & Hagmann-von Arx, 2018), the WJIV (3 to 5 and 20 to 39 years; McGrew, LaForte, & Schrank, 2014), the DOSPERT (Frey et al., 2017; Weber,
#' Blais, & Betz, 2002), the NEO-PI-R (Costa, & McCrae, 1992), and four simulated datasets (baseline, case_1a, case_6b, and case_11b, see \link{test_models} and \link{population_models}) used in Grieder and Steiner (2020).
#'
#'
#' @format A list of 9 containing EFA results for each of the data sets mentioned above. Each of these nine entries is a list of 4 or 8 (see details), of the following structure:
#' \describe{
#'   \item{paf_comm}{(vector) - The final communalities obtained with the FACTOR algorithm with PAF and no rotation. For details, see Grieder and Grob (2019).}
#'   \item{paf_load}{(matrix) - F1 to FN = unrotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names.}
#'   \item{paf_iter}{(numeric) - Number of iterations needed for the principal axis factoring to converge.}
#'   \item{var_load}{(matrix) - F1 to FN = varimax rotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names.}
#'   \item{pro_load}{(matrix) - F1 to FN = promax rotated factor loadings obtained with the FACTOR algorithm with PAF. Rownames are the abbreviated subtest names.}
#'   \item{pro_phi}{(matrix) - F1 to FN = intercorrelations of the promax rotated loadings.}
#'   \item{sl}{(matrix) - g = General / second order factor of the Schmid-Leiman solution. F1 to FN  = First order factors of the Schmid-Leiman solution. h2 = Communalities of the Schmid-Leiman solution. This Schmid-Leiman solution was found using the SPSS Syntax provided by Wolff and Preising (2005).}
#'   \item{L2}{(matrix) - Second order loadings used for the Schmid-Leiman transformation. This Schmid-Leiman solution was found using the SPSS Syntax provided by Wolff and Preising (2005).}
#'  }
#' @details The IDS-2, the two WJIV, the DOSPERT, and the NEO-PI-R contain all the above entries, while the four simulated datasets contain only paf_load, var_load, pro_load, and pro_phi.
#'
#' @source Grieder, S., & Steiner, M.D. (2020). Algorithmic Jingle Jungle:
#' A Comparison of Implementations of Principal Axis Factoring and Promax Rotation
#'  in R and SPSS. Manuscript in Preparation.
#' @source Wolff, H.G., & Preising, K. (2005). Exploring item and higher order factor structure with the Schmid-Leiman solution: Syntax codes for SPSS and SAS. Behavior Research Methods, 37, 48–58. doi: 10.3758/BF03206397
#' @source Grieder, S., & Grob, A. (2019). Exploratory factor analyses of the intelligence and development scales–2: Implications for theory and practice. Assessment. Advance online publication. doi:10.1177/10731911198450
#' @source Grob, A., & Hagmann-von Arx, P. (2018). Intelligence and Development Scales--2 (IDS-2). Intelligenz- und Entwicklungsskalen für Kinder und Jugendliche.
#' [Intelligence and Development Scales for Children and Adolescents.]. Bern, Switzerland: Hogrefe.
#' @source Frey, R., Pedroni, A., Mata, R., Rieskamp, J., & Hertwig, R. (2017). Risk preference shares the psychometric structure of major psychological traits. Science Advances, 3, e1701381.
#' @source McGrew, K. S., LaForte, E. M., & Schrank, F. A. (2014). Technical
#'  Manual. Woodcock-Johnson IV. Rolling Meadows, IL: Riverside.
#' @source Schrank, F. A., McGrew, K. S., & Mather, N. (2014). Woodcock-Johnson IV.
#' Rolling Meadows, IL: Riverside.
#' @source Costa, P. T., & McCrae, R. R. (1992). NEO PI-R professional manual. Odessa, FL: Psychological Assessment Resources, Inc.
"SPSS"
