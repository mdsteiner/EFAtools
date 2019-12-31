#' population_models
#'
#' Pupulation factor models used for the model recovery analyses reported in
#' Grieder and Steiner (2019). All combinations of the loadings and the factor
#' intercorrelations were used in the simulations. Many models are based on cases
#' used in de Winter and Dodou (2012).
#'
#'
#' @format A list of 2 lists "loadings" and "phis".
#' \describe{
#' \code{loadings} contains the following matrices:
#'   \item{baseline}{(matrix) - The loadings of the baseline model. Three factors with six indicators each, all with loadings of .6. Same baseline model as used in de Winter and Dodou (2012).}
#'   \item{case_1a}{(matrix) - Small number of indicators per factor, in this case 2.}
#'   \item{case_1b}{(matrix) - Small number of indicators per factor, in this case 3. Case 5 in de Winter and Dodou (2012).}
#'   \item{case_1c}{(matrix) - Small number of indicators per factor, in this case 4.}
#'   \item{case_1d}{(matrix) - Small number of indicators per factor, in this case 5.}
#'   \item{case_2}{(matrix) - Low factor loadings of .3.}
#'   \item{case_3}{(matrix) - High factor loadings of .9.}
#'   \item{case_4}{(matrix) - Different loadings between factors (one factor with .9, .6, and .3, respectively). Case 7 in de Winter and Dodou (2012).}
#'   \item{case_5}{(matrix) - Different loadings within factors (each factor has two loadings of each .9, .6, and .3). Similar to cases 8/ 9 in de Winter and Dodou (2012).}
#'   \item{case_6a}{(matrix) - 1 crossloading of .4. Similar to case 10 in de Winter and Dodou (2012).}
#'   \item{case_6b}{(matrix) - 3 crossloading of .4 (One factor with 2 and one with 1 crossloading). Similar to case 10 in de Winter and Dodou (2012).}
#'   \item{case_7}{(matrix) - Different number of indicators per factor (2, 4, and 6 respectively). Similar to cases 11/ 12 in de Winter and Dodou (2012).}
#'   \item{case_8}{(matrix) - Random variation in loadings added, drawn from a uniform distribution between [-.2, .2]. Case 13 in de Winter and Dodou (2012).}
#'   \item{case_9a}{(matrix) - Small number of indicators per factor, in this case 2, with different loadings within one of the factors.}
#'   \item{case_9b}{(matrix) - Small number of indicators per factor, in this case 3, with different loadings.}
#'   \item{case_9c}{(matrix) - Small number of indicators per factor, in this case 4, with different loadings.}
#'   \item{case_9d}{(matrix) - Small number of indicators per factor, in this case 5, with different loadings.}
#'   \item{case_10a}{(matrix) - Six factors, all with loadings of .6. Small number of indicators per factor, in this case 2.}
#'   \item{case_10b}{(matrix) - Six factors, all with loadings of .6. Small number of indicators per factor, in this case 3.}
#'   \item{case_10c}{(matrix) - Six factors, all with loadings of .6. Small number of indicators per factor, in this case 4.}
#'   \item{case_10d}{(matrix) - Six factors, all with loadings of .6. Small number of indicators per factor, in this case 5.}
#'   \item{case_10e}{(matrix) - Six factors, all with loadings of .6. Each with 6 indicators.}
#'   \item{case_11a}{(matrix) - Six factors, with different loadings within and between. Small number of indicators per factor, in this case 2.}
#'   \item{case_11b}{(matrix) - Six factors, with different loadings within and between. Small number of indicators per factor, in this case 3.}
#'   \item{case_11c}{(matrix) - Six factors, with different loadings within and between. Small number of indicators per factor, in this case 4.}
#'   \item{case_11d}{(matrix) - Six factors, with different loadings within and between. Small number of indicators per factor, in this case 5.}
#'   \item{case_11e}{(matrix) - Six factors, with different loadings within and between. Each with 6 indicators.}
#'   \code{phis_3} contains the following 3x3 matrices:
#'   \item{zero}{(matrix) - Matrix of factor intercorrelations of 0. Same intercorrelations as used in de Winter and Dodou (2012).}
#'   \item{moderate}{(matrix) - Matrix of moderate factor intercorrelations of .3.}
#'   \item{mixed}{(matrix) - Matrix of mixed (.3, .5, and .7) factor intercorrelations.}
#'   \item{strong}{(matrix) - Matrix of strong factor intercorrelations of .7. Same intercorrelations as used in de Winter and Dodou (2012).}
#'   \code{phis_6} contains the following 6x6 matrices:
#'   \item{zero}{(matrix) - Matrix of factor intercorrelations of 0. Same intercorrelations as used in de Winter and Dodou (2012).}
#'   \item{moderate}{(matrix) - Matrix of moderate factor intercorrelations of .3.}
#'   \item{mixed}{(matrix) - Matrix of mixed (around .3, .5, and .7; smoothing was necessary for the matrix to be positive definite) factor intercorrelations.}
#'   \item{strong}{(matrix) - Matrix of strong factor intercorrelations of .7. Same intercorrelations as used in de Winter and Dodou (2012).}
#'  }
#' @source Grieder, S., & Steiner, M.D.(2019). Algorithmic Jingle Jungle: Comparison of Implementations of an EFA Procedure in R psych Versus SPSS, MacOrtho, and Omega. Submitted Manuscript.
#' @source de Winter, J.C.F., & Dodou, D. (2012). Factor recovery by principal axis factoring and maximum likelihood factor analysis as a function of factor pattern and sample size. Journal of Applied Statistics. 39.
"population_models"
