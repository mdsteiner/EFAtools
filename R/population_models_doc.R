#' population_models
#'
#' Population factor models, some of which (baseline to case_11e) used for the
#' simulation analyses reported in Grieder and Steiner (2019). All combinations
#' of the pattern matrices and the factor
#' intercorrelations were used in the simulations. Many models are based on cases
#' used in de Winter and Dodou (2012).
#'
#'
#' @format A list of 3 lists "loadings", "phis_3", and "phis_6".
#' \describe{
#' \code{loadings} contains the following matrices of pattern coefficients:
#'   \item{baseline}{(matrix) - The pattern coefficients of the baseline model. Three factors with six indicators each, all with pattern coefficients of .6. Same baseline model as used in de Winter and Dodou (2012).}
#'   \item{case_1a}{(matrix) - Three factors with 2 indicators per factor.}
#'   \item{case_1b}{(matrix) - Three factors with 3 indicators per factor. Case 5 in de Winter and Dodou (2012).}
#'   \item{case_1c}{(matrix) - Three factors with 4 indicators per factor.}
#'   \item{case_1d}{(matrix) - Three factors with 5 indicators per factor.}
#'   \item{case_2}{(matrix) - Same as baseline model but with low pattern coefficients of .3.}
#'   \item{case_3}{(matrix) - Same as baseline model but with high pattern coefficients of .9.}
#'   \item{case_4}{(matrix) - Three factors with different pattern coefficients \emph{between} factors (one factor with .9, one with .6, and one with .3, respectively). Case 7 in de Winter and Dodou (2012).}
#'   \item{case_5}{(matrix) - Three factors with different pattern coefficients \emph{within} factors (each factor has two pattern coefficients of each .9, .6, and .3). Similar to cases 8/ 9 in de Winter and Dodou (2012).}
#'   \item{case_6a}{(matrix) - Same as baseline model but with one cross loading of .4. Similar to case 10 in de Winter and Dodou (2012).}
#'   \item{case_6b}{(matrix) - Same as baseline model but with three cross loading of .4 (One factor with 2 and one with 1 crossloading). Similar to case 10 in de Winter and Dodou (2012).}
#'   \item{case_7}{(matrix) - Three factors with different number of indicators per factor (2, 4, and 6 respectively). Similar to cases 11/ 12 in de Winter and Dodou (2012).}
#'   \item{case_8}{(matrix) - Three factors with random variation in pattern coefficients added, drawn from a uniform distribution between [-.2, .2]. Case 13 in de Winter and Dodou (2012).}
#'   \item{case_9a}{(matrix) - Three factors with 2 indicators per factor, with different pattern coefficients within one of the factors.}
#'   \item{case_9b}{(matrix) - Three factors with 3 indicators per factor, with different pattern coefficients.}
#'   \item{case_9c}{(matrix) - Three factors with 4 indicators per factor, with different pattern coefficients.}
#'   \item{case_9d}{(matrix) - Three factors with 5 indicators per factor, with different pattern coefficients.}
#'   \item{case_10a}{(matrix) - Six factors with 2 indicators per factor, all with pattern coefficients of .6.}
#'   \item{case_10b}{(matrix) - Six factors with 3 indicators per factor, all with pattern coefficients of .6.}
#'   \item{case_10c}{(matrix) - Six factors with 4 indicators per factor, all with pattern coefficients of .6.}
#'   \item{case_10d}{(matrix) - Six factors with 5 indicators per factor, all with pattern coefficients of .6.}
#'   \item{case_10e}{(matrix) - Six factors with 6 indicators per factor, all with pattern coefficients of .6.}
#'   \item{case_11a}{(matrix) - Six factors with 2 indicators per factor, with different pattern coefficients within and between factors (.3, .6, and .9).}
#'   \item{case_11b}{(matrix) - Six factors with 3 indicators per factor, with different pattern coefficients within and between factors (.3, .6, and .9).}
#'   \item{case_11c}{(matrix) - Six factors with 4 indicators per factor, with different pattern coefficients within and between factors (.3, .6, and .9).}
#'   \item{case_11d}{(matrix) - Six factors with 5 indicators per factor, with different pattern coefficients within and between factors (.3, .6, and .9).}
#'   \item{case_11e}{(matrix) - Six factors with 6 indicators per factor, with different pattern coefficients within and between factors (.3, .6, and .9).}
#'   \item{case_12a}{(matrix) - One factor, with 2 equal pattern coefficients (.6).}
#'   \item{case_12b}{(matrix) - One factor, with 3 equal pattern coefficients (.6).}
#'   \item{case_12c}{(matrix) - One factor, with 6 equal pattern coefficients (.6).}
#'   \item{case_12d}{(matrix) - One factor, with 10 equal pattern coefficients (.6).}
#'   \item{case_12e}{(matrix) - One factor, with 15 equal pattern coefficients (.6).}
#'   \item{case_13a}{(matrix) - One factor, with 2 different pattern coefficients (.3, and .6).}
#'   \item{case_13b}{(matrix) - One factor, with 3 different pattern coefficients (.3, .6, and .9).}
#'   \item{case_13c}{(matrix) - One factor, with 6 different pattern coefficients (.3, .6, and .9).}
#'   \item{case_13d}{(matrix) - One factor, with 10 different pattern coefficients (.3, .6, and .9).}
#'   \item{case_13e}{(matrix) - One factor, with 15 different pattern coefficients (.3, .6, and .9).}
#'   \item{case_14a}{(matrix) - No factor, 2 variables (0).}
#'   \item{case_14b}{(matrix) - No factor, 3 variables (0).}
#'   \item{case_14c}{(matrix) - No factor, 6 variables (0).}
#'   \item{case_14d}{(matrix) - No factor, 10 variables (0).}
#'   \item{case_14e}{(matrix) - No factor, 15 variables (0).}
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
#' @source Grieder, S., & Steiner, M.D. (2020). Algorithmic Jingle Jungle:
#' A Comparison of Implementations of Principal Axis Factoring and Promax Rotation
#'  in R and SPSS. Manuscript in Preparation.
#' @source de Winter, J.C.F., & Dodou, D. (2012). Factor recovery by principal axis factoring and maximum likelihood factor analysis as a function of factor pattern and sample size. Journal of Applied Statistics. 39.
"population_models"
