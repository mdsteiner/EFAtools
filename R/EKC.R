#' Empirical Kaiser Criterion
#'
#' The empirical Kaiser criterion incorporates random sampling variations of the
#' eigenvalues from the Kaiser-Guttman criterion (\code{\link{KGC}}; see Auerswald & Moshagen
#' , 2019; Braeken & van Assen, 2017). The code is based on Braeken & van Assen, (2017) and on Auerswald and Moshagen
#' (2019).
#'
#' @param x data.frame or matrix. data.frame or matrix of raw data or matrix with
#' correlations.
#' @param N numeric. The number of observations. Only needed if x is a correlation
#'  matrix.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw
#'  data is given as input. Default is  \code{"pairwise.complete.obs"}.
#' @param cor_method character. Passed to \code{\link[stats:cor]{stats::cor}}.
#'  Default is  \code{"pearson"}.
#' @param type character. The calculation of EKC. type \code{"BvA2017"} is the original implementation; type \code{"AM2019"} differs from the original implementation but was used in simulation studies (Auerswald & Moshagen, 2019; Caron, 2025). See details.
#'  Use \code{type = c("BvA2017", "AM2019")} for both implementations. Make sure
#'  to report which version you used.
#'
#' @details The Kaiser-Guttman criterion was defined with the intend that a factor
#'  should only be extracted if it explains at least as much variance as a single
#'  factor (see \code{\link{KGC}}). However, this only applies to population-level
#'  correlation matrices. Due to sampling variation, the KGC strongly overestimates
#'  the number of factors to retrieve (e.g., Zwick & Velicer, 1986). To account
#'  for this and to introduce a factor retention method that performs well with
#'  small number of indicators and correlated factors (cases where the performance
#'  of parallel analysis, see \code{\link{PARALLEL}}, is known to deteriorate)
#'  Braeken and van Assen (2017) introduced the empirical Kaiser criterion in
#'  which a series of reference eigenvalues is created as a function of the
#'  variables-to-sample-size ratio and the observed eigenvalues.
#'
#'  Braeken and van Assen (2017) showed that "(a) EKC performs about as well as
#'  parallel analysis for data arising from the null, 1-factor, or orthogonal
#'  factors model; and (b) clearly outperforms parallel analysis for the specific
#'  case of oblique factors, particularly whenever factor intercorrelation is
#'  moderate to high and the number of variables per factor is small, which is
#'  characteristic of many applications these days" (p.463-464).
#'
#'  In EFAtools version <= 0.5.0 only the implementation of Auerswald and
#'  Moshagen (2019) was implemented (now available with
#'  \code{type = "AM2019"}). However, this implementation, that was probably also used in Caron (2025), differs from the
#'  original implementation by Braeken and van Assen (2017) in that it corrects by the reference values, i.e., without
#'  using the empirical eigenvalues used in the original implementation.
#'  Thanks to Luis Eduardo Garrido for pointing this out and to Johan Braeken for sharing
#'  sample code, based on which the original version is now implemented and used
#'  by default with \code{type = "BvA2017"}.
#'
#'  While the adapted version performed relatively well
#'  in the simulation studies by Auerswald and Moshagen (2019) and Caron (2025),
#'  the theoretical derivations of the EKC as introduced by Braeken and van Assen (2017)
#'  may no longer hold. Currently we are unaware of studies comparing the two implementations,
#'  but based on our own brief comparisons across multiple datasets, the two implementations
#'  appear to often differ substantially regarding the number of factors suggested.
#'
#'  As both implementations exist in different packages and studies, we provide both
#'  versions here. Be sure to state clearly which version you use when reporting your
#'  results to avoid confusion and ensure reproducibility.
#'
#'  The \code{EKC} function can also be called together with other factor
#'   retention criteria in the \code{\link{N_FACTORS}} function.
#'
#' @return A list of class EKC containing
#'
#' \item{eigenvalues}{A vector containing the eigenvalues found on the correlation matrix of the entered data.}
#' \item{n_factors_BvA2017}{The number of factors to retain according to the original empirical Kaiser criterion by Braeken and van Assen (2017).}
#' \item{n_factors_AM2019}{The number of factors to retain according to the adapted empirical Kaiser criterion by Auerswald and Moshagen (2019).}
#' \item{references}{The reference eigenvalues.}
#' \item{settings}{A list with the settings used.}
#'
#' @source Auerswald, M., & Moshagen, M. (2019). How to determine the number of
#' factors to retain in exploratory factor analysis: A comparison of extraction
#' methods under realistic conditions. Psychological Methods, 24(4), 468–491.
#' https://doi.org/10.1037/met0000200
#'
#' @source Braeken, J., & van Assen, M. A. (2017). An empirical Kaiser criterion.
#' Psychological Methods, 22, 450 – 466. http://dx.doi.org/10.1037/met0000074
#'
#' @source Caron, P.-O. (2025). A Comparison of the Next Eigenvalue Sufficiency Test to Other Stopping Rules for the Number of Factors in Factor Analysis.
#' Educational and Psychological Measurement, Online-first publication. https://doi.org/10.1177/00131644241308528
#'
#' @source Zwick, W. R., & Velicer, W. F. (1986). Comparison of five rules for
#' determining the number of components to retain. Psychological Bulletin, 99,
#' 432–442. http://dx.doi.org/10.1037/0033-2909.99.3.432
#'
#' @seealso Other factor retention criteria: \code{\link{CD}},
#'  \code{\link{HULL}}, \code{\link{KGC}}, \code{\link{PARALLEL}},
#'  \code{\link{SMT}}
#'
#'   \code{\link{N_FACTORS}} as a wrapper function for this and all
#'   the above-mentioned factor retention criteria.
#' @export
#'
#' @examples
#' # original implementation
#' EKC(test_models$baseline$cormat, N = 500)
#'
#' # original and adapted implementation
#' EKC(test_models$baseline$cormat, N = 500, type = c("BvA2017", "AM2019"))
EKC <- function(x, N = NA,
                use = c("pairwise.complete.obs", "all.obs",
                           "complete.obs", "everything",
                           "na.or.complete"),
                cor_method = c("pearson", "spearman", "kendall"),
                type = "BvA2017") {

  # Perform argument checks
  if(!inherits(x, c("matrix", "data.frame"))){

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n"))

  }

  checkmate::assert_count(N, na.ok = TRUE)
  use <- match.arg(use)
  cor_method <- match.arg(cor_method)
  type <- match.arg(type, choices = c("BvA2017", "AM2019"), several.ok = TRUE)

  # Check if it is a correlation matrix
  if(.is_cormat(x)){

    R <- x

    if (is.na(N)) {

      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Argument 'N' was NA but correlation matrix was entered. Please either provide N or raw data.\n"))

    }

  } else {

    message(cli::col_cyan(cli::symbol$info, " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n"))

    if (!is.na(N)) {
      warning(crayon::yellow$bold("!"), crayon::yellow(" 'N' was set and data entered. Taking N from data.\n"))
    }

    R <- stats::cor(x, use = use, method = cor_method)
    colnames(R) <- colnames(x)
    N <- nrow(x)

  }

  # Check if correlation matrix is invertable, if it is not, stop with message
  R_i <- try(solve(R), silent = TRUE)

  if (inherits(R_i, "try-error")) {
    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Correlation matrix is singular, no further analyses are performed\n"))
  }

  # Check if correlation matrix is positive definite, if it is not,
  # smooth the matrix (cor.smooth throws a warning)
  if (any(eigen(R, symmetric = TRUE, only.values = TRUE)$values <= .Machine$double.eps^.6)) {

    R <- psych::cor.smooth(R)

  }

  message(cli::col_cyan(cli::symbol$info, " The default implementation of EKC has changed compared to EFAtools version <= 0.5.0 to reflect the original version by Braeken and van Assen (2017). The previous version (which often yields different results from the original) is available with type = 'AM2019'. See details in the help page.\n"))

  n_factors_BvA2017 <- NA
  refs_BvA2017 <- NA
  n_factors_AM2019 <- NA
  refs_AM2019 <- NA


  if ("BvA2017" %in% type) {



    ### implementation in Braeken & van Assen, 2017. An Empirical Kaiser
    ### Criterion. Psychological Methods, 22(3). pp. 450-466
    ### Calculation based on p. 454 and adapted code by Johan Braeken

    # n variables
    J <- ncol(R)

    # eigenvalues
    lambda <- eigen(R, symmetric = TRUE, only.values = TRUE)$values

    # lup: asymptortic max sample eigen value under null model
    #      used as first reference eigenvalue ()
    lup <- (1 + sqrt(J/N))^2

    # correction factor
    correction_factor <- c(J,(J-cumsum(lambda))[-J])/(J:1)

    # Unrestricted EKC reference values
    l_REF <- lup* correction_factor

    # Restricted EKC reference values
    l_EKC <- l_REF
    l_EKC[l_REF < 1] <- 1

    # number of factors to retain
    n_factors_BvA2017 <- sum(lambda > l_EKC)
    refs_BvA2017 <- l_EKC

  }

  if ("AM2019" %in% type) {

    # implementation based on Auerswald and Moshagen 2019:
    # https://osf.io/fnc86?view_only=d03efba1fd0f4c849a87db82e6705668


    p <- ncol(R)

    # eigenvalues
    lambda <- eigen(R, symmetric = TRUE, only.values = TRUE)$values

    # reference values
    refs <- vector("double", p)
    for (i in seq_len(p)) {
      refs[i] <- max( ((1 + sqrt(p / N))^2) * (p - sum(refs))/
                        (p - i + 1), 1)

    }

    n_factors_AM2019 <- which(lambda <= refs)[1] - 1
    refs_AM2019 <- refs

  }


  out <- list(
    eigenvalues = lambda,
    n_factors_BvA2017 = n_factors_BvA2017,
    n_factors_AM2019 = n_factors_AM2019,
    references = data.frame(
      BvA2017 = refs_BvA2017,
      AM2019 = refs_AM2019
    ),
    settings = list(
      use = use,
      cor_method = cor_method,
      N = N,
      type = type
    )
  )

  class(out) <- "EKC"

  return(out)

}
