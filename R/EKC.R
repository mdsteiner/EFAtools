#' Empirical Kaiser Criterion
#'
#' The empirical Kaiser criterion incorporates random sampling variations of the
#' eigenvalues from the Kaiser-Guttman criterion ([KGC()]; see Auerswald & Moshagen
#' , 2019; Braeken & van Assen, 2017). The code is based on Braeken & van Assen, (2017) and on Auerswald and Moshagen
#' (2019).
#'
#' @param x data.frame or matrix. data.frame or matrix of raw data or matrix with
#' correlations.
#' @param N numeric. The number of observations. Only needed if x is a correlation
#'  matrix.
#' @param use character. Passed to [stats::cor()] if raw
#'  data is given as input. Default is  `"pairwise.complete.obs"`.
#' @param cor_method character. Passed to [stats::cor()].
#'  Default is  `"pearson"`.
#' @param type character. The calculation of EKC. type `"BvA2017"` is the original implementation; type `"AM2019"` differs from the original implementation but was used in simulation studies (Auerswald & Moshagen, 2019; Caron, 2025). See details.
#'  Use `type = c("BvA2017", "AM2019")` for both implementations. Make sure
#'  to report which version you used.
#'
#' @details The Kaiser-Guttman criterion was defined with the intend that a factor
#'  should only be extracted if it explains at least as much variance as a single
#'  factor (see [KGC()]). However, this only applies to population-level
#'  correlation matrices. Due to sampling variation, the KGC strongly overestimates
#'  the number of factors to retrieve (e.g., Zwick & Velicer, 1986). To account
#'  for this and to introduce a factor retention method that performs well with
#'  small number of indicators and correlated factors (cases where the performance
#'  of parallel analysis, see [PARALLEL()], is known to deteriorate)
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
#'  `type = "AM2019"`). However, this implementation, that was probably also used in Caron (2025), differs from the
#'  original implementation by Braeken and van Assen (2017) in that it corrects by the reference values, i.e., without
#'  using the empirical eigenvalues used in the original implementation.
#'  Thanks to Luis Eduardo Garrido for pointing this out and to Johan Braeken for sharing
#'  sample code, based on which the original version is now implemented and used
#'  by default with `type = "BvA2017"`.
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
#'  The `EKC` function can also be called together with other factor
#'   retention criteria in the [N_FACTORS()] function.
#'
#' @returns An object of class `efa_retention` (see [print.efa_retention()] and
#'   [plot.efa_retention()] for the print and plot methods). Its main fields are:
#' \item{n_factors}{A named numeric vector with the suggested number of factors
#'   for each requested implementation (`"BvA2017"` and/or `"AM2019"`).}
#' \item{results}{A list with one record per implementation, each holding the
#'   eigenvalues, the reference eigenvalues, and the retained solution used for
#'   printing and plotting.}
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
#' @seealso Other factor retention criteria: [CD()],
#'  [HULL()], [KGC()], [PARALLEL()],
#'  [SMT()]
#'
#'   [N_FACTORS()] as a wrapper function for this and all
#'   the above-mentioned factor retention criteria.
#' @export
#'
#' @examples
#' # original implementation
#' EKC(test_models$baseline$cormat, N = 500)
EKC <- function(x, N = NA,
                use = c("pairwise.complete.obs", "all.obs",
                           "complete.obs", "everything",
                           "na.or.complete"),
                cor_method = c("pearson", "spearman", "kendall"),
                type = "BvA2017") {

  # Perform argument checks
  .assert_cor_input(x)

  checkmate::assert_count(N, na.ok = TRUE)
  use <- match.arg(use)
  cor_method <- match.arg(cor_method)
  type <- match.arg(type, choices = c("BvA2017", "AM2019"), several.ok = TRUE)

  # Detect or compute the correlation matrix, check it, and smooth it if needed
  prep <- .prepare_cor_input(x, N = N, use = use, cor_method = cor_method,
                             N_policy = "required")
  R <- prep$R
  N <- prep$N

  # eigenvalues of the correlation matrix (shared by both implementations)
  lambda <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  J <- ncol(R)

  results <- list()

  if ("BvA2017" %in% type) {

    ### implementation in Braeken & van Assen, 2017. An Empirical Kaiser
    ### Criterion. Psychological Methods, 22(3). pp. 450-466
    ### Calculation based on p. 454 and adapted code by Johan Braeken

    # lup: asymptotic max sample eigenvalue under the null model, used as the
    #      first reference eigenvalue
    lup <- (1 + sqrt(J / N))^2

    # correction factor
    correction_factor <- c(J, (J - cumsum(lambda))[-J]) / (J:1)

    # Unrestricted EKC reference values
    l_REF <- lup * correction_factor

    # Restricted EKC reference values
    l_EKC <- l_REF
    l_EKC[l_REF < 1] <- 1

    # number of factors to retain
    temp <- cumsum(lambda > l_EKC) * cumprod(lambda > l_EKC)
    if (sum(temp) == 0) {
      n_factors_BvA2017 <- 0
    } else {
      n_factors_BvA2017 <- which.max(temp)
    }

    results[["BvA2017"]] <- list(
      name = "BvA2017",
      label = "Original implementation (Braeken & van Assen, 2017)",
      n_factors = n_factors_BvA2017,
      plot_type = "eigen",
      x = seq_along(lambda),
      y = lambda,
      reference = l_EKC,
      threshold = NULL,
      highlight = if (n_factors_BvA2017 >= 1) n_factors_BvA2017 else NULL
    )

  }

  if ("AM2019" %in% type) {

    # implementation based on Auerswald and Moshagen 2019:
    # https://osf.io/fnc86?view_only=d03efba1fd0f4c849a87db82e6705668

    p <- ncol(R)

    # reference values
    refs <- vector("double", p)
    for (i in seq_len(p)) {
      refs[i] <- max(((1 + sqrt(p / N))^2) * (p - sum(refs)) / (p - i + 1), 1)
    }

    n_factors_AM2019 <- which(lambda <= refs)[1] - 1

    results[["AM2019"]] <- list(
      name = "AM2019",
      label = "Adapted implementation (Auerswald & Moshagen, 2019)",
      n_factors = n_factors_AM2019,
      plot_type = "eigen",
      x = seq_along(lambda),
      y = lambda,
      reference = refs,
      threshold = NULL,
      highlight = if (!is.na(n_factors_AM2019) && n_factors_AM2019 >= 1) {
        n_factors_AM2019
      } else {
        NULL
      }
    )

  }

  out <- .new_efa_retention(
    "EKC",
    results = unname(results),
    settings = list(use = use, cor_method = cor_method, N = N, type = type),
    note = paste0("Multiple implementations of EKC exist; make sure to report ",
                  "which one you used (see the EKC help page for details).")
  )

  return(out)

}
