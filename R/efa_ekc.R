#' Empirical Kaiser criterion
#'
#' The empirical Kaiser criterion incorporates random sampling variations of the
#' eigenvalues from the Kaiser-Guttman criterion ([efa_kgc()]; see Auerswald &
#' Moshagen, 2019; Braeken & van Assen, 2017). The implementation follows Braeken
#' and van Assen (2017).
#'
#' @inheritParams efa_kgc
#' @param N numeric. The number of observations. Only needed if x is a correlation
#'  matrix. Must be larger than the number of variables.
#' @param cor_method character. Correlation computed from raw data: `"pearson"`,
#'   `"spearman"`, or `"kendall"` (passed to [stats::cor()]), or `"poly"` /
#'   `"tetra"` for polychoric / tetrachoric correlations of ordinal / binary data
#'   (a two-step estimator). Default is `"pearson"`. Note that the EKC reference
#'   values rest on the Marchenko-Pastur law for the eigenvalues of a sample
#'   correlation matrix of independent variables, which assumes the sampling
#'   behaviour of product-moment correlations; with `"poly"` / `"tetra"` (and, to a
#'   lesser degree, the rank-based methods) the reference series is therefore an
#'   approximation.
#' @param type `r lifecycle::badge("deprecated")` Accepted and ignored. It selected
#'   between two ways to compute the reference values. The `"AM2019"` reference values
#'   do not depend on the observed eigenvalues, so they do not apply the empirical
#'   correction that defines the criterion, and they are no longer computed.
#'
#' @details The Kaiser-Guttman criterion was defined with the intend that a factor
#'  should only be extracted if it explains at least as much variance as a single
#'  factor (see [efa_kgc()]). However, this only applies to population-level
#'  correlation matrices. Due to sampling variation, the KGC strongly overestimates
#'  the number of factors to retrieve (e.g., Zwick & Velicer, 1986). To account
#'  for this and to introduce a factor retention method that performs well with
#'  small number of indicators and correlated factors (cases where the performance
#'  of parallel analysis, see [efa_parallel()], is known to deteriorate)
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
#' @returns An object of class `efa_retention` (see [print.efa_retention()] and
#'   [plot.efa_retention()] for the print and plot methods). Its main fields are:
#' \item{n_factors}{A numeric vector of length one, named `"BvA2017"`, with the
#'   suggested number of factors. The factors up to the first observed eigenvalue
#'   that fails to exceed its reference value are retained. The "all-exceed"
#'   convention of parallel analysis ([efa_parallel()]), which retains all *J* factors
#'   when no such crossing is found, cannot be reached here: the reference values are
#'   never below 1, while the eigenvalues of a correlation matrix sum to *J* and are
#'   sorted downwards, so the last of them is never above 1.}
#' \item{results}{A list with one record, holding the eigenvalues, the reference
#'   eigenvalues, and the retained solution used for printing and plotting.}
#' \item{settings}{A list with the settings used.}
#'
#' @source Auerswald, M., & Moshagen, M. (2019). How to determine the number of
#' factors to retain in exploratory factor analysis: A comparison of extraction
#' methods under realistic conditions. Psychological Methods, 24(4), 468–491.
#' https://doi.org/10.1037/met0000200
#'
#' @source Braeken, J., & van Assen, M. A. (2017). An empirical Kaiser criterion.
#' Psychological Methods, 22, 450 – 466. https://doi.org/10.1037/met0000074
#'
#' @source Zwick, W. R., & Velicer, W. F. (1986). Comparison of five rules for
#' determining the number of components to retain. Psychological Bulletin, 99,
#' 432–442. https://doi.org/10.1037/0033-2909.99.3.432
#'
#' @family factor retention criteria
#'
#' @seealso [efa_retain()] as a wrapper function for this and the other factor
#'   retention criteria.
#' @export
#'
#' @examples
#' efa_ekc(test_models$baseline$cormat, N = 500)
efa_ekc <- function(x, N = NA,
                use = c("pairwise.complete.obs", "all.obs",
                           "complete.obs", "everything",
                           "na.or.complete"),
                cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
                type = lifecycle::deprecated()) {

  # Perform argument checks
  .assert_cor_input(x)

  checkmate::assert_count(N, na.ok = TRUE)
  use <- .match_arg_ci(use)
  cor_method <- .match_arg_ci(cor_method)

  if (lifecycle::is_present(type)) .deprecate_ekc_type("efa_ekc(type)")

  # Detect or compute the correlation matrix, check it, and smooth it if needed
  prep <- .prepare_cor_input(x, N = N, use = use, cor_method = cor_method,
                             N_policy = "required")
  R <- prep$R
  N <- prep$N
  .assert_n_gt_vars(N, ncol(R))

  # eigenvalues of the correlation matrix
  lambda <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  J <- ncol(R)

  ### implementation in Braeken & van Assen, 2017. An Empirical Kaiser
  ### Criterion. Psychological Methods, 22(3). pp. 450-466
  ### Calculation based on p. 454

  # lup: asymptotic max sample eigenvalue under the null model, used as the
  #      first reference eigenvalue
  lup <- (1 + sqrt(J / N))^2

  # correction factor: the variance the preceding eigenvalues do not explain, divided
  # over the factors that remain
  correction_factor <- c(J, (J - cumsum(lambda))[-J]) / (J:1)

  # Unrestricted EKC reference values
  l_REF <- lup * correction_factor

  # Restricted EKC reference values
  l_EKC <- l_REF
  l_EKC[l_REF < 1] <- 1

  # Number of factors to retain: the length of the leading run of eigenvalues above their
  # reference value. cumprod() is 1 across that run and 0 from the first failure onwards,
  # so its sum is the length of the run.
  #
  # cumprod() coerces the logical comparison to double, so the count is cast back to the
  # integer storage mode a count belongs in.
  n_factors <- as.integer(sum(cumprod(lambda > l_EKC)))

  out <- .new_efa_retention(
    "EKC",
    results = list(list(
      name = "BvA2017",
      label = "Braeken & van Assen (2017)",
      n_factors = n_factors,
      plot_type = "eigen",
      x = seq_along(lambda),
      y = lambda,
      reference = l_EKC,
      threshold = NULL,
      highlight = if (n_factors >= 1) n_factors else NULL
    )),
    settings = list(use = use, cor_method = cor_method, N = N)
  )

  return(out)

}
