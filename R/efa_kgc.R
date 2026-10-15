#' Kaiser-Guttman criterion
#'
#' Retains as many factors as there are eigenvalues above a cutoff. For the
#' correlation matrix (`"PCA"`), the cutoff is 1: the classic
#' eigenvalues-greater-than-one rule (Kaiser, 1960). For the reduced correlation
#' matrix, with communality estimates in the diagonal (`"SMC"` and `"EFA"`), the
#' cutoff is 0 (Guttman, 1954).
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#' correlations.
#' @param eigen_type character. On what the eigenvalues should be found. Can be
#'  either "PCA", "SMC", or "EFA", or some combination of them. If using "PCA",
#'  the diagonal values of the correlation matrices are left to be 1. If using
#'  "SMC", the diagonal of the
#'  correlation matrices is replaced by the squared multiple correlations (SMCs)
#'  of the indicators. If using "EFA", eigenvalues are found on the correlation
#'  matrices with the final communalities of an exploratory factor analysis
#'  solution (default is principal axis factoring extracting 1 factor) as
#'  diagonal. Default is `c("PCA", "SMC", "EFA")`, i.e. all three; `"EFA"` is the
#'  only one that fits a model.
#' @param use character. Passed to [stats::cor()] if raw
#'  data is given as input. Default is "pairwise.complete.obs".
#' @param cor_method character. Correlation computed from raw data: `"pearson"`,
#'   `"spearman"`, or `"kendall"` (passed to [stats::cor()]), or `"poly"` /
#'   `"tetra"` for polychoric / tetrachoric correlations of ordinal / binary data
#'   (a two-step estimator).
#' Default is "pearson".
#' @param n_factors numeric. Number of factors to extract if "EFA" is included in
#' `eigen_type`. Default is 1.
#' @param estimate_control an [estimate_control()] object with the estimation settings for the
#'  [efa_fit()] fit that provides the communalities when `"EFA"` is included in `eigen_type`.
#'  `NULL` (default) uses the [efa_fit()] defaults. The fit is unrotated, so no rotation settings
#'  apply.
#' @param ... Additional arguments passed to [efa_fit()]. For example,
#' `estimator`, to change the estimator (PAF is default). The estimation tuning knobs are not
#' passed here; they live in `estimate_control`, and the standard-error arguments (`se`,
#' `b_boot`, `ci`, `seed`) are not accepted because the fit is an internal step that keeps
#' only its communalities.
#'
#' @details The two cutoffs apply to different matrices. The eigenvalues of the
#' correlation matrix sum to the number of variables, and an eigenvalue greater
#' than 1 identifies a component that explains more variance than one standardized
#' variable.
#' The eigenvalues of the reduced correlation matrix sum to the sum of the
#' communalities, so the cutoff 1 does not apply to them. With SMCs in the
#' diagonal, the number of positive eigenvalues is a lower bound to the number of
#' common factors in the population. This bound is the strongest of the three lower
#' bounds of Guttman (1954). `"EFA"` applies the same cutoff to the final communalities of an
#' EFA solution.
#'
#' Both rules are often inaccurate in samples. Because of sampling error, the reduced
#' correlation matrix of a sample usually has many small positive eigenvalues, so the
#' cutoff 0 tends to suggest too many factors. The classic rule also tends to
#' suggest too many factors (e.g., Zwick & Velicer, 1986; see Auerswald & Moshagen,
#' 2019, for a comparison with other criteria). Which criteria are informative
#' depends on the data at hand, so rather than substituting one for another, run
#' several of them together and compare their suggestions.
#'
#' The `efa_kgc` function can also be called together with other factor
#' retention criteria in the [efa_retain()] function.
#'
#' @returns An object of class `efa_retention` (see [print.efa_retention()] and
#'   [plot.efa_retention()] for the print and plot methods). Its main fields are:
#' \item{n_factors}{A named numeric vector with the suggested number of factors
#'   for each requested eigenvalue type (`"PCA"`, `"SMC"`, and/or `"EFA"`).}
#' \item{results}{A list with one record per eigenvalue type, each holding the
#'   eigenvalues and the retained solution used for printing and plotting.}
#' \item{settings}{A list of the settings used.}
#'
#' @source Auerswald, M., & Moshagen, M. (2019). How to determine the number of
#' factors to retain in exploratory factor analysis: A comparison of extraction
#' methods under realistic conditions. Psychological Methods, 24(4), 468–491.
#' https://doi.org/10.1037/met0000200
#'
#' @source Guttman, L. (1954). Some necessary conditions for common-factor analysis.
#' Psychometrika, 19, 149–161. https://doi.org/10.1007/BF02289162
#'
#' @source Kaiser, H. F. (1960). The application of electronic computers to factor
#' analysis. Educational and Psychological Measurement, 20, 141–151.
#' https://doi.org/10.1177/001316446002000116
#'
#' @source Zwick, W. R., & Velicer, W. F. (1986). Comparison of five rules for
#' determining the number of components to retain. Psychological Bulletin, 99,
#' 432–442. https://doi.org/10.1037/0033-2909.99.3.432
#'
#' @family factor retention criteria
#'
#' @seealso [efa_retain()] as a wrapper function for this and the other factor
#'   retention criteria.
#'
#' @export
#'
#' @examples
#' efa_kgc(test_models$baseline$cormat, eigen_type = c("PCA", "SMC"))
efa_kgc <- function(x, eigen_type = c("PCA", "SMC", "EFA"),
                use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                        "everything", "na.or.complete"),
                cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"), n_factors = 1,
                estimate_control = NULL,
                ...){

  # Perform argument checks
  .reject_flat_knobs(...names(), fn = "efa_kgc")
  .reject_unknown_fit_dots(...names(), fn = "efa_kgc", unrotated = TRUE)
  .reject_rotation_dots(list(...), fn = "efa_kgc")
  .assert_cor_input(x)

  eigen_type <- .match_arg_ci(eigen_type, several.ok = TRUE)
  use <- .match_arg_ci(use)
  cor_method <- .match_arg_ci(cor_method)
  checkmate::assert_count(n_factors)
  .assert_estimate_control(estimate_control)

  # Detect or compute the correlation matrix, check it, and smooth it if needed
  prep <- .prepare_cor_input(x, use = use, cor_method = cor_method,
                             N_policy = "none")
  R <- prep$R


  # Calculate the PCA / SMC / EFA eigenvalues for the requested types
  eigen_list <- .three_eigen(R, eigen_type, n_factors = n_factors,
                             estimate_control = estimate_control, ...)

  # Count the eigenvalues above the cutoff: 1 for the correlation matrix (Kaiser, 1960), 0
  # for the reduced matrices (Guttman, 1954)
  cutoff <- c(PCA = 1, SMC = 0, EFA = 0)

  # one record per requested eigenvalue type
  results <- list()
  for (et in c("PCA", "SMC", "EFA")) {
    if (!(et %in% eigen_type)) next
    eig <- eigen_list[[et]]
    # p * eps * max|eig| absorbs eigensolver round-off, not the convergence error of an "EFA" fit
    tol <- length(eig) * .Machine$double.eps * max(abs(eig))
    n_fac <- sum(eig > cutoff[[et]] + tol)
    results[[et]] <- list(
      name = et,
      label = et,
      n_factors = n_fac,
      plot_type = "eigen",
      x = seq_along(eig),
      y = eig,
      reference = NULL,
      threshold = cutoff[[et]],
      highlight = if (n_fac >= 1) n_fac else NULL
    )
  }

  # the report states the cutoff of each requested type
  rule <- if (all(eigen_type == "PCA")) {
    "cutoff 1"
  } else if (!"PCA" %in% eigen_type) {
    "cutoff 0"
  } else {
    paste0("cutoff 1 (PCA), 0 (", paste(setdiff(eigen_type, "PCA"), collapse = ", "), ")")
  }

  output <- .new_efa_retention(
    "KGC",
    results = unname(results),
    settings = list(eigen_type = eigen_type, use = use,
                    cor_method = cor_method, n_factors = n_factors),
    subtitle = .eigen_subtitle(eigen_type, rule)
  )

  return(output)

}
