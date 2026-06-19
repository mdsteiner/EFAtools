#' Kaiser-Guttman Criterion
#'
#' Probably the most popular factor retention criterion. Kaiser and Guttman suggested
#' to retain as many factors as there are sample eigenvalues greater than or equal to 1.
#' This is why the criterion is also known as eigenvalues-greater-than-one rule.
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
#'  diagonal.
#' @param use character. Passed to [stats::cor()] if raw
#'  data is given as input. Default is "pairwise.complete.obs".
#' @param cor_method character. Correlation computed from raw data: `"pearson"`,
#'   `"spearman"`, or `"kendall"` (passed to [stats::cor()]), or `"poly"` /
#'   `"tetra"` for polychoric / tetrachoric correlations of ordinal / binary data
#'   (a two-step estimator with no empty-cell continuity correction).
#' Default is "pearson".
#' @param n_factors numeric. Number of factors to extract if "EFA" is included in
#' `eigen_type`. Default is 1.
#' @param ... Additional arguments passed to [EFA()]. For example,
#' to change the extraction method (PAF is default).
#'
#' @details Originally, the Kaiser-Guttman criterion was intended for the use
#' with prinicpal components, hence with eigenvalues derived from the original
#' correlation matrix. This can be done here by setting `eigen_type` to
#' "PCA". However, it is well-known that this criterion is often inaccurate and
#' that it tends to overestimate the number of factors, especially for unidimensional
#' or orthogonal factor structures (e.g., Zwick & Velicer, 1986).
#'
#' The criterion's inaccuracy in these cases is somewhat addressed if it is
#' applied on the correlation matrix with communalities in the diagonal, either
#' initial communalities estimated from SMCs (done setting `eigen_type` to
#' "SMC") or final communality estimates from an EFA (done setting `eigen_type`
#' to "EFA"; see Auerswald & Moshagen, 2019). However, although this variant
#' of the KGC is more accurate in some cases compared to the traditional KGC, it
#' is at the same time less accurate than the PCA-variant in other cases, and it
#' is still often less accurate than other factor retention methods, for
#' example parallel analysis ([PARALLEL()]), the Hull method
#' [HULL()], or sequential \eqn{chi^2} model tests ([SMT()];
#' see Auerswald & Moshagen, 2019).
#'
#' The `KGC` function can also be called together with other factor
#' retention criteria in the [N_FACTORS()] function.
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
#' Psychometrika, 19, 149 –161. http://dx.doi.org/10.1007/BF02289162
#'
#' @source Kaiser, H. F. (1960). The application of electronic computers to factor
#' analysis. Educational and Psychological Measurement, 20, 141–151.
#' http://dx.doi.org/10.1177/001316446002000116
#'
#' @source Zwick, W. R., & Velicer, W. F. (1986). Comparison of five rules for
#' determining the number of components to retain. Psychological Bulletin, 99,
#' 432–442. http://dx.doi.org/10.1037/0033-2909.99.3.432
#'
#' @seealso Other factor retention criteria: [CD()], [EKC()],
#' [HULL()], [PARALLEL()], [SMT()]
#'
#' [N_FACTORS()] as a wrapper function for this and all the
#' above-mentioned factor retention criteria.
#'
#' @export
#'
#' @examples
#' KGC(test_models$baseline$cormat, eigen_type = c("PCA", "SMC"))
KGC <- function(x, eigen_type = c("PCA", "SMC", "EFA"),
                use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                        "everything", "na.or.complete"),
                cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"), n_factors = 1,
                ...){

  # Perform argument checks
  .assert_cor_input(x)

  eigen_type <- match.arg(eigen_type, several.ok = TRUE)
  use <- match.arg(use)
  cor_method <- match.arg(cor_method)
  checkmate::assert_count(n_factors)

  # Detect or compute the correlation matrix, check it, and smooth it if needed
  prep <- .prepare_cor_input(x, use = use, cor_method = cor_method,
                             N_policy = "none")
  R <- prep$R


  # Calculate the PCA / SMC / EFA eigenvalues for the requested types
  eigen_list <- .three_eigen(R, eigen_type, n_factors = n_factors, ...)

  # Kaiser-Guttman: retain as many factors as there are eigenvalues greater than
  # or equal to 1 (NA for any eigenvalue type that was not requested)
  nfac_list <- lapply(eigen_list, function(eig) {
    if (length(eig) == 1L && is.na(eig)) NA else sum(eig >= 1)
  })

  # one record per requested eigenvalue type (eigenvalues with the >= 1 rule)
  results <- list()
  for (et in c("PCA", "SMC", "EFA")) {
    if (!(et %in% eigen_type)) next
    eig <- eigen_list[[et]]
    n_fac <- nfac_list[[et]]
    results[[et]] <- list(
      name = et,
      label = paste0(et, " eigenvalues"),
      n_factors = n_fac,
      plot_type = "eigen",
      x = seq_along(eig),
      y = eig,
      reference = NULL,
      threshold = 1,
      highlight = if (!is.na(n_fac) && n_fac >= 1) n_fac else NULL
    )
  }

  output <- .new_efa_retention(
    "KGC",
    results = unname(results),
    settings = list(eigen_type = eigen_type, use = use,
                    cor_method = cor_method, n_factors = n_factors)
  )

  return(output)

}
