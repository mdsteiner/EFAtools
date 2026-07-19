#' Scree plot
#'
#' The scree plot was originally introduced by Cattell (1966) to perform the
#' scree test. In a scree plot, the eigenvalues of the factors / components are
#' plotted against the index of the factors / components, ordered from 1 to N
#' factors components, hence from largest to smallest eigenvalue. According to
#' the scree test, the number of factors / components to retain is the number of
#' factors / components to the left of the "elbow" (where the curve starts to
#' level off) in the scree plot.
#'
#' @inheritParams efa_kgc
#'
#' @details As the scree test requires visual examination, the test has been
#' especially criticized for its subjectivity and with this low inter-rater
#' reliability. Moreover, a scree plot can be ambiguous if there are either no
#' clear "elbow" or multiple "elbows", making it difficult to judge just where
#' the eigenvalues do level off. Finally, the scree test has also been found to
#' be less accurate than other factor retention criteria. For all these reasons,
#' the scree test has been recommended against, at least for exclusive use as a
#' factor retention criterion (Zwick & Velicer, 1986)
#'
#' The `efa_scree` function can also be called together with other factor
#' retention criteria in the [efa_retain()] function.
#'
#' @returns An object of class `efa_retention` (see [print.efa_retention()] and
#'   [plot.efa_retention()] for the print and plot methods). The scree plot is a
#'   visual criterion, so it returns no numeric suggestion. Its main fields are:
#' \item{results}{A list with one record per requested eigenvalue type, each
#'   holding the eigenvalues used for the scree plot.}
#' \item{settings}{A list of the settings used.}
#'
#' @source Cattell, R. B. (1966). The scree test for the number of factors.
#' Multivariate Behavioral Research, 1(2), 245–276.
#' https://doi.org/10.1207/s15327906mbr0102_10
#'
#' @source Zwick, W. R., & Velicer, W. F. (1986). Comparison of five rules for
#' determining the number of components to retain. Psychological Bulletin, 99,
#' 432–442. http://dx.doi.org/10.1037/0033-2909.99.3.432
#'
#' @family factor retention criteria
#'
#' @seealso [efa_retain()] as a wrapper function for this and the other factor
#'   retention criteria.
#'
#' @export
#'
#' @examples
#' efa_scree(test_models$baseline$cormat, eigen_type = c("PCA", "SMC"))
efa_scree <- function(x, eigen_type = c("PCA", "SMC", "EFA"),
                  use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                        "everything", "na.or.complete"),
                  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
                  n_factors = 1, estimate_control = NULL, ...){

  # Perform argument checks
  .reject_flat_knobs(...names(), fn = "efa_scree")
  .assert_cor_input(x)

  eigen_type <- .match_arg_ci(eigen_type, several.ok = TRUE)
  use <- match.arg(use)
  cor_method <- match.arg(cor_method)
  checkmate::assert_count(n_factors)
  .assert_estimate_control(estimate_control)

  # Detect or compute the correlation matrix, check it, and smooth it if needed
  prep <- .prepare_cor_input(x, use = use, cor_method = cor_method,
                             N_policy = "none")
  R <- prep$R


  # Calculate the PCA / SMC / EFA eigenvalues for the requested types
  eigen_list <- .three_eigen(R, eigen_type, n_factors = n_factors,
                             estimate_control = estimate_control, ...)

  # one eigenvalue record per requested type; the scree plot is purely visual, so
  # there is no numeric suggestion (n_factors is NA)
  results <- list()
  for (et in c("PCA", "SMC", "EFA")) {
    if (!(et %in% eigen_type)) next
    eig <- eigen_list[[et]]
    results[[et]] <- list(
      name = et,
      label = paste0(et, " eigenvalues"),
      n_factors = NA_real_,
      plot_type = "eigen",
      x = seq_along(eig),
      y = eig
    )
  }

  output <- .new_efa_retention(
    "SCREE",
    results = unname(results),
    settings = list(eigen_type = eigen_type, use = use,
                    cor_method = cor_method, n_factors = n_factors),
    subtitle = paste0("Eigenvalues found using ",
                      cli::ansi_collapse(eigen_type), "."),
    note = "Scree plot is a visual criterion; inspect the plot to identify the elbow."
  )

  return(output)

}
