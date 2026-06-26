#' Velicer's Minimum Average Partial (MAP) Criterion
#'
#' @description
#' Computes Velicer's Minimum Average Partial (MAP) criterion for determining the number of
#' factors/components to retain. The function implements the original MAP criterion
#' (Velicer, 1976), expressed via the \eqn{\mathrm{TR2}} representation, and the revised
#' \eqn{\mathrm{TR4}} variant proposed by Velicer, Eaton, and Fava (2000).
#'
#' @details
#'
#' MAP is based on the idea that systematic common variance is increasingly removed from a
#' correlation matrix \eqn{R} as principal components are partialled out. After removing the
#' first \eqn{m} components, a residual (partial) covariance matrix is obtained as
#' \deqn{C_m = R - A_m A_m',}
#' where \eqn{A_m} contains the first \eqn{m} principal component loading vectors (PCA loadings).
#' This residual matrix is then standardized to a partial correlation matrix
#' \deqn{R^*_m = D_m^{-1/2} \, C_m \, D_m^{-1/2},}
#' with \eqn{D_m = \mathrm{diag}(C_m)}. The MAP criteria summarize the off-diagonal association
#' remaining in \eqn{R^*_m}. The recommended number of factors/components is the \eqn{m} that
#' minimizes the chosen criterion.
#'
#' This function returns two MAP criteria, both summarizing the off-diagonal elements
#' \eqn{r^*_{ij}} of \eqn{R^*_m}:
#' \itemize{
#'   \item **TR2 (original MAP):** \deqn{\mathrm{MAP}_m = \frac{\sum_{i \neq j} (r^*_{ij})^2}{p(p-1)}}
#'   the mean squared off-diagonal partial correlation, corresponding to Velicer's original
#'   MAP procedure.
#'   \item **TR4 (revised MAP):** \deqn{\mathrm{MAP4}_m = \frac{\sum_{i \neq j} (r^*_{ij})^4}{p(p-1)}}
#'   the mean fourth-power off-diagonal partial correlation, which places more weight on dominant
#'   residual association structure and is less affected by small partial correlations.
#' }
#'
#' **Input handling.** `x` can be a correlation matrix or raw data. If `x` is not a
#' correlation matrix, correlations are computed using [stats::cor()] with the requested
#' missing-data handling (`use`) and association measure (`cor_method`).
#'
#' **Matrix conditioning.** The function stops if the correlation matrix is singular (non-invertible),
#' because subsequent computations rely on stable matrix operations. If the correlation matrix is not
#' positive definite (e.g., due to sampling error), it is smoothed using [psych::cor.smooth()].
#'
#' **PCA-based partialing.** The PCA loading matrix \eqn{A} is obtained from the eigen-decomposition
#' of \eqn{R} as \eqn{A = V \Lambda^{1/2}}. For each \eqn{m = 0, \dots, p-1}, the first \eqn{m} columns
#' of \eqn{A} are used to compute \eqn{C_m = R - A_m A_m'}. The residual is re-standardized to the
#' partial correlation matrix \eqn{R^*_m} using \eqn{D_m^{-1/2}} (i.e., dividing by the square roots of
#' residual variances).
#'
#' **Termination.** If any residual variance (the diagonal of \eqn{C_m}) drops to at most a small
#' positive floor (1e-5) or becomes non-finite, the loop terminates early because \eqn{R^*_m}
#' cannot be formed reliably.
#'
#' @param x A numeric `matrix` or `data.frame`. Can be either (a) a correlation matrix, or
#'   (b) raw data (rows = observations, columns = variables) from which correlations are computed.
#' @param use Character string specifying the treatment of missing values when computing correlations.
#'   Passed to [stats::cor()]. Defaults to `"pairwise.complete.obs"`.
#' @param cor_method Character string specifying the correlation coefficient to be computed if raw
#'   data are supplied. One of `"pearson"`, `"spearman"`, or `"kendall"` (passed to
#'   [stats::cor()]), or `"poly"` / `"tetra"` for polychoric / tetrachoric correlations
#'   of ordinal / binary data (a two-step estimator with no empty-cell continuity
#'   correction). Defaults to `"pearson"`.
#'
#' @return An object of class `efa_retention` (see [print.efa_retention()] for the
#'   print method) with the following main elements:
#' \itemize{
#'   \item `n_factors`: A named numeric vector (`"TR2"`, `"TR4"`) with the index
#'   \eqn{m} that minimizes the original (TR2) and revised (TR4) MAP criterion.
#'   \item `results`: A list with one record per criterion, each holding the
#'   criterion values over \eqn{m}.
#'   \item `settings`: A list containing `use` and `cor_method`.
#' }
#'
#' @references
#' Velicer, W. F. (1976). Determining the number of components from the matrix of partial correlations.
#' *Psychometrika, 41*, 321--327.
#'
#' Velicer, W. F., Eaton, C. A., & Fava, J. L. (2000). Construct explication through factor or component analysis: A review and evaluation of alternative procedures for determining the number of factors or components. In Goffin, R. D. & Helmes, E. (Eds.), *Problems and Solutions in
#' Human Assessment: Honoring Douglas N. Jackson at Seventy* (pp. 41--71). Boston: Kluwer.
#'
#'
#' @examples
#' ## Example with raw data
#' res <- MAP(GRiPS_raw)
#' res
#'
#' ## Example with a correlation matrix
#' res2 <- MAP(test_models$baseline$cormat)
#' res2
#'
#' @family factor retention criteria
#'
#' @seealso [N_FACTORS()] as a wrapper function for this and the other factor
#'   retention criteria.
#'
#' @export
MAP <- function(x,
                use = c("pairwise.complete.obs", "all.obs",
                        "complete.obs", "everything",
                        "na.or.complete"),
                cor_method = c("pearson", "spearman", "kendall", "poly", "tetra")) {

  # Perform argument checks
  .assert_cor_input(x)

  use <- match.arg(use)
  cor_method <- match.arg(cor_method)

  # Detect or compute the correlation matrix, check it, and smooth it if needed
  prep <- .prepare_cor_input(x, use = use, cor_method = cor_method,
                             N_policy = "none")
  R <- prep$R

  p <- ncol(R)

  # set up m_grid
  m_max <- p - 1
  ms <- 0:m_max
  criteria <- matrix(NA_real_, nrow = p, ncol = 2)
  colnames(criteria) <- c("TR2 (orig. MAP)", "TR4 (revised MAP)")


  # PCA to get the loadings A
  ed <- eigen(R, symmetric = TRUE)
  vals <- ed$values
  vecs <- ed$vectors
  # A = V Lambda^(1/2): scale each eigenvector column by the square root of its
  # eigenvalue (equivalent to vecs %*% sqrt(diag(vals)) but without forming the
  # p x p diagonal matrix)
  A <- sweep(vecs, 2, sqrt(vals), "*")

  map_from_partials <- function(M) {
    # Velicer's (1976) MAP and the revised criterion of Velicer, Eaton & Fava (2000)
    # summarize the off-diagonal partial correlations of M (which has a unit diagonal):
    #   TR2 = mean squared off-diagonal partial correlation
    #   TR4 = mean fourth-power off-diagonal partial correlation
    # Subtracting p removes the p diagonal ones; p(p - 1) is the number of off-diagonal cells.
    c((sum(M^2) - p) / (p * (p - 1)), # TR2 criterion (original MAP)
      (sum(M^4) - p) / (p * (p - 1))) # TR4 criterion (revised MAP)
  }

  # m=0, Rstar = R
  criteria[1, ] <- map_from_partials(R)

  # run through ms
  for (m in seq_len(m_max)) {

    Am <- A[, 1:m]

    # Partial covariance: C = R - AA'
    Cm <- R - Am %*% t(Am)

    d <- diag(Cm)
    # Guard against zero/negative residual variances (can happen at very high m or numerical issues)
    if (any(!is.finite(d)) || any(d <= 1e-5)) {
      break
    }
    # D^(-1/2) to standardize to correlation matrix
    Dm <- diag(1 / sqrt(d))

    Rstar <- Dm %*% Cm %*% Dm
    criteria[m + 1, ] <- map_from_partials(Rstar)

  }


  n_factors_TR2 <- ms[which.min(criteria[, "TR2 (orig. MAP)"])]
  n_factors_TR4 <- ms[which.min(criteria[, "TR4 (revised MAP)"])]

  # one record per MAP criterion (criterion values over the number of partialled
  # components m)
  results <- list(
    list(name = "TR2", label = "Original implementation (TR2)",
         n_factors = n_factors_TR2, plot_type = "none",
         x = ms, y = criteria[, "TR2 (orig. MAP)"]),
    list(name = "TR4", label = "Revised implementation (TR4)",
         n_factors = n_factors_TR4, plot_type = "none",
         x = ms, y = criteria[, "TR4 (revised MAP)"])
  )

  out <- .new_efa_retention(
    "MAP",
    results = results,
    settings = list(use = use, cor_method = cor_method)
  )

  return(out)

}
