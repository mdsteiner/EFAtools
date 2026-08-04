#' Velicer's minimum average partial (MAP) criterion
#'
#' @description
#' Computes Velicer's Minimum Average Partial (MAP) criterion for determining the number of
#' factors/components to retain. The function implements the original MAP criterion
#' (Velicer, 1976), expressed via the \eqn{\mathrm{TR2}} representation, and the revised
#' \eqn{\mathrm{TR4}} variant proposed by Velicer, Eaton, and Fava (2000).
#'
#' @details
#'
#' MAP partials successive principal components out of the correlation matrix and,
#' after removing \eqn{m} components, summarizes the off-diagonal partial
#' correlations \eqn{r^*_{ij}} that remain in the \eqn{m}-th partial correlation
#' matrix \eqn{M} (which has a unit diagonal); the suggested number of factors is
#' the \eqn{m} that minimizes the criterion. Two criteria are returned, each
#' rescaling the trace of a matrix power of \eqn{M} by the number of off-diagonal
#' cells \eqn{p(p-1)}:
#' \itemize{
#'   \item **TR2 (original MAP; Velicer, 1976):** the average squared off-diagonal
#'   partial correlation,
#'   \deqn{\mathrm{TR2}_m = \frac{\mathrm{tr}(M^2) - p}{p(p-1)} = \frac{\sum_{i \neq j} (r^*_{ij})^2}{p(p-1)},}
#'   where subtracting \eqn{p} removes the \eqn{p} unit diagonal entries.
#'   \item **TR4 (revised MAP; Velicer, Eaton, & Fava, 2000):** the analogous
#'   fourth-power summary, formed from the trace of the fourth matrix power,
#'   \deqn{\mathrm{TR4}_m = \frac{\mathrm{tr}(M^4) - p}{p(p-1)}.}
#'   Moving from the squared to the fourth power downweights the small partial
#'   correlations relative to the large ones, which can sharpen the minimum.
#'   Unlike TR2, \eqn{\mathrm{tr}(M^4)} is *not* the sum of the fourth powers of the
#'   individual partial correlations; the matrix power is intended and is what
#'   Velicer, Eaton, and Fava (2000) describe.
#' }
#'
#' MAP is most dependable when the components are well determined, that is with many
#' indicators per factor and substantial loadings. It has a well-documented tendency
#' to under-extract, particularly with few indicators per factor or weak loadings
#' (Zwick & Velicer, 1986; Auerswald & Moshagen, 2019), so it is best read as a lower
#' bound and paired with a criterion that errs in the other direction, such as the
#' Kaiser-Guttman criterion ([efa_kgc()]).
#'
#' A non-positive-definite input correlation matrix (e.g. from sampling error) is
#' smoothed with [psych::cor.smooth()].
#'
#' @param x A numeric `matrix` or `data.frame`. Can be either (a) a correlation matrix, or
#'   (b) raw data (rows = observations, columns = variables) from which correlations are computed.
#' @param use Character string specifying the treatment of missing values when computing correlations.
#'   Passed to [stats::cor()]. Defaults to `"pairwise.complete.obs"`.
#' @param cor_method Character string specifying the correlation coefficient to be computed if raw
#'   data are supplied. One of `"pearson"`, `"spearman"`, or `"kendall"` (passed to
#'   [stats::cor()]), or `"poly"` / `"tetra"` for polychoric / tetrachoric correlations
#'   of ordinal / binary data (a two-step estimator). Defaults to `"pearson"`.
#'
#' @returns An object of class `efa_retention` (see [print.efa_retention()] for the
#'   print method). MAP has no plot; [plot.efa_retention()] returns `NULL` with a
#'   message for it. Its main elements are:
#' \itemize{
#'   \item `n_factors`: A named numeric vector (`"TR2"`, `"TR4"`) with the index
#'   \eqn{m} that minimizes the original (TR2) and revised (TR4) MAP criterion.
#'   \item `results`: A list with one record per criterion, each holding the
#'   criterion values over \eqn{m}.
#'   \item `settings`: A list containing `use` and `cor_method`.
#' }
#'
#' @source Auerswald, M., & Moshagen, M. (2019). How to determine the number of
#' factors to retain in exploratory factor analysis: A comparison of extraction
#' methods under realistic conditions. *Psychological Methods, 24*(4), 468--491.
#' https://doi.org/10.1037/met0000200
#' @source Velicer, W. F. (1976). Determining the number of components from the matrix of partial correlations.
#' *Psychometrika, 41*, 321--327.
#' @source Velicer, W. F., Eaton, C. A., & Fava, J. L. (2000). Construct explication through factor or component analysis: A review and evaluation of alternative procedures for determining the number of factors or components. In Goffin, R. D. & Helmes, E. (Eds.), *Problems and Solutions in
#' Human Assessment: Honoring Douglas N. Jackson at Seventy* (pp. 41--71). Boston: Kluwer.
#' @source Zwick, W. R., & Velicer, W. F. (1986). Comparison of five rules for
#' determining the number of components to retain. *Psychological Bulletin, 99*,
#' 432--442. http://dx.doi.org/10.1037/0033-2909.99.3.432
#'
#'
#' @examples
#' ## Example with raw data
#' res <- efa_map(GRiPS_raw)
#' res
#'
#' ## Example with a correlation matrix
#' res2 <- efa_map(test_models$baseline$cormat)
#' res2
#'
#' @family factor retention criteria
#'
#' @seealso [efa_retain()] as a wrapper function for this and the other factor
#'   retention criteria.
#'
#' @export
efa_map <- function(x,
                use = c("pairwise.complete.obs", "all.obs",
                        "complete.obs", "everything",
                        "na.or.complete"),
                cor_method = c("pearson", "spearman", "kendall", "poly", "tetra")) {

  # Perform argument checks
  .assert_cor_input(x)

  use <- .match_arg_ci(use)
  cor_method <- .match_arg_ci(cor_method)

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
    M2 <- M %*% M
    M4 <- M2 %*% M2
    c((sum(diag(M2)) - p) / (p * (p - 1)), # TR2 criterion (original MAP)
      (sum(diag(M4)) - p) / (p * (p - 1))) # TR4 criterion (revised MAP)
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
