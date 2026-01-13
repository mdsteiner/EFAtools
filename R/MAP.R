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
#' #' MAP is based on the idea that systematic common variance is increasingly removed from a
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
#' This function returns two MAP criteria:
#' \itemize{
#'   \item \strong{TR2 (original MAP):} \deqn{\mathrm{MAP}_m = \frac{\mathrm{Trace}(R^{*2}_m) - p}{p(p-1)}}
#'   which is algebraically equivalent to the mean squared off-diagonal partial correlations and
#'   corresponds to Velicer's original MAP procedure.
#'   \item \strong{TR4 (revised MAP):} \deqn{\mathrm{MAP4}_m = \frac{\mathrm{Trace}(R^{*4}_m) - p}{p(p-1)}}
#'   a higher-order variant that places more weight on dominant residual association structure.
#' }
#'
#' \strong{Input handling.} \code{x} can be a correlation matrix or raw data. If \code{x} is not a
#' correlation matrix, correlations are computed using \code{\link[stats]{cor}} with the requested
#' missing-data handling (\code{use}) and association measure (\code{cor_method}). If a correlation
#' matrix is supplied, \code{N} must be provided.
#'
#' \strong{Matrix conditioning.} The function stops if the correlation matrix is singular (non-invertible),
#' because subsequent computations rely on stable matrix operations. If the correlation matrix is not
#' positive definite (e.g., due to sampling error), it is smoothed using \code{\link[psych]{cor.smooth}}.
#'
#' \strong{PCA-based partialing.} The PCA loading matrix \eqn{A} is obtained from the eigen-decomposition
#' of \eqn{R} as \eqn{A = V \Lambda^{1/2}}. For each \eqn{m = 0, \dots, p-1}, the first \eqn{m} columns
#' of \eqn{A} are used to compute \eqn{C_m = R - A_m A_m'}. The residual is re-standardized to the
#' partial correlation matrix \eqn{R^*_m} using \eqn{D_m^{-1/2}} (i.e., dividing by the square roots of
#' residual variances).
#'
#' \strong{Termination.} If residual variances (the diagonal of \eqn{C_m}) become non-positive or
#' numerically unstable, the loop terminates early because \eqn{R^*_m} cannot be formed reliably.
#'
#' @param x A numeric \code{matrix} or \code{data.frame}. Can be either (a) a correlation matrix, or
#'   (b) raw data (rows = observations, columns = variables) from which correlations are computed.
#' @param use Character string specifying the treatment of missing values when computing correlations.
#'   Passed to \code{\link[stats]{cor}}. Defaults to \code{"pairwise.complete.obs"}.
#' @param cor_method Character string specifying the correlation coefficient to be computed if raw
#'   data are supplied. Passed to \code{\link[stats]{cor}}. Defaults to \code{"pearson"}.
#'
#' @return An object of class \code{"MAP"} with the following elements:
#' \itemize{
#'   \item \code{eigenvalues}: Eigenvalues of the (possibly smoothed) correlation matrix.
#'   \item \code{n_factors_TR2}: Index \eqn{m} that minimizes the TR2 (original MAP) criterion.
#'   \item \code{n_factors_TR4}: Index \eqn{m} that minimizes the TR4 (revised MAP) criterion.
#'   \item \code{criteria}: A \code{matrix} with columns \code{m}, \code{TR2 (orig. MAP)}, and
#'   \code{TR4 (revised MAP)}.
#'   \item \code{settings}: A list containing \code{use}, \code{cor_method}, and \code{N}.
#' }
#'
#' @references
#' Velicer, W. F. (1976). Determining the number of components from the matrix of partial correlations.
#' \emph{Psychometrika, 41}, 321--327.
#'
#' Velicer, W. F., Eaton, C. A., \& Fava, J. L. (2000). Construct explication through factor or component analysis: A review and evaluation of alternative procedures for determining the number of factors or components. In R. D. Goffin \& E. Helmes (Eds.), \emph{Problems and Solutions in
#' Human Assessment: Honoring Douglas N. Jackson at Seventy} (pp. 41--71). Boston: Kluwer.
#'
#'
#' @examples
#' ## Example with raw data
#' res <- MAP(GRiPS_raw)
#' res
#'
#' ## Example with a correlation matrix
#' res2 <- MAP(test_models$baseline$cormat, N = 500)
#' res2
#'
#' @export
MAP <- function(x,
                use = c("pairwise.complete.obs", "all.obs",
                        "complete.obs", "everything",
                        "na.or.complete"),
                cor_method = c("pearson", "spearman", "kendall")) {

  # Perform argument checks
  if(!inherits(x, c("matrix", "data.frame"))){

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n"))

  }

  use <- match.arg(use)
  cor_method <- match.arg(cor_method)

  # Check if it is a correlation matrix
  if(.is_cormat(x)){

    R <- x

  } else {

    message(cli::col_cyan(cli::symbol$info, " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n"))

    R <- stats::cor(x, use = use, method = cor_method)
    colnames(R) <- colnames(x)

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

  p <- ncol(R)

  # set up m_grid
  m_max <- p - 1
  ms <- 0:m_max
  criteria <- matrix(NA_real_, nrow = p, ncol = 3)
  criteria[, 1] <- ms
  colnames(criteria) <- c("m", "TR2 (orig. MAP)", "TR4 (revised MAP)")


  # PCA to get the loadings A
  ed <- eigen(R, symmetric = TRUE)
  vals <- ed$values
  vecs <- ed$vectors
  A <- vecs %*% sqrt(diag(vals))

  map_from_trace_power <- function(M) {
    # "Trace of the matrix ... to the k-th power" in the VEF (2000) sense (matrix power),
    # then transformed to the MAP scale:
    #   (Trace(M^k) - p) / (p(p-1))  (matches the TR2-to-MAP equivalence described in the chapter)
    M2 <- M %*% M
    M4 <- M2 %*% M2
    c((sum(diag(M2)) - p) / (p * (p - 1)), # TR2 criterion
      (sum(diag(M4)) - p) / (p * (p - 1))) # TR4 criterion
  }

  # m=0, Rstar = R
  criteria[1, c(2, 3)] <- map_from_trace_power(R)

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
    criteria[m + 1, c(2, 3)] <- map_from_trace_power(Rstar)

  }


  out <- list(
    eigenvalues = vals,
    n_factors_TR2 = ms[which.min(criteria[, 2])],
    n_factors_TR4 = ms[which.min(criteria[, 3])],
    criteria = criteria,
    settings = list(
      use = use,
      cor_method = cor_method
    )
  )

  class(out) <- "MAP"

  return(out)

}
