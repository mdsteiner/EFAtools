#' Maximum Likelihood Estimation of Factor Loadings
#'
#' @param x matrix or data.frame. A raw data or correlation matrix.
#' @param n_factors numeric. The number of factors to extract.
#' @param cors logical. Whether x is a correlation matrix.
#' @param N numeric. The number of cases. Only necessary if correlation matrix is
#'  specified. Needed for some fit indices.
#' @param signed_loadings logical. If \code{TRUE} (default), the sign of
#' factors with negative sum of loadings is reflected.
#' @param start_method character. How to specify the starting values for the
#'  optimization prodedure. Default is "factanal" which takes the starting values
#'  specified in the \link{stats}{factanal} function. "psych" takes the starting
#'  values specified in \link{psych}{fa}. Solutions are very similar.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw data
#'  is given as input. Note that in this case \code{cors} must be set to
#'  \code{FALSE}. Default is "pairwise.complete.obs".
#'
#' @return A list of class ML containing the following
#' \item{orig_R}{Original correlation matrix.}
#' \item{h2}{Final communality estimates.}
#' \item{orig_eigen}{Eigen values of the original correlation matrix.}
#' \item{final_eigen}{Eigenvalues of the final correlation matrix with the estimated communalities as diagonal.}
#' \item{unrot_loadings}{Loading matrix containing the final loadings.}
#' \item{vars_accounted}{Matrix of explained variances and sums of squared loadings}
#' \item{fit_indices}{Fit indices as returned by
#'  \code{\link[psych:factor.stats]{psych::factor.stats}}}
#' \item{settings}{list. The settings (arguments) used in ML}
#' @export
#'
#' @examples
#' # call as single function
#' ML(IDS2_R, n_factors = 5)
ML <- function(x, n_factors, cors = TRUE, N = NA, signed_loadings = TRUE,
               start_method = c("factanal", "psych"), use = "pairwise.complete.obs") {

  start_method <- match.arg(start_method)
  # create R correlation matrix object, if from data, using
  # pairwise binary correlations
  if (isTRUE(cors)) {
    R <- x

    # test whether a real correlation matrix is used
    if (nrow(R) != ncol(R)) {
      stop("Entered data is no correlation matrix but cors = TRUE. Either set ",
           "cors = FALSE if you entered raw data, or enter a correlation matrix.")
    }

    if (is.null(N)) {
      stop("Argument 'N' is NULL. Either provide N, N = NA, or raw data.")
    }

  } else {
    R <- stats::cor(x, use = use)
    colnames(R) <- colnames(x)
    N <- nrow(x)
  }

  ml <- .fit_ml(R, n_factors, start_method)

  L <- ml$loadings
  orig_R <- R
  h2 = diag(L %*% t(L))
  diag(R) <- h2

  if (signed_loadings) {
    # reverse the sign of loadings as done in the psych package,
    # and spss
    if (n_factors > 1) {
      signs <- sign(colSums(L))
      signs[signs == 0] <- 1
      L <- L %*% diag(signs)
    } else {
      if (sum(L) < 0) {
        L <- -as.matrix(L)
      } else {
        L <- as.matrix(L)
      }

    }

  }

  if (!is.null(colnames(orig_R))) {
    # name the loading matrix so the variables can be identified
    rownames(L) <- colnames(orig_R)
  }

  colnames(L) <- paste0("ML", 1:n_factors)

  vars_accounted <- .compute_vars(L_unrot = L, L_rot = L)

  colnames(vars_accounted) <- colnames(L)

  # compute fit indices
  fit_ind <- try(psych::factor.stats(orig_R, L, n.obs = N), silent = TRUE)

  if (all(class(fit_ind) == "try-error")) {
    fit_ind <- NA
  }


  # create the output object
  class(L) <- "LOADINGS"

  # store the settings used:

  settings <- list(
    N = N,
    iter = ml$res$counts[1],
    convergence = ml$res$convergence,
    signed_loadings = signed_loadings
  )


  output <- list(
    orig_R = orig_R,
    h2 = diag(L %*% t(L)),
    orig_eigen = eigen(orig_R, symmetric = TRUE)$values,
    final_eigen = eigen(R, symmetric = TRUE)$values,
    unrot_loadings = L,
    vars_accounted = vars_accounted,
    fit_indices = fit_ind,
    settings = settings
  )

  class(output$h2) <- "COMMUNALITIES"
  class(output$orig_eigen) <- "EIGEN"
  class(output$final_eigen) <- "EIGEN"

  class(output) <- "ML"

  output
}

# function to obtain the uls fit; adapted from the psych package
.fit_ml <- function(R, n_fac, start_method) {

  if (start_method == "psych") {
    R.smc <- (1 - 1 / diag(solve(R)))
    if((sum(R.smc) == n_fac) && (n_fac > 1)) {
      start <- rep(.5, n_fac)
    }  else {
      start <- diag(R)- R.smc
    }
  } else if (start_method == "factanal") {
    start <- (1 - 0.5 * n_fac / ncol(R)) / diag(solve(R))
  }


  res <- stats::optim(start, .error_ml, gr = .grad_ml, method = "L-BFGS-B",
                      lower = .005, upper = 1,
                      control = c(list(fnscale = 1,
                                       parscale = rep(0.01, length(start)))),
                      R = R, n_fac = n_fac)

  Lambda <- .FAout(res$par, R, n_fac)

  result <- list(loadings = Lambda, res = res, R = R)
}

# .error_ml2 <- function(psi, R, n_fac)
# {
#   sc <- diag(1/sqrt(psi))
#   Rs <- sc %*% R %*% sc
#   E <- eigen(Rs, symmetric = TRUE, only.values = TRUE)
#   e <- E$values[-(1:n_fac)]
#   e <- sum(log(e) - e) - n_fac + nrow(R)
#   -e
# }
# .grad_ml2 <- function(psi, R, n_fac) {
#   sc <- diag(1 / sqrt(psi))
#   Rs <- sc %*% R %*% sc
#   E <- eigen(Rs, symmetric = TRUE)
#   L <- E$vectors[, 1:n_fac, drop = FALSE]
#   load <- L %*% diag(sqrt(pmax(E$values[1:n_fac] - 1, 0)), n_fac)
#   load <- diag(sqrt(psi)) %*% load
#   g <- load %*% t(load) + diag(psi) - R     # g <- model - data
#   diag(g) / psi^2                             #normalized
# }

# taken from factanal
.FAout <- function(psi, R, n_fac) {
  sc <- diag(1 / sqrt(psi))
  Rs <- sc %*% R %*% sc
  E <- eigen(Rs, symmetric = TRUE)
  L <- E$vectors[, 1L:n_fac, drop = FALSE]
  load <- L %*% diag(sqrt(pmax(E$values[1L:n_fac] - 1, 0)),
                     n_fac)
  diag(sqrt(psi)) %*% load
}
