#' Comparison Data
#'
#' Factor retention method introduced by Ruscio and Roche (2012). The code was
#' adapted from the CD code by Auerswald and Moshagen (2017) available at
#' \url{https://osf.io/x5cz2/?view_only=d03efba1fd0f4c849a87db82e6705668}
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data.
#' @param n_factors_max numeric. The maximum number of factors to test against.
#'  Larger numbers will increase the duration the procedure takes, but test more
#'  possible solutions. If left NA (default) the maximum number of factors for
#'  which the model is still over-identified (df > 0) is used.
#' @param N_pop numeric. Size of finite populations of comparison data. Default
#'  is 10000.
#' @param N_samples numeric. Number of samples drawn from each population.
#'  Default is 500.
#' @param alpha numeric. The alpha level used to test the significance of the
#'  improvement added by an additional factor. Default is .30.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}}. Default
#'  is "pairwise.complete.obs".
#' @param cor_method character. Passed to \code{\link[stats:cor]{stats::cor}}.
#' Default is "pearson".
#' @param max_iter numeric. The maximum number of iterations to perform after
#'  which the iterative PAF procedure is halted. Default is 50.
#'
#' @details "Parallel analysis (PA) is an effective stopping rule that compares
#' the eigenvalues of randomly generated data with those for the actual data.
#' PA takes into account sampling error, and at present it is widely considered
#' the best available method. We introduce a variant of PA that goes even further
#' by reproducing the observed correlation matrix rather than generating random
#' data. Comparison data (CD) with known factorial structure are first generated
#' using 1 factor, and then the number of factors is increased until the
#' reproduction of the observed eigenvalues fails to improve significantly"
#' (Ruscio & Roche, 2012, p. 282).
#'
#' The CD implementation here is based on the code by Ruscio and Roche (2012), but
#' is slightly adapted to increase speed by performing the principal axis factoring
#' using a C++ based function.
#'
#' The \code{CD} function can also be called together with other factor retention
#' criteria in the \code{\link{N_FACTORS}} function.
#'
#' @return A list of class CD containing
#'
#' \item{n_factors}{The number of factors to retain according to comparison data results.}
#' \item{eigenvalues}{A vector containing the eigenvalues of the entered data.}
#' \item{RMSE_eigenvalues}{A matrix containing the RMSEs between the eigenvalues of the generated data and those of the entered data.}
#' \item{settings}{A list of the settings used.}
#'
#' @source Auerswald, M., & Moshagen, M. (2019). How to determine the number of
#' factors to retain in exploratory factor analysis: A comparison of extraction
#' methods under realistic conditions. Psychological Methods, 24(4), 468–491.
#' https://doi.org/10.1037/met0000200
#'
#' @source Ruscio, J., & Roche, B. (2012). Determining the number of factors to
#' retain in an exploratory factor analysis using comparison data of known
#' factorial structure. Psychological Assessment, 24, 282–292.
#' doi: 10.1037/a0025697
#'
#' @seealso Other factor retention criteria: \code{\link{EKC}},
#'  \code{\link{HULL}}, \code{\link{KGC}}, \code{\link{PARALLEL}}, \code{\link{SMT}}
#'
#'   \code{\link{N_FACTORS}} as a wrapper function for this and all
#'   the above-mentioned factor retention criteria.
#'
#' @export
#'
#'
#' @examples
#' \donttest{
#' # determine n factors of the GRiPS
#' CD(GRiPS_raw)
#'
#' # determine n factors of the DOSPERT risk subscale
#' CD(DOSPERT_raw)
#'}
CD <- function(x, n_factors_max = NA, N_pop = 10000, N_samples = 500, alpha = .30,
               use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                    "everything", "na.or.complete"),
               cor_method = c("pearson", "spearman", "kendall"),
               max_iter = 50) {

  # Perform argument checks
  if(!inherits(x, c("matrix", "data.frame"))){

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is neither a matrix nor a dataframe. Provide a dataframe or matrix with raw data.\n"))

  }

  if (.is_cormat(x)) {
    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is a correlation matrix, but CD only works with raw data.\n"))
  }

  if (inherits(x, c("tbl_df", "tbl"))) {
    x <- as.data.frame(x)
  }

  use <- match.arg(use)
  cor_method <- match.arg(cor_method)

  checkmate::assert_count(n_factors_max, na.ok = TRUE)
  checkmate::assert_count(N_pop)
  checkmate::assert_count(N_samples)
  checkmate::assert_number(alpha, lower = 0, upper = 1)
  checkmate::assert_count(max_iter)

  # Create correlation matrix
  R <- stats::cor(x, use = use, method = cor_method)
  colnames(R) <- colnames(x)
  n_cases <- nrow(x)
  k <- ncol(x)

  m_possible <- .det_max_factors(k)

  if (is.na(n_factors_max) || n_factors_max > m_possible) {

    if (!is.na(n_factors_max) & n_factors_max > m_possible) {
      warning(crayon::yellow$bold("!"), crayon::yellow(" n_factors_max was set to",
              n_factors_max, "but maximum possible",
              "factors to extract is", m_possible, ". Setting n_factors_max to",
              m_possible, ".\n"))
    }

    n_factors_max <- m_possible

  }

  eigvals_real <- eigen(R, symmetric = TRUE, only.values = TRUE)$values

  # initialize objects for iterative procedures
  RMSE_eigvals <- matrix(0, nrow = N_samples, ncol = n_factors_max)
  sig <- TRUE
  n_factors <- 1



  while (n_factors <= n_factors_max && isTRUE(sig)) {

    pop <- .gen_data(x, cor_method = cor_method, use, n_factors, N_pop,
                     max_iter = max_iter)

    for (j in 1:N_samples) {

      samp <- pop[sample(1:N_pop, size = n_cases, replace = TRUE),]
      R_samp <- stats::cor(samp, method = cor_method)
      eigvals_samp <- eigen(R_samp, symmetric = TRUE, only.values = TRUE)$values
      RMSE_eigvals[j,n_factors] <- sqrt(sum((eigvals_samp - eigvals_real) *
                                         (eigvals_samp - eigvals_real)) / k)
    }

    if (n_factors > 1) {

      sig <- (stats::wilcox.test(RMSE_eigvals[,n_factors],
                          RMSE_eigvals[,(n_factors - 1)], "less")$p.value < alpha)
    }
    if (isTRUE(sig)) {
      n_factors <- n_factors + 1
    }

  }

  n_factors <- n_factors - 1

  settings <- list(
    n_factors_max = n_factors_max,
    N_pop = N_pop,
    N_samples = N_samples,
    alpha = alpha,
    use = use,
    cor_method = cor_method,
    max_iter = max_iter
  )

  out <- list(
    n_factors = n_factors,
    eigenvalues = eigvals_real,
    RMSE_eigenvalues = RMSE_eigvals,
    settings = settings
  )

  class(out) <- "CD"

  return(out)

}


.gen_data <- function(x, cor_method, use, n_factors, N, max_trials = 5,
                      initial_multiplier = 1, max_iter = 100) {
  # Steps refer to description in the following article:
  # Ruscio, J., & Kaczetow, W. (2008). Simulating multivariate nonnormal data using an iterative algorithm.
  # Multivariate Behavioral Research, 43(3), 355-381.

  # Initialize variables and (if applicable) set random number seed (step 1) -------------------------------------

  k <- ncol(x)
  sim_dat <- matrix(0, nrow = N, ncol = k)         # Matrix to store the simulated data
  dists <- matrix(0, nrow = N, ncol = k)   # Matrix to store each variable's score distribution
  iter <- 0                                   # Iteration counter
  best_RMSE <- 1                                   # Lowest RMSE correlation
  t_no_impr <- 0                  # Trial counter

  # Generate distribution for each variable (step 2) -------------------------------------------------------------

  for (i in 1:k) {
    dists[,i] <- sort(sample(x[,i], size = N, replace = TRUE))
  }


  # Calculate and store a copy of the target correlation matrix (step 3) -----------------------------------------

  R <- stats::cor(x, method = cor_method, use = use)
  R_inter <- R

  # Generate random normal data for shared and unique components, initialize factor loadings (steps 5, 6) --------

  shared_comp <- matrix(stats::rnorm(N * n_factors, 0, 1), nrow = N,
                        ncol = n_factors)
  unique_comp <- matrix(stats::rnorm(N * k, 0, 1), nrow = N, ncol = k)
  shared_load <- matrix(0, nrow = k, ncol = n_factors)
  unique_load <- matrix(0, nrow = k, ncol = 1)

  # Begin loop that ends when specified number of iterations pass without improvement in RMSE correlation --------

  while (t_no_impr < max_trials) {
    iter <- iter + 1

    # Calculate factor loadings and apply to reproduce desired correlations (steps 7, 8) ---------------------------

    L <- suppressWarnings(.paf_iter(rep(1, k), criterion = .001, R = R_inter,
                   n_fac = n_factors, abs_eig = TRUE, crit_type = 2,
                   max_iter = max_iter)$L)

    shared_load[,1:n_factors] <- L

    # get rid of Heywood cases
    shared_load[shared_load > 1] <- 1
    shared_load[shared_load < -1] <- -1

    if (shared_load[1, 1] < 0) {
      shared_load <- shared_load * -1
    }

    for (i in seq_len(k)) {
      if (sum(shared_load[i,] * shared_load[i,]) < 1) {
        unique_load[i, 1] <-
          (1 - sum(shared_load[i,] * shared_load[i,]))
      } else {
        unique_load[i, 1] <- 0
      }
    }

    unique_load <- sqrt(unique_load)

    for (i in seq_len(k)) {
      sim_dat[, i] <- (shared_comp %*% t(shared_load))[, i] +
        unique_comp[, i] * unique_load[i, 1]
    }


      # Replace normal with nonnormal distributions (step 9) ---------------------------------------------------------

      for (i in seq_len(k)) {
        sim_dat <- sim_dat[sort.list(sim_dat[, i]),]
        sim_dat[,i] <- dists[, i]
      }

      # Calculate RMSE correlation, compare to lowest value, take appropriate action (steps 10, 11, 12) --------------

      R_rep <- stats::cor(sim_dat, method = cor_method)
      R_res <- R - R_rep
      # check whether this also works:
      # sqrt(sum(diag(t(R_res) %*% (R_res))) / prod(dim(R)))
      RMSE <- sqrt(sum(R_res[lower.tri(R_res)] * R_res[lower.tri(R_res)]) /
                     (.5 * (k * k - k)))

      if (RMSE < best_RMSE) {
        best_RMSE <- RMSE
        R_best <- R_inter
        res_best <- R_res
        R_inter <- R_inter + initial_multiplier * R_res
        t_no_impr <- 0
      } else {
        t_no_impr <- t_no_impr + 1
        current_multiplier <- initial_multiplier * .5 ^ t_no_impr
        R_inter <- R_best + current_multiplier * res_best
      }
  }

  # Construct the data set with the lowest RMSE correlation (step 13) --------------------------------------------

  L <- suppressWarnings(.paf_iter(rep(1, k), criterion = .001, R = R_best,
                n_fac = n_factors, abs_eig = TRUE, crit_type = 2,
                max_iter = max_iter)$L)
  shared_load[, seq_len(n_factors)] <- L

  shared_load[shared_load > 1] <- 1
  shared_load[shared_load < -1] <- -1
  if (shared_load[1, 1] < 0) {
    shared_load <- shared_load * -1
  }

  for (i in seq_len(k)) {
    if (sum(shared_load[i,] * shared_load[i,]) < 1) {
      unique_load[i, 1] <-
        (1 - sum(shared_load[i,] * shared_load[i,]))
    } else {
      unique_load[i, 1] <- 0
    }
  }

  unique_load <- sqrt(unique_load)

  for (i in seq_len(k)) {
    sim_dat[, i] <- (shared_comp %*% t(shared_load))[, i] +
      unique_comp[, i] * unique_load[i, 1]
  }

  sim_dat <- apply(sim_dat, 2, scale) # standardizes each variable in the matrix
  for (i in seq_len(k)) {
    sim_dat <- sim_dat[sort.list(sim_dat[, i]),]
    sim_dat[,i] <- dists[, i]
  }

  # Return the simulated data set (step 14) ----------------------------------------------------------------------

  return(sim_dat)
}
