# DWLS estimator: diagonally weighted least squares fit to a polychoric matrix. The
# matrix and its asymptotic covariance are covered in test-polychoric.R, the routing and
# bootstrap in test-cor-method-poly.R, and the lavaan parity in test-regression-estimators.R.
# These tests cover the fitter's output contract, the analytic gradient, and the
# fit-index handling.

x_dwls <- DOSPERT_raw[stats::complete.cases(DOSPERT_raw), ]
pc_dwls <- suppressWarnings(
  .polychoric(x_dwls, acov = "diag"),
  classes = "efa_cor_sparse_cells"
)
R_dwls <- pc_dwls$R
W_dwls <- .poly_weight_matrix(pc_dwls$acov, ncol(R_dwls))

DWLS_test   <- .estimate_model(R_dwls, method = "DWLS", n_factors = 3,
                               N = nrow(x_dwls), weights = W_dwls)
DWLS_test_1 <- .estimate_model(R_dwls, method = "DWLS", n_factors = 1,
                               N = nrow(x_dwls), weights = W_dwls)

test_that("output class and dimensions are correct", {
  expect_s3_class(DWLS_test$unrot_loadings, "LOADINGS")
  expect_output(str(DWLS_test), "List of 12")
  expect_s3_class(DWLS_test_1$unrot_loadings, "LOADINGS")
  expect_output(str(DWLS_test_1), "List of 12")
})

test_that("outputs are correct", {
  expect_equal(DWLS_test$orig_R, R_dwls)
  expect_equal(sum(DWLS_test$orig_eigen), ncol(R_dwls))
  expect_lt(sum(DWLS_test$final_eigen), ncol(R_dwls))
  expect_equal(DWLS_test$convergence, 0)

  expect_equal(DWLS_test_1$orig_R, R_dwls)
  expect_equal(sum(DWLS_test_1$orig_eigen), ncol(R_dwls))
  expect_lt(sum(DWLS_test_1$final_eigen), ncol(R_dwls))
  expect_equal(DWLS_test_1$convergence, 0)
})

test_that("the chi-square block is NA but the descriptive indices are reported", {
  # DWLS does not report the ML-discrepancy chi-square or the indices derived from it;
  # the weight-free residual indices are still computed.
  expect_true(is.na(DWLS_test$fit_indices$chi))
  expect_true(is.na(DWLS_test$fit_indices$p_chi))
  expect_true(is.na(DWLS_test$fit_indices$CFI))
  expect_true(is.na(DWLS_test$fit_indices$TLI))
  expect_true(is.na(DWLS_test$fit_indices$RMSEA))
  expect_true(is.na(DWLS_test$fit_indices$AIC))
  expect_true(is.na(DWLS_test$fit_indices$BIC))

  expect_false(is.na(DWLS_test$fit_indices$RMSR))
  expect_false(is.na(DWLS_test$fit_indices$SRMR))
  expect_false(is.na(DWLS_test$fit_indices$CAF))
  expect_type(DWLS_test$fit_indices$Fm, "double")
})

test_that("the DWLS gradient matches finite differences", {
  # the analytic gradient -2 (W o (R - LL')) L of the weighted off-diagonal objective,
  # checked against central differences on a random fixture
  set.seed(1)
  p <- 6L; m <- 2L
  A <- matrix(stats::rnorm(p * p), p, p)
  Rg <- stats::cov2cor(A %*% t(A))
  Wg <- matrix(stats::runif(p * p, 0.5, 5), p, p)
  Wg <- (Wg + t(Wg)) / 2
  diag(Wg) <- 0
  par <- stats::rnorm(p * m, 0, 0.3)

  ga <- as.vector(.grad_dwls(par, Rg, m, Wg))
  h <- 1e-6
  gn <- vapply(seq_along(par), function(k) {
    e <- numeric(length(par)); e[k] <- h
    (.dwls_residuals(par + e, Rg, m, Wg) - .dwls_residuals(par - e, Rg, m, Wg)) / (2 * h)
  }, numeric(1))

  expect_equal(ga, gn, tolerance = 1e-5)
})

test_that("the DWLS objective equals the weighted off-diagonal residual sum", {
  set.seed(2)
  p <- 5L; m <- 2L
  A <- matrix(stats::rnorm(p * p), p, p)
  Rg <- stats::cov2cor(A %*% t(A))
  Wg <- matrix(stats::runif(p * p, 0.5, 5), p, p)
  Wg <- (Wg + t(Wg)) / 2
  diag(Wg) <- 0
  par <- stats::rnorm(p * m, 0, 0.3)

  L <- matrix(par, p, m)
  E <- Rg - tcrossprod(L)
  expect_equal(.dwls_residuals(par, Rg, m, Wg),
               sum(Wg[upper.tri(Wg)] * E[upper.tri(E)]^2))
})

test_that(".DWLS requires a weight matrix", {
  expect_error(.DWLS(R_dwls, n_factors = 2, weights = NULL),
               class = "efa_dwls_no_weights")
})

test_that(".poly_weight_matrix rejects a non-positive asymptotic variance", {
  expect_error(.poly_weight_matrix(c(0.1, 0, 0.2), 3),
               class = "efa_dwls_degenerate_weight")
  expect_error(.poly_weight_matrix(c(0.1, -1, 0.2), 3),
               class = "efa_dwls_degenerate_weight")
})

test_that("DWLS leaves psi NULL so a near-1 communality is not a false boundary Heywood case", {
  # .DWLS optimises free loadings with no lower bound on the uniquenesses, so it must not
  # carry psi into the box-boundary Heywood heuristic that only applies to ML/ULS.
  raw <- .DWLS(R_dwls, n_factors = 3, weights = W_dwls)
  expect_null(raw$psi)

  # a proper solution (all communalities < 1) with one communality near 1 (psi ~ 0.004):
  # not a Heywood case for DWLS, but the box-constrained ML/ULS path treats a pinned
  # uniqueness as a boundary Heywood case.
  L <- matrix(c(0.998, 0.6, 0.5, 0.4), 4, 1)
  h2 <- as.vector(diag(tcrossprod(L)))
  orig_R <- matrix(c(1, .7, .6, .5, .7, 1, .5, .4, .6, .5, 1, .3, .5, .4, .3, 1), 4)
  R_final <- orig_R; diag(R_final) <- h2
  base <- list(L = L, h2 = h2, Fm = 0, iter = 1L, convergence = 0L,
               orig_R = orig_R, R_final = R_final)

  dwls_fin <- suppressWarnings(.finalize_fit(c(base, list(psi = NULL)), N = 500, method = "DWLS"))
  expect_length(dwls_fin$heywood, 0)

  uls_fin <- suppressWarnings(.finalize_fit(c(base, list(psi = 1 - h2)), N = 500, method = "ULS"))
  expect_true(1L %in% uls_fin$heywood)
})

test_that(".fit_dwls_cpp tolerates an indefinite (non-positive-definite) correlation matrix", {
  # bootstrap replicate matrices are fed in unsmoothed and can be indefinite; the warm-start
  # inverse must not throw (as inv_sympd would), so DWLS tolerates the same resamples ULS does.
  R_npd <- matrix(c(1, 0.6, 0.6, 0.6, 1, -0.6, 0.6, -0.6, 1), 3)
  expect_lt(min(eigen(R_npd, symmetric = TRUE, only.values = TRUE)$values), 0)
  W <- matrix(c(0, 2, 3, 2, 0, 1.5, 3, 1.5, 0), 3)
  fit <- .fit_dwls_cpp(R_npd, 1L, W)
  expect_true(all(is.finite(fit$loadings)))
})

test_that("DWLS ignores the weight matrix diagonal", {
  # objective and gradient weight only off-diagonal residuals, regardless of W's diagonal
  set.seed(3)
  p <- 5L; m <- 2L
  A <- matrix(stats::rnorm(p * p), p, p)
  Rg <- stats::cov2cor(A %*% t(A))
  W0 <- matrix(stats::runif(p * p, 0.5, 5), p, p)
  W0 <- (W0 + t(W0)) / 2
  diag(W0) <- 0
  W1 <- W0; diag(W1) <- 10
  par <- stats::rnorm(p * m, 0, 0.3)

  expect_equal(.dwls_residuals(par, Rg, m, W1), .dwls_residuals(par, Rg, m, W0))
  expect_equal(as.vector(.grad_dwls(par, Rg, m, W1)), as.vector(.grad_dwls(par, Rg, m, W0)))
})

rm(x_dwls, pc_dwls, R_dwls, W_dwls, DWLS_test, DWLS_test_1)
