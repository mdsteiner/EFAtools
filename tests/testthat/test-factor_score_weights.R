# Fixtures: an oblique and an orthogonal EFA fit supply Lambda / Phi / R that are
# passed identically to the native engine and to psych::factor.scores(), so the
# two weight matrices must agree regardless of the fit quality.
R_mat <- test_models$baseline$cormat

# The alias check below compares two separate invocations, so it uses expect_equal() with an
# explicit tolerance rather than expect_identical(): a threaded BLAS (Apple's Accelerate, for
# one) is free to vary its GEMM reduction order between calls, so one function on one input
# can differ in the last ulp. waldo still compares S3 classes, names, attributes and structure
# exactly, so only numeric values are given slack -- an alias resolved to a different weight
# formula moves the weights in the second or third decimal and still fails loudly. (A set
# tolerance does relax integer against double, which these weights never are.)
fp_tol <- 1e-8

efa_ob <- suppressMessages(suppressWarnings(
  EFA(R_mat, n_factors = 3, N = 500, type = "EFAtools", method = "PAF",
      rotation = "oblimin")))
Lambda_ob <- unclass(efa_ob$rot_loadings)
Phi_ob <- unclass(efa_ob$Phi)
h2_ob <- diag(Lambda_ob %*% Phi_ob %*% t(Lambda_ob))

efa_or <- suppressMessages(suppressWarnings(
  EFA(R_mat, n_factors = 3, N = 500, type = "EFAtools", method = "PAF",
      rotation = "varimax")))
Lambda_or <- unclass(efa_or$rot_loadings)
Phi_or <- diag(ncol(Lambda_or))
h2_or <- diag(Lambda_or %*% Phi_or %*% t(Lambda_or))


test_that("native weights match psych::factor.scores for every method", {
  skip_on_cran()
  skip_if_not_installed("psych")

  # Oblique fixture: every method except Anderson (which is orthogonal only).
  for (m in c("Thurstone", "Bartlett", "tenBerge", "Harman", "components")) {
    W <- EFAtools:::.factor_score_weights(Lambda_ob, Phi_ob, R_mat, h2_ob,
                                          method = m)
    W_psych <- psych::factor.scores(x = R_mat, f = Lambda_ob, Phi = Phi_ob,
                                    method = m)$weights
    expect_equal(W, W_psych, tolerance = 1e-6, ignore_attr = TRUE,
                 info = paste("method:", m))
  }

  # Anderson-Rubin is defined for orthogonal factors only.
  W_and <- EFAtools:::.factor_score_weights(Lambda_or, Phi_or, R_mat, h2_or,
                                            method = "Anderson")
  W_and_psych <- psych::factor.scores(x = R_mat, f = Lambda_or, Phi = Phi_or,
                                      method = "Anderson")$weights
  expect_equal(W_and, W_and_psych, tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("regression is an alias of Thurstone", {
  W_thu <- EFAtools:::.factor_score_weights(Lambda_ob, Phi_ob, R_mat, h2_ob,
                                            method = "Thurstone")
  W_reg <- EFAtools:::.factor_score_weights(Lambda_ob, Phi_ob, R_mat, h2_ob,
                                            method = "regression")
  expect_equal(W_thu, W_reg, tolerance = fp_tol)
})

test_that("orthogonal regression weights reproduce the SMC validity coefficients", {
  # For Phi = I the regression weights are W = R^-1 Lambda, so diag(Lambda' W)
  # equals diag(Lambda' R^-1 Lambda): each factor's squared multiple correlation
  # with the observed variables (Grice, 2001).
  W <- EFAtools:::.factor_score_weights(Lambda_or, Phi_or, R_mat, h2_or,
                                        method = "regression")
  smc <- diag(t(Lambda_or) %*% solve(R_mat, Lambda_or))
  expect_equal(diag(crossprod(Lambda_or, W)), smc, tolerance = 1e-8)
})

test_that("weights satisfy their estimator-defining properties (psych-free)", {
  # Closed-form identities that pin the weight *values* without psych, so the
  # numeric correctness is still guarded on CRAN (where the psych cross-check is
  # skipped).

  # Bartlett is conditionally unbiased: Lambda' W = I (Bartlett, 1937).
  W_bart <- EFAtools:::.factor_score_weights(Lambda_ob, Phi_ob, R_mat, h2_ob,
                                             method = "Bartlett")
  expect_equal(crossprod(Lambda_ob, W_bart), diag(ncol(Lambda_ob)),
               tolerance = 1e-8, ignore_attr = TRUE)

  # Anderson-Rubin scores are uncorrelated with unit variance: W' R W = I
  # (Anderson & Rubin, 1956).
  W_and <- EFAtools:::.factor_score_weights(Lambda_or, Phi_or, R_mat, h2_or,
                                            method = "Anderson")
  expect_equal(crossprod(W_and, R_mat %*% W_and), diag(ncol(Lambda_or)),
               tolerance = 1e-8, ignore_attr = TRUE)

  # ten Berge scores preserve the factor correlations: W' R W = Phi
  # (ten Berge et al., 1999).
  W_tb <- EFAtools:::.factor_score_weights(Lambda_ob, Phi_ob, R_mat, h2_ob,
                                           method = "tenBerge")
  expect_equal(crossprod(W_tb, R_mat %*% W_tb), Phi_ob,
               tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("Bartlett and Anderson weights lie in the column space of Psi^-1 Lambda", {
  # The identities above pin the right-hand normalizer but not the uniqueness
  # weighting: Lambda' W = I holds for W = D Lambda (Lambda' D Lambda)^-1 with
  # *any* symmetric D, and W' R W = I holds for the analogous Anderson-Rubin form,
  # so both stay at machine precision even if Psi^-1 is replaced by Psi, by
  # diag(1 / h2), or by the identity. What identifies the weighting is the column
  # space: both estimators are Psi^-1 Lambda times an invertible m x m factor, so
  # W must leave no residual when regressed on Psi^-1 Lambda. The residuals are
  # judged against the scale of W itself, so the bound does not silently loosen if
  # the fixture is exchanged for one with larger weights.

  # Bartlett: W = Psi^-1 Lambda (Lambda' Psi^-1 Lambda)^-1 (Bartlett, 1937).
  Lu_ob <- Lambda_ob / (1 - h2_ob)
  W_bart <- EFAtools:::.factor_score_weights(Lambda_ob, Phi_ob, R_mat, h2_ob,
                                             method = "Bartlett")
  expect_lt(max(abs(qr.resid(qr(Lu_ob), W_bart))), 1e-8 * max(abs(W_bart)))

  # Anderson-Rubin: W = Psi^-1 Lambda (Lambda' Psi^-1 R Psi^-1 Lambda)^-1/2
  # (Anderson & Rubin, 1956).
  Lu_or <- Lambda_or / (1 - h2_or)
  W_and <- EFAtools:::.factor_score_weights(Lambda_or, Phi_or, R_mat, h2_or,
                                            method = "Anderson")
  expect_lt(max(abs(qr.resid(qr(Lu_or), W_and))), 1e-8 * max(abs(W_and)))

  # W' R W = I fixes the Anderson-Rubin right factor only up to an orthogonal
  # rotation. The estimator takes the *symmetric* inverse root, which leaves
  # Lu' R W = (Lu' R Lu)^1/2 symmetric; a mis-specified Psi makes the right factor
  # non-symmetric and shows up here even though it cannot disturb W' R W = I.
  M_and <- crossprod(Lu_or, R_mat %*% W_and)
  expect_lt(max(abs(M_and - t(M_and))), 1e-8 * max(abs(M_and)))
})

test_that("the singular-matrix pseudo-inverse fallback works", {
  # .pinv is the Moore-Penrose inverse of a rank-deficient matrix: M Mp M = M
  # and Mp M Mp = Mp.
  set.seed(7)
  B <- matrix(rnorm(12), 4, 3)
  M <- tcrossprod(cbind(B, B[, 1]))  # 4x4, rank 3 (singular)
  Mp <- EFAtools:::.pinv(M)
  expect_equal(M %*% Mp %*% M, M, tolerance = 1e-8)
  expect_equal(Mp %*% M %*% Mp, Mp, tolerance = 1e-8)

  # An exactly singular correlation matrix routes solve() through .pinv instead
  # of erroring, so the engine still returns the pseudo-inverse weights.
  R_sing <- matrix(c(1, 0.5, 0.5,
                     0.5, 1,   1,
                     0.5, 1,   1), nrow = 3)  # rows 2 and 3 identical -> singular
  Lam <- matrix(c(0.7, 0.2,
                  0.5, 0.3,
                  0.4, 0.4), nrow = 3, byrow = TRUE)
  h2s <- diag(Lam %*% t(Lam))
  W_sing <- EFAtools:::.factor_score_weights(Lam, diag(2), R_sing, h2s,
                                             method = "Thurstone")
  expect_true(all(is.finite(W_sing)))
  expect_equal(W_sing, EFAtools:::.pinv(R_sing) %*% Lam,
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("non-positive uniquenesses abort for Bartlett and Anderson", {
  Lambda <- matrix(c(0.9, 0.4,
                     0.5, 0.3,
                     0.4, 0.8), nrow = 3, byrow = TRUE)
  Phi <- diag(2)
  R <- matrix(c(1,    0.30, 0.20,
                0.30, 1,    0.25,
                0.20, 0.25, 1),    nrow = 3)
  h2 <- c(1.02, 0.4, 0.5)  # first communality > 1 -> non-positive uniqueness

  expect_error(
    EFAtools:::.factor_score_weights(Lambda, Phi, R, h2, method = "Bartlett"),
    class = "efa_scores_nonpositive_uniqueness")
  expect_error(
    EFAtools:::.factor_score_weights(Lambda, Phi, R, h2, method = "Anderson"),
    class = "efa_scores_nonpositive_uniqueness")

  # A proper solution (all uniquenesses positive) does not warn.
  Lambda_ok <- matrix(c(0.7, 0.1,
                        0.6, 0.2,
                        0.2, 0.7), nrow = 3, byrow = TRUE)
  h2_ok <- diag(Lambda_ok %*% Phi %*% t(Lambda_ok))
  expect_no_warning(
    EFAtools:::.factor_score_weights(Lambda_ok, Phi, R, h2_ok, method = "Bartlett"))
})

test_that("matrix (inverse) square-root helpers are correct, including 1x1", {
  set.seed(1)
  A <- crossprod(matrix(rnorm(4), 2, 2)) + diag(2)  # 2x2 symmetric positive definite

  expect_equal(EFAtools:::.mat_sqrt(A) %*% EFAtools:::.mat_sqrt(A), A,
               tolerance = 1e-10)
  expect_equal(EFAtools:::.inv_mat_sqrt(A), solve(EFAtools:::.mat_sqrt(A)),
               tolerance = 1e-10)
  expect_equal(EFAtools:::.mat_sqrt(A) %*% EFAtools:::.inv_mat_sqrt(A), diag(2),
               tolerance = 1e-10)

  # 1x1 inputs must return the scalar (inverse) root, not a mis-dimensioned
  # identity from diag(scalar).
  expect_equal(EFAtools:::.mat_sqrt(matrix(4, 1, 1)), matrix(2, 1, 1),
               tolerance = 1e-12)
  expect_equal(EFAtools:::.inv_mat_sqrt(matrix(4, 1, 1)), matrix(0.5, 1, 1),
               tolerance = 1e-12)

  # Indefinite / singular inputs exercise the eigenvalue clamps: without them
  # .mat_sqrt would take the root of a negative eigenvalue (NaN) and
  # .inv_mat_sqrt would divide by a zero eigenvalue (Inf).
  expect_true(all(is.finite(EFAtools:::.mat_sqrt(diag(c(4, 1, -0.25))))))
  expect_true(all(is.finite(EFAtools:::.inv_mat_sqrt(diag(c(4, 1, 0))))))
})

test_that("W has the right shape and carries the loading dimnames for every method", {
  # The psych cross-check uses ignore_attr = TRUE, so dimnames are pinned here
  # (across methods, since the engine sets them in one shared step).
  for (m in c("Thurstone", "Bartlett", "tenBerge", "Harman", "components")) {
    W <- EFAtools:::.factor_score_weights(Lambda_ob, Phi_ob, R_mat, h2_ob,
                                          method = m)
    expect_equal(dim(W), dim(Lambda_ob), info = m)
    expect_identical(dimnames(W), dimnames(Lambda_ob), info = m)
  }
})

rm(R_mat, efa_ob, Lambda_ob, Phi_ob, h2_ob, efa_or, Lambda_or, Phi_or, h2_or)
