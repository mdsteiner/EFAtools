# Synthetic fixtures: a two-factor model with a known pattern and factor correlation.
# The model-implied correlation matrix R = Lambda Phi Lambda' + diag(1 - h2) is used so
# the diagnostics line up with the theoretical determinacy of the regression scores (the
# quantity the psych oracle and the closed forms below report). The r.scores identities
# hold for any PD R, so the same fixtures serve every test.
Lambda <- matrix(c(0.8, 0.1,
                   0.7, 0.0,
                   0.6, 0.1,
                   0.1, 0.8,
                   0.0, 0.7,
                   0.2, 0.6), nrow = 6, byrow = TRUE,
                 dimnames = list(paste0("V", 1:6), c("F1", "F2")))
Phi <- matrix(c(1, 0.3, 0.3, 1), 2, dimnames = list(c("F1", "F2"), c("F1", "F2")))
h2  <- diag(Lambda %*% Phi %*% t(Lambda))
S   <- Lambda %*% Phi
R   <- Lambda %*% Phi %*% t(Lambda) + diag(1 - h2)

# Orthogonal counterpart (Phi = I) for the Anderson-Rubin identity.
Phi_or <- diag(2)
h2_or  <- diag(Lambda %*% t(Lambda))
S_or   <- Lambda
R_or   <- Lambda %*% t(Lambda) + diag(1 - h2_or)


test_that("r.scores reproduces the target factor correlations", {
  # ten Berge scores preserve the factor correlations (W' R W = Phi), so their
  # intercorrelation matrix equals Phi (ten Berge et al., 1999).
  W_tb <- EFAtools:::.factor_score_weights(Lambda, Phi, R, h2, method = "tenBerge")
  d_tb <- EFAtools:::.factor_score_diagnostics(W_tb, R, S)
  expect_equal(d_tb$r.scores, Phi, tolerance = 1e-8, ignore_attr = TRUE)

  # Anderson-Rubin scores are uncorrelated with unit variance (W' R W = I), so their
  # intercorrelation matrix is the identity (Anderson & Rubin, 1956).
  W_an <- EFAtools:::.factor_score_weights(Lambda, Phi_or, R_or, h2_or, method = "Anderson")
  d_an <- EFAtools:::.factor_score_diagnostics(W_an, R_or, S_or)
  expect_equal(d_an$r.scores, diag(2), tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("score_cor has the right shape, carries the factor names, and holds determinacy", {
  W <- EFAtools:::.factor_score_weights(Lambda, Phi, R, h2, method = "Thurstone")
  d <- EFAtools:::.factor_score_diagnostics(W, R, S)

  expect_equal(dim(d$score_cor), c(ncol(Lambda), ncol(Lambda)))
  expect_identical(dimnames(d$score_cor), list(colnames(Lambda), colnames(Lambda)))
  expect_identical(dimnames(d$r.scores), dimnames(d$score_cor))
  # For a proper solution the reported determinacy is exactly the score_cor diagonal.
  expect_equal(d$determinacy, diag(d$score_cor), tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("determinacy follows psych's cap-at-1 / non-positive-to-NA convention", {
  # A single-factor toy where the weight is scaled so the score-factor "correlation"
  # exceeds 1: with R = 1, W'S / sqrt(W'RW) = S / |W| for W, S scalars.
  over <- EFAtools:::.factor_score_diagnostics(W = matrix(1), R = matrix(1),
                                               S = matrix(2))
  expect_equal(unname(over$determinacy), 1)            # capped from 2
  expect_equal(unname(over$guttman), 1)                # 2 * 1^2 - 1

  # A negative score-factor "correlation" is set to NA (and propagates to guttman).
  neg <- EFAtools:::.factor_score_diagnostics(W = matrix(1), R = matrix(1),
                                              S = matrix(-0.5))
  expect_true(is.na(neg$determinacy))
  expect_true(is.na(neg$guttman))

  # A non-PD R (here a 1x1 "correlation" of -1) makes diag(C) < 0, so the score sd and
  # the validity are NaN; the guard must report NA (not a leaked NaN). Warnings from the
  # degenerate input (sqrt/cov2cor) are expected and suppressed.
  nan <- suppressWarnings(
    EFAtools:::.factor_score_diagnostics(W = matrix(1), R = matrix(-1), S = matrix(1)))
  expect_true(is.na(nan$determinacy))
  expect_true(is.na(nan$guttman))
})

test_that("regression determinacy matches psych's factor-score validity", {
  skip_on_cran()
  skip_if_not_installed("psych")

  W  <- EFAtools:::.factor_score_weights(Lambda, Phi, R, h2, method = "Thurstone")
  d  <- EFAtools:::.factor_score_diagnostics(W, R, S)

  # psych's Grice = TRUE value is gricef(), diag(S' R^-1 S) -- the same regression
  # determinacy we report, so the two must agree to machine precision. psych's default
  # (Grice = FALSE) builds the augmented error block from rowSums(f^2) rather than the
  # true communality diag(Lambda Phi Lambda'), and diverges for oblique factors; it is
  # a different quantity, not a looser version of this one.
  fs <- psych::factor.scores(x = R, f = unclass(Lambda), Phi = unname(Phi),
                             method = "regression", Grice = TRUE)

  # psych reports R2 (squared validity); compare against our determinacy (the validity
  # rho). The measured agreement is 1.1e-16; 1e-10 keeps headroom for BLAS variation.
  expect_equal(unname(d$determinacy), sqrt(fs$R2), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("the Guttman index matches its closed form", {
  skip_on_cran()

  W <- EFAtools:::.factor_score_weights(Lambda, Phi, R, h2, method = "Thurstone")
  d <- EFAtools:::.factor_score_diagnostics(W, R, S)

  # Guttman's indeterminacy index is the minimum correlation between two sets of factor
  # scores equally consistent with the model: G = 2 rho^2 - 1, where rho is the multiple
  # correlation between a factor and its regression scores, rho_j^2 = s_j' R^-1 s_j taken
  # on the structure matrix S = Lambda Phi (Guttman, 1955, Psychometrika, 20(3), 173-197,
  # doi:10.1007/BF02289047). The check evaluates that definition directly, independently of
  # the weight matrix the implementation routes through. Cross-validated once against
  # fungible::fsIndeterminacy on this fixture, which agrees to 4.4e-16.
  rho <- sqrt(diag(t(S) %*% solve(R) %*% S))
  expect_equal(unname(d$guttman), unname(2 * rho^2 - 1), tolerance = 1e-6)
})

rm(Lambda, Phi, h2, S, R, Phi_or, h2_or, S_or, R_or)
