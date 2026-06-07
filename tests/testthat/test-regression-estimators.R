# Regression net: pin the current estimator output (PAF / ML / ULS) against the
# reference implementations it is designed to reproduce, so an accidental change is
# caught. Comparisons use quantities that are invariant to the sign, order, and
# orthogonal rotation of the unrotated factors (per-variable uniquenesses and the
# off-diagonal reproduced correlations L %*% t(L)); both are robust yet still pin the
# identified solution. Tolerance is set generously relative to the observed agreement
# (noted per block) so the contract is portable across BLAS implementations.

# off-diagonal of the reproduced correlation matrix from a loading matrix
repro_offdiag <- function(loadings) {
  L <- unclass(loadings)
  M <- L %*% t(L)
  M[upper.tri(M)]
}

# baseline correlation matrix (3 factors) and a single-factor matrix of a different
# size; GRiPS is passed as a correlation matrix so the estimator, not the correlation
# construction, is what is compared.
reg_fixtures <- list(
  list(R = test_models$baseline$cormat, N = 500, k = 3),
  list(R = stats::cor(GRiPS_raw, use = "pairwise.complete.obs"), N = 810, k = 1)
)

test_that("PAF reproduces psych::fa(fm = 'pa')", {
  skip_on_cran()
  skip_if_not_installed("psych")

  for (fx in reg_fixtures) {
    efa <- suppressWarnings(EFA(fx$R, n_factors = fx$k, N = fx$N,
                                type = "psych", method = "PAF"))
    ref <- suppressMessages(suppressWarnings(
      psych::fa(fx$R, nfactors = fx$k, n.obs = fx$N, fm = "pa",
                rotate = "none", SMC = TRUE, max.iter = 50)))

    # observed agreement on the baseline model: ~1e-15
    expect_equal(unname(1 - efa$h2), unname(1 - ref$communality), tolerance = 1e-4)
    expect_equal(repro_offdiag(efa$unrot_loadings),
                 repro_offdiag(ref$loadings), tolerance = 1e-4)
  }
})

test_that("ML reproduces stats::factanal", {
  skip_on_cran()

  for (fx in reg_fixtures) {
    efa <- suppressWarnings(EFA(fx$R, n_factors = fx$k, N = fx$N, method = "ML"))
    ref <- stats::factanal(covmat = fx$R, factors = fx$k, n.obs = fx$N,
                           rotation = "none")

    # observed agreement on the baseline model: ~2e-6
    expect_equal(unname(1 - efa$h2), unname(ref$uniquenesses), tolerance = 1e-4)
    expect_equal(repro_offdiag(efa$unrot_loadings),
                 repro_offdiag(ref$loadings), tolerance = 1e-4)
  }
})

test_that("ULS reproduces psych::fa(fm = 'minres')", {
  skip_on_cran()
  skip_if_not_installed("psych")

  for (fx in reg_fixtures) {
    efa <- suppressWarnings(EFA(fx$R, n_factors = fx$k, N = fx$N, method = "ULS"))
    ref <- suppressMessages(suppressWarnings(
      psych::fa(fx$R, nfactors = fx$k, n.obs = fx$N, fm = "minres",
                rotate = "none")))

    # observed agreement on the baseline model: ~6e-6
    expect_equal(unname(1 - efa$h2), unname(1 - ref$communality), tolerance = 1e-4)
    expect_equal(repro_offdiag(efa$unrot_loadings),
                 repro_offdiag(ref$loadings), tolerance = 1e-4)
  }
})

rm(repro_offdiag, reg_fixtures)
