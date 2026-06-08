# Regression net: pin the current rotation output against the GPArotation / stats engine
# each rotation is designed to wrap, so an accidental change (or a re-implementation of the
# rotation engine) is caught. The rotation is isolated on a single fixed unrotated solution
# and compared against the oracle on that same loading matrix. Sign and column-permutation
# differences are removed with the package's own COMPARE() (Tucker-congruence reorder plus
# sign reflection); the oblique factor correlations are compared via a permutation/sign
# invariant fingerprint. For the GPArotation engines the oracle is given the same
# randomStarts = 100 under a shared seed as EFAtools uses, which makes the local-minimum
# prone criteria (bifactor, simplimax) align deterministically. Tolerance is generous
# relative to the observed agreement (exact for every rotation here) so the contract is
# portable.

# max absolute difference after aligning columns/signs of two loading matrices
aligned_max_diff <- function(x, y) {
  COMPARE(x, y, reorder = "congruence", corres = FALSE,
          plot = FALSE, print_diff = FALSE)$max_abs_diff
}

# permutation/sign-invariant fingerprint of the factor intercorrelations
phi_fingerprint <- function(Phi) sort(abs(Phi[upper.tri(Phi)]))

unrot <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500)
L <- unclass(unrot$unrot_loadings)
seed <- 42

test_that("orthogonal rotations reproduce their GPArotation engine", {
  skip_on_cran()
  skip_if_not_installed("GPArotation")

  # rotation name -> the GPArotation engine EFAtools currently dispatches to. quartimax
  # intentionally maps to bentlerT: that is the engine EFAtools calls for orthogonal
  # quartimax today, so the regression net pins the current behaviour.
  engines <- list(
    equamax   = function() GPArotation::cfT(L, kappa = ncol(L) / (2 * nrow(L)),
                                            normalize = TRUE, eps = 1e-5,
                                            randomStarts = 100),
    bentlerT  = function() GPArotation::bentlerT(L, normalize = TRUE, eps = 1e-5,
                                                 randomStarts = 100),
    quartimax = function() GPArotation::bentlerT(L, normalize = TRUE, eps = 1e-5,
                                                 randomStarts = 100),
    geominT   = function() GPArotation::geominT(L, normalize = TRUE, eps = 1e-5,
                                                randomStarts = 100),
    bifactorT = function() GPArotation::bifactorT(L, normalize = TRUE, eps = 1e-5,
                                                  randomStarts = 100)
  )

  for (rn in names(engines)) {
    set.seed(seed)
    efa <- suppressWarnings(.rotate_model(unrot, rotation = rn, type = "EFAtools"))
    set.seed(seed)
    ref <- suppressWarnings(engines[[rn]]())
    expect_lt(aligned_max_diff(efa$rot_loadings, ref$loadings), 1e-4)
  }
})

test_that("oblique rotations reproduce their GPArotation engine", {
  skip_on_cran()
  skip_if_not_installed("GPArotation")

  engines <- list(
    oblimin   = function() GPArotation::oblimin(L, normalize = TRUE, eps = 1e-5,
                                                randomStarts = 100),
    quartimin = function() GPArotation::quartimin(L, normalize = TRUE, eps = 1e-5,
                                                  randomStarts = 100),
    simplimax = function() GPArotation::simplimax(L, k = nrow(L), normalize = TRUE,
                                                  eps = 1e-5, randomStarts = 100),
    bentlerQ  = function() GPArotation::bentlerQ(L, normalize = TRUE, eps = 1e-5,
                                                 randomStarts = 100),
    geominQ   = function() GPArotation::geominQ(L, normalize = TRUE, eps = 1e-5,
                                                randomStarts = 100),
    bifactorQ = function() GPArotation::bifactorQ(L, normalize = TRUE, eps = 1e-5,
                                                  randomStarts = 100)
  )

  for (rn in names(engines)) {
    set.seed(seed)
    efa <- suppressWarnings(.rotate_model(unrot, rotation = rn, type = "EFAtools"))
    set.seed(seed)
    ref <- suppressWarnings(engines[[rn]]())
    expect_lt(aligned_max_diff(efa$rot_loadings, ref$loadings), 1e-4)
    expect_equal(phi_fingerprint(efa$Phi), phi_fingerprint(ref$Phi), tolerance = 1e-4)
  }
})

test_that("varimax reproduces stats::varimax", {
  skip_on_cran()

  efa <- suppressWarnings(.rotate_model(unrot, rotation = "varimax", type = "psych"))   # varimax_type 'svd'
  ref <- stats::varimax(L, normalize = TRUE, eps = 1e-5)
  expect_lt(aligned_max_diff(efa$rot_loadings, ref$loadings), 1e-4)
})

test_that("promax reproduces stats::promax", {
  skip_on_cran()

  # type 'psych' uses the unnormalized Hendrickson-White target (k = 4) on an svd varimax
  # base, which is what stats::promax(m = 4) computes.
  efa <- suppressWarnings(.rotate_model(unrot, rotation = "promax", type = "psych"))
  ref <- stats::promax(L, m = 4)
  expect_lt(aligned_max_diff(efa$rot_loadings, ref$loadings), 1e-4)
})

rm(aligned_max_diff, phi_fingerprint, unrot, L, seed)
