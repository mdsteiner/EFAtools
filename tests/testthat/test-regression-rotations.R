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

# shared fixture set for the native-engine relaxed-parity tests below (the non-convex CF,
# oblimin, and geomin criteria are exercised across three real correlation matrices of
# different factor counts so the Q-dominance contract is not pinned to a single solution)
fixtures <- list(
  baseline = list(R = test_models$baseline$cormat, N = 500,                nf = 3),
  GRiPS    = list(R = stats::cor(GRiPS_raw),        N = nrow(GRiPS_raw),    nf = 2),
  DOSPERT  = list(R = stats::cor(DOSPERT_raw),      N = nrow(DOSPERT_raw),  nf = 4)
)

test_that("orthogonal rotations reproduce their GPArotation engine", {
  skip_on_cran()
  skip_if_not_installed("GPArotation")

  # rotation name -> the GPArotation engine EFAtools dispatches to for each
  # GPArotation-backed orthogonal rotation. The Crawford-Ferguson criteria (quartimax,
  # equamax) and geomin (geominT) are computed by the native engine and validated separately
  # below.
  engines <- list(
    bentlerT  = function() GPArotation::bentlerT(L, normalize = TRUE, eps = 1e-5,
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

test_that("quartimax and equamax route through the native CF rotation engine", {
  skip_on_cran()
  skip_if_not_installed("GPArotation")

  # quartimax and equamax are computed by the native Crawford-Ferguson engine
  # (kappa = 0 and kappa = k / (2 p) respectively); the remaining orthogonal criteria
  # stay on GPArotation. CF is non-convex, so parity is relaxed rather than byte-exact:
  # the engine must reach an as-good-or-better minimum of the SAME criterion, and the
  # aligned loadings must agree to ~1e-4 rather than the ~1e-6 the convex engines reach.
  # The like-for-like reference is cfT() with the matching kappa -- GPArotation's
  # quartimax() optimises the same rotation but reports a differently scaled criterion
  # value, so cfT(kappa = 0) is the correct Q reference for quartimax.
  for (fx_name in names(fixtures)) {
    fx <- fixtures[[fx_name]]
    unrot_fx <- suppressWarnings(EFA(fx$R, n_factors = fx$nf, N = fx$N))
    Lx <- unclass(unrot_fx$unrot_loadings)
    p <- nrow(Lx)
    k <- ncol(Lx)

    for (rn in c("quartimax", "equamax")) {
      kappa <- if (rn == "quartimax") 0 else k / (2 * p)

      set.seed(seed)
      native <- .rotate_cf_orth(Lx, kappa = kappa, eps = 1e-5, normalize = TRUE,
                                random_starts = 100)
      set.seed(seed)
      ref <- GPArotation::cfT(Lx, kappa = kappa, normalize = TRUE, eps = 1e-5,
                              randomStarts = 100)
      ref_Q <- ref$Table[nrow(ref$Table), 2]

      # (a) as-good-or-better minimum of the same criterion (same Q definition)
      expect_lte(native$value, ref_Q + 1e-6)
      # (b) aligned loadings within relaxed parity
      expect_lt(aligned_max_diff(native$loadings, ref$loadings), 1e-4)

      # the public rotation path routes through the native engine and reproduces it
      set.seed(seed)
      efa <- suppressWarnings(.rotate_model(unrot_fx, rotation = rn, type = "EFAtools"))
      expect_lt(aligned_max_diff(efa$rot_loadings, ref$loadings), 1e-4)
    }
  }
})

test_that("geominT routes through the native geomin GPF engine", {
  skip_on_cran()
  skip_if_not_installed("GPArotation")

  # geominT is computed by the native gradient-projection engine; the remaining orthogonal
  # criteria stay on GPArotation. The geomin criterion is prone to local minima, so parity is
  # relaxed: the PRIMARY contract is that the native engine reaches an as-good-or-better minimum
  # of the SAME criterion. geominT's Table column 2 is the geomin value, computed on the same
  # Kaiser-normalized loadings the native engine reports, so it is directly comparable to
  # native$value. On these fixtures the criterion is well identified at randomStarts = 100, so
  # the aligned loadings also agree to ~1e-4; should an environment ever settle on a different
  # but equal-Q minimum, the Q-dominance check still passes and only the aligned-loadings check
  # would legitimately diverge.

  # the native engine pins the geomin offset at delta = 0.01, chosen to match GPArotation's
  # default; guard that the reference below is taken at the same offset
  expect_equal(formals(GPArotation::geominT)$delta, 0.01)

  for (fx_name in names(fixtures)) {
    fx <- fixtures[[fx_name]]
    unrot_fx <- suppressWarnings(EFA(fx$R, n_factors = fx$nf, N = fx$N))
    Lx <- unclass(unrot_fx$unrot_loadings)

    set.seed(seed)
    native <- .rotate_geomin_orth(Lx, delta = 0.01, eps = 1e-5, normalize = TRUE,
                                  random_starts = 100)
    set.seed(seed)
    ref <- GPArotation::geominT(Lx, delta = 0.01, normalize = TRUE, eps = 1e-5,
                                randomStarts = 100)
    ref_Q <- ref$Table[nrow(ref$Table), 2]

    # (a) as-good-or-better minimum of the same geomin criterion
    expect_lte(native$value, ref_Q + 1e-6)
    # (b) aligned loadings within relaxed parity (well-identified on these fixtures)
    expect_lt(aligned_max_diff(native$loadings, ref$loadings), 1e-4)

    # the public rotation path routes through the native engine and reproduces it exactly (same
    # seed, same compiled entry), so it inherits the Q-dominance guarantee above
    set.seed(seed)
    efa <- suppressWarnings(.rotate_model(unrot_fx, rotation = "geominT", type = "EFAtools"))
    expect_lt(aligned_max_diff(efa$rot_loadings, native$loadings), 1e-6)
  }
})

test_that("oblique rotations reproduce their GPArotation engine", {
  skip_on_cran()
  skip_if_not_installed("GPArotation")

  # rotation name -> the GPArotation engine EFAtools dispatches to for each
  # GPArotation-backed oblique rotation. The oblimin and quartimin criteria and geomin
  # (geominQ) are computed by the native engine and validated separately below.
  engines <- list(
    simplimax = function() GPArotation::simplimax(L, k = nrow(L), normalize = TRUE,
                                                  eps = 1e-5, randomStarts = 100),
    bentlerQ  = function() GPArotation::bentlerQ(L, normalize = TRUE, eps = 1e-5,
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

test_that("oblimin and quartimin route through the native oblique GPF engine", {
  skip_on_cran()
  skip_if_not_installed("GPArotation")

  # oblimin and quartimin are computed by the native gradient-projection engine (oblimin
  # with the GPArotation default gam = 0; quartimin is the same criterion); the remaining
  # oblique criteria stay on GPArotation. The criterion is non-convex, so parity is relaxed
  # rather than byte-exact: the engine must reach an as-good-or-better minimum of the SAME
  # criterion, the aligned loadings must agree to ~1e-4, and the factor correlations must
  # match under the permutation/sign-invariant fingerprint. The single GPArotation::oblimin()
  # reference serves both public names (quartimin is the same gam = 0 criterion).
  for (fx_name in names(fixtures)) {
    fx <- fixtures[[fx_name]]
    unrot_fx <- suppressWarnings(EFA(fx$R, n_factors = fx$nf, N = fx$N))
    Lx <- unclass(unrot_fx$unrot_loadings)

    set.seed(seed)
    native <- .rotate_oblimin(Lx, gam = 0, eps = 1e-5, normalize = TRUE,
                              random_starts = 100)
    set.seed(seed)
    ref <- suppressWarnings(GPArotation::oblimin(Lx, normalize = TRUE, eps = 1e-5,
                                                 randomStarts = 100))
    ref_Q <- ref$Table[nrow(ref$Table), 2]

    # the native engine reaches an as-good-or-better minimum of the same criterion, with
    # aligned loadings and matching factor correlations
    expect_lte(native$value, ref_Q + 1e-6)
    expect_lt(aligned_max_diff(native$loadings, ref$loadings), 1e-4)
    expect_equal(phi_fingerprint(native$Phi), phi_fingerprint(ref$Phi), tolerance = 1e-4)

    # both public rotation names dispatch to the native engine and reproduce the reference
    for (rn in c("oblimin", "quartimin")) {
      set.seed(seed)
      efa <- suppressWarnings(.rotate_model(unrot_fx, rotation = rn, type = "EFAtools"))
      expect_lt(aligned_max_diff(efa$rot_loadings, ref$loadings), 1e-4)
      expect_equal(phi_fingerprint(efa$Phi), phi_fingerprint(ref$Phi), tolerance = 1e-4)
    }
  }
})

test_that("the native oblimin engine generalizes to non-zero gamma", {
  skip_on_cran()
  skip_if_not_installed("GPArotation")

  # The package only ever uses gam = 0 (oblimin == quartimin), but the native engine
  # implements the full oblimin family parameterized by gamma. Exercise the gamma-centering
  # path with gam = 0.5 (biquartimin) against GPArotation's oblimin at the same gamma, under
  # the same relaxed-parity scheme (Q-dominance, aligned loadings, Phi fingerprint).
  set.seed(seed)
  native <- .rotate_oblimin(L, gam = 0.5, eps = 1e-5, normalize = TRUE, random_starts = 100)
  set.seed(seed)
  ref <- suppressWarnings(GPArotation::oblimin(L, gam = 0.5, normalize = TRUE, eps = 1e-5,
                                               randomStarts = 100))

  expect_lte(native$value, ref$Table[nrow(ref$Table), 2] + 1e-6)
  expect_lt(aligned_max_diff(native$loadings, ref$loadings), 1e-4)
  expect_equal(phi_fingerprint(native$Phi), phi_fingerprint(ref$Phi), tolerance = 1e-4)
})

test_that("geominQ routes through the native geomin GPF engine", {
  skip_on_cran()
  skip_if_not_installed("GPArotation")

  # geominQ is computed by the native gradient-projection engine; the remaining oblique criteria
  # stay on GPArotation. As with geominT, the geomin criterion is prone to local minima, so the
  # PRIMARY contract is an as-good-or-better minimum of the SAME criterion (geominQ's Table
  # column 2, directly comparable to native$value) rather than byte-exact reproduction. On these
  # fixtures the criterion is well identified at randomStarts = 100, so the aligned loadings and
  # the factor correlations also agree; where they diverged at equal Q, the Q-dominance check
  # would still hold.

  # the native engine pins the geomin offset at delta = 0.01, chosen to match GPArotation's
  # default; guard that the reference below is taken at the same offset
  expect_equal(formals(GPArotation::geominQ)$delta, 0.01)

  for (fx_name in names(fixtures)) {
    fx <- fixtures[[fx_name]]
    unrot_fx <- suppressWarnings(EFA(fx$R, n_factors = fx$nf, N = fx$N))
    Lx <- unclass(unrot_fx$unrot_loadings)

    set.seed(seed)
    native <- .rotate_geomin_oblq(Lx, delta = 0.01, eps = 1e-5, normalize = TRUE,
                                  random_starts = 100)
    set.seed(seed)
    ref <- suppressWarnings(GPArotation::geominQ(Lx, delta = 0.01, normalize = TRUE,
                                                 eps = 1e-5, randomStarts = 100))
    ref_Q <- ref$Table[nrow(ref$Table), 2]

    # (a) as-good-or-better minimum of the same geomin criterion
    expect_lte(native$value, ref_Q + 1e-6)
    # (b) aligned loadings and factor correlations within relaxed parity (well-identified here)
    expect_lt(aligned_max_diff(native$loadings, ref$loadings), 1e-4)
    expect_equal(phi_fingerprint(native$Phi), phi_fingerprint(ref$Phi), tolerance = 1e-4)

    # the public rotation path routes through the native engine and reproduces it exactly (same
    # seed, same compiled entry), so it inherits the Q-dominance guarantee above
    set.seed(seed)
    efa <- suppressWarnings(.rotate_model(unrot_fx, rotation = "geominQ", type = "EFAtools"))
    expect_lt(aligned_max_diff(efa$rot_loadings, native$loadings), 1e-6)
    expect_equal(phi_fingerprint(efa$Phi), phi_fingerprint(native$Phi), tolerance = 1e-6)
  }
})

test_that("varimax reproduces stats::varimax", {
  skip_on_cran()

  efa <- suppressWarnings(.rotate_model(unrot, rotation = "varimax", type = "psych"))   # varimax_type 'svd'
  ref <- stats::varimax(L, normalize = TRUE, eps = 1e-5)
  expect_lt(aligned_max_diff(efa$rot_loadings, ref$loadings), 1e-4)
})

test_that("kaiser varimax reproduces stats::varimax", {
  skip_on_cran()

  # The svd path is pinned above; this targets the Kaiser successive-pairwise
  # algorithm (varimax_type "kaiser", used by type "EFAtools"/"SPSS"). It optimises
  # the same Kaiser-normalised varimax criterion as stats::varimax's SVD updates but
  # via a different engine, so the loadings agree to ~5e-4 rather than exactly; this
  # pins the rotated loadings only, not the .SV criterion-monitoring path. The bound
  # sits well below a genuine rotation regression (which shifts loadings by >=1e-2)
  # while leaving headroom over the cross-engine difference for BLAS portability.
  ref <- stats::varimax(L, normalize = TRUE, eps = 1e-5)
  for (ty in c("EFAtools", "SPSS")) {
    efa <- suppressWarnings(.rotate_model(unrot, rotation = "varimax", type = ty))
    expect_lt(aligned_max_diff(efa$rot_loadings, ref$loadings), 2e-3)
  }
})

test_that("promax reproduces stats::promax", {
  skip_on_cran()

  # type 'psych' uses the unnormalized Hendrickson-White target (k = 4) on an svd varimax
  # base, which is what stats::promax(m = 4) computes.
  efa <- suppressWarnings(.rotate_model(unrot, rotation = "promax", type = "psych"))
  ref <- stats::promax(L, m = 4)
  expect_lt(aligned_max_diff(efa$rot_loadings, ref$loadings), 1e-4)
})

rm(aligned_max_diff, phi_fingerprint, unrot, L, seed, fixtures)
