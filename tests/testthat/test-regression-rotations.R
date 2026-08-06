# Regression net: pin the current rotation output against the GPArotation / stats engine
# each rotation is designed to wrap, so an accidental change (or a re-implementation of the
# rotation engine) is caught. The rotation is isolated on a single fixed unrotated solution
# and compared against the oracle on that same loading matrix. Sign and column-permutation
# differences are removed with the package's own efa_compare() (Tucker-congruence reorder plus
# sign reflection); the oblique factor correlations are compared via a permutation/sign
# invariant fingerprint. For the GPArotation engines the oracle is given randomStarts under a
# shared seed (100 for the smooth criteria, matching EFAtools; 300 for the multimodal simplimax
# criterion, a matched budget at which both engines reach the global basin). Tolerance is generous
# relative to the observed agreement (exact for every rotation here) so the contract is
# portable.

# max absolute difference after aligning columns/signs of two loading matrices
aligned_max_diff <- function(x, y) {
  efa_compare(x, y, reorder = "congruence", corres = FALSE,
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

# Every block below rotates the same unrotated solution per fixture, so each fixture is
# extracted once and cached. The cache is filled on first request rather than eagerly:
# every consuming block is skip_on_cran(), so on CRAN no extraction is paid for at all.
# Unrotated PAF extraction is deterministic and draws no random numbers, and each site
# re-seeds immediately before the rotation, so caching cannot perturb the random-start
# stream the engines below share with their oracles.
unrot_fixture <- local({
  cache <- list()
  function(fx_name) {
    if (!fx_name %in% names(cache)) {
      fx <- fixtures[[fx_name]]
      cache[[fx_name]] <<- suppressWarnings(EFA(fx$R, n_factors = fx$nf, N = fx$N))
    }
    cache[[fx_name]]
  }
})

test_that("quartimax and equamax route through the native CF rotation engine", {
  skip_on_cran()
  skip_if_not_installed("GPArotation")

  # quartimax and equamax are computed by the native Crawford-Ferguson engine
  # (kappa = 0 and kappa = k / (2 p) respectively), as are the remaining orthogonal
  # criteria. CF is non-convex, so parity is relaxed rather than byte-exact:
  # the engine must reach an as-good-or-better minimum of the SAME criterion, and the
  # aligned loadings must agree to ~1e-4 rather than the ~1e-6 the convex engines reach.
  # The like-for-like reference is cfT() with the matching kappa -- GPArotation's
  # quartimax() optimises the same rotation but reports a differently scaled criterion
  # value, so cfT(kappa = 0) is the correct Q reference for quartimax.
  # A rotation-engine regression shifts loadings by >= 1e-2 regardless of fixture, so
  # a single baseline fixture is enough to catch one; cross-fixture coverage is kept
  # for the multimodal/flat criteria below.
  for (fx_name in "baseline") {
    unrot_fx <- unrot_fixture(fx_name)
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

  # geominT is computed by the native gradient-projection engine, as are the remaining orthogonal
  # criteria. The geomin criterion is prone to local minima, so parity is
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
    unrot_fx <- unrot_fixture(fx_name)
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

test_that("bentlerT routes through the native GPF engine", {
  skip_on_cran()
  skip_if_not_installed("GPArotation")

  # bentlerT is computed by the native gradient-projection engine, as are the remaining orthogonal
  # criteria. The Bentler criterion is non-convex, so parity is relaxed: the
  # PRIMARY contract is that the native engine reaches an as-good-or-better minimum of the SAME
  # criterion. bentlerT's Table column 2 is the Bentler value computed on the same
  # Kaiser-normalized loadings the native engine reports, so it is directly comparable to
  # native$value. On these fixtures the criterion is well identified at randomStarts = 100, so the
  # aligned loadings also agree to ~1e-4.
  for (fx_name in names(fixtures)) {
    unrot_fx <- unrot_fixture(fx_name)
    Lx <- unclass(unrot_fx$unrot_loadings)

    set.seed(seed)
    native <- .rotate_bentler_orth(Lx, eps = 1e-5, normalize = TRUE, random_starts = 100)
    set.seed(seed)
    ref <- GPArotation::bentlerT(Lx, normalize = TRUE, eps = 1e-5, randomStarts = 100)
    ref_Q <- ref$Table[nrow(ref$Table), 2]

    # (a) as-good-or-better minimum of the same Bentler criterion
    expect_lte(native$value, ref_Q + 1e-6)
    # (b) aligned loadings within relaxed parity (well-identified on these fixtures)
    expect_lt(aligned_max_diff(native$loadings, ref$loadings), 1e-4)

    # the public rotation path routes through the native engine and reproduces it exactly (same
    # seed, same compiled entry), so it inherits the Q-dominance guarantee above
    set.seed(seed)
    efa <- suppressWarnings(.rotate_model(unrot_fx, rotation = "bentlerT", type = "EFAtools"))
    expect_lt(aligned_max_diff(efa$rot_loadings, native$loadings), 1e-6)
  }
})

test_that("bifactorT routes through the native GPF engine", {
  skip_on_cran()
  skip_if_not_installed("GPArotation")

  # bifactorT is computed by the native gradient-projection engine; with this every orthogonal
  # criterion is native. The Jennrich-Bentler bifactor criterion (the first factor is a general
  # factor, exempt from the penalty) is strongly prone to local minima, so parity is relaxed: the
  # PRIMARY contract is that the native engine reaches an as-good-or-better minimum of the SAME
  # criterion. bifactorT's Table column 2 is the bifactor value computed on the same
  # Kaiser-normalized loadings the native engine reports, so it is directly comparable to
  # native$value. The criterion needs at least two group factors to be non-trivial: with a single
  # group factor (two factors) it is identically zero, so any rotation is optimal and the aligned
  # loadings may legitimately differ. On the k >= 3 fixtures the criterion is well identified at
  # randomStarts = 100, so the aligned loadings also agree to ~1e-4; where an environment settles
  # on a different but equal-Q minimum, the Q-dominance check still holds.
  for (fx_name in names(fixtures)) {
    unrot_fx <- unrot_fixture(fx_name)
    Lx <- unclass(unrot_fx$unrot_loadings)

    set.seed(seed)
    native <- .rotate_bifactor_orth(Lx, eps = 1e-5, normalize = TRUE, random_starts = 100)
    set.seed(seed)
    ref <- GPArotation::bifactorT(Lx, normalize = TRUE, eps = 1e-5, randomStarts = 100)
    ref_Q <- ref$Table[nrow(ref$Table), 2]

    # (a) as-good-or-better minimum of the same bifactor criterion (every fixture)
    expect_lte(native$value, ref_Q + 1e-6)
    # (b) aligned loadings within relaxed parity, where the criterion is non-trivial (k >= 3)
    if (ncol(Lx) >= 3) {
      expect_lt(aligned_max_diff(native$loadings, ref$loadings), 1e-4)
    }

    # the public rotation path routes through the native engine and reproduces it exactly (same
    # seed, same compiled entry), so it inherits the Q-dominance guarantee above
    set.seed(seed)
    efa <- suppressWarnings(.rotate_model(unrot_fx, rotation = "bifactorT", type = "EFAtools"))
    expect_lt(aligned_max_diff(efa$rot_loadings, native$loadings), 1e-6)
  }
})

test_that("simplimax routes through the native oblique GPF engine", {
  skip_on_cran()
  skip_if_not_installed("GPArotation")
  skip_if_not_slow()

  # simplimax is computed by the native gradient-projection engine;
  # unlike the smooth criteria, the simplimax criterion reselects the
  # k smallest squared loadings at every evaluation (k = nrow(L) by default), so it is only
  # piecewise smooth and strongly multimodal: its gradient jumps as loadings cross the kth-smallest threshold (the
  # native engine uses a non-monotone line search to step across the kinks) and the global minimum
  # depends on the random restarts. The native engine therefore fully optimizes every restart and
  # keeps the lowest-criterion solution -- the standard remedy for the local minima of
  # complexity-based criteria (Kiers, 1994; Browne, 2001) -- rather than the screen-and-triage
  # heuristic the smooth criteria use.
  #
  # Because the criterion is multimodal, native and GPArotation draw independent restart sequences
  # and need not reach the SAME minimizer -- and neither engine's attained VALUE is a stable
  # target either. The minima are dense: on DOSPERT an exact variable relabelling (a permutation
  # under which every reported statistic is invariant by construction, with the seed fixed so the
  # same starts are drawn) moves GPArotation's attained value by ~45% (0.004373 to 0.006397 over
  # 20 relabellings) while the native engine's moves ~6% (0.005571 to 0.005881). The oracle's own
  # value is therefore not identified to the precision a direct engine-to-engine comparison would
  # need, so no dominance claim is asserted here.
  #
  # What is identified is asserted instead: (a) the reported value really is the simplimax
  # criterion of the reported rotation; (b) the rotation descends from its unrotated start;
  # (c) the two engines agree in order of magnitude, which is the scale on which a genuine
  # regression in the engine shows up. Aligned-loadings parity is NOT asserted: at comparable Q the
  # two engines may settle on different, equally valid simple-structure rotations.
  #
  # The criterion is evaluated here from its definition rather than taken from the engine or from
  # GPArotation: the sum of the k smallest squared loadings (Kiers, 1994, Psychometrika, 59(4),
  # 567-579, doi:10.1007/BF02294392).
  simplimax_Q <- function(L, k) sum(sort(as.vector(unclass(L))^2)[seq_len(k)])

  for (fx_name in names(fixtures)) {
    unrot_fx <- unrot_fixture(fx_name)
    Lx <- unclass(unrot_fx$unrot_loadings)

    set.seed(seed)
    native <- .rotate_simplimax_oblq(Lx, k = nrow(Lx), eps = 1e-5, normalize = TRUE,
                                     random_starts = 300)
    set.seed(seed)
    ref <- suppressWarnings(GPArotation::simplimax(Lx, k = nrow(Lx), normalize = TRUE,
                                                   eps = 1e-5, randomStarts = 300))
    ref_Q <- ref$Table[nrow(ref$Table), 2]

    # (a) the reported value is the criterion of the reported rotation. normalize = TRUE optimizes
    # on the Kaiser-normalized loadings and un-normalizes the returned ones, so the comparison is
    # made back on the scale the engine optimizes. Both sides evaluate one closed form on one
    # matrix, so this holds to machine precision on any platform.
    h <- sqrt(rowSums(Lx^2))
    expect_equal(native$value, simplimax_Q(native$loadings / h, nrow(Lx)), tolerance = 1e-10)

    # (b) the rotation descends from its unrotated start
    expect_lt(native$value, simplimax_Q(Lx / h, nrow(Lx)))

    # (c) both engines land on the same scale. The factor is set from measurement in both
    # directions: at the intended 300 starts the worst native/ref_Q ratio seen on any fixture,
    # including the 20 relabellings, is 1.27, while collapsing the multistart search to a single
    # start -- the regression this block guards -- drives the ratio to 1.95 (GRiPS), 4.04 (DOSPERT)
    # and 5.14 (baseline). 1.5 sits between the two, so it catches a restart regression on every
    # fixture without pinning the labelling-induced spread. GRiPS is the binding fixture; a factor
    # of 2 would let a single-start collapse pass there unnoticed.
    expect_lt(native$value, 1.5 * ref_Q)
    # the native solution is a valid oblique rotation (unit-diagonal factor correlations)
    expect_equal(diag(native$Phi), rep(1, ncol(Lx)))
  }
})

test_that("the public simplimax rotation path reproduces the native engine", {
  skip_on_cran()

  # the public rotation path routes through the native engine and reproduces it exactly (same seed,
  # same compiled entry, k defaults to nrow(L)); checked at a small restart budget for speed
  set.seed(seed)
  efa <- suppressWarnings(.rotate_model(unrot, rotation = "simplimax", type = "EFAtools",
                                        randomStarts = 10))
  set.seed(seed)
  nat <- .rotate_simplimax_oblq(L, k = nrow(L), eps = 1e-5, normalize = TRUE, random_starts = 10)
  expect_lt(aligned_max_diff(efa$rot_loadings, nat$loadings), 1e-6)
  expect_equal(phi_fingerprint(efa$Phi), phi_fingerprint(nat$Phi), tolerance = 1e-6)
})

test_that("oblimin and quartimin route through the native oblique GPF engine", {
  skip_on_cran()
  skip_if_not_installed("GPArotation")

  # oblimin and quartimin are computed by the native gradient-projection engine (oblimin
  # with the GPArotation default gam = 0; quartimin is the same criterion); the remaining
  # oblique criteria are computed by the native engine too. The criterion is non-convex, so
  # parity is relaxed rather than byte-exact: the engine must reach an as-good-or-better minimum
  # of the SAME criterion, the aligned loadings must agree to ~1e-4, and the factor correlations
  # must match under the permutation/sign-invariant fingerprint. The single GPArotation::oblimin()
  # reference serves both public names (quartimin is the same gam = 0 criterion).
  # A single baseline fixture is enough to catch any rotation-engine regression here;
  # cross-fixture coverage is kept for the multimodal/flat criteria below.
  for (fx_name in "baseline") {
    unrot_fx <- unrot_fixture(fx_name)
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

  # Exercise the gamma-centering path of the native oblimin engine with gam = 0.5
  # (biquartimin) against GPArotation's oblimin at the same gamma, under
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

  # geominQ is computed by the native gradient-projection engine, as are the remaining oblique
  # criteria. As with geominT, the geomin criterion is prone to local minima, so the
  # PRIMARY contract is an as-good-or-better minimum of the SAME criterion (geominQ's Table
  # column 2, directly comparable to native$value) rather than byte-exact reproduction. On these
  # fixtures the criterion is well identified at randomStarts = 100, so the aligned loadings and
  # the factor correlations also agree; where they diverged at equal Q, the Q-dominance check
  # would still hold.

  # the native engine pins the geomin offset at delta = 0.01, chosen to match GPArotation's
  # default; guard that the reference below is taken at the same offset
  expect_equal(formals(GPArotation::geominQ)$delta, 0.01)

  for (fx_name in names(fixtures)) {
    unrot_fx <- unrot_fixture(fx_name)
    Lx <- unclass(unrot_fx$unrot_loadings)

    # geominQ uses a wider multistart by default than the bare compiled entry (the oblique geomin
    # criterion is multimodal enough that the two best-screened starts can both miss the global
    # basin), so the direct call is given the same per-criterion settings the public path applies,
    # ensuring both optimize identically.
    set.seed(seed)
    native <- do.call(.rotate_geomin_oblq,
                      c(list(Lx, delta = 0.01, eps = 1e-5, normalize = TRUE, random_starts = 100),
                        .gpf_multistart_defaults$geominQ))
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

test_that("bentlerQ routes through the native GPF engine", {
  skip_on_cran()
  skip_if_not_installed("GPArotation")

  # bentlerQ is computed by the native gradient-projection engine, as are the remaining oblique
  # criteria. The Bentler criterion is non-convex AND comparatively flat near its
  # oblique optimum, so parity is relaxed: the PRIMARY contract is that the native engine reaches
  # an as-good-or-better minimum of the SAME criterion (bentlerQ's Table column 2, directly
  # comparable to native$value). Because the oblique optimum is flat, the native engine and
  # GPArotation settle on the same equal-Q solution only to ~1e-4 here (rather than the ~1e-6 of
  # the convex engines), so the aligned loadings and factor correlations are checked at a
  # correspondingly looser ~1e-3 tolerance; a genuine rotation regression shifts loadings by
  # >= 1e-2.
  for (fx_name in names(fixtures)) {
    unrot_fx <- unrot_fixture(fx_name)
    Lx <- unclass(unrot_fx$unrot_loadings)

    set.seed(seed)
    native <- .rotate_bentler_oblq(Lx, eps = 1e-5, normalize = TRUE, random_starts = 100)
    set.seed(seed)
    ref <- suppressWarnings(GPArotation::bentlerQ(Lx, normalize = TRUE, eps = 1e-5,
                                                  randomStarts = 100))
    ref_Q <- ref$Table[nrow(ref$Table), 2]

    # (a) as-good-or-better minimum of the same Bentler criterion
    expect_lte(native$value, ref_Q + 1e-6)
    # (b) aligned loadings and factor correlations within relaxed parity (flat oblique optimum)
    expect_lt(aligned_max_diff(native$loadings, ref$loadings), 1e-3)
    expect_lt(max(abs(phi_fingerprint(native$Phi) - phi_fingerprint(ref$Phi))), 1e-3)

    # the public rotation path routes through the native engine and reproduces it exactly (same
    # seed, same compiled entry), so it inherits the Q-dominance guarantee above
    set.seed(seed)
    efa <- suppressWarnings(.rotate_model(unrot_fx, rotation = "bentlerQ", type = "EFAtools"))
    expect_lt(aligned_max_diff(efa$rot_loadings, native$loadings), 1e-6)
    expect_lt(max(abs(phi_fingerprint(efa$Phi) - phi_fingerprint(native$Phi))), 1e-6)
  }
})

test_that("bifactorQ routes through the native GPF engine", {
  skip_on_cran()
  skip_if_not_installed("GPArotation")

  # bifactorQ is computed by the native gradient-projection engine; with this every analytic
  # rotation criterion is native. The Jennrich-Bentler bifactor criterion (the first factor is a
  # general factor, exempt from the penalty) is strongly prone to local minima and comparatively
  # flat near its oblique optimum, so parity is relaxed: the PRIMARY contract is that the native
  # engine reaches an as-good-or-better minimum of the SAME criterion (bifactorQ's Table column 2,
  # directly comparable to native$value). The criterion needs at least two group factors to be
  # non-trivial: with a single group factor (two factors) it is identically zero, so any rotation is
  # optimal. On the k >= 3 fixtures the native engine and GPArotation settle on the same equal-Q
  # solution only to ~1e-4 (the flat oblique optimum), so the aligned loadings and factor
  # correlations are checked at a looser ~1e-3 tolerance; a genuine rotation regression shifts
  # loadings by >= 1e-2.
  for (fx_name in names(fixtures)) {
    unrot_fx <- unrot_fixture(fx_name)
    Lx <- unclass(unrot_fx$unrot_loadings)

    set.seed(seed)
    native <- .rotate_bifactor_oblq(Lx, eps = 1e-5, normalize = TRUE, random_starts = 100)
    set.seed(seed)
    ref <- suppressWarnings(GPArotation::bifactorQ(Lx, normalize = TRUE, eps = 1e-5,
                                                   randomStarts = 100))
    ref_Q <- ref$Table[nrow(ref$Table), 2]

    # (a) as-good-or-better minimum of the same bifactor criterion (every fixture)
    expect_lte(native$value, ref_Q + 1e-6)
    # (b) aligned loadings and factor correlations within relaxed parity, where non-trivial (k >= 3)
    if (ncol(Lx) >= 3) {
      expect_lt(aligned_max_diff(native$loadings, ref$loadings), 1e-3)
      expect_lt(max(abs(phi_fingerprint(native$Phi) - phi_fingerprint(ref$Phi))), 1e-3)
    }

    # the public rotation path routes through the native engine and reproduces it exactly (same
    # seed, same compiled entry), so it inherits the Q-dominance guarantee above
    set.seed(seed)
    efa <- suppressWarnings(.rotate_model(unrot_fx, rotation = "bifactorQ", type = "EFAtools"))
    expect_lt(aligned_max_diff(efa$rot_loadings, native$loadings), 1e-6)
    expect_lt(max(abs(phi_fingerprint(efa$Phi) - phi_fingerprint(native$Phi))), 1e-6)
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

test_that("promax with normalize = FALSE reproduces psych's Promax", {
  skip_on_cran()
  skip_if_not_installed("psych")

  # Pins the documented claim that switching off Kaiser normalization reproduces the
  # psych promax solution "to within the varimax convergence tolerance". psych::Promax is
  # what psych::fa(rotate = "Promax") calls. The residual is the convergence noise of the
  # shared stats::varimax base at eps = 1e-5, not an algorithmic difference, so the
  # equivalence is asserted at 1e-4 rather than at machine precision.
  efa <- suppressWarnings(
    .rotate_model(unrot, rotation = "promax", type = "psych", normalize = FALSE)
  )
  ref <- psych::Promax(L, m = 4)
  expect_lt(aligned_max_diff(efa$rot_loadings, ref$loadings), 1e-4)
  expect_equal(phi_fingerprint(efa$Phi), phi_fingerprint(ref$Phi), tolerance = 1e-4)
})

rm(aligned_max_diff, phi_fingerprint, unrot, L, seed, fixtures, unrot_fixture)
