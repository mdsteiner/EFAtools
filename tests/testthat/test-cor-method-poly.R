# cor_method = "poly" / "tetra": polychoric and tetrachoric correlations wired into
# the public API. The matrix estimator itself is covered in test-polychoric.R; these
# tests cover the routing through .prepare_cor_input() (EFA and the suitability /
# retention functions), the `use` handling, the EFA() bootstrap recompute, the
# tetrachoric binary assertion, and the clean rejection by the criteria whose
# reference data are continuous (CD, PARALLEL, NEST, HULL).

# Muffle only the named warning classes, letting any other (unexpected) warning
# surface so it is not silently hidden by a blanket suppressWarnings().
.muffle <- function(expr, ...) {
  classes <- c(...)
  withCallingHandlers(expr, warning = function(w) {
    if (inherits(w, classes)) invokeRestart("muffleWarning")
  })
}

# Collect the class strings of every warning raised while evaluating `expr` (several fire on
# the degenerate path), so a single expected class can be asserted without expect_warning
# stopping at the first one.
.warn_classes <- function(expr) {
  classes <- character()
  withCallingHandlers(
    expr,
    warning = function(w) { classes <<- c(classes, class(w)); invokeRestart("muffleWarning") }
  )
  classes
}

# A six-item ordinal set whose first pair is at a Frechet bound (its response table is the
# comonotone coupling of its own margins), so that pair is estimated at the boundary of the
# parameter space and has no asymptotic variance; the other four items are ordinary filler.
# Deterministic: a fixed contingency table expanded to raw data, plus seeded filler.
.degenerate_ordinal <- function() {
  tab <- rbind(c(41, 0), c(13, 0), c(28, 0), c(171, 47))
  pair <- do.call(rbind, lapply(which(tab > 0), function(k) {
    a <- ((k - 1L) %% nrow(tab)) + 1L
    b <- ((k - 1L) %/% nrow(tab)) + 1L
    cbind(rep(a, tab[a, b]), rep(b, tab[a, b]))
  }))
  set.seed(11)
  filler <- matrix(sample(1:4, nrow(pair) * 4L, TRUE), nrow(pair), 4L)
  out <- cbind(pair, filler)
  colnames(out) <- paste0("i", seq_len(ncol(out)))
  out
}

test_that("EFA() runs ordinal factor analysis and matches psych::fa on a polychoric matrix", {
  skip_on_cran()
  skip_if_not_installed("psych")

  x <- DOSPERT_raw[stats::complete.cases(DOSPERT_raw), ]
  N <- nrow(x)
  k <- 3L

  efa <- EFA(x, n_factors = k, method = "ULS", cor_method = "poly",
             rotation = "none")

  # Reference: psych's minres factor analysis (ULS is the same estimator) of psych's
  # own polychoric matrix, with correct = FALSE so psych applies no empty-cell
  # continuity correction, matching .polychoric().
  Rp <- suppressWarnings(psych::polychoric(x, correct = FALSE, global = FALSE)$rho)
  fa_ref <- suppressWarnings(psych::fa(Rp, nfactors = k, fm = "minres",
                                       rotate = "none", n.obs = N))
  L_ref <- unclass(fa_ref$loadings)
  L <- efa$unrot_loadings

  # Compare via the reproduced common-factor correlations (L %*% t(L)), which are
  # invariant to factor sign and column order, so no fragile per-column alignment
  # is needed. The two polychoric matrices agree to ~1e-4 and ULS == minres, so the
  # common parts and communalities agree well within 1e-3.
  expect_lt(max(abs(tcrossprod(L) - tcrossprod(L_ref))), 1e-3)
  expect_lt(max(abs(efa$h2 - fa_ref$communality)), 1e-3)
})

test_that("cor_method = 'tetra' equals 'poly' on binary data", {
  g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), ]
  gb <- apply(g, 2L, function(col) as.integer(col > stats::median(col)))

  # A Heywood case (and any non-PD smoothing) on this small binary example is
  # incidental to the equivalence being checked; muffle only those specific
  # warnings rather than blanket-suppressing every condition.
  e_poly  <- .muffle(EFA(gb, n_factors = 2, method = "ULS",
                         cor_method = "poly", rotation = "none"),
                     "efa_heywood", "efa_cor_smoothed")
  e_tetra <- .muffle(EFA(gb, n_factors = 2, method = "ULS",
                         cor_method = "tetra", rotation = "none"),
                     "efa_heywood", "efa_cor_smoothed")

  expect_equal(e_tetra$unrot_loadings, e_poly$unrot_loadings)
  expect_equal(e_tetra$orig_R, e_poly$orig_R)
})

test_that("cor_method = 'tetra' rejects variables with more than two categories", {
  expect_error(EFA(DOSPERT_raw, n_factors = 3, cor_method = "tetra"),
               class = "efa_cor_not_binary")
})

test_that("a Frechet-bound polychoric pair refuses DWLS and withholds the sandwich SEs", {
  x <- .degenerate_ordinal()

  # The first pair's response table is at a Frechet bound, so it has no asymptotic variance (see
  # test-polychoric.R). DWLS weights each element by the inverse of that variance, so it cannot be
  # formed and the fit is refused rather than run on a weight invented from rounding noise -- which
  # previously either dropped the pair from the solution or let it dominate, depending on the
  # platform.
  expect_error(
    suppressWarnings(
      efa_fit(x, n_factors = 2, estimator = "DWLS", cor_method = "poly", rotation = "none")),
    class = "efa_dwls_degenerate_weight")

  # The ULS sandwich consumes the same covariance. It still returns a fit -- the point estimates
  # need no asymptotic covariance -- but the standard errors are withheld rather than reported from
  # a contaminated meat matrix.
  cl <- character()
  sw <- withCallingHandlers(
    suppressMessages(efa_fit(x, n_factors = 2, estimator = "ULS", cor_method = "poly",
                             se = "sandwich", rotation = "none")),
    warning = function(w) { cl <<- c(cl, class(w)); invokeRestart("muffleWarning") })
  expect_s3_class(sw, "efa")
  expect_true(all(is.na(sw$SE$unrot_loadings)))
  expect_true("efa_se_unreliable" %in% cl)

  # The boundary pair is named, and reporting its correlation at the boundary makes this matrix
  # indefinite, so it is smoothed. Both are asserted because both are consequences of the boundary
  # estimate that a change to it would silently alter: at the old (lower) estimate this matrix was
  # positive definite and no smoothing happened.
  expect_true("efa_cor_boundary" %in% cl)
  expect_true("efa_cor_smoothed" %in% cl)

  # The degeneracy screen fires on this data, and does so on every platform: the variance reaches
  # it as NA by construction rather than as a number whose size decided the branch.
  expect_true("efa_acov_degenerate" %in% cl)

  # Clean ordinal data raises no degeneracy warning.
  cl_clean <- .warn_classes(suppressMessages(
    .prepare_cor_input(DOSPERT_raw[stats::complete.cases(DOSPERT_raw), 1:6],
                       cor_method = "poly", acov = "diag", dwls = TRUE,
                       inform_from_data = FALSE)))
  expect_false("efa_acov_degenerate" %in% cl_clean)
})

test_that("a continuity-corrected binary pair still refuses DWLS", {
  # The correction repairs the point estimate of a binary pair with an empty cell, but not its
  # asymptotic variance: the sandwich would treat the nudged counts as if they were the observed
  # data, and simulated coverage of an interval built that way is not trustworthy. So a corrected
  # pair is withheld exactly like a boundary pair, and DWLS refuses the data rather than weighting
  # by a number it cannot justify. This is the test that pins that decision.
  #
  # Six binary items that are monotone in one another, so every pair's 2x2 table has one empty
  # cell and every pair is corrected.
  x <- vapply(1:6, function(k) as.integer((0:239L) %% 7L >= k), integer(240L))
  colnames(x) <- paste0("i", 1:6)

  poly <- suppressWarnings(.polychoric(x, acov = "diag"))
  expect_true(all(poly$zero_corrected))
  expect_true(all(is.na(poly$acov)))

  expect_error(
    suppressWarnings(
      efa_fit(x, n_factors = 1, estimator = "DWLS", cor_method = "tetra", rotation = "none")),
    class = "efa_dwls_degenerate_weight")

  # ULS needs no asymptotic covariance, so the corrected matrix still fits -- which is the point
  # of the correction: these data are analysable now, where the boundary estimate made every
  # correlation 0.9999.
  cl <- character()
  fit <- withCallingHandlers(
    suppressMessages(efa_fit(x, n_factors = 1, estimator = "ULS", cor_method = "tetra",
                             rotation = "none")),
    warning = function(w) { cl <<- c(cl, class(w)); invokeRestart("muffleWarning") })
  expect_s3_class(fit, "efa")
  expect_true("efa_cor_zero_cell" %in% cl)
  expect_false("efa_cor_boundary" %in% cl)
  expect_true(all(abs(fit$orig_R[upper.tri(fit$orig_R)]) < 0.9999))
})

test_that(".warn_acov_degenerate flags an unusable asymptotic variance and passes clean ones", {
  # A variance above 1 already implies a +/- 1 SE interval as wide as the whole [-1, 1] range
  # of a correlation; non-finite and non-positive values are unusable a fortiori. Asserted on
  # supplied vectors rather than on an estimate, so the gate is exercised exactly rather than
  # through whichever side of it a fitted pair happens to land on.
  expect_warning(.warn_acov_degenerate(c(a = 0.01, b = 12, c = 0.02)),
                 class = "efa_acov_degenerate")
  expect_warning(.warn_acov_degenerate(c(a = 0.01, b = Inf, c = 0.02)),
                 class = "efa_acov_degenerate")
  expect_warning(.warn_acov_degenerate(c(a = 0.01, b = NaN, c = 0.02)),
                 class = "efa_acov_degenerate")
  expect_warning(.warn_acov_degenerate(c(a = 0.01, b = 0, c = 0.02)),
                 class = "efa_acov_degenerate")
  expect_warning(.warn_acov_degenerate(c(a = 0.01, b = -1, c = 0.02)),
                 class = "efa_acov_degenerate")

  # Ordinary variances pass, and an unnamed vector is reported positionally rather than failing
  # on the missing labels.
  expect_no_warning(.warn_acov_degenerate(c(a = 0.01, b = 0.5, c = 0.02)))
  expect_warning(.warn_acov_degenerate(c(0.01, 12, 0.02)), class = "efa_acov_degenerate")
})

test_that("both acov requests screen the same pair-labelled asymptotic variances", {
  # .prepare_cor_input() screens the variances once, before either consumer touches them, so a
  # degenerate pair is reported identically whether it reaches the DWLS weights (acov = "diag")
  # or the sandwich meat (acov = "full", whose diagonal is the same quantity). Captured from the
  # screen itself: the value that decides the warning is undetermined for a boundary pair, so
  # what is asserted is the routing -- that both paths screen, and screen the same labelled
  # diagonal -- on clean data where the screen is silent either way.
  seen <- list()
  local_mocked_bindings(
    .warn_acov_degenerate = function(acov_diag, labels = names(acov_diag)) {
      seen[[length(seen) + 1L]] <<- acov_diag
      invisible(NULL)
    }
  )
  x <- DOSPERT_raw[stats::complete.cases(DOSPERT_raw), 1:5]
  labels <- apply(utils::combn(colnames(x), 2L), 2L, paste, collapse = "-")

  suppressMessages(.prepare_cor_input(x, cor_method = "poly", acov = "diag", dwls = TRUE,
                                      inform_from_data = FALSE))
  suppressMessages(.prepare_cor_input(x, cor_method = "poly", acov = "full", dwls = TRUE,
                                      inform_from_data = FALSE))

  expect_length(seen, 2L)
  expect_named(seen[[1L]], labels)
  expect_named(seen[[2L]], labels)
  # The "full" path reaches the screen through diag(Gamma), which accumulates the same per-cell
  # influences by a different route (scatter then crossprod) than the "diag" path's direct sum.
  expect_equal(seen[[1L]], seen[[2L]], tolerance = 1e-6)
})

test_that("cor_method = 'poly' rejects unordered factor / character columns through efa_fit", {
  set.seed(123)
  N <- 400L
  z1 <- stats::rnorm(N)
  z2 <- 0.7 * z1 + sqrt(1 - 0.7^2) * stats::rnorm(N)
  br <- stats::qnorm(c(.33, .66))
  lab <- c("low", "mid", "high")
  df <- data.frame(a = factor(lab[findInterval(z1, br) + 1L]),   # unordered
                   b = factor(lab[findInterval(z2, br) + 1L]))
  expect_error(
    suppressMessages(efa_fit(df, n_factors = 1, cor_method = "poly")),
    class = "efa_cor_unordered_factor")
})

test_that("method = 'DWLS' requires a polychoric asymptotic covariance", {
  # a correlation matrix carries no raw data to estimate an asymptotic covariance from
  expect_error(EFA(test_models$baseline$cormat, n_factors = 3, method = "DWLS", N = 500),
               class = "efa_dwls_no_acov")
  # a continuous correlation method produces no polychoric asymptotic covariance
  expect_error(EFA(DOSPERT_raw, n_factors = 3, method = "DWLS", cor_method = "pearson"),
               class = "efa_dwls_no_acov")
  # the abort fires before fitting, regardless of bootstrap settings
  expect_error(EFA(test_models$baseline$cormat, n_factors = 3, method = "DWLS", N = 500,
                   se = "np-boot"),
               class = "efa_dwls_no_acov")
})

test_that("KMO() and BARTLETT() honour cor_method = 'poly'", {
  x <- DOSPERT_raw[stats::complete.cases(DOSPERT_raw), ]
  N <- nrow(x)
  Rp <- .polychoric(x)$R

  # cor_method = "poly" must produce exactly the result of passing that polychoric
  # matrix in directly: this pins the routing through .prepare_cor_input(). The raw
  # call infers N = nrow(x), the same N supplied with the matrix.
  expect_equal(KMO(x, cor_method = "poly")$KMO, KMO(Rp)$KMO)
  expect_equal(BARTLETT(x, cor_method = "poly")$chisq,
               BARTLETT(Rp, N = N)$chisq)
})

test_that("KMO() and BARTLETT() with cor_method = 'poly' match an external polychoric reference", {
  skip_on_cran()
  skip_if_not_installed("psych")

  x <- DOSPERT_raw[stats::complete.cases(DOSPERT_raw), ]
  N <- nrow(x)
  Rp <- suppressWarnings(psych::polychoric(x, correct = FALSE, global = FALSE)$rho)

  # The two polychoric matrices agree to ~1e-4, so the suitability statistics agree
  # to well within 1e-2.
  expect_equal(KMO(x, cor_method = "poly")$KMO, psych::KMO(Rp)$MSA,
               tolerance = 1e-2)
  expect_equal(BARTLETT(x, cor_method = "poly")$chisq,
               psych::cortest.bartlett(Rp, n = N)$chisq,
               tolerance = 1e-2, ignore_attr = TRUE)
})

test_that("the simulation-based criteria reject poly / tetra with a classed condition", {
  expect_error(CD(DOSPERT_raw, cor_method = "poly"),
               class = "efa_cor_method_unsupported")
  expect_error(PARALLEL(DOSPERT_raw, cor_method = "tetra"),
               class = "efa_cor_method_unsupported")
  expect_error(NEST(DOSPERT_raw, N = 200, cor_method = "poly"),
               class = "efa_cor_method_unsupported")
  # HULL derives its factor-search bound from an internal parallel analysis, so it
  # rejects poly/tetra for the same reason.
  expect_error(HULL(DOSPERT_raw, N = 200, cor_method = "tetra"),
               class = "efa_cor_method_unsupported")
})

test_that("cor_method = 'poly' honours `use` for missing data", {
  g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), ]
  gm <- g
  gm[1:5, 1] <- NA

  # all.obs and everything do not delete missing data, so they error on it,
  # matching the outcome of stats::cor() for the Pearson path.
  expect_error(KMO(gm, cor_method = "poly", use = "all.obs"),
               class = "efa_cor_na")
  expect_error(KMO(gm, cor_method = "poly", use = "everything"),
               class = "efa_cor_na")

  # complete.obs listwise-deletes first, so the matrix equals the polychoric of
  # the complete-case data.
  R_listwise <- .polychoric(gm[stats::complete.cases(gm), ])$R
  expect_equal(KMO(gm, cor_method = "poly", use = "complete.obs")$KMO,
               KMO(R_listwise)$KMO)
})

test_that("cor_method = 'poly' uses pairwise-complete thresholds under missing data", {
  skip_on_cran()
  skip_if_not_installed("polycor")

  set.seed(202)
  N <- 2000
  L <- chol(matrix(c(1, .6, .6, 1), 2L))
  Z <- matrix(stats::rnorm(2 * N), N, 2L) %*% L
  a <- findInterval(Z[, 1], stats::qnorm(c(.3, .6)))
  b <- findInterval(Z[, 2], stats::qnorm(c(.25, .55, .8)))
  # Missing-at-random: a is dropped more often where b is high, so b's pairwise-complete
  # marginal differs from its full-column marginal. A pairwise polychoric must take its
  # thresholds from each pair's own complete cases (as polycor does); using the full-column
  # marginals instead would bias this estimate by ~0.06.
  am <- a; am[(b >= 2) & (stats::runif(N) < 0.7)] <- NA
  ours <- .polychoric(cbind(am, b))$R[1, 2]
  ref <- suppressWarnings(polycor::polychor(am, b, ML = FALSE))
  expect_equal(ours, ref, tolerance = 1e-3)
})

test_that("an asymptotic covariance over pairwise data reports the listwise override", {
  g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), ]
  gm <- g; gm[1:5, 1] <- NA

  # With missing data, an asymptotic covariance (here the DWLS weights) forces listwise
  # deletion even though `use` asks for pairwise-complete estimation; that override is reported.
  expect_message(
    suppressWarnings(
      .prepare_cor_input(gm, cor_method = "poly", acov = "diag",
                         use = "pairwise.complete.obs", inform_from_data = FALSE),
      classes = "efa_cor_sparse_cells"
    ),
    class = "efa_acov_listwise"
  )
  # Complete data needs no override, so nothing is reported.
  expect_no_message(
    suppressWarnings(
      .prepare_cor_input(g, cor_method = "poly", acov = "diag",
                         use = "pairwise.complete.obs", inform_from_data = FALSE),
      classes = "efa_cor_sparse_cells"
    ),
    class = "efa_acov_listwise"
  )
  # A diagonal (DWLS-weight) covariance is an ordinal construct; requesting it for a Pearson
  # correlation is rejected rather than silently producing no weights.
  expect_error(
    .prepare_cor_input(gm, cor_method = "pearson", acov = "diag",
                       use = "pairwise.complete.obs", inform_from_data = FALSE),
    class = "efa_acov_unsupported"
  )
})

test_that("an acov requested for a correlation matrix is ignored with a warning", {
  set.seed(1)
  R <- stats::cor(matrix(stats::rnorm(200L), ncol = 4L))
  expect_warning(
    .prepare_cor_input(R, N = 50, acov = "full", inform_from_data = FALSE),
    class = "efa_acov_ignored"
  )
})

test_that("an unsupported acov/cor_method is rejected before any listwise notice", {
  set.seed(1)
  dat <- matrix(stats::rnorm(150L), ncol = 3L)
  dat[1:3, 1L] <- NA
  # full + a rank correlation is unsupported; with missing data the abort must precede (and
  # replace) the listwise-override notice rather than follow it.
  expect_error(
    .prepare_cor_input(dat, cor_method = "spearman", acov = "full",
                       use = "pairwise.complete.obs", inform_from_data = FALSE),
    class = "efa_acov_unsupported"
  )
})

test_that("N_FACTORS skips reference-based criteria under poly with an informative note", {
  classes <- character()
  withCallingHandlers(
    N_FACTORS(GRiPS_raw, criteria = c("EKC", "HULL"), cor_method = "poly",
              method = "ULS", suitability = FALSE),
    warning = function(w) {
      classes <<- c(classes, class(w))
      invokeRestart("muffleWarning")
    }
  )
  # HULL is skipped with the dedicated 'skipped' class, not reported as a generic
  # failure; EKC still runs on the polychoric matrix.
  expect_true("efa_criterion_skipped" %in% classes)
  expect_false("efa_criterion_failed" %in% classes)
})

test_that("N_FACTORS skips SMT under poly because its normal-theory test is invalid", {
  classes <- character()
  out <- withCallingHandlers(
    N_FACTORS(GRiPS_raw, criteria = c("EKC", "SMT"), cor_method = "poly",
              method = "ULS", suitability = FALSE),
    warning = function(w) {
      classes <<- c(classes, class(w))
      invokeRestart("muffleWarning")
    }
  )
  # SMT's sequential model tests rest on the normal-theory ML chi-square, which is
  # not valid for polychoric/tetrachoric correlations, so it is skipped (not run on
  # an inappropriate matrix); EKC, a descriptive eigenvalue rule, still runs.
  expect_true("efa_criterion_skipped" %in% classes)
  expect_true("SMT" %in% names(out$not_run))
  expect_false("SMT" %in% names(out$outputs))
  expect_true("EKC" %in% names(out$outputs))
})

test_that("the polychoric bootstrap is reproducible and positive definite under resampling", {
  g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), ]

  b1 <- EFA(g, n_factors = 1, method = "ULS", cor_method = "poly",
            se = "np-boot", b_boot = 6, seed = 42)
  b2 <- EFA(g, n_factors = 1, method = "ULS", cor_method = "poly",
            se = "np-boot", b_boot = 6, seed = 42)

  # Same seed -> identical bootstrap SEs, independent of how many replicates were fit.
  expect_identical(b1$boot$SE$unrot_loadings, b2$boot$SE$unrot_loadings)
  # Replicate matrices are recomputed per resample (degenerate ones are dropped at
  # the fit), so the surviving replicates yield finite SEs.
  expect_true(all(is.finite(b1$boot$SE$unrot_loadings)))
})

test_that("the DWLS polychoric bootstrap is reproducible and positive definite under resampling", {
  g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), ]

  b1 <- suppressWarnings(EFA(g, n_factors = 1, method = "DWLS", cor_method = "poly",
                             se = "np-boot", b_boot = 6, seed = 42))
  b2 <- suppressWarnings(EFA(g, n_factors = 1, method = "DWLS", cor_method = "poly",
                             se = "np-boot", b_boot = 6, seed = 42))

  # Same seed -> identical bootstrap SEs; each replicate recomputes its own polychoric
  # matrix and diagonal-ACOV weights, and the surviving replicates yield finite SEs.
  expect_identical(b1$boot$SE$unrot_loadings, b2$boot$SE$unrot_loadings)
  expect_true(all(is.finite(b1$boot$SE$unrot_loadings)))
})

test_that("the polychoric bootstrap on DOSPERT_raw completes (timed)", {
  skip_on_cran()
  skip_if_not_slow()

  b_boot <- 200L
  t <- system.time(
    fit <- EFA(DOSPERT_raw, n_factors = 3, method = "ULS", cor_method = "poly",
               se = "np-boot", b_boot = b_boot, seed = 1)
  )
  elapsed <- unname(t["elapsed"])
  cli::cli_inform(
    "Polychoric bootstrap: {b_boot} replicates on DOSPERT_raw in {round(elapsed, 1)}s ({round(elapsed / b_boot, 3)}s/replicate).")
  expect_true(all(is.finite(fit$boot$SE$unrot_loadings)))
})
