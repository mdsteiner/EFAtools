# Robust (Godambe sandwich) standard errors and the scaled (Satorra-Bentler) chi-square for the
# ordinal/polychoric path (method ML/ULS/DWLS, cor_method poly/tetra).

test_that("sandwich SEs fill the unrotated SE/CI schema and a scaled chi-square", {
  dat <- DOSPERT_raw[, 1:8]
  fit <- EFA(dat, n_factors = 2, cor_method = "poly", method = "DWLS",
             rotation = "none", se = "sandwich")

  # Unrotated loading and uniqueness SEs are present, finite, and positive.
  expect_equal(dim(fit$SE$unrot_loadings), c(8L, 2L))
  expect_true(all(is.finite(fit$SE$unrot_loadings)))
  expect_true(all(fit$SE$unrot_loadings > 0))
  expect_length(fit$SE$uniquenesses, 8L)
  expect_true(all(is.finite(fit$SE$uniquenesses)))

  # Wald intervals bracket the point estimates.
  ci <- fit$CI$unrot_loadings
  expect_true(all(ci$lower <= unclass(fit$unrot_loadings)))
  expect_true(all(ci$upper >= unclass(fit$unrot_loadings)))

  # The scaled chi-square block is populated (it is NA from .gof() for DWLS).
  fi <- fit$fit_indices
  expect_true(is.finite(fi$chi))
  expect_true(is.finite(fi$p_chi))
  expect_true(all(is.finite(c(fi$CFI, fi$TLI, fi$RMSEA))))
  expect_identical(fi$chi_scaled_type, "scaled.shifted")
  expect_true(all(is.finite(c(fi$chi_scaling, fi$chi_shift, fi$chi_unscaled,
                              fi$chi_mean_adjusted, fi$chi_mean_var, fi$df_mean_var))))
  # AIC/BIC/ECVI are likelihood-ratio chi-square criteria with no meaning for the moment-scaled
  # statistic, so they are deliberately left NA on the sandwich path.
  expect_true(all(is.na(c(fi$AIC, fi$BIC, fi$ECVI))))
})


test_that("the unrotated sandwich CI-provenance note names the robust sandwich, not the information matrix", {
  skip_on_cran()

  dat <- DOSPERT_raw[, 1:8]
  fit <- EFA(dat, n_factors = 2, cor_method = "poly", method = "DWLS",
             rotation = "none", se = "sandwich")

  # The full summary carries the CI-provenance note; for a robust (Godambe sandwich) fit it must
  # describe the sandwich covariance, not the expected information matrix (which is correct only
  # for se = "information").
  body <- cli::ansi_strip(format(summary(fit)))
  note <- body[grepl("Note:", body)]
  expect_true(any(grepl("sandwich", note)))
  expect_false(any(grepl("expected information matrix", note)))
})


test_that(".chi_fit_indices is robust to an undefined (NA) baseline chi-square", {
  # A degenerate scaled baseline yields chi_null = NA while the model chi stays finite; the shared
  # fit-index helper must return NA CFI/TLI (not error) and still compute the model-only RMSEA. The
  # scaled path passes its statistics as the noncentrality inputs (chi_cfi / chi_null_cfi).
  idx <- EFAtools:::.chi_fit_indices(chi = 50, df = 10, chi_null = NA_real_,
                                     df_null = 28, N = 500, m = 8, ci = TRUE,
                                     chi_cfi = 50, chi_null_cfi = NA_real_)
  expect_true(is.na(idx$CFI))
  expect_true(is.na(idx$TLI))
  expect_true(is.na(idx$p_null))
  expect_true(is.finite(idx$RMSEA))
  expect_true(is.finite(idx$p_chi))
})


test_that("the unscaled sandwich statistic equals the DWLS objective Fm", {
  dat <- DOSPERT_raw[, 1:8]
  fit <- EFA(dat, n_factors = 2, cor_method = "poly", method = "DWLS",
             rotation = "none", se = "sandwich")
  # T = N (s - sigma)' V (s - sigma) with V = 1/diag(Gamma) equals the weighted off-diagonal
  # objective the DWLS backend minimised, because the variance-scale weights absorb the N.
  expect_equal(fit$fit_indices$chi_unscaled, fit$fit_indices$Fm, tolerance = 1e-8)
})


test_that("requesting the sandwich does not change the DWLS point estimate", {
  dat <- DOSPERT_raw[, 1:6]
  base <- EFA(dat, n_factors = 1, cor_method = "poly", method = "DWLS",
              rotation = "none", se = "none")
  rob <- EFA(dat, n_factors = 1, cor_method = "poly", method = "DWLS",
             rotation = "none", se = "sandwich")
  # Both fit on the listwise-complete rows, and the full ACOV's diagonal reproduces the DWLS
  # "diag" weights, so the converged loadings are identical (up to the optimiser tolerance).
  expect_equal(unclass(rob$unrot_loadings), unclass(base$unrot_loadings), tolerance = 1e-6)
})


test_that("sandwich SEs propagate through an oblique rotation", {
  dat <- DOSPERT_raw[, 1:8]
  fit <- EFA(dat, n_factors = 2, cor_method = "poly", method = "DWLS",
             rotation = "oblimin", se = "sandwich")

  expect_equal(dim(fit$SE$rot_loadings), c(8L, 2L))
  expect_true(all(is.finite(fit$SE$rot_loadings)))
  # Factor correlations: a symmetric SE matrix with a fixed (zero-variance) unit diagonal.
  expect_equal(dim(fit$SE$Phi), c(2L, 2L))
  expect_equal(diag(fit$SE$Phi), c(0, 0))
  expect_equal(fit$SE$Phi, t(fit$SE$Phi))
  expect_true(is.finite(fit$SE$Phi[1, 2]))
  # Structure coefficients and communalities are reported for the oblique solution.
  expect_equal(dim(fit$SE$Structure), c(8L, 2L))
  expect_length(fit$SE$communalities, 8L)
  # Communalities are rotation-invariant, so their SEs equal the unrotated uniqueness SEs.
  expect_equal(unname(fit$SE$communalities), unname(fit$SE$uniquenesses),
               tolerance = 1e-8)
})


test_that("sandwich SEs are available for ULS and ML on the ordinal path", {
  dat <- DOSPERT_raw[, 1:8]
  for (m in c("ULS", "ML")) {
    fit <- EFA(dat, n_factors = 2, cor_method = "poly", method = m,
               rotation = "none", se = "sandwich")
    expect_true(all(is.finite(fit$SE$unrot_loadings)),
                info = paste("method", m))
    expect_true(is.finite(fit$fit_indices$chi), info = paste("method", m))
  }
})


test_that("sandwich SEs reject the unsupported method/rotation/correlation combinations", {
  dat <- DOSPERT_raw[, 1:6]
  cormat <- test_models$baseline$cormat

  # Continuous correlations have no polychoric ACOV for the robust meat.
  expect_error(
    EFA(dat, n_factors = 2, cor_method = "pearson", method = "DWLS", se = "sandwich"),
    class = "efa_se_unsupported"
  )
  # PAF has no discrepancy-based information for a sandwich.
  expect_error(
    EFA(dat, n_factors = 2, cor_method = "poly", method = "PAF", se = "sandwich"),
    class = "efa_se_unsupported"
  )
  # promax and simplimax have no usable analytic rotation Jacobian.
  expect_error(
    EFA(dat, n_factors = 2, cor_method = "poly", method = "DWLS",
        rotation = "promax", se = "sandwich"),
    class = "efa_se_unsupported"
  )
  expect_error(
    EFA(dat, n_factors = 2, cor_method = "poly", method = "DWLS",
        rotation = "simplimax", se = "sandwich"),
    class = "efa_se_unsupported"
  )
  # A correlation matrix carries no raw data to estimate the polychoric ACOV from. (DWLS rejects a
  # correlation matrix even earlier, via its own weight precondition, so use ML here to exercise
  # the sandwich-specific raw-data guard.)
  expect_error(
    EFA(cormat, n_factors = 3, N = 500, cor_method = "poly", method = "ML", se = "sandwich"),
    class = "efa_se_unsupported"
  )
})


test_that("single-factor sandwich loading SEs match lavaan robust.sem (DWLS and ULS)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  dat <- as.data.frame(DOSPERT_raw[, 1:6])
  vn <- colnames(dat)
  mod <- paste0("f =~ ", paste(vn, collapse = " + "))

  for (pair in list(c("DWLS", "DWLS"), c("ULS", "ULS"))) {
    fit <- EFA(dat, n_factors = 1, cor_method = "poly", method = pair[1],
               rotation = "none", se = "sandwich")

    lfit <- lavaan::cfa(mod, data = dat, ordered = vn, estimator = pair[2],
                        se = "robust.sem", test = "scaled.shifted", std.lv = TRUE)
    pe <- lavaan::parameterEstimates(lfit, standardized = TRUE)
    lam <- pe[pe$op == "=~", ]
    L_lav <- lam$std.all; SE_lav <- lam$se

    L_efa <- as.numeric(fit$unrot_loadings)
    sgn <- sign(sum(L_efa * L_lav))

    # One factor: no rotational indeterminacy, so the loadings (up to sign) and their robust SEs
    # are directly comparable. The residual gap is the small polychoric-correlation difference.
    expect_equal(as.numeric(fit$unrot_loadings) * sgn, L_lav, tolerance = 0.01,
                 info = pair[1])
    expect_equal(as.numeric(fit$SE$unrot_loadings), SE_lav, tolerance = 0.01,
                 info = pair[1])
  }
})


test_that("the scaled chi-square matches lavaan WLSMV", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  dat <- as.data.frame(DOSPERT_raw[, 1:8])
  vn <- colnames(dat)

  fit <- EFA(dat, n_factors = 2, cor_method = "poly", method = "DWLS",
             rotation = "none", se = "sandwich")

  lf <- lavaan::efa(dat, nfactors = 2, ordered = vn, estimator = "DWLS",
                    rotation = "geomin", test = "scaled.shifted", se = "robust.sem")
  lfit <- lf[[1]]
  ss <- lavaan::lavInspect(lfit, "test")[["scaled.shifted"]]

  # The scaled-shifted statistic, its shift, and the degrees of freedom match lavaan WLSMV; the
  # chi-square gap is the small polychoric-correlation difference (the shift, a function of the
  # scale-invariant trace coefficients, matches more tightly).
  expect_equal(fit$fit_indices$df, ss$df)
  expect_equal(fit$fit_indices$chi, ss$stat, tolerance = 0.02)
  expect_equal(fit$fit_indices$chi_shift, ss$shift.parameter, tolerance = 1e-3)
  # lavaan stores the scaling factor as the reciprocal of the multiplier a.
  expect_equal(1 / fit$fit_indices$chi_scaling, ss$scaling.factor, tolerance = 1e-3)
})


test_that("the EFAtools polychoric Gamma matches lavaan's NACOV up to the N scale", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  dat <- as.data.frame(DOSPERT_raw[, 1:6])
  vn <- colnames(dat)

  poly <- EFAtools:::.polychoric(as.matrix(dat), nearest_pd = FALSE,
                                 binary_only = FALSE, acov = "full")
  N <- sum(stats::complete.cases(dat))

  lf <- lavaan::efa(dat, nfactors = 1, ordered = vn, estimator = "DWLS", se = "robust.sem")
  lfit <- lf[[1]]
  Glav <- lavaan::lavInspect(lfit, "gamma")
  if (is.list(Glav)) Glav <- Glav[[1]]
  cor_rows <- grep("~~", rownames(Glav))
  Glav_corr <- unname(Glav[cor_rows, cor_rows])

  # EFAtools' Gamma is on the variance scale Var(rho-hat); lavaan's NACOV is N times that.
  rel_F <- norm(N * poly$acov - Glav_corr, "F") / norm(Glav_corr, "F")
  expect_lt(rel_F, 1e-3)
})


# Continuous (Pearson) path: the fourth-moment (ADF; Browne, 1984) covariance of the sample
# correlations is the robust meat for method ML/ULS, giving MLM/MLR-style robust SEs and a scaled
# (Satorra-Bentler) chi-square without an ordinal correlation method.

# Deterministic non-normal continuous fixture: a k-factor model whose common factors and unique
# errors are standardised chi-square(1) draws (mean 0, variance 1, skew ~ 2.8, excess kurtosis
# ~ 12), so the columns are markedly non-normal and the ADF covariance departs clearly from
# normal theory. Each column has approximately unit variance and the population correlations are
# the simple-structure loading products.
.adf_fixture <- function(N, loadings, seed) {
  set.seed(seed)
  p <- nrow(loadings)
  k <- ncol(loadings)
  skew <- function(n) (stats::rchisq(n, df = 1) - 1) / sqrt(2)
  Fmat <- matrix(skew(N * k), N, k)
  E <- matrix(skew(N * p), N, p)
  uniq <- sqrt(pmax(1 - rowSums(loadings^2), 0))
  X <- Fmat %*% t(loadings) + sweep(E, 2L, uniq, "*")
  colnames(X) <- paste0("V", seq_len(p))
  as.data.frame(X)
}


test_that("continuous Pearson sandwich SEs fill the SE/CI schema and a scaled chi-square", {
  dat <- .adf_fixture(N = 600, loadings = matrix(c(.7, .65, .6, .75, .55, .5, .45, .6), ncol = 1L),
                      seed = 1)

  for (m in c("ML", "ULS")) {
    fit <- EFA(dat, n_factors = 1, cor_method = "pearson", method = m,
               rotation = "none", se = "sandwich")

    # Unrotated loading and uniqueness SEs are present, finite, and positive.
    expect_true(all(is.finite(fit$SE$unrot_loadings)), info = m)
    expect_true(all(fit$SE$unrot_loadings > 0), info = m)
    expect_true(all(is.finite(fit$SE$uniquenesses)), info = m)

    # Wald intervals bracket the point estimates.
    ci <- fit$CI$unrot_loadings
    expect_true(all(ci$lower <= unclass(fit$unrot_loadings)), info = m)
    expect_true(all(ci$upper >= unclass(fit$unrot_loadings)), info = m)

    # The scaled chi-square block is populated (the unscaled ML/ULS discrepancy is not robust to
    # non-normality).
    fi <- fit$fit_indices
    expect_true(is.finite(fi$chi), info = m)
    expect_true(is.finite(fi$p_chi), info = m)
    expect_true(all(is.finite(c(fi$CFI, fi$TLI, fi$RMSEA))), info = m)
    expect_identical(fi$chi_scaled_type, "scaled.shifted")
    expect_true(all(is.na(c(fi$AIC, fi$BIC, fi$ECVI))), info = m)
  }
})


test_that("continuous Pearson sandwich SEs propagate through an oblique rotation", {
  L <- matrix(0, 8L, 2L)
  L[1:4, 1] <- c(.7, .65, .6, .55)
  L[5:8, 2] <- c(.7, .6, .65, .5)
  dat <- .adf_fixture(N = 700, loadings = L, seed = 2)

  fit <- EFA(dat, n_factors = 2, cor_method = "pearson", method = "ML",
             rotation = "oblimin", se = "sandwich")

  expect_equal(dim(fit$SE$rot_loadings), c(8L, 2L))
  expect_true(all(is.finite(fit$SE$rot_loadings)))
  expect_equal(dim(fit$SE$Phi), c(2L, 2L))
  expect_equal(diag(fit$SE$Phi), c(0, 0))
  expect_equal(fit$SE$Phi, t(fit$SE$Phi))
  expect_true(is.finite(fit$SE$Phi[1, 2]))
  expect_equal(dim(fit$SE$Structure), c(8L, 2L))
  expect_length(fit$SE$communalities, 8L)
  # Communalities are rotation-invariant, so their SEs equal the unrotated uniqueness SEs.
  expect_equal(unname(fit$SE$communalities), unname(fit$SE$uniquenesses),
               tolerance = 1e-8)
})


test_that("continuous Pearson sandwich rejects the unsupported method/rotation/data combinations", {
  dat <- .adf_fixture(N = 400, loadings = matrix(c(.7, .6, .65, .55, .5, .6), ncol = 1L), seed = 3)
  cormat <- test_models$baseline$cormat

  # DWLS has no continuous asymptotic covariance (it is an ordinal estimator).
  expect_error(
    EFA(dat, n_factors = 1, cor_method = "pearson", method = "DWLS", se = "sandwich"),
    class = "efa_se_unsupported"
  )
  # promax and simplimax have no usable analytic rotation Jacobian.
  expect_error(
    EFA(dat, n_factors = 2, cor_method = "pearson", method = "ML",
        rotation = "promax", se = "sandwich"),
    class = "efa_se_unsupported"
  )
  expect_error(
    EFA(dat, n_factors = 2, cor_method = "pearson", method = "ML",
        rotation = "simplimax", se = "sandwich"),
    class = "efa_se_unsupported"
  )
  # Spearman/Kendall correlations have no fourth-moment covariance in this implementation.
  expect_error(
    EFA(dat, n_factors = 1, cor_method = "spearman", method = "ML", se = "sandwich"),
    class = "efa_se_unsupported"
  )
  # A correlation matrix carries no raw data to estimate the ADF covariance from.
  expect_error(
    EFA(cormat, n_factors = 3, N = 500, cor_method = "pearson", method = "ML", se = "sandwich"),
    class = "efa_se_unsupported"
  )
})


test_that(".prepare_cor_input rejects a full ADF covariance for non-Pearson correlations", {
  dat <- .adf_fixture(N = 200, loadings = matrix(c(.7, .6, .5, .55), ncol = 1L), seed = 7)
  # The fourth-moment ADF covariance is Pearson-specific; a rank correlation with acov = "full"
  # would mismatch R and Gamma. EFA() gates this earlier, but the helper guards its own contract.
  expect_error(
    EFAtools:::.prepare_cor_input(as.matrix(dat), cor_method = "spearman", acov = "full"),
    class = "efa_acov_unsupported"
  )
})


test_that("the continuous ADF Gamma matches lavaan's correlation NACOV up to the N scale", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if(utils::packageVersion("lavaan") < "0.6.9")

  L <- matrix(c(.7, .65, .6, .75, .55, .5), ncol = 1L)
  dat <- .adf_fixture(N = 600, loadings = L, seed = 4)
  vn <- colnames(dat)
  N <- nrow(dat)

  Gamma <- EFAtools:::.adf_gamma(as.matrix(dat))

  # correlation = TRUE reparameterises lavaan onto the correlation structure: its NACOV omits the
  # variance (Vi~~Vi) rows and orders the off-diagonal ~~ rows exactly as utils::combn(p, 2).
  mod <- paste0("f =~ ", paste(vn, collapse = " + "))
  lfit <- lavaan::cfa(mod, data = dat, estimator = "MLM", std.lv = TRUE,
                      correlation = TRUE, meanstructure = TRUE)
  Glav <- lavaan::lavInspect(lfit, "gamma")
  if (is.list(Glav)) Glav <- Glav[[1]]
  # correlation = TRUE drops the variance (Vi~~Vi) rows, leaving exactly the p(p-1)/2 off-diagonal
  # correlation rows in utils::combn(p, 2) order. `correlation` is consumed via ... (not a formal
  # cfa() argument), so a lavaan build that silently ignored it keeps the variance rows (a
  # covariance-structure NACOV that is not comparable) -- detect that and skip rather than fail on
  # a feature-version guard that does not track when the option landed.
  sides <- strsplit(rownames(Glav), "~~")
  is_var <- vapply(sides, function(s) length(s) == 2L && s[1] == s[2], logical(1))
  skip_if(any(is_var), "this lavaan build did not apply correlation = TRUE")
  cor_rows <- which(vapply(sides, function(s) length(s) == 2L && s[1] != s[2], logical(1)))
  Glav_corr <- unname(Glav[cor_rows, cor_rows])

  # EFAtools' Gamma is on the variance scale Var(rho-hat); lavaan's NACOV is N times that.
  rel_F <- norm(N * unname(Gamma) - Glav_corr, "F") / norm(Glav_corr, "F")
  expect_lt(rel_F, 1e-4)
})


test_that("single-factor continuous sandwich loading SEs and scaled chi match lavaan robust", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if(utils::packageVersion("lavaan") < "0.6.9")

  L <- matrix(c(.7, .65, .6, .75, .55, .5), ncol = 1L)
  dat <- .adf_fixture(N = 700, loadings = L, seed = 5)
  vn <- colnames(dat)
  mod <- paste0("f =~ ", paste(vn, collapse = " + "))

  # Compared against ML, the unambiguous oracle: the MLMV shortcut (= ML + scaled.shifted test +
  # robust.sem SEs) routes through lavaan's working correlation-structure path, whereas the manual
  # se = "robust.sem" override errors there. ULS has no such lavaan oracle (its continuous robust
  # correlation path is unsupported), so it is cross-checked against the bootstrap below instead.
  fit <- EFA(dat, n_factors = 1, cor_method = "pearson", method = "ML",
             rotation = "none", se = "sandwich")

  lfit <- lavaan::cfa(mod, data = dat, estimator = "MLMV", std.lv = TRUE,
                      correlation = TRUE, meanstructure = TRUE)
  pe <- lavaan::parameterEstimates(lfit)
  lam <- pe[pe$op == "=~", ]
  L_lav <- lam$est; SE_lav <- lam$se

  L_efa <- as.numeric(fit$unrot_loadings)
  sgn <- sign(sum(L_efa * L_lav))

  # One factor: no rotational indeterminacy, and the Pearson correlations are identical to
  # lavaan's, so the loadings (up to sign) and their robust SEs match closely. The SEs differ only
  # by the N vs N-1 convention (a uniform ~sqrt(N/(N-1)) ratio), hence the tight 2e-3 band.
  expect_equal(L_efa * sgn, L_lav, tolerance = 0.01)
  expect_equal(as.numeric(fit$SE$unrot_loadings), SE_lav, tolerance = 2e-3)

  ss <- lavaan::lavInspect(lfit, "test")[["scaled.shifted"]]
  expect_equal(fit$fit_indices$df, ss$df)
  # The chi-square carries a small (~1%) gap: EFAtools' statistic is the quadratic-form
  # approximation N (s - sigma)' V (s - sigma), lavaan's is the exact ML likelihood ratio.
  expect_equal(fit$fit_indices$chi, ss$stat, tolerance = 0.02)
  # The shift is invariant to the (arbitrary) scale of the ML weight matrix, so it matches lavaan
  # tightly; the scaling factor a = sqrt(df / c2) is not (EFAtools uses the 1/2 normal-theory ML
  # weight, lavaan a different normalisation), but it cancels in the statistic.
  expect_equal(fit$fit_indices$chi_shift, ss$shift.parameter, tolerance = 1e-3)
})


test_that("continuous sandwich SEs are in the same ballpark as the bootstrap (ML and ULS)", {
  skip_on_cran()

  L <- matrix(c(.7, .65, .6, .75, .55, .5), ncol = 1L)
  dat <- .adf_fixture(N = 500, loadings = L, seed = 6)

  for (m in c("ML", "ULS")) {
    rob <- EFA(dat, n_factors = 1, cor_method = "pearson", method = m,
               rotation = "none", se = "sandwich")
    boot <- EFA(dat, n_factors = 1, cor_method = "pearson", method = m,
                rotation = "none", se = "np-boot", b_boot = 50, seed = 123)

    se_rob <- as.numeric(rob$SE$unrot_loadings)
    se_boot <- as.numeric(boot$SE$unrot_loadings)

    # Robust and bootstrap SEs estimate the same sampling variability (observed median ratio ~1.02
    # for both ML and ULS). This is a coarse cross-check robust to the b_boot = 50 Monte-Carlo
    # noise -- the precise validation is the lavaan oracle above (ML); the band still catches a
    # gross (>~40%) systematic miscalibration of the robust SE, which is the only numeric guard for
    # ULS (lavaan's continuous ULS robust correlation path is unsupported).
    expect_true(all(se_rob > 0 & se_boot > 0), info = m)
    expect_gt(stats::median(se_rob / se_boot), 0.7)
    expect_lt(stats::median(se_rob / se_boot), 1.4)
  }
})


# Unrotated robust loading SEs are gauge-dependent: they must be scaled in the orientation the
# estimator reports its loadings in (Lambda'Lambda diagonal for the eigen-based ULS/DWLS, Lambda'
# Psi^-1 Lambda diagonal for ML; Lawley & Maxwell). Independent reference: the Godambe sandwich with
# the gauge-constraint Jacobian taken by CENTRAL FINITE DIFFERENCES of the orientation function,
# sharing no code with the package's analytic constraint (.se_sandwich_constraint).
.ref_sandwich_loading_se <- function(L, Gamma, N, method, gauge) {
  p <- nrow(L); k <- ncol(L); pk <- p * k; G <- N * Gamma
  pr <- utils::combn(p, 2L); pii <- pr[1, ]; pjj <- pr[2, ]; n <- ncol(pr)
  Delta <- matrix(0, n, pk)
  for (f in seq_len(k)) {
    Delta[cbind(seq_len(n), (f - 1L) * p + pii)] <- L[pjj, f]
    Delta[cbind(seq_len(n), (f - 1L) * p + pjj)] <- L[pii, f]
  }
  if (method == "ULS") {
    VD <- Delta
  } else {
    Sig <- tcrossprod(L); diag(Sig) <- 1; P <- solve(Sig)
    Vmat <- matrix(0, n, n)
    for (s in seq_len(n)) {
      a <- pii[s]; b <- pjj[s]
      Vmat[, s] <- 0.5 * (P[pii, a] * P[pjj, b] + P[pii, b] * P[pjj, a])
    }
    VD <- Vmat %*% Delta
  }
  A <- crossprod(Delta, VD)
  gfun <- function(par) {
    Lm <- matrix(par, p, k)
    M <- if (gauge == "LtL") crossprod(Lm) else {
      ps <- 1 - rowSums(Lm^2); crossprod(Lm, Lm / ps)
    }
    M[upper.tri(M)]
  }
  nc <- k * (k - 1L) / 2L; Cmat <- matrix(0, nc, pk); h <- 1e-6; par0 <- as.vector(L)
  for (j in seq_len(pk)) {
    pp <- par0; pp[j] <- pp[j] + h; pm <- par0; pm[j] <- pm[j] - h
    Cmat[, j] <- (gfun(pp) - gfun(pm)) / (2 * h)
  }
  Aug <- rbind(cbind(A, t(Cmat)), cbind(Cmat, matrix(0, nc, nc)))
  Ab <- solve(Aug)[seq_len(pk), seq_len(pk)]
  V_AA <- (Ab %*% crossprod(VD, G %*% VD) %*% Ab) / (N - 1)
  matrix(sqrt(diag(V_AA)), p, k)
}


test_that("the unrotated sandwich loading SEs match the estimator's reporting gauge", {
  # A two-factor structure with cross-loadings, so the Lambda'Lambda and Lambda'Psi^-1 Lambda
  # orientations differ noticeably (and a wrong gauge would be caught).
  set.seed(2024)
  Lp <- matrix(c(.75, .70, .65, .20, .10, .15,
                 .15, .10, .20, .75, .70, .65), 6, 2)
  Sig <- Lp %*% t(Lp); diag(Sig) <- 1
  X <- matrix(stats::rnorm(800 * 6), 800) %*% chol(Sig)
  colnames(X) <- paste0("v", seq_len(6))
  N <- nrow(X)

  # ML reports loadings in the Lambda'Psi^-1 Lambda orientation, so its robust unrotated loading
  # SEs must be scaled in that gauge -- and are materially different from the Lambda'Lambda gauge.
  fit_ml <- EFA(X, n_factors = 2, cor_method = "pearson", method = "ML",
                rotation = "none", se = "sandwich")
  L_ml <- unclass(fit_ml$unrot_loadings)
  ref_ml_correct <- .ref_sandwich_loading_se(L_ml, fit_ml$Gamma, N, "ML", "LtPiL")
  ref_ml_wrong   <- .ref_sandwich_loading_se(L_ml, fit_ml$Gamma, N, "ML", "LtL")
  expect_equal(unclass(fit_ml$SE$unrot_loadings), ref_ml_correct,
               tolerance = 1e-5, ignore_attr = TRUE)
  expect_gt(max(abs(ref_ml_correct - ref_ml_wrong)) / mean(ref_ml_correct), 0.05)

  # ULS reports loadings in the (eigen-based) Lambda'Lambda orientation, so its robust unrotated
  # loading SEs stay in that gauge.
  fit_uls <- EFA(X, n_factors = 2, cor_method = "pearson", method = "ULS",
                 rotation = "none", se = "sandwich")
  L_uls <- unclass(fit_uls$unrot_loadings)
  ref_uls_correct <- .ref_sandwich_loading_se(L_uls, fit_uls$Gamma, N, "ULS", "LtL")
  expect_equal(unclass(fit_uls$SE$unrot_loadings), ref_uls_correct,
               tolerance = 1e-5, ignore_attr = TRUE)
})


test_that(".is_psd accepts PSD covariances and rejects non-PSD ones with a non-negative diagonal", {
  expect_true(EFAtools:::.is_psd(diag(3)))
  expect_true(EFAtools:::.is_psd(matrix(c(2, 1, 1, 2), 2)))
  # Positive diagonal with a round-off-level negative eigenvalue (eigenvalues ~ 2 and -1e-9): inside
  # the -1e-8 tolerance, so accepted. Exercises the eigenvalue gate, not just the diagonal.
  expect_true(EFAtools:::.is_psd(matrix(c(1, 1 + 1e-9, 1 + 1e-9, 1), 2)))
  # All diagonal entries are positive, but the matrix is indefinite (eigenvalues -1 and 3): the
  # diagonal-only gate would wrongly accept it.
  expect_false(EFAtools:::.is_psd(matrix(c(1, 2, 2, 1), 2)))
  expect_false(EFAtools:::.is_psd(matrix(c(1, NA, NA, 1), 2)))
  expect_false(EFAtools:::.is_psd(matrix(c(-1, 0, 0, 1), 2)))
})


test_that("the sandwich core degrades gracefully when the gauge is undefined or undetermined", {
  # Regression: an NA loading or a Heywood uniqueness (psi <= 0) makes the Lambda'Psi^-1 Lambda
  # orientation undefined. Gauge detection must rule it out and fall back to the Lambda'Lambda gauge
  # without erroring on the comparison, and must not abort the (gauge-invariant) rest of the
  # computation.
  G <- diag(3L)   # p = 3 -> p (p - 1) / 2 = 3 off-diagonal pairs

  # An NA loading propagates to psi and to the off-diagonal comparison; detection must coerce the
  # undefined gauge to the Lambda'Lambda fallback rather than raise "missing value where TRUE/FALSE
  # needed".
  fo_na <- list(unrot_loadings = matrix(c(0.7, NA, 0.6, 0.5, 0.4, 0.3), 3, 2),
                orig_R = diag(3), fit_indices = list(df = 1))
  expect_no_error(out_na <- EFAtools:::.se_sandwich_core(fo_na, N = 200, Gamma = G, method = "ML"))
  expect_false(out_na$reliable)
  expect_true(all(is.na(out_na$loadings_se)))

  # A communality above one drives psi[1] < 0; the core must not error.
  fo_hey <- list(unrot_loadings = matrix(c(0.95, 0.6, 0.5, 0.4, 0.3, 0.2), 3, 2),
                 orig_R = diag(3), fit_indices = list(df = 1))
  expect_no_error(EFAtools:::.se_sandwich_core(fo_hey, N = 200, Gamma = G, method = "ML"))

  # Homogeneous uniquenesses (Psi proportional to I) make BOTH orientations diagonal, so the gauge is
  # undetermined by the loadings and the detection takes the method tie-break branch. Build such a
  # solution -- an equal-norm simple structure rotated to introduce cross-loadings -- and confirm the
  # tie-break runs for either estimator without error (this construction is numerically degenerate at
  # the bordering step, so both return reliable = FALSE rather than a value).
  L0 <- matrix(0, 6L, 2L); L0[1:3, 1] <- 0.6; L0[4:6, 2] <- 0.6
  Q <- matrix(c(cos(pi / 6), -sin(pi / 6), sin(pi / 6), cos(pi / 6)), 2L)
  Sig_tie <- tcrossprod(L0 %*% Q); diag(Sig_tie) <- 1
  fo_tie <- list(unrot_loadings = L0 %*% Q, orig_R = Sig_tie, fit_indices = list(df = 4))
  G_tie <- diag(15L)   # p = 6 -> 15 off-diagonal pairs
  expect_no_error(EFAtools:::.se_sandwich_core(fo_tie, N = 300, Gamma = G_tie, method = "ML"))
  expect_no_error(EFAtools:::.se_sandwich_core(fo_tie, N = 300, Gamma = G_tie, method = "ULS"))
})
