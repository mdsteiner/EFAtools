# cor_method = "fiml": two-stage / full-information ML correlations wired into EFA().
# The EM moments engine itself is covered in test-fiml-moments.R; these tests cover the
# routing through .prepare_cor_input() and EFA(), the point estimates (which match
# psych::corFiml() -> psych::fa() and lavaan(missing = "two.stage"), not lavaan::efa()),
# and the classed conditions guarding the unsupported input/option combinations.

# A clean two-factor population covariance: items 1..p/2 load on factor 1, the rest on
# factor 2, with unit diagonal.
fiml_pop_cov <- function(p = 6) {
  L <- matrix(0, p, 2)
  half <- p / 2
  L[seq_len(half), 1] <- 0.7
  L[(half + 1):p, 2] <- 0.7
  S <- tcrossprod(L)
  diag(S) <- 1
  S
}

# A missing-at-random fixture: column 1 is fully observed and drives the missingness in the
# others, so the mechanism depends only on observed data.
fiml_mar_data <- function(n = 800, seed = 456) {
  set.seed(seed)
  X <- MASS::mvrnorm(n, mu = rep(0, 6), Sigma = fiml_pop_cov(6))
  colnames(X) <- paste0("V", seq_len(6))
  X[X[, 1] >  0.8, 2] <- NA
  X[X[, 1] < -0.8, 3] <- NA
  X[X[, 1] >  1.2, 4] <- NA
  X[X[, 1] < -1.2, 5] <- NA
  X
}

# A single-factor MAR fixture: a clean SE / chi-square oracle, with no rotational indeterminacy.
fiml_mar_data_1f <- function(n = 1000, p = 6, load = 0.65, seed = 11) {
  set.seed(seed)
  L <- matrix(load, p, 1)
  S <- tcrossprod(L); diag(S) <- 1
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = S)
  colnames(X) <- paste0("V", seq_len(p))
  X[X[, 1] >  0.8, 2] <- NA
  X[X[, 1] < -0.8, 3] <- NA
  X[X[, 1] >  1.0, 4] <- NA
  X[X[, 1] < -1.0, 5] <- NA
  X
}

test_that("complete data: FIML loadings equal the Pearson-correlation loadings", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  set.seed(123)
  X <- MASS::mvrnorm(500, mu = rep(0, 6), Sigma = fiml_pop_cov(6))
  colnames(X) <- paste0("V", seq_len(6))

  fiml <- EFA(X, n_factors = 2, method = "ML", cor_method = "fiml")
  pear <- EFA(X, n_factors = 2, method = "ML", cor_method = "pearson")

  # With no missing data the FIML moments reduce to the sample correlation, so the
  # extracted loadings coincide.
  expect_equal(unclass(fiml$unrot_loadings), unclass(pear$unrot_loadings),
               tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("MAR-missing data: FIML loadings match psych::fa(corFiml)", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  skip_if_not_installed("psych")

  X <- fiml_mar_data()
  k <- 2L

  efa <- EFA(X, n_factors = k, method = "ML", cor_method = "fiml")

  # Reference: psych's FIML correlation analysed by psych's ML factor analysis, unrotated.
  R_ref  <- suppressWarnings(psych::corFiml(X))
  fa_ref <- suppressWarnings(psych::fa(R_ref, nfactors = k, fm = "ml",
                                       rotate = "none", n.obs = nrow(X)))
  L_ref <- unclass(fa_ref$loadings)
  L     <- unclass(efa$unrot_loadings)

  # Sign- and column-order-invariant comparison via the reproduced common-factor
  # correlations (L %*% t(L)); the FIML correlation and the ML extraction match the
  # established corFiml -> fa workflow. Both come from the same EM -> ML workflow, so the
  # reproduced common parts agree far tighter than the eyeball 0.01: a tolerance two orders
  # of magnitude looser than they differ would let a real drift pass.
  expect_lt(max(abs(tcrossprod(L) - tcrossprod(L_ref))), 1e-3)
})

test_that("FIML resolves N to the EM case count (rows with any observed value)", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X <- fiml_mar_data()
  Xd <- rbind(X, NA)                                     # an all-missing row carries no info

  efa <- EFA(Xd, n_factors = 2, method = "ML", cor_method = "fiml")

  # N is the EM case count: every row with at least one observed value, independent of `use`.
  # This pins the FIML branch of the N resolution, distinguishing em$n from both the raw row
  # count and the (smaller) complete-case count.
  expect_equal(efa$settings$N, nrow(X))                  # the all-NA row is dropped
  expect_lt(efa$settings$N, nrow(Xd))                    # ... so N is not the raw row count
  expect_gt(efa$settings$N, sum(stats::complete.cases(Xd)))  # ... nor the complete-case count
})

test_that("the FIML correlation matches lavaan two-stage (missing = 'ml')", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")

  X <- fiml_mar_data()
  efa <- EFA(X, n_factors = 2, method = "ML", cor_method = "fiml")

  # EFA() analyses the FIML two-stage correlation, not a pairwise/Pearson one.
  Rml <- lavaan::lavCor(as.data.frame(X), missing = "ml", output = "cor")
  expect_equal(unclass(efa$orig_R), as.matrix(Rml),
               tolerance = 1e-3, ignore_attr = TRUE)
})

test_that("FIML fit indices are the corrected (Satorra-Bentler) two-stage statistic", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X <- fiml_mar_data()
  k <- 2L
  efa <- EFA(X, n_factors = k, method = "ML", cor_method = "fiml")
  fi <- efa$fit_indices

  # The reported chi-square is the scaled-shifted Satorra-Bentler statistic a * T + b, not the
  # plain two-stage likelihood-ratio statistic (which is not asymptotically chi^2(df) under the
  # two-stage estimator). Pin the self-consistency of the scaled-statistic block.
  expect_identical(fi$chi_scaled_type, "scaled.shifted")
  expect_equal(fi$chi, fi$chi_scaling * fi$chi_unscaled + fi$chi_shift, tolerance = 1e-8)
  expect_equal(fi$df, ((6 - k)^2 - (6 + k)) / 2)        # = 4, the EFA correlation-structure df
  expect_equal(fi$df_null, 6 * 5 / 2)                   # = 15
  expect_true(all(is.finite(c(fi$chi, fi$chi_unscaled, fi$chi_mean_adjusted, fi$chi_mean_var,
                              fi$df_mean_var, fi$chi_null, fi$CFI, fi$TLI, fi$RMSEA))))
  expect_gt(fi$chi, 0)

  # AIC/BIC/ECVI are likelihood-ratio chi-square criteria with no meaning for the moment-scaled
  # statistic, so they are NA (as on the sandwich path), not the chi-square-based forms.
  expect_true(all(is.na(c(fi$AIC, fi$BIC, fi$ECVI))))
})

test_that("FIML fit indices match lavaan two-stage (Satorra-Bentler)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")

  # A single-factor model: a clean oracle where EFAtools' quadratic-form statistic is close to
  # lavaan's exact two-stage likelihood ratio (the gap widens with misfit, so a saturated-ish
  # cross-loading model is a looser comparison).
  X <- fiml_mar_data_1f()
  efa <- EFA(X, n_factors = 1, method = "ML", cor_method = "fiml")
  fi  <- efa$fit_indices

  vn <- colnames(X)
  mod <- paste0("f =~ ", paste(vn, collapse = " + "))
  lfit <- lavaan::cfa(mod, data = as.data.frame(X), std.lv = TRUE, missing = "two.stage")
  sb <- lavaan::lavInspect(lfit, "test")[["satorra.bentler"]]

  # lavaan's two-stage estimator reports the Satorra-Bentler (mean-adjusted) corrected statistic;
  # EFAtools keeps it as `chi_mean_adjusted` (the reported `chi` is the scaled-shifted variant of
  # the same correction). The mean-adjusted statistic is the quadratic-form approximation of
  # lavaan's exact two-stage likelihood-ratio statistic, matching within a few percent.
  expect_equal(fi$df, sb$df)                            # EFA correlation-structure df match exactly
  expect_equal(fi$chi_mean_adjusted, sb$stat, tolerance = 0.05)

  # The robust CFI/TLI/RMSEA come from the scaled model and baseline statistics and sit near the
  # perfect-fit boundary for this correctly specified model, so they match lavaan's scaled measures.
  fm <- lavaan::fitMeasures(lfit, c("cfi.scaled", "tli.scaled", "rmsea.scaled"))
  expect_lt(abs(fi$CFI   - unname(fm[["cfi.scaled"]])),   0.02)
  expect_lt(abs(fi$TLI   - unname(fm[["tli.scaled"]])),   0.02)
  expect_lt(abs(fi$RMSEA - unname(fm[["rmsea.scaled"]])), 0.02)
})

test_that("FIML leaves the chi-square NA for PAF but keeps the residual indices", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X <- fiml_mar_data()
  efa <- EFA(X, n_factors = 2, method = "PAF", cor_method = "fiml")
  fi  <- efa$fit_indices

  # PAF has no discrepancy-based chi-square (the package convention), so the chi-square-derived
  # indices are NA under FIML as on every other cor_method; the residual indices are still
  # computed against the EM correlation.
  expect_true(is.na(fi$chi))
  expect_true(is.na(fi$CFI))
  expect_true(is.na(fi$RMSEA))
  expect_true(is.finite(fi$CAF))
  expect_true(is.finite(fi$SRMR))
})

test_that("FIML leaves AIC/BIC/ECVI NA for a just-identified (df = 0) fit", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  # A one-factor model on three indicators is just-identified (df = 0), so the corrected two-stage
  # statistic is not formed and .gof() falls back to the plain likelihood-ratio chi-square. AIC,
  # BIC, and ECVI are likelihood-ratio criteria with no meaning under the two-stage estimator, so
  # they stay NA on this fallback too -- not the finite chi-square-derived values the plain LRT
  # would otherwise produce.
  set.seed(77)
  L <- matrix(0.7, 3, 1)
  S <- tcrossprod(L); diag(S) <- 1
  X <- MASS::mvrnorm(500, mu = rep(0, 3), Sigma = S)
  colnames(X) <- paste0("V", seq_len(3))
  X[X[, 1] >  0.8, 2] <- NA
  X[X[, 1] < -0.8, 3] <- NA

  # The just-identified df = 0 model warns (efa_just_identified); assert it so the fit-index
  # check below runs on a clean condition stack.
  expect_warning(
    fi <- EFA(X, n_factors = 1, method = "ML", cor_method = "fiml")$fit_indices,
    class = "efa_just_identified"
  )
  expect_equal(fi$df, 0)
  expect_true(all(is.na(c(fi$AIC, fi$BIC, fi$ECVI))))
})

test_that("FIML NA's the whole chi-square block for a non-PD model-implied matrix", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X  <- fiml_mar_data()
  em <- .fiml_em_moments(X)

  # A Heywood loading (here > 1) makes the model-implied correlation LL' (unit diagonal) and hence
  # Sigma_model non-positive-definite, so the model deviance is undefined. The model and baseline
  # statistics must NA out together (as on every other cor_method), not report a finite baseline
  # against an undefined model.
  LLt <- tcrossprod(matrix(c(1.1, 1.0, 0.2, 0.1, 0.1, 0.1), ncol = 1))
  fc <- .gof_fiml_chisq(LLt, N = em$n, method = "ML", df = 9, m = 6,
                        fiml = list(data = X, mu = em$mu, sigma = em$sigma, logl = em$logl))

  expect_true(is.na(fc$chi))
  expect_true(is.na(fc$chi_null))
  expect_true(is.na(fc$df_null))
  expect_true(is.na(fc$chi_cfi))
  expect_true(is.na(fc$chi_null_cfi))
})

test_that("FIML fit indices ignore a fully-missing row (point-estimate filter)", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X  <- fiml_mar_data()
  Xd <- rbind(X, NA)                                    # an all-missing row carries no information

  fi  <- EFA(X,  n_factors = 2, method = "ML", cor_method = "fiml")$fit_indices
  fid <- EFA(Xd, n_factors = 2, method = "ML", cor_method = "fiml")$fit_indices

  # .gof_fiml_chisq() re-applies the EM's fully-missing-row filter to fiml$data, so the all-NA row
  # drops out of the model and baseline log-likelihoods exactly as it does from the saturated one:
  # the likelihood-ratio chi-squares are identical with or without it. A regression dropping that
  # filter would compute logl_model/logl_null over a row set wider than logl_sat and diverge here.
  expect_equal(fid$chi, fi$chi, tolerance = 1e-8)
  expect_equal(fid$chi_null, fi$chi_null, tolerance = 1e-8)
})

test_that("FIML np-boot fit-index SEs come from the corrected two-stage chi-square", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X <- fiml_mar_data()
  efa <- suppressWarnings(suppressMessages(
    EFA(X, n_factors = 2, method = "ML", cor_method = "fiml",
        se = "np-boot", b_boot = 20, seed = 42)
  ))

  # Each resample re-estimates the EM moments and forms the same corrected (Satorra-Bentler)
  # two-stage statistic as the point estimate, so the chi-square-derived fit-index SEs are
  # populated and finite (not the NA they would be if the bootstrap path silently dropped the FIML
  # moments). AIC/BIC/ECVI are NA for the scaled statistic, so their bootstrap SE is NA too.
  fi_se <- efa$SE$fit_indices
  expect_true(all(is.finite(fi_se[c("chi", "CFI", "TLI", "RMSEA")])))
  expect_gt(fi_se[["chi"]], 0)
})

test_that("FIML aborts on unsupported input/option combinations", {
  raw <- matrix(c(1, 2, 3, NA, 5, 6,
                  2, 1, 4,  3, 2, 5,
                  3, 5, 2,  4, 1, 6), ncol = 3)
  colnames(raw) <- c("a", "b", "c")
  cmat <- stats::cor(raw[stats::complete.cases(raw), ])

  # A correlation matrix carries no cases to estimate the FIML moments from.
  expect_error(EFA(cmat, n_factors = 1, N = 50, cor_method = "fiml"),
               class = "efa_fiml_needs_raw")
  # DWLS needs a polychoric asymptotic covariance the FIML correlation cannot supply.
  expect_error(EFA(raw, n_factors = 1, method = "DWLS", cor_method = "fiml"),
               class = "efa_fiml_unsupported_method")
  # PAF carries no Stage-2 discrepancy weight for the corrected two-stage sandwich, so analytic
  # standard errors are rejected for it under FIML (the bootstrap stays available for PAF).
  expect_error(EFA(raw, n_factors = 1, method = "PAF", cor_method = "fiml", se = "sandwich"),
               class = "efa_se_unsupported")
  expect_error(EFA(raw, n_factors = 1, method = "PAF", cor_method = "fiml", se = "information"),
               class = "efa_se_unsupported")
})


# Corrected two-stage (Yuan & Bentler, 2000; Savalei & Bentler, 2009) analytic standard errors for
# cor_method = "fiml" (se = "information"/"sandwich" both route to the same two-stage sandwich).

# Independent reference: the two-stage sandwich with the gauge-constraint Jacobian taken by CENTRAL
# FINITE DIFFERENCES of the orientation function, sharing no code with .se_fiml_core()'s analytic
# constraint (.se_sandwich_constraint). H is the normal-theory ML weight; the meat is the FIML
# saturated covariance, used directly (it is already Var(r-hat)).
.ref_fiml_loading_se <- function(L, Omega, gauge, weight = "ML") {
  p <- nrow(L); k <- ncol(L); pk <- p * k
  pr <- utils::combn(p, 2L); pii <- pr[1, ]; pjj <- pr[2, ]; n <- ncol(pr)
  Delta <- matrix(0, n, pk)
  for (f in seq_len(k)) {
    Delta[cbind(seq_len(n), (f - 1L) * p + pii)] <- L[pjj, f]
    Delta[cbind(seq_len(n), (f - 1L) * p + pjj)] <- L[pii, f]
  }
  # Stage-2 weight matched to the estimator: ULS the identity, ML the normal-theory weight.
  if (weight == "ULS") {
    H <- diag(n)
  } else {
    Sig <- tcrossprod(L); diag(Sig) <- 1; P <- solve(Sig)
    H <- matrix(0, n, n)
    for (s in seq_len(n)) { a <- pii[s]; b <- pjj[s]; H[, s] <- 0.5 * (P[pii, a] * P[pjj, b] + P[pii, b] * P[pjj, a]) }
  }
  A <- crossprod(Delta, H %*% Delta)
  gfun <- function(par) {
    Lm <- matrix(par, p, k)
    M <- if (gauge == "LtL") crossprod(Lm) else { ps <- 1 - rowSums(Lm^2); crossprod(Lm, Lm / ps) }
    M[upper.tri(M)]
  }
  nc <- k * (k - 1L) / 2L; Cmat <- matrix(0, nc, pk); h <- 1e-6; par0 <- as.vector(L)
  for (j in seq_len(pk)) { pp <- par0; pp[j] <- pp[j] + h; pm <- par0; pm[j] <- pm[j] - h; Cmat[, j] <- (gfun(pp) - gfun(pm)) / (2 * h) }
  Aug <- if (nc > 0) rbind(cbind(A, t(Cmat)), cbind(Cmat, matrix(0, nc, nc))) else A
  Ab <- solve(Aug)[seq_len(pk), seq_len(pk)]
  HD <- H %*% Delta
  V_AA <- Ab %*% crossprod(HD, Omega %*% HD) %*% Ab
  matrix(sqrt(diag(V_AA)), p, k)
}

test_that("FIML unrotated loading and uniqueness SEs match lavaan two-stage", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")

  X <- fiml_mar_data_1f()
  efa <- EFA(X, n_factors = 1, method = "ML", cor_method = "fiml", se = "sandwich")

  vn <- colnames(X)
  mod <- paste0("f =~ ", paste(vn, collapse = " + "))
  lfit <- lavaan::cfa(mod, data = as.data.frame(X), std.lv = TRUE, missing = "two.stage")
  lam <- lavaan::standardizedSolution(lfit)
  lam <- lam[lam$op == "=~", ]
  L_lav <- lam$est.std; SE_lav <- lam$se
  sgn <- sign(sum(as.numeric(efa$unrot_loadings) * L_lav))

  # One factor: no rotational indeterminacy, so the standardized loadings (up to sign) and their
  # corrected two-stage SEs are directly comparable. The residual gap is the EM-vs-optimiser and
  # N-vs-(N-1) difference.
  expect_equal(as.numeric(efa$unrot_loadings) * sgn, L_lav, tolerance = 0.01)
  expect_equal(as.numeric(efa$SE$unrot_loadings), SE_lav, tolerance = 0.03)
  expect_true(all(is.finite(efa$SE$uniquenesses)) && all(efa$SE$uniquenesses > 0))
})

test_that("FIML two-stage loading SEs match an independent reference sandwich", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X <- fiml_mar_data()
  efa <- EFA(X, n_factors = 2, method = "ML", cor_method = "fiml", se = "sandwich")
  L <- unclass(efa$unrot_loadings)

  em <- .fiml_em_moments(X)
  Omega <- .fiml_saturated_acov(X, em$mu, em$sigma)$cor
  # ML reports loadings in the Lambda' Psi^-1 Lambda orientation; the reference takes its gauge
  # constraint by central differences, sharing no code with .se_fiml_core().
  ref <- .ref_fiml_loading_se(L, Omega, "LtPiL")
  expect_equal(unclass(efa$SE$unrot_loadings), ref, tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("FIML rotated loading SEs match lavaan two-stage under a supported rotation", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")

  X <- fiml_mar_data()
  efa <- EFA(X, n_factors = 2, method = "ML", rotation = "geominQ",
             cor_method = "fiml", se = "sandwich")

  # Match EFAtools' geominQ epsilon (default delta = 0.01) so the two rotations coincide; the
  # two-stage estimator carries the corrected rotated-loading SEs through standardizedSolution().
  lf <- lavaan::efa(as.data.frame(X), nfactors = 2, rotation = "geomin",
                    missing = "two.stage", rotation.args = list(geomin.epsilon = 0.01))
  lfit <- if (is.list(lf) && !inherits(lf, "lavaan")) lf[[1]] else lf
  lam <- lavaan::standardizedSolution(lfit)
  lam <- lam[lam$op == "=~", ]
  SE_lav <- matrix(lam$se, nrow = 6, ncol = 2)
  EST_lav <- matrix(lam$est.std, nrow = 6, ncol = 2)

  SE_efa <- unclass(efa$SE$rot_loadings)
  EST_efa <- unclass(efa$rot_loadings)
  # The rotated column order is arbitrary; match EFAtools columns to lavaan by the salient pattern.
  cmap <- apply(abs(crossprod(EST_efa, EST_lav)), 1L, which.max)
  expect_setequal(cmap, 1:2)
  expect_equal(SE_efa[, cmap], SE_lav, tolerance = 0.05, ignore_attr = TRUE)
})

test_that("FIML sandwich SEs fill the analytic SE/CI schema (oblique)", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X <- fiml_mar_data()
  efa <- EFA(X, n_factors = 2, method = "ML", rotation = "oblimin",
             cor_method = "fiml", se = "sandwich")

  expect_named(efa$SE, c("unrot_loadings", "uniquenesses", "rot_loadings", "communalities",
                         "Phi", "Structure"))
  for (nm in c("unrot_loadings", "rot_loadings", "Structure")) {
    expect_true(all(is.finite(efa$SE[[nm]])), info = nm)
  }
  expect_true(all(efa$SE$unrot_loadings > 0))
  expect_length(efa$SE$uniquenesses, 6L)
  expect_true(all(is.finite(efa$SE$uniquenesses)))
  # Factor correlations: symmetric SE matrix with a fixed (zero-variance) unit diagonal.
  expect_equal(diag(efa$SE$Phi), c(0, 0))
  expect_equal(efa$SE$Phi, t(efa$SE$Phi))
  expect_true(is.finite(efa$SE$Phi[1, 2]))
  # Communalities are rotation-invariant, so their SEs equal the unrotated uniqueness SEs.
  expect_equal(unname(efa$SE$communalities), unname(efa$SE$uniquenesses), tolerance = 1e-8)
  # Wald intervals bracket the point estimates.
  expect_true(all(efa$CI$unrot_loadings$lower <= unclass(efa$unrot_loadings)))
  expect_true(all(efa$CI$unrot_loadings$upper >= unclass(efa$unrot_loadings)))
})

test_that("FIML 'information' and 'sandwich' give the same corrected two-stage SEs", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X <- fiml_mar_data()
  fit_i <- EFA(X, n_factors = 2, method = "ML", cor_method = "fiml", se = "information")
  fit_s <- EFA(X, n_factors = 2, method = "ML", cor_method = "fiml", se = "sandwich")
  # Under FIML both analytic settings return the corrected two-stage sandwich, so the SEs and the
  # persisted loading covariance are identical.
  expect_equal(fit_i$SE$unrot_loadings, fit_s$SE$unrot_loadings)
  expect_equal(fit_i$SE$uniquenesses, fit_s$SE$uniquenesses)
  expect_equal(fit_i$vcov_unrot_loadings, fit_s$vcov_unrot_loadings)
})

test_that("FIML ULS uses the estimator's own (identity) Stage-2 weight, not the ML weight", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X <- fiml_mar_data()
  efa_u <- EFA(X, n_factors = 2, method = "ULS", cor_method = "fiml", se = "sandwich")
  L <- unclass(efa_u$unrot_loadings)
  expect_true(all(is.finite(efa_u$SE$unrot_loadings)) && all(efa_u$SE$unrot_loadings > 0))
  expect_true(all(is.finite(efa_u$SE$uniquenesses)))

  em <- .fiml_em_moments(X)
  Omega <- .fiml_saturated_acov(X, em$mu, em$sigma)$cor
  # ULS reports loadings in the Lambda' Lambda orientation and weights the discrepancy by the
  # identity, so the reference sandwich uses the identity weight in that same gauge.
  ref_uls <- .ref_fiml_loading_se(L, Omega, "LtL", weight = "ULS")
  expect_equal(unclass(efa_u$SE$unrot_loadings), ref_uls, tolerance = 1e-5, ignore_attr = TRUE)
  # In the same gauge, the ML-weighted sandwich differs: the SE follows the estimator's weight,
  # not a hardcoded ML weight.
  ref_ml <- .ref_fiml_loading_se(L, Omega, "LtL", weight = "ML")
  expect_gt(max(abs(unclass(efa_u$SE$unrot_loadings) - ref_ml)), 1e-4)

  # 'information' and 'sandwich' coincide for ULS too (both route to the corrected sandwich).
  efa_ui <- EFA(X, n_factors = 2, method = "ULS", cor_method = "fiml", se = "information")
  expect_equal(efa_ui$SE$unrot_loadings, efa_u$SE$unrot_loadings)
})

test_that("FIML analytic SEs NA-fill with a classed warning when the covariance is unusable", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X <- fiml_mar_data()
  em <- .fiml_em_moments(X)
  fiml <- list(data = X, mu = em$mu, sigma = em$sigma)

  # A non-finite (NA) loading makes the model derivative and the parameter covariance undefined;
  # the unrotated SE schema must NA-fill with the classed efa_se_unreliable warning rather than
  # ship a finite SE next to an unusable covariance. (A mild Heywood, by contrast, the two-stage
  # sandwich handles via the Lambda'Lambda gauge.)
  fo <- list(unrot_loadings = matrix(c(0.7, NA, 0.6, 0.5, 0.4, 0.3,
                                       0.1, 0.2, 0.15, 0.7, 0.6, 0.5), 6, 2),
             orig_R = stats::cov2cor(em$sigma))
  expect_warning(
    res <- .se_fiml(fo, rot_info = NULL, N = em$n, ci = 0.95, fiml = fiml, method = "ML"),
    class = "efa_se_unreliable"
  )
  expect_true(all(is.na(res$SE$unrot_loadings)))
  expect_true(all(is.na(res$vcov_unrot_loadings)))
})

test_that("the FIML CI-provenance note names the corrected two-stage sandwich", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  local_reproducible_output()

  efa <- EFA(fiml_mar_data_1f(), n_factors = 1, method = "ML",
             cor_method = "fiml", se = "sandwich")
  body <- cli::ansi_strip(format(summary(efa)))
  note <- body[grepl("Note:", body)]
  # For a non-pooled analytic fit both "information" and "sandwich" yield the corrected two-stage
  # sandwich, so the note must describe that covariance, not the expected information matrix.
  expect_true(any(grepl("two-stage (FIML) sandwich", note, fixed = TRUE)))
  expect_false(any(grepl("expected information matrix", note)))
})

test_that("print/summary.EFA label FIML correlations in the header", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  local_reproducible_output()

  efa <- EFA(fiml_mar_data(), n_factors = 2, method = "ML", cor_method = "fiml")

  # The header flags the two-stage FIML correlation so the analysed matrix is not read as an
  # ordinary complete-case one; decimals are scrubbed so only layout/labels/wording are pinned.
  expect_snapshot(print(efa), transform = scrub_num)

  # summary() renders the same header, so the label is surfaced there as well.
  expect_match(format(summary(efa)), "FIML (two-stage, missing data)",
               fixed = TRUE, all = FALSE)
})

test_that("use is ignored under cor_method = 'fiml' (classed warning)", {
  set.seed(1)
  n <- 200
  f <- stats::rnorm(n)
  load <- 0.7
  X <- vapply(seq_len(4),
              function(j) load * f + stats::rnorm(n, sd = sqrt(1 - load^2)),
              numeric(n))
  colnames(X) <- paste0("V", seq_len(4))

  expect_warning(
    EFA(X, n_factors = 1, method = "ML", cor_method = "fiml", use = "complete.obs"),
    class = "efa_fiml_use_ignored"
  )
})

test_that("FIML np-boot returns the full SE/CI schema with finite SEs (oblique)", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X <- fiml_mar_data()
  b <- 20L

  efa <- suppressWarnings(suppressMessages(
    EFA(X, n_factors = 2, method = "ML", rotation = "oblimin",
        cor_method = "fiml", se = "np-boot", b_boot = b, seed = 42)
  ))

  expect_s3_class(efa, "EFA")

  # The bootstrap SE/CI schema is the same one the pearson/poly paths produce: an oblique
  # rotation adds the factor-correlation (Phi) and structure-coefficient blocks. The
  # replicate cube's last dimension indexes the bootstrap samples.
  expect_named(efa$SE, c("unrot_loadings", "rot_loadings", "Phi", "Structure",
                         "fit_indices", "residuals", "valid_target_rotations"))
  expect_named(efa$CI, c("unrot_loadings", "rot_loadings", "Phi", "Structure",
                         "fit_indices", "residuals"))
  expect_named(efa$replicates, c("unrot_loadings", "rot_loadings", "Phi", "Structure",
                                 "fit_indices", "residuals"))
  expect_identical(dim(efa$replicates$unrot_loadings)[3], b)

  # SEs are finite and non-negative throughout, the loadings strictly positive (Phi's fixed
  # unit diagonal carries SE 0, so it is only checked >= 0); confidence intervals are ordered.
  for (nm in c("unrot_loadings", "rot_loadings", "Phi", "Structure")) {
    expect_true(all(is.finite(efa$SE[[nm]])), info = nm)
    expect_true(all(efa$SE[[nm]] >= 0), info = nm)
    expect_true(all(efa$CI[[nm]]$lower <= efa$CI[[nm]]$upper), info = nm)
  }
  expect_true(all(efa$SE$unrot_loadings > 0))
})

test_that("FIML np-boot is reproducible given a fixed seed", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X <- fiml_mar_data()
  run <- function() suppressWarnings(suppressMessages(
    EFA(X, n_factors = 2, method = "ML", rotation = "oblimin",
        cor_method = "fiml", se = "np-boot", b_boot = 8, seed = 7)
  ))

  # The `seed` argument reseeds inside EFA(), so two calls reproduce the case resampling
  # (and hence the EM recompute and the replicate fits) independently of the caller's RNG.
  a <- run()
  b <- run()

  expect_equal(a$SE$unrot_loadings, b$SE$unrot_loadings)
  expect_equal(a$SE$Phi, b$SE$Phi)
})

test_that("a degenerate FIML resample is dropped rather than aborting np-boot", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X <- fiml_mar_data()
  b <- 10L

  # Capture the real EM engine, then mock it: succeed on the point-estimate call (the first
  # one, inside .prepare_cor_input()) and force two later bootstrap resamples to abort,
  # standing in for a degenerate resample. The loop must store an all-NA matrix for those so
  # the replicate is dropped at the fit stage rather than aborting the whole bootstrap.
  real_em <- .fiml_em_moments
  n_call <- 0L
  testthat::local_mocked_bindings(
    .fiml_em_moments = function(...) {
      n_call <<- n_call + 1L
      if (n_call %in% c(3L, 6L)) {
        cli::cli_abort("forced degenerate resample", class = "efa_fiml_not_posdef")
      }
      real_em(...)
    }
  )

  efa <- suppressWarnings(suppressMessages(
    EFA(X, n_factors = 2, method = "ML", cor_method = "fiml",
        se = "np-boot", b_boot = b, seed = 1)
  ))

  expect_s3_class(efa, "EFA")
  expect_true(all(is.finite(efa$SE$unrot_loadings)))
  # the two dropped replicates leave NA slices; the surviving ones drive the SEs.
  expect_true(any(is.na(efa$replicates$unrot_loadings)))
})

test_that("FIML np-boot resamples the informative rows when a row is fully missing", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  # A fully-missing row carries no information and is excluded from the EM case count
  # (N = em$n < nrow), so the bootstrap must resample the rows with at least one
  # observed value -- not all row positions -- to keep each replicate's sample size
  # equal to N. The np-boot path runs cleanly and N stays the informative-row count.
  X <- fiml_mar_data()
  Xd <- rbind(X, NA)

  efa <- suppressWarnings(suppressMessages(
    EFA(Xd, n_factors = 2, method = "ML", cor_method = "fiml",
        se = "np-boot", b_boot = 12, seed = 3)
  ))

  expect_s3_class(efa, "EFA")
  expect_equal(efa$settings$N, nrow(X))                 # the all-NA row is excluded
  expect_lt(efa$settings$N, nrow(Xd))
  expect_true(all(is.finite(efa$SE$unrot_loadings)))
})
