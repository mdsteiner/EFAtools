# cor_method = "fiml": two-stage / full-information ML correlations wired into EFA().
# The EM moments engine itself is covered in test-fiml-moments.R; these tests cover the
# routing through .prepare_cor_input() and EFA(), the point estimates (which match
# psych::corFiml() -> psych::fa() and lavaan(missing = "two.stage"), not lavaan::efa()),
# and the classed conditions guarding the unsupported input/option combinations.

# A clean two-factor population covariance: items 1..p/2 load on factor 1 at `loadings[1]`, the
# rest on factor 2 at `loadings[2]`, with unit diagonal.
fiml_pop_cov <- function(p = 6, loadings = c(0.7, 0.7)) {
  L <- matrix(0, p, 2)
  half <- p / 2
  L[seq_len(half), 1] <- loadings[1]
  L[(half + 1):p, 2] <- loadings[2]
  S <- tcrossprod(L)
  diag(S) <- 1
  S
}

# A missing-at-random fixture: column 1 is fully observed and drives the missingness in the
# others, so the mechanism depends only on observed data. The multivariate normal draw goes
# through the Cholesky factor rather than an eigendecomposition: at the default (equal) loadings
# the population covariance has repeated eigenvalues, whose eigenvector basis is mathematically
# undetermined and therefore settled by rounding, while chol() is unique for a positive definite
# matrix, so the same seed yields the same sample on every LAPACK build.
fiml_mar_data <- function(n = 800, seed = 456, loadings = c(0.7, 0.7)) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * 6), n) %*% chol(fiml_pop_cov(6, loadings))
  colnames(X) <- paste0("V", seq_len(6))
  X[X[, 1] >  0.8, 2] <- NA
  X[X[, 1] < -0.8, 3] <- NA
  X[X[, 1] >  1.2, 4] <- NA
  X[X[, 1] < -1.2, 5] <- NA
  X
}

# The same fixture with the two factors at different strengths. At the default equal loadings the
# canonical variances diag(Lambda' Psi^-1 Lambda) that identify the unrotated ML solution coincide
# exactly in the population: the loadings can then be rotated within that plane without leaving the
# identifying constraint, so the unrotated orientation -- and with it the unrotated loading standard
# errors -- has no well-defined sampling distribution, and a fitted solution sits right at the
# detector's floor. Every test of an unrotated analytic standard error needs a solution where that
# orientation is determined, which separating the two strengths provides (canonical variances 4.83
# against 0.90, a transversal an order of magnitude clear of the floor). Rotated quantities are
# gauge-invariant and need no such care.
fiml_mar_data_2f <- function(n = 800, seed = 456) {
  fiml_mar_data(n = n, seed = seed, loadings = c(0.8, 0.5))
}

# A single-factor MAR fixture: a clean SE / chi-square oracle, with no rotational indeterminacy.
fiml_mar_data_1f <- function(n = 1000, p = 6, load = 0.65, seed = 11) {
  set.seed(seed)
  L <- matrix(load, p, 1)
  S <- tcrossprod(L); diag(S) <- 1
  X <- matrix(stats::rnorm(n * p), n) %*% chol(S)
  colnames(X) <- paste0("V", seq_len(p))
  X[X[, 1] >  0.8, 2] <- NA
  X[X[, 1] < -0.8, 3] <- NA
  X[X[, 1] >  1.0, 4] <- NA
  X[X[, 1] < -1.0, 5] <- NA
  X
}

test_that("complete data: FIML loadings equal the Pearson-correlation loadings", {
  skip_on_cran()

  set.seed(123)
  X <- matrix(stats::rnorm(500 * 6), 500) %*% chol(fiml_pop_cov(6))
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

  X <- fiml_mar_data()
  efa <- EFA(X, n_factors = 2, method = "ML", cor_method = "fiml")

  # EFA() analyses the FIML two-stage correlation, not a pairwise/Pearson one.
  Rml <- lavaan::lavCor(as.data.frame(X), missing = "ml", output = "cor")
  expect_equal(unclass(efa$orig_R), as.matrix(Rml),
               tolerance = 1e-3, ignore_attr = TRUE)
})

test_that("FIML fit indices are the corrected (Satorra-Bentler) two-stage statistic", {
  skip_on_cran()

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

# Force the Stage-1 saturated covariance to fail, so the corrected two-stage statistic cannot be
# formed. It is the covariance that is degenerate, not the model, so the plain two-stage
# likelihood-ratio statistic is still available -- which is exactly the situation the fallback
# contract covers. Used by the tests below to reach it deterministically.
.mock_fiml_acov_failure <- function(env = parent.frame()) {
  testthat::local_mocked_bindings(
    .fiml_saturated_acov = function(...) {
      cli::cli_abort("forced degenerate saturated covariance",
                     class = "efa_fiml_singular_information")
    },
    .env = env
  )
}

test_that("a formable correction is tagged as the corrected statistic, with no fallback warning", {
  skip_on_cran()

  expect_no_warning(
    fi <- efa_fit(fiml_mar_data(), n_factors = 2, estimator = "ML",
                  cor_method = "fiml")$fit_indices,
    class = "efa_fiml_uncorrected_chisq"
  )

  # The success side of the tag: `chi` is the corrected statistic, and the components that
  # describe the correction are present alongside it.
  expect_identical(fi$chi_scaled_type, "scaled.shifted")
  expect_true(all(is.finite(c(fi$chi_scaling, fi$chi_shift, fi$chi_unscaled))))
})

test_that("a degenerate saturated covariance falls back to the tagged plain two-stage LRT", {
  skip_on_cran()

  X <- fiml_mar_data()
  .mock_fiml_acov_failure()

  expect_warning(
    fit <- efa_fit(X, n_factors = 2, estimator = "ML", cor_method = "fiml"),
    class = "efa_fiml_uncorrected_chisq"
  )
  fi <- fit$fit_indices

  # Tagged, so no consumer can read the substituted statistic as the Satorra-Bentler one ...
  expect_identical(fi$chi_scaled_type, "uncorrected.lrt")
  # ... and none of the scaling components ships, because nothing was rescaled.
  expect_null(fi$chi_scaling)
  expect_null(fi$chi_shift)
  expect_null(fi$chi_unscaled)
  expect_null(fi$chi_mean_adjusted)

  # The reported statistic is the plain two-stage likelihood ratio, referenced to the
  # chi-square(df) tail, and the block built on it is complete rather than NA-ed away: the
  # statistic is approximate under the two-stage estimator, not meaningless.
  em <- .fiml_em_moments(X)
  ref <- .gof_fiml_chisq(tcrossprod(unclass(fit$unrot_loadings)), N = fit$settings$N,
                         method = "ML", df = fi$df, m = ncol(X),
                         fiml = list(data = X, mu = em$mu, sigma = em$sigma, logl = em$logl))
  expect_equal(fi$chi, ref$chi, tolerance = 1e-8)
  expect_equal(fi$p_chi, stats::pchisq(fi$chi, fi$df, lower.tail = FALSE), tolerance = 1e-10)
  expect_true(all(is.finite(c(fi$CFI, fi$TLI, fi$RMSEA))))

  # AIC/BIC/ECVI stay withheld: they are likelihood-ratio criteria the two-stage estimator
  # cannot support, whether or not the correction was formed.
  expect_true(all(is.na(c(fi$AIC, fi$BIC, fi$ECVI))))
})

test_that("the fallback chi-square prints as uncorrected, never as scaled", {
  skip_on_cran()
  local_reproducible_output()

  .mock_fiml_acov_failure()
  expect_warning(
    fit <- efa_fit(fiml_mar_data(), n_factors = 2, estimator = "ML", cor_method = "fiml"),
    class = "efa_fiml_uncorrected_chisq"
  )

  printed <- cli::ansi_strip(capture.output(print(fit)))
  expect_true(any(grepl("uncorrected χ²(", printed, fixed = TRUE)))
  expect_false(any(grepl("scaled χ²(", printed, fixed = TRUE)))
})

test_that("a just-identified FIML fit is left untagged rather than called uncorrected", {
  skip_on_cran()

  # With df = 0 there is no chi-square test to mislabel, so the fallback tag and its warning
  # must not fire: the tag exists to distinguish a corrected statistic from an uncorrected
  # one, and this fit has neither.
  set.seed(77)
  L <- matrix(0.7, 3, 1)
  S <- tcrossprod(L); diag(S) <- 1
  X <- matrix(stats::rnorm(500 * 3), 500) %*% chol(S)
  colnames(X) <- paste0("V", seq_len(3))
  X[X[, 1] >  0.8, 2] <- NA
  X[X[, 1] < -0.8, 3] <- NA

  # Collect the condition classes rather than suppressing wholesale: the fit legitimately
  # raises efa_just_identified, and a blanket suppressWarnings() would also hide the very
  # condition this test asserts is absent.
  seen <- character()
  fi <- withCallingHandlers(
    efa_fit(X, n_factors = 1, estimator = "ML", cor_method = "fiml")$fit_indices,
    warning = function(w) {
      seen <<- c(seen, class(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true("efa_just_identified" %in% seen)
  expect_false("efa_fiml_uncorrected_chisq" %in% seen)
  expect_equal(fi$df, 0)
  expect_null(fi$chi_scaled_type)
})

test_that("FIML fit indices match lavaan two-stage (Satorra-Bentler)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

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

  # A one-factor model on three indicators is just-identified (df = 0), so the corrected two-stage
  # statistic is not formed and .gof() falls back to the plain likelihood-ratio chi-square. AIC,
  # BIC, and ECVI are likelihood-ratio criteria with no meaning under the two-stage estimator, so
  # they stay NA on this fallback too -- not the finite chi-square-derived values the plain LRT
  # would otherwise produce.
  set.seed(77)
  L <- matrix(0.7, 3, 1)
  S <- tcrossprod(L); diag(S) <- 1
  X <- matrix(stats::rnorm(500 * 3), 500) %*% chol(S)
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

test_that("a converged FIML fit reports its Stage-1 EM diagnostics", {
  skip_on_cran()

  X <- fiml_mar_data()
  fit <- efa_fit(X, n_factors = 2, estimator = "ML", cor_method = "fiml")

  expect_named(fit$fiml, c("converged", "iter", "n_patterns", "n"))
  expect_true(fit$fiml$converged)
  expect_true(fit$fiml$iter >= 1L && fit$fiml$iter < 500L)
  # The pattern count and the informative-row count are the EM's own, so they must match a
  # direct read of the missingness structure rather than the raw data dimensions.
  obs <- !is.na(X)
  keep <- rowSums(obs) > 0L
  expect_equal(fit$fiml$n_patterns, length(.fiml_patterns(obs[keep, , drop = FALSE])))
  expect_equal(fit$fiml$n, sum(keep))
  expect_equal(fit$fiml$n, fit$settings$N)

  # The slot belongs to the FIML route alone; no other correlation method estimates them.
  pear <- efa_fit(X, n_factors = 2, estimator = "ML", cor_method = "pearson")
  expect_null(pear$fiml)
})

test_that("estimate_control's FIML EM knobs reach the moment estimation", {
  skip_on_cran()

  X <- fiml_mar_data()

  # A cap of two iterations is reached on this fixture, so the warning fires and the fit
  # records that the analysed matrix is the last iterate rather than the FIML estimate.
  expect_warning(
    capped <- efa_fit(X, n_factors = 2, estimator = "ML", cor_method = "fiml",
                      estimate_control = estimate_control(fiml_max_iter = 2)),
    class = "efa_fiml_em_nonconvergence"
  )
  expect_false(capped$fiml$converged)
  expect_equal(capped$fiml$iter, 2L)

  # The knobs are not inert: a capped EM stops at a different correlation matrix from the
  # converged one, and a looser tolerance converges in fewer iterations.
  full <- efa_fit(X, n_factors = 2, estimator = "ML", cor_method = "fiml")
  expect_false(isTRUE(all.equal(unclass(capped$orig_R), unclass(full$orig_R))))

  loose <- efa_fit(X, n_factors = 2, estimator = "ML", cor_method = "fiml",
                   estimate_control = estimate_control(fiml_tol = 1e-2))
  expect_true(loose$fiml$converged)
  expect_lt(loose$fiml$iter, full$fiml$iter)
})

test_that("a non-converged FIML EM is visible in the printed output", {
  skip_on_cran()
  local_reproducible_output()

  expect_warning(
    capped <- efa_fit(fiml_mar_data(), n_factors = 2, estimator = "ML", cor_method = "fiml",
                      estimate_control = estimate_control(fiml_max_iter = 2)),
    class = "efa_fiml_em_nonconvergence"
  )

  printed <- cli::ansi_strip(capture.output(print(capped)))
  expect_true(any(grepl("FIML EM stopped after 2 iterations without converging",
                        printed, fixed = TRUE)))
  # summary() renders the same header, so the diagnostic is surfaced there as well.
  expect_match(cli::ansi_strip(format(summary(capped))),
               "without converging", fixed = TRUE, all = FALSE)

  # A converged fit says nothing, so the line marks the exception rather than the rule.
  ok <- cli::ansi_strip(capture.output(print(
    efa_fit(fiml_mar_data(), n_factors = 2, estimator = "ML", cor_method = "fiml"))))
  expect_false(any(grepl("without converging", ok, fixed = TRUE)))
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
# saturated covariance Omega (already on the variance scale Var(r-hat)), put on the unit
# asymptotic-variance scale by N and the covariance divided by N - 1, the package-wide small-sample
# convention shared with the expected-information and polychoric/ADF sandwich paths.
.ref_fiml_loading_se <- function(L, Omega, gauge, N, weight = "ML") {
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
  V_AA <- Ab %*% crossprod(HD, (N * Omega) %*% HD) %*% Ab / (N - 1)
  matrix(sqrt(diag(V_AA)), p, k)
}

test_that("FIML unrotated loading and uniqueness SEs match lavaan two-stage", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

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

  # The unrotated loading SEs are gauge-dependent, so this needs the fixture whose rotational
  # orientation is determined.
  X <- fiml_mar_data_2f()
  efa <- efa_fit(X, n_factors = 2, estimator = "ML", cor_method = "fiml", se = "sandwich")
  L <- unclass(efa$unrot_loadings)

  em <- .fiml_em_moments(X)
  Omega <- .fiml_saturated_acov(X, em$mu, em$sigma)$cor
  # ML reports loadings in the Lambda' Psi^-1 Lambda orientation; the reference takes its gauge
  # constraint by central differences, sharing no code with .se_fiml_core().
  ref <- .ref_fiml_loading_se(L, Omega, "LtPiL", N = nrow(X))
  expect_equal(unclass(efa$SE$unrot_loadings), ref, tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("FIML rotated loading SEs match lavaan two-stage under a supported rotation", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  X <- fiml_mar_data_2f()
  efa <- efa_fit(X, n_factors = 2, estimator = "ML", rotation = "geominQ",
                 cor_method = "fiml", se = "sandwich")

  # Match EFAtools' geominQ epsilon (default delta = 0.01) so the two rotations coincide; the
  # two-stage estimator carries the corrected rotated-loading SEs through standardizedSolution().
  if (utils::packageVersion("lavaan") >= "0.7.0") {
    lf <- lavaan::efa(as.data.frame(X), nfactors = 2,
                      rotation = list("geomin", geomin_epsilon = 0.01),
                      missing = "two.stage")
  } else {
    # lavaan < 0.7 used a separate rotation.args argument.
    lf <- lavaan::efa(as.data.frame(X), nfactors = 2, rotation = "geomin",
                      missing = "two.stage", rotation.args = list(geomin.epsilon = 0.01))
  }
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

  X <- fiml_mar_data_2f()
  efa <- efa_fit(X, n_factors = 2, estimator = "ML", rotation = "oblimin",
                 cor_method = "fiml", se = "sandwich")

  expect_named(efa$SE, c("unrot_loadings", "uniquenesses", "rot_loadings", "communalities",
                         "Phi", "Structure"))
  for (nm in c("unrot_loadings", "rot_loadings", "Structure")) {
    expect_true(all(is.finite(efa$SE[[nm]])), info = nm)
  }
  expect_true(all(efa$SE$unrot_loadings > 0))
  expect_length(efa$SE$uniquenesses, 6L)
  expect_true(all(is.finite(efa$SE$uniquenesses)))
  # Factor correlations: symmetric SE matrix with a fixed (zero-variance) unit diagonal. It carries
  # the factor names, so compare the diagonal on its values alone.
  expect_equal(unname(diag(efa$SE$Phi)), c(0, 0))
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

  X <- fiml_mar_data_2f()
  fit_i <- efa_fit(X, n_factors = 2, estimator = "ML", cor_method = "fiml", se = "information")
  fit_s <- efa_fit(X, n_factors = 2, estimator = "ML", cor_method = "fiml", se = "sandwich")
  # Under FIML both analytic settings return the corrected two-stage sandwich, so the SEs and the
  # persisted loading covariance are identical.
  expect_equal(fit_i$SE$unrot_loadings, fit_s$SE$unrot_loadings)
  expect_equal(fit_i$SE$uniquenesses, fit_s$SE$uniquenesses)
  expect_equal(fit_i$vcov_unrot_loadings, fit_s$vcov_unrot_loadings)
})

test_that("FIML ULS uses the estimator's own (identity) Stage-2 weight, not the ML weight", {
  skip_on_cran()

  X <- fiml_mar_data_2f()
  efa_u <- efa_fit(X, n_factors = 2, estimator = "ULS", cor_method = "fiml", se = "sandwich")
  L <- unclass(efa_u$unrot_loadings)
  expect_true(all(is.finite(efa_u$SE$unrot_loadings)) && all(efa_u$SE$unrot_loadings > 0))
  expect_true(all(is.finite(efa_u$SE$uniquenesses)))

  em <- .fiml_em_moments(X)
  Omega <- .fiml_saturated_acov(X, em$mu, em$sigma)$cor
  # ULS reports loadings in the Lambda' Lambda orientation and weights the discrepancy by the
  # identity, so the reference sandwich uses the identity weight in that same gauge.
  ref_uls <- .ref_fiml_loading_se(L, Omega, "LtL", N = nrow(X), weight = "ULS")
  expect_equal(unclass(efa_u$SE$unrot_loadings), ref_uls, tolerance = 1e-5, ignore_attr = TRUE)
  # In the same gauge, the SE must match the estimator-specific ULS reference more closely than
  # the ML-weighted sandwich reference.
  ref_ml <- .ref_fiml_loading_se(L, Omega, "LtL", N = nrow(X), weight = "ML")
  err_uls <- max(abs(unclass(efa_u$SE$unrot_loadings) - ref_uls))
  err_ml <- max(abs(unclass(efa_u$SE$unrot_loadings) - ref_ml))
  expect_lt(err_uls, err_ml)

  # 'information' and 'sandwich' coincide for ULS too (both route to the corrected sandwich).
  efa_ui <- efa_fit(X, n_factors = 2, estimator = "ULS", cor_method = "fiml", se = "information")
  expect_equal(efa_ui$SE$unrot_loadings, efa_u$SE$unrot_loadings)
})

test_that("FIML analytic SEs NA-fill with a classed warning when the covariance is unusable", {
  skip_on_cran()

  X <- fiml_mar_data()
  em <- .fiml_em_moments(X)
  fiml <- list(data = X, mu = em$mu, sigma = em$sigma)

  # A non-finite (NA) loading makes the model derivative and the parameter covariance undefined;
  # the unrotated SE schema must NA-fill with the classed efa_se_unreliable warning rather than
  # ship a finite SE next to an unusable covariance. This is the unusable-covariance route, which
  # is distinct from the boundary route tested above: there the covariance is fine and it is the
  # Wald approximation that is not.
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

test_that("the FIML sandwich picks its gauge from the loadings, not from the estimator label", {
  skip_on_cran()

  X <- fiml_mar_data_2f()
  em <- .fiml_em_moments(X)
  fiml <- list(data = X, mu = em$mu, sigma = em$sigma)
  Omega <- .fiml_saturated_acov(X, em$mu, em$sigma)$cor

  # ML-fitted loadings are reported in the Lambda' Psi^-1 Lambda orientation. Ask for the ULS
  # Stage-2 weight on those same loadings: the gauge must still be the one the loadings are
  # actually in, because it is detected from the fitted uniquenesses rather than assumed from the
  # estimator name -- the unrotated loading SEs are scaled in whichever gauge is used, so picking
  # the other one would report them in an orientation the solution is not in.
  L <- unclass(efa_fit(X, n_factors = 2, estimator = "ML", cor_method = "fiml")$unrot_loadings)
  fo <- list(unrot_loadings = L, orig_R = stats::cov2cor(em$sigma))
  core <- .se_fiml_core(fo, fiml, N = em$n, method = "ULS")

  ref_ltpil <- .ref_fiml_loading_se(L, Omega, "LtPiL", N = em$n, weight = "ULS")
  ref_ltl <- .ref_fiml_loading_se(L, Omega, "LtL", N = em$n, weight = "ULS")
  expect_equal(core$loadings_se, ref_ltpil, tolerance = 1e-5, ignore_attr = TRUE)
  # The two gauges really do differ here, so the check above is not satisfied by both.
  expect_gt(max(abs(ref_ltpil - ref_ltl)), 1e-3)
})

test_that("the FIML sandwich reports a degenerate rotational gauge", {
  skip_on_cran()

  X <- fiml_mar_data_2f()

  # A weakly determined orientation inflates the unrotated loading SEs without bound while leaving
  # every gauge-invariant quantity intact. This fixture's orientation is well determined, so the
  # transversal is forced below its floor: the mock is then the only cause, and what is under test
  # is that the two-stage path acts on the diagnostic at all, which it previously could not.
  expect_true(all(is.finite(
    efa_fit(X, n_factors = 2, estimator = "ML", cor_method = "fiml",
            se = "sandwich")$SE$unrot_loadings)))

  testthat::local_mocked_bindings(.gauge_transversal = function(...) 0)

  expect_warning(
    fit <- efa_fit(X, n_factors = 2, estimator = "ML", cor_method = "fiml", se = "sandwich"),
    class = "efa_se_unreliable"
  )

  # The unrotated loadings alone are withheld ...
  expect_true(all(is.na(fit$SE$unrot_loadings)))
  # ... while the gauge-invariant uniquenesses and communalities, and the covariance itself, stay:
  # the divergence lives purely in the gauge, so it must not take them down with it.
  expect_true(all(is.finite(fit$SE$uniquenesses)))
  expect_true(all(is.finite(fit$SE$communalities)))
  expect_true(all(is.finite(fit$vcov_unrot_loadings)))
})

test_that("the FIML sandwich withholds standard errors at a uniqueness pinned at the floor", {
  skip_on_cran()

  X <- fiml_mar_data()
  em <- .fiml_em_moments(X)
  fiml <- list(data = X, mu = em$mu, sigma = em$sigma)
  R <- stats::cov2cor(em$sigma)

  # The two-stage sandwich reports the same Wald quantities as the other analytic paths, from a
  # different covariance, so it withholds on the same boundary test: on the parameter-space
  # boundary the Wald approximation fails for every parameter regardless of how the covariance was
  # estimated. The boundary an ML/ULS fitter can actually reach is the uniqueness floor, not zero,
  # so this solution is strictly interior to the gauge detection's own `psi <= 0` test.
  L <- matrix(0, 6, 2)
  L[, 1] <- c(sqrt(1 - .uniqueness_floor), 0.6, 0.5, 0.1, 0.15, 0.2)
  L[, 2] <- c(0.00, 0.05, 0.20, 0.70, 0.60, 0.55)
  psi <- 1 - rowSums(L^2)
  expect_lte(min(psi), .uniqueness_floor + sqrt(.Machine$double.eps))
  expect_gt(min(psi), 0)

  fo <- list(unrot_loadings = L, orig_R = R)

  # Without the gate this covariance is perfectly usable, so what withholds the SEs below is the
  # boundary test and not a degenerate saturated covariance.
  core <- .se_fiml_core(fo, fiml, N = em$n, method = "ML")
  expect_true(core$reliable)
  expect_true(all(is.finite(core$loadings_se)))

  expect_warning(
    out <- .se_fiml(fo, rot_info = NULL, N = em$n, ci = 0.95, fiml = fiml, method = "ML"),
    class = "efa_se_unreliable"
  )
  expect_true(all(is.na(out$SE$unrot_loadings)))
  expect_true(all(is.na(out$SE$uniquenesses)))
  expect_true(all(is.na(out$SE$communalities)))
  expect_true(all(is.na(out$vcov_unrot_loadings)))
  expect_true(all(is.na(unlist(out$CI))))

  # Under a rotation the withholding reaches the rotated quantities too: they are propagated from
  # the same loading covariance, which the gate leaves unusable.
  rot_info <- list(rotation = "oblimin", rotmat = diag(2), rot_loadings = L,
                   Phi = diag(2), normalize = FALSE, crit_args = list(gam = 0, delta = 0.01))
  expect_warning(
    rot <- .se_fiml(fo, rot_info, N = em$n, ci = 0.95, fiml = fiml, method = "ML"),
    class = "efa_se_unreliable"
  )
  expect_true(all(is.na(rot$SE$rot_loadings)))
  expect_true(all(is.na(rot$SE$Phi)))
  expect_true(all(is.na(rot$SE$Structure)))

  # The corrected two-stage chi-square is a discrepancy-function quantity, not a Wald one, so the
  # boundary does not invalidate it: it is built outside the standard-error route and is kept.
  st <- .fiml_scaled_test(L, R, N = em$n, method = "ML", df = 4, m = 6, fiml = fiml)
  expect_false(is.null(st))
  expect_true(is.finite(st$chi))
})

test_that("the FIML fit indices honour ci = FALSE", {
  skip_on_cran()

  X <- fiml_mar_data()
  em <- .fiml_em_moments(X)
  fiml <- list(data = X, mu = em$mu, sigma = em$sigma, logl = em$logl)
  R <- stats::cov2cor(em$sigma)
  L <- unclass(efa_fit(X, n_factors = 2, estimator = "ML", cor_method = "fiml")$unrot_loadings)

  with_ci <- .gof(L, R, em$n, "ML", Fm = NA_real_, ci = TRUE, fiml = fiml)
  no_ci <- .gof(L, R, em$n, "ML", Fm = NA_real_, ci = FALSE, fiml = fiml)

  # The scaled block replaces the whole chi-square-derived block, so it has to honour the request
  # not to solve the RMSEA noncentrality bounds -- which is exactly the per-replicate bootstrap
  # path, where the bounds are dropped again by the aggregation.
  expect_true(any(is.finite(c(with_ci$RMSEA_LB, with_ci$RMSEA_UB))))
  expect_true(is.na(no_ci$RMSEA_LB))
  expect_true(is.na(no_ci$RMSEA_UB))

  # Nothing else changes with it: the flag governs the bounds alone.
  expect_identical(with_ci$chi_scaled_type, "scaled.shifted")
  expect_identical(no_ci$chi_scaled_type, "scaled.shifted")
  expect_equal(no_ci$chi, with_ci$chi)
  expect_equal(no_ci$RMSEA, with_ci$RMSEA)
  expect_equal(no_ci$CFI, with_ci$CFI)
})

test_that("the FIML CI-provenance note names the corrected two-stage sandwich", {
  skip_on_cran()
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

test_that("print/summary.efa label FIML correlations in the header", {
  skip_on_cran()
  local_reproducible_output()

  efa <- EFA(fiml_mar_data(), n_factors = 2, method = "ML", cor_method = "fiml")

  # The header flags the two-stage FIML correlation so the analysed matrix is not read as an
  # ordinary complete-case one. Assert the labels without pinning platform-dependent spacing.
  printed <- cli::ansi_strip(capture.output(print(efa)))
  expect_true(any(grepl(
    "Correlations: FIML (two-stage, missing data)", printed, fixed = TRUE
  )))
  expect_true(any(grepl("Unrotated Loadings", printed, fixed = TRUE)))

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
                         "fit_indices", "residuals", "valid_replicates",
                         "valid_target_rotations"))
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

test_that("the FIML bootstrap carries resample indices, not a matrix per replicate", {
  skip_on_cran()

  X <- fiml_mar_data()
  b <- 8L

  # Intercept what efa_fit() hands the replicate fits, then run them for real: the per-replicate
  # list is a closure global of the parallel map, so whatever it carries is held once per worker
  # on top of the parent.
  real_boot <- .boot_fun
  captured <- NULL
  testthat::local_mocked_bindings(
    .boot_fun = function(..., fiml_list = NULL, fiml_data = NULL) {
      captured <<- list(fiml_list = fiml_list, fiml_data = fiml_data)
      real_boot(..., fiml_list = fiml_list, fiml_data = fiml_data)
    }
  )

  fit <- suppressWarnings(suppressMessages(
    efa_fit(X, n_factors = 2, estimator = "ML", cor_method = "fiml",
            se = "np-boot", b_boot = b, seed = 42)
  ))
  expect_s3_class(fit, "efa")

  reps <- Filter(Negate(is.null), captured$fiml_list)
  expect_length(reps, b)
  for (r in reps) {
    expect_named(r, c("rows", "mu", "sigma", "logl"))
    expect_null(r$data)
    expect_length(r$rows, fit$settings$N)
  }
  # The raw data travels once, alongside the indices that point into it.
  expect_equal(captured$fiml_data, X)

  # The indices identify exactly the resample the stored moments were estimated from, so the
  # replicate fit sees the same data the old per-replicate matrix carried.
  em_1 <- suppressWarnings(.fiml_em_moments(captured$fiml_data[reps[[1]]$rows, , drop = FALSE]))
  expect_equal(em_1$sigma, reps[[1]]$sigma)
  expect_equal(em_1$mu, reps[[1]]$mu)
  expect_equal(em_1$logl, reps[[1]]$logl)

  # ... and it is materially smaller than storing the slices themselves.
  with_data <- lapply(reps, function(r) {
    list(data = X[r$rows, , drop = FALSE], mu = r$mu, sigma = r$sigma, logl = r$logl)
  })
  expect_lt(as.numeric(utils::object.size(captured$fiml_list)),
            as.numeric(utils::object.size(with_data)) / 3)
})

test_that("the FIML bootstrap replicate fit is unchanged by the index representation", {
  skip_on_cran()

  X <- fiml_mar_data()
  rows <- which(rowSums(!is.na(X)) > 0L)
  set.seed(99)
  ind <- sample(rows, size = length(rows), replace = TRUE)
  em_i <- suppressWarnings(.fiml_em_moments(X[ind, , drop = FALSE]))
  R_i <- stats::cov2cor(em_i$sigma)

  # Drive one replicate through .boot_fun() with the index representation ...
  boot <- .boot_fun(array(R_i, c(ncol(X), ncol(X), 1L)), 1L, .estimate_model,
                    method = "ML", n_factors = 2, N = length(rows), type = "EFAtools",
                    start_method = "psych", lean = TRUE,
                    fiml_list = list(list(rows = ind, mu = em_i$mu, sigma = em_i$sigma,
                                          logl = em_i$logl)),
                    fiml_data = X)

  # ... and independently with the resampled matrix supplied directly, as the caller used to
  # store it. The fit -- loadings and the FIML likelihood-ratio fit indices alike -- must be
  # identical: the representation is a storage change, not a numerical one.
  direct <- suppressWarnings(.estimate_model(
    R_i, method = "ML", n_factors = 2, N = length(rows), type = "EFAtools",
    start_method = "psych", lean = TRUE,
    fiml = list(data = X[ind, , drop = FALSE], mu = em_i$mu, sigma = em_i$sigma,
                logl = em_i$logl)))

  expect_equal(boot[[1]]$unrot_loadings, direct$unrot_loadings)
  expect_equal(boot[[1]]$fit_indices, direct$fit_indices)
})

test_that("FIML np-boot is reproducible given a fixed seed", {
  skip_on_cran()

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
