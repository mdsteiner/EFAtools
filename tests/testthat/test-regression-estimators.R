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

# Fit-index parity: the model chi-square is the Bartlett-corrected ML discrepancy at the
# model-implied matrix. For ML this reproduces stats::factanal / psych::fa(fm = "ml"); for
# ULS it reproduces psych::fa(fm = "minres"). lavaan uses the multiplier N (no Bartlett
# correction) but the same fit function, so its statistic differs by exactly bart / N.
test_that("ML chi-square reproduces stats::factanal and psych::fa(fm = 'ml')", {
  skip_on_cran()
  skip_if_not_installed("psych")

  for (fx in reg_fixtures) {
    efa <- suppressWarnings(EFA(fx$R, n_factors = fx$k, N = fx$N, method = "ML"))
    ref_fa <- stats::factanal(covmat = fx$R, factors = fx$k, n.obs = fx$N,
                              rotation = "none")
    ref_ps <- suppressMessages(suppressWarnings(
      psych::fa(fx$R, nfactors = fx$k, n.obs = fx$N, fm = "ml", rotate = "none")))

    # observed agreement on the baseline model: ~1e-9
    expect_equal(efa$fit_indices$chi, unname(ref_fa$STATISTIC), tolerance = 1e-4)
    expect_equal(efa$fit_indices$chi, unname(ref_ps$STATISTIC), tolerance = 1e-4)
    expect_equal(efa$fit_indices$df, ref_fa$dof)
  }
})

test_that("ULS chi-square reproduces psych::fa(fm = 'minres')", {
  skip_on_cran()
  skip_if_not_installed("psych")

  for (fx in reg_fixtures) {
    efa <- suppressWarnings(EFA(fx$R, n_factors = fx$k, N = fx$N, method = "ULS"))
    ref <- suppressMessages(suppressWarnings(
      psych::fa(fx$R, nfactors = fx$k, n.obs = fx$N, fm = "minres", rotate = "none")))

    # observed agreement on the baseline model: ~1e-5
    expect_equal(efa$fit_indices$chi, unname(ref$STATISTIC), tolerance = 1e-3)
  }
})

test_that("ML chi-square relates to lavaan's via the documented Bartlett factor", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  for (fx in reg_fixtures) {
    p <- ncol(fx$R)
    bart <- fx$N - 1 - (2 * p + 5) / 6 - (2 * fx$k) / 3
    efa <- suppressWarnings(EFA(fx$R, n_factors = fx$k, N = fx$N, method = "ML"))

    # lavaan::efa()'s baseline / robust paths can fail on a covariance-matrix input;
    # take the standard test statistic from the fitted model and skip if unavailable.
    lav_stat <- tryCatch({
      fit <- lavaan::efa(sample.cov = fx$R, sample.nobs = fx$N, nfactors = fx$k,
                         rotation = "none", estimator = "ML", se = "none")[[1]]
      lavaan::lavInspect(fit, "test")$standard$stat
    }, error = function(e) NA_real_)
    skip_if(is.na(lav_stat), "lavaan::efa() did not return a standard test statistic")

    expect_equal(efa$fit_indices$chi, lav_stat * bart / fx$N, tolerance = 1e-3)
  }
})

test_that("SRMR reproduces lavaan", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  for (fx in reg_fixtures) {
    efa <- suppressWarnings(EFA(fx$R, n_factors = fx$k, N = fx$N, method = "ML"))

    lav_srmr <- tryCatch({
      fit <- lavaan::efa(sample.cov = fx$R, sample.nobs = fx$N, nfactors = fx$k,
                         rotation = "none", estimator = "ML", se = "none")[[1]]
      unname(lavaan::fitMeasures(fit, "srmr"))
    }, error = function(e) NA_real_)
    skip_if(is.na(lav_srmr), "lavaan::efa() did not return SRMR")

    # SRMR is convention-free, so it matches lavaan exactly
    expect_equal(efa$fit_indices$SRMR, lav_srmr, tolerance = 1e-5)
  }
})

test_that("TLI equals its (N - 1)-scale closed form; ECVI uses the Bartlett-corrected chi", {
  for (fx in reg_fixtures) {
    p <- ncol(fx$R)
    for (mth in c("ML", "ULS")) {
      fi <- suppressWarnings(EFA(fx$R, n_factors = fx$k, N = fx$N,
                                 method = mth))$fit_indices

      # CFI/TLI compare the model and baseline noncentralities on a common (N - 1) scale,
      # so the reported (Bartlett-corrected) model and baseline chi-squares are rescaled
      # from their own multipliers to (N - 1).
      mult_m    <- fx$N - 1 - (2 * p + 5) / 6 - (2 * fx$k) / 3
      mult_null <- fx$N - 1 - (2 * p + 5) / 6
      chi_cfi      <- fi$chi      * (fx$N - 1) / mult_m
      chi_null_cfi <- fi$chi_null * (fx$N - 1) / mult_null
      expect_equal(fi$TLI,
                   ((chi_null_cfi / fi$df_null) - (chi_cfi / fi$df)) /
                     ((chi_null_cfi / fi$df_null) - 1),
                   tolerance = 1e-8)

      # ECVI is still built on the Bartlett-corrected model chi-square.
      n_params <- p * (p + 1) / 2 - fi$df
      expect_equal(fi$ECVI, (fi$chi + 2 * n_params) / (fx$N - 1), tolerance = 1e-8)
    }
  }
})

# The PAF variant code paths reproduce the same psych::fa(fm = "pa") solution as the
# already-tested "psych" type. The EFAtools/SPSS variants (absolute eigenvalues, the
# max-individual convergence criterion) hit the fixed point to ~1e-15; the alternate
# initial-communality paths (mac/unity) reach it within the PAF stopping tolerance.
test_that("PAF variant paths reproduce psych::fa(fm = 'pa')", {
  skip_on_cran()
  skip_if_not_installed("psych")

  R <- test_models$baseline$cormat
  N <- 500
  k <- 3
  ref <- suppressMessages(suppressWarnings(
    psych::fa(R, nfactors = k, n.obs = N, fm = "pa",
              rotate = "none", SMC = TRUE, max.iter = 50)))

  variants <- list(
    list(label = "EFAtools", args = list(type = "EFAtools"),                tol = 1e-4),
    list(label = "SPSS",     args = list(type = "SPSS"),                     tol = 1e-4),
    list(label = "mac",      args = list(type = "psych", init_comm = "mac"),   tol = 5e-3),
    list(label = "unity",    args = list(type = "psych", init_comm = "unity"), tol = 5e-3)
  )

  for (v in variants) {
    efa <- suppressWarnings(do.call(
      EFA, c(list(R, n_factors = k, N = N, method = "PAF"), v$args)))
    expect_equal(unname(1 - efa$h2), unname(1 - ref$communality),
                 tolerance = v$tol, info = v$label)
    expect_equal(repro_offdiag(efa$unrot_loadings),
                 repro_offdiag(ref$loadings), tolerance = v$tol, info = v$label)
  }
})

# AIC and BIC follow the documented chi-square forms, not lavaan's -2logL-based versions.
# CFI and RMSEA place the model (and, for CFI, the baseline) noncentrality on a common
# (N - 1) scale, matching lavaan once its statistics are rescaled from the multiplier N to
# N - 1. The Bartlett correction enters only the reported chi-square test (and the AIC/BIC
# derived from it), not these approximation indices.
test_that("CFI, RMSEA, AIC, and BIC match their references on the baseline", {
  skip_on_cran()

  R <- test_models$baseline$cormat
  N <- 500
  k <- 3
  fi <- suppressWarnings(EFA(R, n_factors = k, N = N, method = "ML"))$fit_indices

  expect_equal(fi$AIC, fi$chi - 2 * fi$df, tolerance = 1e-8)
  expect_equal(fi$BIC, fi$chi - log(N) * fi$df, tolerance = 1e-8)

  skip_if_not_installed("lavaan")

  lav <- tryCatch({
    fit <- lavaan::efa(sample.cov = R, sample.nobs = N, nfactors = k,
                       rotation = "none", estimator = "ML", se = "none")[[1]]
    lavaan::fitMeasures(fit, c("chisq", "df", "baseline.chisq", "baseline.df"))
  }, error = function(e) NULL)
  skip_if(is.null(lav), "lavaan::efa() did not return the baseline fit measures")

  df_m <- unname(lav["df"])
  df_0 <- unname(lav["baseline.df"])
  # lavaan's statistics use the multiplier N; CFI and RMSEA need them on the common (N - 1)
  # noncentrality scale.
  chi_m_n1 <- unname(lav["chisq"])          * (N - 1) / N
  chi_0_n1 <- unname(lav["baseline.chisq"]) * (N - 1) / N

  d_m  <- max(chi_m_n1 - df_m, 0)
  d_0  <- max(chi_0_n1 - df_0, 0)
  cfi_ref   <- 1 - d_m / max(d_m, d_0)
  rmsea_ref <- sqrt(max(0, chi_m_n1 - df_m) / (df_m * (N - 1)))

  # observed agreement on the baseline model: ~1e-13
  expect_equal(fi$CFI, cfi_ref, tolerance = 1e-3)
  expect_equal(fi$RMSEA, rmsea_ref, tolerance = 1e-3)
})

# DWLS reproduces lavaan's diagonally weighted least squares loadings on ordinal data.
# The comparison uses the rotation-invariant reproduced off-diagonal correlations and the
# communalities. Agreement is limited by how closely the two polychoric matrices agree
# (~2e-5), so the observed loading agreement is ~3e-5; the tolerance is set generously at
# 1e-3. The reference is estimator = "DWLS" specifically (not ULS/ULSMV).
test_that("DWLS reproduces lavaan::efa(estimator = 'DWLS') on DOSPERT_raw", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  x <- DOSPERT_raw[stats::complete.cases(DOSPERT_raw), , drop = FALSE]
  k <- 3
  efa <- suppressWarnings(EFA(x, n_factors = k, method = "DWLS",
                              cor_method = "poly", rotation = "none"))
  fit <- suppressWarnings(lavaan::efa(data = as.data.frame(x), nfactors = k,
                                      ordered = colnames(x), estimator = "DWLS",
                                      rotation = "none")[[1]])
  L_ref <- lavaan::lavInspect(fit, "est")$lambda

  expect_equal(repro_offdiag(efa$unrot_loadings), repro_offdiag(L_ref), tolerance = 1e-3)
  expect_equal(unname(efa$h2), unname(diag(tcrossprod(L_ref))), tolerance = 1e-3)
})

test_that("DWLS reproduces lavaan::efa(estimator = 'DWLS') on UPPS_raw", {
  skip_on_cran()
  skip_if_not_slow()
  skip_if_not_installed("lavaan")

  x <- UPPS_raw[stats::complete.cases(UPPS_raw), , drop = FALSE]
  k <- 5
  efa <- suppressWarnings(EFA(x, n_factors = k, method = "DWLS",
                              cor_method = "poly", rotation = "none"))
  fit <- suppressWarnings(lavaan::efa(data = as.data.frame(x), nfactors = k,
                                      ordered = colnames(x), estimator = "DWLS",
                                      rotation = "none")[[1]])
  L_ref <- lavaan::lavInspect(fit, "est")$lambda

  expect_equal(repro_offdiag(efa$unrot_loadings), repro_offdiag(L_ref), tolerance = 1e-3)
  expect_equal(unname(efa$h2), unname(diag(tcrossprod(L_ref))), tolerance = 1e-3)
})

rm(repro_offdiag, reg_fixtures)
