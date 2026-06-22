# MI2S (multiple imputation, two-stage) pooled-inputs route for se = "sandwich"
# in EFA_POOLED(). Pools the correlation matrix and the asymptotic covariance of
# its off-diagonal entries across imputations (Chung & Cai 2019; Sriutaisuk et
# al. 2025), then fits one model on the pooled inputs.

# ---- deterministic fixtures -------------------------------------------------

# m ordinal imputations from a one-factor model, discretised at the quartiles.
.mi2s_ord_imps <- function(m, n = 300L, p = 6L, base_seed = 100L,
                           identical = FALSE) {
  mk <- function(seed) {
    set.seed(seed)
    f <- stats::rnorm(n)
    X <- outer(f, rep(0.6, p)) + matrix(stats::rnorm(n * p, sd = 0.85), n, p)
    d <- apply(X, 2L, function(col) {
      as.integer(cut(col, c(-Inf, stats::quantile(col, c(.25, .5, .75)), Inf)))
    })
    d <- as.data.frame(d)
    names(d) <- paste0("x", seq_len(p))
    d
  }
  if (identical) {
    rep(list(mk(base_seed)), m)
  } else {
    lapply(seq_len(m), function(d) mk(base_seed + d))
  }
}

# m continuous imputations from a one-factor model.
.mi2s_cont_imps <- function(m, n = 400L, p = 6L, base_seed = 200L) {
  lapply(seq_len(m), function(d) {
    set.seed(base_seed + d)
    f <- stats::rnorm(n)
    X <- outer(f, rep(0.6, p)) + matrix(stats::rnorm(n * p, sd = 0.8), n, p)
    out <- as.data.frame(X)
    names(out) <- paste0("v", seq_len(p))
    out
  })
}

# Minimal fit-shaped mocks for the pure pooling helper: a valid correlation
# matrix and a PSD asymptotic covariance over the off-diagonal pairs.
.mi2s_mock_fit <- function(seed, p = 4L, se = "sandwich", cor_method = "poly") {
  set.seed(seed)
  q <- p * (p - 1L) / 2L
  R <- stats::cor(matrix(stats::rnorm(200L * p), 200L, p))
  G <- crossprod(matrix(stats::rnorm(60L * q), 60L, q)) / 60L
  list(orig_R = R, Gamma = G,
       settings = list(se = se, cor_method = cor_method, n_factors = 1L,
                       use = "pairwise.complete.obs", type = "EFAtools"))
}

# Run the MI2S route quietly (the m guidance / "not a correlation matrix"
# messages are not under test here).
.mi2s_fit <- function(imps, ...) {
  suppressMessages(suppressWarnings(
    EFA_POOLED(imps, se = "sandwich", ...)
  ))
}


test_that("identity imputations reproduce the single-fit sandwich solution", {
  imps <- .mi2s_ord_imps(m = 4L, identical = TRUE)

  pooled <- .mi2s_fit(imps, n_factors = 1L, method = "ULS", rotation = "none",
                      cor_method = "poly")
  single <- suppressMessages(
    EFA(imps[[1L]], n_factors = 1L, method = "ULS", rotation = "none",
        se = "sandwich", cor_method = "poly")
  )

  # With identical imputations Gamma_B = 0 and r_bar = r_1, so the pooled fit is
  # the single-imputation fit element-wise.
  expect_equal(unclass(pooled$unrot_loadings), unclass(single$unrot_loadings),
               tolerance = 1e-10)
  expect_equal(pooled$SE$unrot_loadings, single$SE$unrot_loadings,
               tolerance = 1e-10)
  expect_equal(pooled$SE$uniquenesses, single$SE$uniquenesses,
               tolerance = 1e-10)
  expect_equal(pooled$fit_indices$chi, single$fit_indices$chi,
               tolerance = 1e-10)
  expect_equal(pooled$fit_indices$chi_shift, single$fit_indices$chi_shift,
               tolerance = 1e-10)
})


test_that("Gamma_tilde equals the Rubin total ACOV to machine precision", {
  fits <- lapply(seq_len(6L), .mi2s_mock_fit, p = 4L)
  out <- .efa_pooled_mi2s_inputs(fits)

  m <- length(fits)
  p <- ncol(fits[[1L]]$orig_R)
  idx <- utils::combn(p, 2L)
  pair <- cbind(idx[1L, ], idx[2L, ])
  R_list <- lapply(fits, `[[`, "orig_R")
  G_list <- lapply(fits, `[[`, "Gamma")

  V <- t(vapply(R_list, function(R) R[pair], numeric(ncol(idx))))
  Gamma_W <- Reduce(`+`, G_list) / m
  # Recompute the between-imputation ACOV from first principles (centred
  # cross-product with the m - 1 divisor) rather than via stats::cov(), so the
  # m - 1 divisor and the (1 + 1/m) factor are checked independently of the
  # implementation.
  Vc <- sweep(V, 2L, colMeans(V))
  Gamma_B <- crossprod(Vc) / (m - 1L)
  Gamma_ref <- Gamma_W + (1 + 1 / m) * Gamma_B
  Gamma_ref <- (Gamma_ref + t(Gamma_ref)) / 2

  expect_equal(unname(out$Gamma_tilde), unname(Gamma_ref), tolerance = 1e-12)

  r_ref <- Reduce(`+`, R_list) / m
  expect_equal(out$r_bar[lower.tri(out$r_bar)], r_ref[lower.tri(r_ref)],
               tolerance = 1e-12)
  # The off-diagonal stacking follows combn(p, 2): pair k is (idx[1,k], idx[2,k]).
  expect_identical(dim(out$Gamma_tilde), c(ncol(idx), ncol(idx)))
})


test_that("a non-PSD pooled ACOV aborts instead of being projected", {
  fits <- lapply(seq_len(4L), .mi2s_mock_fit, p = 4L)
  q <- ncol(fits[[1L]]$Gamma)

  # PSD inputs pool to a PSD Gamma_tilde -> no error.
  expect_silent(invisible(.efa_pooled_mi2s_inputs(fits)))

  # Inject a clearly indefinite asymptotic covariance into one imputation.
  bad <- diag(q)
  bad[1L, 1L] <- -5
  fits[[1L]]$Gamma <- bad
  expect_error(.efa_pooled_mi2s_inputs(fits),
               class = "efa_pooled_mi2s_acov_not_psd")
})


test_that("MI2S matches the lavaan.mi poolSat input-pooling oracle", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("lavaan.mi")

  m <- 10L
  imps_num <- .mi2s_ord_imps(m = m, n = 400L, p = 6L, base_seed = 1000L)
  vn <- paste0("x", seq_len(6L))
  imps_ord <- lapply(imps_num, function(d) {
    for (j in names(d)) d[[j]] <- ordered(d[[j]])
    d
  })

  pooled <- .mi2s_fit(imps_num, n_factors = 1L, method = "DWLS",
                      rotation = "none", cor_method = "poly")

  # poolSat() pools the inputs the same way (scale.W = FALSE gives exactly the
  # Rubin total ACOV), so a single cfa() on its lavMoments is the matched oracle.
  # lavaan.mi advises attaching itself; that note is irrelevant to the oracle.
  model <- paste("F =~", paste(vn, collapse = " + "))
  lfit <- suppressWarnings(suppressMessages({
    pre <- lavaan.mi::poolSat(imps_ord, ordered = vn, scale.W = FALSE)
    lavaan::cfa(model, data = pre, std.lv = TRUE)
  }))
  pe <- lavaan::parameterEstimates(lfit)
  pe <- pe[pe$op == "=~", ]
  ss <- lavaan::lavInspect(lfit, "test")[["scaled.shifted"]]

  sc <- pre$sample.cov
  if (is.list(sc)) sc <- sc[[1L]]
  rbar <- pooled$orig_R
  expect_equal(rbar[lower.tri(rbar)], sc[lower.tri(sc)], tolerance = 1e-4)

  L_mine <- as.numeric(pooled$unrot_loadings)
  sgn <- sign(sum(L_mine)) * sign(sum(pe$est))
  expect_equal(L_mine, pe$est * sgn, tolerance = 5e-3)
  expect_equal(as.numeric(pooled$SE$unrot_loadings), pe$se, tolerance = 0.02)

  # Native scaled-shifted statistic: df exact, statistic/shift within the
  # polychoric-optimiser gap; lavaan stores the reciprocal of the multiplier.
  expect_equal(pooled$fit_indices$df, ss$df)
  expect_equal(pooled$fit_indices$chi, ss$stat, tolerance = 0.02)
  expect_equal(pooled$fit_indices$chi_shift, ss$shift.parameter, tolerance = 0.02)
  expect_equal(1 / pooled$fit_indices$chi_scaling, ss$scaling.factor,
               tolerance = 0.02)
})


test_that("the pooled fit carries a single native scaled-shifted chi-square", {
  imps <- .mi2s_ord_imps(m = 8L)
  pooled <- .mi2s_fit(imps, n_factors = 1L, method = "DWLS", rotation = "none",
                      cor_method = "poly")
  fi <- pooled$fit_indices

  expect_identical(fi$chi_scaled_type, "scaled.shifted")
  # Scalars from the single fit, not D2 pools of m statistics.
  expect_length(fi$chi_scaling, 1L)
  expect_length(fi$chi_shift, 1L)
  expect_true(is.finite(fi$chi_scaling) && is.finite(fi$chi_shift))
  # Likelihood-ratio information criteria have no meaning on the scaled statistic.
  expect_true(all(is.na(c(fi$AIC, fi$BIC, fi$ECVI))))
})


test_that("the MI2S path invokes no rotation-alignment machinery", {
  imps <- .mi2s_ord_imps(m = 6L)

  called <- new.env(parent = emptyenv())
  called$any <- FALSE
  spy <- function(...) {
    called$any <- TRUE
    NULL
  }
  testthat::local_mocked_bindings(
    PROCRUSTES = spy,
    .gpa_consensus_target = spy,
    .align_solution = spy,
    .efa_pooled_rubin_core = spy
  )

  pooled <- .mi2s_fit(imps, n_factors = 1L, method = "ULS", rotation = "oblimin",
                      cor_method = "poly")

  expect_false(called$any)
  expect_s3_class(pooled, "EFA_POOLED")
})


test_that("the pooled object exposes mi_fit alongside the component fits", {
  imps <- .mi2s_ord_imps(m = 5L)
  pooled <- .mi2s_fit(imps, n_factors = 2L, method = "DWLS", rotation = "oblimin",
                      cor_method = "poly")

  expect_s3_class(pooled, "EFA_POOLED")
  expect_s3_class(pooled$mi_fit, "EFA")
  expect_length(pooled$fits, 5L)

  # Pooled point estimates and fit indices come straight from the single fit.
  expect_equal(pooled$orig_R, pooled$mi_fit$orig_R)
  expect_identical(pooled$fit_indices, pooled$mi_fit$fit_indices)
  expect_identical(pooled$SE, pooled$mi_fit$SE)
  expect_identical(pooled$CI, pooled$mi_fit$CI)

  # Oblique sandwich SE schema flows through to the pooled object.
  expect_setequal(names(pooled$SE),
                  c("unrot_loadings", "uniquenesses", "rot_loadings",
                    "communalities", "Phi", "Structure"))
  # The rotated sandwich SE values must be real numbers, not an all-NA regression
  # (a broken oblique Jacobian or unreliable bordered information would NA them
  # while leaving the schema names intact).
  expect_true(all(is.finite(pooled$SE$rot_loadings)))
  expect_true(all(pooled$SE$rot_loadings > 0))
  expect_true(all(is.finite(pooled$SE$Phi)))
  expect_true(all(is.finite(pooled$SE$Structure)))
  expect_true(all(is.finite(pooled$SE$communalities)))
  # No per-parameter Rubin pooling on this path.
  expect_null(pooled$MI)
  expect_null(pooled$replicates)
  expect_true(isTRUE(pooled$settings$pooled))
  expect_identical(pooled$settings$n_imputations, 5L)
})


test_that("continuous Pearson imputations run through MI2S (ML and ULS)", {
  imps <- .mi2s_cont_imps(m = 6L)

  for (meth in c("ML", "ULS")) {
    pooled <- .mi2s_fit(imps, n_factors = 1L, method = meth, rotation = "none",
                        cor_method = "pearson")
    expect_s3_class(pooled, "EFA_POOLED")
    expect_true(all(is.finite(pooled$SE$unrot_loadings)), info = meth)
    expect_true(all(pooled$SE$unrot_loadings > 0), info = meth)
    expect_identical(pooled$fit_indices$chi_scaled_type, "scaled.shifted")
    expect_s3_class(pooled$mi_fit, "EFA")
  }
})


# Drive .efa_pooled_mi2s() directly with hand-built fit mocks, supplying a
# matching data_list of correlation matrices (so .efa_pooled_get_Ns falls back to
# NA rather than nrow()) unless raw N is provided some other way.
.mi2s_call <- function(fits, efa_args = list()) {
  .efa_pooled_mi2s(
    fits = fits, data_list = lapply(fits, `[[`, "orig_R"), efa_args = efa_args,
    settings = fits[[1L]]$settings, method = "ULS", rotation = "none",
    rotation_type = "none", target_method = "first_target",
    align_unrotated = "signed_tucker_congruence", fit_pool_method = "D2",
    p = 0.05, rmsea_ci_level = 0.90, rmsr_upper = TRUE
  )
}

test_that("inconsistent component inputs fail closed", {
  # Mixed se: one sandwich, one information.
  expect_error(
    .mi2s_call(list(.mi2s_mock_fit(1L, se = "sandwich"),
                    .mi2s_mock_fit(2L, se = "information"))),
    class = "efa_pooled_mi2s_inputs_inconsistent"
  )

  # Mixed cor_method across imputations.
  expect_error(
    .mi2s_call(list(.mi2s_mock_fit(3L, cor_method = "poly"),
                    .mi2s_mock_fit(4L, cor_method = "pearson"))),
    class = "efa_pooled_mi2s_inputs_inconsistent"
  )

  # A consistent but unsupported cor_method (no sandwich asymptotic covariance).
  expect_error(
    .mi2s_call(list(.mi2s_mock_fit(5L, cor_method = "spearman"),
                    .mi2s_mock_fit(6L, cor_method = "spearman"))),
    class = "efa_pooled_mi2s_inputs_inconsistent"
  )

  # A missing per-fit asymptotic covariance.
  no_gamma <- list(.mi2s_mock_fit(7L), .mi2s_mock_fit(8L))
  no_gamma[[2L]]$Gamma <- NULL
  expect_error(.mi2s_call(no_gamma),
               class = "efa_pooled_mi2s_inputs_inconsistent")
})


test_that("MI2S aborts when no imputation supplies N", {
  # Consistent sandwich/poly mocks with valid Gamma but no recoverable N: the
  # data_list holds correlation matrices and efa_args carries no N.
  fits <- lapply(seq_len(3L), .mi2s_mock_fit)
  expect_error(.mi2s_call(fits, efa_args = list()),
               class = "efa_pooled_mi2s_no_n")
})


test_that("MI2S warns on too few imputations and on ignored alignment settings", {
  imps <- .mi2s_ord_imps(m = 5L)

  # Fewer than 20 imputations: the scaled-shifted statistic is under-calibrated.
  expect_warning(
    suppressMessages(EFA_POOLED(imps, n_factors = 1L, method = "ULS",
                                rotation = "none", se = "sandwich",
                                cor_method = "poly")),
    class = "efa_pooled_mi2s_n_too_small"
  )

  # A non-default alignment is inert on the single-fit MI2S path. The m guidance
  # warning also fires here; muffle it so only the alignment note is asserted.
  expect_warning(
    withCallingHandlers(
      suppressMessages(EFA_POOLED(imps, n_factors = 1L, method = "ULS",
                                  rotation = "none", se = "sandwich",
                                  cor_method = "poly",
                                  align_unrotated = "procrustes")),
      efa_pooled_mi2s_n_too_small = function(w) invokeRestart("muffleWarning")
    ),
    class = "efa_pooled_mi2s_alignment_ignored"
  )
})


test_that("tetrachoric (binary) imputations run through MI2S end to end", {
  skip_on_cran()
  # Dichotomise the ordinal generator at the median for a binary fixture.
  imps <- lapply(.mi2s_ord_imps(m = 6L, p = 6L, base_seed = 300L), function(d) {
    as.data.frame(lapply(d, function(col) as.integer(col > stats::median(col))))
  })

  pooled <- .mi2s_fit(imps, n_factors = 1L, method = "DWLS", rotation = "none",
                      cor_method = "tetra")

  expect_s3_class(pooled, "EFA_POOLED")
  expect_true(all(is.finite(pooled$SE$unrot_loadings)))
  expect_true(all(pooled$SE$unrot_loadings > 0))
  expect_identical(pooled$fit_indices$chi_scaled_type, "scaled.shifted")
  expect_identical(pooled$settings$cor_method, "tetra")
})


test_that("MI2S print omits inert pooling settings and labels the CI source", {
  imps <- .mi2s_ord_imps(m = 6L)
  pooled <- .mi2s_fit(imps, n_factors = 1L, method = "DWLS", rotation = "none",
                      cor_method = "poly")

  header <- cli::ansi_strip(format(pooled))
  # The single-fit MI2S path has no per-imputation alignment / chi pooling, so
  # those settings must not appear in the pooled header.
  expect_false(any(grepl("align_unrotated|target_method|fit_pool_method", header)))

  # The full summary carries the CI-provenance note; it must describe the robust
  # sandwich covariance, not the expected information matrix.
  body <- cli::ansi_strip(format(summary(pooled)))
  note <- body[grepl("Note:", body)]
  expect_true(any(grepl("sandwich", note)))
  expect_false(any(grepl("expected information matrix", note)))
})
