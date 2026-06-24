# Rubin's-rules pooling of analytic unrotated-loading and uniqueness SEs in
# EFA_POOLED(). Covers the se = "information" path, where component fits expose
# closed-form expected-information SE matrices and no bootstrap replicate cubes.
#
# Assertion groups:
#   1. Identity-imputation invariance.
#   2. Signed-permutation invariance of marginal SEs after alignment.
#   3. lavaan.mi cross-check (SE only; df conventions differ).
#   4. Rubin (1987) df sanity.
#   5. Uniquenesses pooling shape and invariance.
#   6. Wald CI symmetry and monotonicity.
#   7. FMI / RIV bounds and per-element NA propagation (fail-closed).
#   8. Slot shape, settings preservation, and Procrustes abort.

cormat <- test_models$baseline$cormat
p_vars <- ncol(cormat)
N_id   <- 500L
k      <- 3L
m_id   <- 5L
identical_list <- replicate(m_id, cormat, simplify = FALSE)

# ---- Group 1 ---------------------------------------------------------------

test_that("identity imputations: pooled SE equals the per-fit SE and FMI is zero", {
  pooled <- suppressMessages(EFA_POOLED(
    identical_list, n_factors = k, N = N_id, method = "ML",
    rotation = "none", se = "information"
  ))

  oracle <- EFA(cormat, n_factors = k, N = N_id, method = "ML",
                rotation = "none", se = "information")

  expect_identical(pooled$settings$se, "information")
  expect_identical(pooled$settings$component_se, "information")

  expect_equal(pooled$SE$unrot_loadings, oracle$SE$unrot_loadings, tolerance = 1e-10)
  expect_equal(pooled$SE$uniquenesses, oracle$SE$uniquenesses, tolerance = 1e-10)
  expect_equal(unclass(pooled$unrot_loadings), unclass(oracle$unrot_loadings),
               tolerance = 1e-10)

  fmi <- pooled$MI$unrot_loadings$FMI
  riv <- pooled$MI$unrot_loadings$RIV
  expect_true(all(fmi[is.finite(fmi)] == 0))
  expect_true(all(riv[is.finite(riv)] == 0))

  # Plain Rubin df (nu_com = Inf): identical imputations leave B = 0, so the
  # relative increase in variance r = 0 and nu_old = (m - 1)(1 + 1/r)^2 -> Inf.
  # Every element's df is therefore Inf (the Wald interval reduces to the normal
  # reference), matching lavaan.mi's df -> Inf for asymptotically-normal estimates.
  df <- pooled$MI$unrot_loadings$df
  expect_true(all(is.infinite(df)))
})

# ---- Group 2 ---------------------------------------------------------------

test_that("aligned-pool SEs are invariant to per-imputation sign flips and column permutations", {
  base_fit <- EFA(cormat, n_factors = k, N = N_id, method = "ML",
                  rotation = "none", se = "information")
  L_base  <- unclass(base_fit$unrot_loadings)
  SE_base <- base_fit$SE$unrot_loadings
  psi_se  <- base_fit$SE$uniquenesses

  # The first entry is the anchor (identity gauge); the rest are deliberate
  # column-permutation + sign-flip variants of the same fit. Marginal SEs are
  # non-negative scalars (sign-flip is a no-op) and merely permute under
  # column reordering, so the analytic pool must recover the anchor SE matrix
  # exactly (B = 0 within the identity-gauge orbit).
  gauges <- list(
    list(perm = c(1L, 2L, 3L), signs = c( 1,  1,  1)),
    list(perm = c(2L, 1L, 3L), signs = c(-1,  1,  1)),
    list(perm = c(3L, 2L, 1L), signs = c( 1, -1,  1)),
    list(perm = c(1L, 3L, 2L), signs = c( 1,  1, -1)),
    list(perm = c(2L, 3L, 1L), signs = c(-1, -1,  1))
  )

  permuted_fits <- lapply(gauges, function(g) {
    fit <- base_fit
    L_g  <- L_base[, g$perm, drop = FALSE]
    L_g  <- sweep(L_g, 2, g$signs, `*`)
    fit$unrot_loadings <- structure(L_g,
      class    = class(base_fit$unrot_loadings),
      dimnames = dimnames(base_fit$unrot_loadings)
    )
    fit$SE$unrot_loadings <- SE_base[, g$perm, drop = FALSE]
    fit$SE$uniquenesses   <- psi_se
    fit
  })

  unrot_list <- lapply(permuted_fits, function(f) unclass(f$unrot_loadings))
  aligned <- EFAtools:::.efa_pooled_align_unrotated_list(
    unrot_list,
    align_unrotated = "signed_tucker_congruence",
    return_meta = TRUE
  )

  pool <- EFAtools:::.efa_pooled_analytic_pool(
    fits = permuted_fits,
    unrot_loadings_aligned = aligned$loadings,
    align_meta = aligned$meta,
    ci = 0.95
  )

  expect_equal(unname(pool$SE$unrot_loadings), unname(SE_base), tolerance = 1e-12)
  expect_equal(unname(pool$SE$uniquenesses),   unname(psi_se),  tolerance = 1e-12)
})

# ---- Group 3 ---------------------------------------------------------------

test_that("Rubin pooling of marginal SEs matches lavaan.mi when fed identical per-imp inputs", {
  # This test isolates the pooling math from the fit math. EFAtools' ML
  # extraction and lavaan's ML estimator use different identifications and
  # produce per-imputation loadings that can drift by ~5% in finite samples,
  # so feeding each pipeline its own per-imp fits would mix that drift into
  # the comparison. Instead we fit lavaan per imputation, capture est + se,
  # pool those inputs through .efa_pooled_rubin_core, and compare to
  # lavaan.mi's pooled SE. Pooled SEs are df-agnostic and must agree to ~1e-6.
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("lavaan.mi")
  skip_if_not_installed("MASS")

  set.seed(20260621)
  p_lav <- 6L
  N_lav <- 400L
  m_lav <- 10L
  L_true <- matrix(c(0.70, 0.80, 0.65, 0.55, 0.75, 0.60), p_lav, 1)
  Sigma  <- tcrossprod(L_true)
  diag(Sigma) <- 1

  imps <- lapply(seq_len(m_lav), function(i) {
    X <- MASS::mvrnorm(N_lav, mu = rep(0, p_lav), Sigma = Sigma)
    colnames(X) <- paste0("V", seq_len(p_lav))
    as.data.frame(X)
  })

  vars  <- colnames(imps[[1]])
  model <- paste0("f =~ ", paste(vars, collapse = " + "))

  # Per-imp single-fit estimates and SEs from lavaan, in canonical loading order.
  est_lav <- matrix(NA_real_, m_lav, p_lav)
  se_lav  <- matrix(NA_real_, m_lav, p_lav)
  for (d in seq_len(m_lav)) {
    fit_d <- suppressWarnings(suppressMessages(
      lavaan::cfa(model, data = imps[[d]], std.lv = TRUE)
    ))
    pe <- lavaan::parameterEstimates(fit_d)
    ld <- pe[pe$op == "=~" & pe$lhs == "f", ]
    ld <- ld[match(vars, ld$rhs), ]
    est_lav[d, ] <- ld$est
    se_lav[d, ]  <- ld$se
  }
  if (anyNA(est_lav) || anyNA(se_lav)) skip("lavaan per-imp fits returned NAs")

  Qbar <- colMeans(est_lav)
  Ubar <- colMeans(se_lav^2)
  B    <- apply(est_lav, 2, stats::var)

  efa_core <- EFAtools:::.efa_pooled_rubin_core(
    Qbar, Ubar, B, m = m_lav, N_pool = Inf, alpha = 0.05
  )

  fit_mi <- tryCatch(
    suppressWarnings(suppressMessages(
      lavaan.mi::cfa.mi(model, data = imps, std.lv = TRUE)
    )),
    error = function(e) e
  )
  if (inherits(fit_mi, "error")) {
    skip(paste0("lavaan.mi::cfa.mi failed: ", conditionMessage(fit_mi)))
  }

  # asymptotic = TRUE selects the textbook Rubin (1987) formula
  # T = W + (1 + 1/m) B with normal-z CIs. lavaan.mi's default
  # (asymptotic = FALSE) instead rescales W by (1 + 1/m), which is a
  # lavaan-specific adjustment and not what EFA_POOLED implements.
  pe_mi <- tryCatch(
    lavaan.mi::parameterEstimates.mi(fit_mi, ci = TRUE, level = 0.95,
                                     asymptotic = TRUE),
    error = function(e) NULL
  )
  if (is.null(pe_mi)) skip("lavaan.mi::parameterEstimates.mi failed")
  ld_mi <- pe_mi[pe_mi$op == "=~" & pe_mi$lhs == "f", ]
  ld_mi <- ld_mi[match(vars, ld_mi$rhs), ]
  if (anyNA(ld_mi$se)) skip("lavaan.mi returned NA standard errors")

  # Same per-imp inputs -> same pooled SE under the textbook Rubin formula.
  expect_equal(efa_core$se, ld_mi$se, tolerance = 1e-6)
})

# ---- Group 4 ---------------------------------------------------------------

test_that(".efa_pooled_rubin_core: BR df inflates the critical, is monotone in N_pool, and limits to Rubin", {
  Qbar <- 0.5
  Ubar <- 0.01
  B    <- 0.005
  m    <- 5L

  core_500 <- EFAtools:::.efa_pooled_rubin_core(Qbar, Ubar, B, m,
                                                N_pool = 500, alpha = 0.05)
  core_50  <- EFAtools:::.efa_pooled_rubin_core(Qbar, Ubar, B, m,
                                                N_pool = 50,  alpha = 0.05)
  core_inf <- EFAtools:::.efa_pooled_rubin_core(Qbar, Ubar, B, m,
                                                N_pool = Inf, alpha = 0.05)

  expect_true(is.finite(core_500$df))
  expect_gt(stats::qt(0.975, core_500$df), stats::qnorm(0.975))

  rubin_df <- (m - 1) * (1 + Ubar / ((1 + 1 / m) * B))^2
  expect_equal(core_inf$df, rubin_df, tolerance = 1e-9)

  # Shrinking N_pool shrinks the BR df.
  expect_lt(core_50$df, core_500$df)
  expect_lt(core_500$df, core_inf$df)

  # SE is df-agnostic; identical across all three calls.
  expect_equal(core_50$se,  core_inf$se, tolerance = 1e-12)
  expect_equal(core_500$se, core_inf$se, tolerance = 1e-12)
})

test_that(".efa_pooled_rubin_core: BR fmi = 1 (Ubar = 0, B > 0) falls back to plain-Rubin df", {
  # Within-imputation analytic SE = 0 with non-zero between-imputation variance
  # makes fmi -> 1 and nu_obs -> 0 in the BR formula, which would otherwise
  # collapse df to 0 and qt(., 0) to NaN. The core falls back to df_old.
  core <- EFAtools:::.efa_pooled_rubin_core(
    Qbar = 0.5, Ubar = 0, B = 0.01, m = 5L, N_pool = 200, alpha = 0.05
  )
  expect_true(is.finite(core$df))
  expect_equal(core$df, 5L - 1L, tolerance = 1e-12)
  expect_true(is.finite(core$ci$lower))
  expect_true(is.finite(core$ci$upper))
  expect_equal(core$FMI, 1, tolerance = 1e-12)
})

# ---- Group 5 ---------------------------------------------------------------

test_that("uniqueness pooling has correct length, is non-negative and identity-invariant", {
  pooled <- suppressMessages(EFA_POOLED(
    identical_list, n_factors = k, N = N_id, method = "ML",
    rotation = "none", se = "information"
  ))
  oracle <- EFA(cormat, n_factors = k, N = N_id, method = "ML",
                rotation = "none", se = "information")

  expect_length(pooled$SE$uniquenesses, p_vars)
  expect_true(all(pooled$SE$uniquenesses[!is.na(pooled$SE$uniquenesses)] >= 0))
  expect_equal(unname(pooled$SE$uniquenesses), unname(oracle$SE$uniquenesses),
               tolerance = 1e-10)
})

# ---- Group 6 ---------------------------------------------------------------

test_that("pooled Wald CIs are symmetric around the pooled estimate", {
  set.seed(42)
  imps <- lapply(seq_len(5), function(i) {
    GRiPS_raw[sample.int(nrow(GRiPS_raw), replace = TRUE), , drop = FALSE]
  })
  pooled <- suppressWarnings(suppressMessages(EFA_POOLED(
    imps, n_factors = 1, method = "ML", rotation = "none", se = "information"
  )))

  est <- unclass(pooled$unrot_loadings)
  lo  <- pooled$CI$unrot_loadings$lower
  hi  <- pooled$CI$unrot_loadings$upper

  finite <- is.finite(lo) & is.finite(hi) & is.finite(est)
  expect_true(any(finite))
  expect_true(all(lo[finite] < est[finite]))
  expect_true(all(est[finite] < hi[finite]))
  expect_equal(unname((lo[finite] + hi[finite]) / 2),
               unname(est[finite]), tolerance = 1e-10)
})

# ---- Group 7 ---------------------------------------------------------------

test_that("FMI is in [0,1], RIV >= 0; per-element NA in any imputation propagates", {
  set.seed(42)
  imps <- lapply(seq_len(5), function(i) {
    GRiPS_raw[sample.int(nrow(GRiPS_raw), replace = TRUE), , drop = FALSE]
  })
  pooled <- suppressWarnings(suppressMessages(EFA_POOLED(
    imps, n_factors = 1, method = "ML", rotation = "none", se = "information"
  )))

  fmi <- pooled$MI$unrot_loadings$FMI
  riv <- pooled$MI$unrot_loadings$RIV
  expect_true(all(fmi[is.finite(fmi)] >= 0))
  expect_true(all(fmi[is.finite(fmi)] <= 1))
  expect_true(all(riv[is.finite(riv)] >= 0))

  # Inject NA into one imputation's SE for a single element and re-pool via
  # the helper directly. Expect NA at that position in SE, both CI bounds, and
  # all three MI fields; surrounding elements must remain finite.
  fits <- pooled$fits
  fits[[3]]$SE$unrot_loadings[2, 1] <- NA_real_

  unrot_list <- lapply(fits, function(f) unclass(f$unrot_loadings))
  aligned <- EFAtools:::.efa_pooled_align_unrotated_list(
    unrot_list, align_unrotated = "signed_tucker_congruence", return_meta = TRUE
  )
  na_pool <- EFAtools:::.efa_pooled_analytic_pool(
    fits = fits,
    unrot_loadings_aligned = aligned$loadings,
    align_meta = aligned$meta,
    ci = 0.95
  )

  expect_true(is.na(na_pool$SE$unrot_loadings[2, 1]))
  expect_true(is.na(na_pool$CI$unrot_loadings$lower[2, 1]))
  expect_true(is.na(na_pool$CI$unrot_loadings$upper[2, 1]))
  expect_true(is.na(na_pool$MI$unrot_loadings$FMI[2, 1]))
  expect_true(is.na(na_pool$MI$unrot_loadings$RIV[2, 1]))
  expect_true(is.na(na_pool$MI$unrot_loadings$df[2, 1]))

  # Other positions in the same column unaffected.
  expect_false(is.na(na_pool$SE$unrot_loadings[1, 1]))
  expect_false(is.na(na_pool$SE$unrot_loadings[3, 1]))
})

# ---- Group 8 ---------------------------------------------------------------

test_that("analytic-pool output has correct slot shapes and preserves se = 'information'", {
  pooled <- suppressMessages(EFA_POOLED(
    identical_list, n_factors = k, N = N_id, method = "ML",
    rotation = "none", se = "information"
  ))

  loading_dn <- dimnames(unclass(pooled$unrot_loadings))

  expect_identical(dim(pooled$SE$unrot_loadings),         c(p_vars, k))
  expect_identical(dimnames(pooled$SE$unrot_loadings),    loading_dn)
  expect_identical(dim(pooled$CI$unrot_loadings$lower),   c(p_vars, k))
  expect_identical(dim(pooled$CI$unrot_loadings$upper),   c(p_vars, k))
  expect_identical(dimnames(pooled$CI$unrot_loadings$lower), loading_dn)

  expect_setequal(names(pooled$MI$unrot_loadings), c("RIV", "FMI", "df"))
  expect_identical(dim(pooled$MI$unrot_loadings$RIV), c(p_vars, k))
  expect_identical(dim(pooled$MI$unrot_loadings$FMI), c(p_vars, k))
  expect_identical(dim(pooled$MI$unrot_loadings$df),  c(p_vars, k))

  expect_length(pooled$SE$uniquenesses,             p_vars)
  expect_length(pooled$CI$uniquenesses$lower,       p_vars)
  expect_length(pooled$CI$uniquenesses$upper,       p_vars)
  expect_setequal(names(pooled$MI$uniquenesses),    c("RIV", "FMI", "df"))
  expect_length(pooled$MI$uniquenesses$RIV,         p_vars)
  expect_length(pooled$MI$uniquenesses$FMI,         p_vars)
  expect_length(pooled$MI$uniquenesses$df,          p_vars)

  # Analytic path: replicates slot is present-but-NULL (matches the EFA() schema
  # contract: `expect_true("replicates" %in% names(fit));
  # expect_null(fit$replicates)`), no rotation-failure counters.
  expect_true("replicates" %in% names(pooled))
  expect_null(pooled$replicates)
  expect_null(pooled$MI$unrot_loadings$bootstrap_rotation_failures)
  expect_null(pooled$MI$unrot_loadings$bootstrap_rotation_valid)
})

test_that("analytic-pooled summary note advertises Wald-from-information CIs, not bootstrap", {
  pooled <- suppressMessages(EFA_POOLED(
    identical_list, n_factors = k, N = N_id, method = "ML",
    rotation = "none", se = "information"
  ))
  # summary() exercises the full view that emits the provenance note; the brief
  # print() view suppresses it (passes ci = "none" to .print_efa_bootstrap_note).
  out <- testthat::capture_output(print(summary(pooled)))
  expect_match(out, "Wald CIs from the expected information matrix",
               fixed = TRUE)
  expect_false(grepl("Bootstrap/MI CIs", out, fixed = TRUE))
})

