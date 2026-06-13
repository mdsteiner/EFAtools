# Tests for EFA_POOLED(): structure of the returned object, pooling math,
# classed conditions, bootstrap/MI pooling, and the print/format methods.

cormat <- test_models$baseline$cormat
p_vars <- ncol(cormat)
cormat_list <- list(cormat, cormat, cormat)

pooled_obl <- EFA_POOLED(cormat_list, n_factors = 3, N = 500, method = "PAF",
                         rotation = "promax")
pooled_orth <- EFA_POOLED(cormat_list, n_factors = 3, N = 500, method = "PAF",
                          rotation = "varimax")

set.seed(42)
grips_list <- lapply(1:3, function(i) {
  GRiPS_raw[sample(seq_len(nrow(GRiPS_raw)), replace = TRUE), ]
})
pooled_none <- suppressMessages(
  EFA_POOLED(grips_list, n_factors = 1, method = "ML")
)

test_that("EFA_POOLED returns a well-formed pooled object", {
  expect_s3_class(pooled_obl, c("EFA_POOLED", "EFA"), exact = TRUE)

  core <- c("h2", "unrot_loadings", "vars_accounted", "fit_indices",
            "model_implied_R", "residuals", "orig_R", "settings", "fits",
            "alignment", "mi_diagnostics", "rot_loadings",
            "vars_accounted_rot", "Phi", "Structure")
  expect_true(all(core %in% names(pooled_obl)))

  expect_type(pooled_obl$h2, "double")
  expect_length(pooled_obl$h2, p_vars)
  expect_named(pooled_obl$h2)

  expect_s3_class(pooled_obl$unrot_loadings, "LOADINGS")
  expect_s3_class(pooled_obl$rot_loadings, "LOADINGS")
  expect_s3_class(pooled_obl$Structure, "LOADINGS")
  expect_identical(dim(unclass(pooled_obl$rot_loadings)), c(p_vars, 3L))

  expect_identical(dim(pooled_obl$Phi), c(3L, 3L))
  expect_equal(pooled_obl$Phi, t(pooled_obl$Phi))

  expect_equal(unname(diag(pooled_obl$model_implied_R)), rep(1, p_vars))
  expect_equal(unname(diag(pooled_obl$residuals)), rep(0, p_vars))

  expect_length(pooled_obl$fits, 3)
  expect_true(all(vapply(pooled_obl$fits, inherits, logical(1), "EFA")))
  expect_true(isTRUE(pooled_obl$alignment$converged))

  expect_type(pooled_obl$fit_indices, "list")
  expect_true(is.finite(pooled_obl$fit_indices$CAF))
  expect_true(is.finite(pooled_obl$fit_indices$RMSR))
})

test_that("pooled CAF reproduces the single-solution CAF on identical imputations", {
  # Pooling identical correlation matrices reproduces the single EFA solution, so
  # the pooled CAF (computed on the residual matrix with a unit diagonal) must
  # equal the CAF that EFA() reports for that solution.
  single <- EFA(cormat, n_factors = 3, N = 500, method = "PAF", rotation = "promax")
  expect_equal(pooled_obl$fit_indices$CAF, single$fit_indices$CAF, tolerance = 1e-6)
  expect_gt(pooled_obl$fit_indices$CAF, 0)
  expect_lt(pooled_obl$fit_indices$CAF, 1)
})

test_that("EFA_POOLED records the pooling settings", {
  s <- pooled_obl$settings
  expect_true(s$pooled)
  expect_identical(s$n_imputations, 3L)
  expect_identical(s$target_method, "consensus")
  expect_identical(s$align_unrotated, "signed_tucker_congruence")
  expect_identical(s$fit_pool_method, "D2")
  expect_equal(s$p, 0.05)
  expect_equal(s$ci, 0.95)
  # no bootstrap arrays were available, so the pooled object must not claim SEs
  expect_identical(s$se, "none")
  expect_false(any(c("boot.SE", "boot.CI", "boot.arrays", "boot.MI") %in%
                     names(pooled_obl)))
})

test_that("rotation variants include exactly the applicable components", {
  # orthogonal: rotated loadings but no factor intercorrelations
  expect_s3_class(pooled_orth$rot_loadings, "LOADINGS")
  expect_true(!is.null(pooled_orth$vars_accounted_rot))
  expect_null(pooled_orth$Phi)
  expect_null(pooled_orth$Structure)

  # unrotated: no rotated components at all
  expect_null(pooled_none$rot_loadings)
  expect_null(pooled_none$vars_accounted_rot)
  expect_null(pooled_none$Phi)
  expect_null(pooled_none$Structure)
  expect_null(pooled_none$alignment)
})

test_that("alignment variants produce well-formed pooled objects", {
  # first_target: align every imputation to the first rotated solution rather
  # than to an iteratively updated consensus target
  ft <- EFA_POOLED(cormat_list, n_factors = 3, N = 500, method = "PAF",
                   rotation = "promax", target_method = "first_target")
  expect_s3_class(ft, c("EFA_POOLED", "EFA"), exact = TRUE)
  expect_identical(ft$settings$target_method, "first_target")
  expect_identical(ft$alignment$method, "first_target")
  expect_s3_class(ft$rot_loadings, "LOADINGS")
  expect_identical(dim(ft$Phi), c(3L, 3L))

  # align_unrotated = "procrustes": orthogonal Procrustes of the unrotated axes
  pr <- EFA_POOLED(cormat_list, n_factors = 3, N = 500, method = "PAF",
                   rotation = "none", align_unrotated = "procrustes")
  expect_identical(pr$settings$align_unrotated, "procrustes")
  expect_s3_class(pr$unrot_loadings, "LOADINGS")

  # align_unrotated = "none": average the unrotated loadings as returned. With
  # identical imputations the pooled result must equal the single fit.
  nn <- EFA_POOLED(cormat_list, n_factors = 3, N = 500, method = "PAF",
                   rotation = "none", align_unrotated = "none")
  expect_identical(nn$settings$align_unrotated, "none")
  single <- EFA(cormat, n_factors = 3, N = 500, method = "PAF",
                rotation = "none")
  expect_equal(unclass(nn$unrot_loadings), unclass(single$unrot_loadings),
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("pooling identical imputations reproduces the single fit", {
  single <- EFA(cormat, n_factors = 3, N = 500, method = "PAF",
                rotation = "promax")

  # unrotated alignment is a pure sign/permutation step, so it must be exact
  expect_equal(unclass(pooled_obl$unrot_loadings),
               unclass(single$unrot_loadings),
               tolerance = 1e-10, ignore_attr = TRUE)
  # rotated solutions are re-derived by oblique Procrustes alignment, which
  # recovers the promax solution up to the solver tolerance
  expect_equal(unclass(pooled_obl$rot_loadings),
               unclass(single$rot_loadings),
               tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(pooled_obl$Phi, single$Phi, tolerance = 1e-4,
               ignore_attr = TRUE)
  expect_equal(pooled_obl$orig_R, cormat, ignore_attr = TRUE)
})

test_that("pooled components are internally consistent", {
  L <- unclass(pooled_obl$rot_loadings)
  Phi <- pooled_obl$Phi

  expect_equal(unclass(pooled_obl$Structure), L %*% Phi, ignore_attr = TRUE)
  expect_equal(unname(pooled_obl$h2), unname(diag(L %*% Phi %*% t(L))))

  res <- pooled_obl$orig_R - pooled_obl$model_implied_R
  diag(res) <- 0
  expect_equal(pooled_obl$residuals, res, ignore_attr = TRUE)
})

test_that("pooled fit indices D2-pool the imputation chi-squares", {
  fi <- pooled_none$fit_indices
  expect_identical(fi$pool_method, "D2")
  expect_equal(fi$df, pooled_none$fits[[1]]$fit_indices$df)
  expect_true(is.finite(fi$chi) && fi$chi >= 0)
  expect_true(is.finite(fi$p_chi))
  expect_true(is.finite(fi$RMSR))

  md <- pooled_none$mi_diagnostics
  expect_identical(md$m, 3L)
  expect_gte(md$ARIV, 0)
  expect_gte(md$FMI, 0)
  expect_lte(md$FMI, 1)

  # the pooled observed correlation matrix is the mean across imputations
  R_mean <- Reduce(`+`, lapply(pooled_none$fits, function(f) f$orig_R)) / 3
  expect_equal(pooled_none$orig_R, R_mean, ignore_attr = TRUE)
})

test_that("EFA_POOLED validates its arguments with classed conditions", {
  expect_error(
    EFA_POOLED(cormat_list, p = 0, n_factors = 3, N = 500, method = "PAF",
               rotation = "none"),
    class = "efa_pooled_bad_p"
  )
  expect_error(
    EFA_POOLED(cormat_list, rmsea_ci_level = 1, n_factors = 3, N = 500,
               method = "PAF", rotation = "none"),
    class = "efa_pooled_bad_ci_level"
  )
  expect_warning(
    EFA_POOLED(cormat_list, ci = .8, n_factors = 3, N = 500, method = "PAF",
               rotation = "none"),
    class = "efa_pooled_ci_ignored"
  )
})

test_that("EFA_POOLED rejects non-conformable imputations", {
  expect_error(
    EFA_POOLED(list(cormat, cormat[1:(p_vars - 1), 1:(p_vars - 1)]),
               n_factors = 3, N = 500, method = "PAF", rotation = "none"),
    class = "efa_pooled_dim_mismatch"
  )

  renamed <- cormat
  dimnames(renamed) <- list(paste0("X", seq_len(p_vars)),
                            paste0("X", seq_len(p_vars)))
  expect_error(
    EFA_POOLED(list(cormat, renamed), n_factors = 3, N = 500, method = "PAF",
               rotation = "none"),
    class = "efa_pooled_var_mismatch"
  )
})

test_that("EFA_POOLED warns when imputations have different N", {
  expect_warning(
    suppressMessages(
      EFA_POOLED(list(GRiPS_raw[1:300, ], GRiPS_raw[1:400, ]), n_factors = 1,
                 method = "PAF", rotation = "none")
    ),
    class = "efa_pooled_unequal_n"
  )
})

test_that("EFA_POOLED warns when N cannot be recovered", {
  # correlation-matrix input carries no N; a method that needs N for chi-square
  # fit cannot compute chi-square-based indices. The component ML fits also warn
  # about the missing N; muffle those so only the pooled class is asserted.
  suppressWarnings(
    expect_warning(
      suppressMessages(
        EFA_POOLED(cormat_list, n_factors = 3, method = "ML", rotation = "none")
      ),
      class = "efa_pooled_no_n"
    )
  )

  # N recoverable for the raw-data imputation but not the correlation matrix
  expect_warning(
    suppressMessages(
      EFA_POOLED(list(GRiPS_raw, stats::cor(GRiPS_raw)), n_factors = 1,
                 method = "PAF", rotation = "none")
    ),
    class = "efa_pooled_partial_n"
  )
})

test_that("bootstrap arrays are pooled into MI SEs and CIs", {
  skip_on_cran()
  local_reproducible_output()

  set.seed(1)
  boot_list <- lapply(1:2, function(i) {
    GRiPS_raw[sample(seq_len(nrow(GRiPS_raw)), 250, replace = TRUE), ]
  })
  pooled_boot <- suppressMessages(
    EFA_POOLED(boot_list, n_factors = 1, method = "ML",
               se = "np-boot", b_boot = 16)
  )

  expect_true(all(c("boot.SE", "boot.CI", "boot.arrays", "boot.MI",
                    "standardized_residuals") %in% names(pooled_boot)))

  L <- unclass(pooled_boot$unrot_loadings)
  se <- pooled_boot$boot.SE$unrot_loadings
  expect_identical(dim(se), dim(L))
  # the SEs must actually be computed (finite), not silently all-NA, and the
  # bootstrap must have produced real variation (positive SE somewhere)
  expect_true(all(is.finite(se)))
  expect_true(all(se >= 0))
  expect_true(any(se > 0))

  # Wald-type MI intervals: finite, correctly ordered, and centred on the
  # pooled point estimate. The centring check is the real test (it fails if the
  # interval is built around the wrong quantity or the bounds are swapped);
  # asserting only lower <= L <= upper would hold for any symmetric interval.
  ci <- pooled_boot$boot.CI$unrot_loadings
  expect_true(all(is.finite(ci$lower) & is.finite(ci$upper)))
  expect_true(all(ci$upper >= ci$lower))
  expect_equal((ci$lower + ci$upper) / 2, L, ignore_attr = TRUE)

  expect_equal(pooled_boot$settings$b_boot, 16)
  expect_identical(pooled_boot$settings$se, "np-boot")

  # FMIs must be computed (at least some finite) and in [0, 1]; an all-NA vector
  # would pass a bare all(..., na.rm = TRUE) range check vacuously.
  fmi <- pooled_boot$boot.MI$unrot_loadings$FMI
  expect_true(any(is.finite(fmi)))
  expect_true(all(fmi[is.finite(fmi)] >= 0 & fmi[is.finite(fmi)] <= 1))

  expect_identical(dim(pooled_boot$standardized_residuals),
                   dim(pooled_boot$residuals))

  # summary() additionally shows the MI uncertainty summary
  expect_snapshot(print(summary(pooled_boot)), transform = scrub_num)
})

test_that("oblique bootstrap pooling produces rotated SEs, CIs, and Phi", {
  skip_on_cran()

  set.seed(2)
  boot_list <- lapply(1:2, function(i) {
    GRiPS_raw[sample(seq_len(nrow(GRiPS_raw)), 250, replace = TRUE), ]
  })
  pooled_boot <- suppressWarnings(suppressMessages(
    EFA_POOLED(boot_list, n_factors = 2, method = "PAF", rotation = "promax",
               se = "np-boot", b_boot = 12)
  ))

  # the oblique branch pools rotated loadings, factor correlations, and
  # structure coefficients in addition to the unrotated quantities
  for (comp in c("rot_loadings", "Phi", "Structure")) {
    se <- pooled_boot$boot.SE[[comp]]
    ci <- pooled_boot$boot.CI[[comp]]
    expect_false(is.null(se))
    expect_true(all(is.finite(se)))
    expect_true(all(se >= 0))
    expect_true(all(is.finite(ci$lower) & is.finite(ci$upper)))
    expect_true(all(ci$upper >= ci$lower))
  }

  expect_identical(dim(pooled_boot$boot.SE$Phi), dim(as.matrix(pooled_boot$Phi)))
  expect_identical(dim(pooled_boot$boot.SE$Structure),
                   dim(unclass(pooled_boot$Structure)))

  fmi <- pooled_boot$boot.MI$Phi$FMI
  expect_true(any(is.finite(fmi)))
  expect_true(all(fmi[is.finite(fmi)] >= 0 & fmi[is.finite(fmi)] <= 1))
})

test_that("print.EFA_POOLED output is stable (PAF, promax)", {
  local_reproducible_output()

  expect_snapshot(print(pooled_obl), transform = scrub_num)
})

test_that("print.EFA_POOLED output is stable (ML, unrotated)", {
  local_reproducible_output()

  expect_snapshot(print(pooled_none), transform = scrub_num)
})

test_that("summary.EFA_POOLED output is stable (PAF, promax)", {
  local_reproducible_output()

  expect_snapshot(print(summary(pooled_obl)), transform = scrub_num)
})

test_that("summary.EFA_POOLED output is stable (ML, unrotated)", {
  local_reproducible_output()

  expect_snapshot(print(summary(pooled_none)), transform = scrub_num)
})

test_that("format.EFA_POOLED matches the printed output", {
  local_reproducible_output()

  expect_identical(format(pooled_obl),
                   utils::capture.output(print(pooled_obl)))
})
