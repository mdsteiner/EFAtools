# Tests for efa_mi(): structure of the returned object, pooling math,
# classed conditions, bootstrap/MI pooling, and the print/format methods.

cormat <- test_models$baseline$cormat
p_vars <- ncol(cormat)
cormat_list <- list(cormat, cormat, cormat)

pooled_obl <- efa_mi(cormat_list, n_factors = 3, N = 500, estimator = "PAF",
                         rotation = "promax")
pooled_orth <- efa_mi(cormat_list, n_factors = 3, N = 500, estimator = "PAF",
                          rotation = "varimax")

# The single-solution counterpart of `pooled_obl`: pooling identical imputations must
# reproduce it, so several blocks below compare against the same fit. PAF extraction and
# promax rotation are deterministic, so one fit serves all of them.
single_obl <- EFA(cormat, n_factors = 3, N = 500, method = "PAF", rotation = "promax")

set.seed(42)
grips_list <- lapply(1:3, function(i) {
  GRiPS_raw[sample(seq_len(nrow(GRiPS_raw)), replace = TRUE), ]
})
pooled_none <- suppressMessages(
  efa_mi(grips_list, n_factors = 1, estimator = "ML")
)

test_that("efa_mi returns a well-formed pooled object", {
  expect_s3_class(pooled_obl, c("efa_mi", "EFA_POOLED", "efa", "EFA"), exact = TRUE)

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
  expect_equal(pooled_obl$fit_indices$CAF, single_obl$fit_indices$CAF, tolerance = 1e-6)
  expect_gt(pooled_obl$fit_indices$CAF, 0)
  expect_lt(pooled_obl$fit_indices$CAF, 1)
})

test_that("efa_mi records the pooling settings", {
  s <- pooled_obl$settings
  expect_true(s$pooled)
  expect_identical(s$n_imputations, 3L)
  expect_identical(s$target_method, "first_target")
  expect_identical(s$align_unrotated, "signed_tucker_congruence")
  expect_identical(s$fit_pool_method, "D2")
  expect_equal(s$p, 0.05)
  expect_equal(s$ci, 0.95)
  # no bootstrap arrays were available, so the pooled object must not claim SEs;
  # neither the current names nor the historic flattened ones may appear.
  expect_identical(s$se, "none")
  expect_false(any(c("SE", "CI", "replicates", "MI",
                     "vcov_unrot_loadings", "Gamma",
                     "boot.SE", "boot.CI", "boot.arrays", "boot.MI") %in%
                     names(pooled_obl)))
})

test_that("efa_mi records the component fits' admissibility", {
  adm <- pooled_obl$mi_admissibility
  expect_identical(adm$m, 3L)
  expect_identical(adm$n_heywood_items, c(0L, 0L, 0L))
  expect_identical(adm$heywood_imputations, integer(0))
  expect_identical(adm$nonconverged, integer(0))
  expect_true(all(is.finite(adm$iter)))

  # The record is read off the component fits, so an inadmissible or non-converged
  # component is captured whatever the pooled matrix looks like.
  affected <- pooled_obl$fits
  affected[[1]]$heywood <- c(V1 = 1L, V5 = 5L)
  affected[[3]]$heywood <- c(V3 = 3L)
  affected[[3]]$convergence <- 1L
  adm_affected <- .efa_pooled_admissibility(affected)
  expect_identical(adm_affected$n_heywood_items, c(2L, 0L, 1L))
  expect_identical(adm_affected$heywood_imputations, c(1L, 3L))
  expect_identical(adm_affected$nonconverged, 3L)
})

test_that("summary() reports component Heywood cases behind a clean pooled count", {
  # Averaging aligned solutions pulls boundary communalities back inside the admissible
  # range, so a pool of improper solutions can have a proper pooled matrix. The pooled
  # count alone then reads as an all-clear for the analysis.
  local_reproducible_output()

  affected <- pooled_obl
  affected$fits[[1]]$heywood <- c(V1 = 1L, V5 = 5L)
  affected$fits[[3]]$heywood <- c(V3 = 3L)
  affected$mi_admissibility <- .efa_pooled_admissibility(affected$fits)

  heywood_line <- function(x) {
    grep("^Heywood cases:", cli::ansi_strip(format(summary(x))), value = TRUE)
  }
  # Three flags, not three distinct variables: a variable flagged in two imputations
  # counts once per imputation.
  expect_identical(heywood_line(affected),
                   "Heywood cases: 0 pooled (3 flags across 2 of 3 imputations)")

  # A pool whose component fits were all proper keeps the unqualified count.
  expect_identical(heywood_line(pooled_obl), "Heywood cases: 0")
})

test_that("a pool of improper component fits does not report zero Heywood cases", {
  skip_on_cran()
  local_reproducible_output()

  # One small imputation among larger ones: several component fits hit the uniqueness
  # boundary while their average stays well inside the admissible range.
  withr::local_seed(3)
  imps <- lapply(1:6, function(i) {
    n <- if (i == 1L) 45L else 300L
    UPPS_raw[sample.int(nrow(UPPS_raw), n, replace = TRUE), 1:12]
  })
  pooled <- suppressWarnings(suppressMessages(
    efa_mi(imps, n_factors = 4, estimator = "ML", rotation = "varimax")
  ))

  # precondition: the components really are improper and the pooled matrix is not
  expect_gt(sum(pooled$mi_admissibility$n_heywood_items), 0)
  expect_lt(max(pooled$h2), 1)

  expect_match(
    grep("^Heywood cases:", cli::ansi_strip(format(summary(pooled))), value = TRUE),
    "imputations)", fixed = TRUE
  )
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
  ft <- efa_mi(cormat_list, n_factors = 3, N = 500, estimator = "PAF",
                   rotation = "promax", target_method = "first_target")
  expect_s3_class(ft, c("efa_mi", "EFA_POOLED", "efa", "EFA"), exact = TRUE)
  expect_identical(ft$settings$target_method, "first_target")
  expect_identical(ft$alignment$method, "first_target")
  expect_s3_class(ft$rot_loadings, "LOADINGS")
  expect_identical(dim(ft$Phi), c(3L, 3L))

  # align_unrotated = "procrustes": orthogonal Procrustes of the unrotated axes
  pr <- efa_mi(cormat_list, n_factors = 3, N = 500, estimator = "PAF",
                   rotation = "none", align_unrotated = "procrustes")
  expect_identical(pr$settings$align_unrotated, "procrustes")
  expect_s3_class(pr$unrot_loadings, "LOADINGS")

  # align_unrotated = "none": average the unrotated loadings as returned. With
  # identical imputations the pooled result must equal the single fit.
  nn <- efa_mi(cormat_list, n_factors = 3, N = 500, estimator = "PAF",
                   rotation = "none", align_unrotated = "none")
  expect_identical(nn$settings$align_unrotated, "none")
  single <- EFA(cormat, n_factors = 3, N = 500, method = "PAF",
                rotation = "none")
  expect_equal(unclass(nn$unrot_loadings), unclass(single$unrot_loadings),
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("target_method = 'consensus' pools a solution that does not depend on the imputation order", {
  skip_on_cran()

  # One atypical imputation among larger ones. The GPA iteration keeps the gauge of the
  # solution it starts from, so this is exactly the setting in which the anchor decides
  # the pooled orientation.
  withr::local_seed(11)
  imps <- lapply(1:4, function(i) {
    n <- if (i == 1L) 60L else 250L
    GRiPS_raw[sample.int(nrow(GRiPS_raw), n, replace = TRUE), ]
  })
  fit_mi <- function(dl, tm, ...) {
    suppressWarnings(suppressMessages(
      efa_mi(dl, n_factors = 2, estimator = "ML", rotation = "varimax",
             target_method = tm, ...)
    ))
  }
  perm <- c(3L, 1L, 4L, 2L)

  consensus <- fit_mi(imps, "consensus")
  consensus_perm <- fit_mi(imps[perm], "consensus")
  expect_equal(unclass(consensus_perm$rot_loadings), unclass(consensus$rot_loadings),
               tolerance = 1e-8)
  expect_equal(as.numeric(consensus_perm$h2), as.numeric(consensus$h2),
               tolerance = 1e-8)

  # The GPA record is returned and the iteration is started at the medoid rotated
  # solution, which is a property of the set rather than of the list order.
  expect_true(isTRUE(consensus$alignment$converged))
  expect_length(consensus$alignment$aligned_loadings, length(imps))
  expect_equal(
    consensus$alignment$start,
    .efa_pooled_medoid_anchor(lapply(consensus$fits,
                                     function(f) unclass(f$rot_loadings)))
  )

  # The fixture is discriminating: the first-imputation anchor does move under the same
  # permutation, and by more than a relabelling of the columns.
  first_target <- fit_mi(imps, "first_target")
  first_target_perm <- fit_mi(imps[perm], "first_target")
  matched <- .align_solution(L_target = unclass(first_target$rot_loadings),
                             L = unclass(first_target_perm$rot_loadings))$loadings
  expect_gt(max(abs(unclass(matched) - unclass(first_target$rot_loadings))), 1e-3)

  # ... because "first_target" is anchored on the first imputation by definition.
  expect_equal(unclass(first_target$alignment$target),
               unclass(first_target$fits[[1]]$rot_loadings), tolerance = 1e-12)

  # An explicit start overrides the medoid anchor, and is honoured.
  from_first <- fit_mi(imps, "consensus", consensus_args = list(start = 1))
  expect_equal(from_first$alignment$start, 1)
})

test_that("procrustes_args forwards algorithm controls and nothing else", {
  target <- matrix(0, p_vars, 3)
  reserved <- list(list(A = target), list(Target = target),
                   list(rotation = "orthogonal"), list(S = diag(3)),
                   # partially spelled names are matched the way do.call() matches them
                   list(Tar = target))
  for (args in reserved) {
    expect_error(
      efa_mi(cormat_list, n_factors = 3, N = 500, estimator = "PAF",
             rotation = "promax", procrustes_args = args),
      class = "efa_pooled_bad_procrustes_args"
    )
  }

  # An unknown name is rejected here rather than surfacing as an unused-argument error
  # from inside the alignment, and an unnamed element cannot be matched at all.
  expect_error(
    efa_mi(cormat_list, n_factors = 3, N = 500, estimator = "PAF",
           rotation = "promax", procrustes_args = list(oblique_maxitt = 5)),
    class = "efa_pooled_bad_procrustes_args"
  )
  expect_error(
    efa_mi(cormat_list, n_factors = 3, N = 500, estimator = "PAF",
           rotation = "promax", procrustes_args = list(50)),
    class = "efa_pooled_bad_procrustes_args"
  )

  # A genuine control reaches efa_procrustes(): capping the oblique solver at one
  # iteration stops the alignment exactly there.
  tuned <- suppressWarnings(efa_mi(
    cormat_list, n_factors = 3, N = 500, estimator = "PAF", rotation = "promax",
    procrustes_args = list(oblique_maxit = 1L)
  ))
  expect_identical(
    vapply(tuned$alignment$target_rotations[-1L],
           function(x) as.numeric(x$iterations), numeric(1)),
    rep(1, length(cormat_list) - 1L)
  )
  expect_false(tuned$alignment$converged)
})

test_that("pooling identical imputations reproduces the single fit", {
  # unrotated alignment is a pure sign/permutation step, so it must be exact
  expect_equal(unclass(pooled_obl$unrot_loadings),
               unclass(single_obl$unrot_loadings),
               tolerance = 1e-10, ignore_attr = TRUE)
  # rotated solutions are re-derived by oblique Procrustes alignment, which
  # recovers the promax solution up to the solver tolerance
  expect_equal(unclass(pooled_obl$rot_loadings),
               unclass(single_obl$rot_loadings),
               tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(pooled_obl$Phi, single_obl$Phi, tolerance = 1e-4,
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

test_that("RMSR is the root mean square of the unique off-diagonal residuals", {
  for (pooled in list(pooled_obl, pooled_orth, pooled_none)) {
    res <- pooled$residuals
    expect_equal(pooled$fit_indices$RMSR,
                 sqrt(mean(res[upper.tri(res)]^2)))
    # The pooled residual matrix is symmetric, which is why counting each residual
    # pair twice instead of once cannot change that mean square.
    expect_equal(res, t(res), ignore_attr = TRUE)
    expect_equal(pooled$fit_indices$RMSR,
                 sqrt(mean(res[row(res) != col(res)]^2)))
  }
})

test_that("the deprecated rmsr_upper is accepted, warned about, and ignored", {
  expect_warning(
    deprecated <- efa_mi(cormat_list, n_factors = 3, N = 500, estimator = "PAF",
                         rotation = "promax", rmsr_upper = TRUE),
    class = "lifecycle_warning_deprecated"
  )
  expect_equal(deprecated$fit_indices$RMSR, pooled_obl$fit_indices$RMSR)

  expect_warning(
    flipped <- efa_mi(cormat_list, n_factors = 3, N = 500, estimator = "PAF",
                      rotation = "promax", rmsr_upper = FALSE),
    class = "lifecycle_warning_deprecated"
  )
  expect_equal(flipped$fit_indices$RMSR, pooled_obl$fit_indices$RMSR)

  # It no longer records a setting, and the frozen wrapper still accepts it without
  # passing it on, so legacy calls stay silent.
  expect_null(pooled_obl$settings$rmsr_upper)
  expect_no_warning(
    legacy <- EFA_POOLED(cormat_list, n_factors = 3, N = 500, method = "PAF",
                         rotation = "promax", rmsr_upper = TRUE),
    class = "lifecycle_warning_deprecated"
  )
  expect_equal(legacy$fit_indices$RMSR, pooled_obl$fit_indices$RMSR)
})

test_that("pooled fit indices D2-pool the imputation chi-squares", {
  fi <- pooled_none$fit_indices
  md <- pooled_none$mi_diagnostics
  expect_identical(fi$pool_method, "D2")
  # The reported names and their order are pinned by one shared definition, so a new
  # index cannot be added to the assembled list without also being placed there.
  expect_identical(names(fi), EFAtools:::.efa_pooled_fit_index_order)
  expect_equal(fi$df, pooled_none$fits[[1]]$fit_indices$df)
  expect_true(is.finite(fi$chi) && fi$chi >= 0)
  expect_true(is.finite(fi$p_chi))
  expect_true(is.finite(fi$RMSR))

  # The pooled set reports CFI/TLI as the average of the per-imputation indices
  # (kept in range and consistent with the component fits; verified in
  # test-EFA_POOLED-cfi-scale.R). The separately pooled (N - 1)-scale model and
  # baseline noncentralities remain exposed via mi_diagnostics chi_cfi /
  # chi_null_cfi for reconciliation against lavaan.mi, while ECVI keeps the
  # reported pooled chi-square. No mislabeled Fm.
  expect_true(is.finite(fi$TLI))
  expect_true(is.finite(fi$ECVI))
  expect_null(fi$Fm)
  expect_true(is.finite(md$chi_cfi))
  expect_true(is.finite(md$chi_null_cfi))
  N_used <- pooled_none$settings$N
  n_vars <- ncol(pooled_none$orig_R)
  n_params <- n_vars * (n_vars + 1) / 2 - fi$df
  expect_equal(fi$ECVI, (fi$chi + 2 * n_params) / (N_used - 1))

  # The diagnostics live in the top-level slot only; `fit_indices` holds scalars, the
  # same shape a single efa_fit() returns.
  expect_null(fi$mi_diagnostics)
  expect_true(all(vapply(fi, function(v) {
    (is.numeric(v) || is.character(v)) && length(v) == 1L
  }, logical(1))))

  expect_identical(md$m, 3L)
  expect_gte(md$ARIV, 0)
  expect_gte(md$FMI, 0)
  expect_lte(md$FMI, 1)

  # the pooled observed correlation matrix is the mean across imputations
  R_mean <- Reduce(`+`, lapply(pooled_none$fits, function(f) f$orig_R)) / 3
  expect_equal(pooled_none$orig_R, R_mean, ignore_attr = TRUE)
})

test_that(".efa_pooled_D2 ARIV matches the Li et al. (1991) formula", {
  # The average relative increase in variance is the between-imputation variance
  # of the sqrt-transformed statistics, i.e. (1 + 1/M) * var(sqrt(chi^2)).
  chis <- c(30, 35, 50)
  d2 <- .efa_pooled_D2(chis, df = 20)
  expect_equal(d2$ARIV, (1 + 1 / length(chis)) * stats::var(sqrt(chis)))
  expect_equal(d2$FMI, d2$ARIV / (1 + d2$ARIV))
})

test_that("a negative D2 statistic is reported as such, and the pooled fit stays in range", {
  # The pooling statistic is returned unfloored, so it can go negative when the
  # between-imputation variability of the component statistics exceeds the pooled
  # discrepancy; that is a diagnostic of the pool. The fit it feeds must not inherit
  # it: the pooled chi-square is floored at zero and its p-value is 1.
  d2 <- .efa_pooled_D2(c(10, 20, 40, 60, 80), df = 33)
  expect_lt(d2$F, 0)
  expect_identical(d2$chi, 0)
  expect_identical(d2$p, 1)
  expect_gt(d2$ARIV, 0)

  # An ordinary pool keeps a positive statistic and a chi-square derived from it.
  ok <- .efa_pooled_D2(c(30, 35, 50), df = 20)
  expect_gt(ok$F, 0)
  expect_equal(ok$chi, ok$df1 * ok$F)
})

test_that("seed and b_boot in the dots govern the whole pooled call", {
  skip_on_cran()

  withr::local_seed(20260812)
  boot_list <- lapply(1:2, function(i) {
    GRiPS_raw[sample(seq_len(nrow(GRiPS_raw)), 200, replace = TRUE), 1:5]
  })
  run <- function(seed) {
    suppressWarnings(suppressMessages(
      efa_mi(boot_list, n_factors = 1, estimator = "ML", se = "np-boot",
             b_boot = 4, seed = seed)
    ))
  }

  a <- run(99)
  b <- run(99)
  expect_identical(a$SE$unrot_loadings, b$SE$unrot_loadings)
  # b_boot is the per-imputation replicate count and is recorded as such
  expect_equal(a$settings$b_boot, 4)
  expect_identical(dim(a$replicates$unrot_loadings[[1]])[1], 4L)

  # A different seed moves the bootstrap SEs, so the reproducibility above is not
  # vacuous ...
  expect_false(isTRUE(all.equal(a$SE$unrot_loadings, run(100)$SE$unrot_loadings)))

  # ... and the caller's own random stream is left where it was.
  set.seed(7)
  expected <- stats::runif(1)
  set.seed(7)
  invisible(run(99))
  expect_identical(stats::runif(1), expected)
})

test_that("efa_mi validates its arguments with classed conditions", {
  expect_error(
    efa_mi(cormat_list, p = 0, n_factors = 3, N = 500, estimator = "PAF",
               rotation = "none"),
    class = "efa_pooled_bad_p"
  )
  expect_error(
    efa_mi(cormat_list, rmsea_ci_level = 1, n_factors = 3, N = 500,
               estimator = "PAF", rotation = "none"),
    class = "efa_pooled_bad_ci_level"
  )
  expect_warning(
    efa_mi(cormat_list, ci = .8, n_factors = 3, N = 500, estimator = "PAF",
               rotation = "none"),
    class = "efa_pooled_ci_ignored"
  )
  # Too few imputations signals the documented condition class, not a bare
  # assertion error, and does so before any component fit is run.
  expect_error(
    efa_mi(list(cormat), n_factors = 3, N = 500, estimator = "PAF",
               rotation = "none"),
    class = "efa_pooled_min_fits"
  )
  expect_error(
    efa_mi(list(), n_factors = 3, N = 500, estimator = "PAF", rotation = "none"),
    class = "efa_pooled_min_fits"
  )
  # A mids object is itself a list, so without the dedicated check it would fail deep
  # inside the per-dataset assertion instead of naming the conversion.
  expect_error(
    efa_mi(structure(list(data = cormat, m = 2L), class = "mids"),
               n_factors = 3, N = 500, estimator = "PAF", rotation = "none"),
    class = "efa_pooled_mids_input"
  )
})

test_that("a failing component fit names its imputation and keeps the original condition", {
  local_reproducible_output()

  mk <- function(seed) {
    set.seed(seed)
    f <- stats::rnorm(60)
    X <- outer(f, rep(0.7, 4)) + matrix(stats::rnorm(240, sd = 0.7), 60, 4)
    colnames(X) <- paste0("V", seq_len(4))
    X
  }
  raw <- lapply(1:3, mk)
  # A constant column in one imputation only: the pool is all-or-nothing, so the whole
  # call ends, and without the index the user has to bisect `data_list` by hand.
  raw[[2]][, 3] <- 1

  err <- expect_error(
    efa_mi(raw, n_factors = 1, estimator = "PAF", rotation = "none"),
    class = "efa_pooled_fit_failed"
  )
  # the component condition is chained, not replaced, so its own diagnosis survives
  expect_s3_class(err$parent, "efa_cor_uncomputable")

  expect_snapshot(
    error = TRUE,
    efa_mi(raw, n_factors = 1, estimator = "PAF", rotation = "none")
  )
})

test_that("efa_mi rejects non-conformable imputations", {
  expect_error(
    efa_mi(list(cormat, cormat[1:(p_vars - 1), 1:(p_vars - 1)]),
               n_factors = 3, N = 500, estimator = "PAF", rotation = "none"),
    class = "efa_pooled_dim_mismatch"
  )

  renamed <- cormat
  dimnames(renamed) <- list(paste0("X", seq_len(p_vars)),
                            paste0("X", seq_len(p_vars)))
  expect_error(
    efa_mi(list(cormat, renamed), n_factors = 3, N = 500, estimator = "PAF",
               rotation = "none"),
    class = "efa_pooled_var_mismatch"
  )
})

test_that("efa_mi warns when imputations have different N", {
  expect_warning(
    suppressMessages(
      efa_mi(list(GRiPS_raw[1:300, ], GRiPS_raw[1:400, ]), n_factors = 1,
                 estimator = "PAF", rotation = "none")
    ),
    class = "efa_pooled_unequal_n"
  )
})

test_that("efa_mi warns when N cannot be recovered", {
  # correlation-matrix input carries no N; a method that needs N for chi-square
  # fit cannot compute chi-square-based indices. The component ML fits also warn
  # about the missing N; muffle those so only the pooled class is asserted.
  suppressWarnings(
    expect_warning(
      suppressMessages(
        efa_mi(cormat_list, n_factors = 3, estimator = "ML", rotation = "none")
      ),
      class = "efa_pooled_no_n"
    )
  )

  # N recoverable for the raw-data imputation but not the correlation matrix
  expect_warning(
    suppressMessages(
      efa_mi(list(GRiPS_raw, stats::cor(GRiPS_raw)), n_factors = 1,
                 estimator = "PAF", rotation = "none")
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
    efa_mi(boot_list, n_factors = 1, estimator = "ML",
               se = "np-boot", b_boot = 6)
  )

  expect_true(all(c("SE", "CI", "replicates", "MI",
                    "standardized_residuals") %in% names(pooled_boot)))
  # the historic flattened slot names must not leak in alongside the new ones
  expect_false(any(c("boot.SE", "boot.CI", "boot.arrays", "boot.MI") %in%
                     names(pooled_boot)))

  L <- unclass(pooled_boot$unrot_loadings)
  se <- pooled_boot$SE$unrot_loadings
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
  ci <- pooled_boot$CI$unrot_loadings
  expect_true(all(is.finite(ci$lower) & is.finite(ci$upper)))
  expect_true(all(ci$upper >= ci$lower))
  expect_equal((ci$lower + ci$upper) / 2, L, ignore_attr = TRUE)

  expect_equal(pooled_boot$settings$b_boot, 6)
  expect_identical(pooled_boot$settings$se, "np-boot")

  # Communality SEs/CIs use efa_fit()'s name on this route too, with `h2` kept as an
  # alias of it; the MI diagnostics carry the canonical name alone, so each family
  # enters the printed FMI/RIV summary once.
  expect_false(is.null(pooled_boot$SE$communalities))
  expect_identical(pooled_boot$SE$h2, pooled_boot$SE$communalities)
  expect_identical(pooled_boot$CI$h2, pooled_boot$CI$communalities)
  expect_false(is.null(pooled_boot$MI$communalities))
  expect_null(pooled_boot$MI$h2)

  # FMIs must be computed (at least some finite) and in [0, 1]; an all-NA vector
  # would pass a bare all(..., na.rm = TRUE) range check vacuously.
  fmi <- pooled_boot$MI$unrot_loadings$FMI
  expect_true(any(is.finite(fmi)))
  expect_true(all(fmi[is.finite(fmi)] >= 0 & fmi[is.finite(fmi)] <= 1))

  expect_identical(dim(pooled_boot$standardized_residuals),
                   dim(pooled_boot$residuals))

  fit_ci <- pooled_boot$CI$fit_indices_descriptive
  expect_true(all(is.finite(fit_ci$lower[c("SRMR", "TLI", "ECVI")])))
  expect_true(all(is.finite(fit_ci$upper[c("SRMR", "TLI", "ECVI")])))
  expect_true("RMSR" %in% names(pooled_boot$fit_indices))

  summary_lines <- cli::ansi_strip(format(summary(pooled_boot)))
  expect_false(any(grepl("^RMSR\\b", summary_lines)))
  expect_true(any(grepl("^SRMR \\[95% bootstrap/MI-CI\\]:", summary_lines)))
  # The incremental indices carry the averaged-over-imputations label ahead of the CI tag
  expect_true(any(grepl("^TLI \\(avg\\. over imputations\\) \\[95% bootstrap/MI-CI\\]:",
                        summary_lines)))
  expect_true(any(grepl("^ECVI \\[95% bootstrap/MI-CI\\]:", summary_lines)))

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
    efa_mi(boot_list, n_factors = 2, estimator = "PAF", rotation = "promax",
               se = "np-boot", b_boot = 6)
  ))

  # the oblique branch pools rotated loadings, factor correlations, and
  # structure coefficients in addition to the unrotated quantities
  for (comp in c("rot_loadings", "Phi", "Structure")) {
    se <- pooled_boot$SE[[comp]]
    ci <- pooled_boot$CI[[comp]]
    expect_false(is.null(se))
    expect_true(all(is.finite(se)))
    expect_true(all(se >= 0))
    expect_true(all(is.finite(ci$lower) & is.finite(ci$upper)))
    expect_true(all(ci$upper >= ci$lower))
  }

  expect_identical(dim(pooled_boot$SE$Phi), dim(as.matrix(pooled_boot$Phi)))
  expect_identical(dim(pooled_boot$SE$Structure),
                   dim(unclass(pooled_boot$Structure)))

  fmi <- pooled_boot$MI$Phi$FMI
  expect_true(any(is.finite(fmi)))
  expect_true(all(fmi[is.finite(fmi)] >= 0 & fmi[is.finite(fmi)] <= 1))
})

test_that("a failed bootstrap replicate is skipped, not fatal", {
  # A component EFA NA-fills a bootstrap replicate it could not fit; the pooled
  # bootstrap must skip that replicate (classed warning) and still produce finite
  # SEs from the valid ones, rather than aborting on the NA in the alignment step.
  withr::local_seed(1)
  p <- 3L; k <- 1L; B <- 4L; m <- 2L
  mk_arr <- function(na = NULL) {
    a <- array(stats::rnorm(p * k * B, 0.6, 0.05), dim = c(p, k, B))
    if (!is.null(na)) a[, , na] <- NA_real_
    a
  }
  mk_res <- function(na = NULL) {
    a <- array(stats::rnorm(p * p * B, 0, 0.02), dim = c(p, p, B))
    if (!is.null(na)) a[, , na] <- NA_real_
    a
  }
  # bootstrap fit-index arrays so the Rubin-Wald descriptive fit path also runs
  mk_fit <- function(na = NULL) {
    f <- matrix(stats::rnorm(B * 2, 1, 0.1), nrow = B,
                dimnames = list(NULL, c("chi", "CFI")))
    if (!is.null(na)) f[na, ] <- NA_real_
    f
  }
  L <- matrix(0.6, p, k)
  fits <- list(
    list(fit_indices = list(chi = 1, CFI = 1),
         replicates = list(unrot_loadings = mk_arr(na = 2), residuals = mk_res(na = 2),
                           fit_indices = mk_fit(na = 2))),
    list(fit_indices = list(chi = 1, CFI = 1),
         replicates = list(unrot_loadings = mk_arr(),       residuals = mk_res(),
                           fit_indices = mk_fit()))
  )
  orig_R <- replicate(m, { R <- matrix(0.4, p, p); diag(R) <- 1; R }, simplify = FALSE)
  args <- list(
    fits = fits, orig_R_list = orig_R,
    unrot_loadings_aligned = replicate(m, L, simplify = FALSE),
    mean_unrot_loadings = L, rotation_type = "none",
    align_unrotated = "signed_tucker_congruence",
    h2 = rep(0.36, p), residuals = matrix(0, p, p)
  )

  expect_warning(
    pooled <- do.call(.efa_pooled_bootstrap_pool, args),
    class = "efa_pooled_boot_failed"
  )
  expect_true(all(is.finite(pooled$SE$unrot_loadings)))
  # the Rubin-Wald descriptive fit path also ran (its failed replicate skipped,
  # not fatal) and produced finite fit-index SEs
  expect_false(is.null(pooled$SE$fit_indices_descriptive))
  expect_true(all(is.finite(pooled$SE$fit_indices_descriptive)))

  # The skipped (NA-filled) replicate is tallied as a source failure, and the
  # valid-rotation count subtracts it (not only rotation failures), so it never
  # overstates the replicates that entered the pool.
  expect_identical(pooled$MI$bootstrap_source_failures, c(1L, 0L))
  expect_identical(pooled$MI$bootstrap_rotation_failures, c(0L, 0L))
  expect_identical(pooled$MI$bootstrap_rotation_valid, c(B - 1L, B))

  # if an imputation is left with fewer than two valid replicates, no SEs can be
  # computed and the pooled bootstrap returns NULL (the existing "no SEs" path)
  fits_fail <- fits
  fits_fail[[1]]$replicates$unrot_loadings <- mk_arr(na = seq_len(B))
  args_fail <- args
  args_fail$fits <- fits_fail
  expect_warning(
    pooled_fail <- do.call(.efa_pooled_bootstrap_pool, args_fail),
    class = "efa_pooled_boot_insufficient"
  )
  expect_null(pooled_fail)
})

test_that("the pooled bootstrap sample count reports the usable replicates when they differ", {
  # A pooled fit records per-imputation FAILURE counts rather than a survivor count, so without a
  # fallback the printed sample count states the requested b_boot whatever the survival rate -- the
  # same claim of an unearned precision that a single fit used to make.
  spec <- list(b_boot = 4L, is_pooled = TRUE, np_boot = TRUE)
  x <- list(MI = list(bootstrap_source_failures = c(1L, 0L)),
            settings = list(b_boot = 4L))

  expect_identical(EFAtools:::.efa_valid_replicate_counts(x), c(3L, 4L))
  expect_identical(EFAtools:::.efa_usable_replicate_text(x, spec),
                   " (between 3 and 4 usable per imputation)")

  # A run in which every replicate survived, in either shape, says nothing extra.
  x_clean <- list(MI = list(bootstrap_source_failures = c(0L, 0L)),
                  settings = list(b_boot = 4L))
  expect_identical(EFAtools:::.efa_usable_replicate_text(x_clean, spec), "")
  expect_identical(
    EFAtools:::.efa_usable_replicate_text(list(SE = list(valid_replicates = 4L)), spec), "")

  # A single fit keeps the unqualified wording.
  expect_identical(
    EFAtools:::.efa_usable_replicate_text(list(SE = list(valid_replicates = 3L)),
                                          list(b_boot = 4L)),
    " (3 usable)")
})

test_that("print.efa_mi output is stable (PAF, promax)", {
  local_reproducible_output()

  expect_snapshot(print(pooled_obl), transform = scrub_num)
})

test_that("print.efa_mi output is stable (ML, unrotated)", {
  local_reproducible_output()

  expect_snapshot(print(pooled_none), transform = scrub_num)
})

test_that("summary.efa_mi output is stable (PAF, promax)", {
  local_reproducible_output()

  expect_snapshot(print(summary(pooled_obl)), transform = scrub_num)
})

test_that("summary.efa_mi output is stable (ML, unrotated)", {
  local_reproducible_output()

  expect_snapshot(print(summary(pooled_none)), transform = scrub_num)
})

test_that("the pooled header and pooling settings follow the console width", {
  # The pooled header names the analysis, and its settings line is the longest line of the
  # report; both are emitted verbatim so a "setting = 'value'" token is never split, and both
  # therefore have to be packed to the console at the separators between those tokens.

  # The settings an item reports, independent of where its lines happen to break.
  tokens <- function(lines) {
    txt <- sub("^Pooling( settings)?: ", "", paste(trimws(lines), collapse = " "))
    strsplit(sub("\\.$", "", txt), ", ", fixed = TRUE)[[1L]]
  }
  at_width <- function(x, w) withr::with_options(list(cli.width = w), cli::ansi_strip(format(x)))

  reference <- tokens(wrapped_item(at_width(pooled_obl, 120L), "^Pooling settings:"))
  expect_length(reference, 3L)

  for (w in c(120L, 80L, 60L)) {
    out <- at_width(pooled_obl, w)
    lines <- c(wrapped_item(out, "^Pooled EFA"), wrapped_item(out, "^Pooling settings:"))
    expect_true(all(cli::ansi_nchar(lines, type = "width") <= w))
    # wrapping only moves the line breaks: the same settings, in the same order ...
    expect_identical(tokens(wrapped_item(out, "^Pooling settings:")), reference)
    # ... and no line ends inside a "setting = 'value'" token
    expect_false(any(grepl("= '[^']*$", lines)))
    expect_true(any(grepl("estimator = 'PAF'", lines, fixed = TRUE)))
    expect_true(any(grepl("rotation = 'promax'", lines, fixed = TRUE)))

    # the summary's diagnostics entry reports the same settings and is packed the same way
    pooling <- wrapped_item(at_width(summary(pooled_obl), w), "^Pooling:")
    expect_true(all(cli::ansi_nchar(pooling, type = "width") <= w))
    expect_identical(tokens(pooling), reference)
  }
})

test_that("the wrapped pooled header keeps its styling when colours are on", {
  withr::local_options(cli.num_colors = 256, cli.width = 60)

  out <- format(pooled_obl)
  # the estimator/rotation values stay emphasised after the line is packed ...
  expect_true(cli::ansi_has_any(out[2L]) || cli::ansi_has_any(out[3L]))
  # ... and with colours off the same lines are plain text, wrapped identically
  withr::local_options(cli.num_colors = 1)
  plain <- format(pooled_obl)
  expect_false(cli::ansi_has_any(paste(plain[1:4], collapse = "")))
  expect_identical(cli::ansi_strip(out[1:4]), plain[1:4])
})

test_that("format.efa_mi matches the printed output", {
  local_reproducible_output()

  expect_identical(format(pooled_obl),
                   utils::capture.output(print(pooled_obl)))
})

test_that("the averaged-index flag names only the indices that were reported", {
  local_reproducible_output()

  inc_lines <- function(x) {
    lines <- cli::ansi_strip(format(x))
    lines[grepl("^CFI|^TLI|averaged over the imputations", lines)]
  }

  # Both incremental indices defined: both are labelled and the note is plural.
  both <- inc_lines(pooled_none)
  expect_true(any(grepl("^CFI \\(avg\\. over imputations\\):", both)))
  expect_true(any(grepl("^TLI \\(avg\\. over imputations\\):", both)))
  expect_true(any(grepl("CFI and TLI are averaged over the imputations", both)))

  # A degenerate baseline leaves TLI undefined: its line is dropped, so the note must
  # name CFI alone and in the singular rather than claiming a TLI that was not printed.
  one <- pooled_none
  one$fit_indices$TLI <- NA_real_
  one_lines <- inc_lines(one)
  expect_true(any(grepl("^CFI \\(avg\\. over imputations\\):", one_lines)))
  expect_false(any(grepl("^TLI", one_lines)))
  expect_true(any(grepl("CFI is averaged over the imputations", one_lines)))

  # Neither index defined: no label on the bare CFI line and no note at all.
  none <- pooled_none
  none$fit_indices$CFI <- NA_real_
  none$fit_indices$TLI <- NA_real_
  none_lines <- cli::ansi_strip(format(none))
  expect_true(any(grepl("^CFI: ", none_lines)))
  expect_false(any(grepl("avg\\. over imputations", none_lines)))
  expect_false(any(grepl("averaged over the imputations", none_lines)))

  # The D2 chi-square note is independent of the incremental-index note.
  expect_true(any(grepl("the pooled .* is the D2 statistic", none_lines)))
})

# ---- Default unrotated alignment: medoid anchor and canonical gauge ---------

# These use resampled imputations on purpose. A `list(cormat, cormat, cormat)`
# fixture has no between-imputation variance, so the aligned solutions are
# identical, their average is still in the extraction's gauge, and neither the
# anchor choice nor the re-gauging can be detected.

make_resamples <- function(dat, n = 5, seed = 42) {
  set.seed(seed)
  lapply(seq_len(n), function(i) {
    dat[sample.int(nrow(dat), replace = TRUE), , drop = FALSE]
  })
}

test_that("pooled unrotated loadings do not depend on the order of data_list", {
  # A weakly determined second factor makes the column matching ambiguous, which
  # is precisely when the anchor choice matters: with a first-imputation anchor
  # the pooled loadings move by up to 0.84 across orderings of this fixture.
  # Where the factors are well separated the matching is unambiguous and every
  # anchor already agrees, so such a fixture could not detect the difference.
  imps <- make_resamples(GRiPS_raw)
  base <- suppressWarnings(suppressMessages(
    efa_mi(imps, n_factors = 2, estimator = "ML", rotation = "none")))

  set.seed(7)
  for (i in seq_len(3)) {
    permuted <- suppressWarnings(suppressMessages(
      efa_mi(imps[sample(length(imps))], n_factors = 2, estimator = "ML",
             rotation = "none")))
    expect_equal(unclass(permuted$unrot_loadings),
                 unclass(base$unrot_loadings), tolerance = 1e-8)
    expect_equal(as.numeric(permuted$h2), as.numeric(base$h2), tolerance = 1e-8)
    expect_equal(permuted$fit_indices$RMSR, base$fit_indices$RMSR,
                 tolerance = 1e-8)
  }
})

test_that("pooled unrotated loadings are returned in the extraction's gauge", {
  imps <- make_resamples(DOSPERT_raw)

  # ML identifies the unrotated solution by diagonal L' Psi^-1 L ...
  ml <- suppressWarnings(suppressMessages(
    efa_mi(imps, n_factors = 3, estimator = "ML", rotation = "none")))
  L <- unclass(ml$unrot_loadings)
  A <- crossprod(L, L / pmax(1 - rowSums(L^2), 1e-6))
  expect_lt(sum(abs(A[upper.tri(A)])) / sum(abs(diag(A))), 1e-8)

  # ... a principal-axis extraction by diagonal L'L, and the pooled solution
  # must follow whichever its component fits use.
  paf <- suppressWarnings(suppressMessages(
    efa_mi(imps, n_factors = 3, estimator = "PAF", rotation = "none")))
  Lp <- unclass(paf$unrot_loadings)
  G <- crossprod(Lp)
  expect_lt(sum(abs(G[upper.tri(G)])) / sum(abs(diag(G))), 1e-8)
})

test_that("the canonical re-gauging leaves gauge-invariant quantities alone", {
  imps <- make_resamples(DOSPERT_raw)
  fits <- lapply(imps, function(d) suppressWarnings(suppressMessages(
    efa_fit(d, n_factors = 3, estimator = "ML", rotation = NULL))))
  Ls <- lapply(fits, function(f) unclass(f$unrot_loadings))

  aligned <- .efa_pooled_align_unrotated_list(
    Ls, align_unrotated = "signed_tucker_congruence")
  C <- aligned$meta[[1]]$C
  expect_false(is.null(C))
  expect_equal(crossprod(C), diag(ncol(C)), tolerance = 1e-10)

  # Undo the common rotation to recover the pre-canonical average, then check
  # that the quantities a rotation cannot change are indeed unchanged.
  M <- .average_matrices(aligned$loadings)
  raw <- M %*% t(C)
  expect_equal(rowSums(M^2), rowSums(raw^2), tolerance = 1e-10)
  expect_equal(M %*% t(M), raw %*% t(raw), tolerance = 1e-10)
})

test_that("analytic SEs survive the column-mixing gauge rotation", {
  imps <- make_resamples(DOSPERT_raw)
  pooled <- suppressWarnings(suppressMessages(
    efa_mi(imps, n_factors = 3, estimator = "ML", rotation = "none",
           se = "information")))

  se <- pooled$SE$unrot_loadings
  expect_false(is.null(se))
  expect_true(all(is.finite(as.matrix(se))))
  expect_true(all(as.matrix(se) > 0))

  # Wald intervals stay centred on the pooled estimate they belong to.
  est <- unclass(pooled$unrot_loadings)
  lo <- pooled$CI$unrot_loadings$lower
  hi <- pooled$CI$unrot_loadings$upper
  expect_equal(unname((lo + hi) / 2), unname(est), tolerance = 1e-8)
})

test_that("an unreliable per-imputation SE still blanks the pooled element", {
  # The re-gauging reads the full covariance block rather than the marginal SEs,
  # so the marginal NA mask has to be applied explicitly. Because the rotation
  # mixes columns, one NA contaminates its whole row.
  imps <- make_resamples(DOSPERT_raw)
  pooled <- suppressWarnings(suppressMessages(
    efa_mi(imps, n_factors = 3, estimator = "ML", rotation = "none",
           se = "information")))

  fits <- pooled$fits
  fits[[3]]$SE$unrot_loadings[2, 1] <- NA_real_
  Ls <- lapply(fits, function(f) unclass(f$unrot_loadings))
  aligned <- .efa_pooled_align_unrotated_list(
    Ls, align_unrotated = "signed_tucker_congruence")
  na_pool <- .efa_pooled_analytic_pool(
    fits = fits, unrot_loadings_aligned = aligned$loadings,
    align_meta = aligned$meta, ci = 0.95)

  expect_true(all(is.na(na_pool$SE$unrot_loadings[2, ])))
  expect_true(all(is.finite(na_pool$SE$unrot_loadings[3, ])))
})

test_that("a missing gauge covariance is reported once, not by both unrotated warnings", {
  # Dropping the covariance block leaves that imputation's row all-NA, which is also
  # the shape the NA-filled-SE warning keys on. The two must not both fire for one
  # downgrade: the missing-covariance branch is the specific diagnosis.
  imps <- make_resamples(DOSPERT_raw)
  pooled <- suppressWarnings(suppressMessages(
    efa_mi(imps, n_factors = 3, estimator = "ML", rotation = "none",
           se = "information")))

  fits <- pooled$fits
  fits[[2]]$vcov_unrot_loadings <- NULL
  Ls <- lapply(fits, function(f) unclass(f$unrot_loadings))
  aligned <- .efa_pooled_align_unrotated_list(
    Ls, align_unrotated = "signed_tucker_congruence")
  # precondition: the gauge rotation is active, so the covariance block is required
  expect_false(is.null(aligned$meta[[1]]$C))

  expect_warning(
    gauge_pool <- withCallingHandlers(
      .efa_pooled_analytic_pool(fits = fits,
                                unrot_loadings_aligned = aligned$loadings,
                                align_meta = aligned$meta, ci = 0.95),
      efa_pooled_unrotated_se_unreliable = function(w) {
        testthat::fail("the missing-covariance downgrade was reported twice")
      }
    ),
    class = "efa_pooled_gauge_vcov_missing"
  )
  expect_true(all(is.na(gauge_pool$SE$unrot_loadings)))
})

test_that("the medoid anchor is order-invariant with only two imputations", {
  # Two imputations are tied by construction -- each is exactly as far from the
  # other -- so the anchor has to be settled on the candidates' own content
  # rather than on list position, or the tie puts the order dependence back.
  set.seed(42)
  imps <- lapply(seq_len(2), function(i) {
    GRiPS_raw[sample.int(nrow(GRiPS_raw), replace = TRUE), , drop = FALSE]
  })
  a <- suppressWarnings(suppressMessages(
    efa_mi(imps, n_factors = 2, estimator = "ML", rotation = "none")))
  b <- suppressWarnings(suppressMessages(
    efa_mi(rev(imps), n_factors = 2, estimator = "ML", rotation = "none")))

  expect_equal(unclass(a$unrot_loadings), unclass(b$unrot_loadings),
               tolerance = 1e-10)
  expect_equal(as.numeric(a$h2), as.numeric(b$h2), tolerance = 1e-10)

  Ls <- lapply(a$fits, function(x) unclass(x$unrot_loadings))
  expect_identical(.efa_pooled_medoid_anchor(Ls), 1L)
  expect_identical(.efa_pooled_medoid_anchor(rev(Ls)), 2L)
})

test_that("re-gauging stands down when the gauge is not identified", {
  Ls <- lapply(1:3, function(i) unclass(
    efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
            estimator = "ML", rotation = NULL)$unrot_loadings))

  # An improper average would enter Psi^-1 at the clamp and let one variable
  # decide the gauge, so it is left alone.
  heywood <- Ls[[1]]
  heywood[1, ] <- heywood[1, ] / sqrt(sum(heywood[1, ]^2)) * 1.05
  expect_identical(.efa_pooled_canonical_gauge(heywood, Ls[[1]]),
                   diag(ncol(heywood)))

  # So is one whose canonical values are tied, where any rotation inside the
  # tied subspace diagonalises equally well.
  tied <- matrix(0, 6, 2)
  tied[1:3, 1] <- 0.6
  tied[4:6, 2] <- 0.6
  expect_identical(.efa_pooled_canonical_gauge(tied, tied), diag(2))
})

test_that("re-gauging still fires when a uniqueness sits on the estimation floor", {
  # Psi is reconstructed from the pooled loadings, and maximum likelihood pins a
  # uniqueness at .uniqueness_floor. There 1 - rowSums(L^2) lands a whisker below
  # the value that actually diagonalised the fit, and Psi^-1 multiplies that gap
  # by ~1/floor. Reconstructing with a smaller floor pushes the measured defect
  # past the detection tolerance, which would silently switch the re-gauging off
  # for exactly the boundary solutions.
  set.seed(5)
  imps <- lapply(seq_len(5), function(i) {
    UPPS_raw[sample.int(nrow(UPPS_raw), 60, replace = TRUE), 1:12, drop = FALSE]
  })
  fits <- lapply(imps, function(d) suppressWarnings(suppressMessages(
    efa_fit(d, n_factors = 4, estimator = "ML", rotation = NULL))))

  # precondition: at least one fit really is on the floor, or the test is vacuous
  min_psi <- min(vapply(fits, function(f) min(1 - as.numeric(f$h2)), numeric(1)))
  expect_lte(min_psi, .uniqueness_floor + 1e-9)

  Ls <- lapply(fits, function(f) unclass(f$unrot_loadings))
  anchor <- Ls[[.efa_pooled_medoid_anchor(Ls)]]
  defect <- function(A) sum(abs(A[upper.tri(A)])) / sum(abs(diag(A)))
  wgram <- function(L, floor) crossprod(L, L / pmax(1 - rowSums(L^2), floor))

  # too small a floor makes the anchor look like it is in no recognised gauge ...
  expect_gt(defect(wgram(anchor, 1e-6)), 1e-4)
  # ... while the estimation floor recovers the constraint the fit satisfies
  expect_lt(defect(wgram(anchor, .uniqueness_floor)), 1e-4)

  # so detection must fire and the pooled matrix must come back canonical
  aligned <- .efa_pooled_align_unrotated_list(
    Ls, align_unrotated = "signed_tucker_congruence")
  expect_false(is.null(aligned$meta[[1]]$C))

  pooled <- suppressWarnings(suppressMessages(
    efa_mi(imps, n_factors = 4, estimator = "ML", rotation = "none")))
  L <- unclass(pooled$unrot_loadings)
  expect_lt(defect(wgram(L, .uniqueness_floor)), 1e-8)
})

test_that("pooled AIC/BIC/ECVI are withheld when the component fits are scaled", {
  # A missing-at-random fixture: column 1 is fully observed and drives the
  # missingness in the others, so cor_method = "fiml" fits a genuine two-stage
  # model and each component carries the corrected (scaled-shifted) statistic.
  imps <- lapply(1:3, function(i) {
    set.seed(100 + i)
    L <- matrix(0, 6, 2)
    L[1:3, 1] <- 0.7
    L[4:6, 2] <- 0.7
    S <- tcrossprod(L)
    diag(S) <- 1
    # Colour the deviates with the Cholesky factor, not an eigendecomposition: this population
    # covariance has repeated eigenvalues, whose eigenvector basis is undetermined and settled
    # by rounding, whereas chol() is unique for a positive definite matrix, so the seed above
    # reproduces the same imputation on every LAPACK build.
    X <- matrix(stats::rnorm(400 * 6), 400) %*% chol(S)
    colnames(X) <- paste0("V", seq_len(6))
    X[X[, 1] > 0.8, 2] <- NA
    X[X[, 1] < -0.8, 3] <- NA
    X[X[, 1] > 1.2, 4] <- NA
    X
  })

  pooled <- suppressMessages(suppressWarnings(
    efa_mi(imps, n_factors = 2, estimator = "ML", rotation = "none",
           cor_method = "fiml")
  ))

  # precondition: the components really are scaled and withhold the three indices
  expect_true(all(vapply(pooled$fits,
                         function(f) !is.null(f$fit_indices$chi_scaled_type),
                         logical(1))))
  expect_true(all(vapply(pooled$fits,
                         function(f) is.na(f$fit_indices$AIC), logical(1))))

  # so the D2 pool of those statistics must withhold them as well ...
  expect_true(is.na(pooled$fit_indices$AIC))
  expect_true(is.na(pooled$fit_indices$BIC))
  expect_true(is.na(pooled$fit_indices$ECVI))

  # ... while the pooled chi-square and the indices that remain interpretable stay
  expect_true(is.finite(pooled$fit_indices$chi))
  expect_true(is.finite(pooled$fit_indices$CFI))
  expect_true(is.finite(pooled$fit_indices$RMSEA))
})

test_that("pooled AIC/BIC/ECVI are withheld when the components fell back to the plain LRT", {
  # Same fixture as above, but with the Stage-1 saturated covariance forced to fail, so no
  # component can form its correction and each reports the tagged plain two-stage
  # likelihood-ratio statistic instead. The pooler takes the components' withholding
  # decision, so the three criteria must stay NA on this route too -- the tag tells a reader
  # which statistic was pooled, and it must never say "corrected" here.
  imps <- lapply(1:3, function(i) {
    set.seed(100 + i)
    L <- matrix(0, 6, 2)
    L[1:3, 1] <- 0.7
    L[4:6, 2] <- 0.7
    S <- tcrossprod(L)
    diag(S) <- 1
    X <- matrix(stats::rnorm(400 * 6), 400) %*% chol(S)
    colnames(X) <- paste0("V", seq_len(6))
    X[X[, 1] > 0.8, 2] <- NA
    X[X[, 1] < -0.8, 3] <- NA
    X[X[, 1] > 1.2, 4] <- NA
    X
  })

  testthat::local_mocked_bindings(
    .fiml_saturated_acov = function(...) {
      cli::cli_abort("forced degenerate saturated covariance",
                     class = "efa_fiml_singular_information")
    }
  )

  pooled <- suppressMessages(suppressWarnings(
    efa_mi(imps, n_factors = 2, estimator = "ML", rotation = "none",
           cor_method = "fiml")
  ))

  expect_true(all(vapply(
    pooled$fits,
    function(f) identical(f$fit_indices$chi_scaled_type, "uncorrected.lrt"),
    logical(1))))
  expect_true(all(vapply(pooled$fits,
                         function(f) is.na(f$fit_indices$AIC), logical(1))))

  expect_true(is.na(pooled$fit_indices$AIC))
  expect_true(is.na(pooled$fit_indices$BIC))
  expect_true(is.na(pooled$fit_indices$ECVI))
  expect_true(is.finite(pooled$fit_indices$chi))
})

test_that("pooled AIC/BIC/ECVI are reported on an unscaled ML pool", {
  # The guard is conditional, not a blanket suppression: with plain (unscaled)
  # component statistics the descriptive criteria are still returned.
  pooled_ml <- efa_mi(cormat_list, n_factors = 3, N = 500, estimator = "ML",
                      rotation = "none")

  expect_true(all(vapply(pooled_ml$fits,
                         function(f) is.null(f$fit_indices$chi_scaled_type),
                         logical(1))))
  expect_true(is.finite(pooled_ml$fit_indices$AIC))
  expect_true(is.finite(pooled_ml$fit_indices$BIC))
  expect_true(is.finite(pooled_ml$fit_indices$ECVI))
})

test_that("a pooled RMSEA interval that misses its point estimate is reported as undefined", {
  # The same non-convergent noncentrality solve the single-fit path guards: at a chi in the
  # millions stats::pchisq() stops converging and both bounds collapse below the point
  # estimate. An interval that does not contain the RMSEA reported beside it is not an
  # interval, so it is withheld rather than shipped as a usable range.
  point_big <- .rmsea_point(4101811, 135, 5e6)
  ci_big <- .efa_pooled_rmsea_ci(4101811, df = 135, N = 5e6, point = point_big)
  expect_true(is.na(ci_big[["lower"]]))
  expect_true(is.na(ci_big[["upper"]]))

  # An ordinary pooled statistic keeps its interval, and that interval brackets the point
  # estimate the caller reports.
  point_ok <- .rmsea_point(200, 102, 500)
  ci_ok <- .efa_pooled_rmsea_ci(200, df = 102, N = 500, point = point_ok)
  expect_true(is.finite(ci_ok[["lower"]]))
  expect_true(is.finite(ci_ok[["upper"]]))
  expect_lte(ci_ok[["lower"]], point_ok)
  expect_lte(point_ok, ci_ok[["upper"]])

  # Containment is judged against the caller's point estimate, not one re-derived here from
  # chi and N: the two callers report estimates built on different statistics and scales.
  outside <- .efa_pooled_rmsea_ci(200, df = 102, N = 500,
                                  point = ci_ok[["upper"]] + .05)
  expect_true(is.na(outside[["lower"]]))
  expect_true(is.na(outside[["upper"]]))
})
