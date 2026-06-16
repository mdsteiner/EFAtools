# Tests for the non-parametric bootstrap standard error path of EFA()
# (se = "np-boot"): end-to-end coverage, output structure and validity,
# reproducibility, and graceful handling of degenerate bootstrap replicates.

# A clean oblique (promax) bootstrap fit reused by several structure/validity
# tests. GRiPS_raw is well conditioned, so all replicates should succeed.
set.seed(42)
boot_promax <- suppressWarnings(suppressMessages(
  EFA(GRiPS_raw, n_factors = 2, method = "PAF", rotation = "promax",
      se = "np-boot", b_boot = 12)
))

test_that("np-boot runs end to end for all methods and rotation families", {
  skip_on_cran()

  combos <- expand.grid(
    method = c("PAF", "ML", "ULS"),
    rotation = c("none", "varimax", "promax", "oblimin"),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(combos))) {
    method <- combos$method[i]
    rotation <- combos$rotation[i]
    label <- paste(method, rotation)

    set.seed(100 + i)
    res <- suppressWarnings(suppressMessages(
      EFA(GRiPS_raw, n_factors = 2, method = method, rotation = rotation,
          se = "np-boot", b_boot = 8)
    ))

    expect_s3_class(res, "EFA")
    expect_false(is.null(res$boot.SE), info = label)
    expect_false(is.null(res$boot.CI), info = label)
    expect_false(is.null(res$boot.arrays), info = label)
    # the third array dimension is the number of bootstrap replicates
    expect_identical(dim(res$boot.arrays$unrot_loadings)[3], 8L, info = label)
  }
})

test_that("np-boot resamples only the complete cases under listwise deletion", {
  skip_on_cran()

  # Build raw data whose complete cases sit at the tail: the first 35 rows each
  # carry a missing value and the 25 complete cases are rows 36-60. Under
  # use = "complete.obs" the correlation matrix - and hence N - rests on those
  # 25 complete cases, so the bootstrap must resample them. Resampling row
  # positions 1:N (= 1:25) would instead draw only all-missing rows, and
  # stats::cor(use = "complete.obs") then errors with "no complete element
  # pairs"; resampling the complete cases yields NA-free bootstrap correlations
  # and finite standard errors.
  set.seed(123)
  f <- rnorm(60)
  dat <- vapply(1:6, function(j) 0.6 * f + 0.8 * rnorm(60), numeric(60))
  colnames(dat) <- paste0("V", 1:6)
  dat[1:35, 1] <- NA

  expect_equal(sum(stats::complete.cases(dat)), 25)        # complete cases are rows 36-60
  expect_true(all(!stats::complete.cases(dat[1:25, ])))    # the old 1:N pool was all-missing

  set.seed(123)
  res <- suppressWarnings(suppressMessages(
    EFA(dat, n_factors = 1, method = "PAF", se = "np-boot", b_boot = 30,
        use = "complete.obs")
  ))

  expect_s3_class(res, "EFA")
  expect_true(all(is.finite(res$boot.SE$unrot_loadings)))
})

test_that("oblique np-boot output has the expected structure", {
  se <- boot_promax$boot.SE
  ci <- boot_promax$boot.CI
  arr <- boot_promax$boot.arrays

  expect_named(se, c("unrot_loadings", "rot_loadings", "Phi", "Structure",
                     "fit_indices", "residuals", "valid_target_rotations"))
  expect_named(ci, c("unrot_loadings", "rot_loadings", "Phi", "Structure",
                     "fit_indices", "residuals"))
  expect_named(arr, c("unrot_loadings", "rot_loadings", "Phi", "Structure",
                      "fit_indices", "residuals"))

  L <- boot_promax$rot_loadings
  expect_equal(dim(se$rot_loadings), dim(L))
  expect_equal(dim(se$unrot_loadings), dim(L))
  expect_equal(dim(se$Phi), c(ncol(L), ncol(L)))
  expect_identical(dim(arr$rot_loadings)[3], 12L)

  # the number of usable target rotations is reported and within bounds
  expect_true(se$valid_target_rotations >= 1 &&
                se$valid_target_rotations <= 12)
})

test_that("np-boot standard errors and confidence intervals are valid", {
  se <- boot_promax$boot.SE
  ci <- boot_promax$boot.CI

  for (nm in c("unrot_loadings", "rot_loadings", "Phi", "Structure")) {
    expect_true(all(se[[nm]] >= 0), info = nm)            # SEs are non-negative
    expect_true(all(is.finite(se[[nm]])), info = nm)
    expect_true(all(ci[[nm]]$lower <= ci[[nm]]$upper), info = nm)  # ordered CIs
  }

  # standardized residuals are added from the bootstrap residual SEs; the
  # off-diagonal entries (the ones of interest) are finite
  sr <- boot_promax$standardized_residuals
  expect_equal(dim(sr), dim(boot_promax$residuals))
  expect_true(all(is.finite(sr[upper.tri(sr)])))
})

test_that("np-boot is reproducible with a fixed seed", {
  skip_on_cran()

  run <- function() suppressWarnings(suppressMessages(
    EFA(GRiPS_raw, n_factors = 2, method = "PAF", rotation = "promax",
        se = "np-boot", b_boot = 8)
  ))

  set.seed(7); a <- run()
  set.seed(7); b <- run()

  expect_equal(a$boot.SE$unrot_loadings, b$boot.SE$unrot_loadings)
  expect_equal(a$boot.SE$rot_loadings, b$boot.SE$rot_loadings)
  expect_equal(a$boot.SE$Phi, b$boot.SE$Phi)
})

test_that("np-boot on a correlation matrix warns and disables the bootstrap", {
  expect_warning(
    res <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
               method = "PAF", rotation = "promax", se = "np-boot"),
    "Cannot compute bootstrap standard errors"
  )
  expect_s3_class(res, "EFA")
  expect_identical(res$settings$se, "none")
  expect_null(res$boot.SE)
})

test_that("a failed bootstrap replicate is skipped with a warning", {
  set.seed(11)
  x <- GRiPS_raw
  R <- stats::cor(x)
  N <- nrow(x)
  m <- ncol(R)
  b <- 6

  R_boot <- array(NA_real_, c(m, m, b))
  for (i in seq_len(b)) {
    ind <- sample(N, size = N, replace = TRUE)
    R_boot[, , i] <- stats::cor(x[ind, ])
  }

  fit_target <- suppressWarnings(
    .estimate_model(R, method = "PAF", n_factors = 2, N = N, type = "EFAtools"))
  boot_fit <- suppressWarnings(
    .boot_fun(R_boot, b, .estimate_model, method = "PAF", n_factors = 2,
              N = N, type = "EFAtools"))

  # inject a failed replicate without dropping the list element
  boot_fit[2] <- list(NULL)

  expect_warning(
    res <- .boot_se_ci(fit_target, L_rot = NULL, boot_fit,
                       boot_rot = "none", ci = 0.95, b = b),
    class = "efa_boot_replicate_failed"
  )

  # the failed replicate's slice stays NA; SEs from the rest are finite
  expect_true(all(is.na(res$arrays$unrot_loadings[, , 2])))
  expect_true(all(is.finite(res$SE$unrot_loadings)))
})

test_that(".boot_fun records a real failed replicate as NULL without dropping it", {
  # drive an actual fit failure (degenerate, non-finite correlation matrix)
  # through .boot_fun: a failed replicate must be recorded as a length-preserving
  # NULL, including when it is the LAST replicate, so .boot_se_ci can skip it.
  set.seed(13)
  x <- GRiPS_raw[1:200, ]
  R <- stats::cor(x)
  N <- nrow(x)
  m <- ncol(R)
  b <- 5

  bad <- R
  bad[1, ] <- NaN
  bad[, 1] <- NaN
  diag(bad) <- 1

  for (fail_at in c(2L, b)) {            # mid replicate and the last replicate
    R_boot <- array(NA_real_, c(m, m, b))
    for (i in seq_len(b)) {
      ind <- sample(N, size = N, replace = TRUE)
      R_boot[, , i] <- stats::cor(x[ind, ])
    }
    R_boot[, , fail_at] <- bad

    boot_fit <- suppressWarnings(
      .boot_fun(R_boot, b, .estimate_model, method = "PAF", n_factors = 2,
                N = N, type = "EFAtools"))

    # length preserved and exactly the degenerate replicate is NULL
    expect_length(boot_fit, b)
    expect_true(is.null(boot_fit[[fail_at]]))
    expect_equal(sum(vapply(boot_fit, is.null, logical(1))), 1L)

    # .boot_se_ci skips it gracefully (warns, does not error)
    fit_target <- suppressWarnings(
      .estimate_model(R, method = "PAF", n_factors = 2, N = N, type = "EFAtools"))
    expect_warning(
      res <- .boot_se_ci(fit_target, L_rot = NULL, boot_fit,
                         boot_rot = "none", ci = 0.95, b = b),
      class = "efa_boot_replicate_failed"
    )
    expect_true(all(is.finite(res$SE$unrot_loadings)))
  }
})

test_that("np-boot aborts when every replicate fails", {
  x <- GRiPS_raw[1:100, ]
  R <- stats::cor(x)
  b <- 4

  fit_target <- suppressWarnings(
    .estimate_model(R, method = "PAF", n_factors = 2, N = nrow(x),
                    type = "EFAtools"))
  boot_fit <- rep(list(NULL), b)

  expect_error(
    .boot_se_ci(fit_target, L_rot = NULL, boot_fit, boot_rot = "none",
                ci = 0.95, b = b),
    class = "efa_boot_all_failed"
  )
})

# Count how many signalled warnings carry a given condition class, muffling all
# warnings so the wrapped expression runs to completion.
count_warning_class <- function(expr, cls) {
  n <- 0L
  withCallingHandlers(
    force(expr),
    warning = function(w) {
      if (inherits(w, cls)) n <<- n + 1L
      invokeRestart("muffleWarning")
    }
  )
  n
}

test_that("a pinned argument re-warns once, not once per bootstrap replicate", {
  skip_on_cran()
  # `efa_type_override` reflects the (type, pinned-argument) combination, which is
  # identical for every replicate and already surfaced once by the point-estimate
  # fit. The bootstrap loop must not repeat it b_boot times.
  set.seed(42)
  n_override <- count_warning_class(
    EFA(GRiPS_raw, n_factors = 2, method = "PAF", type = "EFAtools",
        max_iter = 500, se = "np-boot", b_boot = 6),
    "efa_type_override"
  )
  expect_equal(n_override, 1L)
})

test_that("bootstrap non-convergence is summarized in a single classed warning", {
  skip_on_cran()
  # A tiny iteration cap forces every replicate to hit the maximum-iteration limit;
  # the per-replicate fitter warnings must be suppressed and replaced by a single
  # classed summary rather than one warning per replicate.
  set.seed(1)
  n_summary <- count_warning_class(
    EFA(GRiPS_raw, n_factors = 2, method = "PAF", type = "none",
        init_comm = "smc", criterion = 1e-3, criterion_type = "sum",
        abs_eigen = TRUE, max_iter = 1, se = "np-boot", b_boot = 5),
    "efa_boot_nonconvergence"
  )
  expect_equal(n_summary, 1L)

  # a cleanly converging bootstrap emits no non-convergence summary
  set.seed(7)
  n_clean <- count_warning_class(
    EFA(GRiPS_raw, n_factors = 2, method = "PAF", rotation = "promax",
        se = "np-boot", b_boot = 8),
    "efa_boot_nonconvergence"
  )
  expect_equal(n_clean, 0L)
})

test_that("eigendecomposition guards turn degenerate matrices into errors", {
  # a non-finite (constant-column-style) correlation matrix makes the symmetric
  # eigendecomposition fail; the guarded fitters must error, not crash R
  m <- 6
  R_bad <- diag(m)
  R_bad[1, ] <- NaN
  R_bad[, 1] <- NaN
  diag(R_bad) <- 1
  psi <- rep(0.5, m)

  expect_error(.paf_iter(psi, 0.001, R_bad, 2L, TRUE, 2L, 100L),
               "Eigendecomposition failed")
  expect_error(.grad_ml(psi, R_bad, 2L), "Eigendecomposition failed")
  expect_error(.error_ml(psi, R_bad, 2L), "Eigendecomposition failed")
  expect_error(.grad_uls(psi, R_bad, 2L), "Eigendecomposition failed")
  expect_error(.uls_residuals(psi, R_bad, 2L), "Eigendecomposition failed")
})

test_that("over-extraction guards turn n_fac >= ncol into errors", {
  # the eigenvalue-based extraction reads the largest n_fac eigenpairs; with
  # n_fac >= ncol(R) it would index past the available eigenvalues (undefined
  # behaviour in an unchecked build). The guarded fitters must error, not crash R
  m <- 6L
  R <- diag(m)
  psi <- rep(0.5, m)

  expect_error(.paf_iter(rep(1, m), 0.001, R, m, TRUE, 2L, 10L),
               "smaller than the number of variables")
  expect_error(.grad_ml(psi, R, m), "smaller than the number of variables")
  expect_error(.error_ml(psi, R, m), "smaller than the number of variables")
  expect_error(.grad_uls(psi, R, m), "smaller than the number of variables")
  expect_error(.uls_residuals(psi, R, m), "smaller than the number of variables")
})

test_that(".estimate_model(lean = TRUE) returns only the bootstrap-aggregated quantities", {
  # The bootstrap replicate fitter computes just the loadings, fit indices, and
  # residuals that .boot_se_ci() aggregates. Those must match the full fit
  # exactly, except for the analytic RMSEA bounds, which the lean fit does not
  # solve (they are not meaningfully bootstrapped).
  R <- test_models$baseline$cormat
  N <- 500
  bounds <- c("RMSEA_LB", "RMSEA_UB")

  method_args <- list(
    PAF = list(type = "EFAtools"),
    ML  = list(start_method = "psych"),
    ULS = list()
  )

  for (method in names(method_args)) {
    common <- c(list(R, method = method, n_factors = 3, N = N),
                method_args[[method]])
    full <- suppressWarnings(do.call(.estimate_model, common))
    lean <- suppressWarnings(do.call(.estimate_model, c(common, list(lean = TRUE))))

    expect_named(lean, c("unrot_loadings", "fit_indices", "residuals", "convergence"))

    expect_equal(as.vector(lean$unrot_loadings), as.vector(full$unrot_loadings),
                 info = method)
    expect_equal(as.vector(lean$residuals), as.vector(full$residuals), info = method)
    expect_identical(lean$convergence, full$convergence, info = method)

    keep <- setdiff(names(full$fit_indices), bounds)
    expect_equal(lean$fit_indices[keep], full$fit_indices[keep], info = method)
    expect_true(all(is.na(unlist(lean$fit_indices[bounds]))), info = method)
  }
})

test_that("ML np-boot drops only the analytic RMSEA bounds from the fit-index SEs", {
  skip_on_cran()
  set.seed(202)
  res <- suppressWarnings(suppressMessages(
    EFA(GRiPS_raw, n_factors = 2, method = "ML", rotation = "none",
        se = "np-boot", b_boot = 10)
  ))

  se_fit <- res$boot.SE$fit_indices
  # the per-replicate analytic RMSEA bounds are not bootstrapped
  expect_true(all(is.na(se_fit[c("RMSEA_LB", "RMSEA_UB")])))
  # the bootstrapped fit indices that are aggregated stay finite
  expect_true(all(is.finite(se_fit[c("CAF", "RMSR", "SRMR")])))
  # the point estimate keeps its full analytic RMSEA confidence interval
  expect_true(is.finite(res$fit_indices$RMSEA_LB))
  expect_true(is.finite(res$fit_indices$RMSEA_UB))
})
