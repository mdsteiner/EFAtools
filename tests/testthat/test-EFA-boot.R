# Tests for the non-parametric bootstrap standard error path of EFA()
# (se = "np-boot"): end-to-end coverage, output structure and validity,
# reproducibility, and graceful handling of degenerate bootstrap replicates.

# The reproducibility assertions below compare two separate invocations, so they use
# expect_equal() with an explicit tolerance rather than expect_identical(): a threaded BLAS
# (Apple's Accelerate, for one) is free to vary its GEMM reduction order between calls, so
# one function on one input can differ in the last ulp, and an iterative rotation carries
# that into the reported loadings. waldo still compares S3 classes, names, attributes and
# structure exactly, so only numeric values are given slack -- a seed that was ignored moves
# the compared values by orders of magnitude more than this and still fails loudly. (A set
# tolerance does relax integer against double, so where the storage mode is itself part of
# the contract it is asserted separately.)
fp_tol <- 1e-8

# A clean oblique (promax) bootstrap fit reused by several structure/validity
# tests. GRiPS_raw is well conditioned, so all replicates should succeed.
set.seed(42)
boot_promax <- suppressWarnings(suppressMessages(
  EFA(GRiPS_raw, n_factors = 2, method = "PAF", rotation = "promax",
      se = "np-boot", b_boot = 12)
))

test_that("np-boot runs end to end for all methods and rotation families", {
  skip_on_cran()

  # Three (method, rotation) representatives -- one per method and one per rotation
  # family (none / orthogonal / oblique) -- guard the method x np-boot dispatch
  # without redoing the full method x rotation cross product. The rotation-family
  # slot schema is independently pinned across se = information / sandwich / np-boot
  # by test-EFA-fields.R.
  combos <- data.frame(
    method   = c("PAF", "ML", "ULS"),
    rotation = c("promax", "varimax", "none"),
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
    expect_false(is.null(res$SE), info = label)
    expect_false(is.null(res$CI), info = label)
    expect_false(is.null(res$replicates), info = label)
    # the third array dimension is the number of bootstrap replicates
    expect_identical(dim(res$replicates$unrot_loadings)[3], 8L, info = label)
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
  expect_true(all(is.finite(res$SE$unrot_loadings)))
})

test_that("oblique np-boot output has the expected structure", {
  se <- boot_promax$SE
  ci <- boot_promax$CI
  arr <- boot_promax$replicates

  expect_named(se, c("unrot_loadings", "rot_loadings", "Phi", "Structure",
                     "fit_indices", "residuals", "valid_replicates",
                     "valid_target_rotations"))
  expect_named(ci, c("unrot_loadings", "rot_loadings", "Phi", "Structure",
                     "fit_indices", "residuals"))
  expect_named(arr, c("unrot_loadings", "rot_loadings", "Phi", "Structure",
                      "fit_indices", "residuals"))

  L <- boot_promax$rot_loadings
  expect_equal(dim(se$rot_loadings), dim(L))
  expect_equal(dim(se$unrot_loadings), dim(L))
  expect_equal(dim(se$Phi), c(ncol(L), ncol(L)))
  expect_identical(dim(arr$rot_loadings)[3], 12L)

  # the effective B, and the number of usable target rotations, are reported and in bounds
  expect_true(se$valid_replicates >= 1 && se$valid_replicates <= 12)
  expect_true(se$valid_target_rotations >= 1 &&
                se$valid_target_rotations <= se$valid_replicates)
})

test_that("the printed bootstrap sample count reports the usable replicates when they differ", {
  # The effective B is what the unrotated-loading, residual and fit-index intervals rest on. Without
  # it the output asserts the requested b_boot whatever the survival rate, so an interval built from
  # 4 order statistics reads as one built from 20 -- and states a number that is wrong.
  degraded <- boot_promax
  degraded$SE$valid_replicates <- 4L

  body <- cli::ansi_strip(format(summary(degraded)))
  expect_true(any(grepl("Bootstrap samples: 12 (4 usable)", body, fixed = TRUE)))
  expect_true(any(grepl("12 bootstrap samples (4 usable)", body, fixed = TRUE)))

  # A run in which every replicate survived says nothing extra.
  clean <- cli::ansi_strip(format(summary(boot_promax)))
  expect_false(any(grepl("usable", clean, fixed = TRUE)))
})

test_that("np-boot standard errors and confidence intervals are valid", {
  se <- boot_promax$SE
  ci <- boot_promax$CI

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

  expect_equal(a$SE$unrot_loadings, b$SE$unrot_loadings)
  expect_equal(a$SE$rot_loadings, b$SE$rot_loadings)
  expect_equal(a$SE$Phi, b$SE$Phi)
})

test_that("np-boot with a fixed seed is reproducible at 1 vs 2 workers", {
  skip_on_cran()
  skip_if_not_slow()
  # The replicate fits are parallelised across workers with future.apply, and
  # future.seed = TRUE binds each replicate's RNG stream to its index. With a fixed
  # `seed` the bootstrap must therefore return the same result regardless of the
  # number of workers. The comparison uses a small tolerance rather than bit-for-bit
  # equality: the worker fits run in separate processes whose BLAS/LAPACK may sum in a
  # different order, so results can differ in the last bit or two while remaining
  # numerically equivalent. The multisession workers are fresh R processes that load
  # the installed package, so run this under devtools::check() / after
  # devtools::install() for the worker code to match the main process. (multicore is
  # unavailable on Windows.)
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  run <- function() suppressWarnings(suppressMessages(
    EFA(GRiPS_raw, n_factors = 2, method = "PAF", rotation = "promax",
        se = "np-boot", b_boot = 12, seed = 2024)
  ))

  future::plan(future::sequential)
  one <- run()

  future::plan(future::multisession, workers = 2)
  two <- run()

  expect_equal(one$replicates, two$replicates, tolerance = 1e-10)
  expect_equal(one$SE, two$SE, tolerance = 1e-10)
  expect_equal(one$CI, two$CI, tolerance = 1e-10)
})

test_that("a supplied seed leaves the caller's RNG stream unchanged", {
  # Passing `seed` must not have a lasting side effect on the global RNG: the stream
  # is restored on exit, so a draw taken after a seeded bootstrap is identical to one
  # taken without the call having happened.
  set.seed(1)
  state_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)

  # An oblique rotation exercises every RNG consumer the restore must neutralize:
  # the case resampling, the point-estimate rotation random starts, and the oblique
  # Procrustes random starts that run after the replicate fits.
  invisible(suppressWarnings(suppressMessages(
    EFA(GRiPS_raw, n_factors = 2, method = "PAF", rotation = "promax",
        se = "np-boot", b_boot = 6, seed = 99)
  )))

  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   state_before)
})

test_that("seed makes a criterion rotation reproducible without a bootstrap", {
  # `seed` governs the whole fit, not just the bootstrap: the criterion-based rotations
  # draw random starts from the RNG, so two seeded calls at the default se = "none" must
  # agree. simplimax is used because its criterion is only piecewise smooth, making it
  # the most prone to landing in different optima from different starts -- unseeded, two
  # calls disagree by a wide margin, so this would fail if `seed` were ignored.
  fit <- function(s) suppressWarnings(
    efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500, estimator = "PAF",
            rotation = "simplimax", seed = s))

  a <- fit(42)
  b <- fit(42)
  expect_equal(unclass(a$rot_loadings), unclass(b$rot_loadings), tolerance = fp_tol)
  expect_equal(a$Phi, b$Phi, tolerance = fp_tol)

  # A different seed explores different starts, so the argument is doing something
  # rather than the fit being deterministic to begin with.
  set.seed(1)
  unseeded_1 <- suppressWarnings(
    efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500, estimator = "PAF",
            rotation = "simplimax"))
  set.seed(2)
  unseeded_2 <- suppressWarnings(
    efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500, estimator = "PAF",
            rotation = "simplimax"))
  expect_false(identical(unclass(unseeded_1$rot_loadings),
                         unclass(unseeded_2$rot_loadings)))
})

test_that("seed leaves the caller's RNG stream unchanged without a bootstrap", {
  # The restore-on-exit contract holds on the non-bootstrap path too.
  set.seed(1)
  state_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)

  invisible(suppressWarnings(
    efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500, estimator = "PAF",
            rotation = "simplimax", seed = 99)))

  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   state_before)
})

test_that("np-boot on a correlation matrix is an error", {
  expect_error(
    efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500, se = "np-boot"),
    class = "efa_boot_cormat"
  )
})

# A healthy PAF point-estimate fit on GRiPS_raw together with the matching list of `b`
# replicate fits. The failure-mode blocks below all start from this pair and damage it in
# different ways, so it is built by one helper rather than retyped. The seed lives inside the
# helper so every call draws the same resamples from the same stream whatever ran before it.
make_boot_pair <- function(b, n_factors) {
  set.seed(11)
  x <- GRiPS_raw
  R <- stats::cor(x)
  N <- nrow(x)
  m <- ncol(R)

  R_boot <- array(NA_real_, c(m, m, b))
  for (i in seq_len(b)) {
    ind <- sample(N, size = N, replace = TRUE)
    R_boot[, , i] <- stats::cor(x[ind, ])
  }

  list(
    fit_target = suppressWarnings(
      .estimate_model(R, method = "PAF", n_factors = n_factors, N = N, type = "EFAtools")),
    boot_fit = suppressWarnings(
      .boot_fun(function(i) list(R = R_boot[, , i]), seq_len(b), .estimate_model,
                method = "PAF", n_factors = n_factors, N = N, type = "EFAtools"))
  )
}

test_that("a failed bootstrap replicate is skipped with a warning", {
  b <- 6
  pair <- make_boot_pair(b, n_factors = 2)
  boot_fit <- pair$boot_fit

  # inject a failed replicate without dropping the list element
  boot_fit[2] <- list(NULL)

  expect_warning(
    res <- .boot_se_ci(pair$fit_target, L_rot = NULL, boot_fit,
                       boot_rot = "none", ci = 0.95, b = b),
    class = "efa_boot_replicate_failed"
  )

  # the failed replicate's slice stays NA; SEs from the rest are finite
  expect_true(all(is.na(res$replicates$unrot_loadings[, , 2])))
  expect_true(all(is.finite(res$SE$unrot_loadings)))
  # the effective B is recorded on the object, not only in the warning
  expect_equal(res$SE$valid_replicates, b - 1)
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
      .boot_fun(function(i) list(R = R_boot[, , i]), seq_len(b), .estimate_model,
                method = "PAF", n_factors = 2, N = N, type = "EFAtools"))

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

test_that("a replicate that cannot be target-rotated warns with its own class", {
  # The rotated block has a survival count of its own: a replicate can fit and align to the
  # unrotated point estimate and still fail the target rotation. A non-finite target forces that on
  # the orthogonal branch; the oblique branch signals the same class from its batched aligner.
  b <- 4
  pair <- make_boot_pair(b, n_factors = 2)

  bad_target <- unclass(pair$fit_target$unrot_loadings)
  bad_target[1, 1] <- NaN

  expect_warning(
    res <- .boot_se_ci(pair$fit_target, L_rot = bad_target, pair$boot_fit,
                       boot_rot = "orthogonal", ci = 0.95, b = b),
    class = "efa_boot_rotation_failed"
  )
  expect_equal(res$SE$valid_target_rotations, 0)
  expect_true(all(is.na(res$SE$rot_loadings)))
  # the unrotated block is unaffected: only the target rotation failed
  expect_equal(res$SE$valid_replicates, b)
  expect_true(all(is.finite(res$SE$unrot_loadings)))
})

test_that("b_boot below two is rejected", {
  # A bootstrap standard error is the dispersion across replicates: at b_boot = 1 every SE is the
  # sd() of a single value and comes back NA, and the interval collapses onto that replicate. That
  # used to happen silently -- the only route on which an SE returned NA with no condition at all.
  expect_error(
    efa_fit(GRiPS_raw, n_factors = 1, estimator = "PAF", rotation = "none",
            se = "np-boot", b_boot = 1),
    class = "efa_b_boot_too_small"
  )
  expect_error(
    efa_fit(GRiPS_raw, n_factors = 1, estimator = "PAF", rotation = "none",
            se = "np-boot", b_boot = 0),
    class = "efa_b_boot_too_small"
  )
  expect_error(
    efa_fit(GRiPS_raw, n_factors = 1, estimator = "PAF", rotation = "none",
            se = "np-boot", b_boot = -5),
    class = "efa_b_boot_too_small"
  )
})

test_that("a ci at either end of the unit interval is rejected", {
  # The bound is open at both ends, on every SE path. The analytic paths read the level as
  # z = qnorm(1 - (1 - ci) / 2): at ci = 1 that is infinite and each Wald bound came back
  # -Inf/Inf, at ci = 0 it is zero and each interval collapsed onto the point estimate. Both
  # passed silently, because the assertion took the closed interval the documentation described.
  cm <- test_models$baseline$cormat

  expect_error(
    efa_fit(cm, n_factors = 3, N = 500, estimator = "ML", se = "information", ci = 1),
    class = "efa_ci_out_of_bounds"
  )
  expect_error(
    efa_fit(cm, n_factors = 3, N = 500, estimator = "ML", se = "information", ci = 0),
    class = "efa_ci_out_of_bounds"
  )
  # The guard sits in the argument checks, so it fires whatever the SE path -- including none.
  expect_error(
    efa_fit(cm, n_factors = 3, N = 500, estimator = "PAF", se = "np-boot", ci = 1),
    class = "efa_ci_out_of_bounds"
  )
  # It owns the whole range, so a level outside the unit interval reports the same rule.
  expect_error(
    efa_fit(cm, n_factors = 3, N = 500, estimator = "ML", se = "information", ci = 1.5),
    class = "efa_ci_out_of_bounds"
  )
  # A non-number stays with the shared argument assertion.
  expect_error(
    efa_fit(cm, n_factors = 3, N = 500, estimator = "ML", se = "information", ci = NA),
    class = "efa_invalid_argument"
  )

  # A level inside the interval still passes.
  expect_no_error(
    efa_fit(cm, n_factors = 3, N = 500, estimator = "ML", rotation = "none",
            se = "information", ci = .99)
  )
})

test_that("fewer than two surviving replicates is flagged as unreliable", {
  # The input bound cannot see this one: b_boot is legal, but the replicates fail at run time and
  # leave a single usable draw behind. Without a condition the object returns all-NA SEs and a
  # collapsed interval that print as though the requested b_boot had stood.
  b <- 6
  pair <- make_boot_pair(b, n_factors = 1)
  fit_target <- pair$fit_target
  boot_fit <- pair$boot_fit

  one_left <- boot_fit
  one_left[2:b] <- rep(list(NULL), b - 1L)     # a single survivor

  expect_warning(
    expect_warning(
      res <- .boot_se_ci(fit_target, L_rot = NULL, one_left, boot_rot = "none",
                         ci = 0.95, b = b),
      class = "efa_boot_replicate_failed"
    ),
    class = "efa_se_unreliable"
  )
  expect_equal(res$SE$valid_replicates, 1L)
  expect_true(all(is.na(res$SE$unrot_loadings)))
  expect_equal(res$CI$unrot_loadings$lower, res$CI$unrot_loadings$upper)

  # Two survivors are enough for a defined (if very wide) standard error: the failed-replicate
  # warning still fires there, but the unreliability warning must not.
  two_left <- boot_fit
  two_left[3:b] <- rep(list(NULL), b - 2L)
  seen <- character(0)
  withCallingHandlers(
    res2 <- .boot_se_ci(fit_target, L_rot = NULL, two_left, boot_rot = "none",
                        ci = 0.95, b = b),
    warning = function(w) {
      seen <<- c(seen, class(w)[1])
      invokeRestart("muffleWarning")
    }
  )
  expect_false("efa_se_unreliable" %in% seen)
  expect_equal(res2$SE$valid_replicates, 2L)
  expect_true(all(is.finite(res2$SE$unrot_loadings)))
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

  expect_eigendecomp_error <- function(expr) {
    # Armadillo writes a diagnostic to stderr before returning failure; capture it
    # so the test asserts the classed R error without noisy logs.
    invisible(utils::capture.output(
      expect_error(force(expr), "Eigendecomposition failed"),
      type = "message"
    ))
  }

  expect_eigendecomp_error(.paf_iter(psi, 0.001, R_bad, 2L, TRUE, 2L, 100L))
  expect_eigendecomp_error(.grad_ml(psi, R_bad, 2L))
  expect_eigendecomp_error(.error_ml(psi, R_bad, 2L))
  expect_eigendecomp_error(.grad_uls(psi, R_bad, 2L))
  expect_eigendecomp_error(.uls_residuals(psi, R_bad, 2L))
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

test_that("the PAF kernel rejects a non-positive iteration budget", {
  # With max_iter < 1 the iteration never runs and the kernel would return an empty
  # loading matrix, which .finalize_fit() would read as a zero-factor solution. The
  # R-side control validation is the normal user-facing path and keeps its condition
  # class; this asserts the backstop on the other side of the boundary.
  m <- 6L
  R <- matrix(0.3, m, m)
  diag(R) <- 1
  h2 <- rep(0.5, m)

  expect_error(.paf_iter(h2, 0.001, R, 2L, TRUE, 2L, 0L), "at least 1")
  expect_error(.paf_iter(h2, 0.001, R, 2L, TRUE, 2L, -5L), "at least 1")

  expect_error(estimate_control(max_iter = 0), class = "efa_control_input")
  expect_error(estimate_control(max_iter = -5), class = "efa_control_input")
})

test_that(".array_se_ci matches a per-probability sweep at every shape it is called with", {
  # The helper takes both percentile bounds from one pass over the replicate array. The bounds
  # must equal what a separate sweep per probability produced -- same values, same shapes, same
  # dimnames -- at every rank its call sites use: a replicate matrix (the fit indices, M = 2),
  # a 3-D cube (loadings, Phi, residuals, M = c(1, 2)), and the group bootstrap's 4-D cube.
  ref <- function(x, probs, M) {
    ci <- lapply(probs, function(p) apply(x, M, stats::quantile, probs = p, na.rm = TRUE))
    stats::setNames(ci, c("lower", "upper"))
  }
  probs <- c(0.025, 0.975)

  set.seed(20)
  cases <- list(
    matrix_M2 = list(
      x = matrix(stats::rnorm(40 * 5), nrow = 40,
                 dimnames = list(NULL, paste0("idx", 1:5))),
      M = 2),
    cube_M12 = list(
      x = array(stats::rnorm(6 * 3 * 40), c(6, 3, 40),
                dimnames = list(paste0("V", 1:6), paste0("F", 1:3), NULL)),
      M = c(1, 2)),
    cube4d_M123 = list(
      x = array(stats::rnorm(2 * 2 * 4 * 30), c(2, 2, 4, 30),
                dimnames = list(c("g1", "g2"), c("g1", "g2"), paste0("V", 1:4), NULL)),
      M = c(1L, 2L, 3L))
  )

  for (nm in names(cases)) {
    x <- cases[[nm]]$x
    M <- cases[[nm]]$M
    got <- .array_se_ci(x, probs, M = M)
    want <- ref(x, probs, M)
    expect_equal(got$ci$lower, want$lower, info = nm)
    expect_equal(got$ci$upper, want$upper, info = nm)
    expect_identical(dim(got$ci$lower), dim(got$se), info = nm)
    expect_identical(dimnames(got$ci$lower), dimnames(got$se), info = nm)
  }

  # NA handling is unchanged too: a cell whose replicates are all missing must come back NA
  # rather than error, and a partially missing cell must use its observed replicates only.
  y <- array(stats::rnorm(3 * 2 * 20), c(3, 2, 20))
  y[1, 1, ] <- NA_real_
  y[2, 1, 1:5] <- NA_real_
  got <- .array_se_ci(y, probs)
  want <- ref(y, probs, c(1, 2))
  expect_equal(got$ci$lower, want$lower)
  expect_equal(got$ci$upper, want$upper)
  expect_true(is.na(got$ci$lower[1, 1]))
  expect_false(is.na(got$ci$lower[2, 1]))
})

test_that(".estimate_model(lean = TRUE) returns only the bootstrap-aggregated quantities", {
  # The bootstrap replicate fitter computes just the loadings and residuals that
  # .boot_se_ci() aggregates and the fit indices that .yhy_fit_boot() aggregates.
  # Those must match the full fit exactly, except for the analytic RMSEA bounds,
  # which the lean fit does not solve (they are not meaningfully bootstrapped).
  # `gof = FALSE` drops the fit indices only.
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

    no_gof <- suppressWarnings(do.call(.estimate_model,
                                       c(common, list(lean = TRUE, gof = FALSE))))
    expect_null(no_gof$fit_indices)
    expect_identical(no_gof[c("unrot_loadings", "residuals", "convergence")],
                     lean[c("unrot_loadings", "residuals", "convergence")])
  }
})

test_that("np-boot gives bootstrap SEs and intervals of RMSEA, CFI, and TLI only", {
  skip_on_cran()
  res <- suppressWarnings(
    efa_fit(GRiPS_raw, n_factors = 2, estimator = "ML", rotation = "none",
            se = "np-boot", b_boot = 10, seed = 202)
  )

  keep <- c("RMSEA", "CFI", "TLI")
  se_fit <- res$SE$fit_indices
  ci_fit <- res$CI$fit_indices
  expect_named(se_fit, keep)
  expect_named(ci_fit$lower, keep)
  expect_named(ci_fit$upper, keep)
  expect_true(all(is.finite(se_fit)))
  expect_true(all(ci_fit$lower <= ci_fit$upper))
  # the replicate matrix labels its columns, and the replicate is its FIRST dimension
  expect_identical(colnames(res$replicates$fit_indices), keep)
  expect_identical(nrow(res$replicates$fit_indices), 10L)
  # the point estimate keeps its full analytic RMSEA confidence interval
  expect_true(is.finite(res$fit_indices$RMSEA_LB))
  expect_true(is.finite(res$fit_indices$RMSEA_UB))
})

test_that("print() and summary() show the bootstrap intervals of CFI and TLI", {
  ci_fit <- boot_promax$CI$fit_indices
  expect_true(all(is.finite(c(ci_fit$lower, ci_fit$upper))))
  for (out in list(format(boot_promax), format(summary(boot_promax)))) {
    out <- cli::ansi_strip(out)
    expect_true(any(startsWith(out, "CFI [95% bootstrap CI]: ")))
    expect_true(any(startsWith(out, "TLI [95% bootstrap CI]: ")))
    # RMSEA keeps its analytic interval
    expect_true(any(startsWith(out, "RMSEA [90% CI]: ")))
    expect_true(any(grepl("^SRMR: ", out)))
    expect_false(any(grepl("^(RMSR|CAF)\\b", out)))
  }
  # summary(ci = "none") shows no interval at all
  out <- cli::ansi_strip(format(summary(boot_promax, ci = "none")))
  expect_false(any(grepl("bootstrap CI", out, fixed = TRUE)))
  expect_true(any(startsWith(out, "CFI: ")))
})

# A 2-factor model fitted to a 3-factor population with a weak third factor: RMSEA .057 and
# CFI .955 in the population. Shared by the fit-index bootstrap tests below.
yhy_pop <- local({
  L <- matrix(0, 12, 3)
  L[1:5, 1] <- c(.7, .65, .6, .7, .55)
  L[6:10, 2] <- c(.6, .7, .65, .55, .6)
  L[c(3, 5, 8, 11, 12), 3] <- c(.3, .3, .3, .5, .5)
  Phi <- matrix(.3, 3, 3)
  diag(Phi) <- 1
  R <- L %*% Phi %*% t(L)
  diag(R) <- 1
  dimnames(R) <- list(paste0("V", 1:12), paste0("V", 1:12))
  R
})
yhy_sample <- function(N, seed) {
  set.seed(seed)
  x <- matrix(stats::rnorm(N * ncol(yhy_pop)), N) %*% chol(yhy_pop)
  colnames(x) <- colnames(yhy_pop)
  x
}

test_that(".yhy_transform hits the target misfit and keeps the fitted solution", {
  x <- yhy_sample(300, 1)
  R <- stats::cor(x)
  N <- nrow(x)
  for (method in c("ML", "ULS", "PAF")) {
    fit_args <- list(method = method, n_factors = 2, N = N, type = "EFAtools",
                     max_iter = 1000, init_comm = "smc", criterion = 1e-8,
                     criterion_type = "sum", abs_eigen = TRUE, start_method = "psych")
    fit <- do.call(.estimate_model, c(list(R), fit_args))
    R_hat <- fit$model_implied_R
    df <- fit$fit_indices$df
    tau <- .ml_discrepancy(R_hat, R) - df / (N - 1)
    expect_gt(tau, 0)

    yhy <- .yhy_transform(scale(x), R_hat, N, df)
    R_a <- yhy$a * R + (1 - yhy$a) * R_hat
    expect_true(yhy$a > 0 && yhy$a < 1, info = method)
    # The transformed cases have the correlation matrix R_a, whose ML discrepancy from the
    # model-implied matrix is the estimated misfit.
    expect_equal(stats::cor(yhy$data), R_a, tolerance = 1e-10, ignore_attr = TRUE,
                 info = method)
    expect_equal(.ml_discrepancy(R_hat, R_a), tau, tolerance = 1e-8, info = method)
    # The fitted loadings also fit R_a best (Yuan, Hayashi, & Yanagihara, 2007, Theorem 1).
    fit_a <- do.call(.estimate_model, c(list(R_a), fit_args))
    expect_equal(unclass(fit_a$unrot_loadings), unclass(fit$unrot_loadings),
                 tolerance = 1e-4, ignore_attr = TRUE, info = method)
  }

  # No estimated misfit (tau = 0): a = 0, the transformed cases fit the model exactly.
  fit <- .estimate_model(R, method = "ML", n_factors = 2, N = N, start_method = "psych")
  R_hat <- fit$model_implied_R
  yhy <- .yhy_transform(scale(x), R_hat, N = 3, df = fit$fit_indices$df)
  expect_identical(yhy$a, 0)
  expect_equal(stats::cor(yhy$data), R_hat, tolerance = 1e-10, ignore_attr = TRUE)

  # A model-implied matrix that is not positive definite has no transformation.
  R_bad <- R_hat
  R_bad[1, 2] <- R_bad[2, 1] <- 1.5
  expect_null(.yhy_transform(scale(x), R_bad, N, 1))
  # Nor do cases with a singular correlation matrix (a duplicated variable).
  expect_null(.yhy_transform(scale(cbind(x, x[, 1])), diag(13), N, 1))
})

test_that("the fit-index bootstrap changes no other result and no random-number stream", {
  skip_on_cran()
  # Its fits draw no random numbers, so without `seed` the caller's stream after the fit, and
  # all other results, are the same as without it.
  x <- yhy_sample(200, 2)
  run <- function() {
    set.seed(3)
    fit <- efa_fit(x, n_factors = 2, se = "np-boot", b_boot = 20)
    list(fit = fit, next_draw = stats::runif(1))
  }
  with_yhy <- run()
  .efa_state$no_fit_boot <- TRUE
  on.exit(.efa_state$no_fit_boot <- FALSE, add = TRUE)
  without <- run()

  expect_identical(with_yhy$next_draw, without$next_draw)
  expect_null(without$fit$SE$fit_indices)
  expect_named(with_yhy$fit$SE$fit_indices, c("RMSEA", "CFI", "TLI"))
  with_yhy$fit$SE$fit_indices <- NULL
  with_yhy$fit$CI$fit_indices <- NULL
  with_yhy$fit$replicates$fit_indices <- NULL
  expect_equal(with_yhy$fit, without$fit, tolerance = fp_tol)
})

test_that("np-boot gives no fit-index intervals outside the supported data", {
  skip_on_cran()
  no_fit_boot <- function(fit) {
    expect_null(fit$SE$fit_indices)
    expect_null(fit$CI$fit_indices)
    expect_null(fit$replicates$fit_indices)
  }
  x_na <- GRiPS_raw
  x_na[1:5, 1] <- NA
  # pairwise deletion of missing values
  no_fit_boot(suppressWarnings(
    efa_fit(x_na, n_factors = 1, se = "np-boot", b_boot = 4, seed = 1)))
  # a rank correlation
  no_fit_boot(efa_fit(GRiPS_raw, n_factors = 1, cor_method = "spearman", se = "np-boot",
                      b_boot = 4, seed = 1))
  # polychoric correlations
  no_fit_boot(suppressWarnings(
    efa_fit(DOSPERT_raw[, 1:5], n_factors = 1, cor_method = "poly", se = "np-boot",
            b_boot = 4, seed = 1)))
  # a singular correlation matrix (a duplicated variable, which the psych preset accepts)
  no_fit_boot(suppressMessages(suppressWarnings(
    efa_fit(cbind(GRiPS_raw, copy = GRiPS_raw[, 1]), n_factors = 1, estimator = "PAF",
            estimate_control = estimate_control(type = "psych"), se = "np-boot",
            b_boot = 4, seed = 1))))
  # a just-identified model (df = 0) has no misfit to transform to
  no_fit_boot(suppressWarnings(
    efa_fit(GRiPS_raw[, 1:3], n_factors = 1, se = "np-boot", b_boot = 4, seed = 1)))
  # listwise deletion of the same missing values is supported
  listwise <- efa_fit(x_na, n_factors = 1, use = "complete.obs", se = "np-boot",
                      b_boot = 4, seed = 1)
  expect_named(listwise$CI$fit_indices$lower, c("RMSEA", "CFI", "TLI"))
})

test_that("the fit-index bootstrap intervals have approximately correct coverage", {
  skip_on_cran()
  skip_if_not_slow()
  # The population values of the 2-factor model, from its fit to the population matrix.
  L_pop <- .estimate_model(yhy_pop, method = "ML", n_factors = 2,
                           start_method = "psych")$unrot_loadings
  pop <- .efa_population_fit(unclass(L_pop), yhy_pop)
  covered <- vapply(1:60, function(i) {
    fit <- efa_fit(yhy_sample(300, 100 + i), n_factors = 2, estimator = "ML",
                   rotation = "none", se = "np-boot", b_boot = 100, seed = i)
    ci_fit <- fit$CI$fit_indices
    c(RMSEA = ci_fit$lower[["RMSEA"]] <= pop$rmsea && pop$rmsea <= ci_fit$upper[["RMSEA"]],
      CFI = ci_fit$lower[["CFI"]] <= pop$cfi && pop$cfi <= ci_fit$upper[["CFI"]])
  }, logical(2))
  # Nominal 95%; about 94% in a 400-sample simulation of this design. Below 85% of 60
  # samples has a probability of about .2% if the true coverage is 94%.
  expect_true(all(rowMeans(covered) >= .85))
})

test_that(".oblique_procrustes_batch isolates an unalignable replicate", {
  # A single non-finite slice must not abort the whole batch: it is reported
  # valid = FALSE with NA loadings/Phi/diagnostics, while the OTHER (distinct)
  # slices align independently and correctly. This mirrors the per-replicate
  # failure isolation of the per-replicate efa_procrustes() loop the batch replaces.
  set.seed(311)
  efa <- suppressWarnings(suppressMessages(
    EFA(GRiPS_raw, n_factors = 2, method = "PAF", rotation = "promax")))
  L_rot <- efa$rot_loadings
  p <- nrow(L_rot); m <- ncol(L_rot); N <- nrow(GRiPS_raw)

  # three DISTINCT replicate loadings, so a cross-slice contamination bug cannot
  # hide behind identical input slices
  slice_fit <- function(seed) {
    set.seed(seed)
    ind <- sample(N, N, replace = TRUE)
    suppressWarnings(.estimate_model(stats::cor(GRiPS_raw[ind, ]), method = "PAF",
        n_factors = m, N = N, type = "EFAtools", lean = TRUE))$unrot_loadings
  }
  cube <- array(NA_real_, c(p, m, 3))
  cube[, , 1] <- slice_fit(11)
  cube[, , 2] <- slice_fit(12)
  cube[, , 3] <- slice_fit(13)
  cube[1, 1, 2] <- Inf                       # make the middle replicate unalignable

  set.seed(1)
  res <- .oblique_procrustes_batch(cube, L_rot, random_starts = 5)

  expect_identical(as.logical(res$valid), c(TRUE, FALSE, TRUE))
  expect_true(all(is.na(res$loadings[, , 2])))           # invalid slice -> NA, not garbage
  expect_true(all(is.na(res$Phi[, , 2])))
  # invalid slice: every diagnostic is NA, only `valid` is a definite FALSE
  expect_true(all(is.na(c(res$value[2], res$iterations[2],
                          res$convergence[2], res$line_search_failed[2]))))
  expect_true(all(is.finite(res$loadings[, , c(1, 3)])))
  expect_true(all(is.finite(res$Phi[, , c(1, 3)])))

  # the surviving slices are aligned independently: distinct inputs give distinct
  # outputs (a contamination bug would tie them together)...
  expect_gt(max(abs(res$loadings[, , 1] - res$loadings[, , 3])), 1e-3)
  # ...and the first surviving slice reproduces a stand-alone single-slice
  # alignment of the same input under the same RNG offset (uncontaminated result)
  set.seed(1)
  solo <- .oblique_procrustes_batch(cube[, , 1, drop = FALSE], L_rot, random_starts = 5)
  expect_equal(res$loadings[, , 1], solo$loadings[, , 1], tolerance = 1e-10)
})

test_that(".oblique_procrustes_batch matches the per-replicate efa_procrustes alignment", {
  skip_on_cran()
  # The batched alignment must reproduce the single-matrix oblique efa_procrustes path
  # it replaces (same warm start, same R::rnorm random-start stream) to tolerance.
  set.seed(321)
  efa <- suppressWarnings(suppressMessages(
    EFA(GRiPS_raw, n_factors = 3, method = "PAF", rotation = "promax")))
  L_rot <- efa$rot_loadings
  p <- nrow(L_rot); m <- ncol(L_rot)
  N <- nrow(GRiPS_raw)

  b <- 10L
  cube <- array(NA_real_, c(p, m, b))
  for (i in seq_len(b)) {
    ind <- sample(N, N, replace = TRUE)
    fit <- suppressWarnings(.estimate_model(stats::cor(GRiPS_raw[ind, ]),
              method = "PAF", n_factors = m, N = N, type = "EFAtools", lean = TRUE))
    cube[, , i] <- fit$unrot_loadings
  }

  set.seed(99)
  loop_L <- array(NA_real_, c(p, m, b))
  for (j in seq_len(b)) {
    loop_L[, , j] <- efa_procrustes(cube[, , j], Target = L_rot, rotation = "oblique",
                                    oblique_random_starts = 5)$loadings
  }
  set.seed(99)
  batch <- .oblique_procrustes_batch(cube, L_rot, random_starts = 5)

  expect_equal(batch$loadings, loop_L, tolerance = 1e-8)
})

test_that("oblique np-boot runs end to end for a single-factor model", {
  skip_on_cran()
  # A one-factor model with an oblique-family rotation reaches the batch's m == 1
  # closed-form path; the bootstrap must run end to end and return finite SEs.
  set.seed(404)
  res <- suppressWarnings(suppressMessages(
    EFA(GRiPS_raw, n_factors = 1, method = "PAF", rotation = "promax",
        se = "np-boot", b_boot = 8)))

  expect_s3_class(res, "EFA")
  expect_false(is.null(res$SE$rot_loadings))
  expect_true(all(is.finite(res$SE$rot_loadings)))
  expect_identical(dim(res$replicates$rot_loadings)[3], 8L)
})
