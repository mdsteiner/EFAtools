# Automatic routing of EFA_POOLED()'s standard-error pooling by the SE method of
# its component fits, and the classed conditions that replace the former silent
# `settings$se <- "none"` downgrade.
#
#   - .efa_pooled_route() maps the (identical) component se onto one of
#     "none" / "information" / "sandwich" / "np-boot", or "mixed".
#   - Mixed se aborts (efa_pooled_mixed_se).
#   - A route that cannot produce pooled SEs warns (efa_pooled_se_unavailable)
#     and falls back to point-estimate-only pooling, downgrading settings$se.
#   - The sandwich/MI2S route aborts directly on a structural failure
#     (efa_pooled_mi2s_acov_not_psd), with no se-unavailable fallback.

cormat <- test_models$baseline$cormat
p_vars <- ncol(cormat)
k      <- 3L
N_id   <- 500L

# Fail if the umbrella se-unavailable warning fires; return the expression value.
expect_no_se_unavailable <- function(expr) {
  withCallingHandlers(
    expr,
    efa_pooled_se_unavailable = function(w) {
      testthat::fail("unexpected efa_pooled_se_unavailable warning")
    }
  )
}

# Minimal EFA-shaped fit mock: enough for .efa_pooled_check_fits() and
# .efa_pooled_route() to run (loadings/orig_R dims + names, matching method /
# rotation / n_factors settings, and the se under test).
.route_mock_fit <- function(se, p = 4L, kk = 1L) {
  R <- diag(p)
  dimnames(R) <- list(paste0("v", seq_len(p)), paste0("v", seq_len(p)))
  L <- matrix(0.5, p, kk)
  list(unrot_loadings = L, orig_R = R,
       settings = list(se = se, method = "ML", rotation = "none",
                       n_factors = kk))
}

# m ordinal imputations from a one-factor model, discretised at the quartiles.
.route_ord_imps <- function(m, n = 300L, p = 6L, base_seed = 700L) {
  lapply(seq_len(m), function(d) {
    set.seed(base_seed + d)
    f <- stats::rnorm(n)
    X <- outer(f, rep(0.6, p)) + matrix(stats::rnorm(n * p, sd = 0.85), n, p)
    out <- as.data.frame(apply(X, 2L, function(col) {
      as.integer(cut(col, c(-Inf, stats::quantile(col, c(.25, .5, .75)), Inf)))
    }))
    names(out) <- paste0("x", seq_len(p))
    out
  })
}

# ---- 1. .efa_pooled_route() -------------------------------------------------

test_that(".efa_pooled_route maps identical component se and flags mixed se", {
  for (s in c("none", "information", "sandwich", "np-boot")) {
    fits <- replicate(3L, .route_mock_fit(s), simplify = FALSE)
    expect_identical(EFAtools:::.efa_pooled_route(fits), s)
  }

  mixed <- list(.route_mock_fit("information"), .route_mock_fit("sandwich"))
  expect_identical(EFAtools:::.efa_pooled_route(mixed), "mixed")

  # A NULL/absent se is treated as non-conformable (mixed), not silently coerced.
  missing_se <- list(.route_mock_fit("information"),
                     list(settings = list(method = "ML")))
  expect_identical(EFAtools:::.efa_pooled_route(missing_se), "mixed")
})

# ---- 2. Information happy path ----------------------------------------------

test_that("information route pools SEs without an se-unavailable warning", {
  imps <- replicate(5L, cormat, simplify = FALSE)

  pooled <- expect_no_se_unavailable(suppressMessages(EFA_POOLED(
    imps, n_factors = k, N = N_id, method = "ML", rotation = "none",
    se = "information"
  )))

  expect_identical(pooled$settings$se, "information")
  expect_identical(pooled$settings$component_se, "information")
  expect_false(is.null(pooled$SE$unrot_loadings))
  expect_true(all(is.finite(pooled$SE$unrot_loadings)))
})

# ---- 3. Sandwich (MI2S) happy path ------------------------------------------

test_that("sandwich route fits a single MI2S model exposed as mi_fit", {
  skip_on_cran()
  imps <- .route_ord_imps(m = 5L)

  pooled <- expect_no_se_unavailable(suppressMessages(suppressWarnings(
    EFA_POOLED(imps, n_factors = 1L, method = "DWLS", rotation = "none",
               se = "sandwich", cor_method = "poly")
  )))

  expect_identical(pooled$settings$se, "sandwich")
  expect_identical(pooled$settings$component_se, "sandwich")
  expect_s3_class(pooled$mi_fit, "EFA")
  expect_true(all(is.finite(pooled$SE$unrot_loadings)))

  # The print labels the SE source as the robust sandwich covariance.
  note <- cli::ansi_strip(format(summary(pooled)))
  note <- note[grepl("Note:", note)]
  expect_true(any(grepl("sandwich", note)))
})

# ---- 4. np-boot happy path --------------------------------------------------

test_that("np-boot route pools bootstrap replicates without downgrade", {
  skip_on_cran()
  set.seed(42)
  imps <- lapply(1:3, function(i) {
    GRiPS_raw[sample(seq_len(nrow(GRiPS_raw)), replace = TRUE), ]
  })

  pooled <- expect_no_se_unavailable(suppressMessages(suppressWarnings(
    EFA_POOLED(imps, n_factors = 1L, method = "ML", rotation = "none",
               se = "np-boot", b_boot = 6L)
  )))

  expect_identical(pooled$settings$se, "np-boot")
  expect_false(is.null(pooled$SE$unrot_loadings))
  expect_false(is.null(pooled$replicates))

  # The print labels the SE source as the bootstrap/MI intervals.
  body_np <- cli::ansi_strip(format(summary(pooled)))
  expect_true(any(grepl("Bootstrap", body_np)))
})

# ---- 5. None route ----------------------------------------------------------

test_that("none route pools point estimates with no SE machinery or warning", {
  imps <- replicate(3L, cormat, simplify = FALSE)

  pooled <- expect_no_se_unavailable(
    EFA_POOLED(imps, n_factors = k, N = N_id, method = "PAF", rotation = "none")
  )

  expect_identical(pooled$settings$se, "none")
  expect_false(any(c("SE", "CI", "MI") %in% names(pooled)))
})

# ---- 6. Mixed se aborts -----------------------------------------------------

test_that("mixed component se aborts with efa_pooled_mixed_se", {
  i <- 0L
  fake_efa <- function(x, ...) {
    i <<- i + 1L
    .route_mock_fit(if (i == 1L) "information" else "sandwich",
                    p = p_vars, kk = k)
  }
  testthat::local_mocked_bindings(EFA = fake_efa)

  expect_error(
    EFA_POOLED(list(cormat, cormat), n_factors = k, N = N_id, method = "ML",
               rotation = "none"),
    class = "efa_pooled_mixed_se"
  )
})

# ---- 7. Information route: pooling failure -> se-unavailable + fallback ------

test_that("an analytic-pool abort falls back with efa_pooled_se_unavailable", {
  imps <- replicate(5L, cormat, simplify = FALSE)

  # Simulate an unreliable analytic covariance (Heywood / singular information)
  # surfacing from the pool, regardless of the alignment branch.
  testthat::local_mocked_bindings(
    .efa_pooled_analytic_pool = function(...) {
      cli::cli_abort("unreliable", class = "efa_pooled_unreliable_vcov")
    }
  )

  expect_warning(
    pooled <- suppressMessages(EFA_POOLED(
      imps, n_factors = k, N = N_id, method = "ML", rotation = "none",
      se = "information"
    )),
    class = "efa_pooled_se_unavailable"
  )

  # The point-estimate pooling still completed; only the SEs were dropped.
  expect_false(is.null(pooled$unrot_loadings))
  expect_false(is.null(pooled$residuals))
  expect_false(is.null(pooled$fit_indices))
  expect_false("SE" %in% names(pooled))
  expect_identical(pooled$settings$se, "none")
  expect_identical(pooled$settings$component_se, "information")
})

test_that("the reliability gate aborts on an NA-filled vcov under any alignment", {
  base_fit <- EFA(cormat, n_factors = k, N = N_id, method = "ML",
                  rotation = "none", se = "information")
  fits <- replicate(3L, base_fit, simplify = FALSE)
  fits[[2]]$vcov_unrot_loadings[] <- NA_real_

  unrot_list <- lapply(fits, function(f) unclass(f$unrot_loadings))
  aligned <- EFAtools:::.efa_pooled_align_unrotated_list(
    unrot_list, align_unrotated = "signed_tucker_congruence", return_meta = TRUE
  )

  # The default alignment now aborts too (previously only "procrustes" did), so
  # the dispatcher can fall back cleanly instead of emitting all-NA pooled SEs.
  expect_error(
    EFAtools:::.efa_pooled_analytic_pool(
      fits = fits,
      unrot_loadings_aligned = aligned$loadings,
      align_meta = aligned$meta,
      ci = 0.95,
      align_unrotated = "signed_tucker_congruence"
    ),
    class = "efa_pooled_unreliable_vcov"
  )
})

test_that("the reliability gate signals a classed abort (not a bare error) for >= 2 bad imputations", {
  # The classed abort's message uses cli pluralizers; with two or more flagged
  # imputations a {?...} pluralizer without an explicit cli::qty() count would
  # throw a bare formatting error that escapes the route's tryCatch. Assert the
  # abort is the documented classed condition for both one and several bad fits.
  base_fit <- EFA(cormat, n_factors = k, N = N_id, method = "ML",
                  rotation = "none", se = "information")
  unrot_list <- lapply(replicate(4L, base_fit, simplify = FALSE),
                       function(f) unclass(f$unrot_loadings))
  aligned <- EFAtools:::.efa_pooled_align_unrotated_list(
    unrot_list, align_unrotated = "signed_tucker_congruence", return_meta = TRUE
  )

  for (bad in list(2L, c(2L, 3L), c(1L, 2L, 4L))) {
    fits <- replicate(4L, base_fit, simplify = FALSE)
    for (d in bad) fits[[d]]$vcov_unrot_loadings[] <- NA_real_
    expect_error(
      EFAtools:::.efa_pooled_analytic_pool(
        fits = fits,
        unrot_loadings_aligned = aligned$loadings,
        align_meta = aligned$meta,
        ci = 0.95,
        align_unrotated = "signed_tucker_congruence"
      ),
      class = "efa_pooled_unreliable_vcov"
    )
  }
})

# ---- 8. np-boot route: soft failure -> se-unavailable + fallback ------------

test_that("a NULL bootstrap pool falls back with efa_pooled_se_unavailable", {
  skip_on_cran()
  # Raw data (not a correlation matrix), so the component fits actually carry
  # se = "np-boot" and the np-boot route is taken; EFA() downgrades np-boot to
  # "none" for correlation-matrix input, which would otherwise route to "none".
  set.seed(7)
  imps <- lapply(1:3, function(i) {
    GRiPS_raw[sample(seq_len(nrow(GRiPS_raw)), 150L, replace = TRUE), ]
  })

  testthat::local_mocked_bindings(
    .efa_pooled_bootstrap_pool = function(...) NULL
  )

  expect_warning(
    pooled <- suppressMessages(EFA_POOLED(
      imps, n_factors = 1L, method = "ML", rotation = "none",
      se = "np-boot", b_boot = 6L
    )),
    class = "efa_pooled_se_unavailable"
  )

  expect_identical(pooled$settings$se, "none")
  expect_identical(pooled$settings$component_se, "np-boot")
  expect_false("SE" %in% names(pooled))
  expect_false(is.null(pooled$unrot_loadings))
})

# ---- 9. Sandwich route: structural failure aborts directly ------------------

test_that("a non-PD pooled ACOV aborts directly, with no se-unavailable fallback", {
  skip_on_cran()
  imps <- .route_ord_imps(m = 4L)

  testthat::local_mocked_bindings(
    .efa_pooled_mi2s_inputs = function(...) {
      cli::cli_abort("not psd", class = "efa_pooled_mi2s_acov_not_psd")
    }
  )

  # Propagates as an error (not converted to a warning + fallback).
  expect_error(
    suppressMessages(suppressWarnings(EFA_POOLED(
      imps, n_factors = 1L, method = "DWLS", rotation = "none",
      se = "sandwich", cor_method = "poly"
    ))),
    class = "efa_pooled_mi2s_acov_not_psd"
  )
})

# ---- 10. Print labels the SE source per route -------------------------------

test_that("print labels the SE source for the information and none routes", {
  info <- suppressMessages(EFA_POOLED(
    replicate(5L, cormat, simplify = FALSE), n_factors = k, N = N_id,
    method = "ML", rotation = "none", se = "information"
  ))
  note <- cli::ansi_strip(format(summary(info)))
  note <- note[grepl("Note:", note)]
  expect_true(any(grepl("expected information matrix", note)))

  none <- EFA_POOLED(replicate(3L, cormat, simplify = FALSE), n_factors = k,
                     N = N_id, method = "PAF", rotation = "none")
  body_none <- cli::ansi_strip(format(summary(none)))
  expect_false(any(grepl("information matrix|Bootstrap", body_none)))
})
