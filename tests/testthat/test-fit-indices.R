# Unit tests for the fit-index helpers in R/fit-indices.R: the variance-accounted
# table (.compute_vars), the goodness-of-fit block (.gof, including the CFI
# noncentrality truncation, the small-N Bartlett guard, and an unlocatable RMSEA
# bound), and the common-part accounted-for index (.compute_caf). The ML/ULS/PAF
# fixtures are fitted once at the top from a seeded random matrix so every .gof
# branch is exercised on the same solutions.

efa_temp <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500)
efa_pro <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                rotation = "promax")
test_that(".compute_vars works", {
  checkmate::expect_matrix(
    .compute_vars(efa_temp$unrot_loadings, efa_temp$unrot_loadings)
  )
  checkmate::expect_matrix(
    .compute_vars(efa_pro$rot_loadings, efa_pro$unrot_loadings, efa_pro$Phi)
  )
})

set.seed(42)
efa_ml <- suppressWarnings(EFA(cbind(rnorm(100), rnorm(100), rnorm(100), rnorm(100),
                                     rnorm(100), rnorm(100)), 3, N = 500,
                               method = "ML"))
efa_uls <- suppressWarnings(EFA(cbind(rnorm(100), rnorm(100), rnorm(100), rnorm(100),
                    rnorm(100), rnorm(100)), 3, method = "ULS"))
efa_paf <- suppressWarnings(EFA(cbind(rnorm(100), rnorm(100), rnorm(100), rnorm(100),
                    rnorm(100), rnorm(100)), 3, method = "PAF"))

# .gof() derives CAF from a KMO on the residual matrix and signals the classed
# efa_caf_failed warning when that solve is unusable, in which case CAF falls back to 0
# instead of the ~.5 these near-null fixtures otherwise give. Which of the two branches ran
# is therefore part of what the expectations below need to know, so record it alongside the
# value. withCallingHandlers reads the condition without unwinding, so one evaluation yields
# both; tryCatch(warning = ) aborts the call and would force .gof() to be run a second time.
# Keying the handler on the class rather than on "warning" means an unrelated warning cannot
# silently divert the CAF expectation to the wrong branch.
.gof_caf <- function(...) {
  caf_failed <- FALSE
  gof <- withCallingHandlers(
    .gof(...),
    efa_caf_failed = function(w) {
      caf_failed <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  list(gof = gof, caf_failed = caf_failed)
}

res_ml <- .gof_caf(efa_ml$unrot_loadings, efa_ml$orig_R, efa_ml$settings$N,
                   "ML", efa_ml$fit_indices$Fm)
gof_ml <- res_ml$gof
res_uls <- .gof_caf(efa_uls$unrot_loadings, efa_uls$orig_R, efa_uls$settings$N,
                    "ULS", efa_uls$fit_indices$Fm)
gof_uls <- res_uls$gof
res_paf <- .gof_caf(efa_paf$unrot_loadings, efa_paf$orig_R, efa_paf$settings$N,
                    "PAF", NA)
gof_paf <- res_paf$gof
m <- 6 # n variables
q <- 3 # n factors
test_that(".gof works", {
  expect_type(gof_ml, "list")
  expect_named(gof_ml,
               c("chi", "df", "p_chi", "CAF", "RMSR", "SRMR", "CFI", "TLI",
                 "RMSEA", "RMSEA_LB", "RMSEA_UB", "AIC", "BIC", "ECVI", "Fm",
                 "chi_null", "df_null", "p_null"))
  expect_true(is.na(gof_ml$p_chi))    # df == 0: the chi-square test is undefined
  expect_equal(gof_ml$CFI, 1)
  expect_equal(gof_ml$RMSEA, 0)
  if (res_ml$caf_failed) {
    expect_equal(gof_ml$CAF, 0, tolerance = .01)
  } else {
    expect_equal(gof_ml$CAF, .5, tolerance = .1)
  }
  expect_equal(gof_ml$df, ((m - q)**2 - (m + q)) / 2)

  expect_type(gof_uls, "list")
  expect_named(gof_uls,
               c("chi", "df", "p_chi", "CAF", "RMSR", "SRMR", "CFI", "TLI",
                 "RMSEA", "RMSEA_LB", "RMSEA_UB", "AIC", "BIC", "ECVI", "Fm",
                 "chi_null", "df_null", "p_null"))
  expect_true(is.na(gof_uls$p_chi))   # df == 0: the chi-square test is undefined
  expect_equal(gof_uls$CFI, 1)
  expect_equal(gof_uls$RMSEA, 0)
  if (res_uls$caf_failed) {
    expect_equal(gof_uls$CAF, 0, tolerance = .01)
  } else {
    expect_equal(gof_uls$CAF, .5, tolerance = .1)
  }

  expect_equal(gof_uls$df, ((m - q)**2 - (m + q)) / 2)

  expect_type(gof_paf, "list")
  expect_named(gof_paf,
               c("chi", "df", "p_chi", "CAF", "RMSR", "SRMR", "CFI", "TLI",
                 "RMSEA", "RMSEA_LB", "RMSEA_UB", "AIC", "BIC", "ECVI", "Fm",
                 "chi_null", "df_null", "p_null"))
  expect_equal(gof_paf$chi, NA_real_)
  expect_equal(gof_paf$p_chi, NA_real_)
  expect_equal(gof_paf$CFI, NA_real_)
  expect_equal(gof_paf$RMSEA, NA_real_)
  if (res_paf$caf_failed) {
    expect_equal(gof_paf$CAF, 0, tolerance = .01)
  } else {
    expect_equal(gof_paf$CAF, .5, tolerance = .1)
  }
  expect_equal(gof_paf$df, ((m - q)**2 - (m + q)) / 2)
  expect_equal(gof_paf$chi_null, NA_real_)
  expect_equal(gof_paf$df_null, NA_real_)
  expect_equal(gof_paf$p_null, NA_real_)
})


test_that(".gof CFI uses the Bentler noncentrality truncation (matches lavaan)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  # Random near-uncorrelated data: a one-factor model barely improves on the
  # baseline, so the model noncentrality (chi - df) is negative. The Bentler
  # convention floors it at 0 (CFI = 1); the previous unfloored ratio deflated
  # CFI toward 0 (here it returned about 0.21).
  set.seed(1)
  X <- matrix(stats::rnorm(100 * 7), 100, 7)
  R <- stats::cor(X)
  colnames(R) <- rownames(R) <- paste0("V", seq_len(7))

  cfi_efa <- suppressWarnings(
    EFA(R, n_factors = 1, N = 100, method = "ML")
  )$fit_indices$CFI

  mod <- paste0("f =~ ", paste(colnames(R), collapse = " + "))
  fit_l <- lavaan::cfa(mod, sample.cov = R, sample.nobs = 100, std.lv = TRUE)
  cfi_lav <- unname(lavaan::fitMeasures(fit_l, "cfi"))

  expect_equal(cfi_efa, cfi_lav, tolerance = 1e-3)
  expect_equal(cfi_efa, 1)
})


test_that(".gof guards the Bartlett multiplier against tiny N (no negative chi)", {
  # For small N relative to the number of variables the model multiplier
  # N - 1 - (2m + 5)/6 - (2q)/3 turns non-positive (here N = 5, m = 6, q = 3).
  # The chi-square must then be NA, not a negative statistic masquerading as
  # perfect fit (p_chi -> 1, RMSEA floored to 0, CFI from a meaningless number).
  gof_tiny <- .gof(efa_ml$unrot_loadings, efa_ml$orig_R, N = 5, "ML",
                   efa_ml$fit_indices$Fm)
  expect_true(is.na(gof_tiny$chi))
  expect_true(is.na(gof_tiny$p_chi))
  expect_true(is.na(gof_tiny$CFI))
  expect_true(is.na(gof_tiny$RMSEA))

  # The null (baseline) multiplier N - 1 - (2p + 5)/6 is guarded the same way.
  expect_true(is.na(.null_chisq(efa_ml$orig_R, N = 3)))

  # A missing N must propagate NA, not crash the `if (mult <= 0)` guard with
  # if(NA), and an undefined chi-square has no RMSEA noncentrality bound.
  expect_true(is.na(.null_chisq(efa_ml$orig_R, N = NA_real_)))
  expect_true(is.na(.rmsea_lambda(NA_real_, df = 10, goal = .95)))
  # An undefined df has no bound either, and must not reach the `if (pchisq(NA) >= goal)`
  # comparison.
  expect_true(is.na(.rmsea_lambda(20, df = NA_real_, goal = .95)))
})


test_that("an unlocatable RMSEA noncentrality bound is NA, not an error", {
  # A statistic whose tail probability never crosses the target quantile leaves the
  # noncentrality root unbracketed even after extendInt: the bound is undefined, so it is
  # reported as NA rather than aborting the fit. chi = Inf makes goal - pchisq(Inf, df, ncp)
  # negative for every ncp, which is exactly that case.
  expect_true(is.na(.rmsea_lambda(Inf, df = 10, goal = .95)))

  # An undefined common-scale statistic must survive the caps in .chi_fit_indices() as NA
  # rather than failing an `if (NA > 1)` comparison, for the point estimate and both
  # bounds. The reported chi-square stays finite here, so only the RMSEA block is undefined.
  idx <- .chi_fit_indices(chi = 30, df = 10, chi_null = 300, df_null = 45, N = 200,
                          m = 10, ci = TRUE, chi_cfi = NA_real_, chi_null_cfi = 300)
  expect_true(is.na(idx$RMSEA))
  expect_true(is.na(idx$RMSEA_LB))
  expect_true(is.na(idx$RMSEA_UB))
  expect_true(is.finite(idx$p_chi))
})


test_that(".compute_caf returns 0 (with warning) when KMO is not computable", {
  # Hollow residual matrix: solve() succeeds (not a try-error) but the inverse
  # has a negative diagonal, so KMO is NaN. CAF must fall back to 0, not NaN.
  delta_hat <- matrix(c(0, .5, .5, .5, 0, .5, .5, .5, 0), 3)
  expect_warning(caf <- .compute_caf(delta_hat), class = "efa_caf_failed")
  expect_equal(caf, 0)
})


rm(efa_pro, efa_temp, efa_ml, efa_uls, efa_paf, gof_ml, gof_uls, gof_paf,
   res_ml, res_uls, res_paf, .gof_caf, m, q)
