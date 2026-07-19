## Tests for the exported efa_reliability() front-end: adapter dispatch ->
## .reliability_core -> tidy long result, the coefficient selector, and errors.

efa_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
sl_mod <- efa_schmid_leiman(efa_mod, estimate_control = estimate_control(type = "EFAtools"), estimator = "PAF")
fc <- sl_mod$sl[, c("F1", "F2", "F3")] >= .2

test_that("efa_reliability wires the SL adapter, core, and result builder together", {
  spec <- .rel_adapt_SL(sl_mod, factor_corres = fc, type = "EFAtools")
  core <- .reliability_core(spec, "correlation", add_ind = TRUE, add_rel = TRUE)
  expected <- .reliability_result(core, settings = list(variance = "correlation"))

  res <- efa_reliability(sl_mod, factor_map = fc)
  expect_s3_class(res, "efa_reliability")
  expect_equal(res, expected)
})

test_that("efa_reliability surfaces the full coefficient menu; CR/AVE never appear", {
  res <- efa_reliability(sl_mod, factor_map = fc)
  expect_setequal(res$coefficient,
                  c("omega_total", "omega_hierarchical", "omega_subscale",
                    "alpha", "H", "ECV", "PUC"))
  expect_false(any(c("CR", "AVE", "GLB", "beta") %in% res$coefficient))
})

test_that("the numbers tie to the validated OMEGA path", {
  res <- efa_reliability(sl_mod, factor_map = fc)
  om <- OMEGA(sl_mod, type = "EFAtools", factor_corres = fc)

  get <- function(coef, fac) res$value[res$coefficient == coef & res$factor == fac]
  expect_equal(get("omega_total", "g"), unname(om["g", "tot"]))
  expect_equal(get("omega_hierarchical", "F2"), unname(om["F2", "hier"]))
  expect_equal(get("H", "F3"), unname(om["F3", "H"]))
  expect_equal(get("ECV", "g"), unname(om["g", "ECV"]))
})

test_that("the coefficients selector returns only the requested coefficients", {
  res <- efa_reliability(sl_mod, factor_map = fc,
                         coefficients = c("omega_total", "alpha"))
  expect_s3_class(res, "efa_reliability")
  expect_setequal(res$coefficient, c("omega_total", "alpha"))
  # The kind attribute is subset to match, and the values are unchanged.
  expect_setequal(names(attr(res, "kind")), c("omega_total", "alpha"))
  full <- efa_reliability(sl_mod, factor_map = fc)
  expect_equal(res$value[res$coefficient == "alpha"],
               full$value[full$coefficient == "alpha"])
})

test_that("selecting only common-variance indices drops the reliability rows", {
  res <- efa_reliability(sl_mod, factor_map = fc,
                         coefficients = c("ECV", "PUC"))
  expect_setequal(res$coefficient, c("ECV", "PUC"))
  expect_true(all(res$level == "general"))
})

test_that("an unknown coefficient is a classed error", {
  expect_error(
    efa_reliability(sl_mod, factor_map = fc, coefficients = "omega"),
    class = "efa_reliability_bad_coefficients")
})

test_that("invalid model and missing/malformed manual components are classed errors", {
  expect_error(efa_reliability(model = 1:5), class = "efa_reliability_invalid_model")
  expect_error(
    efa_reliability(model = NULL, var_names = rownames(sl_mod$sl),
                    s_load = sl_mod$sl[, c("F1", "F2", "F3")], u2 = sl_mod$sl[, "u2"],
                    factor_map = fc),
    class = "efa_reliability_missing_args")
  # A vector s_load must abort, not be silently coerced to a one-column matrix.
  expect_error(
    efa_reliability(model = NULL, var_names = rownames(sl_mod$sl),
                    g_load = sl_mod$sl[, "g"], s_load = as.numeric(sl_mod$sl[, "F1"]),
                    u2 = sl_mod$sl[, "u2"]),
    class = "efa_reliability_bad_s_load")
  # A non-numeric matrix model is rejected rather than failing deep in the algebra.
  expect_error(efa_reliability(matrix("a", 4, 2)), class = "efa_reliability_bad_matrix")
})

test_that("correlation mode without an available correlation matrix aborts (not silent NA)", {
  s <- matrix(0, 6, 2); s[1:3, 1] <- 0.4; s[4:6, 2] <- 0.4
  # Manual components with no cormat (and no pattern/Phi): correlation mode has no
  # denominator, so it aborts with an actionable message instead of returning NA omegas.
  expect_error(
    efa_reliability(model = NULL, var_names = paste0("V", 1:6), g_load = rep(0.5, 6),
                    s_load = s, u2 = rep(0.5, 6)),
    class = "efa_reliability_need_cormat")
  # The same input is fine under variance = "sums_load" (no correlation matrix needed).
  expect_s3_class(
    efa_reliability(model = NULL, var_names = paste0("V", 1:6), g_load = rep(0.5, 6),
                    s_load = s, u2 = rep(0.5, 6), variance = "sums_load"),
    "efa_reliability")
})

test_that("manual components of mismatched length abort with a classed error", {
  s <- matrix(0, 5, 2); s[1:3, 1] <- 0.4; s[4:5, 2] <- 0.4
  expect_error(
    efa_reliability(model = NULL, var_names = paste0("V", 1:6), g_load = rep(0.5, 6),
                    s_load = s, u2 = rep(0.5, 6)),
    class = "efa_reliability_length_mismatch")
})

test_that("a correlated-factors manual input needs a correlation matrix (sums_load too)", {
  s <- matrix(0, 6, 2); s[1:3, 1] <- 0.4; s[4:6, 2] <- 0.4
  # g_load all zero = a manual no-general (correlated-factors) input. sums_load would ignore
  # the factor correlations, so a correlation matrix is required either way -- otherwise the
  # whole-scale omega total would be silently wrong while the result is treated as no-general.
  expect_error(
    efa_reliability(model = NULL, var_names = paste0("V", 1:6), g_load = rep(0, 6),
                    s_load = s, u2 = 1 - rowSums(s^2), variance = "sums_load"),
    class = "efa_reliability_need_cormat")
  expect_error(
    efa_reliability(model = NULL, var_names = paste0("V", 1:6), g_load = rep(0, 6),
                    s_load = s, u2 = 1 - rowSums(s^2)),
    class = "efa_reliability_need_cormat")
})

test_that("omitting factor_map auto-maps each item to its highest-loading factor", {
  # No factor_map: each item is assigned to its highest |group loading|, so the
  # result matches supplying that same 0/1 assignment explicitly.
  s <- sl_mod$sl[, c("F1", "F2", "F3")]
  auto <- matrix(0, nrow(s), ncol(s))
  for (i in seq_len(nrow(s))) auto[i, which.max(abs(s[i, ]))] <- 1
  expect_equal(efa_reliability(sl_mod), efa_reliability(sl_mod, factor_map = auto))
})

test_that("an oblique EFA is scored as a correlated-factors model (no general factor)", {
  res <- efa_reliability(efa_mod)
  # No general factor: the bifactor indices omega hierarchical, ECV, and PUC are not
  # meaningful and are omitted.
  expect_false(any(res$coefficient == "omega_hierarchical"))
  expect_false(any(res$coefficient == "ECV"))
  expect_false(any(res$coefficient == "PUC"))
  # The synthetic all-zero general-factor column also leaves the whole-scale (g) row's
  # omega subscale and H index undefined; they are dropped (H of an all-zero loading
  # vector would otherwise print a misleading 0). Group factors keep their own omega
  # subscale and H.
  expect_false(any(res$coefficient == "omega_subscale" & res$factor == "g"))
  expect_false(any(res$coefficient == "H" & res$factor == "g"))
  expect_true(any(res$coefficient == "omega_subscale" & res$factor != "g"))
  expect_true(any(res$coefficient == "H" & res$factor != "g"))
  # The whole-scale omega total and alpha, and each group factor's omega/H, remain.
  expect_true(any(res$coefficient == "omega_total" & res$factor == "g"))
  expect_true(any(res$coefficient == "alpha" & res$factor == "g"))
  spec <- .rel_adapt_efa(efa_mod, type = "psych")
  expect_equal(res$value[res$coefficient == "omega_total" & res$factor == "g"],
               1 - sum(spec$u2) / sum(spec$cormat))
})

test_that("sums_load on an oblique EFA warns and falls back to correlation", {
  # The correlation-based total is correct for correlated factors; sums_load would ignore
  # Phi. The default is correlation, so a plain call does not warn.
  expect_no_warning(efa_reliability(efa_mod))
  expect_warning(res <- efa_reliability(efa_mod, variance = "sums_load"),
                 class = "efa_reliability_variance_override")
  expect_identical(attr(res, "settings"), list(variance = "correlation"))
  # The fallback yields the same numbers as an explicit correlation request.
  expect_equal(res$value, suppressWarnings(efa_reliability(efa_mod))$value)
})

test_that("an explicitly requested but undefined coefficient is reported, not silent", {
  # omega hierarchical is not defined for a correlated-factors EFA.
  expect_message(efa_reliability(efa_mod, coefficients = "omega_hierarchical"),
                 class = "efa_reliability_coef_undefined")
  res <- suppressMessages(
    efa_reliability(efa_mod, coefficients = c("omega_total", "omega_hierarchical")))
  # The defined coefficient is still returned; the undefined one is dropped.
  expect_setequal(res$coefficient, "omega_total")
})

test_that("a raw bifactor loading matrix is accepted", {
  L <- cbind(g = rep(0.5, 6), matrix(c(rep(0.4, 3), rep(0, 6), rep(0.4, 3)), 6, 2))
  colnames(L) <- c("g", "F1", "F2")
  rownames(L) <- paste0("V", 1:6)
  res <- efa_reliability(L)
  expect_s3_class(res, "efa_reliability")
  expect_setequal(unique(res$factor), c("g", "F1", "F2"))
})

## lavaan paths --------------------------------------------------------------

test_that("a lavaan bifactor fit matches OMEGA and yields one block", {
  skip_if_not_installed("lavaan")
  mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
          F2 =~ V7 + V8 + V9 + V10 + V11 + V12
          F3 =~ V13 + V14 + V15 + V16 + V17 + V18
          g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
               V13 + V14 + V15 + V16 + V17 + V18'
  fit <- suppressWarnings(lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
                                      sample.nobs = 500, estimator = "ml",
                                      orthogonal = TRUE))
  res <- efa_reliability(fit, g_name = "g")
  om <- OMEGA(fit, g_name = "g")

  expect_true(all(is.na(res$group)))       # single ungrouped fit
  expect_identical(attr(res, "settings"), list(variance = "sums_load"))
  get <- function(coef, fac) res$value[res$coefficient == coef & res$factor == fac]
  expect_equal(get("omega_total", "g"), unname(om["g", "tot"]))
  expect_equal(get("H", "F1"), unname(om["F1", "H"]))
})

test_that("a lavaan single-factor fit surfaces only omega total and H, and informs", {
  skip_if_not_installed("lavaan")
  fit <- suppressWarnings(lavaan::cfa(
    'g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
          V13 + V14 + V15 + V16 + V17 + V18',
    sample.cov = test_models$baseline$cormat, sample.nobs = 500,
    estimator = "ml", orthogonal = TRUE))

  expect_message(efa_reliability(fit, g_name = "g"),
                 class = "efa_reliability_single_factor")
  res <- suppressMessages(efa_reliability(fit, g_name = "g"))
  expect_setequal(res$coefficient, c("omega_total", "H"))
  expect_true(all(res$factor == "g"))
})

test_that("a lavaan second-order fit informs about the Schmid-Leiman transform", {
  skip_if_not_installed("lavaan")
  fit <- suppressWarnings(lavaan::cfa(
    'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
     F2 =~ V7 + V8 + V9 + V10 + V11 + V12
     F3 =~ V13 + V14 + V15 + V16 + V17 + V18
     g =~ F1 + F2 + F3',
    sample.cov = test_models$baseline$cormat, sample.nobs = 500, estimator = "ml"))
  expect_message(efa_reliability(fit, g_name = "g"),
                 class = "efa_reliability_g_second_order")
})

test_that("a lavaan multiple-group fit yields one labelled block per group", {
  skip_if_not_installed("lavaan")
  mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
          F2 =~ V7 + V8 + V9 + V10 + V11 + V12
          F3 =~ V13 + V14 + V15 + V16 + V17 + V18
          g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
               V13 + V14 + V15 + V16 + V17 + V18'
  fit <- suppressWarnings(lavaan::cfa(mod,
    sample.cov = list(test_models$baseline$cormat, test_models$baseline$cormat),
    sample.nobs = c(500, 500), estimator = "ml", orthogonal = TRUE))
  res <- efa_reliability(fit, g_name = "g", group_names = c("Young", "Old"))
  expect_identical(unique(res$group), c("Young", "Old"))

  # Duplicate labels would merge groups in the tidy output, so they abort.
  expect_error(efa_reliability(fit, g_name = "g", group_names = c("Dup", "Dup")),
               class = "efa_reliability_group_names")
})

test_that("a bifactor with an under-loaded item informs about too few loadings", {
  skip_if_not_installed("lavaan")
  # g loads on all items except V11/V12, so those load on a single factor.
  mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
          F2 =~ V7 + V8 + V9 + V10 + V11 + V12
          F3 =~ V13 + V14 + V15 + V16 + V17 + V18
          g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
               V13 + V14 + V15 + V16 + V17 + V18'
  fit <- suppressWarnings(lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
                                      sample.nobs = 500, estimator = "ml",
                                      orthogonal = TRUE))
  expect_message(efa_reliability(fit, g_name = "g"),
                 class = "efa_reliability_few_loadings")
})

test_that("print output is stable", {
  local_reproducible_output()

  # full menu, single group (reliability and common-variance blocks)
  expect_snapshot(print(efa_reliability(sl_mod, factor_map = fc)),
                  transform = scrub_num)

  # coefficient selector (only the requested columns)
  expect_snapshot(print(efa_reliability(sl_mod, factor_map = fc,
                                        coefficients = c("omega_total", "alpha"))),
                  transform = scrub_num)

  # oblique EFA (no general factor): the g row omits omega subscale and H, which print
  # blank rather than as "NA".
  expect_snapshot(print(efa_reliability(efa_mod)), transform = scrub_num)
})

rm(efa_mod, sl_mod, fc)
