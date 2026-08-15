## Tests for the exported efa_reliability() front-end: adapter dispatch ->
## .reliability_core -> tidy long result, the coefficient selector, and errors.

efa_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
sl_mod <- efa_schmid_leiman(efa_mod, estimate_control = estimate_control(type = "EFAtools"), estimator = "PAF")
fc <- sl_mod$sl[, c("F1", "F2", "F3")] >= .2

# The map efa_reliability derives when factor_map is omitted: each item on its highest
# |group loading|.
auto <- local({
  s <- sl_mod$sl[, c("F1", "F2", "F3")]
  m <- matrix(0, nrow(s), ncol(s))
  for (i in seq_len(nrow(s))) m[i, which.max(abs(s[i, ]))] <- 1
  m
})

test_that("efa_reliability wires the SL adapter, core, and result builder together", {
  spec <- .rel_adapt_SL(sl_mod, factor_corres = fc, type = "EFAtools")
  core <- .reliability_core(spec, "correlation", add_ind = TRUE, add_rel = TRUE)
  expected <- .reliability_result(core, settings = list(variance = "correlation",
                                                        no_general = FALSE))

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
  g <- rep(0.5, 6)
  u2 <- 1 - g^2 - rowSums(s^2)
  # Manual components with no cormat (and no pattern/Phi): correlation mode has no
  # denominator, so it aborts with an actionable message instead of returning NA omegas.
  expect_error(
    efa_reliability(model = NULL, var_names = paste0("V", 1:6), g_load = g,
                    s_load = s, u2 = u2),
    class = "efa_reliability_need_cormat")
  # The same input is fine under variance = "sums_load" (no correlation matrix needed).
  expect_s3_class(
    efa_reliability(model = NULL, var_names = paste0("V", 1:6), g_load = g,
                    s_load = s, u2 = u2, variance = "sums_load"),
    "efa_reliability")
})

test_that("manual components of mismatched length abort with a classed error", {
  s <- matrix(0, 5, 2); s[1:3, 1] <- 0.4; s[4:5, 2] <- 0.4
  expect_error(
    efa_reliability(model = NULL, var_names = paste0("V", 1:6), g_load = rep(0.5, 6),
                    s_load = s, u2 = rep(0.5, 6)),
    class = "efa_reliability_length_mismatch")
})

test_that("only correlation mode needs a correlation matrix, general factor or not", {
  s <- matrix(0, 6, 2); s[1:3, 1] <- 0.4; s[4:6, 2] <- 0.4
  u2 <- 1 - rowSums(s^2)
  args <- list(model = NULL, var_names = paste0("V", 1:6), s_load = s, u2 = u2)

  # Correlation mode has no denominator without one, with or without a general factor.
  # The uniquenesses follow each general column so the components stay standardized.
  for (g in list(rep(0, 6), rep(0.3, 6))) {
    expect_error(
      do.call(efa_reliability,
              c(args[setdiff(names(args), "u2")],
                list(g_load = g, u2 = 1 - g^2 - rowSums(s^2)))),
      class = "efa_reliability_need_cormat")
  }

  # The model-implied convention needs none: its denominator is the composite's common
  # variance plus the uniquenesses. An all-zero general column says the solution has no
  # general factor, not that a correlation matrix is required.
  res <- efa_reliability(model = NULL, var_names = paste0("V", 1:6), g_load = rep(0, 6),
                         s_load = s, u2 = u2, variance = "sums_load")
  tot <- res$value[res$coefficient == "omega_total" & res$factor == "g"]
  expect_equal(tot, sum(tcrossprod(s)) / (sum(tcrossprod(s)) + sum(u2)))
  expect_identical(attr(res, "settings"),
                   list(variance = "sums_load", no_general = TRUE))
})

test_that("Phi without pattern makes manual components a correlated-factors solution", {
  # The defect this closes: the same solution scored through the components and through the
  # fitted object must return the same coefficients. Without Phi in the spec the numerator is
  # 1'LL'1, which describes an orthogonal model and falls short of 1'L Phi L'1.
  L <- matrix(0, 6, 2)
  L[1:3, 1] <- c(.75, .70, .65)
  L[4:6, 2] <- c(.60, .65, .70)
  rownames(L) <- paste0("V", 1:6)
  Phi <- matrix(c(1, .5, .5, 1), 2, 2)
  common <- L %*% Phi %*% t(L)
  u2 <- 1 - diag(common)
  R <- common + diag(u2)
  dimnames(R) <- list(rownames(L), rownames(L))
  manual <- function(...) efa_reliability(NULL, var_names = rownames(L),
                                          g_load = rep(0, 6), s_load = L, u2 = u2, ...)
  get <- function(x) x$value[x$coefficient == "omega_total" & x$factor == "g"]

  expect_equal(get(manual(cormat = R, Phi = Phi)), sum(common) / sum(R))
  # Without it, the orthogonal reading -- materially different, so the test is not vacuous.
  expect_equal(get(manual(cormat = R)), sum(tcrossprod(L)) / sum(R))
  expect_gt(get(manual(cormat = R, Phi = Phi)), get(manual(cormat = R)))

  # It needs no correlation matrix in the model-implied convention.
  expect_equal(get(manual(Phi = Phi, variance = "sums_load")),
               sum(common) / (sum(common) + sum(u2)))

  # A pattern says Phi belongs to that separate solution, which is stated rather than assumed.
  expect_warning(res <- manual(Phi = Phi, pattern = L),
                 class = "efa_reliability_phi_pattern")
  expect_equal(get(res), get(manual(cormat = R)))
})

test_that("missing and infinite components are rejected, not scored", {
  # Every check that might name them is blind to exactly these values -- the Heywood report
  # tests `u2 < eps` and skips NA, and the keying and standardization checks compare totals
  # that are themselves non-finite -- while the table that comes back looks ordinary: an
  # infinite uniqueness drives its composite's omega total to a flat 0, an infinite general
  # loading leaves a plausible coefficient, and a missing one silently drops rows. The
  # assertions come from checkmate, so unlike the other component checks they carry no
  # package class and can only be asserted as errors.
  L <- matrix(c(.7, .7, .7, 0, 0, 0, 0, 0, 0, .6, .6, .6), 6, 2)
  u2 <- 1 - rowSums(L^2)
  manual <- function(g = rep(0, 6), uniq = u2) {
    efa_reliability(NULL, var_names = paste0("V", 1:6), g_load = g, s_load = L,
                    u2 = uniq, variance = "sums_load")
  }

  s_bad <- function(bad) {
    efa_reliability(NULL, var_names = paste0("V", 1:6), g_load = rep(0, 6),
                    s_load = replace(L, 2, bad), u2 = u2, variance = "sums_load")
  }
  for (bad in list(NA_real_, Inf, -Inf)) {
    expect_error(manual(uniq = replace(u2, 2, bad)))
    expect_error(manual(g = replace(rep(0.3, 6), 2, bad)))
    expect_error(s_bad(bad))
  }
  expect_error(manual(uniq = rep(Inf, 6)))
  expect_error(manual(g = rep(NA_real_, 6)))

  # Checking the group loadings numerically also rejects a matrix that is not numeric,
  # which the class test alone admits.
  expect_error(efa_reliability(NULL, var_names = paste0("V", 1:6), g_load = rep(0, 6),
                               s_load = matrix(as.character(L), 6, 2), u2 = u2,
                               variance = "sums_load"))

  # The same components, complete and finite, are scored.
  expect_s3_class(manual(), "efa_reliability")
  expect_s3_class(manual(g = rep(0.3, 6), uniq = 1 - 0.09 - rowSums(L^2)),
                  "efa_reliability")

  # A bare bifactor loading matrix takes `u2` as an override of the uniquenesses it would
  # otherwise derive, so that route needs the same guard on its own account.
  Lb <- cbind(g = rep(.5, 6), L)
  rownames(Lb) <- paste0("V", 1:6)
  u2b <- 1 - rowSums(Lb^2)
  for (bad in list(NA_real_, Inf, -Inf)) {
    expect_error(efa_reliability(Lb, u2 = replace(u2b, 2, bad)))
  }
  # Omitting it, so the adapter derives them, is untouched.
  expect_s3_class(efa_reliability(Lb, u2 = u2b), "efa_reliability")
  expect_s3_class(efa_reliability(Lb), "efa_reliability")
})

test_that("supplied uniquenesses are checked against the loadings they accompany", {
  # The printed Schmid-Leiman table ends in h2 and u2, and published tables often report
  # only h2, so passing communalities where uniquenesses belong is a one-column slip. It is
  # silent otherwise: u2 does not enter the coefficients under variance = "correlation"
  # with a correlation matrix, only in the model-implied convention a user without one is
  # sent to.
  sl <- unclass(sl_mod$sl)
  g <- sl[, "g"]
  s <- sl[, c("F1", "F2", "F3")]
  h2 <- sl[, "h2"]
  manual <- function(u2, ...) {
    efa_reliability(NULL, var_names = rownames(sl), g_load = g, s_load = s, u2 = u2,
                    variance = "sums_load", ...)
  }

  expect_warning(manual(u2 = h2), class = "efa_reliability_u2_not_standardized")
  expect_no_warning(manual(u2 = sl[, "u2"]),
                    class = "efa_reliability_u2_not_standardized")

  # A published table rounded to two decimals is a legitimate transcription and must not
  # trip it, so the tolerance has to clear ordinary rounding by a wide margin.
  expect_no_warning(manual(u2 = round(sl[, "u2"], 2)),
                    class = "efa_reliability_u2_not_standardized")

  # The bare bifactor matrix takes its own uniquenesses, and derives them when it is not
  # given any -- derived ones complete by construction and must stay silent.
  L <- cbind(g = g, s)
  expect_warning(efa_reliability(L, u2 = h2, cormat = test_models$baseline$cormat),
                 class = "efa_reliability_u2_not_standardized")
  expect_no_warning(efa_reliability(L, cormat = test_models$baseline$cormat),
                    class = "efa_reliability_u2_not_standardized")

  # It reads the factor correlations where the solution has them. That only shows on a
  # solution with cross-loadings: under exact simple structure the off-diagonal terms of
  # L Phi L' never reach the diagonal, so the Phi-aware and Phi-blind uniquenesses coincide
  # and no fixture built that way could tell the two apart.
  Lo <- matrix(0, 8, 2)
  Lo[1:4, 1] <- .70; Lo[1:4, 2] <- .30
  Lo[5:8, 2] <- .70; Lo[5:8, 1] <- .30
  Phi <- matrix(c(1, .5, .5, 1), 2, 2)
  aware <- 1 - diag(Lo %*% Phi %*% t(Lo))
  blind <- 1 - rowSums(Lo^2)
  expect_gt(max(abs(aware - blind)), .05)

  obl <- function(u2) efa_reliability(NULL, var_names = paste0("V", 1:8),
                                      g_load = rep(0, 8), s_load = Lo, u2 = u2,
                                      Phi = Phi, variance = "sums_load")
  expect_no_warning(obl(aware), class = "efa_reliability_u2_not_standardized")
  expect_warning(obl(blind), class = "efa_reliability_u2_not_standardized")
})

test_that("a correlation matrix and a pattern to rebuild one are an either/or", {
  # Two ways to give one matrix. Taking both is ambiguous rather than redundant: the
  # reconstruction never runs when a cormat is in hand, so `pattern` is read by nothing --
  # yet its presence alone would decide whether `Phi` describes this solution or another.
  L <- matrix(0, 6, 2); L[1:3, 1] <- .7; L[4:6, 2] <- .6
  rownames(L) <- paste0("V", 1:6)
  Phi <- matrix(c(1, .5, .5, 1), 2, 2)
  u2 <- 1 - diag(L %*% Phi %*% t(L))
  R <- L %*% Phi %*% t(L) + diag(u2)
  dimnames(R) <- list(rownames(L), rownames(L))
  manual <- function(...) efa_reliability(NULL, var_names = rownames(L),
                                          g_load = rep(0, 6), s_load = L, u2 = u2, ...)

  expect_error(manual(cormat = R, pattern = L, Phi = Phi),
               class = "efa_reliability_cormat_and_pattern")
  # A pattern with no Phi cannot reconstruct anything either, so it is the same clash.
  expect_error(manual(cormat = R, pattern = L),
               class = "efa_reliability_cormat_and_pattern")

  # Each on its own is fine, and Phi with a cormat is the correlated-factors call.
  expect_no_error(manual(cormat = R, Phi = Phi))
  expect_no_error(suppressWarnings(manual(pattern = L, Phi = Phi)))
})

test_that("Phi is rejected or reported where it cannot describe the group factors", {
  L <- matrix(c(.7, .7, .7, 0, 0, 0, 0, 0, 0, .6, .6, .6), 6, 2)
  u2 <- 1 - rowSums(L^2)
  manual <- function(...) efa_reliability(NULL, var_names = paste0("V", 1:6),
                                          s_load = L, u2 = u2, variance = "sums_load", ...)
  ok <- matrix(c(1, .5, .5, 1), 2, 2)

  # Correlated group factors and a general factor are not the decomposition these
  # coefficients report.
  expect_error(manual(g_load = rep(.5, 6), Phi = ok),
               class = "efa_reliability_phi_with_general")

  # A Phi that is not a correlation matrix of the solution's group factors.
  for (bad in list(diag(3), diag(2, 2), matrix(c(1, .9, .1, 1), 2, 2),
                   matrix(c(1, NA, NA, 1), 2, 2), matrix(c(1, 1.4, 1.4, 1), 2, 2),
                   as.data.frame(ok), "nonsense")) {
    expect_error(manual(g_load = rep(0, 6), Phi = bad),
                 class = "efa_reliability_bad_phi")
  }

  # A Phi paired with a pattern builds the correlation matrix every coefficient is divided
  # by, so it is checked on that reading too, against the pattern's own factor count.
  expect_error(manual(g_load = rep(0, 6), pattern = L, Phi = diag(3)),
               class = "efa_reliability_bad_phi")
  expect_warning(manual(g_load = rep(0, 6), pattern = L, Phi = ok),
                 class = "efa_reliability_phi_pattern")

  # A single group factor takes a 1 x 1 Phi, which must stay reachable. It is one factor,
  # so the single-factor note is expected here and is muffled rather than left to print.
  expect_no_error(suppressMessages(
    efa_reliability(NULL, var_names = paste0("V", 1:6), g_load = rep(0, 6),
                    s_load = matrix(.6, 6, 1), u2 = rep(.64, 6),
                    Phi = matrix(1, 1, 1), variance = "sums_load")))

  # A fitted solution carries its own components, so every argument that describes manually
  # supplied ones is reported rather than dropped -- not just the two that changed meaning.
  cls <- "efa_reliability_args_ignored"
  expect_warning(efa_reliability(sl_mod, factor_map = fc, Phi = ok), class = cls)
  expect_warning(efa_reliability(sl_mod, factor_map = fc, pattern = L), class = cls)
  expect_warning(efa_reliability(sl_mod, factor_map = fc, u2 = rep(.5, 18)), class = cls)
  expect_warning(efa_reliability(sl_mod, factor_map = fc, g_load = rep(.5, 18)), class = cls)
  expect_warning(efa_reliability(sl_mod, factor_map = fc, var_names = paste0("z", 1:18)),
                 class = cls)
  # and the coefficients are untouched by any of them
  expect_equal(suppressWarnings(efa_reliability(sl_mod, factor_map = fc, Phi = ok)),
               efa_reliability(sl_mod, factor_map = fc))
  expect_equal(suppressWarnings(efa_reliability(sl_mod, factor_map = fc,
                                                u2 = rep(.5, 18))),
               efa_reliability(sl_mod, factor_map = fc))

  # A bare loading matrix is the one `model` that does read `u2`, so it must not be
  # reported there -- but it must line up with the loadings.
  Lb <- cbind(g = rep(.5, 6), matrix(c(rep(.4, 3), rep(0, 6), rep(.4, 3)), 6, 2))
  rownames(Lb) <- paste0("V", 1:6)
  expect_no_warning(efa_reliability(Lb, u2 = 1 - rowSums(Lb^2)), class = cls)
  expect_error(efa_reliability(Lb, u2 = rep(.5, 4)))

  # A lavaan fit reads neither the correspondences nor a correlation matrix, so those are
  # reported there too, and the remedy must not point the user back at them.
  skip_if_not_installed("lavaan")
  fit <- suppressWarnings(lavaan::cfa(
    paste('F1 =~ V1 + V2 + V3 + V4 + V5 + V6
           F2 =~ V7 + V8 + V9 + V10 + V11 + V12
           F3 =~ V13 + V14 + V15 + V16 + V17 + V18
           g =~', paste(paste0("V", 1:18), collapse = " + ")),
    sample.cov = test_models$baseline$cormat, sample.nobs = 500, estimator = "ml",
    orthogonal = TRUE))
  expect_warning(efa_reliability(fit, g_name = "g",
                                 cormat = test_models$baseline$cormat), class = cls)
  expect_warning(efa_reliability(fit, g_name = "g", factor_map = matrix(1, 18, 3)),
                 class = cls)
  # and they never moved a coefficient
  expect_equal(suppressWarnings(
    efa_reliability(fit, g_name = "g", cormat = test_models$baseline$cormat)),
    efa_reliability(fit, g_name = "g"))
})

test_that("a labelling argument is reported for an input that does not read it", {
  # `fac_names` and `group_names` label the result rather than describe components, so each
  # is reported where it does not apply rather than dropped: `fac_names` for a lavaan fit,
  # whose factor labels come from the model syntax, and `group_names` for anything but a
  # lavaan fit, every other input being a single ungrouped solution.
  cls <- "efa_reliability_arg_not_applicable"

  expect_warning(efa_reliability(sl_mod, factor_map = fc, group_names = c("a", "b")),
                 class = cls)
  Lb <- cbind(g = rep(.5, 6), F1 = c(rep(.4, 3), rep(0, 3)),
              F2 = c(rep(0, 3), rep(.4, 3)))
  rownames(Lb) <- paste0("V", 1:6)
  expect_warning(efa_reliability(Lb, group_names = c("a", "b"), variance = "sums_load"),
                 class = cls)
  # Manually supplied components read no `group_names` either, and are not covered by the
  # report for a fitted `model`, there being none.
  expect_warning(efa_reliability(NULL, var_names = rownames(Lb), g_load = Lb[, 1],
                                 s_load = Lb[, -1], u2 = 1 - rowSums(Lb^2),
                                 group_names = c("a", "b"), variance = "sums_load"),
                 class = cls)
  # and it never moved a coefficient
  expect_equal(suppressWarnings(
    efa_reliability(sl_mod, factor_map = fc, group_names = c("a", "b"))),
    efa_reliability(sl_mod, factor_map = fc))

  # `fac_names` is read on every non-lavaan route, so it stays silent there.
  expect_no_warning(efa_reliability(sl_mod, factor_map = fc,
                                    fac_names = c("A", "B", "C")), class = cls)

  skip_if_not_installed("lavaan")
  fit <- suppressWarnings(lavaan::cfa(
    paste('F1 =~ V1 + V2 + V3 + V4 + V5 + V6
           F2 =~ V7 + V8 + V9 + V10 + V11 + V12
           F3 =~ V13 + V14 + V15 + V16 + V17 + V18
           g =~', paste(paste0("V", 1:18), collapse = " + ")),
    sample.cov = test_models$baseline$cormat, sample.nobs = 500, estimator = "ml",
    orthogonal = TRUE))
  expect_warning(efa_reliability(fit, g_name = "g", fac_names = c("A", "B", "C")),
                 class = cls)
  # Including one of the wrong length, which the length check -- belonging to the routes
  # that read it -- does not refuse here.
  expect_warning(efa_reliability(fit, g_name = "g", fac_names = c("A", "B")), class = cls)
  expect_equal(suppressWarnings(
    efa_reliability(fit, g_name = "g", fac_names = c("A", "B", "C"))),
    efa_reliability(fit, g_name = "g"))

  # A single-group lavaan fit is scored as one ungrouped solution and drops the labels
  # unread, exactly as the inputs above do, so it is reported too -- the report keys on the
  # group count, not on the class.
  expect_warning(efa_reliability(fit, g_name = "g", group_names = "one"), class = cls)
  expect_true(all(is.na(suppressWarnings(
    efa_reliability(fit, g_name = "g", group_names = "one"))$group)))

  # The mirror, and the one input that does read them: a multiple-group fit, whose blocks
  # the labels key. Without it an inverted guard would warn here and pass the suite.
  fit_gr <- suppressWarnings(lavaan::cfa(
    paste('F1 =~ V1 + V2 + V3 + V4 + V5 + V6
           F2 =~ V7 + V8 + V9 + V10 + V11 + V12
           F3 =~ V13 + V14 + V15 + V16 + V17 + V18
           g =~', paste(paste0("V", 1:18), collapse = " + ")),
    sample.cov = list(test_models$baseline$cormat, test_models$baseline$cormat),
    sample.nobs = c(500, 500), estimator = "ml", orthogonal = TRUE))
  expect_no_warning(efa_reliability(fit_gr, g_name = "g",
                                    group_names = c("Young", "Old")), class = cls)
})

test_that("omitting factor_map auto-maps each item to its highest-loading factor", {
  # No factor_map: each item is assigned to its highest |group loading|, so the
  # result matches supplying that same 0/1 assignment explicitly.
  expect_equal(efa_reliability(sl_mod), efa_reliability(sl_mod, factor_map = auto))
})

test_that("a factor_map whose columns are in the wrong order is flagged", {
  # The columns of factor_map are matched to the group factors by position, so a
  # permuted map assigns each factor the items of another one -- well-formed, but
  # meaningless. The items then hardly load on the factor they are mapped to while
  # loading substantially on another, which is what the check picks up.
  expect_warning(efa_reliability(sl_mod, factor_map = auto[, c(2, 3, 1)]),
                 class = "efa_reliability_implausible_map")
  # A sound map -- explicit or auto-derived -- stays silent.
  expect_no_warning(efa_reliability(sl_mod, factor_map = fc))
  expect_no_warning(efa_reliability(sl_mod, factor_map = auto))
  expect_no_warning(efa_reliability(sl_mod))

  # The check belongs to the efa_reliability front-end; OMEGA()'s behaviour is
  # unchanged (it derives and validates its map in its own worker).
  expect_no_warning(OMEGA(sl_mod, type = "EFAtools",
                          factor_corres = auto[, c(2, 3, 1)]),
                    class = "efa_reliability_implausible_map")
})

test_that("a bifactor loading matrix validates and checks a supplied factor_map", {
  L <- cbind(g = rep(0.5, 6),
             F1 = c(rep(0.5, 3), rep(0, 3)),
             F2 = c(rep(0, 3), rep(0.5, 3)))
  rownames(L) <- paste0("V", 1:6)
  map <- cbind(c(rep(1, 3), rep(0, 3)), c(rep(0, 3), rep(1, 3)))

  expect_no_warning(efa_reliability(L, factor_map = map))
  # Swapping the columns points each factor at the other one's items.
  expect_warning(efa_reliability(L, factor_map = map[, c(2, 1)]),
                 class = "efa_reliability_implausible_map")
  # A map that does not conform to the group loadings is an error, not coefficients
  # above 1 from a recycled map.
  expect_error(efa_reliability(L, factor_map = map[1:3, , drop = FALSE]),
               class = "efa_reliability_map_dim")
  expect_error(efa_reliability(L, factor_map = cbind(map, 0)),
               class = "efa_reliability_map_dim")
  # a map of the right shape that is not a matrix is rejected by the same check
  expect_error(efa_reliability(L, factor_map = as.data.frame(map)),
               class = "efa_reliability_map_dim")
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
               sum(spec$s_load %*% spec$Phi %*% t(spec$s_load)) / sum(spec$cormat))
})

test_that("the correlated-factors group omega totals match psych::omega", {
  skip_on_cran()
  skip_if_not_installed("psych")
  # An external check on the composite common variance. psych scores the same data
  # through its own Schmid-Leiman solution and adds the general and own-group terms of
  # that solution, so its group totals sit a hair below ours, which count the
  # cross-loading terms too. Compared as sorted values, because the two packages order
  # their factors independently.
  res <- efa_reliability(efa_mod)
  ours <- res$value[res$coefficient == "omega_total" & res$level == "group"]
  po <- suppressWarnings(suppressMessages(
    psych::omega(test_models$baseline$cormat, nfactors = 3, n.obs = 500, fm = "pa",
                 rotate = "Promax", plot = FALSE)))
  theirs <- sort(unname(as.matrix(po$omega.group)[-1, "total"]))
  expect_equal(sort(ours), theirs, tolerance = 5e-3)
  # Agreement to three decimals is not the whole claim. The values must belong to the
  # factors they are reported for: each factor's total is its own composite's, not a
  # sorted set. Checked against the composite common variance of that factor's assigned
  # items -- the assertion that fails if the cross terms are dropped again, whatever the
  # numerical agreement with psych. The cross terms are also shown to matter here, so
  # the identity is not vacuous; their sign is not fixed in general, because
  # cs' Phi cs is a quadratic form and can fall either side of cs_j^2 depending on how
  # the other columns' sums line up with the factor correlations.
  spec <- .rel_adapt_efa(efa_mod, type = "psych")
  common <- spec$s_load %*% spec$Phi %*% t(spec$s_load)
  for (j in seq_len(ncol(spec$s_load))) {
    mem <- which(spec$map[, j] == 1)
    Vgr <- sum(spec$cormat[mem, mem])
    expect_equal(res$value[res$coefficient == "omega_total" &
                             res$factor == spec$fac_names[j]],
                 sum(common[mem, mem]) / Vgr)
    expect_gt(sum(common[mem, mem]) / Vgr, sum(spec$s_load[mem, j])^2 / Vgr)
  }
})

test_that("sums_load on an oblique EFA is scored, and accounts for the factor correlations", {
  # Both conventions are defined for a correlated-factors solution: they share the
  # composite's model-implied common variance, which reads the factor correlations, and
  # differ only in the denominator -- the observed total variance against the
  # model-implied one.
  expect_no_warning(efa_reliability(efa_mod))
  res <- efa_reliability(efa_mod, variance = "sums_load")
  expect_identical(attr(res, "settings"), list(variance = "sums_load",
                                               no_general = TRUE))

  spec <- .rel_adapt_efa(efa_mod, type = "psych")
  common <- spec$s_load %*% spec$Phi %*% t(spec$s_load)
  tot <- res$value[res$coefficient == "omega_total" & res$factor == "g"]
  expect_equal(tot, sum(common) / (sum(common) + sum(spec$u2)))

  # The factor correlations are genuinely in there: ignoring them, as summing the squared
  # loading-column sums alone would, gives a materially different value.
  s <- colSums(spec$s_load)
  expect_false(isTRUE(all.equal(sum(common), sum(s^2))))

  # This solution reproduces its correlation matrix closely, so the two denominators --
  # and hence the two conventions -- nearly agree; they are not the same computation.
  expect_equal(tot, suppressWarnings(efa_reliability(efa_mod))$value[1],
               tolerance = 1e-4)
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

test_that("a correlation matrix is rejected, a square bifactor matrix is scored", {
  # Read as a bifactor matrix, a correlation matrix takes the first variable's
  # correlations for general-factor loadings and returns plausible but meaningless
  # coefficients. It carries no loading class, so nothing but its shape identifies it.
  expect_error(efa_reliability(test_models$baseline$cormat),
               class = "efa_reliability_cormat_as_model")
  expect_error(efa_reliability(stats::cov2cor(stats::cov(GRiPS_raw))),
               class = "efa_reliability_cormat_as_model")

  # The guard keys on the unit diagonal, not on squareness: a bifactor matrix of p
  # variables with a general factor and p - 1 group factors is square too, and is still
  # scored, because its diagonal holds loadings rather than ones.
  L <- matrix(c(0.60, 0.60, 0.60, 0.60,
                0.40, 0.40, 0.00, 0.00,
                0.00, 0.00, 0.40, 0.40,
                0.30, 0.00, 0.00, 0.30), nrow = 4)
  expect_s3_class(efa_reliability(L), "efa_reliability")

  # The manual route takes the group loadings separately, where a correlation matrix
  # passes the class and numeric checks and has one row per variable, so nothing else
  # would catch it.
  cm <- test_models$baseline$cormat
  p <- ncol(cm)
  expect_error(
    efa_reliability(var_names = colnames(cm), g_load = rep(0.5, p), s_load = cm,
                    u2 = rep(0.5, p), variance = "sums_load"),
    class = "efa_reliability_cormat_as_s_load")

  # The group loadings of a genuine bifactor solution are still accepted there.
  s <- L[, 2:3]
  expect_s3_class(
    efa_reliability(var_names = paste0("V", 1:4), g_load = L[, 1], s_load = s,
                    u2 = 1 - L[, 1]^2 - rowSums(s^2), variance = "sums_load"),
    "efa_reliability")

  # A missing value is reported as such: the correlation-matrix test aborts on one in
  # terms of its own argument, so it must not be the check that speaks first.
  s_na <- test_models$baseline$cormat
  s_na[2, 3] <- s_na[3, 2] <- NA
  expect_error(
    efa_reliability(var_names = colnames(s_na), g_load = rep(0.5, ncol(s_na)),
                    s_load = s_na, u2 = rep(0.5, ncol(s_na)), variance = "sums_load"),
    class = "efa_invalid_argument")
  # Likewise on the matrix route, where the guard sees an incomplete matrix at all.
  expect_error(efa_reliability(s_na), class = "efa_invalid_argument")
})

test_that("the whole Schmid-Leiman table is refused as s_load, its group columns are not", {
  # Passed whole, its h2 and u2 columns would be scored as two more group factors and its
  # general column as a third, on top of the g_load named separately.
  expect_error(
    efa_reliability(var_names = rownames(sl_mod$sl), g_load = sl_mod$sl[, "g"],
                    s_load = sl_mod$sl, u2 = sl_mod$sl[, "u2"], variance = "sums_load"),
    class = "efa_reliability_s_load_display_cols")

  # Either column alone is caught too, which also renders the singular message.
  for (nm in c("h2", "u2")) {
    s1 <- sl_mod$sl[, c(nm, "F2")]
    expect_error(
      efa_reliability(var_names = rownames(sl_mod$sl), g_load = sl_mod$sl[, "g"],
                      s_load = s1, u2 = sl_mod$sl[, "u2"], variance = "sums_load"),
      class = "efa_reliability_s_load_display_cols")
  }

  # The group-factor columns of that same table are the documented input and still score.
  expect_s3_class(
    efa_reliability(var_names = rownames(sl_mod$sl), g_load = sl_mod$sl[, "g"],
                    s_load = sl_mod$sl[, c("F1", "F2", "F3")], u2 = sl_mod$sl[, "u2"],
                    variance = "sums_load"),
    "efa_reliability")
})

test_that("a bare loading matrix with a missing or infinite loading is refused", {
  # With no u2 the adapter derives it as 1 - rowSums(loadings^2), so an unchecked loading
  # carries a missing or infinite uniqueness into the coefficients -- and the table it
  # produces looks ordinary. The manual route refuses the same components.
  L <- matrix(c(0.6, 0.6, 0.6, 0.6,
                0.5, 0.5, 0.0, 0.0,
                0.0, 0.0, 0.5, 0.5), nrow = 4)
  for (bad in list(NA_real_, Inf)) {
    L_bad <- L
    L_bad[2, 2] <- bad
    expect_error(efa_reliability(L_bad, variance = "sums_load"),
                 class = "efa_invalid_argument")
  }
  expect_s3_class(efa_reliability(L, variance = "sums_load"), "efa_reliability")
})

test_that("a factor_map holding values other than 0 and 1 is refused", {
  # The map states membership, and is read both as the set of assigned variables and as a
  # multiplier on the loadings and uniquenesses. A map of 2s names the same composites but
  # empties every one of them, returning NA coefficients under a plausible whole-scale row.
  L <- matrix(c(0.6, 0.6, 0.6, 0.6,
                0.5, 0.5, 0.0, 0.0,
                0.0, 0.0, 0.5, 0.5), nrow = 4)
  map <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), nrow = 4)
  expect_error(efa_reliability(L, factor_map = map * 2, variance = "sums_load"),
               class = "efa_reliability_map_values")
  # A missing entry is no more a membership statement than a 2 is.
  map_na <- map; map_na[1, 1] <- NA
  expect_error(efa_reliability(L, factor_map = map_na, variance = "sums_load"),
               class = "efa_reliability_map_values")
  # Several distinct offending values, so the message lists more than one of them.
  map_mixed <- matrix(c(1, 2, 3, 0, 0, 0, 1, 1), nrow = 4)
  expect_error(efa_reliability(L, factor_map = map_mixed, variance = "sums_load"),
               class = "efa_reliability_map_values")

  # The check sits where every route meets, so OMEGA()'s factor_corres gets it too --
  # that front end validates only the dimensions of its own argument.
  expect_error(OMEGA(sl_mod, factor_corres = matrix(2, nrow(sl_mod$sl), 3)),
               class = "efa_reliability_map_values")

  # 0/1 and the logical map that stands for it are both accepted, and agree.
  res <- efa_reliability(L, factor_map = map, variance = "sums_load")
  expect_s3_class(res, "efa_reliability")
  expect_equal(res, efa_reliability(L, factor_map = map == 1, variance = "sums_load"))
})

test_that("an ordinary loading matrix is rejected, a Schmid-Leiman table is scored", {
  # A pattern matrix read as a bifactor matrix takes its first factor for a general
  # factor and returns plausible but meaningless coefficients, so it is refused.
  expect_error(efa_reliability(efa_mod$rot_loadings),
               class = "efa_reliability_pattern_matrix")
  expect_error(efa_reliability(efa_mod$unrot_loadings),
               class = "efa_reliability_pattern_matrix")
  # A genuine bifactor matrix is still read as one once it carries no loading class.
  expect_s3_class(efa_reliability(unclass(efa_mod$rot_loadings)), "efa_reliability")

  # The Schmid-Leiman loading table is a bifactor matrix plus its h2/u2 display
  # columns: it is scored without them, and with its own uniquenesses, so it gives
  # exactly what the Schmid-Leiman object gives.
  res <- efa_reliability(sl_mod$sl, factor_map = fc,
                         cormat = test_models$baseline$cormat)
  expect_setequal(unique(res$factor), c("g", "F1", "F2", "F3"))
  expect_equal(res, efa_reliability(sl_mod, factor_map = fc))

  # With no factor_map the table gets the same auto-derived correspondences as the
  # object. The nonzero-pattern default of a raw bifactor matrix would instead assign
  # every variable to every group factor -- the table is estimated and has no structural
  # zeros -- which collapses PUC to zero and repeats the whole-scale alpha on every row.
  expect_equal(efa_reliability(sl_mod$sl, cormat = test_models$baseline$cormat),
               efa_reliability(sl_mod))

  # The h2/u2 columns are always the last two, so an unlabelled table is read off that
  # layout instead of scoring them as two more group factors (which puts omega total
  # above 1 and names a factor "5").
  unlabelled <- sl_mod$sl
  colnames(unlabelled) <- NULL
  res_bare <- efa_reliability(unlabelled, cormat = test_models$baseline$cormat)
  expect_length(unique(res_bare$factor), 4L)
  expect_true(all(res_bare$value[res_bare$coefficient == "omega_total"] <= 1))
  # Only the row labels differ from the labelled table, not the coefficients.
  expect_equal(res_bare$value, res$value)
})

test_that("an ordinary loading matrix with its Phi is scored as a correlated-factors solution", {
  # The same solution in three spellings -- the fitted object, its loading matrix with the
  # factor correlations, and the components one by one -- must give one result. The matrix
  # route is the object's own components read off the matrix, so this is an identity and
  # not an approximation.
  L <- efa_mod$rot_loadings
  Phi <- efa_mod$Phi
  cm <- test_models$baseline$cormat
  u2 <- 1 - diag(L %*% Phi %*% t(L))

  res <- efa_reliability(L, Phi = Phi, cormat = cm)
  expect_s3_class(res, "efa_reliability")
  expect_equal(res, efa_reliability(efa_mod, cormat = cm))
  expect_equal(res, efa_reliability(NULL, var_names = rownames(L), g_load = rep(0, nrow(L)),
                                    s_load = unclass(L), u2 = u2, Phi = Phi, cormat = cm,
                                    fac_names = colnames(L)))

  # It is the correlated-factors menu, not the bifactor one: the coefficients a solution
  # without a general factor does not define are omitted rather than returned as zeros.
  expect_setequal(res$coefficient, c("omega_total", "omega_subscale", "alpha", "H"))
  expect_identical(attr(res, "settings"), list(variance = "correlation",
                                               no_general = TRUE))
  # Labels come from the matrix itself.
  expect_setequal(unique(res$factor), c("g", "F1", "F2", "F3"))

  # Without Phi nothing says the columns are correlated factors, so the matrix is still
  # refused rather than read as a bifactor one.
  expect_error(efa_reliability(L, cormat = cm), class = "efa_reliability_pattern_matrix")
  # And Phi is no longer reported as ignored on this route, while it still is on the
  # others -- an unclassed matrix is documented as a bifactor one, whose group factors are
  # orthogonal by construction, and reading its first column as an ordinary factor because
  # a second argument was given would be the silent reinterpretation the refusal prevents.
  expect_no_warning(efa_reliability(L, Phi = Phi, cormat = cm),
                    class = "efa_reliability_args_ignored")
  expect_warning(efa_reliability(unclass(L), Phi = Phi, cormat = cm),
                 class = "efa_reliability_args_ignored")

  # The uniquenesses are derived from the loadings under Phi, and a supplied set overrides
  # them, as on the bare bifactor matrix route. Shown in the model-implied convention,
  # which is the one they enter (with a correlation matrix in hand they do not).
  derived <- efa_reliability(L, Phi = Phi, variance = "sums_load")
  expect_equal(efa_reliability(L, Phi = Phi, variance = "sums_load", u2 = u2), derived)
  own <- suppressWarnings(efa_reliability(L, Phi = Phi, variance = "sums_load",
                                          u2 = u2 / 2))
  expect_false(isTRUE(all.equal(own$value, derived$value)))
  # A malformed one is refused there too, rather than scored into a plausible table.
  expect_error(efa_reliability(L, Phi = Phi, cormat = cm, u2 = replace(u2, 2, NA_real_)),
               class = "efa_invalid_argument")
  # A Phi that cannot be the correlation matrix of these columns is named as such.
  expect_error(efa_reliability(L, Phi = diag(2), cormat = cm),
               class = "efa_reliability_bad_phi")
})

test_that("the matrix route carries the factor correlations into the coefficients", {
  # Not a wiring test: the whole-scale omega total is the Phi-aware common variance over
  # the observed variance, and the Phi-blind reading is materially different.
  L <- efa_mod$rot_loadings
  Phi <- efa_mod$Phi
  cm <- test_models$baseline$cormat
  res <- efa_reliability(L, Phi = Phi, cormat = cm)
  tot <- res$value[res$coefficient == "omega_total" & res$factor == "g"]
  common <- unclass(L) %*% Phi %*% t(unclass(L))

  expect_equal(tot, sum(common) / sum(cm))
  expect_false(isTRUE(all.equal(tot, sum(tcrossprod(unclass(L))) / sum(cm))))

  # Model-implied variances need no correlation matrix; correlation mode does, and says so
  # rather than dividing by nothing.
  expect_s3_class(efa_reliability(L, Phi = Phi, variance = "sums_load"),
                  "efa_reliability")
  expect_error(efa_reliability(L, Phi = Phi), class = "efa_reliability_need_cormat")
})

test_that("a supplied cormat must hold the solution's variables", {
  cm <- test_models$baseline$cormat

  # The same variables in another order are reordered to the solution, so the
  # coefficients are the ones of the correctly ordered matrix rather than silently
  # wrong subscale values.
  expect_equal(efa_reliability(sl_mod, factor_map = fc, cormat = cm[18:1, 18:1]),
               efa_reliability(sl_mod, factor_map = fc, cormat = cm))
  # A different number of variables, or different ones, cannot be matched.
  expect_error(efa_reliability(sl_mod, factor_map = fc, cormat = cm[1:6, 1:6]),
               class = "efa_reliability_cormat_dim")
  other <- cm
  dimnames(other) <- list(paste0("X", 1:18), paste0("X", 1:18))
  expect_error(efa_reliability(sl_mod, factor_map = fc, cormat = other),
               class = "efa_reliability_cormat_names")
  # An unlabelled matrix of the right size is still matched by position.
  bare <- unname(cm)
  expect_equal(efa_reliability(sl_mod, factor_map = fc, cormat = bare),
               efa_reliability(sl_mod, factor_map = fc, cormat = cm))
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
  expect_identical(attr(res, "settings"), list(variance = "sums_load",
                                               no_general = FALSE))
  get <- function(coef, fac) res$value[res$coefficient == coef & res$factor == fac]
  expect_equal(get("omega_total", "g"), unname(om["g", "tot"]))
  expect_equal(get("H", "F1"), unname(om["F1", "H"]))
})

test_that("a lavaan correlated-factors fit is scored as the oblique solution it is", {
  skip_if_not_installed("lavaan")
  mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
          F2 =~ V7 + V8 + V9 + V10 + V11 + V12
          F3 =~ V13 + V14 + V15 + V16 + V17 + V18'
  fit <- suppressWarnings(lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
                                      sample.nobs = 500, estimator = "ml"))

  # How the fit was read is stated rather than left to be inferred from a missing column.
  expect_message(efa_reliability(fit), class = "efa_reliability_no_general_factor")
  res <- suppressMessages(efa_reliability(fit))

  # The correlated-factors menu: the coefficients a solution with no general factor does
  # not define are omitted, not returned as the zeros the arithmetic would give.
  expect_setequal(res$coefficient, c("omega_total", "omega_subscale", "alpha", "H"))
  expect_false(any(res$coefficient == "omega_subscale" & res$factor == "g"))
  expect_false(any(res$coefficient == "H" & res$factor == "g"))
  expect_true(all(is.na(res$group)))
  expect_identical(attr(res, "settings"), list(variance = "sums_load",
                                               no_general = TRUE))

  # The values are the model-implied ones, as for every lavaan input: the composite's
  # common variance Lambda Psi Lambda' over that variance plus the residual variances.
  std <- suppressWarnings(lavaan::lavInspect(fit, "std",
                                             drop.list.single.group = FALSE))[[1]]
  common <- std$lambda %*% std$psi %*% t(std$lambda)
  expect_equal(res$value[res$coefficient == "omega_total" & res$factor == "g"],
               sum(common) / (sum(common) + sum(diag(std$theta))))
  # Ignoring the factor correlations would give a materially different value, so the
  # identity above is not one the orthogonal reading would also satisfy.
  expect_false(isTRUE(all.equal(sum(common), sum(tcrossprod(std$lambda)))))

  # `g_name` names a general factor, which this fit has none of: it is read from the
  # structure, so naming one of the ordinary factors changes nothing.
  expect_equal(suppressMessages(efa_reliability(fit, g_name = "F1")), res)

  # A multiple-group fit is scored per group, as the supported shapes are.
  fit_gr <- suppressWarnings(lavaan::cfa(mod, sample.cov =
    list(test_models$baseline$cormat, test_models$baseline$cormat),
    sample.nobs = c(500, 500), estimator = "ml"))
  res_gr <- suppressMessages(efa_reliability(fit_gr, group_names = c("Young", "Old")))
  expect_identical(unique(res_gr$group), c("Young", "Old"))
  expect_setequal(res_gr$coefficient, c("omega_total", "omega_subscale", "alpha", "H"))

  # The printed explanation of the two omega columns has to describe the table it stands
  # above for a grouped result as well. Every factor label repeats once per group here, so
  # a comparison keyed on the label alone matches only each factor's first group and can
  # never find the columns equal -- and under exact simple structure they are equal.
  tot <- res_gr$value[res_gr$coefficient == "omega_total" & res_gr$level == "group"]
  sub_om <- res_gr$value[res_gr$coefficient == "omega_subscale" & res_gr$level == "group"]
  expect_equal(tot, sub_om)
  testthat::local_reproducible_output()
  txt <- paste(cli::ansi_strip(format(res_gr)), collapse = " ")
  expect_match(txt, "subscale omega equals its total omega")

  # OMEGA() keeps its own surface: its wide per-factor output reports the general
  # factor's omega hierarchical, ECV, and PUC, which this fit does not define.
  expect_error(OMEGA(fit), class = "efa_reliability_invalid_lavaan")
  expect_error(OMEGA(fit, g_name = "F1"), class = "efa_reliability_invalid_lavaan")
})

test_that("a fit that is not a correlated-factors one still needs its general factor named", {
  skip_if_not_installed("lavaan")
  # A variable on two or more factors is the bifactor structure, and a fit with one is
  # read as a bifactor solution -- so its general factor still has to be named, and a
  # correlated-factors model with a cross-loading is not routed round that.
  mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6 + V7
          F2 =~ V7 + V8 + V9 + V10 + V11 + V12
          F3 =~ V13 + V14 + V15 + V16 + V17 + V18'
  fit <- suppressWarnings(lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
                                      sample.nobs = 500, estimator = "ml"))
  expect_error(efa_reliability(fit), class = "efa_reliability_g_name")

  # And a genuine bifactor whose factors were left correlated stays refused: the general
  # and group parts of such a model do not partition a composite's variance, which is a
  # different problem from having no general factor and is named as one.
  mod_bi <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
             F2 =~ V7 + V8 + V9 + V10 + V11 + V12
             F3 =~ V13 + V14 + V15 + V16 + V17 + V18
             g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
                  V13 + V14 + V15 + V16 + V17 + V18'
  fit_bi <- suppressWarnings(lavaan::cfa(mod_bi, sample.cov = test_models$baseline$cormat,
                                         sample.nobs = 500, estimator = "ml",
                                         orthogonal = FALSE))
  expect_error(efa_reliability(fit_bi, g_name = "g"),
               class = "efa_reliability_correlated_factors")
})

test_that("a lavaan single-factor fit surfaces omega total, alpha, and H, and informs", {
  skip_if_not_installed("lavaan")
  fit <- suppressWarnings(lavaan::cfa(
    'g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
          V13 + V14 + V15 + V16 + V17 + V18',
    sample.cov = test_models$baseline$cormat, sample.nobs = 500,
    estimator = "ml", orthogonal = TRUE))

  expect_message(efa_reliability(fit, g_name = "g"),
                 class = "efa_reliability_single_factor")
  res <- suppressMessages(efa_reliability(fit, g_name = "g"))
  expect_setequal(res$coefficient, c("omega_total", "alpha", "H"))
  expect_true(all(res$factor == "g"))

  # The values are the model-implied ones every lavaan input takes: omega total from the
  # loading sum over that sum plus the residual variances, standardized alpha from the
  # model-implied correlation matrix, and the H index over the loadings.
  std <- suppressWarnings(lavaan::lavInspect(fit, "std",
                                             drop.list.single.group = FALSE))[[1]]
  lambda <- std$lambda[, 1]
  theta <- diag(std$theta)
  k <- length(lambda)
  R_imp <- stats::cov2cor(tcrossprod(lambda) + diag(theta))
  get <- function(coef) res$value[res$coefficient == coef]
  expect_equal(get("omega_total"), sum(lambda)^2 / (sum(lambda)^2 + sum(theta)))
  expect_equal(get("alpha"), k / (k - 1) * (1 - sum(diag(R_imp)) / sum(R_imp)))
  expect_equal(get("H"), 1 / (1 + 1 / sum(lambda^2 / (1 - lambda^2))))

  # Both expectations below are about the undefined-coefficient note alone, and every call
  # they make also raises the single-factor one, which would otherwise print outside any
  # expectation. Muffled by class rather than with suppressMessages(), which would hide the
  # note the expectations are about.
  quiet <- function(expr) {
    withCallingHandlers(
      expr,
      efa_reliability_single_factor = function(m) invokeRestart("muffleMessage"))
  }

  # Alpha assumes essential tau-equivalence, which is nested in a one-factor model, so it
  # is defined here and requesting it returns it instead of reporting it undefined.
  expect_no_message(quiet(efa_reliability(fit, g_name = "g", coefficients = "alpha")),
                    class = "efa_reliability_coef_undefined")
  sel <- suppressMessages(efa_reliability(fit, g_name = "g", coefficients = "alpha"))
  expect_identical(nrow(sel), 1L)
  expect_equal(sel$value, get("alpha"))

  # The coefficients a single factor genuinely does not define still report as undefined.
  expect_message(
    quiet(efa_reliability(fit, g_name = "g", coefficients = "omega_subscale")),
    class = "efa_reliability_coef_undefined")
})

test_that("every input route scores one and the same one-factor solution alike", {
  skip_if_not_installed("lavaan")

  # One solution, handed in four ways. It is scored against its own model-implied
  # correlation matrix, which a one-factor model reproduces exactly: lavaan then recovers
  # the loadings the EFA already has, and the observed and model-implied total variances
  # coincide, so the input route is the only difference left between the four results.
  cm <- test_models$baseline$cormat
  vn <- rownames(cm)
  fit_efa <- efa_fit(cm, N = 500, n_factors = 1, estimator = "PAF")
  L <- as.numeric(unclass(fit_efa$unrot_loadings))
  u2 <- 1 - L^2
  L_col <- matrix(L, ncol = 1, dimnames = list(vn, NULL))
  R_imp <- tcrossprod(L) + diag(u2)
  dimnames(R_imp) <- dimnames(cm)

  fit_lav <- suppressWarnings(lavaan::cfa(
    paste("g =~", paste(vn, collapse = " + ")),
    sample.cov = R_imp, sample.nobs = 500, estimator = "ml"))

  routes <- suppressMessages(list(
    efa_fit = efa_reliability(fit_efa, cormat = R_imp),
    matrix  = efa_reliability(L_col, cormat = R_imp),
    manual  = efa_reliability(NULL, var_names = vn, g_load = rep(0, length(L)),
                              s_load = L_col, u2 = u2, Phi = matrix(1, 1, 1),
                              cormat = R_imp),
    lavaan  = efa_reliability(fit_lav, g_name = "g")))

  # Every route reaches the same menu -- what one factor defines, and nothing else -- on the
  # one row a single-factor solution has, which is its general-factor row whatever the input
  # calls the factor.
  for (route in routes) {
    expect_setequal(route$coefficient, c("omega_total", "alpha", "H"))
    expect_true(all(route$level == "general"))
    expect_identical(length(unique(route$factor)), 1L)
  }

  # The label is the name the input gives the factor: the EFA's own column label, the name
  # lavaan's model syntax gives it, and -- for the two routes that name it nothing -- the
  # default first-factor label, which the EFA's happens to be as well. Not "g", which would
  # claim a general factor over group factors this solution has none of.
  expect_identical(unique(routes$efa_fit$factor),
                   colnames(unclass(fit_efa$unrot_loadings)))
  expect_identical(unique(routes$efa_fit$factor), "F1")
  expect_identical(unique(routes$matrix$factor), "F1")
  expect_identical(unique(routes$manual$factor), "F1")
  expect_identical(unique(routes$lavaan$factor), "g")   # the model calls it "g"

  # ... and the same values, not only the same coefficient names.
  vals <- function(x) stats::setNames(x$value, x$coefficient)
  expect_equal(vals(routes$matrix), vals(routes$efa_fit))
  expect_equal(vals(routes$manual), vals(routes$efa_fit))
  # lavaan re-estimates the solution rather than being handed it, so it agrees to the
  # precision of its own ML fit of a model this matrix satisfies exactly.
  expect_equal(vals(routes$lavaan), vals(routes$efa_fit), tolerance = 1e-6)

  # The agreement is on the right numbers: each coefficient is what its definition gives
  # for this solution, so three routes matching a fourth is not three routes sharing a bug.
  k <- length(L)
  expect_equal(vals(routes$efa_fit)[["omega_total"]], sum(L)^2 / (sum(L)^2 + sum(u2)))
  expect_equal(vals(routes$efa_fit)[["alpha"]],
               k / (k - 1) * (1 - sum(diag(R_imp)) / sum(R_imp)))
  expect_equal(vals(routes$efa_fit)[["H"]], 1 / (1 + 1 / sum(L^2 / (1 - L^2))))
})

test_that("a one-column loading matrix is scored, not refused as a pattern matrix", {
  # The refusal of an ordinary loading matrix asks the caller either for the factor
  # intercorrelations or for a hierarchy, and a solution with one factor has neither to
  # give: one factor has no correlations with other factors, and one column is no
  # hierarchy for it to be mistaken for. So a one-column matrix is read as the single
  # factor it is, whether or not it carries the loading class of the solution it came
  # from -- which is what a one-factor `efa_fit()` object's loadings do carry.
  cm <- test_models$baseline$cormat
  fit <- efa_fit(cm, N = 500, n_factors = 1, estimator = "PAF")
  L <- fit$unrot_loadings
  expect_s3_class(L, "efa_loadings")

  expect_message(efa_reliability(L, cormat = cm),
                 class = "efa_reliability_single_factor")
  res <- suppressMessages(efa_reliability(L, cormat = cm))
  expect_setequal(res$coefficient, c("omega_total", "alpha", "H"))
  expect_identical(unique(res$factor), "F1")

  # The class diverts nothing: the same matrix stripped of it takes the same route, and
  # the object the matrix came from reaches the same coefficients through its own adapter.
  expect_identical(res, suppressMessages(efa_reliability(unclass(L), cormat = cm)))
  expect_equal(res, suppressMessages(efa_reliability(fit)))
  # The superseded loading class is read the same way.
  L_old <- unclass(L)
  class(L_old) <- "LOADINGS"
  expect_identical(suppressMessages(efa_reliability(L_old, cormat = cm)), res)

  # Nor does it divert the model-implied correlation matrix the bare-matrix route falls
  # back to when `cormat` is absent, or either variance convention.
  expect_identical(suppressMessages(efa_reliability(L)),
                   suppressMessages(efa_reliability(unclass(L))))
  expect_identical(suppressMessages(efa_reliability(L, variance = "sums_load")),
                   suppressMessages(efa_reliability(unclass(L), variance = "sums_load")))

  # A single factor's Phi can only be the 1 x 1 identity, which says nothing about the
  # solution, so the route that reads the matrix as a correlated-factors one lands on the
  # same single factor and the same coefficients.
  expect_equal(suppressMessages(efa_reliability(L, Phi = matrix(1, 1, 1), cormat = cm)),
               res)

  # The exemption stops at one column: with two or more, nothing still says whether they
  # are a hierarchy or correlated factors, and the matrix is refused as before.
  two_col <- unclass(efa_mod$rot_loadings)[, 1:2, drop = FALSE]
  class(two_col) <- c("efa_loadings", "LOADINGS")
  expect_error(efa_reliability(two_col, cormat = cm),
               class = "efa_reliability_pattern_matrix")
})

test_that("a fac_names with the wrong number of names is refused", {
  # One label per factor. Too few or too many otherwise reach rownames() on the coefficient
  # matrix and fail there in base R's terms -- an error that names neither the argument nor
  # the count it needed, and that is localized on a non-English machine.
  expect_error(efa_reliability(sl_mod, factor_map = fc, fac_names = c("A", "B")),
               class = "efa_reliability_fac_names_length")
  expect_error(efa_reliability(sl_mod, factor_map = fc,
                               fac_names = c("A", "B", "C", "D")),
               class = "efa_reliability_fac_names_length")
  expect_error(efa_reliability(efa_mod, cormat = test_models$baseline$cormat,
                               fac_names = c("A", "B")),
               class = "efa_reliability_fac_names_length")
  # One name per group factor is what the solution has, so it is accepted and used.
  named <- efa_reliability(sl_mod, factor_map = fc, fac_names = c("A", "B", "C"))
  expect_setequal(unique(named$factor), c("g", "A", "B", "C"))

  # A single-factor solution has one factor to name, and the general-factor row is not one
  # of them -- so two names are as wrong here as they are above, and are refused rather
  # than dropped in favour of the default label.
  cm <- test_models$baseline$cormat
  fit1 <- efa_fit(cm, N = 500, n_factors = 1, estimator = "PAF")
  expect_error(suppressMessages(efa_reliability(fit1, fac_names = c("A", "B"))),
               class = "efa_reliability_fac_names_length")
  expect_error(suppressMessages(efa_reliability(fit1$unrot_loadings, cormat = cm,
                                                fac_names = character(0))),
               class = "efa_reliability_fac_names_length")
  expect_identical(unique(suppressMessages(
    efa_reliability(fit1, fac_names = "Impulsivity"))$factor), "Impulsivity")
})

test_that("fac_names labels the one factor of a single-factor solution", {
  cm <- test_models$baseline$cormat
  vn <- rownames(cm)
  L <- as.numeric(unclass(efa_fit(cm, N = 500, n_factors = 1,
                                  estimator = "PAF")$unrot_loadings))
  L_col <- matrix(L, ncol = 1, dimnames = list(vn, NULL))
  manual <- function(...) suppressMessages(efa_reliability(
    NULL, var_names = vn, g_load = rep(0, length(L)), s_load = L_col,
    u2 = 1 - L^2, cormat = cm, ...))

  named <- manual(fac_names = "Impulsivity")
  expect_identical(unique(named$factor), "Impulsivity")
  # Still the whole-scale row: there are no group factors for it to be distinguished from,
  # which is why omega subscale is not among the coefficients.
  expect_true(all(named$level == "general"))
  expect_setequal(named$coefficient, c("omega_total", "alpha", "H"))
  # Only the label changes with it.
  expect_equal(named$value, manual()$value)

  # Named nothing, the factor takes the default first-factor label rather than "g", which
  # would claim a general factor over group factors a single-factor solution has none of.
  # The manual route's own default is the bare column position, which is not a name.
  expect_identical(unique(manual()$factor), "F1")
  expect_identical(unique(suppressMessages(efa_reliability(L_col, cormat = cm))$factor), "F1")
  # A name that happens to be a digit is still a name.
  expect_identical(unique(manual(fac_names = "1")$factor), "1")

  # A one-column matrix names its factor on the only column it has, which is read as the
  # general one; that name reaches the label, so the matrix and the components agree.
  named_col <- L_col
  colnames(named_col) <- "Impulsivity"
  expect_identical(
    unique(suppressMessages(efa_reliability(named_col, cormat = cm))$factor),
    unique(named$factor))
})

test_that("a general factor with one group factor is two factors, not one", {
  # The single-factor reading keys on how many factors a solution has, not on how many
  # columns any one slot holds: a bifactor solution with a single group factor still has
  # two, and keeps the coefficients that decompose a composite into them.
  L <- cbind(g = rep(0.5, 6), F1 = c(rep(0.4, 3), rep(0, 3)))
  rownames(L) <- paste0("V", 1:6)
  expect_no_message(efa_reliability(L), class = "efa_reliability_single_factor")
  res <- efa_reliability(L)
  expect_setequal(unique(res$factor), c("g", "F1"))
  expect_true(any(res$coefficient == "omega_hierarchical"))
  expect_true(any(res$coefficient == "omega_subscale"))
})

test_that("a one-factor efa_fit() is scored, whether or not a rotation was asked for", {
  cm <- test_models$baseline$cormat
  fit <- efa_fit(cm, N = 500, n_factors = 1, estimator = "PAF")

  expect_message(efa_reliability(fit), class = "efa_reliability_single_factor")
  res <- suppressMessages(efa_reliability(fit))
  expect_setequal(res$coefficient, c("omega_total", "alpha", "H"))

  # A single factor cannot be rotated, so asking for a rotation changes neither the
  # solution nor its coefficients -- the unrotated loadings are read either way.
  rot <- suppressWarnings(suppressMessages(
    efa_fit(cm, N = 500, n_factors = 1, estimator = "PAF", rotation = "promax")))
  expect_equal(suppressMessages(efa_reliability(rot)), res)

  # Scored against the object's own correlation matrix, as an oblique solution is.
  L <- as.numeric(unclass(fit$unrot_loadings))
  expect_equal(res$value[res$coefficient == "omega_total"],
               sum(tcrossprod(L)) / sum(cm))

  # A solution with more than one factor still needs an oblique rotation to be scored,
  # which is advice a one-factor solution could not follow.
  orth <- suppressWarnings(suppressMessages(
    efa_fit(cm, N = 500, n_factors = 3, estimator = "PAF", rotation = "varimax")))
  expect_error(efa_reliability(orth), class = "efa_reliability_not_oblique")
})

test_that("the one-factor menu does not depend on the variance convention", {
  # The two conventions are two total-variance denominators, not two coefficient menus:
  # they move omega total and leave the set of coefficients a single factor defines alone.
  fit <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 1, estimator = "PAF")
  corr <- suppressMessages(efa_reliability(fit, variance = "correlation"))
  sums <- suppressMessages(efa_reliability(fit, variance = "sums_load"))

  expect_setequal(corr$coefficient, sums$coefficient)
  expect_setequal(corr$coefficient, c("omega_total", "alpha", "H"))
  # Not vacuous: the denominators genuinely differ, the solution not reproducing the
  # observed correlations exactly.
  expect_false(isTRUE(all.equal(corr$value[corr$coefficient == "omega_total"],
                                sums$value[sums$coefficient == "omega_total"])))
})

test_that("the single-factor alpha is standardized alpha", {
  # An external check on the one coefficient the single-factor menu gained. Unlike the
  # psych::omega comparison above this needs no factor solution from psych -- alpha is a
  # closed-form function of the correlation matrix -- so it is cheap and stable enough to
  # run everywhere.
  skip_if_not_installed("psych")
  cm <- test_models$baseline$cormat
  fit <- efa_fit(cm, N = 500, n_factors = 1, estimator = "PAF")
  res <- suppressMessages(efa_reliability(fit))
  expect_equal(res$value[res$coefficient == "alpha"],
               suppressWarnings(psych::alpha(cm, n.obs = 500))$total$std.alpha,
               tolerance = 1e-8)
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

test_that("a multiple-group single-factor fit yields one single-factor block per group", {
  skip_if_not_installed("lavaan")
  # The only path that assembles several single-factor matrices into one result: each block
  # has to be labelled, blanked, and scored on its own, and the note stated once for all.
  cm <- test_models$baseline$cormat
  fit <- suppressWarnings(lavaan::cfa(
    paste("ability =~", paste(rownames(cm), collapse = " + ")),
    sample.cov = list(cm, cm), sample.nobs = c(500, 500), estimator = "ml"))

  expect_message(efa_reliability(fit, group_names = c("Young", "Old")),
                 class = "efa_reliability_single_factor")
  res <- suppressMessages(efa_reliability(fit, group_names = c("Young", "Old")))

  expect_identical(unique(res$group), c("Young", "Old"))
  for (grp in c("Young", "Old")) {
    block <- res[res$group == grp, , drop = FALSE]
    expect_setequal(block$coefficient, c("omega_total", "alpha", "H"))
    # The factor's name in the model, on the whole-scale row.
    expect_true(all(block$factor == "ability"))
    expect_true(all(block$level == "general"))
  }
  # Identical groups, so the two blocks must agree coefficient by coefficient.
  expect_equal(res$value[res$group == "Young"], res$value[res$group == "Old"])
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

test_that("efa_reliability is invariant across the oblique rotations", {
  skip_if_not_installed("GPArotation")

  # The factor columns are reordered by the number in their labels. A rotation whose
  # labels are absent, or carry no number, must therefore fall back to the columns'
  # own order: an unguarded reorder subscripts the loadings down to zero columns, and
  # the whole-scale omega total then silently becomes a function of the correlation
  # matrix alone while every per-factor row disappears. The whole-scale omega total is
  # rotation-invariant, so all oblique rotations must reach the same one.
  ref <- efa_reliability(efa_mod)
  g_tot <- function(x) x$value[x$coefficient == "omega_total" & x$factor == "g"]

  res_list <- lapply(.oblq_rotations, function(rot) {
    fit <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
                   estimator = "PAF", rotation = rot)
    # bifactorQ leaves a group factor without salient items on this fixture, which is
    # reported by the classed efa_omega_empty_factor warning and is not under test here.
    suppressWarnings(efa_reliability(fit))
  })
  names(res_list) <- .oblq_rotations

  expect_true(all(vapply(res_list, nrow, integer(1)) > 2))
  expect_equal(vapply(res_list, g_tot, numeric(1)),
               stats::setNames(rep(g_tot(ref), length(res_list)), names(res_list)),
               tolerance = 1e-8)

  # An oblimin fit additionally keeps the full per-factor structure of a promax fit
  # of the same extraction.
  res_obl <- res_list[["oblimin"]]
  expect_setequal(res_obl$factor, c("g", "F1", "F2", "F3"))
  for (cf in c("omega_subscale", "H")) {
    expect_setequal(res_obl$factor[res_obl$coefficient == cf], c("F1", "F2", "F3"))
  }
  expect_setequal(res_obl$coefficient, unique(ref$coefficient))
})

## keying -------------------------------------------------------------------

# One factor's block of variables reverse-coded, which is what a scale whose
# reverse-worded items were never reverse-coded looks like. It exercises both shapes the
# check has to cover at once: in a Schmid-Leiman solution the block's general-factor
# loadings turn negative, whereas the oblique solution reflects the factor column instead
# and comes back with positive loadings throughout and a negative factor correlation.
flip <- rep(1, 18)
flip[7:12] <- -1
cormat_flipped <- diag(flip) %*% test_models$baseline$cormat %*% diag(flip)
dimnames(cormat_flipped) <- dimnames(test_models$baseline$cormat)
efa_flipped <- suppressWarnings(EFA(cormat_flipped, N = 500, n_factors = 3,
                                    type = "EFAtools", method = "PAF",
                                    rotation = "promax"))
sl_flipped <- suppressWarnings(efa_schmid_leiman(
  efa_flipped, estimate_control = estimate_control(type = "EFAtools"),
  estimator = "PAF"))

test_that("the fixture really is mis-keyed in both of the shapes the check must cover", {
  # Guards the tests below: were the fixture to lose either shape, they would still pass
  # while covering only one of the two routes.
  expect_true(any(sl_flipped$sl[, 1] < 0))
  expect_true(any(efa_flipped$Phi[upper.tri(efa_flipped$Phi)] < 0))
})

test_that("variables keyed against each other warn on every non-lavaan input route", {
  cls <- "efa_reliability_mixed_keying"
  sl <- unclass(sl_flipped$sl)
  g <- sl[, 1]
  s <- sl[, 2:4, drop = FALSE]
  u2 <- sl[, "u2"]

  expect_warning(efa_reliability(sl_flipped), class = cls)
  # The reverse-keyed block also leaves one of this fixture's subscale composites with a
  # smaller observed variance than the loadings imply, which trips the independent
  # out-of-range guard; muffle it so this expectation is about the keying alone.
  expect_warning(
    withCallingHandlers(efa_reliability(efa_flipped),
                        efa_omega_out_of_range = function(w) invokeRestart("muffleWarning")),
    class = cls)
  expect_warning(efa_reliability(cbind(g = g, s), u2 = u2, cormat = cormat_flipped),
                 class = cls)
  expect_warning(efa_reliability(NULL, var_names = rownames(sl), g_load = g,
                                 s_load = s, u2 = u2, cormat = cormat_flipped),
                 class = cls)

  skip_if_not_installed("psych")
  sm <- suppressWarnings(psych::schmid(cormat_flipped, nfactors = 3, fm = "minres",
                                       rotate = "promax"))
  expect_warning(efa_reliability(sm), class = cls)
})

test_that("a lavaan fit of mis-keyed variables warns, single-factor path included", {
  skip_if_not_installed("lavaan")
  cls <- "efa_reliability_mixed_keying"
  items <- paste0("V", 1:18)

  bi <- suppressWarnings(lavaan::cfa(
    paste('F1 =~ V1 + V2 + V3 + V4 + V5 + V6
           F2 =~ V7 + V8 + V9 + V10 + V11 + V12
           F3 =~ V13 + V14 + V15 + V16 + V17 + V18
           g =~', paste(items, collapse = " + ")),
    sample.cov = cormat_flipped, sample.nobs = 500, estimator = "ml",
    orthogonal = TRUE))
  expect_warning(efa_reliability(bi, g_name = "g"), class = cls)

  # The single-factor path reaches the core on a spec with no group factors, which no
  # other input to the check has, so it needs its own coverage.
  sf <- suppressWarnings(lavaan::cfa(paste("g =~", paste(items, collapse = " + ")),
                                     sample.cov = cormat_flipped, sample.nobs = 500,
                                     estimator = "ml"))
  expect_warning(suppressMessages(efa_reliability(sf, g_name = "g")), class = cls)
})

test_that("a correctly keyed solution raises no keying warning", {
  cls <- "efa_reliability_mixed_keying"
  expect_no_warning(efa_reliability(sl_mod, factor_map = fc), class = cls)
  expect_no_warning(efa_reliability(efa_mod), class = cls)
  expect_no_warning(efa_reliability(sl_mod$sl), class = cls)

  skip_if_not_installed("lavaan")
  items <- paste0("V", 1:18)
  sf <- suppressWarnings(lavaan::cfa(paste("g =~", paste(items, collapse = " + ")),
                                     sample.cov = test_models$baseline$cormat,
                                     sample.nobs = 500, estimator = "ml"))
  expect_no_warning(suppressMessages(efa_reliability(sf, g_name = "g")), class = cls)
})

test_that("the keying warning changes no coefficient", {
  # The check is a diagnostic only; the mis-keyed solution's coefficients are still the
  # ones the definitions give.
  res <- suppressWarnings(efa_reliability(sl_flipped))
  get <- function(x, coef, fac) x$value[x$coefficient == coef & x$factor == fac]

  # The whole-scale omega total is still the composite's model-implied common variance
  # over its observed variance -- collapsed by the keying, not by the warning.
  spec <- .rel_adapt_SL(sl_flipped, type = "psych")
  L <- cbind(spec$g_load, spec$s_load)
  expect_equal(get(res, "omega_total", "g"), sum(tcrossprod(L)) / sum(spec$cormat))

  # The same solution scored through the core directly, which is where the check sits, is
  # identical to what the front end returns.
  core <- .reliability_result(
    suppressWarnings(.reliability_core(spec, "correlation", add_ind = TRUE,
                                       add_rel = TRUE, arg = "factor_map")),
    settings = list(variance = "correlation", no_general = FALSE))
  expect_equal(res, core)
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
  # blank rather than as "NA", and the context line separates tot from sub.
  expect_snapshot(print(efa_reliability(efa_mod)), transform = scrub_num)
})

test_that("the correlated-factors context line follows the printed columns", {
  # It relates the tot and sub columns, so it is shown only when both are printed --
  # not for a selection that drops omega subscale, nor for a solution that has a
  # general factor.
  expect_match(format(efa_reliability(efa_mod)), "Correlated-factors", all = FALSE)
  sub_dropped <- suppressMessages(
    efa_reliability(efa_mod, coefficients = c("omega_total", "alpha")))
  expect_false(any(grepl("Correlated-factors", format(sub_dropped))))
  expect_false(any(grepl("Correlated-factors",
                         format(efa_reliability(sl_mod, factor_map = fc)))))

  # Which relation it states follows the printed values. An oblique fit has factor
  # correlations, so its totals exceed its subscale omegas.
  obl <- format(efa_reliability(efa_mod))
  expect_match(obl, "receives from every factor", all = FALSE)
  expect_false(any(grepl("equals its total omega", obl)))

  # Components supplied with a zero general factor, no factor correlations and simple
  # structure are a correlated-factors solution too, but nothing lets a composite draw
  # on another factor, so the two columns coincide and the line must say so instead.
  mk <- function(s) {
    u <- 1 - rowSums(s^2)
    efa_reliability(model = NULL, var_names = paste0("V", 1:6),
                    g_load = rep(0, 6), s_load = s, u2 = u,
                    cormat = s %*% t(s) + diag(u))
  }
  s <- matrix(0, 6, 2)
  s[1:3, 1] <- c(0.7, 0.6, 0.5)
  s[4:6, 2] <- c(0.7, 0.6, 0.5)
  no_phi <- mk(s)
  tot <- no_phi$value[no_phi$coefficient == "omega_total" & no_phi$level == "group"]
  sub_om <- no_phi$value[no_phi$coefficient == "omega_subscale"]
  expect_equal(tot, sub_om)
  txt <- format(no_phi)
  expect_match(txt, "equals its total omega", all = FALSE)
  expect_false(any(grepl("receives from every factor", txt)))

  # Cross-loadings alone are enough to separate the two columns, even with no factor
  # correlations at all: a composite still receives true score variance from the other
  # factor through them. The line must then take the general wording.
  s_cross <- s
  s_cross[3, 2] <- 0.3
  s_cross[6, 1] <- 0.25
  cross <- mk(s_cross)
  tot_c <- cross$value[cross$coefficient == "omega_total" & cross$level == "group"]
  sub_c <- cross$value[cross$coefficient == "omega_subscale"]
  expect_true(all(tot_c > sub_c))
  txt_c <- format(cross)
  expect_match(txt_c, "receives from every factor", all = FALSE)
  expect_false(any(grepl("equals its total omega", txt_c)))
})

rm(efa_mod, sl_mod, fc, auto, flip, cormat_flipped, efa_flipped, sl_flipped)
