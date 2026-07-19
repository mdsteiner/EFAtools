## Tests for .reliability_core() -- the reliability engine behind OMEGA().

test_that(".reliability_core reproduces the OMEGA coefficient menu (regression)", {
  # Reference coefficients captured from OMEGA's flexible engine (both variance
  # modes), pinned to a tight tolerance. A regression in the lifted math would move
  # these numbers far beyond it; the comparison is not bit-exact because the last bits
  # differ across BLAS/platforms.
  ref_corr <- structure(c(0.8828685009816174, 0.76916616970747365,
    0.76340460593899595, 0.74403137565562616, 0.73987835131636559,
    0.49973769215022423, 0.4935927140365538, 0.51896424927096363,
    0.12494360637362421, 0.26942847755724941, 0.26981189190244209,
    0.22506712638466256, 0.8423884892261212, 0.4630872563777037,
    0.47192399221710424, 0.40799477228556474, 0.65197317077447814, NA, NA,
    NA, 0.70588235294117641, NA, NA, NA), dim = c(4L, 6L),
    dimnames = list(c("g", "1", "2", "3"),
                    c("tot", "hier", "sub", "H", "ECV", "PUC")), class = "OMEGA")
  ref_sums <- structure(c(0.88286138027838279, 0.7694136764855295,
    0.76539916799807128, 0.74540130311749742, 0.73992333028627189,
    0.49989850066069758, 0.49488233331890324, 0.51991977802965228,
    0.14293804999211088, 0.26951517582483187, 0.27051683467916798,
    0.22548152508784508, 0.8423884892261212, 0.4630872563777037,
    0.47192399221710424, 0.40799477228556474, 0.65197317077447814, NA, NA,
    NA, 0.70588235294117641, NA, NA, NA), dim = c(4L, 6L),
    dimnames = list(c("g", "1", "2", "3"),
                    c("tot", "hier", "sub", "H", "ECV", "PUC")), class = "OMEGA")
  # The SL adapter labels the group rows with the factor names; same numbers.
  ref_sl <- ref_corr
  rownames(ref_sl) <- c("g", "F1", "F2", "F3")

  efa_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
                 type = "EFAtools", method = "PAF", rotation = "promax")
  sl_mod <- efa_schmid_leiman(efa_mod, estimate_control = estimate_control(type = "EFAtools"), estimator = "PAF")
  fc <- sl_mod$sl[, c("F1", "F2", "F3")] >= .2
  spec <- list(g_load = sl_mod$sl[, "g"], s_load = sl_mod$sl[, c("F1", "F2", "F3")],
               u2 = sl_mod$sl[, "u2"], map = fc, cormat = test_models$baseline$cormat,
               var_names = rownames(sl_mod$sl), fac_names = seq_len(3))

  # Both variance modes, driven directly through the core.
  expect_equal(.reliability_core(spec, "correlation", TRUE), ref_corr, tolerance = 1e-6)
  spec_sums <- spec
  spec_sums$cormat <- NULL
  expect_equal(.reliability_core(spec_sums, "sums_load", TRUE), ref_sums, tolerance = 1e-6)

  # add_ind = FALSE drops the H/ECV/PUC columns but leaves the omegas unchanged.
  noadd <- .reliability_core(spec, "correlation", FALSE)
  expect_equal(unclass(noadd), unclass(ref_corr)[, c("tot", "hier", "sub")],
               tolerance = 1e-6)
  expect_s3_class(noadd, "OMEGA")

  # The .OMEGA_FLEX adapter must normalise an SL object to a spec that yields the
  # same numbers through the core.
  expect_equal(.OMEGA_FLEX(sl_mod, type = "EFAtools", factor_corres = fc,
                           variance = "correlation"), ref_sl, tolerance = 1e-6)
})

test_that("a loading of absolute value >= 1 returns that construct's H index as NA and notes it", {
  s <- matrix(0, 6, 2)
  s[1:3, 1] <- 0.4
  s[4:6, 2] <- 0.4
  u2 <- rep(0.5, 6)                          # proper unique variances: not a Heywood case
  map <- matrix(0, 6, 2)
  map[1:3, 1] <- 1
  map[4:6, 2] <- 1

  # A loading >= 1 on the general factor (item 6) NAs only the general-factor H. With
  # proper unique variances this is not a Heywood case (e.g. an oblique pattern
  # coefficient above 1), so it is reported as an undefined H, not a Heywood case.
  spec_g <- list(g_load = c(rep(0.5, 5), 1.05), s_load = s, u2 = u2, map = map,
                 cormat = NULL, var_names = paste0("V", 1:6),
                 fac_names = c("F1", "F2"))
  expect_warning(om_g <- .reliability_core(spec_g, "sums_load", TRUE),
                 class = "efa_reliability_h_undefined")
  expect_true(is.na(om_g["g", "H"]))
  expect_false(is.na(om_g["F1", "H"]))
  expect_false(is.na(om_g["F2", "H"]))

  # A loading >= 1 on a single group factor (item 4 on F2) NAs only that H.
  s2 <- s
  s2[4, 2] <- 1.1
  spec_s <- list(g_load = rep(0.5, 6), s_load = s2, u2 = u2, map = map, cormat = NULL,
                 var_names = paste0("V", 1:6), fac_names = c("F1", "F2"))
  expect_warning(om_s <- .reliability_core(spec_s, "sums_load", TRUE),
                 class = "efa_reliability_h_undefined")
  expect_true(is.na(om_s["F2", "H"]))
  expect_false(is.na(om_s["F1", "H"]))
  expect_false(is.na(om_s["g", "H"]))
})

test_that("a Heywood case (a unique variance at or below zero) is flagged from the variances", {
  s <- matrix(0, 6, 2)
  s[1:3, 1] <- 0.4
  s[4:6, 2] <- 0.4
  map <- matrix(0, 6, 2)
  map[1:3, 1] <- 1
  map[4:6, 2] <- 1

  # Item 1 has a non-positive unique variance (a communality at or above 1): an improper
  # (Heywood) solution, detected from the uniquenesses rather than the loadings.
  u2_bad <- c(-0.05, rep(0.5, 5))
  spec_bad <- list(g_load = rep(0.5, 6), s_load = s, u2 = u2_bad, map = map,
                   cormat = NULL, var_names = paste0("V", 1:6),
                   fac_names = c("F1", "F2"))
  expect_warning(.reliability_core(spec_bad, "sums_load", TRUE),
                 class = "efa_reliability_heywood")

  # Proper loadings (|lambda| < 1) with proper variances warn about neither.
  spec_ok <- list(g_load = rep(0.5, 6), s_load = s, u2 = rep(0.5, 6), map = map,
                  cormat = NULL, var_names = paste0("V", 1:6),
                  fac_names = c("F1", "F2"))
  expect_no_warning(.reliability_core(spec_ok, "sums_load", TRUE))
})

test_that("an admissible spec yields a finite H index and no Heywood warning", {
  spec <- list(g_load = rep(0.6, 6), s_load = matrix(0.3, 6, 1), u2 = rep(0.4, 6),
               map = matrix(1, 6, 1), cormat = NULL,
               var_names = paste0("V", 1:6), fac_names = "F1")
  expect_no_warning(om <- .reliability_core(spec, "sums_load", TRUE))
  expect_true(all(is.finite(om[, "H"])))
})

test_that("the h-undefined note is per-construct, not globally suppressed by a Heywood case", {
  collect_classes <- function(spec) {
    cls <- character(0)
    withCallingHandlers(
      .reliability_core(spec, "sums_load", TRUE),
      warning = function(w) { cls <<- c(cls, class(w)[1]); invokeRestart("muffleWarning") })
    cls
  }

  # Item 1 is a Heywood case (u2 < 0). Group factor F2's item 4 has a proper oblique pattern
  # loading >= 1 (u2 > 0), so F2's H is NA for a non-Heywood reason. Both are reported: the
  # Heywood case (naming V1) AND the undefined H for F2.
  s <- matrix(0, 6, 2); s[1:3, 1] <- 0.4; s[4:6, 2] <- 0.4; s[4, 2] <- 1.2
  map <- matrix(0, 6, 2); map[1:3, 1] <- 1; map[4:6, 2] <- 1
  spec <- list(g_load = rep(0.5, 6), s_load = s, u2 = c(-0.05, rep(0.5, 5)), map = map,
               cormat = NULL, var_names = paste0("V", 1:6), fac_names = c("F1", "F2"))
  cls <- collect_classes(spec)
  expect_true("efa_reliability_heywood" %in% cls)
  expect_true("efa_reliability_h_undefined" %in% cls)

  # A loading >= 1 that is ITSELF the Heywood case (negative u2 on the same item) is reported
  # only as a Heywood case, not additionally as an undefined H (no redundant double warning).
  s2 <- matrix(0, 6, 2); s2[1:3, 1] <- 0.4; s2[4:6, 2] <- 0.4
  spec_hw <- list(g_load = c(rep(0.5, 5), 1.2), s_load = s2, u2 = c(rep(0.5, 5), -0.1),
                  map = map, cormat = NULL, var_names = paste0("V", 1:6),
                  fac_names = c("F1", "F2"))
  cls_hw <- collect_classes(spec_hw)
  expect_true("efa_reliability_heywood" %in% cls_hw)
  expect_false("efa_reliability_h_undefined" %in% cls_hw)
})

test_that("an NA uniqueness is not misread as a Heywood case", {
  # An NA in u2 (a malformed component) must not fire a false Heywood warning naming the
  # variable "NA", nor suppress the legitimate h-undefined note. which(u2 < eps) drops NA.
  s <- matrix(0, 6, 2); s[1:3, 1] <- 0.4; s[4:6, 2] <- 0.4
  map <- matrix(0, 6, 2); map[1:3, 1] <- 1; map[4:6, 2] <- 1
  spec <- list(g_load = rep(0.5, 6), s_load = s, u2 = c(0.5, NA, rep(0.5, 4)), map = map,
               cormat = NULL, var_names = paste0("V", 1:6), fac_names = c("F1", "F2"))
  expect_no_warning(.reliability_core(spec, "sums_load", TRUE))

  # A genuine |loading| >= 1 still triggers the h-undefined note even when an unrelated
  # uniqueness is NA (the NA does not count as a Heywood variable that would suppress it).
  spec_h <- spec
  spec_h$g_load <- c(rep(0.5, 5), 1.05)
  expect_warning(.reliability_core(spec_h, "sums_load", TRUE),
                 class = "efa_reliability_h_undefined")
})

test_that("a missing (NA) loading yields an NA H index rather than an error", {
  # A missing loading is not a Heywood case: it is excluded from the >= 1 test and
  # flows through the sum to NA, so H must be NA (not an aborted call). This guards
  # against `if (any(abs(NA) >= 1))` erroring with "missing value where TRUE/FALSE".
  s <- matrix(0, 6, 2)
  s[1:3, 1] <- 0.4
  s[4:6, 2] <- 0.4
  map <- matrix(0, 6, 2)
  map[1:3, 1] <- 1
  map[4:6, 2] <- 1
  spec <- list(g_load = c(0.5, NA, rep(0.5, 4)), s_load = s, u2 = rep(0.5, 6),
               map = map, var_names = paste0("V", 1:6), fac_names = c("F1", "F2"))
  expect_no_warning(om <- .reliability_core(spec, "sums_load", TRUE))
  expect_true(is.na(om["g", "H"]))
  expect_true(is.finite(om["F1", "H"]))
  expect_true(is.finite(om["F2", "H"]))
})

## add_rel: Cronbach's alpha, composite reliability (CR), and AVE -----------------

test_that("add_rel = FALSE leaves the OMEGA coefficient columns unchanged", {
  # The OMEGA path never sets add_rel, so its coefficient menu must be untouched;
  # add_rel = TRUE appends exactly alpha, CR, and AVE.
  spec <- list(g_load = rep(0.4, 6), s_load = matrix(0.5, 6, 1), u2 = rep(0.59, 6),
               map = matrix(1, 6, 1), cormat = NULL, var_names = paste0("V", 1:6),
               fac_names = "F1")
  expect_identical(colnames(.reliability_core(spec, "sums_load", TRUE)),
                   c("tot", "hier", "sub", "H", "ECV", "PUC"))
  expect_identical(colnames(.reliability_core(spec, "sums_load", FALSE)),
                   c("tot", "hier", "sub"))
  expect_identical(colnames(.reliability_core(spec, "sums_load", TRUE, add_rel = TRUE)),
                   c("tot", "hier", "sub", "H", "ECV", "PUC", "alpha", "CR", "AVE"))
})

test_that("standardized alpha matches psych::alpha for the whole scale and subscales", {
  skip_if_not_installed("psych")
  R <- test_models$baseline$cormat            # 18 items, three 6-item blocks
  s <- matrix(0, 18, 3)
  s[1:6, 1] <- 0.5; s[7:12, 2] <- 0.5; s[13:18, 3] <- 0.5
  map <- matrix(0, 18, 3)
  map[1:6, 1] <- 1; map[7:12, 2] <- 1; map[13:18, 3] <- 1
  # variance = "sums_load" so the omega denominators ignore the (arbitrary) loadings;
  # only the appended alpha uses the supplied correlation matrix R.
  spec <- list(g_load = rep(0.4, 18), s_load = s, u2 = rep(0.2, 18), map = map,
               cormat = R, var_names = rownames(R), fac_names = c("F1", "F2", "F3"))
  om <- .reliability_core(spec, "sums_load", add_ind = FALSE, add_rel = TRUE)

  # Our alpha is standardized (computed from the correlation matrix), matching psych's
  # std.alpha; raw (covariance-based) alpha would need the item standard deviations.
  a_whole <- suppressWarnings(
    psych::alpha(R, check.keys = FALSE, warnings = FALSE))$total[["std.alpha"]]
  expect_equal(unname(om["g", "alpha"]), a_whole, tolerance = 1e-8)
  for (j in 1:3) {
    mem <- which(map[, j] == 1)
    a_sub <- suppressWarnings(
      psych::alpha(R[mem, mem], check.keys = FALSE, warnings = FALSE))$total[["std.alpha"]]
    expect_equal(unname(om[c("F1", "F2", "F3")[j], "alpha"]), a_sub, tolerance = 1e-8)
  }
})

test_that("alpha falls back to the model-implied correlation matrix when cormat is NULL", {
  g <- rep(0.4, 6)
  s <- matrix(0, 6, 2); s[1:3, 1] <- 0.5; s[4:6, 2] <- 0.5
  # Uniquenesses that do NOT complete to unit item variance (item 1), so the
  # model-implied Lambda Lambda' + diag(u2) has a non-unit, non-constant diagonal:
  # standardized alpha must standardize it (cov2cor) rather than read it as a
  # correlation matrix as-is.
  u2 <- c(0.30, rep(1 - 0.4^2 - 0.5^2, 5))    # item 1 diagonal = 0.16 + 0.25 + 0.30 != 1
  map <- matrix(0, 6, 2); map[1:3, 1] <- 1; map[4:6, 2] <- 1
  spec <- list(g_load = g, s_load = s, u2 = u2, map = map, cormat = NULL,
               var_names = paste0("V", 1:6), fac_names = c("F1", "F2"))
  om <- .reliability_core(spec, "sums_load", add_ind = FALSE, add_rel = TRUE)

  # Reference: standardize the model-implied matrix to a correlation matrix first.
  R <- stats::cov2cor(cbind(g, s) %*% t(cbind(g, s)) + diag(u2))
  expect_equal(unname(om["g", "alpha"]), 6 / 5 * (1 - sum(diag(R)) / sum(R)))
  Rk <- R[1:3, 1:3]
  expect_equal(unname(om["F1", "alpha"]), 3 / 2 * (1 - sum(diag(Rk)) / sum(Rk)))
})

test_that("CR and AVE match hand-computed values; the general-factor row is NA", {
  g <- rep(0.4, 6)
  s <- matrix(0, 6, 2); s[1:3, 1] <- 0.5; s[4:6, 2] <- 0.5
  u2 <- rep(1 - 0.4^2 - 0.5^2, 6)             # 0.59
  map <- matrix(0, 6, 2); map[1:3, 1] <- 1; map[4:6, 2] <- 1
  spec <- list(g_load = g, s_load = s, u2 = u2, map = map, cormat = NULL,
               var_names = paste0("V", 1:6), fac_names = c("F1", "F2"))
  om <- .reliability_core(spec, "sums_load", add_ind = FALSE, add_rel = TRUE)

  eF1 <- 3 * 0.59
  expect_equal(unname(om["F1", "CR"]),  1.5^2 / (1.5^2 + eF1))
  expect_equal(unname(om["F1", "AVE"]), (3 * 0.25) / (3 * 0.25 + eF1))
  expect_equal(unname(om["F2", "CR"]),  unname(om["F1", "CR"]))   # symmetry
  expect_equal(unname(om["F2", "AVE"]), unname(om["F1", "AVE"]))

  # CR and AVE are per-subscale only; the general-factor row is NA. Alpha is not.
  expect_true(is.na(om["g", "CR"]))
  expect_true(is.na(om["g", "AVE"]))
  expect_false(is.na(om["g", "alpha"]))
})

test_that("CR uses the assigned composite: a foreign item's cross-loading does not change it", {
  s0 <- matrix(0, 6, 2); s0[1:3, 1] <- 0.5; s0[4:6, 2] <- 0.5
  map <- matrix(0, 6, 2); map[1:3, 1] <- 1; map[4:6, 2] <- 1
  base <- list(g_load = rep(0.3, 6), s_load = s0, u2 = 1 - 0.3^2 - rowSums(s0^2),
               map = map, cormat = NULL, var_names = paste0("V", 1:6),
               fac_names = c("F1", "F2"))
  om0 <- .reliability_core(base, "sums_load", add_ind = FALSE, add_rel = TRUE)

  # Item 4 (assigned to F2) gains a 0.3 cross-loading on F1. F1's assigned composite
  # (items 1-3) is unchanged, so CR_F1 and AVE_F1 must not move.
  sX <- s0; sX[4, 1] <- 0.3
  cross <- base; cross$s_load <- sX; cross$u2 <- 1 - 0.3^2 - rowSums(sX^2)
  omX <- .reliability_core(cross, "sums_load", add_ind = FALSE, add_rel = TRUE)

  expect_equal(unname(omX["F1", "CR"]),  unname(om0["F1", "CR"]))
  expect_equal(unname(omX["F1", "AVE"]), unname(om0["F1", "AVE"]))

  # It equals the assigned-only formula, not a full-column one that would fold in the 0.3.
  eF1 <- 3 * (1 - 0.3^2 - 0.25)
  expect_equal(unname(omX["F1", "CR"]), 1.5^2 / (1.5^2 + eF1))
  expect_false(isTRUE(all.equal(unname(omX["F1", "CR"]), 1.8^2 / (1.8^2 + eF1))))
})

test_that("CR equals the lavaan single-factor omega for the k = 1 case", {
  skip_if_not_installed("lavaan")
  mod <- 'g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
               V13 + V14 + V15 + V16 + V17 + V18'
  fit <- suppressWarnings(lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
                                      sample.nobs = 500, estimator = "ml",
                                      orthogonal = TRUE))
  ref <- suppressMessages(.OMEGA_LAVAAN(fit))          # c(Omega, H)

  # Represent the single factor as one group factor (no general factor) and read CR.
  sf <- .rel_adapt_lavaan(fit)$groups[[1]]
  p <- length(sf$g_load)
  spec <- list(g_load = rep(0, p), s_load = matrix(sf$g_load, ncol = 1), u2 = sf$u2,
               map = matrix(1, p, 1), cormat = NULL, var_names = sf$var_names,
               fac_names = "F1")
  om <- .reliability_core(spec, "sums_load", add_ind = FALSE, add_rel = TRUE)
  expect_equal(unname(om["F1", "CR"]), unname(ref[["Omega"]]), tolerance = 1e-8)
})

test_that("empty and single-item subscales get NA alpha/CR/AVE as appropriate", {
  # Empty F2: all items assigned to F1, so F2 has an all-zero map column.
  s <- matrix(0, 6, 2); s[, 1] <- 0.5
  map <- matrix(0, 6, 2); map[, 1] <- 1
  spec <- list(g_load = rep(0.3, 6), s_load = s, u2 = 1 - 0.3^2 - rowSums(s^2),
               map = map, cormat = NULL, var_names = paste0("V", 1:6),
               fac_names = c("F1", "F2"))
  expect_warning(
    om <- .reliability_core(spec, "sums_load", add_ind = FALSE, add_rel = TRUE),
    class = "efa_omega_empty_factor")
  expect_true(all(is.na(om["F2", c("alpha", "CR", "AVE")])))
  expect_false(any(is.na(om["F1", c("alpha", "CR", "AVE")])))

  # Single-item subscale F1: alpha is undefined (NA), CR and AVE are finite.
  s2 <- matrix(0, 4, 2); s2[1, 1] <- 0.6; s2[2:4, 2] <- 0.5
  map2 <- matrix(0, 4, 2); map2[1, 1] <- 1; map2[2:4, 2] <- 1
  spec2 <- list(g_load = rep(0.3, 4), s_load = s2, u2 = 1 - 0.3^2 - rowSums(s2^2),
                map = map2, cormat = NULL, var_names = paste0("V", 1:4),
                fac_names = c("F1", "F2"))
  om2 <- .reliability_core(spec2, "sums_load", add_ind = FALSE, add_rel = TRUE)
  expect_true(is.na(om2["F1", "alpha"]))
  expect_true(is.finite(om2["F1", "CR"]))
  expect_true(is.finite(om2["F1", "AVE"]))
})
