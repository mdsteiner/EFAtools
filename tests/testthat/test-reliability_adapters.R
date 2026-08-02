## Tests for the .rel_adapt_* front-end adapters. Each normalizes one input source
## to the spec consumed by .reliability_core(); the spec must drive the core to the
## same numbers the current OMEGA path produces for that input.

efa_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
sl_mod <- efa_schmid_leiman(efa_mod, estimate_control = estimate_control(type = "EFAtools"), estimator = "PAF")
schmid_mod <- psych::schmid(test_models$baseline$cormat, nfactors = 3,
                            n.obs = 500, fm = "pa", rotate = "Promax")
fc <- sl_mod$sl[, c("F1", "F2", "F3")] >= .2

test_that(".rel_adapt_SL reproduces the SL OMEGA path", {
  spec <- .rel_adapt_SL(sl_mod, factor_corres = fc, type = "EFAtools")

  expect_identical(
    .reliability_core(spec, "correlation"),
    .OMEGA_FLEX(sl_mod, type = "EFAtools", factor_corres = fc,
                variance = "correlation"))
  expect_identical(
    .reliability_core(spec, "sums_load"),
    .OMEGA_FLEX(sl_mod, type = "EFAtools", factor_corres = fc,
                variance = "sums_load"))
  expect_identical(
    .reliability_core(spec, "correlation", add_ind = FALSE),
    .OMEGA_FLEX(sl_mod, type = "EFAtools", factor_corres = fc,
                variance = "correlation", add_ind = FALSE))

  # The stored correlation matrix is picked up from the SL object.
  expect_identical(spec$cormat, sl_mod$orig_R)
})

test_that(".rel_adapt_schmid reproduces the psych::schmid OMEGA path", {
  # Cormat reconstructed from the schmid pattern and Phi.
  expect_equal(
    .reliability_core(.rel_adapt_schmid(schmid_mod, type = "psych"), "correlation"),
    .OMEGA_FLEX(schmid_mod, type = "psych", variance = "correlation"))

  # User-supplied cormat.
  expect_equal(
    .reliability_core(
      .rel_adapt_schmid(schmid_mod, type = "psych",
                        cormat = test_models$baseline$cormat), "correlation"),
    .OMEGA_FLEX(schmid_mod, type = "psych", cormat = test_models$baseline$cormat,
                variance = "correlation"))
})

test_that(".rel_adapt_manual reproduces the model = NULL OMEGA path", {
  spec <- .rel_adapt_manual(g_load = sl_mod$sl[, "g"],
                            s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                            u2 = sl_mod$sl[, "u2"],
                            var_names = rownames(sl_mod$sl),
                            factor_corres = fc, type = "EFAtools",
                            cormat = test_models$baseline$cormat)

  expect_identical(
    .reliability_core(spec, "correlation"),
    .OMEGA_FLEX(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
                g_load = sl_mod$sl[, "g"], s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                u2 = sl_mod$sl[, "u2"], cormat = test_models$baseline$cormat,
                factor_corres = fc, variance = "correlation"))
  expect_identical(
    .reliability_core(spec, "sums_load"),
    .OMEGA_FLEX(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
                g_load = sl_mod$sl[, "g"], s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                u2 = sl_mod$sl[, "u2"], cormat = test_models$baseline$cormat,
                factor_corres = fc, variance = "sums_load"))
})

test_that(".rel_adapt_efa gives model-appropriate correlated-factors coefficients", {
  spec <- .rel_adapt_efa(efa_mod)

  # No general factor: g_load is zero, so omega hierarchical and ECV vanish.
  expect_true(all(spec$g_load == 0))
  expect_no_warning(om <- .reliability_core(spec, "correlation"))
  expect_true(all(om[, "hier"] == 0))
  expect_equal(unname(om["g", "ECV"]), 0)

  # Whole-scale omega total is 1 - sum(u2) / 1'R1 (all common variance over total).
  expect_equal(unname(om["g", "tot"]),
               1 - sum(spec$u2) / sum(spec$cormat))

  # Each group row is that factor's congeneric omega from its own indicators.
  for (j in seq_len(ncol(spec$s_load))) {
    mem <- which(spec$map[, j] == 1)
    cong <- sum(spec$s_load[mem, j])^2 / sum(spec$cormat[mem, mem])
    expect_equal(unname(om[spec$fac_names[j], "sub"]), cong)
    expect_equal(unname(om[spec$fac_names[j], "hier"]), 0)
  }

  # Works for a two-factor oblique EFA too (no Schmid-Leiman, so no identification
  # floor at three factors).
  efa2 <- EFA(test_models$baseline$cormat, N = 500, n_factors = 2,
              type = "EFAtools", method = "PAF", rotation = "promax")
  expect_no_error(.reliability_core(.rel_adapt_efa(efa2), "correlation"))
})

test_that(".rel_adapt_efa handles rotations that leave the loading columns unlabelled", {
  skip_if_not_installed("GPArotation")

  # efa_fit() labels the factor columns for every rotation, but the adapter also sees
  # loading matrices from other sources (a bare pattern matrix, an externally built fit),
  # where the labels can be missing. Unlabelled columns must not cost the adapter its
  # group factors, so construct that state explicitly.
  efa_obl <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
                     estimator = "PAF", rotation = "oblimin")
  colnames(efa_obl$rot_loadings) <- NULL
  dimnames(efa_obl$Phi) <- NULL

  spec <- .rel_adapt_efa(efa_obl)

  # All three factors survive, and every group factor still gets a row label.
  expect_identical(dim(spec$s_load), c(18L, 3L))
  expect_identical(spec$fac_names, c("F1", "F2", "F3"))
  expect_identical(spec$s_load, unclass(efa_obl$rot_loadings))

  # Communalities come from the real pattern, not from an empty loading matrix.
  expect_true(all(spec$u2 > 0 & spec$u2 < 1))

  om <- .reliability_core(spec, "correlation")
  expect_identical(rownames(om), c("g", "F1", "F2", "F3"))
  # Each group factor scores a real congeneric omega. Tested on the group rows
  # only: the g row carries a "sub" value whether or not the group factors
  # survived, so it cannot discriminate a lost factor.
  expect_true(all(om[c("F1", "F2", "F3"), "sub"] > 0 &
                    om[c("F1", "F2", "F3"), "sub"] < 1))

  # Whole-scale omega total is 1 - sum(u2) / 1'R1, and L Phi L' -- hence u2 -- is
  # invariant to the choice of oblique rotation, so the promax fit of the same
  # extraction must agree on it.
  om_pro <- .reliability_core(.rel_adapt_efa(efa_mod), "correlation")
  expect_equal(unname(om["g", "tot"]), unname(om_pro["g", "tot"]),
               tolerance = 1e-8)
})

test_that(".rel_adapt_bifactor reproduces the manual OMEGA path", {
  # A raw orthogonal-bifactor loading matrix: g in the first column.
  g_load <- rep(0.5, 6)
  s_load <- matrix(0, 6, 2)
  s_load[1:3, 1] <- 0.4
  s_load[4:6, 2] <- 0.4
  L <- cbind(g = g_load, s_load)
  colnames(L) <- c("g", "F1", "F2")
  rownames(L) <- paste0("V", 1:6)
  u2 <- 1 - rowSums(L^2)
  Rm <- L %*% t(L) + diag(u2)
  map <- abs(s_load) > 0

  spec <- .rel_adapt_bifactor(L, fac_names = c("F1", "F2"))
  expect_identical(spec$u2, u2)
  expect_identical(spec$cormat, Rm)

  ref_corr <- .OMEGA_FLEX(model = NULL, type = "EFAtools", var_names = rownames(L),
                          g_load = g_load, s_load = s_load, u2 = u2, cormat = Rm,
                          factor_corres = map, fac_names = c("F1", "F2"),
                          variance = "correlation")
  ref_sums <- .OMEGA_FLEX(model = NULL, type = "EFAtools", var_names = rownames(L),
                          g_load = g_load, s_load = s_load, u2 = u2,
                          factor_corres = map, fac_names = c("F1", "F2"),
                          variance = "sums_load")

  expect_identical(.reliability_core(spec, "correlation"), ref_corr)
  expect_identical(.reliability_core(spec, "sums_load"), ref_sums)
})

test_that("adapters reject a non-correlation cormat via .is_cormat", {
  bad <- matrix(rnorm(36), 6, 6)
  expect_error(.rel_adapt_SL(sl_mod, factor_corres = fc, cormat = bad),
               class = "efa_reliability_not_cormat")
  expect_error(.rel_adapt_bifactor(cbind(rep(0.5, 6), matrix(0.4, 6, 1)),
                                   cormat = bad),
               class = "efa_reliability_not_cormat")
  # the guard runs before the coercion, so a data frame that is not a correlation matrix is
  # still rejected rather than silently coerced
  expect_error(.rel_adapt_SL(sl_mod, factor_corres = fc, cormat = as.data.frame(bad)),
               class = "efa_reliability_not_cormat")
})

test_that("a cormat supplied as a data frame is accepted and handed on as a matrix", {
  # read.csv(file, row.names = 1) returns a data frame, and this guard is the only check a
  # user-supplied cormat passes on the reliability route -- it does not go through
  # .prepare_cor_input(). Without the coercion the omega columns still come out right (they
  # only sum and subset the matrix) and only the standardized-alpha branch fails, on diag();
  # efa_reliability() always requests that branch, so the whole call breaks.
  cm <- test_models$baseline$cormat
  cdf <- as.data.frame(cm)

  expect_true(is.matrix(.rel_check_cormat(cdf)))
  expect_identical(.rel_check_cormat(cdf), cm)

  # the adapters, and the exported entry point, give the matrix route's answer
  expect_equal(.rel_adapt_SL(sl_mod, factor_corres = fc, cormat = cdf),
               .rel_adapt_SL(sl_mod, factor_corres = fc, cormat = cm))
  expect_equal(efa_reliability(efa_mod, cormat = cdf),
               efa_reliability(efa_mod, cormat = cm))
})

test_that(".rel_adapt_efa requires an oblique solution", {
  efa_orth <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
                  type = "EFAtools", method = "PAF", rotation = "varimax")
  expect_error(.rel_adapt_efa(efa_orth), class = "efa_reliability_not_oblique")
})

## lavaan adapters -----------------------------------------------------------

test_that(".rel_adapt_lavaan reproduces the bifactor OMEGA path", {
  skip_if_not_installed("lavaan")
  mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
          F2 =~ V7 + V8 + V9 + V10 + V11 + V12
          F3 =~ V13 + V14 + V15 + V16 + V17 + V18
          g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
               V13 + V14 + V15 + V16 + V17 + V18'
  fit <- suppressWarnings(lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
                                      sample.nobs = 500, estimator = "ml",
                                      orthogonal = TRUE))

  adapt <- .rel_adapt_lavaan(fit, g_name = "g")
  expect_length(adapt$groups, 1L)
  expect_false(adapt$groups[[1]]$single)
  expect_identical(adapt$variance, "sums_load")

  expect_equal(.reliability_core(adapt$groups[[1]], "sums_load"),
               .OMEGA_LAVAAN(fit, g_name = "g"))
})

test_that(".rel_adapt_lavaan Schmid-Leiman transforms a second-order model", {
  skip_if_not_installed("lavaan")
  mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
          F2 =~ V7 + V8 + V9 + V10 + V11 + V12
          F3 =~ V13 + V14 + V15 + V16 + V17 + V18
          g =~ F1 + F2 + F3'
  fit <- suppressWarnings(lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
                                      sample.nobs = 500, estimator = "ml"))

  spec <- .rel_adapt_lavaan(fit, g_name = "g")$groups[[1]]
  expect_equal(.reliability_core(spec, "sums_load"),
               suppressMessages(.OMEGA_LAVAAN(fit, g_name = "g")))
})

test_that(".rel_adapt_lavaan handles a single factor as a special path", {
  skip_if_not_installed("lavaan")
  mod <- 'g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
               V13 + V14 + V15 + V16 + V17 + V18'
  fit <- suppressWarnings(lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
                                      sample.nobs = 500, estimator = "ml",
                                      orthogonal = TRUE))

  sf <- .rel_adapt_lavaan(fit)$groups[[1]]
  expect_true(sf$single)

  ref <- suppressMessages(.OMEGA_LAVAAN(fit))  # c(Omega, H)
  omega <- sum(sf$g_load)^2 / (sum(sf$g_load)^2 + sum(sf$u2))
  expect_equal(unname(omega), unname(ref["Omega"]))
  expect_equal(unname(.h_index(sf$g_load)), unname(ref["H"]))
})

test_that(".rel_adapt_lavaan returns one spec per group", {
  skip_if_not_installed("lavaan")
  mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
          F2 =~ V7 + V8 + V9 + V10 + V11 + V12
          F3 =~ V13 + V14 + V15 + V16 + V17 + V18
          g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
               V13 + V14 + V15 + V16 + V17 + V18'
  fit <- suppressWarnings(lavaan::cfa(mod, sample.cov =
                                        list(test_models$baseline$cormat,
                                             test_models$baseline$cormat),
                                      sample.nobs = c(500, 500),
                                      estimator = "ml", orthogonal = TRUE))

  adapt <- .rel_adapt_lavaan(fit, g_name = "g", group_names = c("Some", "Others"))
  expect_identical(adapt$group_names, c("Some", "Others"))
  expect_length(adapt$groups, 2L)

  ref <- .OMEGA_LAVAAN(fit, g_name = "g", group_names = c("Some", "Others"))
  expect_equal(.reliability_core(adapt$groups[[1]], "sums_load"), ref$Some)
  expect_equal(.reliability_core(adapt$groups[[2]], "sums_load"), ref$Others)
})

rm(efa_mod, sl_mod, schmid_mod, fc)
