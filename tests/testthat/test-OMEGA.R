## Use with a lavaan output
if (requireNamespace("lavaan", quietly = TRUE)) {
lav_mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
        F2 =~ V7 + V8 + V9 + V10 + V11 + V12
        F3 =~ V13 + V14 + V15 + V16 + V17 + V18
        g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
             V13 + V14 + V15 + V16 + V17 + V18'
lav_fit <- lavaan::cfa(lav_mod, sample.cov = test_models$baseline$cormat,
                   sample.nobs = 500, estimator = "ml", orthogonal = TRUE)
om_lav <- OMEGA(lav_fit, g_name = "g")

# Second-order model (g loads on the first-order factors): exercises the
# Schmid-Leiman transform of a higher-order solution.
lav_mod_ho <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
               F2 =~ V7 + V8 + V9 + V10 + V11 + V12
               F3 =~ V13 + V14 + V15 + V16 + V17 + V18
               g =~ F1 + F2 + F3'
lav_fit_ho <- suppressWarnings(lavaan::cfa(lav_mod_ho,
                                           sample.cov = test_models$baseline$cormat,
                                           sample.nobs = 500, estimator = "ml"))
om_lav_ho <- suppressMessages(OMEGA(lav_fit_ho, g_name = "g"))
}

## Use with an output from the SL function, with type EFAtools
efa_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
sl_mod <- SL(efa_mod, type = "EFAtools", method = "PAF")
om_sl <- OMEGA(sl_mod, type = "EFAtools",
               factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2)

## Use with an output from the psych::schmid function, with type psych
schmid_mod <- psych::schmid(test_models$baseline$cormat, nfactors = 3,
                            n.obs = 500, fm = "pa", rotate = "Promax")
# Find correlation matrix from phi and pattern matrix from psych::schmid outpu
om_schmid <- OMEGA(schmid_mod, type = "psych")

## Manually specify components
om_man <- OMEGA(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
                g_load = sl_mod$sl[, "g"], s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                u2 = sl_mod$sl[, "u2"], cormat = test_models$baseline$cormat,
                factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2)

test_that("output class and dimensions are correct", {
  skip_if_not_installed("lavaan")
  expect_s3_class(om_lav, "OMEGA")
  expect_s3_class(om_sl, "OMEGA")
  expect_s3_class(om_schmid, "OMEGA")
  expect_s3_class(om_man, "OMEGA")

  expect_output(str(om_lav), "List of 2")
  expect_output(str(om_sl), "List of 2")
  expect_output(str(om_schmid), "List of 2")
  expect_output(str(om_man), "List of 2")
})

test_that("a second-order lavaan model is Schmid-Leiman transformed", {
  skip_if_not_installed("lavaan")

  # The second-order general factor triggers the SL transform of the higher-order
  # solution; OMEGA informs the user and returns the full coefficient matrix.
  expect_message(OMEGA(lav_fit_ho, g_name = "g"),
                 class = "efa_omega_g_second_order")
  expect_s3_class(om_lav_ho, "OMEGA")
  expect_identical(dim(om_lav_ho), c(4L, 6L))
  expect_identical(rownames(om_lav_ho), c("g", "F1", "F2", "F3"))
  expect_identical(colnames(om_lav_ho),
                   c("tot", "hier", "sub", "H", "ECV", "PUC"))

  # The transformed group-factor loadings yield finite omega coefficients within
  # [0, 1], and omega hierarchical never exceeds omega total.
  coefs <- unclass(om_lav_ho)[, c("tot", "hier", "sub")]
  expect_true(all(is.finite(coefs)))
  expect_true(all(coefs >= 0 & coefs <= 1))
  expect_true(all(om_lav_ho[, "tot"] >= om_lav_ho[, "hier"]))
})

test_that("OMEGA lavaan numbers are unchanged (regression)", {
  skip_if_not_installed("lavaan")

  # Reference coefficient outputs captured from OMEGA()'s lavaan path, which scores
  # the fitted solution through the shared reliability core. These pin the coefficients
  # (non-Heywood fixtures) to a tight tolerance so any change to that path is caught;
  # the comparison is not bit-exact because the last bits differ across BLAS/platforms.
  ref_bi <- structure(c(0.88301174174790287, 0.74751625677111733, 0.766678791034203,
    0.7685498467792552, 0.76548235222658689, 0.5703940426328894,
    0.50146589411449449, 0.49415766808792488, 0.11752938952131604,
    0.17712221413822785, 0.26521289691970851, 0.2743921786913302,
    0.84884383208407344, 0.37901193583533627, 0.48206335211458812,
    0.47292561296446434, 0.67182974097583426, NA, NA, NA, 0.70588235294117641,
    NA, NA, NA), dim = c(4L, 6L), dimnames = list(c("g", "F1", "F2", "F3"),
    c("tot", "hier", "sub", "H", "ECV", "PUC")), class = "OMEGA")
  expect_equal(om_lav, ref_bi, tolerance = 1e-6)

  ref_ho <- structure(c(0.88158341868254952, 0.74400564775764777, 0.76359153637480304,
    0.76786220097481273, 0.76381778912616816, 0.54868574924161106,
    0.50423868000554961, 0.50487682924077948, 0.11776562955638133,
    0.19531989851603676, 0.25935285636925343, 0.26298537173403319,
    0.84773834907453449, 0.36084418350415226, 0.44810509596935005,
    0.45435071956970075, 0.68379181434818159, NA, NA, NA, 0.70588235294117641,
    NA, NA, NA), dim = c(4L, 6L), dimnames = list(c("g", "F1", "F2", "F3"),
    c("tot", "hier", "sub", "H", "ECV", "PUC")), class = "OMEGA")
  expect_equal(om_lav_ho, ref_ho, tolerance = 1e-6)

  # A single factor returns a named c(Omega, H) vector.
  ref_sf <- structure(c(Omega = 0.86802682659256869, H = 0.86873279471693676),
                      class = "OMEGA")
  fit_sf <- suppressWarnings(lavaan::cfa(
    'g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
          V13 + V14 + V15 + V16 + V17 + V18',
    sample.cov = test_models$baseline$cormat, sample.nobs = 500,
    estimator = "ml", orthogonal = TRUE))
  expect_equal(suppressMessages(OMEGA(fit_sf, g_name = "g")), ref_sf, tolerance = 1e-6)

  # Multiple groups return a named list of (unclassed) coefficient matrices.
  ref_mg <- structure(list(A = structure(c(0.88301175860369074, 0.74751650041121109,
    0.7666785723706202, 0.76854962402215088, 0.7654826788236242, 0.57038939781173303,
    0.50147399385562386, 0.49415517812176823, 0.11752907978006645, 0.17712710259947803,
    0.26520457851499635, 0.27439444590038276, 0.84884398691526919, 0.37901393509779763,
    0.48206053742597404, 0.47292785665908216, 0.67183010049739167, NA, NA, NA,
    0.70588235294117641, NA, NA, NA), dim = c(4L, 6L), dimnames = list(
    c("g", "F1", "F2", "F3"), c("tot", "hier", "sub", "H", "ECV", "PUC"))),
    B = structure(c(0.88301177654016794, 0.74751702091877537, 0.76667849734450322,
    0.76854951437076557, 0.76548249348611708, 0.57038915952811098, 0.50147394437097303,
    0.49415464215133981, 0.1175292830540509, 0.17712786139066439, 0.26520455297353013,
    0.27439487221942582, 0.84884395119146006, 0.37901520165768471, 0.48206079699355442,
    0.47292821672600055, 0.6718296371036121, NA, NA, NA, 0.70588235294117641,
    NA, NA, NA), dim = c(4L, 6L), dimnames = list(c("g", "F1", "F2", "F3"),
    c("tot", "hier", "sub", "H", "ECV", "PUC")))), class = "OMEGA")
  fit_mg <- suppressWarnings(lavaan::cfa(lav_mod,
    sample.cov = list(test_models$baseline$cormat, test_models$baseline$cormat),
    sample.nobs = c(500, 500), estimator = "ml", orthogonal = TRUE))
  # Looser than the single-group pins above: these coefficients come out of lavaan's
  # multigroup ML optimizer, which settles on a slightly different point under a
  # different BLAS and moves them by about 1e-6. Drift in the omega arithmetic itself
  # would move them by orders of magnitude more.
  expect_equal(OMEGA(fit_mg, g_name = "g", group_names = c("A", "B")), ref_mg,
               tolerance = 1e-4)
})

test_that("the general-factor row is labelled with g_name, not a hardcoded 'g'", {
  skip_if_not_installed("lavaan")
  mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
          F2 =~ V7 + V8 + V9 + V10 + V11 + V12
          F3 =~ V13 + V14 + V15 + V16 + V17 + V18
          gen =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
                 V13 + V14 + V15 + V16 + V17 + V18'
  fit <- suppressWarnings(lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
                                      sample.nobs = 500, estimator = "ml",
                                      orthogonal = TRUE))
  expect_identical(rownames(OMEGA(fit, g_name = "gen")),
                   c("gen", "F1", "F2", "F3"))
})

test_that("errors are thrown correctly", {
  expect_error(OMEGA(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
                     g_load = sl_mod$sl[, "g"], s_load = 1:7,
                     u2 = sl_mod$sl[, "u2"], cormat = test_models$baseline$cormat,
                     factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2), class = "efa_omega_bad_s_load")
  expect_error(OMEGA(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
                     g_load = sl_mod$sl[, "g"],
                     s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                     u2 = sl_mod$sl[, "u2"],
                     factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                     pattern = 1:5), class = "efa_omega_bad_pattern")
  expect_warning(OMEGA(sl_mod, type = "EFAtools", g_load = sl_mod$sl[, "g"],
                       s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                       u2 = sl_mod$sl[, "u2"],
                       factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2), class = "efa_omega_model_args_ignored")
  expect_error(OMEGA(model = 1:4), class = "efa_omega_invalid_model")
  expect_error(OMEGA(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
                     s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                     u2 = sl_mod$sl[, "u2"],
                     factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2), class = "efa_omega_missing_args")
})

test_that("a user cormat is honoured for a flexible-input SL object", {
  # A flexible-input SL object stores orig_R = NA; OMEGA must use the supplied
  # cormat rather than overwriting it (and must not crash).
  sl_flex <- SL(efa_mod$rot_loadings, Phi = efa_mod$Phi, type = "EFAtools",
                method = "PAF")
  expect_true(identical(sl_flex$orig_R, NA))

  fc <- sl_flex$sl[, c("F1", "F2", "F3")] >= .2
  om_flex <- OMEGA(sl_flex, type = "EFAtools",
                   cormat = test_models$baseline$cormat, factor_corres = fc)
  expect_s3_class(om_flex, "OMEGA")

  # The coefficients match the equivalent model = NULL manual path with the same
  # loadings, uniquenesses, and correlation matrix (the SL object carries factor
  # names, so only the row labels differ).
  om_flex_ref <- OMEGA(model = NULL, type = "EFAtools",
                       var_names = rownames(sl_flex$sl),
                       g_load = sl_flex$sl[, "g"],
                       s_load = sl_flex$sl[, c("F1", "F2", "F3")],
                       u2 = sl_flex$sl[, "u2"],
                       cormat = test_models$baseline$cormat, factor_corres = fc)
  expect_equal(unname(unclass(om_flex)), unname(unclass(om_flex_ref)))
})

test_that("a solution whose variables are keyed against each other is flagged", {
  # The coefficients describe the unit-weighted sum of the variables as supplied, so a
  # scale whose reverse-worded items were never reverse-coded is signalled rather than
  # returned as a merely poor one. Shared with efa_reliability() through the reliability
  # engine, so it reaches this surface too.
  flip <- rep(1, 18)
  flip[7:12] <- -1
  cormat_flipped <- diag(flip) %*% test_models$baseline$cormat %*% diag(flip)
  dimnames(cormat_flipped) <- dimnames(test_models$baseline$cormat)
  efa_flipped <- suppressWarnings(EFA(cormat_flipped, N = 500, n_factors = 3,
                                      type = "EFAtools", method = "PAF",
                                      rotation = "promax"))
  sl_flipped <- suppressWarnings(SL(efa_flipped, type = "EFAtools", method = "PAF"))
  fc <- sl_flipped$sl[, c("F1", "F2", "F3")] >= .2

  expect_warning(OMEGA(sl_flipped, type = "EFAtools", factor_corres = fc),
                 class = "efa_reliability_mixed_keying")
  expect_no_warning(OMEGA(sl_mod, type = "EFAtools",
                          factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2),
                    class = "efa_reliability_mixed_keying")
})

test_that("print output is stable", {
  local_reproducible_output()

  # single group, full coefficient matrix (general factor plus group factors)
  expect_snapshot(print(om_sl), transform = scrub_num)

  # single group, single factor
  expect_snapshot(print(structure(0.85, class = "OMEGA")), transform = scrub_num)
  expect_snapshot(print(structure(c(0.85, 2.10), class = "OMEGA")), transform = scrub_num)

  # multiple groups, full coefficient matrix (long header with ECV and PUC)
  om_mg6 <- structure(list(GroupA = unclass(om_sl), GroupB = unclass(om_sl)),
                      class = "OMEGA")
  expect_snapshot(print(om_mg6), transform = scrub_num)

  # multiple groups, three-column matrix (omega total, hierarchical, subscale)
  om_mg3 <- structure(list(GroupA = unclass(om_sl)[, 1:3],
                           GroupB = unclass(om_sl)[, 1:3]), class = "OMEGA")
  expect_snapshot(print(om_mg3), transform = scrub_num)

  # multiple groups, single factor
  expect_snapshot(print(structure(list(GroupA = 0.85, GroupB = 0.80), class = "OMEGA")),
                  transform = scrub_num)
  expect_snapshot(print(structure(list(GroupA = c(0.85, 2.10), GroupB = c(0.80, 1.90)),
                                  class = "OMEGA")), transform = scrub_num)
})

rm(efa_mod, sl_mod, om_sl, schmid_mod, om_schmid, om_man)
if (requireNamespace("lavaan", quietly = TRUE)) {
  rm(lav_mod, lav_fit, om_lav, lav_mod_ho, lav_fit_ho, om_lav_ho)
}
