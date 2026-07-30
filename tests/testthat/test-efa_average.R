# efa_average runs its grid of EFAs through future.apply with future.seed = TRUE,
# so the random starts of the oblique rotation engines draw from the session RNG.
# Seed once here so the fixtures below - and the print snapshot built from them -
# are reproducible regardless of what ran earlier in the test run.
#
# The vector-valued default tuning settings expand to ~144 inner EFA fits per
# fixture, so the un-collapsed `efa_def` / `efa_all_none` / `efa_all_md`
# fixtures dominate the file (~25s combined). Skipped by default; opt in with
# `Sys.setenv(EFATOOLS_TEST_SLOW = "true")`. See helper-slow.R.
if (is_slow_test()) {
set.seed(42)

efa_def <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                       show_progress = FALSE)
efa_ml <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                      estimator = "ML", show_progress = FALSE)
efa_uls <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                      estimator = "ULS", show_progress = FALSE)

# Scalar overrides on the six default-vector tuning axes (init_comm,
# criterion_type, varimax_type, k_promax, P_type, start_method) for the
# fixtures whose downstream tests don't depend on the multi-value defaults.
# `efa_def` (above) still pins the vector defaults end-to-end; `efa_all_none`
# and `efa_all_md` (below) keep their full-grid expansion because they feed
# print-snapshot tests. Each of these scalar fixtures drops from ~1728 to ~12
# inner EFA fits.
.scalar_axes <- list(
  init_comm     = "smc",
  criterion_type = "sum",
  varimax_type  = "svd",
  k_promax      = 4,
  P_type        = "norm",
  start_method  = "psych"
)
efa_all <- do.call(efa_average, c(
  list(test_models$baseline$cormat, n_factors = 3, N = 500,
       estimator = c("PAF", "ML", "ULS"),
       type = c("none", "EFAtools", "psych", "SPSS"),
       salience_threshold = .2, show_progress = FALSE),
  .scalar_axes))
efa_all_oblq <- do.call(efa_average, c(
  list(test_models$baseline$cormat, n_factors = 3, N = 500,
       estimator = c("PAF", "ML", "ULS"),
       type = c("none", "EFAtools", "psych", "SPSS"),
       rotation = "oblique", show_progress = FALSE),
  .scalar_axes))
efa_all_orth <- do.call(efa_average, c(
  list(test_models$baseline$cormat, n_factors = 3, N = 500,
       estimator = c("PAF", "ML", "ULS"),
       type = c("none", "EFAtools", "psych", "SPSS"),
       rotation = "orthogonal", show_progress = FALSE),
  .scalar_axes))
efa_all_none <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                            estimator = c("PAF", "ML", "ULS"),
                            type = c("none", "EFAtools", "psych", "SPSS"),
                            rotation = "none", show_progress = FALSE)

efa_all_md <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                            estimator = c("PAF", "ML", "ULS"),
                            type = c("none", "EFAtools", "psych", "SPSS"),
                            rotation = "oblique", averaging = "median", show_progress = FALSE)
efa_all_tm <- do.call(efa_average, c(
  list(test_models$baseline$cormat, n_factors = 3, N = 500,
       estimator = c("PAF", "ML", "ULS"),
       type = c("none", "EFAtools", "psych", "SPSS"),
       rotation = "oblique", averaging = "mean",
       trim = .2, show_progress = FALSE),
  .scalar_axes))
efa_raw <- efa_average(GRiPS_raw, n_factors = 1, rotation = "none", show_progress = FALSE)
efa_raw_p <- efa_average(GRiPS_raw, n_factors = 2, rotation = "promax", show_progress = FALSE)
}  # is_slow_test()

test_that("efa_average runs end to end on a cheap grid (smoke test)", {
  # The rest of this file is gated behind EFATOOLS_TEST_SLOW, so without an
  # always-on test CRAN would never run efa_average() end to end. A two-row grid
  # (PAF and ML at the EFAtools preset, promax) fits in well under a second.
  local_reproducible_output()
  set.seed(42)
  res <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                     estimator = c("PAF", "ML"), type = "EFAtools",
                     rotation = "promax", start_method = "psych",
                     show_progress = FALSE)

  expect_s3_class(res, "efa_average")
  # PAF and ML at one preset and one rotation collapse to exactly two rows.
  expect_identical(nrow(res$implementations_grid), 2L)

  # Slot shapes: variables x factors loadings (carrying the LOADINGS class), a
  # square factor-correlation matrix from the oblique rotation, one communality
  # per variable, and the same-shaped indicator-to-factor correspondence matrix.
  m <- ncol(test_models$baseline$cormat)
  expect_s3_class(res$loadings$average, "LOADINGS")
  expect_identical(dim(res$loadings$average), c(m, 3L))
  expect_identical(dim(res$Phi$average), c(3L, 3L))
  expect_length(res$h2$average, m)
  expect_identical(dim(res$ind_fac_corres), c(m, 3L))

  # The averaging-rate sentence must never contain NaN (the empty-denominator guard).
  expect_no_match(paste(cli::ansi_strip(format(res)), collapse = " "), "NaN",
                  fixed = TRUE)
})

test_that("output class and dimensions are correct", {
  skip_if_not_slow()
  expect_s3_class(efa_def, "efa_average")
  expect_s3_class(efa_ml, "efa_average")
  expect_s3_class(efa_uls, "efa_average")
  expect_s3_class(efa_all, "efa_average")
  expect_s3_class(efa_all_oblq, "efa_average")
  expect_s3_class(efa_all_orth, "efa_average")
  expect_s3_class(efa_all_none, "efa_average")
  expect_s3_class(efa_all_md, "efa_average")
  expect_s3_class(efa_all_tm, "efa_average")
  expect_s3_class(efa_raw, "efa_average")
  expect_s3_class(efa_raw_p, "efa_average")

  # the old class string is kept alongside the new one, so `inherits(x, "EFA_AVERAGE")`
  # in existing code still resolves
  expect_identical(class(efa_def), c("efa_average", "EFA_AVERAGE"))

  expect_s3_class(efa_def$loadings$average, "LOADINGS")
  expect_s3_class(efa_ml$loadings$average, "LOADINGS")
  expect_s3_class(efa_uls$loadings$average, "LOADINGS")
  expect_s3_class(efa_all$loadings$average, "LOADINGS")
  expect_s3_class(efa_all_oblq$loadings$average, "LOADINGS")
  expect_s3_class(efa_all_orth$loadings$average, "LOADINGS")
  expect_s3_class(efa_all_none$loadings$average, "LOADINGS")
  expect_s3_class(efa_all_md$loadings$average, "LOADINGS")
  expect_s3_class(efa_all_tm$loadings$average, "LOADINGS")
  expect_s3_class(efa_raw$loadings$average, "LOADINGS")
  expect_s3_class(efa_raw_p$loadings$average, "LOADINGS")

  expect_named(efa_def, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_ml, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_uls, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_all, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_all_oblq, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_all_orth, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_all_none, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_all_md, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_all_tm, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_raw, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_raw_p, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
})

test_that("settings are returned correctly", {
  skip_if_not_slow()
  expect_named(efa_def$settings, c("estimator", "method", "rotation", "type", "n_factors", "N",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen", "varimax_type", "normalize",
                                   "k_promax", "k_simplimax", "P_type",
                                   "precision", "start_method", "use",
                                   "cor_method", "max_iter", "averaging",
                                   "trim", "salience_threshold"))
  expect_named(efa_ml$settings, c("estimator", "method", "rotation", "type", "n_factors", "N",
                                  "init_comm", "criterion", "criterion_type",
                                  "abs_eigen", "varimax_type", "normalize",
                                  "k_promax", "k_simplimax", "P_type",
                                  "precision", "start_method", "use",
                                  "cor_method", "max_iter", "averaging",
                                  "trim", "salience_threshold"))
  expect_named(efa_uls$settings, c("estimator", "method", "rotation", "type", "n_factors", "N",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen", "varimax_type", "normalize",
                                   "k_promax", "k_simplimax", "P_type",
                                   "precision", "start_method", "use",
                                   "cor_method", "max_iter", "averaging",
                                   "trim", "salience_threshold"))
  expect_named(efa_all$settings, c("estimator", "method", "rotation", "type", "n_factors", "N",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen", "varimax_type", "normalize",
                                   "k_promax", "k_simplimax", "P_type",
                                   "precision", "start_method", "use",
                                   "cor_method", "max_iter", "averaging",
                                   "trim", "salience_threshold"))
  expect_named(efa_all_oblq$settings, c("estimator", "method", "rotation", "type", "n_factors", "N",
                                        "init_comm", "criterion", "criterion_type",
                                        "abs_eigen", "varimax_type", "normalize",
                                        "k_promax", "k_simplimax", "P_type",
                                        "precision", "start_method", "use",
                                        "cor_method", "max_iter", "averaging",
                                        "trim", "salience_threshold"))
  expect_named(efa_all_orth$settings, c("estimator", "method", "rotation", "type", "n_factors", "N",
                                        "init_comm", "criterion", "criterion_type",
                                        "abs_eigen", "varimax_type", "normalize",
                                        "k_promax", "k_simplimax", "P_type",
                                        "precision", "start_method", "use",
                                        "cor_method", "max_iter", "averaging",
                                        "trim", "salience_threshold"))
  expect_named(efa_all_none$settings, c("estimator", "method", "rotation", "type", "n_factors", "N",
                                        "init_comm", "criterion", "criterion_type",
                                        "abs_eigen", "varimax_type", "normalize",
                                        "k_promax", "k_simplimax", "P_type",
                                        "precision", "start_method", "use",
                                        "cor_method", "max_iter", "averaging",
                                        "trim", "salience_threshold"))
  expect_named(efa_all_md$settings, c("estimator", "method", "rotation", "type", "n_factors", "N",
                                      "init_comm", "criterion", "criterion_type",
                                      "abs_eigen", "varimax_type", "normalize",
                                      "k_promax", "k_simplimax", "P_type",
                                      "precision", "start_method", "use",
                                      "cor_method", "max_iter", "averaging",
                                      "trim", "salience_threshold"))
  expect_named(efa_all_tm$settings, c("estimator", "method", "rotation", "type", "n_factors", "N",
                                      "init_comm", "criterion", "criterion_type",
                                      "abs_eigen", "varimax_type", "normalize",
                                      "k_promax", "k_simplimax", "P_type",
                                      "precision", "start_method", "use",
                                      "cor_method", "max_iter", "averaging",
                                      "trim", "salience_threshold"))
  expect_named(efa_raw$settings, c("estimator", "method", "rotation", "type", "n_factors", "N",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen", "varimax_type", "normalize",
                                   "k_promax", "k_simplimax", "P_type",
                                   "precision", "start_method", "use",
                                   "cor_method", "max_iter", "averaging",
                                   "trim", "salience_threshold"))
  expect_named(efa_raw_p$settings, c("estimator", "method", "rotation", "type", "n_factors", "N",
                                     "init_comm", "criterion", "criterion_type",
                                     "abs_eigen", "varimax_type", "normalize",
                                     "k_promax", "k_simplimax", "P_type",
                                     "precision", "start_method", "use",
                                     "cor_method", "max_iter", "averaging",
                                     "trim", "salience_threshold"))


  expect_equal(efa_def$settings$estimator, "PAF")
  expect_equal(efa_ml$settings$estimator, "ML")
  expect_equal(efa_uls$settings$estimator, "ULS")
  expect_equal(efa_all$settings$estimator, c("PAF", "ML", "ULS"))
  expect_equal(efa_all_oblq$settings$estimator, c("PAF", "ML", "ULS"))
  expect_equal(efa_all_orth$settings$estimator, c("PAF", "ML", "ULS"))
  expect_equal(efa_all_none$settings$estimator, c("PAF", "ML", "ULS"))
  expect_equal(efa_all_md$settings$estimator, c("PAF", "ML", "ULS"))
  expect_equal(efa_all_tm$settings$estimator, c("PAF", "ML", "ULS"))
  expect_equal(efa_raw$settings$estimator, "PAF")
  expect_equal(efa_raw_p$settings$estimator, "PAF")

  expect_equal(efa_def$settings$rotation, "promax")
  expect_equal(efa_ml$settings$rotation, "promax")
  expect_equal(efa_uls$settings$rotation, "promax")
  expect_equal(efa_all$settings$rotation, "promax")
  expect_equal(efa_all_oblq$settings$rotation, "oblique")
  expect_equal(efa_all_orth$settings$rotation, "orthogonal")
  expect_equal(efa_all_none$settings$rotation, "none")
  expect_equal(efa_all_md$settings$rotation, "oblique")
  expect_equal(efa_all_tm$settings$rotation, "oblique")
  expect_equal(efa_raw$settings$rotation, "none")
  expect_equal(efa_raw_p$settings$rotation, "promax")

  expect_equal(efa_def$settings$type, "none")
  expect_equal(efa_ml$settings$type, "none")
  expect_equal(efa_uls$settings$type, "none")
  expect_equal(efa_all$settings$type, c("none", "EFAtools", "psych", "SPSS"))
  expect_equal(efa_all_oblq$settings$type, c("none", "EFAtools", "psych", "SPSS"))
  expect_equal(efa_all_orth$settings$type, c("none", "EFAtools", "psych", "SPSS"))
  expect_equal(efa_all_none$settings$type, c("none", "EFAtools", "psych", "SPSS"))
  expect_equal(efa_all_md$settings$type, c("none", "EFAtools", "psych", "SPSS"))
  expect_equal(efa_all_tm$settings$type, c("none", "EFAtools", "psych", "SPSS"))
  expect_equal(efa_raw$settings$type, "none")
  expect_equal(efa_raw_p$settings$type, "none")

  expect_equal(efa_def$settings$n_factors, 3)
  expect_equal(efa_ml$settings$n_factors, 3)
  expect_equal(efa_uls$settings$n_factors, 3)
  expect_equal(efa_all$settings$n_factors, 3)
  expect_equal(efa_all_oblq$settings$n_factors, 3)
  expect_equal(efa_all_orth$settings$n_factors, 3)
  expect_equal(efa_all_none$settings$n_factors, 3)
  expect_equal(efa_all_md$settings$n_factors, 3)
  expect_equal(efa_all_tm$settings$n_factors, 3)
  expect_equal(efa_raw$settings$n_factors, 1)
  expect_equal(efa_raw_p$settings$n_factors, 2)

  expect_equal(efa_def$settings$N, 500)
  expect_equal(efa_ml$settings$N, 500)
  expect_equal(efa_uls$settings$N, 500)
  expect_equal(efa_all$settings$N, 500)
  expect_equal(efa_all_oblq$settings$N, 500)
  expect_equal(efa_all_orth$settings$N, 500)
  expect_equal(efa_all_none$settings$N, 500)
  expect_equal(efa_all_md$settings$N, 500)
  expect_equal(efa_all_tm$settings$N, 500)
  expect_equal(efa_raw$settings$N, 810)
  expect_equal(efa_raw_p$settings$N, 810)

  expect_equal(efa_def$settings$init_comm, c("smc", "mac", "unity"))
  expect_equal(efa_ml$settings$init_comm, c("smc", "mac", "unity"))
  expect_equal(efa_uls$settings$init_comm, c("smc", "mac", "unity"))
  expect_equal(efa_all$settings$init_comm, "smc")
  expect_equal(efa_all_oblq$settings$init_comm, "smc")
  expect_equal(efa_all_orth$settings$init_comm, "smc")
  expect_equal(efa_all_none$settings$init_comm, c("smc", "mac", "unity"))
  expect_equal(efa_all_md$settings$init_comm, c("smc", "mac", "unity"))
  expect_equal(efa_all_tm$settings$init_comm, "smc")
  expect_equal(efa_raw$settings$init_comm, c("smc", "mac", "unity"))
  expect_equal(efa_raw_p$settings$init_comm, c("smc", "mac", "unity"))

  expect_equal(efa_def$settings$criterion, 0.001)
  expect_equal(efa_ml$settings$criterion, 0.001)
  expect_equal(efa_uls$settings$criterion, 0.001)
  expect_equal(efa_all$settings$criterion, 0.001)
  expect_equal(efa_all_oblq$settings$criterion, 0.001)
  expect_equal(efa_all_orth$settings$criterion, 0.001)
  expect_equal(efa_all_none$settings$criterion, 0.001)
  expect_equal(efa_all_md$settings$criterion, 0.001)
  expect_equal(efa_all_tm$settings$criterion, 0.001)
  expect_equal(efa_raw$settings$criterion, 0.001)
  expect_equal(efa_raw_p$settings$criterion, 0.001)

  expect_equal(efa_def$settings$criterion_type, c("sum", "max_individual"))
  expect_equal(efa_ml$settings$criterion_type, c("sum", "max_individual"))
  expect_equal(efa_uls$settings$criterion_type, c("sum", "max_individual"))
  expect_equal(efa_all$settings$criterion_type, "sum")
  expect_equal(efa_all_oblq$settings$criterion_type, "sum")
  expect_equal(efa_all_orth$settings$criterion_type, "sum")
  expect_equal(efa_all_none$settings$criterion_type, c("sum", "max_individual"))
  expect_equal(efa_all_md$settings$criterion_type, c("sum", "max_individual"))
  expect_equal(efa_all_tm$settings$criterion_type, "sum")
  expect_equal(efa_raw$settings$criterion_type, c("sum", "max_individual"))
  expect_equal(efa_raw_p$settings$criterion_type, c("sum", "max_individual"))

  expect_equal(efa_def$settings$abs_eigen, TRUE)
  expect_equal(efa_ml$settings$abs_eigen, TRUE)
  expect_equal(efa_uls$settings$abs_eigen, TRUE)
  expect_equal(efa_all$settings$abs_eigen, TRUE)
  expect_equal(efa_all_oblq$settings$abs_eigen, TRUE)
  expect_equal(efa_all_orth$settings$abs_eigen, TRUE)
  expect_equal(efa_all_none$settings$abs_eigen, TRUE)
  expect_equal(efa_all_md$settings$abs_eigen, TRUE)
  expect_equal(efa_all_tm$settings$abs_eigen, TRUE)
  expect_equal(efa_raw$settings$abs_eigen, TRUE)
  expect_equal(efa_raw_p$settings$abs_eigen, TRUE)

  expect_equal(efa_def$settings$abs_eigen, TRUE)
  expect_equal(efa_ml$settings$abs_eigen, TRUE)
  expect_equal(efa_uls$settings$abs_eigen, TRUE)
  expect_equal(efa_all$settings$abs_eigen, TRUE)
  expect_equal(efa_all_oblq$settings$abs_eigen, TRUE)
  expect_equal(efa_all_orth$settings$abs_eigen, TRUE)
  expect_equal(efa_all_none$settings$abs_eigen, TRUE)
  expect_equal(efa_all_md$settings$abs_eigen, TRUE)
  expect_equal(efa_all_tm$settings$abs_eigen, TRUE)
  expect_equal(efa_raw$settings$abs_eigen, TRUE)
  expect_equal(efa_raw_p$settings$abs_eigen, TRUE)

  expect_equal(efa_def$settings$varimax_type, c("svd", "kaiser"))
  expect_equal(efa_ml$settings$varimax_type, c("svd", "kaiser"))
  expect_equal(efa_uls$settings$varimax_type, c("svd", "kaiser"))
  expect_equal(efa_all$settings$varimax_type, "svd")
  expect_equal(efa_all_oblq$settings$varimax_type, "svd")
  expect_equal(efa_all_orth$settings$varimax_type, "svd")
  expect_equal(efa_all_none$settings$varimax_type, c("svd", "kaiser"))
  expect_equal(efa_all_md$settings$varimax_type, c("svd", "kaiser"))
  expect_equal(efa_all_tm$settings$varimax_type, "svd")
  expect_equal(efa_raw$settings$varimax_type, c("svd", "kaiser"))
  expect_equal(efa_raw_p$settings$varimax_type, c("svd", "kaiser"))

  expect_equal(efa_def$settings$normalize, TRUE)
  expect_equal(efa_ml$settings$normalize, TRUE)
  expect_equal(efa_uls$settings$normalize, TRUE)
  expect_equal(efa_all$settings$normalize, TRUE)
  expect_equal(efa_all_oblq$settings$normalize, TRUE)
  expect_equal(efa_all_orth$settings$normalize, TRUE)
  expect_equal(efa_all_none$settings$normalize, TRUE)
  expect_equal(efa_all_md$settings$normalize, TRUE)
  expect_equal(efa_all_tm$settings$normalize, TRUE)
  expect_equal(efa_raw$settings$normalize, TRUE)
  expect_equal(efa_raw_p$settings$normalize, TRUE)

  expect_equal(efa_def$settings$k_promax, 2:4)
  expect_equal(efa_ml$settings$k_promax, 2:4)
  expect_equal(efa_uls$settings$k_promax, 2:4)
  expect_equal(efa_all$settings$k_promax, 4)
  expect_equal(efa_all_oblq$settings$k_promax, 4)
  expect_equal(efa_all_orth$settings$k_promax, 4)
  expect_equal(efa_all_none$settings$k_promax, 2:4)
  expect_equal(efa_all_md$settings$k_promax, 2:4)
  expect_equal(efa_all_tm$settings$k_promax, 4)
  expect_equal(efa_raw$settings$k_promax, 2:4)
  expect_equal(efa_raw_p$settings$k_promax, 2:4)

  expect_equal(efa_def$settings$k_simplimax, 18)
  expect_equal(efa_ml$settings$k_simplimax, 18)
  expect_equal(efa_uls$settings$k_simplimax, 18)
  expect_equal(efa_all$settings$k_simplimax, 18)
  expect_equal(efa_all_oblq$settings$k_simplimax, 18)
  expect_equal(efa_all_orth$settings$k_simplimax, 18)
  expect_equal(efa_all_none$settings$k_simplimax, 18)
  expect_equal(efa_all_md$settings$k_simplimax, 18)
  expect_equal(efa_all_tm$settings$k_simplimax, 18)
  expect_equal(efa_raw$settings$k_simplimax, 8)
  expect_equal(efa_raw_p$settings$k_simplimax, 8)

  expect_equal(efa_def$settings$P_type, c("norm", "unnorm"))
  expect_equal(efa_ml$settings$P_type, c("norm", "unnorm"))
  expect_equal(efa_uls$settings$P_type, c("norm", "unnorm"))
  expect_equal(efa_all$settings$P_type, "norm")
  expect_equal(efa_all_oblq$settings$P_type, "norm")
  expect_equal(efa_all_orth$settings$P_type, "norm")
  expect_equal(efa_all_none$settings$P_type, c("norm", "unnorm"))
  expect_equal(efa_all_md$settings$P_type, c("norm", "unnorm"))
  expect_equal(efa_all_tm$settings$P_type, "norm")
  expect_equal(efa_raw$settings$P_type, c("norm", "unnorm"))
  expect_equal(efa_raw_p$settings$P_type, c("norm", "unnorm"))

  expect_equal(efa_def$settings$precision, 1e-5)
  expect_equal(efa_ml$settings$precision, 1e-5)
  expect_equal(efa_uls$settings$precision, 1e-5)
  expect_equal(efa_all$settings$precision, 1e-5)
  expect_equal(efa_all_oblq$settings$precision, 1e-5)
  expect_equal(efa_all_orth$settings$precision, 1e-5)
  expect_equal(efa_all_none$settings$precision, 1e-5)
  expect_equal(efa_all_md$settings$precision, 1e-5)
  expect_equal(efa_all_tm$settings$precision, 1e-5)
  expect_equal(efa_raw$settings$precision, 1e-5)
  expect_equal(efa_raw_p$settings$precision, 1e-5)

  expect_equal(efa_def$settings$start_method, c("psych", "factanal"))
  expect_equal(efa_ml$settings$start_method, c("psych", "factanal"))
  expect_equal(efa_uls$settings$start_method, c("psych", "factanal"))
  expect_equal(efa_all$settings$start_method, "psych")
  expect_equal(efa_all_oblq$settings$start_method, "psych")
  expect_equal(efa_all_orth$settings$start_method, "psych")
  expect_equal(efa_all_none$settings$start_method, c("psych", "factanal"))
  expect_equal(efa_all_md$settings$start_method, c("psych", "factanal"))
  expect_equal(efa_all_tm$settings$start_method, "psych")
  expect_equal(efa_raw$settings$start_method, c("psych", "factanal"))
  expect_equal(efa_raw_p$settings$start_method, c("psych", "factanal"))

  expect_equal(efa_def$settings$use, "pairwise.complete.obs")
  expect_equal(efa_ml$settings$use, "pairwise.complete.obs")
  expect_equal(efa_uls$settings$use, "pairwise.complete.obs")
  expect_equal(efa_all$settings$use, "pairwise.complete.obs")
  expect_equal(efa_all_oblq$settings$use, "pairwise.complete.obs")
  expect_equal(efa_all_orth$settings$use, "pairwise.complete.obs")
  expect_equal(efa_all_none$settings$use, "pairwise.complete.obs")
  expect_equal(efa_all_md$settings$use, "pairwise.complete.obs")
  expect_equal(efa_all_tm$settings$use, "pairwise.complete.obs")
  expect_equal(efa_raw$settings$use, "pairwise.complete.obs")
  expect_equal(efa_raw_p$settings$use, "pairwise.complete.obs")

  expect_equal(efa_def$settings$cor_method, "pearson")
  expect_equal(efa_ml$settings$cor_method, "pearson")
  expect_equal(efa_uls$settings$cor_method, "pearson")
  expect_equal(efa_all$settings$cor_method, "pearson")
  expect_equal(efa_all_oblq$settings$cor_method, "pearson")
  expect_equal(efa_all_orth$settings$cor_method, "pearson")
  expect_equal(efa_all_none$settings$cor_method, "pearson")
  expect_equal(efa_all_md$settings$cor_method, "pearson")
  expect_equal(efa_all_tm$settings$cor_method, "pearson")
  expect_equal(efa_raw$settings$cor_method, "pearson")
  expect_equal(efa_raw_p$settings$cor_method, "pearson")

  expect_equal(efa_def$settings$max_iter, 10000)
  expect_equal(efa_ml$settings$max_iter, 10000)
  expect_equal(efa_uls$settings$max_iter, 10000)
  expect_equal(efa_all$settings$max_iter, 10000)
  expect_equal(efa_all_oblq$settings$max_iter, 10000)
  expect_equal(efa_all_orth$settings$max_iter, 10000)
  expect_equal(efa_all_none$settings$max_iter, 10000)
  expect_equal(efa_all_md$settings$max_iter, 10000)
  expect_equal(efa_all_tm$settings$max_iter, 10000)
  expect_equal(efa_raw$settings$max_iter, 10000)
  expect_equal(efa_raw_p$settings$max_iter, 10000)

  expect_equal(efa_def$settings$averaging, "mean")
  expect_equal(efa_ml$settings$averaging, "mean")
  expect_equal(efa_uls$settings$averaging, "mean")
  expect_equal(efa_all$settings$averaging, "mean")
  expect_equal(efa_all_oblq$settings$averaging, "mean")
  expect_equal(efa_all_orth$settings$averaging, "mean")
  expect_equal(efa_all_none$settings$averaging, "mean")
  expect_equal(efa_all_md$settings$averaging, "median")
  expect_equal(efa_all_tm$settings$averaging, "mean")
  expect_equal(efa_raw$settings$averaging, "mean")
  expect_equal(efa_raw_p$settings$averaging, "mean")

  expect_equal(efa_def$settings$trim, 0)
  expect_equal(efa_ml$settings$trim, 0)
  expect_equal(efa_uls$settings$trim, 0)
  expect_equal(efa_all$settings$trim, 0)
  expect_equal(efa_all_oblq$settings$trim, 0)
  expect_equal(efa_all_orth$settings$trim, 0)
  expect_equal(efa_all_none$settings$trim, 0)
  expect_equal(efa_all_md$settings$trim, 0)
  expect_equal(efa_all_tm$settings$trim, 0.2)
  expect_equal(efa_raw$settings$trim, 0)
  expect_equal(efa_raw_p$settings$trim, 0)

  expect_equal(efa_def$settings$salience_threshold, 0.3)
  expect_equal(efa_ml$settings$salience_threshold, 0.3)
  expect_equal(efa_uls$settings$salience_threshold, 0.3)
  expect_equal(efa_all$settings$salience_threshold, 0.2)
  expect_equal(efa_all_oblq$settings$salience_threshold, 0.3)
  expect_equal(efa_all_orth$settings$salience_threshold, 0.3)
  expect_equal(efa_all_none$settings$salience_threshold, 0.3)
  expect_equal(efa_all_md$settings$salience_threshold, 0.3)
  expect_equal(efa_all_tm$settings$salience_threshold, 0.3)
  expect_equal(efa_raw$settings$salience_threshold, 0.3)
  expect_equal(efa_raw_p$settings$salience_threshold, 0.3)


})


# Create singular correlation matrix for tests
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z), ncol = 3)
cor_sing <- stats::cor(dat_sing)

cor_nposdef <- matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), ncol = 3)

test_that("errors are thrown correctly", {
  skip_if_not_slow()
  expect_error(efa_average(1:5, show_progress = FALSE), class = "efa_input_not_matrix")
  # Extracting two factors from the unidimensional GRiPS overshoots, and psych's
  # 50-iteration cap now stops that solution short of convergence, so it is
  # excluded with a warning; suppress that incidental warning to assert the target
  # condition (as the n_factors = 1 cases below already do).
  expect_message(suppressWarnings(
    efa_average(GRiPS_raw, n_factors = 2, estimator = "PAF", type = c("EFAtools", "psych"), show_progress = FALSE),
    classes = "efa_avg_excluded_solutions"),
                 class = "efa_cor_from_data")
  expect_warning(suppressWarnings(
    efa_average(GRiPS_raw, n_factors = 2, estimator = "PAF", type = c("EFAtools", "psych"),
                N = 20, show_progress = FALSE),
    classes = "efa_avg_excluded_solutions"),
                 class = "efa_n_from_data")
  expect_error(efa_average(dat_sing, n_factors = 1, show_progress = FALSE),
               class = "efa_cor_singular")
  expect_error(efa_average(cor_sing, N = 10, n_factors = 1, show_progress = FALSE),
               class = "efa_cor_singular")
  expect_error(efa_average(matrix(rnorm(30), ncol = 3), n_factors = 2, show_progress = FALSE),
               class = "efa_underidentified")
  expect_warning(suppressWarnings(efa_average(matrix(rnorm(30), ncol = 3), n_factors = 1,
                             estimator = "PAF", type = c("EFAtools", "psych"), show_progress = FALSE),
                             class = c("efa_avg_excluded_solutions")),
                 class = "efa_just_identified")
  expect_warning(suppressWarnings(efa_average(cor_nposdef, n_factors = 1, N = 10, estimator = "PAF",
                     type = c("EFAtools", "psych"), show_progress = FALSE), classes = c("efa_just_identified", "efa_avg_excluded_solutions")), class = "efa_cor_smoothed")
  expect_message(efa_average(GRiPS_raw, n_factors = 1, estimator = "PAF", type = c("EFAtools", "psych"), show_progress = FALSE),
                 class = "efa_avg_single_factor_rotation")
  expect_warning(efa_average(GRiPS_raw, n_factors = 1, estimator = "PAF", type = c("EFAtools"),
                             rotation = "none", show_progress = FALSE),
                 class = "efa_avg_single_combination")
})

test_that("the grid progress bar declares exactly the steps it signals", {
  # progressr is inert in a non-interactive session and every other test in this file
  # runs with show_progress = FALSE, so this is the only place the progress accounting
  # is exercised. Sending more updates than the progressor declares makes
  # with_progress() warn that it is "no longer listening to this progressor"; the
  # warning then also styles the returned object's auto-printed output as warning
  # output. Small grids are the exposed case: they report on every EFA, so the per-EFA
  # updates alone complete the progressor.
  withr::local_options(progressr.enable = TRUE)

  # 3 EFAs (one per estimator at the EFAtools preset), the grid from @examples. Grids
  # of ten or fewer report on every EFA, so stepsize is 1.
  expect_no_warning(
    efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                estimator = c("PAF", "ULS", "ML"), type = "EFAtools",
                start_method = "psych"))

  # 24 unrotated PAF EFAs, enough for stepsize to exceed 1 (round(24/10) = 2). That is
  # the branch in which the declared steps stop being a whole number of updates, so it
  # needs its own case; the cheap unrotated PAF grid keeps it fast.
  expect_no_warning(
    efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                estimator = "PAF", type = "none", rotation = "none",
                init_comm = c("smc", "mac", "unity"),
                criterion = c(1e-3, 1e-4),
                criterion_type = c("sum", "max_individual"),
                abs_eigen = c(TRUE, FALSE)))
})

test_that("the Model Fit block reports the effective n per index class", {
  # The print snapshot masks these counts (they follow the per-solution convergence and
  # Heywood outcomes, which move across BLAS implementations), so the wiring - which grid
  # column feeds which printed line - is pinned here instead.
  local_reproducible_output()
  set.seed(42)

  # PAF leaves the chi-square-based indices NA but still contributes CAF, RMSR, and SRMR,
  # so the two index classes are averaged over different numbers of solutions.
  mixed <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                       estimator = c("PAF", "ML"), type = "EFAtools",
                       rotation = "promax", start_method = "psych",
                       show_progress = FALSE)
  grid <- mixed$implementations_grid
  n_chisq <- sum(!is.na(grid$chisq))
  n_resid <- sum(!is.na(grid$caf))
  expect_gt(n_resid, n_chisq)

  out <- paste(cli::ansi_strip(format(mixed)), collapse = " ")
  expect_match(out, paste0("Chi-square-based indices averaged over ", n_chisq,
                           " of ", nrow(grid), " solution"), fixed = TRUE)
  expect_match(out, paste0("CAF, RMSR, and SRMR averaged over ", n_resid,
                           " of ", nrow(grid), " solution"), fixed = TRUE)

  # With PAF alone no chi-square-based index is reported at all, so that line is absent.
  paf <- suppressWarnings(
    efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                estimator = "PAF", type = c("EFAtools", "psych"),
                rotation = "promax", show_progress = FALSE))
  out_paf <- paste(cli::ansi_strip(format(paf)), collapse = " ")
  expect_no_match(out_paf, "Chi-square-based indices", fixed = TRUE)
  expect_match(out_paf, paste0("CAF, RMSR, and SRMR averaged over ",
                               sum(!is.na(paf$implementations_grid$caf)), " of "),
               fixed = TRUE)
})

test_that("the Model Fit block names the aggregation it reports", {
  # The sentence introducing the block must agree with the M / Md header below it.
  local_reproducible_output()
  set.seed(42)
  args <- list(test_models$baseline$cormat, n_factors = 3, N = 500,
               estimator = c("PAF", "ML"), type = "EFAtools", rotation = "promax",
               start_method = "psych", show_progress = FALSE)

  mn <- paste(cli::ansi_strip(format(do.call(efa_average, args))), collapse = " ")
  expect_match(mn, "the mean of the per-solution fit indices", fixed = TRUE)

  md <- paste(cli::ansi_strip(format(do.call(efa_average, c(args, averaging = "median")))),
              collapse = " ")
  expect_match(md, "the median of the per-solution fit indices", fixed = TRUE)
  expect_no_match(md, "the average of the per-solution fit indices", fixed = TRUE)
})

test_that("a trim that the median cannot use is flagged", {
  # trim only reaches mean(), so a non-zero trim under averaging = "median" changes
  # nothing; it must not be accepted silently. The message is raised during argument
  # checking, before any EFA is fitted, so this needs no grid.
  expect_message(
    suppressWarnings(efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                                 estimator = "PAF", type = "EFAtools", rotation = "none",
                                 averaging = "median", trim = 0.2,
                                 show_progress = FALSE)),
    class = "efa_avg_trim_ignored"
  )
  # ... and not raised when the trim is actually used.
  expect_no_message(
    suppressWarnings(efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                                 estimator = "PAF", type = "EFAtools", rotation = "none",
                                 averaging = "mean", trim = 0.2,
                                 show_progress = FALSE)),
    class = "efa_avg_trim_ignored"
  )
})

test_that("an all-failed averaging grid returns an empty (NA) result", {
  skip_if_not_slow()
  # When every solution fails (here all runs hit max_iter and do not converge),
  # the averaged result is NA rather than an error or an average over an empty set.
  res <- suppressWarnings(efa_average(test_models$baseline$cormat, n_factors = 3,
                                      N = 500, max_iter = 1, estimator = "PAF",
                                      type = "none", rotation = "none",
                                      show_progress = FALSE))
  expect_s3_class(res, "efa_average")
  expect_true(all((res$implementations_grid$converged != 0) %in% TRUE))
  expect_true(all(is.na(res$h2)))
  expect_true(all(is.na(res$loadings)))
})

test_that("an all-Heywood averaging grid returns an empty (NA) result", {
  skip_if_not_slow()
  # A 3-variable matrix whose single-factor solution implies a communality > 1, so
  # every fitted solution is a Heywood case and is excluded; with nothing left to
  # average, the communalities and loadings come back as NA rather than an error.
  m <- matrix(c(1, .8, .8,  .8, 1, .5,  .8, .5, 1), 3, 3)
  res <- suppressWarnings(efa_average(m, n_factors = 1, N = 200, estimator = "PAF",
                                      type = "none", rotation = "none",
                                      show_progress = FALSE))
  expect_s3_class(res, "efa_average")
  expect_true(all(res$implementations_grid$heywood))
  expect_true(all(is.na(res$h2)))
  expect_true(all(is.na(res$loadings)))
})

test_that("print output is stable", {
  skip_if_not_slow()
  skip_on_cran()
  local_reproducible_output()

  # default view (oblique -> Phi section, average + range)
  expect_snapshot(print(efa_def, plot = FALSE), transform = scrub_num_pct)

  # all statistics, no Phi (orthogonal-free rotation = "none")
  expect_snapshot(print(efa_all_none, stat = c("average", "sd", "min", "max"),
                        plot = FALSE), transform = scrub_num_pct)

  # median averaging (Md / Median labels)
  expect_snapshot(print(efa_all_md, plot = FALSE), transform = scrub_num_pct)
})

test_that("format.efa_average is the source of truth", {
  skip_on_cran()
  skip_if_not_slow()
  local_reproducible_output()

  # print() is exactly cat(format(x), sep = "\n") (with plot = FALSE, the default, so
  # nothing is drawn), so the two agree line for line.
  expect_identical(utils::capture.output(print(efa_def)), format(efa_def))
})

test_that("plot returns a ggplot", {
  skip_on_cran()
  skip_if_not_slow()
  expect_s3_class(plot(efa_def), "ggplot")
})

test_that("the average-loading plot is visually stable", {
  skip_if_not_installed("vdiffr")

  # A literal fixture: plot.efa_average() reads only the average/min/max loading
  # matrices and the two settings below, so hand-chosen values keep the baseline
  # deterministic (the fitted fixtures above are grids of many EFAs).
  L <- matrix(c( 0.82, -0.11,
                 0.45,  0.60,
                -0.30,  0.71,
                 0.04, -0.05),
              nrow = 4, byrow = TRUE,
              dimnames = list(paste0("V", 1:4), c("F1", "F2")))
  avg <- structure(
    list(loadings = list(average = L, min = L - 0.08, max = L + 0.12),
         settings = list(averaging = "mean", salience_threshold = 0.2)),
    class = "efa_average")

  # the axis title, the caption naming the four marks, and the marks themselves
  vdiffr::expect_doppelganger("efa_average loading plot", plot(avg))
})

test_that("a vector-valued precision is recycled across the grid", {
  skip_if_not_slow()
  # A vector precision is expanded into the grid; each EFA must receive its own
  # scalar value rather than the whole vector (which would fail every fit and
  # return an all-NA result).
  res <- suppressWarnings(
    efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                estimator = "ML", rotation = "promax", type = "none",
                precision = c(1e-5, 1e-3), show_progress = FALSE))
  expect_s3_class(res, "efa_average")
  expect_setequal(unique(res$implementations_grid$precision), c(1e-5, 1e-3))
  expect_false(any(res$implementations_grid$errors))
  expect_false(all(is.na(res$loadings$average)))
  expect_true(all(is.finite(res$loadings$average)))
})

test_that("problematic solutions are summarised in a single warning", {
  skip_if_not_slow()
  # All runs hit max_iter and are excluded, so one summary warning is raised
  # rather than one per model.
  expect_warning(
    efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                max_iter = 1, estimator = "PAF", type = "none", rotation = "none",
                show_progress = FALSE),
    class = "efa_avg_excluded_solutions"
  )
})

# The rates summary reads only the implementations_grid, so the denominator states below can
# be exercised without fitting anything.
rates_text <- function(grid) {
  paste(cli::ansi_strip(cli::cli_format_method(.efa_emit_average_rates(grid))),
        collapse = " ")
}

test_that("each rate names the denominator it is computed over", {
  local_reproducible_output()
  # Convergence is assessed for every solution that did not error, but Heywood cases and
  # admissibility only for those that then converged. Here: 2 errored, 2 did not converge,
  # and 6 converged (3 Heywood, 1 proper but not salient, 2 admissible). The Heywood and
  # admissibility denominators are therefore 6, not the 8 that did not error - over 8 the
  # Heywood rate would read 38% rather than 50%.
  grid <- data.frame(
    errors     = c(rep(TRUE, 2), rep(FALSE, 8)),
    converged  = c(rep(NA, 2), rep(1, 2), rep(0, 6)),
    heywood    = c(rep(NA, 4), rep(TRUE, 3), rep(FALSE, 3)),
    admissible = c(rep(NA, 4), rep(FALSE, 4), rep(TRUE, 2))
  )

  out <- rates_text(grid)
  expect_match(out, "The error rate is at 20%.", fixed = TRUE)
  expect_match(out, "Of the solutions that did not result in an error, 75% converged.",
               fixed = TRUE)
  expect_match(out,
               "Of the solutions that converged, 50% contained Heywood cases and 33% were admissible.",
               fixed = TRUE)
})

test_that("Heywood cases and admissibility are not reported when no solution converged", {
  local_reproducible_output()
  # Nothing reached the stage at which either is assessed, so both rates are undefined
  # rather than 0%.
  grid <- data.frame(errors = FALSE, converged = c(1, 1), heywood = NA, admissible = NA)

  out <- rates_text(grid)
  expect_no_match(out, "NaN", fixed = TRUE)
  expect_match(out, "Of the solutions that did not result in an error, 0% converged.",
               fixed = TRUE)
  expect_match(out,
               "No solution converged, so Heywood cases and admissibility could not be assessed.",
               fixed = TRUE)
})

test_that("only the error rate is reported when every EFA errored", {
  local_reproducible_output()
  # No model was fitted at all, so convergence is undefined too.
  grid <- data.frame(errors = TRUE, converged = c(NA, NA), heywood = NA, admissible = NA)

  out <- rates_text(grid)
  expect_no_match(out, "NaN", fixed = TRUE)
  expect_match(out, "The error rate is at 100%.", fixed = TRUE)
  expect_match(out,
               "No solution could be fitted, so convergence, Heywood cases, and admissibility could not be assessed.",
               fixed = TRUE)
})

test_that("printing a run in which no EFA converged reports no NaN rate", {
  skip_if_not_slow()
  local_reproducible_output()
  # Every run hits max_iter, so the Heywood and admissibility rates have an empty
  # denominator on the way through the public API.
  res <- suppressWarnings(
    efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                max_iter = 1, estimator = "PAF", type = "none", rotation = "none",
                show_progress = FALSE))

  expect_no_match(paste(cli::ansi_strip(format(res)), collapse = " "), "NaN", fixed = TRUE)
})

test_that("printing does not plot by default", {
  expect_false(formals(getS3method("print", "efa_average"))$plot)
})

test_that("admissibility is reported as an outcome, not a varied setting", {
  skip_on_cran()
  skip_if_not_slow()
  local_reproducible_output()
  obj <- efa_def
  # Force mixed admissibility so the column would surface in the varied-settings
  # list if it were (incorrectly) treated as a setting.
  obj$implementations_grid$admissible[1] <- !obj$implementations_grid$admissible[1]
  out <- paste(cli::ansi_strip(utils::capture.output(print(obj, plot = FALSE))),
               collapse = " ")
  # "admissible" appears only in the convergence summary, never in the settings list.
  expect_equal(unname(lengths(gregexpr("admissible", out, fixed = TRUE))), 1L)
})

test_that("settings varied equally across the grid are still listed", {
  local_reproducible_output()
  # A grid whose settings columns each hold the same number of distinct values must
  # still report every varied setting. Every model is a Heywood case, so the report
  # stops after the summary and no averaged solution is needed.
  obj <- structure(list(
    settings = list(estimator = c("PAF", "ML"), rotation = c("promax", "oblimin"),
                    averaging = "mean", trim = 0),
    implementations_grid = data.frame(
      estimator  = c("PAF", "ML"),
      rotation   = c("promax", "oblimin"),
      errors     = FALSE,
      converged  = 0,
      heywood    = TRUE,
      admissible = FALSE
    )), class = "efa_average")

  out <- paste(cli::ansi_strip(format(obj)), collapse = " ")
  expect_match(out, "settings: estimator and rotation")
})

test_that("averaged fit indices match the per-model grid means", {
  skip_if_not_slow()
  # efa_def uses the untrimmed mean, so each averaged index equals the column
  # mean over the included models (and the new residual indices are averaged too).
  fi <- efa_def$fit_indices
  grid <- efa_def$implementations_grid
  expect_equal(fi$average[fi$index == "caf"], mean(grid$caf, na.rm = TRUE))
  expect_equal(fi$average[fi$index == "srmr"], mean(grid$srmr, na.rm = TRUE))

  # df is constant across the averaged solutions, so its dispersion is zero: the
  # sd and range cells of the df row are 0 (average/min/max hold df itself).
  expect_equal(fi$sd[fi$index == "df"], 0)
  expect_equal(fi$range[fi$index == "df"], 0)
  expect_equal(fi$average[fi$index == "df"], fi$min[fi$index == "df"])
})

if (is_slow_test()) {
  rm(efa_def, efa_ml, efa_uls, efa_all, efa_all_oblq, efa_all_orth, efa_all_none,
     efa_all_md, efa_all_tm, efa_raw, efa_raw_p)
}



