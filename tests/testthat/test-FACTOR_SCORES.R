EFA_raw <- suppressMessages(EFA(DOSPERT_raw, n_factors = 10, type = "EFAtools",
                                method = "PAF", rotation = "oblimin"))
fac_scores_raw <- FACTOR_SCORES(DOSPERT_raw, f = EFA_raw, method = "Bartlett")

EFA_cor <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
               type = "EFAtools", method = "PAF", rotation = "oblimin")
fac_scores_cor <- suppressMessages(FACTOR_SCORES(test_models$baseline$cormat,
                                                 f = EFA_cor))

EFA_unrot <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                 type = "EFAtools", method = "PAF", rotation = "none")
fac_scores_unrot <- suppressMessages(FACTOR_SCORES(test_models$baseline$cormat,
                                                   f = EFA_unrot))

fac_scores_man <- suppressMessages(FACTOR_SCORES(test_models$baseline$cormat,
                                                 f = unclass(EFA_unrot$unrot_loadings)))


test_that("output is correct", {
  expect_s3_class(fac_scores_raw, "FACTOR_SCORES")
  expect_s3_class(fac_scores_cor, "FACTOR_SCORES")
  expect_s3_class(fac_scores_unrot, "FACTOR_SCORES")
  expect_s3_class(fac_scores_man, "FACTOR_SCORES")

  expect_named(fac_scores_raw, c("scores", "weights", "r.scores", "missing",
                                 "R2", "settings"))
  expect_named(fac_scores_cor, c("scores", "weights", "r.scores",
                                 "R2", "settings"))
  expect_named(fac_scores_unrot, c("scores", "weights", "r.scores", "R2",
                                   "settings"))
  expect_named(fac_scores_man, c("scores", "weights", "r.scores",
                                 "R2", "settings"))
})

test_that("factor scores have the expected shape and recover a known factor", {
  skip_on_cran()

  # raw-data scores: one row per observation, one column per factor (guards against
  # a transposed or weight-matrix mix-up in the returned object)
  expect_equal(dim(fac_scores_raw$scores), c(nrow(DOSPERT_raw), 10L))
  expect_false(anyNA(fac_scores_raw$scores))

  # On a strongly unidimensional scale the single estimated factor score must track
  # the factor it represents; the mean item response is a faithful proxy for that
  # latent factor, so the two correlate almost perfectly.
  efa_1f <- suppressMessages(EFA(GRiPS_raw, n_factors = 1, type = "EFAtools",
                                 method = "PAF"))
  fs_1f <- suppressMessages(FACTOR_SCORES(GRiPS_raw, f = efa_1f, method = "Bartlett"))
  expect_equal(dim(fs_1f$scores), c(nrow(GRiPS_raw), 1L))
  expect_gt(abs(stats::cor(fs_1f$scores[, 1], rowMeans(GRiPS_raw))), 0.95)
})

test_that("settings are returned correctly", {
  expect_named(fac_scores_raw$settings, "method")
  expect_named(fac_scores_cor$settings, "method")
  expect_named(fac_scores_unrot$settings, "method")
  expect_named(fac_scores_man$settings, "method")

  expect_equal(fac_scores_raw$settings$method, "Bartlett")
  expect_equal(fac_scores_cor$settings$method, "Thurstone")
  expect_equal(fac_scores_unrot$settings$method, "Thurstone")
  expect_equal(fac_scores_man$settings$method, "Thurstone")
})


test_that("warnings and errors are thrown correctly", {
  expect_error(FACTOR_SCORES(1:5), class = "efa_input_not_matrix")
  expect_error(FACTOR_SCORES(DOSPERT_raw, f = 1:5), class = "efa_scores_bad_f")
  expect_message(FACTOR_SCORES(test_models$baseline$cormat, f = EFA_cor), class = "efa_scores_needs_raw")
  expect_message(FACTOR_SCORES(test_models$baseline$cormat,
                               f = unclass(EFA_unrot$unrot_loadings)), class = "efa_scores_phi_null")

})


test_that("a non-Pearson EFA with rho = NULL on raw data warns", {
  skip_on_cran()

  # When rho is NULL, psych derives the score machinery from the Pearson cor(x):
  # for regression-based methods the weights/scores, and for every method the
  # score intercorrelations (r.scores), become inconsistent with a non-Pearson
  # fit; warn unless the matching matrix is supplied via rho.
  efa_sp <- suppressMessages(EFA(GRiPS_raw, n_factors = 1, type = "EFAtools",
                                 method = "PAF", cor_method = "spearman"))
  expect_warning(FACTOR_SCORES(GRiPS_raw, f = efa_sp, method = "Thurstone"),
                 class = "efa_scores_cor_method")

  # Supplying the matching correlation via rho silences the warning.
  rho_sp <- stats::cor(GRiPS_raw, method = "spearman")
  expect_no_warning(FACTOR_SCORES(GRiPS_raw, f = efa_sp, rho = rho_sp,
                                  method = "Thurstone"),
                    class = "efa_scores_cor_method")

  # Bartlett weights/scores are correlation-free, but psych still computes
  # r.scores from the Pearson cor(x), so the warning still fires for a
  # non-Pearson fit with rho = NULL.
  expect_warning(FACTOR_SCORES(GRiPS_raw, f = efa_sp, method = "Bartlett"),
                 class = "efa_scores_cor_method")
})

test_that("components scores work on a raw data.frame and Anderson guards single factor", {
  skip_on_cran()

  # psych's components path does `x %*% w` without coercing, so a data.frame raw
  # input must be coerced to a matrix by the wrapper; this previously errored.
  # EFA_raw (file top) is the matching 10-factor DOSPERT fit and is still in scope.
  fs_comp <- FACTOR_SCORES(DOSPERT_raw, f = EFA_raw, method = "components")
  expect_s3_class(fs_comp, "FACTOR_SCORES")
  expect_equal(dim(fs_comp$scores), c(nrow(DOSPERT_raw), 10L))
  expect_false(anyNA(fs_comp$scores))

  # Anderson-Rubin orthogonalisation is undefined for a single factor (and psych
  # aborts opaquely); the wrapper raises a classed error instead.
  efa_1f <- suppressMessages(EFA(GRiPS_raw, n_factors = 1, type = "EFAtools",
                                 method = "PAF"))
  expect_error(FACTOR_SCORES(GRiPS_raw, f = efa_1f, method = "Anderson"),
               class = "efa_scores_anderson_single")
})


rm(EFA_raw, fac_scores_raw, EFA_cor, fac_scores_cor, EFA_unrot, fac_scores_unrot,
   fac_scores_man)
