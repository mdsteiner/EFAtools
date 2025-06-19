nf_grips <- suppressMessages(suppressWarnings(N_FACTORS(GRiPS_raw)))

test_that("output class and dimensions are correct", {
  expect_is(nf_grips, "N_FACTORS")
  expect_is(nf_grips$outputs, "list")
  expect_is(nf_grips$settings, "list")
  expect_is(nf_grips$n_factors, "numeric")

  expect_named(nf_grips, c("outputs", "n_factors", "settings"))
  expect_named(nf_grips$outputs, c("bart_out", "kmo_out", "cd_out", "ekc_out",
                                   "hull_out", "kgc_out", "parallel_out",
                                  "nest_out", "scree_out", "smt_out"))
  expect_named(nf_grips$n_factors, c("nfac_CD",
                                     "nfac_EKC_BvA2017", "nfac_EKC_AM2019",
                                     "nfac_HULL_CAF",
                                     "nfac_HULL_CFI", "nfac_HULL_RMSEA",
                                     "nfac_KGC_PCA", "nfac_KGC_SMC",
                                     "nfac_KGC_EFA",
                                     "nfac_PA_PCA", "nfac_PA_SMC", "nfac_PA_EFA",
                                     "nfac_NEST",
                                     "nfac_SMT_chi", "nfac_RMSEA", "nfac_AIC"))
  expect_named(nf_grips$settings, c("criteria", "suitability", "N", "use",
                                    "n_factors_max", "N_pop", "N_samples", "alpha",
                                    "cor_method", "max_iter_CD", "n_fac_theor",
                                    "method", "gof", "eigen_type_HULL",
                                    "eigen_type_other", "n_factors", "n_datasets",
                                    "percent", "decision_rule",
                                    "ekc_type", "n_datasets_nest",
                                    "alpha_nest"))
})

x <- rnorm(100)
y <- rnorm(100)
z <- x + y



burt <- matrix(c(1.00,  0.83,  0.81,  0.80,   0.71, 0.70, 0.54, 0.53,  0.59,  0.24, 0.13,
                 0.83,  1.00,  0.87,  0.62,   0.59, 0.44, 0.58, 0.44,  0.23,  0.45,  0.21,
                 0.81,  0.87,  1.00,  0.63,   0.37, 0.31, 0.30, 0.12,  0.33,  0.33,  0.36,
                 0.80,  0.62,  0.63,  1.00,   0.49, 0.54, 0.30, 0.28,  0.42,  0.29, -0.06,
                 0.71,  0.59,  0.37,  0.49,   1.00, 0.54, 0.34, 0.55,  0.40,  0.19, -0.10,
                 0.70,  0.44,  0.31,  0.54,   0.54, 1.00, 0.50, 0.51,  0.31,  0.11,  0.10,
                 0.54,  0.58,  0.30,  0.30,   0.34, 0.50, 1.00, 0.38,  0.29,  0.21,  0.08,
                 0.53,  0.44,  0.12,  0.28,   0.55, 0.51, 0.38, 1.00,  0.53,  0.10, -0.16,
                 0.59,  0.23,  0.33,  0.42,   0.40, 0.31, 0.29, 0.53,  1.00, -0.09, -0.10,
                 0.24,  0.45,  0.33,  0.29,   0.19, 0.11, 0.21, 0.10, -0.09,  1.00,  0.41,
                 0.13,  0.21,  0.36, -0.06,  -0.10, 0.10, 0.08, -0.16, -0.10, 0.41,  1.00),
               nrow = 11, ncol = 11)

test_that("errors etc. are thrown correctly", {
  expect_error(N_FACTORS(1:10), " 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n")
  expect_warning(N_FACTORS(GRiPS_raw, N = 10), " 'N' was set and data entered. Taking N from data.\n")
  expect_error(N_FACTORS(cbind(x, y, z, z + 1, y + 1, x + 1)), " Correlation matrix is singular, no further analyses are performed\n")
  expect_warning(N_FACTORS(test_models$baseline$cormat, N = 500), " 'x' was a correlation matrix but CD needs raw data. Skipping CD.\n")
  expect_warning(N_FACTORS(burt, N = 170, criteria = c("PARALLEL", "EKC")), "Matrix was not positive definite, smoothing was done")
})

rm(nf_grips, x, y, z, burt)
