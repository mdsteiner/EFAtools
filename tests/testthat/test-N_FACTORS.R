nf_grips <- suppressMessages(suppressWarnings(N_FACTORS(GRiPS_raw)))

test_that("output class and dimensions are correct", {
  expect_is(nf_grips, "N_FACTORS")
  expect_is(nf_grips$outputs, "list")
  expect_is(nf_grips$settings, "list")
  expect_is(nf_grips$n_factors, "numeric")

  expect_named(nf_grips, c("outputs", "n_factors", "settings"))
  expect_named(nf_grips$outputs, c("bart_out", "kmo_out", "cd_out", "ekc_out",
                                   "hull_out", "kgc_out", "parallel_out",
                                   "scree_out", "smt_out"))
  expect_named(nf_grips$n_factors, c("nfac_CD", "nfac_EKC", "nfac_HULL_CAF",
                                     "nfac_HULL_CFI", "nfac_HULL_RMSEA",
                                     "nfac_KGC_PCA", "nfac_KGC_SMC",
                                     "nfac_KGC_EFA",
                                     "nfac_PA_PCA", "nfac_PA_SMC", "nfac_PA_EFA",
                                     "nfac_SMT_chi", "nfac_RMSEA", "nfac_AIC"))
  expect_named(nf_grips$settings, c("criteria", "suitability", "N", "use",
                                    "n_factors_max", "N_pop", "N_samples", "alpha",
                                    "cor_method", "max_iter_CD", "n_fac_theor",
                                    "method", "gof", "eigen_type_HULL",
                                    "eigen_type_other", "n_factors", "n_datasets",
                                    "percent", "decision_rule"))
})

x <- rnorm(100)
y <- rnorm(100)
z <- x + y

# Example from nearPD function from Matrix package
cor_nposdef <- matrix(c(1,     0.477, 0.644, 0.578, 0.651, 0.826,
                        0.477, 1,     0.516, 0.233, 0.682, 0.75,
                        0.644, 0.516, 1,     0.599, 0.581, 0.742,
                        0.478, 0.233, 0.599, 1,     0.741, 0.8,
                        0.651, 0.682, 0.581, 0.741, 1,     0.798,
                        0.826, 0.75,  0.742, 0.8,   0.798, 1),
                      nrow = 6, ncol = 6)

test_that("errors etc. are thrown correctly", {
  expect_error(N_FACTORS(1:10), " 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n")
  expect_warning(N_FACTORS(GRiPS_raw, N = 10), " 'N' was set and data entered. Taking N from data.\n")
  expect_error(N_FACTORS(cbind(x, y, z, z + 1, y + 1, x + 1)), " Correlation matrix is singular, no further analyses are performed\n")
  expect_warning(N_FACTORS(test_models$baseline$cormat, N = 500), " 'x' was a correlation matrix but CD needs raw data. Skipping CD.\n")
  expect_warning(N_FACTORS(cor_nposdef, N = 20), "Matrix was not positive definite, smoothing was done")
})

rm(nf_grips, x, y, z, cor_nposdef)
