kgc_cor <- KGC(test_models$baseline$cormat)
kgc_cor_smc <- KGC(test_models$baseline$cormat, eigen_type = "SMC")
kgc_raw <- KGC(GRiPS_raw)

test_that("output class and dimensions are correct", {
  expect_is(kgc_cor, "KGC")
  expect_output(str(kgc_cor), "List of 7")
  expect_is(kgc_raw, "KGC")
  expect_output(str(kgc_raw), "List of 7")
})

test_that("found eigenvalues are correct", {
  expect_equal(sum(kgc_cor$eigen_PCA), ncol(test_models$baseline$cormat))
  expect_lt(sum(kgc_cor$eigen_SMC), ncol(test_models$baseline$cormat))
  expect_lt(sum(kgc_cor$eigen_EFA), ncol(test_models$baseline$cormat))

  expect_equal(sum(kgc_raw$eigen_PCA), ncol(GRiPS_raw))
  expect_lt(sum(kgc_raw$eigen_SMC), ncol(GRiPS_raw))
  expect_lt(sum(kgc_raw$eigen_EFA), ncol(GRiPS_raw))

  expect_equal(c(kgc_cor_smc$eigen_PCA, kgc_cor_smc$eigen_EFA), c(NA, NA))

})

test_that("identified number of factors is correct", {
  expect_equal(kgc_cor$n_fac_PCA, 3)
  expect_equal(kgc_cor$n_fac_SMC, 1)
  expect_equal(kgc_cor$n_fac_EFA, 1)

  expect_equal(kgc_raw$n_fac_PCA, 1)
  expect_equal(kgc_raw$n_fac_SMC, 1)
  expect_equal(kgc_raw$n_fac_EFA, 1)

  expect_equal(kgc_cor_smc$n_fac_SMC, 1)
  expect_equal(c(kgc_cor_smc$n_fac_PCA, kgc_cor_smc$n_fac_EFA), c(NA, NA))
})
