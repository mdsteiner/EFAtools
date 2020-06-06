hull_cor_paf <- HULL(test_models$baseline$cormat, N = 500)
hull_cor_ml <- HULL(test_models$baseline$cormat, N = 500, method = "ML")
hull_cor_uls <- HULL(test_models$baseline$cormat, N = 500, method = "ULS")
hull_cor_uls_CFI <- HULL(test_models$baseline$cormat, N = 500, method = "ULS",
                         gof = "CFI")
hull_cor_ml_nf <- HULL(test_models$baseline$cormat, N = 500, method = "ML",
                       n_fac_theor = 12)

hull_raw_paf <- HULL(GRiPS_raw)
hull_raw_ml <- HULL(GRiPS_raw, method = "ML")
hull_raw_uls <- HULL(GRiPS_raw, method = "ULS")
hull_raw_uls_CFI <- HULL(GRiPS_raw, method = "ULS", gof = "CFI")
hull_raw_ml_nf <- HULL(GRiPS_raw, N = 500, method = "ML",
                       n_fac_theor = 7)



test_that("output class and dimensions are correct", {
  expect_is(hull_cor_paf, "HULL")
  expect_output(str(hull_cor_paf), "List of 8")
  expect_is(hull_cor_ml, "HULL")
  expect_output(str(hull_cor_ml), "List of 8")
  expect_is(hull_cor_uls, "HULL")
  expect_output(str(hull_cor_uls), "List of 8")
  expect_is(hull_cor_uls_CFI, "HULL")
  expect_output(str(hull_cor_uls_CFI), "List of 8")

  expect_is(hull_raw_paf, "HULL")
  expect_output(str(hull_raw_paf), "List of 8")
  expect_is(hull_raw_ml, "HULL")
  expect_output(str(hull_raw_ml), "List of 8")
  expect_is(hull_raw_uls, "HULL")
  expect_output(str(hull_raw_uls), "List of 8")
  expect_is(hull_raw_uls_CFI, "HULL")
  expect_output(str(hull_raw_uls_CFI), "List of 8")
})

test_that("n_fac_max is correctly specified", {
  expect_lte(hull_cor_paf$n_fac_max,
             floor(ncol(test_models$baseline$cormat) / 2))
  expect_lte(hull_cor_ml$n_fac_max,
             floor(ncol(test_models$baseline$cormat) / 2))
  expect_lte(hull_cor_uls$n_fac_max,
             floor(ncol(test_models$baseline$cormat) / 2))
  expect_lte(hull_cor_uls_CFI$n_fac_max,
             floor(ncol(test_models$baseline$cormat) / 2))

  expect_lte(hull_raw_paf$n_fac_max,
             floor(ncol(GRiPS_raw) / 2))
  expect_lte(hull_raw_ml$n_fac_max,
             floor(ncol(GRiPS_raw) / 2))
  expect_lte(hull_raw_uls$n_fac_max,
             floor(ncol(GRiPS_raw) / 2))
  expect_lte(hull_raw_uls_CFI$n_fac_max,
             floor(ncol(GRiPS_raw) / 2))

  expect_equal(hull_cor_ml_nf$n_fac_max,
             floor(ncol(test_models$baseline$cormat) / 2))
  expect_equal(hull_raw_ml_nf$n_fac_max,
               floor(ncol(GRiPS_raw) / 2))
})

test_that("solution matrices are correctly returned", {
  expect_is(hull_cor_paf$solutions_CAF, "matrix")
  expect_equal(hull_cor_paf$solutions_CFI, NA)
  expect_equal(hull_cor_paf$solutions_RMSEA, NA)

  expect_is(hull_cor_uls$solutions_CAF, "matrix")
  expect_is(hull_cor_uls$solutions_CFI, "matrix")
  expect_is(hull_cor_uls$solutions_RMSEA, "matrix")

  expect_is(hull_cor_uls_CFI$solutions_CFI, "matrix")
  expect_equal(hull_cor_uls_CFI$solutions_CAF, NA)
  expect_equal(hull_cor_uls_CFI$solutions_RMSEA, NA)

  expect_is(hull_raw_paf$solutions_CAF, "matrix")
  expect_equal(hull_raw_paf$solutions_CFI, NA)
  expect_equal(hull_raw_paf$solutions_RMSEA, NA)

  expect_is(hull_raw_uls$solutions_CAF, "matrix")
  expect_is(hull_raw_uls$solutions_CFI, "matrix")
  expect_is(hull_raw_uls$solutions_RMSEA, "matrix")

  expect_is(hull_raw_uls_CFI$solutions_CFI, "matrix")
  expect_equal(hull_raw_uls_CFI$solutions_CAF, NA)
  expect_equal(hull_raw_uls_CFI$solutions_RMSEA, NA)
})

test_that("n_factors are correctly returned", {
  expect_equal(hull_cor_paf$n_fac_CAF, 3)
  expect_equal(hull_cor_paf$n_fac_CFI, NA)
  expect_equal(hull_cor_paf$n_fac_RMSEA, NA)

  expect_equal(hull_cor_uls$n_fac_CAF, 3)
  expect_equal(hull_cor_uls$n_fac_CFI, 3)
  expect_equal(hull_cor_uls$n_fac_RMSEA, 3)

  expect_equal(hull_cor_uls_CFI$n_fac_CAF, NA)
  expect_equal(hull_cor_uls_CFI$n_fac_CFI, 3)
  expect_equal(hull_cor_uls_CFI$n_fac_RMSEA, NA)

  expect_equal(hull_raw_paf$n_fac_CAF, "matrix")
  expect_equal(hull_raw_paf$n_fac_CFI, NA)
  expect_equal(hull_raw_paf$n_fac_RMSEA, NA)

  expect_equal(hull_raw_uls$n_fac_CAF, "matrix")
  expect_equal(hull_raw_uls$n_fac_CFI, "matrix")
  expect_equal(hull_raw_uls$n_fac_RMSEA, "matrix")

  expect_equal(hull_raw_uls_CFI$n_fac_CAF, NA)
  expect_equal(hull_raw_uls_CFI$n_fac_CFI, "matrix")
  expect_equal(hull_raw_uls_CFI$n_fac_RMSEA, NA)
})


FURTHER HULL CHECKS -> N_FACTORS, WARNINGS AND ERRORS, SETTINGS
WARNING IF N_FAC MAX IS TOO LARGE

rm(hull_cor_paf, hull_cor_ml, hull_cor_uls, hull_cor_uls_CFI, hull_cor_paf,
   hull_cor_ml, hull_cor_uls, hull_cor_uls_CFI)
