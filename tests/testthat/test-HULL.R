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

  expect_output(str(hull_cor_paf$settings), "List of 7")
  expect_output(str(hull_cor_ml$settings), "List of 7")
  expect_output(str(hull_cor_uls$settings), "List of 7")
  expect_output(str(hull_cor_uls_CFI$settings), "List of 7")

  expect_output(str(hull_raw_paf$settings), "List of 7")
  expect_output(str(hull_raw_ml$settings), "List of 7")
  expect_output(str(hull_raw_uls$settings), "List of 7")
  expect_output(str(hull_raw_uls_CFI$settings), "List of 7")
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


# Create singular correlation matrix for tests
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z, rnorm(10), rnorm(10), rnorm(10)), ncol = 6)
cor_sing <- stats::cor(dat_sing)

test_that("errors etc are thrown correctly", {
  expect_error(HULL(test_models$baseline$cormat), ' "N" is not specified but is needed for computation of some of the fit indices.')
  expect_error(HULL(test_models$baseline$cormat, method = "ML"), ' "N" is not specified but is needed for computation of some of the fit indices.')
  expect_error(HULL(test_models$baseline$cormat, method = "ULS"), ' "N" is not specified but is needed for computation of some of the fit indices.')

  expect_error(HULL(dat_sing, method = "ML"), ' Correlation matrix is singular, the HULL method cannot be exectued.')
  expect_error(HULL(cor_sing, N = 10, method = "ULS"), ' Correlation matrix is singular, the HULL method cannot be exectued.')

  PROCEED HERE


  expect_warning(PARALLEL(GRiPS_raw, N = 20, eigen_type = "PCA"), '"N" is not specified but is needed for computation of some of the fit indices.')
  expect_error(PARALLEL(test_models$baseline$cormat, eigen_type = "PCA"), '"N" was not set and could not be taken from data. Please specify N and try again.')
  expect_warning(PARALLEL(test_models$baseline$cormat, N = 500,
                          eigen_type = "PCA", decision_rule = "Crawford",
                          percent = 80), "decision_rule == 'Crawford' is specified, but 95 percentile was not used. Using Means instead. To use 'Crawford', make sure to specify percent = 95.")
  expect_error(PARALLEL(dat_sing), " Correlation matrix is singular, parallel analysis is not possible.")
  expect_error(PARALLEL(cor_sing, N = 10), " Correlation matrix is singular, parallel analysis is not possible.")
})

FURTHER HULL CHECKS -> N_FACTORS, WARNINGS AND ERRORS
WARNING IF N_FAC MAX IS TOO LARGE

ERRORS
'"N" is not specified but is needed for computation of some of the fit indices.'
'Correlation matrix is singular, the HULL method cannot be exectued.'
'Data has fewer than 6 indicators. Hull method needs at least 6.'

WARNINGS
' "n_fac_theor" was larger than number of variables / 2. Setting maximum number of factors to number of variables / 2.'
" Suggested maximum number of factors was 2 but must be at least 3 for hull method to work. Setting it to 3."
" Less than three solutions located on the hull have been identified when using", gof_t, "as goodness of fit index. Proceeding by taking the value with the maximum",
gof_t, "as heuristic. You may want to consider additional indices or methods as robustness check."



rm(hull_cor_paf, hull_cor_ml, hull_cor_uls, hull_cor_uls_CFI, hull_cor_paf,
   hull_cor_ml, hull_cor_uls, hull_cor_uls_CFI)
