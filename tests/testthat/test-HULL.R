hull_cor_paf <- suppressMessages(HULL(test_models$baseline$cormat, N = 500))
hull_cor_ml <- HULL(test_models$baseline$cormat, N = 500, method = "ML")
hull_cor_uls <- HULL(test_models$baseline$cormat, N = 500, method = "ULS")
hull_cor_uls_CFI <- HULL(test_models$baseline$cormat, N = 500, method = "ULS",
                         gof = "CFI")
hull_cor_ml_nf <- suppressWarnings(HULL(test_models$baseline$cormat, N = 500,
                                        method = "ML", n_fac_theor = 12))

hull_PCA <- HULL(test_models$baseline$cormat, N = 500, method = "ML",
                 eigen_type = "PCA")
hull_EFA <- HULL(test_models$baseline$cormat, N = 500, method = "ML",
                 eigen_type = "EFA")

hull_raw_paf <- suppressMessages(suppressWarnings(HULL(GRiPS_raw)))
hull_raw_ml <- suppressMessages(suppressWarnings(HULL(GRiPS_raw, method = "ML")))
hull_raw_uls <- suppressMessages(suppressWarnings(HULL(GRiPS_raw, method = "ULS")))
hull_raw_uls_CFI <- suppressMessages(suppressWarnings(HULL(GRiPS_raw,
                                                           method = "ULS",
                                                           gof = "CFI")))
hull_raw_ml_nf <- suppressMessages(suppressWarnings(HULL(GRiPS_raw, N = 500,
                                                         method = "ML",
                                                         n_fac_theor = 7)))


test_that("output class and dimensions are correct", {
  expect_is(hull_cor_paf, "HULL")
  expect_output(str(hull_cor_paf), "List of 8")
  expect_is(hull_cor_ml, "HULL")
  expect_output(str(hull_cor_ml), "List of 8")
  expect_is(hull_cor_uls, "HULL")
  expect_output(str(hull_cor_uls), "List of 8")
  expect_is(hull_cor_uls_CFI, "HULL")
  expect_output(str(hull_cor_uls_CFI), "List of 8")
  expect_is(hull_PCA, "HULL")
  expect_output(str(hull_PCA), "List of 8")
  expect_is(hull_EFA, "HULL")
  expect_output(str(hull_EFA), "List of 8")

  expect_is(hull_raw_paf, "HULL")
  expect_output(str(hull_raw_paf), "List of 8")
  expect_is(hull_raw_ml, "HULL")
  expect_output(str(hull_raw_ml), "List of 8")
  expect_is(hull_raw_uls, "HULL")
  expect_output(str(hull_raw_uls), "List of 8")
  expect_is(hull_raw_uls_CFI, "HULL")
  expect_output(str(hull_raw_uls_CFI), "List of 8")
  expect_is(hull_PCA, "HULL")
  expect_output(str(hull_PCA), "List of 8")
  expect_is(hull_EFA, "HULL")
  expect_output(str(hull_EFA), "List of 8")

  expect_output(str(hull_cor_paf$settings), "List of 7")
  expect_output(str(hull_cor_ml$settings), "List of 7")
  expect_output(str(hull_cor_uls$settings), "List of 7")
  expect_output(str(hull_cor_uls_CFI$settings), "List of 7")
  expect_output(str(hull_PCA$settings), "List of 7")
  expect_output(str(hull_EFA$settings), "List of 7")

  expect_output(str(hull_raw_paf$settings), "List of 7")
  expect_output(str(hull_raw_ml$settings), "List of 7")
  expect_output(str(hull_raw_uls$settings), "List of 7")
  expect_output(str(hull_raw_uls_CFI$settings), "List of 7")
  expect_output(str(hull_PCA$settings), "List of 7")
  expect_output(str(hull_EFA$settings), "List of 7")
})

test_that("n_fac_max is correctly specified", {
  expect_lte(hull_cor_paf$n_fac_max,
             .det_max_factors(ncol(test_models$baseline$cormat)))
  expect_lte(hull_cor_ml$n_fac_max,
             .det_max_factors(ncol(test_models$baseline$cormat)))
  expect_lte(hull_cor_uls$n_fac_max,
             .det_max_factors(ncol(test_models$baseline$cormat)))
  expect_lte(hull_cor_uls_CFI$n_fac_max,
             .det_max_factors(ncol(test_models$baseline$cormat)))
  expect_lte(hull_PCA$n_fac_max,
             .det_max_factors(ncol(test_models$baseline$cormat)))
  expect_lte(hull_EFA$n_fac_max,
             .det_max_factors(ncol(test_models$baseline$cormat)))

  expect_lte(hull_raw_paf$n_fac_max,
             .det_max_factors(ncol(GRiPS_raw)))
  expect_lte(hull_raw_ml$n_fac_max,
             .det_max_factors(ncol(GRiPS_raw)))
  expect_lte(hull_raw_uls$n_fac_max,
             .det_max_factors(ncol(GRiPS_raw)))
  expect_lte(hull_raw_uls_CFI$n_fac_max,
             .det_max_factors(ncol(GRiPS_raw)))

  expect_equal(hull_cor_ml_nf$n_fac_max,
             .det_max_factors(ncol(test_models$baseline$cormat)))
  expect_equal(hull_raw_ml_nf$n_fac_max,
               .det_max_factors(ncol(GRiPS_raw)))
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

  expect_equal(hull_raw_paf$n_fac_CAF, 1)
  expect_equal(hull_raw_paf$n_fac_CFI, NA)
  expect_equal(hull_raw_paf$n_fac_RMSEA, NA)

  expect_equal(hull_raw_uls$n_fac_CAF, 1)
  expect_equal(hull_raw_uls$n_fac_CFI, 1)
  expect_equal(hull_raw_uls$n_fac_RMSEA, 1)

  expect_equal(hull_raw_uls_CFI$n_fac_CAF, NA)
  expect_equal(hull_raw_uls_CFI$n_fac_CFI, 1)
  expect_equal(hull_raw_uls_CFI$n_fac_RMSEA, NA)
})


# Create singular correlation matrix for tests
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z, rnorm(10), rnorm(10), rnorm(10)), ncol = 6)
cor_sing <- stats::cor(dat_sing)

# Example from nearPD function from Matrix package
cor_nposdef <- matrix(c(1,     0.477, 0.644, 0.578, 0.651, 0.826,
                        0.477, 1,     0.516, 0.233, 0.682, 0.75,
                        0.644, 0.516, 1,     0.599, 0.581, 0.742,
                        0.478, 0.233, 0.599, 1,     0.741, 0.8,
                        0.651, 0.682, 0.581, 0.741, 1,     0.798,
                        0.826, 0.75,  0.742, 0.8,   0.798, 1),
                      nrow = 6, ncol = 6)

test_that("errors etc are thrown correctly", {
  expect_error(HULL(1:5), " 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n")
  expect_message(suppressWarnings(HULL(GRiPS_raw)), " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n")
  expect_warning(suppressMessages(HULL(GRiPS_raw)), " Suggested maximum number of factors was 2 but must be at least 3 for hull method to work. Setting it to 3.\n")
  expect_warning(suppressMessages(HULL(IDS2_R, N = 20)), " Less than three solutions located on the hull have been identified when using CAF as goodness of fit index. Proceeding by taking the value with the maximum CAF as heuristic. You may want to consider additional indices or methods as robustness check.\n")
  expect_warning(suppressMessages(HULL(GRiPS_raw, N = 20)), " 'N' was set and data entered. Taking N from data.\n")
  expect_error(HULL(test_models$baseline$cormat), ' "N" is not specified but is needed for computation of some of the fit indices.\n')
  expect_error(HULL(test_models$baseline$cormat, method = "ML"), ' "N" is not specified but is needed for computation of some of the fit indices.\n')
  expect_error(HULL(test_models$baseline$cormat, method = "ULS"), ' "N" is not specified but is needed for computation of some of the fit indices.\n')

  expect_error(HULL(dat_sing, method = "ML"), ' Correlation matrix is singular, the HULL method cannot be exectued.\n')
  expect_error(HULL(cor_sing, N = 20, method = "ML"), ' Correlation matrix is singular, the HULL method cannot be exectued.\n')

  expect_error(HULL(matrix(rnorm(50), ncol = 5)), "Data has fewer than 6 indicators. Hull method needs at least 6.\n")

  expect_message(suppressWarnings(HULL(GRiPS_raw)), 'Only CAF can be used as gof if method "PAF" is used. Setting gof to "CAF"\n')

  expect_warning(HULL(test_models$baseline$cormat, n_fac_theor = 13, N = 500), ' Setting maximum number of factors to 12 to ensure overidentified models.\n')
  expect_warning(HULL(cor_nposdef, N = 20, method = "ML"), "Matrix was not positive definite, smoothing was done")
  expect_warning(HULL(cor_nposdef, N = 20, method = "ML", n_fac_theor = 1), "Matrix was not positive definite, smoothing was done")
})

rm(hull_cor_paf, hull_cor_ml, hull_cor_uls, hull_cor_uls_CFI, hull_raw_paf,
   hull_raw_ml, hull_raw_uls, hull_raw_uls_CFI, hull_raw_ml_nf, hull_cor_ml_nf,
   hull_PCA, hull_EFA, x, y, z, dat_sing, cor_sing, cor_nposdef)
