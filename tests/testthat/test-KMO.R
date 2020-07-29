kmo_cor <- KMO(test_models$baseline$cormat)
kmo_raw <- KMO(GRiPS_raw)
dat_nonames <- test_models$baseline$cormat
colnames(dat_nonames) <- NULL
kmo_nona <- KMO(dat_nonames)

test_that("output class and dimensions are correct", {
  expect_is(kmo_cor, "KMO")
  expect_output(str(kmo_cor), "List of 3")
  expect_is(kmo_raw, "KMO")
  expect_output(str(kmo_raw), "List of 3")
  expect_is(kmo_nona, "KMO")
  expect_output(str(kmo_nona), "List of 3")
})

test_that("KMO values are correct", {
  expect_equal(kmo_cor$KMO, 0.916, tolerance = 1e-3)
  expect_equal(kmo_raw$KMO, 0.955, tolerance = 1e-3)
  expect_equal(kmo_nona$KMO, 0.916, tolerance = 1e-3)

  expect_length(kmo_cor$KMO_i, ncol(test_models$baseline$cormat))
  expect_length(kmo_raw$KMO_i, ncol(GRiPS_raw))
  expect_length(kmo_nona$KMO_i, ncol(dat_nonames))
})

test_that("settings are returned correctly", {
  expect_named(kmo_cor$settings, c("use", "cor_method"))
  expect_named(kmo_raw$settings, c("use", "cor_method"))
  expect_named(kmo_nona$settings, c("use", "cor_method"))

  expect_equal(kmo_cor$settings$use, "pairwise.complete.obs")
  expect_equal(kmo_raw$settings$use, "pairwise.complete.obs")
  expect_equal(kmo_nona$settings$use, "pairwise.complete.obs")

  expect_equal(kmo_cor$settings$cor_method, "pearson")
  expect_equal(kmo_raw$settings$cor_method, "pearson")
  expect_equal(kmo_nona$settings$cor_method, "pearson")
})

# Create singular correlation matrix for tests
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z), ncol = 3)
cor_sing <- stats::cor(dat_sing)

cor_nposdef <- matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), ncol = 3)

test_that("errors are thrown correctly", {
  expect_error(KMO(1:5), " 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n")
  expect_message(KMO(GRiPS_raw), " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n")
  expect_error(KMO(dat_sing), " Correlation matrix is singular, no further analyses are performed.\n")
  expect_error(KMO(cor_sing), " Correlation matrix is singular, no further analyses are performed.\n")
  expect_warning(KMO(cor_nposdef), "Matrix was not positive definite, smoothing was done")
})

rm(kmo_cor, kmo_raw, dat_nonames, kmo_nona, x, y, z, dat_sing, cor_sing, cor_nposdef)
