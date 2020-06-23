test_that(".numformat works", {
  expect_equal(.numformat(0.23), " .23")
  expect_equal(.numformat(0.2345, digits = 3), " .234")
  expect_equal(.numformat(0.2345, digits = 3, print_zero = TRUE), " 0.234")
})


efa_temp <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500)
efa_pro <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                rotation = "promax")
test_that(".compute_vars works", {
  expect_is(.compute_vars(efa_temp$unrot_loadings,
                             efa_temp$unrot_loadings), "matrix")
  expect_is(.compute_vars(efa_pro$rot_loadings,
                             efa_pro$unrot_loadings,
                             efa_pro$Phi), "matrix")
})

x_base <- population_models$loadings$baseline
y_base <- x_base[, c(3,2,1)]

test_that(".factor_congruence works", {
  expect_is(.factor_congruence(x_base, y_base), "matrix")
  expect_equal(sum(.factor_congruence(x_base, y_base)), 3)
})


efa_ml <- suppressWarnings(EFA(cbind(rnorm(100), rnorm(100), rnorm(100), rnorm(100),
                    rnorm(100), rnorm(100)), 3, method = "ML"))
efa_uls <- suppressWarnings(EFA(cbind(rnorm(100), rnorm(100), rnorm(100), rnorm(100),
                    rnorm(100), rnorm(100)), 3, method = "ULS"))
efa_paf <- suppressWarnings(EFA(cbind(rnorm(100), rnorm(100), rnorm(100), rnorm(100),
                    rnorm(100), rnorm(100)), 3, method = "PAF"))

gof_ml <- .gof(efa_ml$unrot_loadings, efa_ml$orig_R, efa_ml$settings$N,
               "ML", efa_ml$fit_indices$Fm)
gof_uls <- .gof(efa_uls$unrot_loadings, efa_uls$orig_R, efa_uls$settings$N,
               "ULS", efa_uls$fit_indices$Fm)
gof_paf <- .gof(efa_paf$unrot_loadings, efa_paf$orig_R, efa_paf$settings$N,
               "PAF", NA)
m <- 6 # n variables
q <- 3 # n factors
test_that(".gof works", {
  expect_is(gof_ml, "list")
  expect_named(gof_ml,
               c("chi", "df", "CAF", "CFI", "RMSEA", "RMSEA_LB", "RMSEA_UB",
                 "AIC", "BIC", "Fm", "chi_null", "df_null"))
  expect_equal(gof_ml$CFI, 1)
  expect_equal(gof_ml$RMSEA, 0)
  expect_equal(gof_ml$CAF, .5, tolerance = .1)
  expect_equal(gof_ml$df, ((m - q)**2 - (m + q)) / 2)

  expect_is(gof_uls, "list")
  expect_named(gof_uls,
               c("chi", "df", "CAF", "CFI", "RMSEA", "RMSEA_LB", "RMSEA_UB",
                 "AIC", "BIC", "Fm", "chi_null", "df_null"))
  expect_equal(gof_uls$CFI, 1)
  expect_equal(gof_uls$RMSEA, 0)
  expect_equal(gof_uls$CAF, .5, tolerance = .1)
  expect_equal(gof_uls$df, ((m - q)**2 - (m + q)) / 2)

  expect_is(gof_paf, "list")
  expect_named(gof_paf,
               c("chi", "df", "CAF", "CFI", "RMSEA", "RMSEA_LB", "RMSEA_UB",
                 "AIC", "BIC", "Fm", "chi_null", "df_null"))
  expect_equal(gof_paf$CFI, NA)
  expect_equal(gof_paf$RMSEA, NA)
  expect_equal(gof_paf$CAF, .5, tolerance = .1)
  expect_equal(gof_paf$df, ((m - q)**2 - (m + q)) / 2)
})




test_that(".is_cormat works", {
  expect_equal(.is_cormat(cor(cbind(rnorm(100), rnorm(100)))), TRUE)
  expect_equal(.is_cormat(cbind(rnorm(100), rnorm(100))), FALSE)
  expect_equal(.is_cormat(cbind(rnorm(2), rnorm(2))), FALSE)
  expect_equal(.is_cormat(cbind(c(1, NA), rnorm(2))), FALSE)
  expect_equal(.is_cormat(matrix(c(1, .1, .3, 1), ncol = 2)), TRUE)
  expect_error(.is_cormat(matrix(c(1, NA, .3, 1), ncol = 2)), ' "x" is likely a correlation matrix but contains missing values. Please check the entered data.\n')
})

q_p <- .det_max_factors(8) + 1
test_that(".det_max_factors works", {
  expect_is(.det_max_factors(8), "numeric")
  expect_lte(((8 - q_p)**2 - (8 + q_p)) / 2, 0)
  expect_equal(.det_max_factors(0), 0)
  expect_equal(.det_max_factors(1), 0)
  expect_equal(.det_max_factors(2), 0)
  expect_equal(.det_max_factors(3), 0)
  expect_gt(.det_max_factors(4), 0)
})

rm(efa_pro, efa_temp, x_base, y_base, efa_ml, efa_uls, efa_paf, gof_ml, gof_uls,
   gof_paf, m, q, q_p)
