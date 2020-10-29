unrot <- EFA(test_models$baseline$cormat, 3, N = 500)
prom <- .PROMAX(unrot, type = "EFAtools")
prom_psych <- .PROMAX(unrot, type = "psych")
prom_spss <- .PROMAX(unrot, type = "SPSS")

unrot_1 <- EFA(test_models$baseline$cormat, 1, N = 500)
prom_1 <- suppressWarnings(.PROMAX(unrot_1, type = "EFAtools"))

test_that("output class and dimensions are correct", {
  expect_is(prom, "list")
  expect_is(prom_1, "list")
  expect_named(prom, c("rot_loadings", "Phi", "Structure", "rotmat",
                       "vars_accounted_rot", "settings"))
  expect_named(prom_1, c("rot_loadings", "Phi", "Structure", "rotmat",
                       "vars_accounted_rot", "settings"))

  expect_is(prom$rot_loadings, "LOADINGS")
  expect_is(prom$Phi, "matrix")
  expect_is(prom$Structure, "matrix")
  expect_is(prom$rotmat, "matrix")
  expect_is(prom$vars_accounted_rot, "matrix")
  expect_is(prom$settings, "list")

  expect_is(prom_1$rot_loadings, "LOADINGS")
  expect_equal(prom_1$Phi, NA)
  expect_equal(prom_1$Structure, NA)
  expect_equal(prom_1$rotmat, NA)
  expect_equal(prom_1$vars_accounted_rot, NA)
  expect_is(prom_1$settings, "list")
})

test_that("settings are returned correctly", {
  expect_named(prom$settings, c("normalize", "P_type", "precision",
                                "order_type", "varimax_type", "k"))
  expect_named(prom_psych$settings, c("normalize", "P_type", "precision",
                                      "order_type", "varimax_type", "k"))
  expect_named(prom_spss$settings, c("normalize", "P_type", "precision",
                                     "order_type", "varimax_type", "k"))
  expect_named(prom_1$settings, c("normalize", "P_type", "precision",
                                  "order_type", "varimax_type", "k"))

  expect_equal(prom$settings$normalize, TRUE)
  expect_equal(prom_psych$settings$normalize, TRUE)
  expect_equal(prom_spss$settings$normalize, TRUE)
  expect_equal(prom_1$settings$normalize, TRUE)

  expect_equal(prom$settings$P_type, "norm")
  expect_equal(prom_psych$settings$P_type, "unnorm")
  expect_equal(prom_spss$settings$P_type, "norm")
  expect_equal(prom_1$settings$P_type, "norm")

  expect_equal(prom$settings$precision, 1e-05)
  expect_equal(prom_psych$settings$precision, 1e-05)
  expect_equal(prom_spss$settings$precision, 1e-05)
  expect_equal(prom_1$settings$precision, 1e-05)

  expect_equal(prom$settings$order_type, "eigen")
  expect_equal(prom_psych$settings$order_type, "eigen")
  expect_equal(prom_spss$settings$order_type, "ss_factors")
  expect_equal(prom_1$settings$order_type, "eigen")

  expect_equal(prom$settings$varimax_type, "svd")
  expect_equal(prom_psych$settings$varimax_type, "svd")
  expect_equal(prom_spss$settings$varimax_type, "kaiser")
  expect_equal(prom_1$settings$varimax_type, "svd")

  expect_equal(prom$settings$k, 4)
  expect_equal(prom_psych$settings$k, 4)
  expect_equal(prom_spss$settings$k, 4)
  expect_equal(prom_1$settings$k, 4)

})

test_that("errors etc. are thrown correctly", {

  expect_error(.PROMAX(unrot, type = "none"), ' One of "P_type", "order_type", "varimax_type", or "k" was NA and no valid "type" was specified. Either use one of "EFAtools", "psych", or "SPSS" for type, or specify all other arguments\n')

  expect_warning(.PROMAX(unrot, type = "EFAtools", normalize = FALSE), " Type and normalize is specified. normalize is used with value ' FALSE '. Results may differ from the specified type\n")
  expect_warning(.PROMAX(unrot, type = "EFAtools", P_type = "norm"), " Type and P_type is specified. P_type is used with value ' norm '. Results may differ from the specified type\n")
  expect_warning(.PROMAX(unrot, type = "EFAtools", order_type = "ss_factors"), " Type and order_type is specified. order_type is used with value ' ss_factors '. Results may differ from the specified type\n")
  expect_warning(.PROMAX(unrot, type = "EFAtools", k = 2), " Type and k is specified. k is used with value ' 2 '. Results may differ from the specified type\n")

  expect_warning(.PROMAX(unrot, type = "psych", normalize = FALSE), " Type and normalize is specified. normalize is used with value ' FALSE '. Results may differ from the specified type.\n")
  expect_warning(.PROMAX(unrot, type = "psych", P_type = "norm"), " Type and P_type is specified. P_type is used with value ' norm '. Results may differ from the specified type\n")
  expect_warning(.PROMAX(unrot, type = "psych", order_type = "ss_factors"), " Type and order_type is specified. order_type is used with value ' ss_factors '. Results may differ from the specified type\n")
  expect_warning(.PROMAX(unrot, type = "psych", k = 2), " Type and k is specified. k is used with value ' 2 '. Results may differ from the specified type\n")

  expect_warning(.PROMAX(unrot, type = "SPSS", normalize = FALSE), " Type and normalize is specified. normalize is used with value ' FALSE '. Results may differ from the specified type.\n")
  expect_warning(.PROMAX(unrot, type = "SPSS", P_type = "unnorm"), " Type and P_type is specified. P_type is used with value ' unnorm '. Results may differ from the specified type\n")
  expect_warning(.PROMAX(unrot, type = "SPSS", order_type = "eigen"), " Type and order_type is specified. order_type is used with value ' eigen '. Results may differ from the specified type\n")
  expect_warning(.PROMAX(unrot, type = "SPSS", k = 2), " Type and k is specified. k is used with value ' 2 '. Results may differ from the specified type\n")
  expect_warning(.PROMAX(unrot, type = "SPSS", varimax_type = "svd"), " Type and varimax_type is specified. varimax_type is used with value ' svd '. Results may differ from the specified type\n")
  expect_warning(.PROMAX(unrot_1, type = "EFAtools"), " Cannot rotate single factor. Unrotated loadings returned.\n")
})

rm(unrot, prom, unrot_1, prom_1, prom_psych, prom_spss)
