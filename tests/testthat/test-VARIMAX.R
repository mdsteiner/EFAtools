unrot <- EFA(test_models$baseline$cormat, 3, N = 500)
vari <- .VARIMAX(unrot, type = "EFAtools")
vari_psych <- .VARIMAX(unrot, type = "psych")
vari_spss <- .VARIMAX(unrot, type = "SPSS")

unrot_1 <- EFA(test_models$baseline$cormat, 1, N = 500)
vari_1 <- suppressWarnings(.VARIMAX(unrot_1, type = "EFAtools"))

test_that("output class and dimensions are correct", {
  expect_is(vari, "list")
  expect_is(vari_1, "list")
  expect_named(vari, c("rot_loadings", "rotmat", "vars_accounted_rot", "settings"))
  expect_named(vari_1, c("rot_loadings", "rotmat", "vars_accounted_rot", "settings"))

  expect_is(vari$rot_loadings, "LOADINGS")
  expect_is(vari$rotmat, "matrix")
  expect_is(vari$vars_accounted_rot, "matrix")
  expect_is(vari$settings, "list")

  expect_is(vari_1$rot_loadings, "LOADINGS")
  expect_equal(vari_1$rotmat, NA)
  expect_equal(vari_1$vars_accounted_rot, NA)
  expect_is(vari_1$settings, "list")
})

test_that("settings are returned correctly", {
  expect_named(vari$settings, c("kaiser", "precision", "order_type"))
  expect_named(vari_psych$settings, c("kaiser", "precision", "order_type"))
  expect_named(vari_spss$settings, c("kaiser", "precision", "order_type"))
  expect_named(vari_1$settings, c("kaiser", "precision", "order_type"))

  expect_equal(vari$settings$kaiser, TRUE)
  expect_equal(vari_psych$settings$kaiser, TRUE)
  expect_equal(vari_spss$settings$kaiser, TRUE)
  expect_equal(vari_1$settings$kaiser, TRUE)

  expect_equal(vari$settings$precision, 1e-05)
  expect_equal(vari_psych$settings$precision, 1e-05)
  expect_equal(vari_spss$settings$precision, 1e-10)
  expect_equal(vari_1$settings$precision, 1e-05)

  expect_equal(vari$settings$order_type, "eigen")
  expect_equal(vari_psych$settings$order_type, "eigen")
  expect_equal(vari_spss$settings$order_type, "ss_factors")
  expect_equal(vari_1$settings$order_type, "eigen")

})

test_that("errors etc. are thrown correctly", {

  expect_error(.VARIMAX(unrot, type = "none"), ' One of "precision", or "order_type" was NULL and no valid "type" was specified. Either use one of "EFAtools", "psych", or "SPSS" for type, or specify all other arguments\n')

  expect_warning(.VARIMAX(unrot, type = "EFAtools", kaiser = FALSE), " Type and kaiser is specified. kaiser is used with value ' FALSE '. Results may differ from the specified type\n")
  expect_warning(.VARIMAX(unrot, type = "EFAtools", precision = 1e-9), " Type and precision is specified. precision is used with value ' 1e-09 '. Results may differ from the specified type\n")
  expect_warning(.VARIMAX(unrot, type = "EFAtools", order_type = "ss_factors"), " Type and order_type is specified. order_type is used with value ' ss_factors '. Results may differ from the specified type\n")

  expect_warning(.VARIMAX(unrot, type = "psych", kaiser = FALSE), " Type and kaiser is specified. kaiser is used with value ' FALSE '. Results may differ from the specified type\n")
  expect_warning(.VARIMAX(unrot, type = "psych", precision = 1e-9), " Type and precision is specified. precision is used with value ' 1e-09 '. Results may differ from the specified type\n")
  expect_warning(.VARIMAX(unrot, type = "psych", order_type = "ss_factors"), " Type and order_type is specified. order_type is used with value ' ss_factors '. Results may differ from the specified type\n")

  expect_warning(.VARIMAX(unrot, type = "SPSS", kaiser = FALSE), " Type and kaiser is specified. kaiser is used with value ' FALSE '. Results may differ from the specified type\n")
  expect_warning(.VARIMAX(unrot, type = "SPSS", precision = 1e-9), " Type and precision is specified. precision is used with value ' 1e-09 '. Results may differ from the specified type\n")
  expect_warning(.VARIMAX(unrot, type = "SPSS", order_type = "eigen"), " Type and order_type is specified. order_type is used with value ' eigen '. Results may differ from the specified type\n")

  expect_warning(.VARIMAX(unrot_1, type = "EFAtools"), " Cannot rotate single factor. Unrotated loadings returned.\n")
})

rm(unrot, vari, unrot_1, vari_1, vari_psych, vari_spss)

