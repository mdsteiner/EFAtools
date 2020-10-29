unrot <- EFA(test_models$baseline$cormat, 3, N = 500)
equa <- .ROTATE_ORTH(unrot, rotation = "equamax", type = "EFAtools")

unrot_1 <- EFA(test_models$baseline$cormat, 1, N = 500)
equa_1 <- suppressWarnings(.ROTATE_ORTH(unrot_1, rotation = "equamax",
                                        type = "EFAtools"))

quarti <- .ROTATE_ORTH(unrot, rotation = "quartimax", type = "psych")
bentT <- .ROTATE_ORTH(unrot, rotation = "bentlerT", type = "none",
                      order_type = "eigen")
geoT <- .ROTATE_ORTH(unrot, rotation = "geominT", type = "SPSS")
bifacT <- .ROTATE_ORTH(unrot, rotation = "bifactorT", type = "EFAtools")

test_that("output class and dimensions are correct", {
  expect_is(equa$rot_loadings, "LOADINGS")
  expect_is(equa_1$rot_loadings, "LOADINGS") # The unrotated loadings here!
  expect_is(quarti$rot_loadings, "LOADINGS")
  expect_is(bentT$rot_loadings, "LOADINGS")
  expect_is(geoT$rot_loadings, "LOADINGS")
  expect_is(bifacT$rot_loadings, "LOADINGS")

  expect_output(str(equa), "List of 4")
  expect_output(str(equa_1), "List of 4")
  expect_output(str(quarti), "List of 4")
  expect_output(str(bentT), "List of 4")
  expect_output(str(geoT), "List of 4")
  expect_output(str(bifacT), "List of 4")

  expect_is(equa$rotmat, "matrix")
  expect_is(equa$vars_accounted_rot, "matrix")

  expect_is(quarti$rotmat, "matrix")
  expect_is(quarti$vars_accounted_rot, "matrix")

  expect_is(bentT$rotmat, "matrix")
  expect_is(bentT$vars_accounted_rot, "matrix")

  expect_is(geoT$rotmat, "matrix")
  expect_is(geoT$vars_accounted_rot, "matrix")

  expect_is(bifacT$rotmat, "matrix")
  expect_is(bifacT$vars_accounted_rot, "matrix")

  expect_equal(equa_1$rotmat, NA)
  expect_equal(equa_1$vars_accounted_rot, NA)
})

test_that("settings are returned correctly", {
  expect_named(equa$settings, c("normalize", "precision", "order_type"))
  expect_named(equa_1$settings,c("normalize", "precision", "order_type"))
  expect_named(quarti$settings, c("normalize", "precision", "order_type"))
  expect_named(bentT$settings, c("normalize", "precision", "order_type"))
  expect_named(geoT$settings, c("normalize", "precision", "order_type"))
  expect_named(bifacT$settings, c("normalize", "precision", "order_type"))

  expect_equal(equa$settings$normalize, TRUE)
  expect_equal(equa_1$settings$normalize, TRUE)
  expect_equal(quarti$settings$normalize, TRUE)
  expect_equal(bentT$settings$normalize, TRUE)
  expect_equal(geoT$settings$normalize, TRUE)
  expect_equal(bifacT$settings$normalize, TRUE)

  expect_equal(equa$settings$precision, 1e-5)
  expect_equal(equa_1$settings$precision, 1e-5)
  expect_equal(quarti$settings$precision, 1e-5)
  expect_equal(bentT$settings$precision, 1e-5)
  expect_equal(geoT$settings$precision, 1e-5)
  expect_equal(bifacT$settings$precision, 1e-5)

  expect_equal(equa$settings$order_type, "eigen")
  expect_equal(equa_1$settings$order_type, "eigen")
  expect_equal(quarti$settings$order_type, "eigen")
  expect_equal(bentT$settings$order_type, "eigen")
  expect_equal(geoT$settings$order_type, "ss_factors")
  expect_equal(bifacT$settings$order_type, "eigen")
})

test_that("errors etc. are thrown correctly", {
  expect_error(.ROTATE_ORTH(unrot, rotation = "equamax", type = "none"), ' "order_type" was NA and no valid "type" was specified. Either use one of "EFAtools", "psych", or "SPSS" for type, or specify the "order_type" argument\n')

  expect_warning(.ROTATE_ORTH(unrot, rotation = "equamax", type = "EFAtools",
                              normalize = FALSE), " Type and normalize is specified. normalize is used with value ' FALSE '. Results may differ from the specified type\n")
  expect_warning(.ROTATE_ORTH(unrot, rotation = "equamax", type = "EFAtools",
                              order_type = "ss_factors"), " Type and order_type is specified. order_type is used with value ' ss_factors '. Results may differ from the specified type\n")

  expect_warning(.ROTATE_ORTH(unrot, rotation = "equamax", type = "psych",
                              normalize = FALSE), " Type and normalize is specified. normalize is used with value ' FALSE '. Results may differ from the specified type\n")
  expect_warning(.ROTATE_ORTH(unrot, rotation = "equamax", type = "psych",
                              order_type = "ss_factors"), " Type and order_type is specified. order_type is used with value ' ss_factors '. Results may differ from the specified type\n")

  expect_warning(.ROTATE_ORTH(unrot, rotation = "equamax", type = "SPSS",
                              normalize = FALSE), " Type and normalize is specified. normalize is used with value ' FALSE '. Results may differ from the specified type\n")
  expect_warning(.ROTATE_ORTH(unrot, rotation = "equamax", type = "SPSS",
                              order_type = "ss_factors"), " Type and order_type is specified. order_type is used with value ' ss_factors '. Results may differ from the specified type\n")

  expect_warning(.ROTATE_ORTH(unrot_1, rotation = "equamax", type = "EFAtools"), " Cannot rotate single factor. Unrotated loadings returned.\n")
})

rm(unrot, equa, unrot_1, equa_1, quarti, bentT, geoT, bifacT)
