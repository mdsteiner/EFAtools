unrot <- EFA(test_models$baseline$cormat, 3, N = 500)
equa <- .ROTATE_OBLQ(unrot, rotation = "equamax", type = "EFAtools")

unrot_1 <- EFA(test_models$baseline$cormat, 1, N = 500)
equa_1 <- .ROTATE_OBLQ(unrot_1, rotation = "equamax", type = "EFAtools")

quarti <- .ROTATE_OBLQ(unrot, rotation = "quartimax", type = "psych")
bentT <- .ROTATE_OBLQ(unrot, rotation = "bentlerT", type = "none",
                      order_type = "eigen")
geoT <- .ROTATE_OBLQ(unrot, rotation = "geominT", type = "EFAtools")
bifacT <- .ROTATE_OBLQ(unrot, rotation = "bifactorT", type = "EFAtools")

# HIER STEHENGEBLIEBEN

test_that("output class and dimensions are correct", {
  expect_is(equa$rot_loadings, "LOADINGS")
  expect_is(equa_1$rot_loadings, "LOADINGS") # The unrotated loadings here!
  expect_is(quarti$rot_loadings, "LOADINGS")
  expect_is(bentT$rot_loadings, "LOADINGS")
  expect_is(geoT$rot_loadings, "LOADINGS")
  expect_is(bifacT$rot_loadings, "LOADINGS")

  expect_output(str(equa), "List of 6")
  expect_output(str(equa_1), "List of 6")
  expect_output(str(quarti), "List of 6")
  expect_output(str(bentT), "List of 6")
  expect_output(str(geoT), "List of 6")
  expect_output(str(bifacT), "List of 6")

  expect_is(obli$rotmat, "matrix")
  expect_is(obli$vars_accounted_rot, "matrix")

  expect_is(quarti$rotmat, "matrix")
  expect_is(quarti$vars_accounted_rot, "matrix")

  expect_is(simpli$Phi, "matrix")
  expect_is(simpli$Structure, "matrix")
  expect_is(simpli$rotmat, "matrix")
  expect_is(simpli$vars_accounted_rot, "matrix")

  expect_is(bentQ$Phi, "matrix")
  expect_is(bentQ$Structure, "matrix")
  expect_is(bentQ$rotmat, "matrix")
  expect_is(bentQ$vars_accounted_rot, "matrix")

  expect_is(geoQ$Phi, "matrix")
  expect_is(geoQ$Structure, "matrix")
  expect_is(geoQ$rotmat, "matrix")
  expect_is(geoQ$vars_accounted_rot, "matrix")

  expect_is(bifacQ$Phi, "matrix")
  expect_is(bifacQ$Structure, "matrix")
  expect_is(bifacQ$rotmat, "matrix")
  expect_is(bifacQ$vars_accounted_rot, "matrix")

  expect_equal(obli_1$Phi, NA)
  expect_equal(obli_1$Structure, NA)
  expect_equal(obli_1$rotmat, NA)
  expect_equal(obli_1$vars_accounted_rot, NA)
})

test_that("settings are returned correctly", {
  expect_named(obli$settings, c("kaiser", "precision", "order_type", "k"))
  expect_named(obli_1$settings, c("kaiser", "precision", "order_type", "k"))
  expect_named(quarti$settings, c("kaiser", "precision", "order_type", "k"))
  expect_named(simpli$settings, c("kaiser", "precision", "order_type", "k"))
  expect_named(bentQ$settings, c("kaiser", "precision", "order_type", "k"))
  expect_named(geoQ$settings, c("kaiser", "precision", "order_type", "k"))
  expect_named(bifacQ$settings, c("kaiser", "precision", "order_type", "k"))

  expect_equal(obli$settings$kaiser, TRUE)
  expect_equal(obli_1$settings$kaiser, TRUE)
  expect_equal(quarti$settings$kaiser, TRUE)
  expect_equal(simpli$settings$kaiser, TRUE)
  expect_equal(bentQ$settings$kaiser, TRUE)
  expect_equal(geoQ$settings$kaiser, TRUE)
  expect_equal(bifacQ$settings$kaiser, TRUE)

  expect_equal(obli$settings$precision, 1e-5)
  expect_equal(obli_1$settings$precision, 1e-5)
  expect_equal(quarti$settings$precision, 1e-5)
  expect_equal(simpli$settings$precision, 1e-5)
  expect_equal(bentQ$settings$precision, 1e-5)
  expect_equal(geoQ$settings$precision, 1e-5)
  expect_equal(bifacQ$settings$precision, 1e-5)

  expect_equal(obli$settings$order_type, "eigen")
  expect_equal(obli_1$settings$order_type, "eigen")
  expect_equal(quarti$settings$order_type, "eigen")
  expect_equal(simpli$settings$order_type, "ss_factors")
  expect_equal(bentQ$settings$order_type, "eigen")
  expect_equal(geoQ$settings$order_type, "eigen")
  expect_equal(bifacQ$settings$order_type, "eigen")

  expect_null(obli$settings$k)
  expect_null(obli_1$settings$k)
  expect_null(quarti$settings$k)
  expect_equal(simpli$settings$k, 18)
  expect_null(bentQ$settings$k)
  expect_null(geoQ$settings$k)
  expect_null(bifacQ$settings$k)
})

test_that("errors etc. are thrown correctly", {
  expect_error(.ROTATE_OBLQ(unrot, rotation = "oblimin", type = "none"), ' "order_type" was NULL and no valid "type" was specified. Either use one of "EFAtools", "psych", or "SPSS" for type, or specify the "order_type" argument')

  expect_warning(.ROTATE_OBLQ(unrot, rotation = "oblimin", type = "EFAtools",
                              kaiser = FALSE), " Type and kaiser is specified. kaiser is used with value ' FALSE '. Results may differ from the specified type")
  expect_warning(.ROTATE_OBLQ(unrot, rotation = "oblimin", type = "EFAtools",
                              order_type = "ss_factors"), " Type and order_type is specified. order_type is used with value ' ss_factors '. Results may differ from the specified type")

  expect_warning(.ROTATE_OBLQ(unrot, rotation = "oblimin", type = "psych",
                              kaiser = FALSE), " Type and kaiser is specified. kaiser is used with value ' FALSE '. Results may differ from the specified type")
  expect_warning(.ROTATE_OBLQ(unrot, rotation = "oblimin", type = "psych",
                              order_type = "ss_factors"), " Type and order_type is specified. order_type is used with value ' ss_factors '. Results may differ from the specified type")

  expect_warning(.ROTATE_OBLQ(unrot, rotation = "oblimin", type = "SPSS",
                              kaiser = FALSE), " Type and kaiser is specified. kaiser is used with value ' FALSE '. Results may differ from the specified type")
  expect_warning(.ROTATE_OBLQ(unrot, rotation = "oblimin", type = "SPSS",
                              order_type = "ss_factors"), " Type and order_type is specified. order_type is used with value ' ss_factors '. Results may differ from the specified type")

  expect_warning(.ROTATE_OBLQ(unrot_1, rotation = "oblimin", type = "EFAtools"), " Cannot rotate single factor. Unrotated loadings returned.")
})

rm(unrot, obli, unrot_1, obli_1, quarti, simpli, bentQ, geoQ, bifacQ)
