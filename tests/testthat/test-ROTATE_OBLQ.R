unrot <- EFA(test_models$baseline$cormat, 3, N = 500)
obli <- .ROTATE_OBLQ(unrot, rotation = "oblimin", type = "EFAtools")

unrot_1 <- EFA(test_models$baseline$cormat, 1, N = 500)
obli_1 <- suppressWarnings(.ROTATE_OBLQ(unrot_1, rotation = "oblimin",
                                        type = "EFAtools"))

quarti <- .ROTATE_OBLQ(unrot, rotation = "quartimin", type = "psych")
simpli <- .ROTATE_OBLQ(unrot, rotation = "simplimax", type = "SPSS",
                       maxit = 2000)
bentQ <- .ROTATE_OBLQ(unrot, rotation = "bentlerQ", type = "none",
                       order_type = "eigen")
geoQ <- .ROTATE_OBLQ(unrot, rotation = "geominQ", type = "EFAtools")
bifacQ <- .ROTATE_OBLQ(unrot, rotation = "bifactorQ", type = "EFAtools")

test_that("output class and dimensions are correct", {
  expect_is(obli$rot_loadings, "LOADINGS")
  expect_is(obli_1$rot_loadings, "LOADINGS") # The unrotated loadings here!
  expect_is(quarti$rot_loadings, "LOADINGS")
  expect_is(simpli$rot_loadings, "LOADINGS")
  expect_is(bentQ$rot_loadings, "LOADINGS")
  expect_is(geoQ$rot_loadings, "LOADINGS")
  expect_is(bifacQ$rot_loadings, "LOADINGS")

  expect_output(str(obli), "List of 6")
  expect_output(str(obli_1), "List of 6")
  expect_output(str(quarti), "List of 6")
  expect_output(str(simpli), "List of 6")
  expect_output(str(bentQ), "List of 6")
  expect_output(str(geoQ), "List of 6")
  expect_output(str(bifacQ), "List of 6")

  expect_is(obli$Phi, "matrix")
  expect_is(obli$Structure, "matrix")
  expect_is(obli$rotmat, "matrix")
  expect_is(obli$vars_accounted_rot, "matrix")

  expect_is(quarti$Phi, "matrix")
  expect_is(quarti$Structure, "matrix")
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
  expect_named(obli$settings, c("normalize", "precision", "order_type", "k"))
  expect_named(obli_1$settings, c("normalize", "precision", "order_type", "k"))
  expect_named(quarti$settings, c("normalize", "precision", "order_type", "k"))
  expect_named(simpli$settings, c("normalize", "precision", "order_type", "k"))
  expect_named(bentQ$settings, c("normalize", "precision", "order_type", "k"))
  expect_named(geoQ$settings, c("normalize", "precision", "order_type", "k"))
  expect_named(bifacQ$settings, c("normalize", "precision", "order_type", "k"))

  expect_equal(obli$settings$normalize, TRUE)
  expect_equal(obli_1$settings$normalize, TRUE)
  expect_equal(quarti$settings$normalize, TRUE)
  expect_equal(simpli$settings$normalize, TRUE)
  expect_equal(bentQ$settings$normalize, TRUE)
  expect_equal(geoQ$settings$normalize, TRUE)
  expect_equal(bifacQ$settings$normalize, TRUE)

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

  expect_equal(obli$settings$k, NA)
  expect_equal(obli_1$settings$k, NA)
  expect_equal(quarti$settings$k, NA)
  expect_equal(simpli$settings$k, 18)
  expect_equal(bentQ$settings$k, NA)
  expect_equal(geoQ$settings$k, NA)
  expect_equal(bifacQ$settings$k, NA)
})

test_that("errors etc. are thrown correctly", {
  expect_error(.ROTATE_OBLQ(unrot, rotation = "oblimin", type = "none"), ' "order_type" was NA and no valid "type" was specified. Either use one of "EFAtools", "psych", or "SPSS" for type, or specify the "order_type" argument\n')

  expect_warning(.ROTATE_OBLQ(unrot, rotation = "oblimin", type = "EFAtools",
                              normalize = FALSE), " Type and normalize is specified. normalize is used with value ' FALSE '. Results may differ from the specified type\n")
  expect_warning(.ROTATE_OBLQ(unrot, rotation = "oblimin", type = "EFAtools",
                              order_type = "ss_factors"), " Type and order_type is specified. order_type is used with value ' ss_factors '. Results may differ from the specified type\n")

  expect_warning(.ROTATE_OBLQ(unrot, rotation = "oblimin", type = "psych",
                              normalize = FALSE), " Type and normalize is specified. normalize is used with value ' FALSE '. Results may differ from the specified type\n")
  expect_warning(.ROTATE_OBLQ(unrot, rotation = "oblimin", type = "psych",
                              order_type = "ss_factors"), " Type and order_type is specified. order_type is used with value ' ss_factors '. Results may differ from the specified type\n")

  expect_warning(.ROTATE_OBLQ(unrot, rotation = "oblimin", type = "SPSS",
                              normalize = FALSE), " Type and normalize is specified. normalize is used with value ' FALSE '. Results may differ from the specified type\n")
  expect_warning(.ROTATE_OBLQ(unrot, rotation = "oblimin", type = "SPSS",
                              order_type = "ss_factors"), " Type and order_type is specified. order_type is used with value ' ss_factors '. Results may differ from the specified type\n")

  expect_warning(.ROTATE_OBLQ(unrot_1, rotation = "oblimin", type = "EFAtools"), " Cannot rotate single factor. Unrotated loadings returned.\n")
})

rm(unrot, obli, unrot_1, obli_1, quarti, simpli, bentQ, geoQ, bifacQ)
