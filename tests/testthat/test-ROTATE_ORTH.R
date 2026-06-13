unrot <- EFA(test_models$baseline$cormat, 3, N = 500)
equa <- .rotate_model(unrot, rotation = "equamax", type = "EFAtools")

unrot_1 <- EFA(test_models$baseline$cormat, 1, N = 500)
equa_1 <- suppressWarnings(.rotate_model(unrot_1, rotation = "equamax",
                                        type = "EFAtools"))

quarti <- .rotate_model(unrot, rotation = "quartimax", type = "psych")
bentT <- .rotate_model(unrot, rotation = "bentlerT", type = "none",
                      order_type = "eigen")
geoT <- .rotate_model(unrot, rotation = "geominT", type = "SPSS")
bifacT <- .rotate_model(unrot, rotation = "bifactorT", type = "EFAtools")

test_that("output class and dimensions are correct", {
  expect_s3_class(equa$rot_loadings, "LOADINGS")
  expect_s3_class(equa_1$rot_loadings, "LOADINGS") # The unrotated loadings here!
  expect_s3_class(quarti$rot_loadings, "LOADINGS")
  expect_s3_class(bentT$rot_loadings, "LOADINGS")
  expect_s3_class(geoT$rot_loadings, "LOADINGS")
  expect_s3_class(bifacT$rot_loadings, "LOADINGS")

  expect_output(str(equa), "List of 4")
  expect_output(str(equa_1), "List of 4")
  expect_output(str(quarti), "List of 4")
  expect_output(str(bentT), "List of 4")
  expect_output(str(geoT), "List of 4")
  expect_output(str(bifacT), "List of 4")

  checkmate::expect_matrix(equa$rotmat)
  checkmate::expect_matrix(equa$vars_accounted_rot)

  checkmate::expect_matrix(quarti$rotmat)
  checkmate::expect_matrix(quarti$vars_accounted_rot)

  checkmate::expect_matrix(bentT$rotmat)
  checkmate::expect_matrix(bentT$vars_accounted_rot)

  checkmate::expect_matrix(geoT$rotmat)
  checkmate::expect_matrix(geoT$vars_accounted_rot)

  checkmate::expect_matrix(bifacT$rotmat)
  checkmate::expect_matrix(bifacT$vars_accounted_rot)

  expect_equal(equa_1$rotmat, NA)
  expect_equal(equa_1$vars_accounted_rot, NA)
})

test_that("settings are returned correctly", {
  expect_named(equa$settings, c("normalize", "precision", "order_type", "randomStarts"))
  expect_named(equa_1$settings,c("normalize", "precision", "order_type", "randomStarts"))
  expect_named(quarti$settings, c("normalize", "precision", "order_type", "randomStarts"))
  expect_named(bentT$settings, c("normalize", "precision", "order_type", "randomStarts"))
  expect_named(geoT$settings, c("normalize", "precision", "order_type", "randomStarts"))
  expect_named(bifacT$settings, c("normalize", "precision", "order_type", "randomStarts"))

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
  expect_error(.rotate_model(unrot, rotation = "equamax", type = "none"), class = "efa_type_none")

  expect_warning(.rotate_model(unrot, rotation = "equamax", type = "EFAtools",
                              normalize = FALSE), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "equamax", type = "EFAtools",
                              order_type = "ss_factors"), class = "efa_type_override")

  expect_warning(.rotate_model(unrot, rotation = "equamax", type = "psych",
                              normalize = FALSE), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "equamax", type = "psych",
                              order_type = "ss_factors"), class = "efa_type_override")

  expect_warning(.rotate_model(unrot, rotation = "equamax", type = "SPSS",
                              normalize = FALSE), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "equamax", type = "SPSS",
                              order_type = "ss_factors"), class = "efa_type_override")

  expect_warning(.rotate_model(unrot_1, rotation = "equamax", type = "EFAtools"), class = "efa_single_factor")
})

test_that("orthogonal rotation matrix reproduces the rotated loadings", {
  skip_on_cran()
  skip_if_not_installed("GPArotation")

  # the engine reflects two factors on this solution, so the rotation matrix must carry
  # both the sign reflection and the factor reordering. Reflect and reorder the engine's
  # own solution the documented way, check the package's rotation matrix matches that
  # (signed, reordered) target, then confirm the self-contained reproduction identity
  # L_unrot %*% rotmat == rot_loadings.
  L <- unrot$unrot_loadings
  set.seed(11)
  ref <- suppressWarnings(GPArotation::bentlerT(unclass(L), eps = 1e-5,
                                                normalize = TRUE, randomStarts = 100))
  signs <- sign(colSums(ref$loadings))
  signs[signs == 0] <- 1
  ord <- order(colSums((ref$loadings %*% diag(signs))^2), decreasing = TRUE)
  exp_rotmat <- (ref$Th %*% diag(signs))[, ord]

  set.seed(11)
  o <- suppressWarnings(.rotate_model(unrot, rotation = "bentlerT", type = "psych"))

  expect_false(all(signs == 1))  # this fixture exercises the reflection path
  expect_equal(o$rotmat, exp_rotmat, ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(unclass(L) %*% o$rotmat, unclass(o$rot_loadings),
               ignore_attr = TRUE, tolerance = 1e-6)
})

rm(unrot, equa, unrot_1, equa_1, quarti, bentT, geoT, bifacT)
