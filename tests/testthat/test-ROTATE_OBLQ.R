unrot <- EFA(test_models$baseline$cormat, 3, N = 500)
obli <- .rotate_model(unrot, rotation = "oblimin", type = "EFAtools")

unrot_1 <- EFA(test_models$baseline$cormat, 1, N = 500)
obli_1 <- suppressWarnings(.rotate_model(unrot_1, rotation = "oblimin",
                                        type = "EFAtools"))

quarti <- suppressWarnings(.rotate_model(unrot, rotation = "quartimin", type = "psych"))
simpli <- suppressWarnings(.rotate_model(unrot, rotation = "simplimax", type = "SPSS",
                       maxit = 2000))
bentQ <- suppressWarnings(.rotate_model(unrot, rotation = "bentlerQ", type = "none",
                       order_type = "eigen"))
geoQ <- suppressWarnings(.rotate_model(unrot, rotation = "geominQ", type = "EFAtools"))
bifacQ <- suppressWarnings(.rotate_model(unrot, rotation = "bifactorQ", type = "EFAtools"))

test_that("output class and dimensions are correct", {
  expect_s3_class(obli$rot_loadings, "LOADINGS")
  expect_s3_class(obli_1$rot_loadings, "LOADINGS") # The unrotated loadings here!
  expect_s3_class(quarti$rot_loadings, "LOADINGS")
  expect_s3_class(simpli$rot_loadings, "LOADINGS")
  expect_s3_class(bentQ$rot_loadings, "LOADINGS")
  expect_s3_class(geoQ$rot_loadings, "LOADINGS")
  expect_s3_class(bifacQ$rot_loadings, "LOADINGS")

  expect_output(str(obli), "List of 6")
  expect_output(str(obli_1), "List of 6")
  expect_output(str(quarti), "List of 6")
  expect_output(str(simpli), "List of 6")
  expect_output(str(bentQ), "List of 6")
  expect_output(str(geoQ), "List of 6")
  expect_output(str(bifacQ), "List of 6")

  checkmate::expect_matrix(obli$Phi)
  expect_s3_class(obli$Structure, "LOADINGS")
  checkmate::expect_matrix(obli$rotmat)
  checkmate::expect_matrix(obli$vars_accounted_rot)

  checkmate::expect_matrix(quarti$Phi)
  expect_s3_class(quarti$Structure, "LOADINGS")
  checkmate::expect_matrix(quarti$rotmat)
  checkmate::expect_matrix(quarti$vars_accounted_rot)

  checkmate::expect_matrix(simpli$Phi)
  expect_s3_class(simpli$Structure, "LOADINGS")
  checkmate::expect_matrix(simpli$rotmat)
  checkmate::expect_matrix(simpli$vars_accounted_rot)

  checkmate::expect_matrix(bentQ$Phi)
  expect_s3_class(bentQ$Structure, "LOADINGS")
  checkmate::expect_matrix(bentQ$rotmat)
  checkmate::expect_matrix(bentQ$vars_accounted_rot)

  checkmate::expect_matrix(geoQ$Phi)
  expect_s3_class(geoQ$Structure, "LOADINGS")
  checkmate::expect_matrix(geoQ$rotmat)
  checkmate::expect_matrix(geoQ$vars_accounted_rot)

  checkmate::expect_matrix(bifacQ$Phi)
  expect_s3_class(bifacQ$Structure, "LOADINGS")
  checkmate::expect_matrix(bifacQ$rotmat)
  checkmate::expect_matrix(bifacQ$vars_accounted_rot)

  expect_null(obli_1$Phi)
  expect_null(obli_1$Structure)
  expect_equal(obli_1$rotmat, NA)
  expect_null(obli_1$vars_accounted_rot)
})

test_that("settings are returned correctly", {
  expect_named(obli$settings, c("normalize", "precision", "order_type", "k", "randomStarts"))
  expect_named(obli_1$settings, c("normalize", "precision", "order_type", "k", "randomStarts"))
  expect_named(quarti$settings, c("normalize", "precision", "order_type", "k", "randomStarts"))
  expect_named(simpli$settings, c("normalize", "precision", "order_type", "k", "randomStarts"))
  expect_named(bentQ$settings, c("normalize", "precision", "order_type", "k", "randomStarts"))
  expect_named(geoQ$settings, c("normalize", "precision", "order_type", "k", "randomStarts"))
  expect_named(bifacQ$settings, c("normalize", "precision", "order_type", "k", "randomStarts"))

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
  expect_error(.rotate_model(unrot, rotation = "oblimin", type = "none"), class = "efa_type_none")

  expect_warning(.rotate_model(unrot, rotation = "oblimin", type = "EFAtools",
                              normalize = FALSE), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "oblimin", type = "EFAtools",
                              order_type = "ss_factors"), class = "efa_type_override")

  expect_warning(.rotate_model(unrot, rotation = "oblimin", type = "psych",
                              normalize = FALSE), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "oblimin", type = "psych",
                              order_type = "ss_factors"), class = "efa_type_override")

  expect_warning(.rotate_model(unrot, rotation = "oblimin", type = "SPSS",
                              normalize = FALSE), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "oblimin", type = "SPSS",
                              order_type = "ss_factors"), class = "efa_type_override")

  expect_warning(.rotate_model(unrot_1, rotation = "oblimin", type = "EFAtools"), class = "efa_single_factor")
})

test_that("factors, Phi, structure, and rotmat are reordered consistently", {
  # ss_factors and eigen use the same ordering criterion, so the oblique solution
  # (loadings, Phi, structure, and rotation matrix) is identical, and the structure
  # matrix equals pattern %*% Phi.
  set.seed(42)
  eig <- suppressWarnings(.rotate_model(unrot, rotation = "oblimin", type = "EFAtools"))
  set.seed(42)
  ss <- suppressWarnings(.rotate_model(unrot, rotation = "oblimin", type = "SPSS"))

  expect_equal(unclass(eig$Structure), unclass(eig$rot_loadings) %*% eig$Phi,
               ignore_attr = TRUE)
  expect_equal(unclass(ss$rot_loadings), unclass(eig$rot_loadings), ignore_attr = TRUE)
  expect_equal(ss$Phi, eig$Phi, ignore_attr = TRUE)
  expect_equal(unclass(ss$Structure), unclass(eig$Structure), ignore_attr = TRUE)
  expect_equal(ss$rotmat, eig$rotmat, ignore_attr = TRUE)
})

test_that("oblique Phi, structure, and rotmat are reflected/reordered with the loadings", {
  skip_on_cran()
  skip_if_not_installed("GPArotation")

  # Reflect and reorder the engine's own solution the documented way, then check the
  # rotated output matches it exactly. This pins the sign reflection of Phi (a factor
  # in this solution has a negative column sum) and the consistent reflection and
  # reordering of Phi and the rotation matrix, which a structure == pattern %*% Phi
  # check alone misses.
  L <- unrot$unrot_loadings
  set.seed(11)
  ref <- suppressWarnings(GPArotation::oblimin(L, eps = 1e-5, normalize = TRUE,
                                               randomStarts = 100))
  signs <- sign(colSums(ref$loadings))
  signs[signs == 0] <- 1
  ord <- order(colSums((ref$loadings %*% diag(signs))^2), decreasing = TRUE)
  exp_load <- (ref$loadings %*% diag(signs))[, ord]
  exp_phi <- (diag(signs) %*% ref$Phi %*% diag(signs))[ord, ord]
  exp_rotmat <- (ref$Th %*% diag(signs))[, ord]

  set.seed(11)
  o <- suppressWarnings(.rotate_model(unrot, rotation = "oblimin", type = "EFAtools"))

  expect_false(all(signs == 1))  # the reflection path is exercised
  expect_equal(unclass(o$rot_loadings), exp_load, ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(o$Phi, exp_phi, ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(o$rotmat, exp_rotmat, ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(unclass(o$Structure), exp_load %*% exp_phi, ignore_attr = TRUE,
               tolerance = 1e-6)
  # the rotation matrix reproduces the rotated pattern: L_unrot %*% t(solve(rotmat))
  expect_equal(unclass(L) %*% t(solve(o$rotmat)), unclass(o$rot_loadings),
               ignore_attr = TRUE, tolerance = 1e-6)
})

rm(unrot, obli, unrot_1, obli_1, quarti, simpli, bentQ, geoQ, bifacQ)
