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

  # Feed .reflect_and_order a raw oblique solution whose first factor has a negative column
  # sum, and check that the documented sign reflection is applied consistently to the
  # loadings, the factor intercorrelations, and the rotation matrix (a
  # structure == pattern %*% Phi check alone misses the Phi/rotmat reflection). This
  # exercises the package's own reflection branch directly and deterministically, rather
  # than relying on a rotation engine to return a negative-column-sum factor (an orientation
  # that varies across GPArotation versions). The raw input is a real oblimin solution with
  # its first factor negated -- still a valid oblique solution (loadings = L %*% t(solve(Th)),
  # Phi = t(Th) %*% Th) -- so reflecting it back must recover the canonical solution.
  L <- unrot$unrot_loadings
  o <- suppressWarnings(.rotate_model(unrot, rotation = "oblimin", type = "EFAtools"))
  k <- ncol(o$rot_loadings)

  neg <- replace(rep(1, k), 1L, -1)
  raw_loadings <- unclass(o$rot_loadings) %*% diag(neg)
  raw_Phi      <- diag(neg) %*% o$Phi %*% diag(neg)
  raw_rotmat   <- o$rotmat %*% diag(neg)

  res <- .reflect_and_order(raw_loadings, Phi = raw_Phi, rotmat = raw_rotmat,
                            L_unrot = L, name_factors = FALSE)

  # the negative factor is reflected back to a non-negative column sum
  expect_true(all(colSums(unclass(res$rot_loadings)) >= 0))
  # loadings, Phi, and the rotation matrix are reflected consistently, recovering the
  # canonical solution
  expect_equal(unclass(res$rot_loadings), unclass(o$rot_loadings), ignore_attr = TRUE)
  expect_equal(res$Phi, o$Phi, ignore_attr = TRUE)
  expect_equal(res$rotmat, o$rotmat, ignore_attr = TRUE)
  # the structure matrix stays pattern %*% Phi, and the rotation matrix reproduces the
  # rotated pattern via L_unrot %*% t(solve(rotmat))
  expect_equal(unclass(res$Structure), unclass(res$rot_loadings) %*% res$Phi,
               ignore_attr = TRUE)
  expect_equal(unclass(L) %*% t(solve(res$rotmat)), unclass(res$rot_loadings),
               ignore_attr = TRUE, tolerance = 1e-6)
})

test_that("the bentlerQ solution satisfies the oblique structure invariants", {
  skip_on_cran()

  # bentlerQ is computed by the native gradient-projection engine. These invariants hold for any
  # valid oblique solution and need no reference package (so coverage survives the criterion
  # moving off GPArotation): the structure matrix equals pattern %*% Phi, the rotation matrix
  # reproduces the rotated pattern via the documented identity L_unrot %*% t(solve(Th)), and the
  # factor correlation matrix has a unit diagonal.
  L <- unrot$unrot_loadings

  expect_equal(unclass(bentQ$Structure), unclass(bentQ$rot_loadings) %*% bentQ$Phi,
               ignore_attr = TRUE)
  expect_equal(unclass(L) %*% t(solve(bentQ$rotmat)), unclass(bentQ$rot_loadings),
               ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(diag(bentQ$Phi), rep(1, ncol(bentQ$rot_loadings)))
})

test_that("the bifactorQ solution satisfies the oblique structure invariants", {
  skip_on_cran()

  # bifactorQ is computed by the native gradient-projection engine. These invariants hold for any
  # valid oblique solution and need no reference package (so coverage survives the criterion moving
  # off GPArotation): the structure matrix equals pattern %*% Phi, the rotation matrix reproduces
  # the rotated pattern via the documented identity L_unrot %*% t(solve(Th)), and the factor
  # correlation matrix has a unit diagonal.
  L <- unrot$unrot_loadings

  expect_equal(unclass(bifacQ$Structure), unclass(bifacQ$rot_loadings) %*% bifacQ$Phi,
               ignore_attr = TRUE)
  expect_equal(unclass(L) %*% t(solve(bifacQ$rotmat)), unclass(bifacQ$rot_loadings),
               ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(diag(bifacQ$Phi), rep(1, ncol(bifacQ$rot_loadings)))
})

rm(unrot, obli, unrot_1, obli_1, quarti, simpli, bentQ, geoQ, bifacQ)
