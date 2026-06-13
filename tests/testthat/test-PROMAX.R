unrot <- EFA(test_models$baseline$cormat, 3, N = 500)
prom <- .rotate_model(unrot, rotation = "promax", type = "EFAtools")
prom_psych <- .rotate_model(unrot, rotation = "promax", type = "psych")
prom_spss <- .rotate_model(unrot, rotation = "promax", type = "SPSS")

unrot_1 <- EFA(test_models$baseline$cormat, 1, N = 500)
prom_1 <- suppressWarnings(.rotate_model(unrot_1, rotation = "promax", type = "EFAtools"))

test_that("output class and dimensions are correct", {
  expect_type(prom, "list")
  expect_type(prom_1, "list")
  expect_named(prom, c("rot_loadings", "Phi", "Structure", "rotmat",
                       "vars_accounted_rot", "settings"))
  expect_named(prom_1, c("rot_loadings", "Phi", "Structure", "rotmat",
                       "vars_accounted_rot", "settings"))

  expect_s3_class(prom$rot_loadings, "LOADINGS")
  checkmate::expect_matrix(prom$Phi)
  expect_s3_class(prom$Structure, "LOADINGS")
  checkmate::expect_matrix(prom$rotmat)
  checkmate::expect_matrix(prom$vars_accounted_rot)
  expect_type(prom$settings, "list")

  expect_s3_class(prom_1$rot_loadings, "LOADINGS")
  expect_equal(prom_1$Phi, NA)
  expect_equal(prom_1$Structure, NA)
  expect_equal(prom_1$rotmat, NA)
  expect_equal(prom_1$vars_accounted_rot, NA)
  expect_type(prom_1$settings, "list")
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

  expect_equal(prom$settings$varimax_type, "kaiser")
  expect_equal(prom_psych$settings$varimax_type, "svd")
  expect_equal(prom_spss$settings$varimax_type, "kaiser")
  expect_equal(prom_1$settings$varimax_type, "kaiser")

  expect_equal(prom$settings$k, 4)
  expect_equal(prom_psych$settings$k, 4)
  expect_equal(prom_spss$settings$k, 4)
  expect_equal(prom_1$settings$k, 4)

})

test_that("errors etc. are thrown correctly", {

  expect_error(.rotate_model(unrot, rotation = "promax", type = "none"), class = "efa_type_none")

  expect_warning(.rotate_model(unrot, rotation = "promax", type = "EFAtools", normalize = FALSE), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "promax", type = "EFAtools", P_type = "norm"), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "promax", type = "EFAtools", order_type = "ss_factors"), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "promax", type = "EFAtools", k = 2), class = "efa_type_override")

  expect_warning(.rotate_model(unrot, rotation = "promax", type = "psych", normalize = FALSE), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "promax", type = "psych", P_type = "norm"), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "promax", type = "psych", order_type = "ss_factors"), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "promax", type = "psych", k = 2), class = "efa_type_override")

  expect_warning(.rotate_model(unrot, rotation = "promax", type = "SPSS", normalize = FALSE), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "promax", type = "SPSS", P_type = "unnorm"), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "promax", type = "SPSS", order_type = "eigen"), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "promax", type = "SPSS", k = 2), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "promax", type = "SPSS", varimax_type = "svd"), class = "efa_type_override")
  expect_warning(.rotate_model(unrot_1, rotation = "promax", type = "EFAtools"), class = "efa_single_factor")
})

test_that("promax rotation matrix reproduces the rotated loadings", {
  # the rotation matrix must carry the sign reflection and factor reordering so that
  # L_unrot %*% rotmat == rot_loadings, for both ordering branches (eigen reorders the
  # pattern after the fit; ss_factors reorders the varimax base before it)
  L <- unclass(unrot$unrot_loadings)
  expect_equal(L %*% prom$rotmat, unclass(prom$rot_loadings),       # order_type "eigen"
               ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(L %*% prom_spss$rotmat, unclass(prom_spss$rot_loadings),  # "ss_factors"
               ignore_attr = TRUE, tolerance = 1e-6)
})

rm(unrot, prom, unrot_1, prom_1, prom_psych, prom_spss)
