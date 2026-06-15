unrot <- EFA(test_models$baseline$cormat, 3, N = 500)
vari <- .rotate_model(unrot, rotation = "varimax", type = "EFAtools")
vari_psych <- .rotate_model(unrot, rotation = "varimax", type = "psych")
vari_spss <- .rotate_model(unrot, rotation = "varimax", type = "SPSS")

unrot_1 <- EFA(test_models$baseline$cormat, 1, N = 500)
vari_1 <- suppressWarnings(.rotate_model(unrot_1, rotation = "varimax", type = "EFAtools"))

test_that("output class and dimensions are correct", {
  expect_type(vari, "list")
  expect_type(vari_1, "list")
  expect_named(vari, c("rot_loadings", "rotmat", "vars_accounted_rot", "settings"))
  expect_named(vari_1, c("rot_loadings", "rotmat", "vars_accounted_rot", "settings"))

  expect_s3_class(vari$rot_loadings, "LOADINGS")
  checkmate::expect_matrix(vari$rotmat)
  checkmate::expect_matrix(vari$vars_accounted_rot)
  expect_type(vari$settings, "list")

  expect_s3_class(vari_1$rot_loadings, "LOADINGS")
  expect_equal(vari_1$rotmat, NA)
  expect_null(vari_1$vars_accounted_rot)
  expect_type(vari_1$settings, "list")
})

test_that("settings are returned correctly", {
  expect_named(vari$settings, c("normalize", "precision", "order_type",
                                "varimax_type"))
  expect_named(vari_psych$settings, c("normalize", "precision", "order_type",
                                      "varimax_type"))
  expect_named(vari_spss$settings, c("normalize", "precision", "order_type",
                                     "varimax_type"))
  expect_named(vari_1$settings, c("normalize", "precision", "order_type",
                                  "varimax_type"))

  expect_equal(vari$settings$normalize, TRUE)
  expect_equal(vari_psych$settings$normalize, TRUE)
  expect_equal(vari_spss$settings$normalize, TRUE)
  expect_equal(vari_1$settings$normalize, TRUE)

  expect_equal(vari$settings$precision, 1e-05)
  expect_equal(vari_psych$settings$precision, 1e-05)
  expect_equal(vari_spss$settings$precision, 1e-5)
  expect_equal(vari_1$settings$precision, 1e-05)

  expect_equal(vari$settings$order_type, "eigen")
  expect_equal(vari_psych$settings$order_type, "eigen")
  expect_equal(vari_spss$settings$order_type, "ss_factors")
  expect_equal(vari_1$settings$order_type, "eigen")

  expect_equal(vari$settings$varimax_type, "kaiser")
  expect_equal(vari_psych$settings$varimax_type, "svd")
  expect_equal(vari_spss$settings$varimax_type, "kaiser")
  expect_equal(vari_1$settings$varimax_type, "kaiser")

})

test_that("errors etc. are thrown correctly", {

  expect_error(.rotate_model(unrot, rotation = "varimax", type = "none"), class = "efa_type_none")

  expect_warning(.rotate_model(unrot, rotation = "varimax", type = "EFAtools", normalize = FALSE), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "varimax", type = "EFAtools", order_type = "ss_factors"), class = "efa_type_override")

  expect_warning(.rotate_model(unrot, rotation = "varimax", type = "psych", normalize = FALSE), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "varimax", type = "psych", order_type = "ss_factors"), class = "efa_type_override")

  expect_warning(.rotate_model(unrot, rotation = "varimax", type = "SPSS", normalize = FALSE), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "varimax", type = "SPSS", order_type = "eigen"), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "varimax", type = "SPSS", varimax_type = "svd"), class = "efa_type_override")
  expect_warning(.rotate_model(unrot_1, rotation = "varimax", type = "EFAtools"), class = "efa_single_factor")
})

test_that(".VARIMAX_SPSS tolerates a zero-communality row (Kaiser normalisation)", {
  # A zero-communality row makes 1/sqrt(h2) = Inf; without floored weights Kaiser
  # normalisation injects NaN that corrupts the rotation or triggers an
  # uninformative if(NA) error. The floored row must stay zero throughout.
  L <- matrix(c(0.8, 0.1,
                0.7, 0.2,
                0.0, 0.0,
                0.1, 0.9,
                0.2, 0.8), ncol = 2, byrow = TRUE)
  res <- .VARIMAX_SPSS(L, normalize = TRUE)
  expect_false(anyNA(res$loadings))
  expect_equal(res$loadings[3, ], c(0, 0))
})

rm(unrot, vari, unrot_1, vari_1, vari_psych, vari_spss)

