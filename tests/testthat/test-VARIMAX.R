unrot <- EFA(test_models$baseline$cormat, 3, N = 500)
vari <- .VARIMAX(unrot, type = "EFAtools")

unrot_1 <- EFA(test_models$baseline$cormat, 1, N = 500)
vari_1 <- .VARIMAX(unrot_1, type = "EFAtools")

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

test_that("errors etc. are thrown correctly", {
  expect_warning(.VARIMAX(unrot_1, type = "EFAtools"), " Cannot rotate single factor. Unrotated loadings returned.")
})

rm(unrot, vari, unrot_1, vari_1)

