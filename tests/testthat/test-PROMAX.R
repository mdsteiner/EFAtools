unrot <- EFA(test_models$baseline$cormat, 3, N = 500)
prom <- .PROMAX(unrot, type = "EFAtools")

unrot_1 <- EFA(test_models$baseline$cormat, 1, N = 500)
prom_1 <- .PROMAX(unrot_1, type = "EFAtools")

test_that("output class and dimensions are correct", {
  expect_is(prom, "list")
  expect_is(prom_1, "list")
  expect_named(prom, c("rot_loadings", "Phi", "Structure", "rotmat",
                       "vars_accounted_rot", "settings"))
  expect_named(prom_1, c("rot_loadings", "Phi", "Structure", "rotmat",
                       "vars_accounted_rot", "settings"))

  expect_is(prom$rot_loadings, "LOADINGS")
  expect_is(prom$Phi, "matrix")
  expect_is(prom$Structure, "matrix")
  expect_is(prom$rotmat, "matrix")
  expect_is(prom$vars_accounted_rot, "matrix")
  expect_is(prom$settings, "list")

  expect_is(prom_1$rot_loadings, "LOADINGS")
  expect_equal(prom_1$Phi, NA)
  expect_equal(prom_1$Structure, NA)
  expect_equal(prom_1$rotmat, NA)
  expect_equal(prom_1$vars_accounted_rot, NA)
  expect_is(prom_1$settings, "list")
})

test_that("errors etc. are thrown correctly", {
  expect_warning(.PROMAX(unrot_1, type = "EFAtools"), " Cannot rotate single factor. Unrotated loadings returned.")
})

rm(unrot, prom, unrot_1, prom_1)
