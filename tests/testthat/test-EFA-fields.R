# Schema contract for the top-level inference slots on EFA: SE, CI, replicates.
# Pins which keys live in each slot under each (se method x rotation) combination,
# and the absence of accidental dual storage under the old `boot.` names.

dat <- GRiPS_raw[, 1:6]
R   <- stats::cor(dat)
N   <- nrow(dat)

# Assert the slot contract for one fit. Old `boot.` names must be absent
# everywhere, the expected SE/CI keys must be present, and `replicates` is a
# present-but-NULL slot under analytic SE methods (a list under "np-boot").
expect_field_contract <- function(fit, expected_SE, expected_CI,
                                  expected_replicates, info = "") {
  expect_false(any(c("boot.SE", "boot.CI", "boot.arrays", "boot.MI")
                   %in% names(fit)), info = info)
  expect_setequal(names(fit$SE), expected_SE)
  expect_setequal(names(fit$CI), expected_CI)
  if (is.null(expected_replicates)) {
    expect_true("replicates" %in% names(fit), info = info)
    expect_null(fit$replicates, info = info)
  } else {
    expect_true(is.list(fit$replicates), info = info)
    expect_setequal(names(fit$replicates), expected_replicates)
  }
}


test_that("se = 'none' leaves no SE/CI/replicates/MI slot on the fit", {
  fit <- EFA(R, n_factors = 2, N = N, method = "ML", rotation = "none",
             se = "none")
  expect_false(any(c("SE", "CI", "replicates", "MI",
                     "boot.SE", "boot.CI", "boot.arrays", "boot.MI")
                   %in% names(fit)))
})


test_that("information SEs fill the schema across all rotation families", {
  # Unrotated: unrot_loadings + uniquenesses; no rotated-loading components.
  fit_none <- EFA(R, n_factors = 2, N = N, method = "ML", rotation = "none",
                  se = "information")
  expect_field_contract(
    fit_none,
    expected_SE = c("unrot_loadings", "uniquenesses"),
    expected_CI = c("unrot_loadings", "uniquenesses"),
    expected_replicates = NULL,
    info = "information / none"
  )

  # Orthogonal: adds rot_loadings + communalities; no Phi/Structure.
  fit_orth <- EFA(R, n_factors = 2, N = N, method = "ML", rotation = "varimax",
                  se = "information")
  expect_field_contract(
    fit_orth,
    expected_SE = c("unrot_loadings", "uniquenesses",
                    "rot_loadings", "communalities"),
    expected_CI = c("unrot_loadings", "uniquenesses",
                    "rot_loadings", "communalities"),
    expected_replicates = NULL,
    info = "information / varimax"
  )

  # Oblique: adds Phi and Structure on top of the orthogonal schema.
  fit_oblq <- EFA(R, n_factors = 2, N = N, method = "ML", rotation = "oblimin",
                  se = "information")
  expect_field_contract(
    fit_oblq,
    expected_SE = c("unrot_loadings", "uniquenesses",
                    "rot_loadings", "communalities", "Phi", "Structure"),
    expected_CI = c("unrot_loadings", "uniquenesses",
                    "rot_loadings", "communalities", "Phi", "Structure"),
    expected_replicates = NULL,
    info = "information / oblimin"
  )
})


test_that("sandwich SEs fill the same schema as information across rotations", {
  skip_on_cran()
  # The Browne (1984) ADF Gamma for Pearson correlations needs the raw data, so
  # these fits use `dat` instead of the precomputed cormat.
  fit_none <- EFA(dat, n_factors = 2, method = "ML", rotation = "none",
                  se = "sandwich")
  expect_field_contract(
    fit_none,
    expected_SE = c("unrot_loadings", "uniquenesses"),
    expected_CI = c("unrot_loadings", "uniquenesses"),
    expected_replicates = NULL,
    info = "sandwich / none"
  )

  fit_orth <- EFA(dat, n_factors = 2, method = "ML", rotation = "varimax",
                  se = "sandwich")
  expect_field_contract(
    fit_orth,
    expected_SE = c("unrot_loadings", "uniquenesses",
                    "rot_loadings", "communalities"),
    expected_CI = c("unrot_loadings", "uniquenesses",
                    "rot_loadings", "communalities"),
    expected_replicates = NULL,
    info = "sandwich / varimax"
  )

  fit_oblq <- EFA(dat, n_factors = 2, method = "ML", rotation = "oblimin",
                  se = "sandwich")
  expect_field_contract(
    fit_oblq,
    expected_SE = c("unrot_loadings", "uniquenesses",
                    "rot_loadings", "communalities", "Phi", "Structure"),
    expected_CI = c("unrot_loadings", "uniquenesses",
                    "rot_loadings", "communalities", "Phi", "Structure"),
    expected_replicates = NULL,
    info = "sandwich / oblimin"
  )
})


test_that("np-boot fills `replicates` and adds fit_indices/residuals SEs", {
  skip_on_cran()
  # Small `b` -- this test pins schema/slot names; the replicate count enters only
  # via dim(...)[3L] == b, not via any numeric agreement.
  b <- 4L

  # Unrotated: SE/CI cover the unrotated loadings, fit indices, and residuals;
  # `replicates` is a list of three 3D cubes whose last dimension is b.
  set.seed(42)
  fit_none <- suppressWarnings(suppressMessages(
    EFA(dat, n_factors = 2, method = "ML", rotation = "none",
        se = "np-boot", b_boot = b)
  ))
  expect_field_contract(
    fit_none,
    expected_SE = c("unrot_loadings", "fit_indices", "residuals"),
    expected_CI = c("unrot_loadings", "fit_indices", "residuals"),
    expected_replicates = c("unrot_loadings", "fit_indices", "residuals"),
    info = "np-boot / none"
  )
  expect_identical(dim(fit_none$replicates$unrot_loadings)[3L], b)

  # Orthogonal: adds rot_loadings; SE also reports valid_target_rotations (the
  # CI list does not). `replicates` mirrors SE minus that integer counter.
  set.seed(42)
  fit_orth <- suppressWarnings(suppressMessages(
    EFA(dat, n_factors = 2, method = "ML", rotation = "varimax",
        se = "np-boot", b_boot = b)
  ))
  expect_field_contract(
    fit_orth,
    expected_SE = c("unrot_loadings", "rot_loadings",
                    "fit_indices", "residuals", "valid_target_rotations"),
    expected_CI = c("unrot_loadings", "rot_loadings",
                    "fit_indices", "residuals"),
    expected_replicates = c("unrot_loadings", "rot_loadings",
                            "fit_indices", "residuals"),
    info = "np-boot / varimax"
  )
  expect_identical(dim(fit_orth$replicates$rot_loadings)[3L], b)

  # Oblique: adds Phi and Structure on top of the orthogonal schema.
  set.seed(42)
  fit_oblq <- suppressWarnings(suppressMessages(
    EFA(dat, n_factors = 2, method = "ML", rotation = "oblimin",
        se = "np-boot", b_boot = b)
  ))
  expect_field_contract(
    fit_oblq,
    expected_SE = c("unrot_loadings", "rot_loadings", "Phi", "Structure",
                    "fit_indices", "residuals", "valid_target_rotations"),
    expected_CI = c("unrot_loadings", "rot_loadings", "Phi", "Structure",
                    "fit_indices", "residuals"),
    expected_replicates = c("unrot_loadings", "rot_loadings", "Phi",
                            "Structure", "fit_indices", "residuals"),
    info = "np-boot / oblimin"
  )
  expect_identical(dim(fit_oblq$replicates$Phi)[3L], b)
})
