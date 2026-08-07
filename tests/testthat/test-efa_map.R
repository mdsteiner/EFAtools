map_cor <- efa_map(test_models$baseline$cormat)
map_raw <- efa_map(GRiPS_raw)

test_that("output class and dimensions are correct", {
  expect_s3_class(map_cor, "efa_retention")
  expect_length(map_cor, 6)
  expect_s3_class(map_raw, "efa_retention")
  expect_length(map_raw, 6)

  expect_named(map_cor$n_factors, c("TR2", "TR4"))
  expect_equal(.retention_record(map_cor, "TR2")$plot_type, "none")
})

test_that("identified number of factors is correct", {
  expect_equal(map_cor$n_factors[["TR2"]], 1)
  expect_equal(map_cor$n_factors[["TR4"]], 3)
})

test_that("criterion series are returned", {
  tr2 <- .retention_record(map_cor, "TR2")
  expect_length(tr2$y, ncol(test_models$baseline$cormat))
  expect_length(tr2$x, ncol(test_models$baseline$cormat))
  expect_equal(tr2$x[1], 0)
  # the suggested m minimizes the criterion
  expect_equal(tr2$x[which.min(tr2$y)], map_cor$n_factors[["TR2"]])
})

test_that("the TR2 criterion series matches psych::vss", {
  skip_if_not_installed("psych")
  # Velicer's original MAP is TR2, which psych::vss returns as $map. EFAtools
  # additionally reports m = 0 (part of Velicer's procedure, omitted by psych), and
  # returns NA once the partialling degenerates at the largest m, where psych
  # substitutes 1; compare the m >= 1 values both report. `n.obs` feeds psych's
  # other statistics only -- MAP is a function of the correlation matrix alone --
  # so one nominal value covers both matrices.
  for (cmat in list(test_models$baseline$cormat, stats::cor(GRiPS_raw))) {
    ours <- .retention_record(efa_map(cmat), "TR2")$y[-1]
    theirs <- psych::vss(cmat, n = ncol(cmat) - 1, n.obs = 500,
                         plot = FALSE)$map

    compared <- !is.na(ours)
    # guard against a vacuously passing comparison
    expect_gt(sum(compared), 0.5 * length(ours))
    # The two sides compute the same quantity through different operation sequences --
    # ours partials with A = V sqrt(Lambda) and takes the trace of a matrix product,
    # psych's sums the squared residuals elementwise off its own principal() fit -- so
    # their rounding does not cancel. The partial matrix is standardised by
    # 1 / sqrt(d_i d_j), which at the largest m amplifies those last-digit differences
    # before they are squared and summed. A definitional divergence (a wrong partialling
    # order, a missing offset) moves the criterion by ~1e-2, so 1e-8 still catches one.
    expect_equal(ours[compared], theirs[compared], tolerance = 1e-8)
  }
})

test_that("TR4 is the trace of the fourth matrix power", {
  # For a three-variable equicorrelation matrix with r = .2, the eigenvalues
  # are 1.4, .8, .8. Hence tr(R^4) = 1.4^4 + 2 * .8^4 and the m = 0
  # revised-MAP value has the following closed form. It is far from the
  # element-wise fourth-power average (.2^4), so this guards that distinction.
  R <- matrix(.2, 3, 3)
  diag(R) <- 1
  tr4_expected <- (1.4^4 + 2 * .8^4 - 3) / (3 * 2)

  tr4 <- .retention_record(efa_map(R), "TR4")$y[1]
  expect_equal(tr4, tr4_expected, tolerance = 1e-14)
  expect_gt(abs(tr4 - .2^4), .1)
})

test_that("settings are returned correctly", {
  expect_named(map_cor$settings, c("use", "cor_method"))
  expect_equal(map_cor$settings$use, "pairwise.complete.obs")
  expect_equal(map_cor$settings$cor_method, "pearson")
})

test_that("errors are thrown correctly", {
  expect_error(efa_map(1:5), class = "efa_input_not_matrix")
  expect_message(efa_map(GRiPS_raw), class = "efa_cor_from_data")
})

rm(map_cor, map_raw)
