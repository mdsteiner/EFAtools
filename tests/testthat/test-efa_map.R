map_cor <- efa_map(test_models$baseline$cormat)
map_raw <- efa_map(GRiPS_raw)

# Three uncorrelated triplets whose within-triplet correlations are near-unit but
# *different* from each other. The leading eigenvalue (2.99998) is then simple and
# well separated from the next two (2.80, 2.60), so its eigenvector is determined up
# to sign on every LAPACK driver and concentrates on the first triplet. That
# triplet's residual variance is 6.7e-06 after one component, below the 1e-5 guard,
# so the partialling degenerates at m = 1 on any platform. Do not equalise the three
# correlations: that makes the leading eigenvalues numerically tied, and any
# orthonormal basis of the resulting three-dimensional eigenspace is then a valid
# result -- a balanced eigenvector leaves every residual variance at 0.667 and the
# search does not truncate at all.
map_trunc_cor <- local({
  R <- diag(9)
  for (b in 0:2) {
    idx <- (b * 3 + 1):(b * 3 + 3)
    R[idx, idx] <- c(.99999, .90, .80)[b + 1]
  }
  diag(R) <- 1
  R
})

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
  expect_named(tr2, c("name", "label", "n_factors", "plot_type", "x", "y",
                      "m_last"))
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

test_that("an ordinary run searches the whole grid and stays silent", {
  # Partialling out all but one component leaves a rank-one residual, so the final
  # point m = p - 1 is undefined on this matrix; that is the expected end of the
  # grid and must not raise the truncation warning.
  expect_no_warning(res <- efa_map(test_models$baseline$cormat),
                    class = "efa_map_truncated")
  expect_equal(.retention_record(res, "TR2")$m_last,
               ncol(test_models$baseline$cormat) - 2)
  expect_null(res$note)

  # GRiPS reaches the last point, so nothing is NA there
  expect_equal(.retention_record(map_raw, "TR2")$m_last, ncol(GRiPS_raw) - 1)
  expect_false(any(is.na(.retention_record(map_raw, "TR2")$y)))
})

test_that("a truncated grid warns, keeps the computed prefix, and records how far it got", {
  expect_warning(res <- efa_map(map_trunc_cor), class = "efa_map_truncated")

  for (nm in c("TR2", "TR4")) {
    rec <- .retention_record(res, nm)
    # the search stopped well short of the grid, which is what the warning reports
    expect_equal(rec$m_last, 0)
    expect_lt(rec$m_last, ncol(map_trunc_cor) - 2)
    # the prefix that could be computed is kept, the rest stays NA -- asserted
    # against m_last itself, so the contract holds wherever the stop lands
    kept <- seq_len(rec$m_last + 1)
    expect_false(any(is.na(rec$y[kept])))
    expect_true(all(is.na(rec$y[-kept])))
    # the suggestion is the minimum over the evaluated range
    expect_equal(res$n_factors[[nm]], rec$x[which.min(rec$y)])
  }

  # the print says so too, rather than presenting a fully searched grid
  expect_length(res$note, 1)
  expect_match(res$note, "not over the full grid", fixed = TRUE)
})

test_that("a truncated MAP does not take the criterion out of efa_retain()", {
  expect_warning(
    nf <- efa_retain(map_trunc_cor, N = 500, suitability = FALSE,
                     criteria = "MAP"),
    class = "efa_map_truncated"
  )
  expect_named(nf$outputs, "MAP")
  expect_null(nf$not_run)
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

rm(map_cor, map_raw, map_trunc_cor)
