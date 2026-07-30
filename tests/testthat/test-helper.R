test_that(".efa_num formats numbers, normalises negative zero, returns plain character", {
  # strips the leading zero of values in (-1, 1) and left-pads to a common width
  expect_equal(.efa_num(0.23, digits = 2), " .23")
  expect_equal(.efa_num(0.2345, digits = 3), " .234")
  expect_equal(.efa_num(0.2345, digits = 3, print_zero = TRUE), " 0.234")

  # a value rounding to zero from below is rendered without a sign, like its
  # positive counterpart (no negative zero)
  expect_equal(.efa_num(-0.0001, digits = 3), .efa_num(0.0001, digits = 3))
  expect_equal(.efa_num(-0.0001, digits = 3, print_zero = TRUE),
               .efa_num(0.0001, digits = 3, print_zero = TRUE))
  expect_false(any(grepl("-", .efa_num(-0.0001, digits = 3), fixed = TRUE)))

  # the padded result is always plain character, never a cli ansi_string
  padded <- .efa_num(0.5, digits = 3, pad = TRUE)
  expect_type(padded, "character")
  expect_false(inherits(padded, "cli_ansi_string"))
})


x_base <- population_models$loadings$baseline
x_NA <- population_models$loadings$baseline
x_NA[1, 3] <- NA
y_base <- x_base[, c(3,2,1)]
y_NA <- y_base
y_NA[2, 2] <- NA

test_that(".factor_congruence works", {
  checkmate::expect_matrix(.factor_congruence(x_base, y_base))
  expect_equal(sum(.factor_congruence(x_base, y_base)), 3)
  expect_warning(.factor_congruence(x_NA, y_NA, na.rm = FALSE), class = "efa_missing_check")
  expect_warning(.factor_congruence(x_NA, y_NA), class = "efa_missing_complete")
})

test_that(".factor_congruence has (i, j) = cos(x_i, y_j) orientation", {
  # Asymmetric inputs so the congruence matrix differs from its transpose; this
  # pins the orientation the efa_compare / averaging reordering depends on.
  A <- matrix(c(2, 0, 1,
                0, 1, 0,
                0, 0, 3,
                1, 2, 0), nrow = 4, byrow = TRUE)
  B <- matrix(c(0, 1, 0,
                3, 0, 1,
                0, 0, 2,
                1, 1, 0), nrow = 4, byrow = TRUE)
  # Reference computed entrywise from column dot products, independent of the
  # matrix-algebra implementation.
  ref <- sapply(seq_len(ncol(B)), function(j)
    sapply(seq_len(ncol(A)), function(i)
      sum(A[, i] * B[, j]) / sqrt(sum(A[, i]^2) * sum(B[, j]^2))))
  expect_false(isTRUE(all.equal(ref, t(ref))))
  expect_equal(unname(.factor_congruence(A, B)), ref)
})

test_that(".factor_congruence aborts when too few complete cases remain", {
  x_few <- population_models$loadings$baseline
  x_few[-1, 1] <- NA   # only the first row has complete data
  expect_error(
    suppressWarnings(.factor_congruence(x_few, x_base)),
    class = "efa_too_few_complete"
  )
})

test_that(".is_cormat works", {
  expect_equal(.is_cormat(cor(cbind(rnorm(100), rnorm(100)))), TRUE)
  expect_equal(.is_cormat(cbind(rnorm(100), rnorm(100))), FALSE)
  expect_equal(.is_cormat(cbind(rnorm(2), rnorm(2))), FALSE)
  expect_equal(.is_cormat(cbind(c(1, NA, .57, .85))), FALSE)
  expect_equal(.is_cormat(matrix(c(1, .1, .3, 1), ncol = 2)), FALSE)
  expect_error(.is_cormat(matrix(c(1, NA, NA, 1), ncol = 2)),
               class = "efa_cormat_has_na")
  # A symmetric matrix with entries in [-1, 1] but a non-unit diagonal (here trace 2.6, which
  # rounds to the row count) is not a correlation matrix and must not be detected as one.
  expect_equal(
    .is_cormat(matrix(c(0.6, 0.2, 0.2, 0.2, 1, 0.2, 0.2, 0.2, 1), ncol = 3)),
    FALSE
  )
})

q_p <- .det_max_factors(8) + 1
test_that(".det_max_factors works", {
  expect_type(.det_max_factors(8), "double")
  expect_lte(((8 - q_p)**2 - (8 + q_p)) / 2, 0)
  expect_equal(.det_max_factors(0), 0)
  expect_equal(.det_max_factors(1), 0)
  expect_equal(.det_max_factors(2), 0)
  expect_equal(.det_max_factors(3), 0)
  expect_gt(.det_max_factors(4), 0)
})


test_that(".decimals works", {
  expect_type(.decimals(8), "double")
  expect_equal(.decimals(8), 0)
  expect_type(.decimals(8), "double")
  expect_error(.decimals("a"), class = "efa_not_numeric")

  # NA elements must not trigger "missing value where TRUE/FALSE needed"; they
  # contribute 0 decimals so the rest of the input still drives the result.
  expect_equal(.decimals(NA_real_), 0)
  expect_equal(.decimals(c(1.23, NA, 4.5)), 2)
  expect_equal(.decimals(matrix(c(1.23, NA, 4.5, 6), ncol = 2)), 2)

  # Small magnitudes that R renders in scientific notation (e.g. "3e-05") used to
  # abort with "subscript out of bounds"; the decimal count must be taken from the
  # fixed-notation representation instead.
  expect_equal(.decimals(3e-5), 5)
  expect_equal(.decimals(1e-4), 4)
  expect_equal(.decimals(1e-6), 6)
  expect_equal(.decimals(1.23e-5), 7)
  expect_equal(.decimals(c(0.5, 3e-5)), 5)
  expect_equal(.decimals(matrix(c(0.5, 3e-5, 0.4, 0.2), ncol = 2)), 5)

  # A comma-decimal OutDec must not collapse the count to 0: the separator is
  # pinned to "." so the split still finds the fractional part.
  withr::with_options(list(OutDec = ","), {
    expect_equal(.decimals(0.123), 3)
    expect_equal(.decimals(3e-5), 5)
  })
})

### test .rotation_family

test_that(".rotation_family works", {
  expect_identical(.rotation_family("none"), "none")

  # representative names hardcoded so the test catches a rotation placed in the
  # wrong family vector (including the special-cased varimax/promax and the
  # geomin/bentler/bifactor T-vs-Q pairs)
  expect_identical(.rotation_family("varimax"), "orthogonal")
  expect_identical(.rotation_family("geominT"), "orthogonal")
  expect_identical(.rotation_family("bifactorT"), "orthogonal")
  expect_identical(.rotation_family("promax"), "oblique")
  expect_identical(.rotation_family("geominQ"), "oblique")
  expect_identical(.rotation_family("bifactorQ"), "oblique")

  # and every canonical name dispatches to its family
  for (rot in .orth_rotations) {
    expect_identical(.rotation_family(rot), "orthogonal")
  }

  for (rot in .oblq_rotations) {
    expect_identical(.rotation_family(rot), "oblique")
  }

  expect_error(.rotation_family("bogus"), class = "efa_unknown_rotation")
})



rm(x_base, y_base, q_p)


test_that(".smc_start returns the squared multiple correlations", {
  R <- test_models$baseline$cormat
  expect_equal(.smc_start(R), 1 - 1 / diag(solve(R)))
})

test_that(".three_eigen returns eigenvalues per requested diagonal convention", {
  R <- test_models$baseline$cormat

  res <- .three_eigen(R, c("PCA", "SMC"))

  # PCA: eigenvalues on the unit-diagonal correlation matrix
  expect_equal(res$PCA, eigen(R, symmetric = TRUE, only.values = TRUE)$values)

  # SMC: eigenvalues with the squared multiple correlations on the diagonal
  R_smc <- R
  diag(R_smc) <- .smc_start(R)
  expect_equal(res$SMC, eigen(R_smc, symmetric = TRUE, only.values = TRUE)$values)

  # a convention that was not requested is NA
  expect_true(is.na(res$EFA))

  # EFA: eigenvalues with the communalities of an n-factor EFA on the diagonal
  res_efa <- .three_eigen(R, "EFA", n_factors = 2)
  R_efa <- R
  diag(R_efa) <- suppressMessages(suppressWarnings(EFA(R, n_factors = 2)$h2))
  expect_equal(res_efa$EFA, eigen(R_efa, symmetric = TRUE, only.values = TRUE)$values)
  expect_true(is.na(res_efa$PCA))
})
