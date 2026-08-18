
ekc_cor <- efa_ekc(test_models$baseline$cormat, N = 500)
ekc_raw <- efa_ekc(GRiPS_raw)

test_that("output class and dimensions are correct", {
  expect_s3_class(ekc_cor, "efa_retention")
  expect_length(ekc_cor, 6)
  expect_s3_class(ekc_raw, "efa_retention")
  expect_length(ekc_raw, 6)

  expect_equal(unname(ekc_cor$criterion["id"]), "EKC")
  expect_length(ekc_cor$results, 1)
  expect_named(ekc_cor$n_factors, "BvA2017")
  expect_equal(ekc_cor$results[[1]]$plot_type, "eigen")
})


test_that("found eigenvalues are correct", {
  expect_equal(sum(ekc_cor$results[[1]]$y),
               ncol(test_models$baseline$cormat))
  expect_equal(sum(ekc_raw$results[[1]]$y), ncol(GRiPS_raw))
  expect_length(ekc_cor$results[[1]]$y,
                ncol(test_models$baseline$cormat))
  expect_length(ekc_raw$results[[1]]$y, ncol(GRiPS_raw))
})

test_that("reference eigenvalues are correct", {
  expect_equal(ekc_cor$results[[1]]$reference[floor(ncol(test_models$baseline$cormat) / 2)], 1)
  expect_equal(ekc_raw$results[[1]]$reference[floor(ncol(GRiPS_raw) / 2)], 1)
  expect_length(ekc_cor$results[[1]]$reference,
                length(ekc_cor$results[[1]]$y))
  expect_length(ekc_raw$results[[1]]$reference, length(ekc_raw$results[[1]]$y))
})

test_that("identified number of factors is correct", {
  expect_equal(ekc_cor$n_factors[["BvA2017"]], 3)
  expect_equal(ekc_raw$n_factors[["BvA2017"]], 1)
})

test_that("the BvA2017 reference series reproduces the published equation", {
  # Braeken and van Assen (2017, p. 454) define the j-th reference eigenvalue as
  # max(((J - sum_{i < j} lambda_i) / (J - j + 1)) * (1 + sqrt(J / N))^2, 1),
  # transcribed here as an explicit loop -- an independent expression of the
  # equation the vectorised cumsum form in efa_ekc() implements. The equation holds
  # for any N, so one nominal value is used for both matrices.
  N <- 500
  for (cmat in list(test_models$baseline$cormat, stats::cor(GRiPS_raw))) {
    lambda <- eigen(cmat, symmetric = TRUE, only.values = TRUE)$values
    J <- ncol(cmat)

    ref <- numeric(J)
    for (j in seq_len(J)) {
      used <- if (j == 1) 0 else sum(lambda[seq_len(j - 1)])
      ref[j] <- max((J - used) / (J - j + 1) * (1 + sqrt(J / N))^2, 1)
    }

    out <- efa_ekc(cmat, N = N)
    expect_equal(.retention_record(out, "BvA2017")$reference, ref)
    # the paper's stopping rule: retain the leading run of eigenvalues above their
    # reference value
    expect_equal(out$n_factors[["BvA2017"]], sum(cumprod(lambda > ref)))
  }
})

test_that("the retained count is always an integer", {
  # The count is derived with cumprod(), which coerces the comparison to double, so pin
  # that it comes back in the integer storage mode a count belongs in. An empty solution
  # is checked alongside a retained one because a single expression has to cover both.
  # (How the count is rendered no longer depends on this: the report and the plot label
  # both go through .retention_count(); see test-efa_retention.R.)
  some <- efa_ekc(test_models$baseline$cormat, N = 500)
  expect_identical(some$results[[1]]$n_factors, 3L)
  expect_identical(some$results[[1]]$highlight, 3L)

  # a sample size just above the number of variables lifts the reference series above every
  # eigenvalue, so nothing is retained
  none <- efa_ekc(test_models$case_1a$cormat, N = 11)
  expect_identical(none$results[[1]]$n_factors, 0L)
  expect_null(none$results[[1]]$highlight)
})

test_that("the deprecated type argument warns and is ignored", {
  cmat <- test_models$baseline$cormat

  expect_warning(efa_ekc(cmat, N = 500, type = "AM2019"),
                 class = "lifecycle_warning_deprecated")
  expect_warning(efa_ekc(cmat, N = 500, type = c("BvA2017", "AM2019")),
                 class = "lifecycle_warning_deprecated")
  expect_warning(efa_retain(cmat, N = 500, criteria = "EKC", suitability = FALSE,
                            ekc_type = "AM2019"),
                 class = "lifecycle_warning_deprecated")

  # ignored: whatever is passed, the result is the one the argument-free call gives
  deprecated <- suppressWarnings(efa_ekc(cmat, N = 500, type = "AM2019"))
  expect_equal(deprecated, ekc_cor)

  # the criteria themselves never trip the deprecation, and the argument-free calls
  # stay silent
  expect_no_condition(efa_ekc(cmat, N = 500))
  expect_no_condition(efa_retain(cmat, N = 500, criteria = "EKC", suitability = FALSE))

  # the superseded wrappers keep both arguments, and a call that does not name one stays
  # silent
  expect_no_condition(EKC(cmat, N = 500))
  expect_warning(EKC(cmat, N = 500, type = "AM2019"),
                 class = "lifecycle_warning_deprecated")
  expect_no_condition(N_FACTORS(cmat, N = 500, criteria = "EKC", suitability = FALSE))
  expect_warning(N_FACTORS(cmat, N = 500, criteria = "EKC", suitability = FALSE,
                           ekc_type = "AM2019"),
                 class = "lifecycle_warning_deprecated")

  # a wrapper that warns but still forwards a value would pass every assertion above, so
  # pin what the wrappers return as well
  expect_equal(suppressWarnings(EKC(cmat, N = 500, type = "AM2019")),
               EKC(cmat, N = 500))
  expect_equal(
    suppressWarnings(N_FACTORS(cmat, N = 500, criteria = "EKC", suitability = FALSE,
                               ekc_type = "AM2019"))$n_factors,
    N_FACTORS(cmat, N = 500, criteria = "EKC", suitability = FALSE)$n_factors
  )
})

test_that("the deprecation warning repeats and is not attributed to the package", {
  skip_if_not_installed("withr")
  # lifecycle decides from the calling environment whether a user named the deprecated
  # argument, or whether a package named it in its own internals. Only the second case is
  # rate limited, and only it gets a footer that asks the reader to report a bug against
  # the package. Every route into the warning enters from this package, the direct call as
  # much as the one through a superseded wrapper, so without the caller EFAtools declares
  # all of them would be read as the second case.
  #
  # lifecycle exempts the package under test from that rule, which hides the difference
  # from the rest of the suite. Clearing TESTTHAT_PKG removes the exemption, so the calls
  # below take the path a user's session takes. The repeat is what separates the two: a
  # rate-limited warning is signalled at most once.
  withr::local_envvar(c(TESTTHAT_PKG = ""))
  cmat <- test_models$baseline$cormat

  count_two <- function(call) {
    fired <- 0L
    for (i in 1:2) {
      withCallingHandlers(call(), lifecycle_warning_deprecated = function(w) {
        fired <<- fired + 1L
        invokeRestart("muffleWarning")
      })
    }
    fired
  }

  # both routes, because both resolve to this package
  expect_identical(count_two(function() efa_ekc(cmat, N = 500, type = "AM2019")), 2L)
  expect_identical(count_two(function() EKC(cmat, N = 500, type = "AM2019")), 2L)
})

test_that("errors are thrown correctly", {
  expect_error(efa_ekc(1:5), class = "efa_input_not_matrix")
  expect_error(efa_ekc(test_models$baseline$cormat), class = "efa_n_required")
  expect_message(efa_ekc(GRiPS_raw), class = "efa_cor_from_data")
  expect_warning(efa_ekc(GRiPS_raw, N = 20), class = "efa_n_from_data")
  expect_error(efa_ekc(sing_raw), class = "efa_cor_singular")
  expect_error(efa_ekc(sing_cor, N = sing_N), class = "efa_cor_singular")
  expect_warning(efa_ekc(cor_nposdef, N = 20), class = "efa_cor_smoothed")
})

test_that("EKC requires more observations than variables", {
  # the EKC reference series (1 + sqrt(J / N))^2 is only defined above that boundary
  R8 <- stats::cor(GRiPS_raw)  # 8 variables
  expect_error(efa_ekc(R8, N = 8), class = "efa_n_too_small")
  expect_s3_class(efa_ekc(R8, N = 9), "efa_retention")
})

test_that("settings are returned correctly", {
  expect_named(ekc_cor$settings, c("use", "cor_method", "N"))
  expect_named(ekc_raw$settings, c("use", "cor_method", "N"))

  expect_equal(ekc_cor$settings$N, 500)
  expect_equal(ekc_raw$settings$N, 810)

  expect_equal(ekc_cor$settings$use, "pairwise.complete.obs")
  expect_equal(ekc_raw$settings$use, "pairwise.complete.obs")

  expect_equal(ekc_cor$settings$cor_method, "pearson")
  expect_equal(ekc_raw$settings$cor_method, "pearson")
})

rm(ekc_cor, ekc_raw)
