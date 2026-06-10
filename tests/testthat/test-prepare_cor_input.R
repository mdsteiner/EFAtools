# Unit tests for the shared correlation-input preparer: correlation-matrix
# detection, the N-handling policy, and the singular / non-positive-definite
# branches.

x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z), ncol = 3)               # raw, collinear -> singular
cor_sing <- stats::cor(dat_sing)                       # correlation matrix, singular
cor_nposdef <- matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), ncol = 3)  # invertible but non-PD
cormat <- test_models$baseline$cormat                  # well-behaved correlation matrix

test_that("a correlation matrix is detected and returned unchanged", {
  prep <- .prepare_cor_input(cormat, N = 500)

  expect_named(prep, c("R", "N", "is_cormat"))
  expect_true(prep$is_cormat)
  expect_equal(prep$R, cormat)
  expect_equal(prep$N, 500)

  # a clean correlation matrix triggers no message, warning, or error
  expect_silent(.prepare_cor_input(cormat, N = 500))
})

test_that("raw data is converted to a correlation matrix", {
  expect_message(prep <- .prepare_cor_input(GRiPS_raw), class = "efa_cor_from_data")

  expect_false(prep$is_cormat)
  expect_equal(prep$R,
               stats::cor(GRiPS_raw, use = "pairwise.complete.obs",
                          method = "pearson"))
  expect_equal(colnames(prep$R), colnames(GRiPS_raw))
  expect_equal(prep$N, nrow(GRiPS_raw))

  # inform_from_data = FALSE silences the "computed from raw data" message
  expect_no_message(.prepare_cor_input(GRiPS_raw, inform_from_data = FALSE))
})

test_that("N_policy governs how N is handled", {
  # required: a correlation matrix without N aborts, but a supplied N is kept
  expect_error(.prepare_cor_input(cormat, N_policy = "required"),
               class = "efa_n_required")
  prep_req <- .prepare_cor_input(cormat, N = 200, N_policy = "required")
  expect_equal(prep_req$N, 200)

  # none: N is passed through untouched and the both-N-and-raw warning is muted
  expect_no_warning(
    prep_none <- suppressMessages(
      .prepare_cor_input(GRiPS_raw, N = 999, N_policy = "none")))
  expect_equal(prep_none$N, 999)

  # optional (default): supplying N with raw data warns and N is taken from data
  expect_warning(
    prep_opt <- suppressMessages(.prepare_cor_input(GRiPS_raw, N = 999)),
    class = "efa_n_from_data")
  expect_equal(prep_opt$N, nrow(GRiPS_raw))
})

test_that("a non-positive-definite matrix is smoothed, or aborts under posdef_abort", {
  # default: the matrix is smoothed to positive-definiteness (psych warns)
  prep_sm <- suppressWarnings(.prepare_cor_input(cor_nposdef, N = 10))
  expect_true(prep_sm$is_cormat)
  expect_true(all(eigen(prep_sm$R, symmetric = TRUE, only.values = TRUE)$values > 0))

  # posdef_abort = TRUE turns the same input into an error instead
  expect_error(.prepare_cor_input(cor_nposdef, N = 10, posdef_abort = TRUE),
               class = "efa_cor_not_posdef")
})

test_that("a singular matrix aborts unless the check is disabled", {
  expect_error(.prepare_cor_input(cor_sing, N = 10), class = "efa_cor_singular")
  expect_error(suppressMessages(.prepare_cor_input(dat_sing)),
               class = "efa_cor_singular")

  # check_singular = FALSE skips the abort; the matrix is smoothed instead
  prep_ns <- suppressWarnings(
    .prepare_cor_input(cor_sing, N = 10, check_singular = FALSE))
  expect_true(all(eigen(prep_ns$R, symmetric = TRUE, only.values = TRUE)$values > 0))
})
