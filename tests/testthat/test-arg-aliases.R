# The tuning arguments `P_type` and `randomStarts` were renamed to the snake_case
# `p_type` and `random_starts`, and the estimator argument `method` was renamed to
# `estimator`. The old tuning spellings must keep working, silently, and the returned
# object must keep its old settings field names (for the estimator, as a `method` entry
# duplicating `estimator`).

cormat <- test_models$baseline$cormat

test_that("EFA() maps the old randomStarts onto random_starts identically", {
  # The rotation random starts draw from the ambient RNG, so fix it identically before
  # each call (passing `seed` to both calls would do the same).
  set.seed(1)
  new <- suppressWarnings(
    EFA(cormat, n_factors = 3, N = 500, rotation = "oblimin", type = "EFAtools",
        random_starts = 20))
  set.seed(1)
  old <- suppressWarnings(
    EFA(cormat, n_factors = 3, N = 500, rotation = "oblimin", type = "EFAtools",
        randomStarts = 20))
  expect_equal(new$rot_loadings, old$rot_loadings)
  expect_equal(new$settings$randomStarts, old$settings$randomStarts)
  expect_equal(new$settings$randomStarts, 20)
})

test_that("EFA() maps the old P_type onto p_type identically", {
  new <- suppressWarnings(
    EFA(cormat, n_factors = 3, N = 500, rotation = "promax", type = "none",
        init_comm = "smc", criterion = 1e-3, criterion_type = "sum",
        max_iter = 300, abs_eigen = TRUE, normalize = TRUE, order_type = "eigen",
        varimax_type = "kaiser", k = 4, p_type = "norm"))
  old <- suppressWarnings(
    EFA(cormat, n_factors = 3, N = 500, rotation = "promax", type = "none",
        init_comm = "smc", criterion = 1e-3, criterion_type = "sum",
        max_iter = 300, abs_eigen = TRUE, normalize = TRUE, order_type = "eigen",
        varimax_type = "kaiser", k = 4, P_type = "norm"))
  expect_equal(new$rot_loadings, old$rot_loadings)
  expect_equal(new$settings$P_type, old$settings$P_type)
  expect_equal(new$settings$P_type, "norm")
})

test_that("the old spellings add no warning of their own", {
  # Comparing the emitted warnings isolates any deprecation warning from the
  # data-dependent warnings both calls share.
  set.seed(1)
  w_new <- capture_warnings(
    EFA(cormat, n_factors = 3, N = 500, rotation = "oblimin", type = "EFAtools",
        random_starts = 15))
  set.seed(1)
  w_old <- capture_warnings(
    EFA(cormat, n_factors = 3, N = 500, rotation = "oblimin", type = "EFAtools",
        randomStarts = 15))
  expect_equal(w_old, w_new)
})

test_that("EFA() keeps the old settings field names in its output", {
  m_promax <- suppressWarnings(
    EFA(cormat, n_factors = 3, N = 500, rotation = "promax", type = "EFAtools"))
  expect_true("P_type" %in% names(m_promax$settings))
  # the estimator keeps its former field name alongside the current one
  expect_identical(m_promax$settings$method, m_promax$settings$estimator)
  expect_identical(m_promax$settings$method, "PAF")

  m_oblq <- suppressWarnings(
    EFA(cormat, n_factors = 3, N = 500, rotation = "oblimin", type = "EFAtools"))
  expect_true("randomStarts" %in% names(m_oblq$settings))
})

test_that("efa_average() and EFA_AVERAGE() accept p_type / P_type identically", {
  # A single-combination grid returns a plain EFA output, keeping the test fast. The
  # frozen wrapper still selects the estimator with `method`; the successor takes
  # `estimator`.
  args <- list(cormat, n_factors = 3, N = 500, rotation = "promax",
               type = "none", init_comm = "smc", criterion = 1e-3,
               criterion_type = "sum", abs_eigen = TRUE, varimax_type = "kaiser",
               normalize = TRUE, k_promax = 4, start_method = "psych")
  args_new <- c(args, list(estimator = "PAF"))
  args_old <- c(args, list(method = "PAF"))

  avg_new <- suppressWarnings(do.call(efa_average, c(args_new, list(p_type = "norm"))))
  avg_old <- suppressWarnings(do.call(efa_average, c(args_new, list(P_type = "norm"))))
  ea_old  <- suppressWarnings(do.call(EFA_AVERAGE, c(args_old, list(P_type = "norm"))))

  expect_equal(avg_new$rot_loadings, avg_old$rot_loadings)
  expect_equal(avg_new$rot_loadings, ea_old$rot_loadings)
  # the output settings keep the old P_type field name.
  expect_true("P_type" %in% names(avg_new$settings))
  expect_equal(avg_new$settings$P_type, "norm")
})
