# A single-factor solution cannot be rotated and has no factor correlations, so
# an oblique rotation request on one factor must degenerate to the no-Phi case on
# every standard-error route (matching a single-fit EFA(), which returns
# Phi = NULL for one factor) rather than crashing on L %*% Phi with Phi = NULL.

cormat <- test_models$baseline$cormat
N_id   <- 500L
m_id   <- 3L
identical_list <- replicate(m_id, cormat, simplify = FALSE)

# A single-fit EFA() with one factor and an oblique rotation: the reference
# behaviour the pooled object must reproduce (no $Phi / $Structure).
single_k1 <- suppressWarnings(suppressMessages(
  EFA(cormat, n_factors = 1L, N = N_id, method = "ML", rotation = "oblimin",
      se = "none")
))


test_that("a single-fit EFA() with one factor and an oblique rotation carries no Phi/Structure", {
  expect_null(single_k1$Phi)
  expect_null(single_k1$Structure)
  expect_false(is.null(single_k1$rot_loadings))
})


test_that("k = 1 oblique pooling does not crash and returns no Phi/Structure (se = none, information)", {
  for (se in c("none", "information")) {
    pooled <- suppressWarnings(suppressMessages(
      EFA_POOLED(identical_list, n_factors = 1L, N = N_id, method = "ML",
                 rotation = "oblimin", se = se)
    ))
    expect_s3_class(pooled, "EFA_POOLED")
    expect_null(pooled$Phi)
    expect_null(pooled$Structure)
    expect_false(is.null(pooled$rot_loadings))
    # Communalities/loadings match the single-fit reference (identical imputations).
    expect_equal(unname(as.numeric(pooled$h2)), unname(as.numeric(single_k1$h2)),
                 tolerance = 1e-6, info = se)
  }
})


test_that("k = 1 oblique pooling fills rotated SE without Phi/Structure SE (se = information)", {
  pooled <- suppressWarnings(suppressMessages(
    EFA_POOLED(identical_list, n_factors = 1L, N = N_id, method = "ML",
               rotation = "oblimin", se = "information")
  ))
  # The rotated-loading SE family is populated (one column), but the factor
  # correlation / structure SE families do not exist for a single factor.
  expect_false(is.null(pooled$SE$rot_loadings))
  expect_identical(dim(as.matrix(pooled$SE$rot_loadings)), c(ncol(cormat), 1L))
  expect_null(pooled$SE$Phi)
  expect_null(pooled$SE$Structure)
  expect_null(pooled$MI$Phi)
})


test_that("k = 1 oblique pooling does not crash on the bootstrap route (se = np-boot)", {
  skip_on_cran()
  set.seed(11)
  raw <- lapply(1:3, function(i) {
    GRiPS_raw[sample(seq_len(nrow(GRiPS_raw)), 200L, replace = TRUE), ]
  })
  pooled <- suppressWarnings(suppressMessages(
    EFA_POOLED(raw, n_factors = 1L, method = "ML", rotation = "oblimin",
               se = "np-boot", b_boot = 20L)
  ))
  expect_s3_class(pooled, "EFA_POOLED")
  expect_null(pooled$Phi)
  expect_null(pooled$Structure)
  expect_false(is.null(pooled$rot_loadings))
  # The bootstrap MI families exist for the loadings, none for Phi/Structure.
  expect_false(is.null(pooled$MI$rot_loadings))
  expect_null(pooled$MI$Phi)
  expect_null(pooled$MI$Structure)
})


test_that("k = 1 oblique pooling does not crash on the sandwich route (MI2S)", {
  skip_on_cran()
  # A few ordinal imputations from a one-factor model.
  mk <- function(seed) {
    set.seed(seed)
    f <- stats::rnorm(300L)
    X <- outer(f, rep(0.6, 6L)) + matrix(stats::rnorm(300L * 6L, sd = 0.85), 300L, 6L)
    d <- as.data.frame(apply(X, 2L, function(col) {
      as.integer(cut(col, c(-Inf, stats::quantile(col, c(.25, .5, .75)), Inf)))
    }))
    names(d) <- paste0("x", seq_len(6L))
    d
  }
  imps <- lapply(1:4, function(i) mk(500L + i))

  pooled <- suppressWarnings(suppressMessages(
    EFA_POOLED(imps, n_factors = 1L, method = "ULS", rotation = "oblimin",
               se = "sandwich", cor_method = "poly")
  ))
  expect_s3_class(pooled, "EFA_POOLED")
  expect_s3_class(pooled$mi_fit, "EFA")
  expect_null(pooled$Phi)
  expect_null(pooled$Structure)
  expect_false(is.null(pooled$rot_loadings))
})
