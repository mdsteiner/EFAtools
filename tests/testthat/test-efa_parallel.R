# Each efa_parallel() call simulates ~100 datasets and eigen-decomposes them, so this
# fixture block dominates the file (~50s). Skipped by default; opt in with
# `Sys.setenv(EFATOOLS_TEST_SLOW = "true")` to run. See helper-slow.R.

# The kernel-routing check below compares two separate eigen-decomposition runs, so it uses
# expect_equal() with an explicit tolerance rather than expect_identical(): a threaded BLAS
# (Apple's Accelerate, for one) is free to vary its GEMM reduction order between calls, so
# one function on one input can differ in the last ulp. waldo still compares S3 classes,
# names, attributes and structure exactly, so only numeric values are given slack -- a draw
# that did not route through the shared kernel consumes a different random-number stream and
# produces wholly different eigenvalues, so it still fails loudly. (A set tolerance does relax
# integer against double, which these eigenvalues never are.)
fp_tol <- 1e-8

if (is_slow_test()) {
# seed the parallel-analysis simulation so the retained-factor counts are
# reproducible (future_lapply uses future.seed = TRUE under the sequential plan)
set.seed(42)
pa_cor <- efa_parallel(test_models$baseline$cormat, N = 500)
pa_cor_pca <- efa_parallel(test_models$baseline$cormat, N = 500, eigen_type = "PCA")
pa_raw <- efa_parallel(GRiPS_raw)
pa_nodat <- efa_parallel(N = 20, n_vars = 5)
pa_craw <- efa_parallel(test_models$baseline$cormat, N = 500, eigen_type = "PCA",
                    decision_rule = "crawford")
pa_perc <- efa_parallel(test_models$baseline$cormat, N = 500, eigen_type = "PCA",
                    decision_rule = "percentile")
}  # is_slow_test()

test_that("output class and dimensions are correct", {
  skip_if_not_slow()
  expect_s3_class(pa_cor, "efa_retention")
  expect_length(pa_cor, 6)
  expect_s3_class(pa_raw, "efa_retention")
  expect_length(pa_raw, 6)
  expect_s3_class(pa_cor_pca, "efa_retention")
  expect_s3_class(pa_nodat, "efa_retention")
  expect_s3_class(pa_craw, "efa_retention")
  expect_s3_class(pa_perc, "efa_retention")

  expect_named(pa_cor$n_factors, c("PCA", "SMC", "EFA"))
  expect_named(pa_cor_pca$n_factors, "PCA")
  expect_named(pa_nodat$n_factors, c("PCA", "SMC", "EFA"))
  expect_equal(.retention_record(pa_cor, "PCA")$plot_type, "eigen")
})


test_that("found eigenvalues are correct", {
  skip_if_not_slow()
  # real eigenvalues form the solid line (record $y)
  expect_equal(sum(.retention_record(pa_cor, "PCA")$y),
               ncol(test_models$baseline$cormat))
  expect_lt(sum(.retention_record(pa_cor, "SMC")$y),
            ncol(test_models$baseline$cormat))
  expect_lt(sum(.retention_record(pa_cor, "EFA")$y),
            ncol(test_models$baseline$cormat))

  # simulated eigenvalues form the dashed reference series (record $references)
  expect_equal(sum(.retention_record(pa_cor, "PCA")$references$Means),
               ncol(test_models$baseline$cormat))
  expect_gt(sum(.retention_record(pa_cor, "PCA")$references[["95 Percentile"]]),
            ncol(test_models$baseline$cormat))
  expect_lt(sum(.retention_record(pa_cor, "SMC")$references$Means),
            ncol(test_models$baseline$cormat))
  expect_lt(sum(.retention_record(pa_cor, "EFA")$references$Means),
            ncol(test_models$baseline$cormat))

  expect_equal(sum(.retention_record(pa_raw, "PCA")$y), ncol(GRiPS_raw))
  expect_named(.retention_record(pa_raw, "PCA")$references,
               c("Means", "95 Percentile"))
  expect_equal(sum(.retention_record(pa_raw, "PCA")$references$Means),
               ncol(GRiPS_raw))
  expect_gt(sum(.retention_record(pa_raw, "PCA")$references[["95 Percentile"]]),
            ncol(GRiPS_raw))
  expect_lt(sum(.retention_record(pa_raw, "SMC")$y), ncol(GRiPS_raw))
  expect_lt(sum(.retention_record(pa_raw, "EFA")$y), ncol(GRiPS_raw))

  # references is a named two-series list (means + percentile)
  refs <- .retention_record(pa_cor, "PCA")$references
  checkmate::expect_list(refs, len = 2)
  expect_named(refs, c("Means", "95 Percentile"))

  # only the requested eigenvalue type produces a record
  expect_null(.retention_record(pa_cor_pca, "SMC"))
  expect_null(.retention_record(pa_cor_pca, "EFA"))

  # no real data: no real-eigenvalue series, but the simulated references remain
  expect_null(.retention_record(pa_nodat, "PCA")$y)
  expect_named(.retention_record(pa_nodat, "PCA")$references,
               c("Means", "95 Percentile"))
})


test_that("identified number of factors is correct", {
  skip_if_not_slow()
  expect_equal(pa_cor$n_factors[["PCA"]], 3)
  expect_equal(pa_cor$n_factors[["SMC"]], 3)
  expect_equal(pa_cor$n_factors[["EFA"]], 7, tolerance = 2)

  expect_equal(pa_raw$n_factors[["PCA"]], 1)
  expect_equal(pa_raw$n_factors[["SMC"]], 1)
  expect_equal(pa_raw$n_factors[["EFA"]], 3, tolerance = 2)

  expect_true(all(is.na(pa_nodat$n_factors)))

  expect_equal(pa_cor_pca$n_factors[["PCA"]], 3)
  expect_equal(pa_craw$n_factors[["PCA"]], 3)
  expect_equal(pa_perc$n_factors[["PCA"]], 3)
})


test_that("settings are returned correctly", {
  skip_if_not_slow()
  expect_named(pa_cor$settings, c("x_dat", "N", "n_vars", "n_datasets", "percent",
                                  "eigen_type", "use", "cor_method", "decision_rule",
                                  "n_factors"))
  expect_named(pa_raw$settings, c("x_dat", "N", "n_vars", "n_datasets", "percent",
                                  "eigen_type", "use", "cor_method", "decision_rule",
                                  "n_factors"))
  expect_named(pa_cor_pca$settings, c("x_dat", "N", "n_vars", "n_datasets", "percent",
                                      "eigen_type", "use", "cor_method",
                                      "decision_rule", "n_factors"))
  expect_named(pa_nodat$settings, c("x_dat", "N", "n_vars", "n_datasets", "percent",
                                      "eigen_type", "use", "cor_method",
                                      "decision_rule", "n_factors"))
  expect_named(pa_craw$settings, c("x_dat", "N", "n_vars", "n_datasets", "percent",
                                      "eigen_type", "use", "cor_method",
                                      "decision_rule", "n_factors"))
  expect_named(pa_perc$settings, c("x_dat", "N", "n_vars", "n_datasets", "percent",
                                      "eigen_type", "use", "cor_method",
                                      "decision_rule", "n_factors"))

  expect_true(pa_cor$settings$x_dat)
  expect_true(pa_raw$settings$x_dat)
  expect_true(pa_cor_pca$settings$x_dat)
  expect_false(pa_nodat$settings$x_dat)
  expect_true(pa_craw$settings$x_dat)
  expect_true(pa_perc$settings$x_dat)

  expect_equal(pa_cor$settings$N, 500)
  expect_equal(pa_raw$settings$N, 810)
  expect_equal(pa_cor_pca$settings$N, 500)
  expect_equal(pa_nodat$settings$N, 20)
  expect_equal(pa_craw$settings$N, 500)
  expect_equal(pa_perc$settings$N, 500)

  expect_equal(pa_cor$settings$n_vars, 18)
  expect_equal(pa_raw$settings$n_vars, 8)
  expect_equal(pa_cor_pca$settings$n_vars, 18)
  expect_equal(pa_nodat$settings$n_vars, 5)
  expect_equal(pa_craw$settings$n_vars, 18)
  expect_equal(pa_perc$settings$n_vars, 18)

  expect_equal(pa_cor$settings$n_datasets, 1000)
  expect_equal(pa_raw$settings$n_datasets, 1000)
  expect_equal(pa_cor_pca$settings$n_datasets, 1000)
  expect_equal(pa_nodat$settings$n_datasets, 1000)
  expect_equal(pa_craw$settings$n_datasets, 1000)
  expect_equal(pa_perc$settings$n_datasets, 1000)

  expect_equal(pa_cor$settings$percent, 95)
  expect_equal(pa_raw$settings$percent, 95)
  expect_equal(pa_cor_pca$settings$percent, 95)
  expect_equal(pa_nodat$settings$percent, 95)
  expect_equal(pa_craw$settings$percent, 95)
  expect_equal(pa_perc$settings$percent, 95)

  expect_equal(pa_cor$settings$eigen_type, c("PCA", "SMC", "EFA"))
  expect_equal(pa_raw$settings$eigen_type, c("PCA", "SMC", "EFA"))
  expect_equal(pa_cor_pca$settings$eigen_type, "PCA")
  expect_equal(pa_nodat$settings$eigen_type, c("PCA", "SMC", "EFA"))
  expect_equal(pa_craw$settings$eigen_type, "PCA")
  expect_equal(pa_perc$settings$eigen_type, "PCA")

  expect_equal(pa_cor$settings$use, "pairwise.complete.obs")
  expect_equal(pa_raw$settings$use, "pairwise.complete.obs")
  expect_equal(pa_cor_pca$settings$use, "pairwise.complete.obs")
  expect_equal(pa_nodat$settings$use, "pairwise.complete.obs")
  expect_equal(pa_craw$settings$use, "pairwise.complete.obs")
  expect_equal(pa_perc$settings$use, "pairwise.complete.obs")

  expect_equal(pa_cor$settings$cor_method, "pearson")
  expect_equal(pa_raw$settings$cor_method, "pearson")
  expect_equal(pa_cor_pca$settings$cor_method, "pearson")
  expect_equal(pa_nodat$settings$cor_method, "pearson")
  expect_equal(pa_craw$settings$cor_method, "pearson")
  expect_equal(pa_perc$settings$cor_method, "pearson")

  expect_equal(pa_cor$settings$decision_rule, "means")
  expect_equal(pa_raw$settings$decision_rule, "means")
  expect_equal(pa_cor_pca$settings$decision_rule, "means")
  expect_equal(pa_nodat$settings$decision_rule, "means")
  expect_equal(pa_craw$settings$decision_rule, "crawford")
  expect_equal(pa_perc$settings$decision_rule, "percentile")

  expect_equal(pa_cor$settings$n_factors, 1)
  expect_equal(pa_raw$settings$n_factors, 1)
  expect_equal(pa_cor_pca$settings$n_factors, 1)
  expect_equal(pa_nodat$settings$n_factors, 1)
  expect_equal(pa_craw$settings$n_factors, 1)
  expect_equal(pa_perc$settings$n_factors, 1)

})

burt <- .burt_cormat()


#sim_NA <- data.frame(rnorm(30), rnorm(30), rnorm(30), rep("a", 30))

test_that("errors are thrown correctly", {
  skip_if_not_slow()
  expect_error(efa_parallel(1:5), class = "efa_input_not_matrix")
  expect_warning(suppressMessages(efa_parallel(GRiPS_raw, n_vars = 5)), class = "efa_nvars_from_data")
  expect_warning(suppressMessages(efa_parallel(GRiPS_raw, N = 20, eigen_type = "PCA")), class = "efa_n_from_data")
  expect_error(suppressMessages(efa_parallel(N = 500)), class = "efa_nvars_required")
  expect_message(efa_parallel(GRiPS_raw, eigen_type = "PCA"), class = "efa_cor_from_data")
  expect_error(efa_parallel(test_models$baseline$cormat, eigen_type = "PCA"), class = "efa_n_required")
  expect_warning(efa_parallel(test_models$baseline$cormat, N = 500,
                          eigen_type = "PCA", decision_rule = "crawford",
                          percent = 80), class = "efa_parallel_crawford")
  expect_error(efa_parallel(sing_raw), class = "efa_cor_singular")
  expect_error(efa_parallel(sing_cor, N = sing_N), class = "efa_cor_singular")
  expect_warning(efa_parallel(burt, N = 100, eigen_type = "PCA"), class = "efa_cor_smoothed")
  expect_error(efa_parallel(test_models$baseline$cormat, N = 15), class = "efa_n_too_small")
  expect_error(efa_parallel(test_models$baseline$cormat, N = 18), class = "efa_n_too_small")
})

test_that("parallel-analysis simulation counts must be positive", {
  # all four are refused at the R boundary with the package's argument class,
  # rather than reaching the simulation and answering from an empty reference
  expect_error(efa_parallel(N = 0, n_vars = 4, n_datasets = 1),
               class = "efa_invalid_argument")
  expect_error(efa_parallel(N = 20, n_vars = 0, n_datasets = 1),
               class = "efa_invalid_argument")
  expect_error(efa_parallel(N = 20, n_vars = 4, n_datasets = 0),
               class = "efa_invalid_argument")
  expect_error(efa_parallel(N = 20, n_vars = 4, n_factors = 0),
               class = "efa_invalid_argument")
})

test_that("the simulated PCA reference matches psych::fa.parallel", {
  skip_if_not_slow()
  skip_if_not_installed("psych")
  # Both draw their reference eigenvalues from independent random normal data, so
  # they agree only up to Monte-Carlo error; 2000 datasets bring that well inside
  # the tolerance (observed max absolute difference ~0.0013).
  R <- test_models$baseline$cormat

  set.seed(42)
  ours <- efa_parallel(R, N = 500, n_datasets = 2000, eigen_type = "PCA")

  # fa.parallel reports its suggestion with cat(), so capture that output
  invisible(utils::capture.output(
    theirs <- suppressWarnings(psych::fa.parallel(R, n.obs = 500, fa = "pc",
                                                  n.iter = 2000, plot = FALSE))
  ))

  expect_equal(.retention_record(ours, "PCA")$references$Means, theirs$pc.sim,
               tolerance = 0.01)
  expect_equal(ours$n_factors[["PCA"]], theirs$ncomp)
})

test_that(".parallel_chunks splits exactly into non-negative integer chunks", {
  # the documented degenerate case: 11 datasets across 7 workers must not yield a
  # negative chunk (the old round-and-backfill produced c(2, 2, 2, 2, 2, 2, -1))
  chunks_11_7 <- .parallel_chunks(11, 7)
  expect_length(chunks_11_7, 7)
  expect_equal(sum(chunks_11_7), 11)
  expect_true(all(chunks_11_7 >= 0))

  # a range of awkward (n_datasets, n_cores) pairs: always exact, never negative
  for (n_cores in 1:8) {
    for (n_datasets in c(1, 3, 11, 100, 1000)) {
      chunks <- .parallel_chunks(n_datasets, n_cores)
      expect_length(chunks, n_cores)
      expect_equal(sum(chunks), n_datasets)
      expect_true(all(chunks >= 0))
    }
  }

  # a single chunk takes all datasets
  expect_equal(.parallel_chunks(1000, 1), 1000)
})


test_that("efa_parallel does not size its chunks from the worker count", {
  # The simulated datasets are split into chunks, one future per chunk, and
  # future.seed = TRUE assigns one random-number stream per chunk. If the chunk count
  # tracked the worker count, the streams -- and hence the reference eigenvalues -- would
  # differ between parallel plans for the same set.seed().
  #
  # Assert that through efa_parallel() itself rather than by re-computing the chunk
  # formula here: a test that only restates the arithmetic would still pass if the call
  # site were reverted to nbrOfWorkers(). Faking the reported worker count exercises the
  # call site while keeping the run sequential, so no workers are started and the check
  # runs on CRAN -- unlike the end-to-end test below.
  skip_if_not_installed("testthat", "3.2.0")

  run_reporting <- function(n_workers) {
    testthat::local_mocked_bindings(nbrOfWorkers = function(...) n_workers,
                                    .package = "future")
    set.seed(2024)
    suppressWarnings(suppressMessages(
      efa_parallel(N = 250, n_vars = 6, n_datasets = 24, eigen_type = "PCA")))
  }

  expect_equal(run_reporting(1L)$results[[1]]$references,
               run_reporting(4L)$results[[1]]$references)
})


test_that(".parallel_chunks is exact at the sizes efa_parallel asks for", {
  # fewer datasets than the fixed chunk count: one chunk each, never an empty chunk
  expect_equal(.parallel_chunks(5, max(1L, min(5, 20L))), rep(1, 5))

  # the degenerate n_datasets = 0 must not divide by a zero chunk count
  expect_equal(.parallel_chunks(0, max(1L, min(0, 20L))), 0)
})


test_that("only the simulation chunk that failed is redrawn", {
  # Deterministic chunk failure: the mocked simulation fails one chunk evaluation and only
  # that one, so the redraw meets a working simulation again. Counting the evaluations shows
  # what a failure costs -- redrawing the whole batch would repeat every chunk, which for a
  # failure in the last of twenty chunks is nineteen sound chunks thrown away.
  skip_if_not_installed("testthat", "3.2.0")

  # The mock replaces the simulation in this process only, so the chunks have to resolve
  # here; an ambient parallel plan would run the real one in a worker and count nothing.
  old_plan <- future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)

  real <- .parallel_sim_eig
  n_chunks <- length(.parallel_chunks(40, max(1L, min(40, 20L))))

  count <- function(fail_on) {
    calls <- 0L
    local_mocked_bindings(
      .parallel_sim_eig = function(n_datasets, ...) {
        calls <<- calls + 1L
        if (isTRUE(calls == fail_on)) {
          cli::cli_abort("simulated chunk failure", class = "efa_test_chunk_failure")
        }
        real(n_datasets, ...)
      }
    )
    set.seed(11)
    suppressWarnings(suppressMessages(
      efa_parallel(N = 200, n_vars = 5, n_datasets = 40, eigen_type = "PCA")))
    calls
  }

  expect_identical(count(NA_integer_), n_chunks)
  expect_identical(count(3L), n_chunks + 1L)
  expect_identical(count(n_chunks), n_chunks + 1L)
})


test_that("a redrawn chunk leaves the other chunks' draws untouched", {
  # The retry runs on the failing chunk alone, so the chunks that succeeded keep the exact
  # draws of the first attempt; only the redrawn one gets a fresh random-number stream.
  skip_if_not_installed("testthat", "3.2.0")

  old_plan <- future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)

  real <- .parallel_sim_eig
  size_vec <- .parallel_chunks(40, max(1L, min(40, 20L)))
  failed_rows <- (sum(size_vec[1:2]) + 1L):sum(size_vec[1:3])

  run <- function(fail_on) {
    calls <- 0L
    local_mocked_bindings(
      .parallel_sim_eig = function(n_datasets, ...) {
        calls <<- calls + 1L
        if (isTRUE(calls == fail_on)) {
          cli::cli_abort("simulated chunk failure", class = "efa_test_chunk_failure")
        }
        real(n_datasets, ...)
      }
    )
    set.seed(11)
    .parallel_sim_chunks(size_vec, label = "PCA", N = 200, n_vars = 5,
                         eigen_type = 1, cor_method = "pearson")
  }

  clean <- run(NA_integer_)
  broken <- run(3L)

  expect_identical(broken[-failed_rows, ], clean[-failed_rows, ])
  expect_false(identical(broken[failed_rows, ], clean[failed_rows, ]))
})


test_that("the parallel-analysis abort carries the failure that caused it", {
  skip_if_not_installed("testthat", "3.2.0")

  old_plan <- future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)

  local_mocked_bindings(
    .parallel_sim_eig = function(n_datasets, ...) {
      cli::cli_abort("no usable draw", class = "efa_test_sim_failure")
    }
  )

  cnd <- expect_error(
    .parallel_sim_chunks(.parallel_chunks(4, 2), label = "PCA", N = 200, n_vars = 5,
                         eigen_type = 1, cor_method = "pearson", max_tries = 3L),
    class = "efa_parallel_sim_failed")
  expect_s3_class(cnd$parent, "efa_test_sim_failure")
})


test_that("a failing parallel backend is retried and reported the same way", {
  # A chunk reports its own failure, so what the batch-level retry covers is the backend
  # itself -- a lost worker, say. It must reach the same classed abort, not leak out raw.
  skip_if_not_installed("testthat", "3.2.0")

  local_mocked_bindings(
    future_lapply = function(...) {
      cli::cli_abort("the workers are gone", class = "efa_test_backend_failure")
    },
    .package = "future.apply")

  cnd <- expect_error(
    .parallel_sim_chunks(.parallel_chunks(4, 2), label = "SMCs", N = 200, n_vars = 5,
                         eigen_type = 2, cor_method = "pearson", max_tries = 2L),
    class = "efa_parallel_sim_failed")
  expect_s3_class(cnd$parent, "efa_test_backend_failure")
})


test_that("a seeded efa_parallel run is invariant to the number of workers", {
  skip_on_cran()
  skip_if_not_slow()
  # End-to-end counterpart to the chunk-vector test above: the same set.seed() must give
  # the same reference eigenvalues sequentially and on a two-worker multisession plan.
  # The multisession workers are fresh R processes that load the *installed* package, so
  # run this under devtools::check() / after devtools::install() for the worker code to
  # match the main process (multicore is unavailable on Windows).
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  run <- function() {
    set.seed(2024)
    suppressWarnings(suppressMessages(
      efa_parallel(N = 300, n_vars = 8, n_datasets = 60, eigen_type = "PCA")))
  }

  future::plan(future::sequential)
  one <- run()

  future::plan(future::multisession, workers = 2)
  two <- run()

  # `results[[1]]$references` holds the simulated reference eigenvalues (means and
  # percentiles), which is what the chunk streams determine. Bit-for-bit equality is not
  # asserted: the workers are separate processes whose BLAS/LAPACK may sum in a different
  # order.
  expect_equal(one$results[[1]]$references, two$results[[1]]$references,
               tolerance = 1e-10)
  expect_equal(one$n_factors, two$n_factors)
})


test_that(".parallel_summarise uses stats::quantile for the percentile series", {
  set.seed(42)
  eig_vals <- matrix(stats::rnorm(1000 * 4), nrow = 1000, ncol = 4)

  res <- .parallel_summarise(eig_vals, percent = 95, n_vars = 4)
  expect_equal(dim(res), c(4, 2))
  expect_equal(res[, 1], colMeans(eig_vals))
  for (root in 1:4) {
    expect_equal(res[root, 2],
                 stats::quantile(eig_vals[, root], probs = 0.95, names = FALSE))
  }

  # multiple percentiles produce one column each
  res2 <- .parallel_summarise(eig_vals, percent = c(50, 95), n_vars = 4)
  expect_equal(dim(res2), c(4, 3))
  for (root in 1:4) {
    expect_equal(res2[root, 2:3],
                 stats::quantile(eig_vals[, root], probs = c(0.5, 0.95),
                                 names = FALSE))
  }

  # a NA/NaN simulated eigenvalue is tolerated, not turned into a hard error
  eig_na <- eig_vals
  eig_na[1, 1] <- NA
  expect_no_error(res_na <- .parallel_summarise(eig_na, percent = 95, n_vars = 4))
  expect_equal(res_na[1, 2],
               stats::quantile(eig_na[, 1], probs = 0.95, names = FALSE,
                               na.rm = TRUE))
})

test_that(".determine_factors retains all components and warns when no crossing", {
  # every real eigenvalue exceeds its reference: the rule finds no crossing, so
  # instead of a silent NA it retains all tested components and flags the boundary
  eigvals_real <- matrix(c(3, 2, 1.5), ncol = 1)
  results <- cbind(c(2, 1, 0.5), c(2.5, 1.2, 0.6))
  colnames(results) <- c("Means", "95 Percentile")

  for (rule in c("means", "percentile", "crawford")) {
    expect_warning(
      n_fac <- .determine_factors(rule, eigvals_real, results, percent = 95),
      class = "efa_parallel_no_crossing"
    )
    expect_equal(n_fac, nrow(eigvals_real), info = rule)
  }
})

test_that(".determine_factors returns a normal crossing and 0 without the warning", {
  results <- cbind(c(1, 1, 1), c(1.2, 1.2, 1.2))
  colnames(results) <- c("Means", "95 Percentile")

  # first eigenvalue exceeds the reference, the second does not -> retain 1
  expect_no_warning(
    n1 <- .determine_factors("means", matrix(c(3, 0.5, 0.2), ncol = 1),
                             results, percent = 95),
    class = "efa_parallel_no_crossing"
  )
  expect_equal(n1, 1)

  # the first eigenvalue is already at/below the reference -> retain 0
  expect_no_warning(
    n0 <- .determine_factors("means", matrix(c(0.5, 0.4, 0.3), ncol = 1),
                             results, percent = 95),
    class = "efa_parallel_no_crossing"
  )
  expect_equal(n0, 0)
})

test_that("the null-model reference draw matches the compiled fast path", {
  # The rank-method and EFA paths draw the null model with the shared MVN kernel on the
  # identity correlation, while the Pearson PCA/SMC path uses the compiled .parallel_sim().
  # This pins that the two draw the same underlying stream: the kernel's eigenvalues agree
  # with the compiled fast path up to the compiled cor() vs stats::cor() rounding (~1e-15)
  # under a fixed seed. (That PARALLEL's rank/EFA paths actually route through the kernel is
  # guarded by the two tests below, which drive .parallel_sim_eig() and .parallel_EFA_sim().)
  N <- 60; p <- 5; nd <- 3

  set.seed(101)
  e_fast <- .parallel_sim(nd, p, N, 1L)          # compiled PCA fast path

  set.seed(101)
  e_kernel <- t(vapply(seq_len(nd), function(i) {
    eigen(stats::cor(.simulate_cfm_mvn(diag(p), N)), symmetric = TRUE,
          only.values = TRUE)$values
  }, numeric(p)))

  expect_equal(e_kernel, e_fast, tolerance = 1e-10)
})

test_that(".parallel_sim_eig draws its rank-method reference from the shared kernel", {
  # The rank-method PCA path never rejects a draw, so its reference eigenvalues must match a
  # direct per-draw recompute through the shared kernel under the same random-number stream
  # -- a guard that the draw routes through .simulate_cfm_mvn().
  N <- 60; p <- 5; nd <- 4

  set.seed(202)
  got <- .parallel_sim_eig(nd, n_vars = p, N = N, eigen_type = 1L,
                           cor_method = "spearman")

  set.seed(202)
  want <- t(vapply(seq_len(nd), function(i) {
    eigen(stats::cor(.simulate_cfm_mvn(diag(p), N), method = "spearman"),
          symmetric = TRUE, only.values = TRUE)$values
  }, numeric(p)))

  expect_equal(got, want, tolerance = fp_tol)
})

test_that("the PCA and SMC series are taken from one set of simulated datasets", {
  # PCA and SMC differ only in the diagonal put into the same simulated correlation matrix,
  # so asking for both must simulate n_datasets datasets, not 2 * n_datasets. Counted at the
  # kernel rather than by restating the branch, so a call site that reverted to two
  # independent simulations would fail here.
  skip_if_not_installed("testthat", "3.2.0")

  # The mock counts draws in this process, so the chunks have to resolve here.
  old_plan <- future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)

  real <- .parallel_sim
  n_datasets <- 20
  n_chunks <- length(.parallel_chunks(n_datasets, max(1L, min(n_datasets, 20L))))

  traced <- function(eigen_type) {
    drawn <- 0L
    types <- numeric(0)
    local_mocked_bindings(
      .parallel_sim = function(n_datasets, n_vars, N, eigen_type, maxit) {
        drawn <<- drawn + n_datasets
        types <<- c(types, eigen_type)
        real(n_datasets, n_vars, N, eigen_type, maxit)
      })
    set.seed(3)
    suppressWarnings(suppressMessages(
      efa_parallel(N = 200, n_vars = 5, n_datasets = n_datasets,
                   eigen_type = eigen_type)))
    list(drawn = drawn, types = types)
  }

  pca <- traced("PCA")
  expect_identical(pca$drawn, n_datasets)
  expect_equal(pca$types, rep(1, n_chunks))

  smc <- traced("SMC")
  expect_identical(smc$drawn, n_datasets)
  expect_equal(smc$types, rep(2, n_chunks))

  both <- traced(c("PCA", "SMC"))
  expect_identical(both$drawn, n_datasets)
  expect_equal(both$types, rep(3, n_chunks))

  # the default set adds "EFA", which fits a model per dataset and draws its own
  all_three <- traced(c("PCA", "SMC", "EFA"))
  expect_identical(all_three$drawn, n_datasets)
  expect_equal(all_three$types, rep(3, n_chunks))
})

test_that("the shared draw feeds both series with the eigenvalues each would get alone", {
  # The shared kernel run must return exactly the two series the single-type runs produce
  # from the same random-number stream: the PCA eigenvalues in the first n_vars columns and
  # the SMC eigenvalues in the next, replicate by replicate.
  N <- 200L; p <- 6L; nd <- 25L

  set.seed(4)
  shared <- .parallel_sim(nd, p, N, 3L, nd * 10L)
  set.seed(4)
  pca_only <- .parallel_sim(nd, p, N, 1L, nd * 10L)
  set.seed(4)
  smc_only <- .parallel_sim(nd, p, N, 2L, nd * 10L)

  expect_identical(dim(shared), c(nd, 2L * p))
  expect_identical(shared[, seq_len(p)], pca_only)
  expect_identical(shared[, p + seq_len(p)], smc_only)
})

test_that("a draw the SMC series cannot use is discarded for the PCA series too", {
  # A simulated correlation matrix with no inverse yields no squared multiple correlations,
  # and the draw is then dropped from both series so that they stay paired. Admissible input
  # essentially never produces such a matrix, so the rejection is forced here, on the
  # rank-based route where it is taken in R.
  skip_if_not_installed("testthat", "3.2.0")

  N <- 60; p <- 5
  real <- .smc_start

  set.seed(606)
  clean <- .parallel_sim_eig(4, n_vars = p, N = N, eigen_type = 3L,
                             cor_method = "spearman")

  calls <- 0L
  local_mocked_bindings(
    .smc_start = function(R, ...) {
      calls <<- calls + 1L
      if (calls == 2L) stop("singular simulated matrix")
      real(R, ...)
    })
  set.seed(606)
  rejected <- .parallel_sim_eig(3, n_vars = p, N = N, eigen_type = 3L,
                                cor_method = "spearman")

  # the second draw is gone from both halves, and the third took its place in both
  expect_identical(rejected[1L, ], clean[1L, ])
  expect_identical(rejected[2L, ], clean[3L, ])
  expect_identical(rejected[3L, ], clean[4L, ])
})

test_that(".parallel_EFA_sim draws its reference from the shared kernel", {
  # The EFA path likewise draws the null-model data with the shared kernel; under a
  # fixed seed its eigenvalues equal a direct recompute that draws the same way.
  N <- 60; p <- 5; nd <- 2

  set.seed(303)
  got <- .parallel_EFA_sim(nd, n_vars = p, N = N, n_factors = 1,
                           cor_method = "pearson")

  set.seed(303)
  want <- t(vapply(seq_len(nd), function(i) {
    R <- stats::cor(.simulate_cfm_mvn(diag(p), N))
    suppressWarnings(suppressMessages(EFA(R, n_factors = 1, N = N)$final_eigen))
  }, numeric(p)))

  expect_equal(got, want)
})

if (is_slow_test()) {
  rm(pa_cor, pa_cor_pca, pa_raw, pa_nodat, pa_craw, pa_perc)
}
rm(burt)
