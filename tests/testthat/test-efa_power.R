# Regression values are the RMSEA power-analysis figures tabulated by MacCallum,
# Browne & Sugawara (1996) for a 100-df model, reproduced with base R's noncentral
# chi-square functions (no external dependency).

test_that("RMSEA close-fit power reproduces MacCallum, Browne & Sugawara (1996)", {
  cf <- efa_power(df = 100, N = 200, type = "close", eps0 = 0.05, eps1 = 0.08)
  expect_equal(cf$power, 0.9549, tolerance = 1e-4)
  expect_equal(cf$solve_for, "power")
  expect_equal(cf$N, 200L)
  # ncp = (N - 1) * df * eps^2: 199 * 100 * .0025 and 199 * 100 * .0064
  expect_equal(unname(cf$ncp), c(49.75, 127.36))
  expect_named(cf$ncp, c("H0", "H1"))
})

test_that("RMSEA not-close-fit power reproduces the tabled value", {
  nc <- efa_power(df = 100, N = 200, type = "notclose", eps0 = 0.05, eps1 = 0.01)
  expect_equal(nc$power, 0.8701, tolerance = 1e-4)
  expect_equal(nc$solve_for, "power")
})

test_that("required sample size is the smallest N reaching the target power", {
  rn <- efa_power(df = 100, power = 0.80, type = "close")
  expect_equal(rn$N, 132L)
  expect_equal(rn$solve_for, "N")

  # 132 is minimal: one fewer falls short of .80, 132 clears it.
  expect_lt(efa_power(df = 100, N = 131)$power, 0.80)
  expect_gte(efa_power(df = 100, N = 132)$power, 0.80)

  # With neither N nor power supplied the target defaults to .80.
  expect_equal(efa_power(df = 100)$N, 132L)
  expect_equal(efa_power(df = 100)$settings$power, 0.80)
})

test_that("df is taken from n_vars and n_factors via the EFA df formula", {
  # n_vars = 20, n_factors = 5 -> ((20 - 5)^2 - (20 + 5)) / 2 = 100
  from_dims <- efa_power(n_vars = 20, n_factors = 5, N = 200)
  from_df <- efa_power(df = 100, N = 200)
  expect_equal(from_dims$settings$df, 100)
  expect_equal(from_dims$power, from_df$power)
})

test_that("p and k are deprecated aliases of n_vars and n_factors", {
  new <- efa_power(n_vars = 20, n_factors = 5, N = 200)
  expect_warning(old_p <- efa_power(p = 20, n_factors = 5, N = 200),
                 class = "lifecycle_warning_deprecated")
  expect_warning(old_k <- efa_power(n_vars = 20, k = 5, N = 200),
                 class = "lifecycle_warning_deprecated")
  expect_equal(old_p, new)
  expect_equal(old_k, new)
})

test_that("the noncentrality uses Steiger's (N - group) / group factor", {
  # Steiger (1998), eqs. 15-23: n. = N - G, RMSEA = sqrt(G) * sqrt(F / df), so
  # ncp = (N - G) * df * eps^2 / G. At N = 200, G = 2: 99 * 100 * eps^2.
  g2 <- efa_power(df = 100, N = 200, group = 2)
  expect_equal(unname(g2$ncp), c(99 * 100 * 0.05^2, 99 * 100 * 0.08^2))
})

test_that("n_vars and n_factors give the df of a separate EFA in each group", {
  expect_equal(efa_power(n_vars = 18, n_factors = 3, N = 200, group = 2)$settings$df, 204)
  # The joint test sums the two groups' chi-squares, so it needs fewer cases per group.
  two <- efa_power(n_vars = 18, n_factors = 3, group = 2)
  expect_lt(two$N_per_group, efa_power(n_vars = 18, n_factors = 3)$N)
})

test_that("an explicit NULL of the other mode is not an error", {
  expect_s3_class(efa_power(df = 100, N = 200, Lambda = NULL, seed = NULL), "efa_power")
})

test_that("with equal groups, power depends only on the group size", {
  # (N - G) / G = m - 1 for G groups of m cases: the single-group noncentrality at m.
  for (g in 2:3) {
    expect_equal(efa_power(df = 102, N = 130 * g, group = g)$power,
                 efa_power(df = 102, N = 130)$power)
  }
  g1 <- efa_power(df = 102, power = 0.80)
  expect_equal(g1$N, 130L)
  for (g in 1:4) {
    x <- efa_power(df = 102, power = 0.80, group = g)
    expect_equal(x$N, 130L * g)
    expect_equal(x$N_per_group, 130)
  }

  # A target below `alpha` is reached at any size: one case per group.
  floor_g2 <- efa_power(df = 102, power = 0.01, group = 2)
  expect_equal(floor_g2$N, 2L)
  expect_equal(floor_g2$N_per_group, 1)

  # A given N is the caller's own total and is kept unchanged.
  odd <- efa_power(df = 102, N = 259L, group = 2)
  expect_equal(odd$N, 259L)
  expect_equal(odd$N_per_group, 129.5)
})

test_that("the single-group solved N matches semTools::findRMSEAsamplesize", {
  skip_if_not_installed("semTools")

  # semTools uses (N - 1) / group for every group count, which agrees with Steiger's
  # (N - group) / group only for one group; group > 1 is pinned by the invariant above.
  # findRMSEAsamplesize() reports the first N with power strictly above the target.
  ref <- semTools::findRMSEAsamplesize(rmsea0 = 0.05, rmseaA = 0.08, df = 102,
                                       power = 0.80, group = 1)
  expect_equal(efa_power(df = 102, power = 0.80)$N, ref)
})

test_that("output has the expected class and structure", {
  x <- efa_power(df = 100, N = 200)
  expect_s3_class(x, "efa_power")
  expect_named(x, c("power", "N", "N_per_group", "crit", "ncp", "solve_for",
                    "settings"))
  expect_named(x$settings, c("mode", "type", "eps0", "eps1", "df", "n_vars",
                             "n_factors", "alpha", "group", "power"))
})

test_that("classed conditions fire as expected", {
  # Equal null and alternative RMSEA: nothing to detect.
  expect_error(efa_power(df = 100, N = 200, eps0 = 0.05, eps1 = 0.05),
               class = "efa_power_unreachable")
  # eps0/eps1 on the wrong side for the chosen test: the result would not be a power.
  expect_error(efa_power(df = 100, N = 200, type = "close", eps0 = 0.08, eps1 = 0.05),
               class = "efa_power_wrong_side")
  expect_error(efa_power(df = 100, N = 200, type = "notclose", eps0 = 0.05, eps1 = 0.08),
               class = "efa_power_wrong_side")
  # Non-positive degrees of freedom, given directly or from an unidentified model.
  expect_error(efa_power(df = 0, N = 200), class = "efa_power_bad_df")
  expect_error(efa_power(n_vars = 3, n_factors = 3, N = 200), class = "efa_power_bad_df")
  # No way to obtain df.
  expect_error(efa_power(N = 200), class = "efa_power_missing_df")
  # Both N and power supplied: nothing to solve for.
  expect_error(efa_power(df = 100, N = 200, power = 0.80),
               class = "efa_power_overdetermined")
  # Out-of-range alpha and target power.
  expect_error(efa_power(df = 100, N = 200, alpha = 0),
               class = "efa_power_bad_alpha")
  expect_error(efa_power(df = 100, power = 1.2), class = "efa_power_bad_power")
  # Fewer cases than groups.
  expect_error(efa_power(df = 100, N = 2, group = 3), class = "efa_power_bad_n")
})

test_that("a sample size not larger than n_vars gives a warning", {
  # With 1480 df the test needs only about 30 cases, fewer than the 60 variables
  # (MacCallum et al., 1996, p. 144).
  expect_warning(small <- efa_power(n_vars = 60, n_factors = 5),
                 class = "efa_power_small_n")
  expect_lte(small$N, 60L)
  expect_warning(efa_power(n_vars = 60, n_factors = 5, N = 60),
                 class = "efa_power_small_n")
  # The check is per group.
  expect_warning(efa_power(n_vars = 60, n_factors = 5, N = 100, group = 2),
                 class = "efa_power_small_n")
  expect_no_warning(efa_power(n_vars = 60, n_factors = 5, N = 61))
})

test_that("an argument of the other mode is an error, and a population selects simulation", {
  expect_error(efa_power(df = 100, N = 200, mode = "rmsea", Lambda = diag(3)),
               class = "efa_power_wrong_mode")
  expect_error(efa_power(df = 100, N = 200, seed = 1), class = "efa_power_wrong_mode")
  expect_error(efa_power("simulation", Lambda = diag(3), N = 200, alpha = 0.01),
               class = "efa_power_wrong_mode")
  # Lambda without `mode` selects simulation mode, so `df` is then the stray argument.
  expect_error(efa_power(df = 100, N = 200, Lambda = diag(3)),
               class = "efa_power_wrong_mode")
})

test_that("power stays correct at a very large noncentrality", {
  # stats::qchisq() fails above a noncentrality of about 1e5 and stats::pchisq() above
  # about 1e7; Pearson's (1959) three-moment approximation takes over there. An
  # alternative this near the null has power near alpha at any N.
  expect_equal(efa_power(df = 100, N = 1e7, eps1 = 0.0500001)$power, 0.05,
               tolerance = 0.01)
  # Above 1e5 the approximation agrees with the exact values wherever those are still
  # accurate: pchisq() up to about 1e7, qchisq() at the switch.
  for (ncp in c(1e5 + 1, 1e6)) {
    q <- 100 + ncp + c(-2, 2) * sqrt(2 * (100 + 2 * ncp))
    # Absolute errors of about 1e-7 at tail probabilities of .02 to .98 (relative tolerance).
    expect_equal(.pchisq_nc(q, 100, ncp, TRUE), stats::pchisq(q, 100, ncp), tolerance = 1e-5)
    expect_equal(stats::pchisq(.qchisq_nc(0.05, 100, ncp, FALSE), 100, ncp, lower.tail = FALSE),
                 0.05, tolerance = 1e-5)
  }
  expect_equal(.qchisq_nc(0.05, 100, 1e5 + 1, FALSE),
               suppressWarnings(stats::qchisq(0.05, 100, 1e5 + 1, lower.tail = FALSE)),
               tolerance = 1e-7)
})

test_that("an unreachable target aborts cleanly without warning noise", {
  # An alternative this near the null cannot reach .80 within the search bound.
  expect_error(efa_power(df = 100, power = 0.80, eps0 = 0.05, eps1 = 0.0500001),
               class = "efa_power_unreached")

  # The search evaluates the noncentral chi-square routines at extreme
  # noncentrality; those base-R precision/convergence notes must not leak.
  expect_no_warning(tryCatch(
    efa_power(df = 100, power = 0.80, eps0 = 0.05, eps1 = 0.0500001),
    error = function(e) NULL))
})

test_that("printed output is stable", {
  local_reproducible_output()

  # Solve for power (test of close fit).
  expect_snapshot(print(efa_power(df = 100, N = 200)))
  # Solve for the required sample size.
  expect_snapshot(print(efa_power(df = 100, power = 0.80)))
  # Test of not-close fit with more than one group.
  expect_snapshot(print(efa_power(df = 100, N = 200, type = "notclose", group = 2)))
})

test_that("plot() returns a ggplot for the various curves", {
  # Base curve (solved for power), and one solved for the required N.
  expect_s3_class(plot(efa_power(df = 100, N = 200)), "ggplot")
  expect_s3_class(plot(efa_power(df = 100, power = 0.80)), "ggplot")
  # Sweep df, sweep eps1, and a caller-supplied sample-size sequence.
  expect_s3_class(plot(efa_power(df = 100, N = 200), df = c(50, 100, 200)), "ggplot")
  expect_s3_class(plot(efa_power(df = 100, N = 200), eps1 = c(0.06, 0.08, 0.10)),
                  "ggplot")
  expect_s3_class(plot(efa_power(df = 100, N = 200), n = seq(50, 400, by = 50)),
                  "ggplot")
  # Not-close fit with more than one group still plots, with finite power on every point.
  pg <- plot(efa_power(df = 100, N = 200, type = "notclose", group = 2))
  expect_s3_class(pg, "ggplot")
  expect_true(all(is.finite(pg$data$power)))

  # An alternative near the null cannot reach 0.99 power within the search bound, so
  # the N axis falls back to a multiple of the object's N rather than erroring out.
  near <- efa_power(df = 100, N = 200, eps1 = 0.0500001)
  expect_s3_class(plot(near), "ggplot")
})

test_that("the object's point is marked only for its own curve", {
  pw <- efa_power(df = 100, N = 200)
  # Own curve: reference lines and the point are drawn (hline, vline, line, point).
  expect_length(plot(pw)$layers, 4L)
  # A scalar df/eps1 override redraws the curve, so the object's point no longer lies
  # on it and the marks are dropped -- only the line remains.
  expect_length(plot(pw, df = 50)$layers, 1L)
  expect_length(plot(pw, eps1 = 0.10)$layers, 1L)
  # A supplied N range that excludes the object's N would place the point off-axis, so
  # the marks are dropped there too.
  expect_length(plot(pw, n = seq(10, 50, by = 5))$layers, 1L)
})

test_that("the subtitle comparator follows the test direction", {
  # Close fit bounds the null from above, not-close fit from below.
  expect_match(plot(efa_power(df = 100, N = 200))$labels$subtitle, "H0 RMSEA \u2264")
  expect_match(plot(efa_power(df = 100, N = 200, type = "notclose"))$labels$subtitle,
               "H0 RMSEA \u2265")
})

test_that("overrides uphold the constructor's invariants", {
  pw <- efa_power(df = 100, N = 200)
  # A non-positive df has no identified model, exactly as efa_power() rejects it.
  expect_error(plot(pw, df = 0), class = "efa_power_bad_df")
  expect_error(plot(pw, df = c(0, 100)), class = "efa_power_bad_df")
  # An alternative RMSEA equal to the null leaves nothing to detect.
  expect_error(plot(pw, eps1 = 0.05), class = "efa_power_unreachable")
  # An alternative on the null's side of a close-fit test is not a power curve.
  expect_error(plot(pw, eps1 = c(0.03, 0.08)), class = "efa_power_wrong_side")
  expect_error(plot(efa_power(df = 100, N = 200, type = "notclose"), eps1 = 0.08),
               class = "efa_power_wrong_side")
})

test_that("sweeping both df and eps1 is refused", {
  expect_error(
    plot(efa_power(df = 100, N = 200), df = c(50, 100), eps1 = c(0.06, 0.08)),
    class = "efa_power_plot_grid")
})

test_that("the power curve is visually stable", {
  skip_if_not_installed("vdiffr")

  vdiffr::expect_doppelganger("efa_power RMSEA power curve",
                              plot(efa_power(df = 100, N = 200)))
  vdiffr::expect_doppelganger("efa_power power curve swept df",
                              plot(efa_power(df = 100, N = 200), df = c(50, 100, 200)))
})


# --- Simulation mode -----------------------------------------------------------
# A clean, strong three-factor population fixture reused across the simulation-mode
# tests: the retention criteria and the loadings recovery should both succeed on it.
sim_Lambda <- population_models$loadings$baseline
sim_Phi <- population_models$phis_3$moderate

test_that("simulation mode is reproducible and leaves the RNG stream unchanged", {
  a <- efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 150,
                 n_datasets = 8, criteria = "EKC", seed = 42)
  b <- efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 150,
                 n_datasets = 8, criteria = "EKC", seed = 42)
  expect_equal(a, b)

  # A fixed seed restores the caller's random-number stream (no global side effect).
  set.seed(99)
  before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 150,
            n_datasets = 4, criteria = "EKC", seed = 7)
  after <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  expect_identical(before, after)
})

test_that("simulation mode is reproducible at 1 vs 2 workers", {
  skip_on_cran()
  skip_if_not_slow()
  # Each replicate is analysed in its own future, and future.seed = TRUE binds its RNG
  # stream (data-independent criteria and the rotation's random starts alike) to its
  # index, so the result must be identical regardless of the number of workers. The
  # multisession workers are fresh R processes that load the installed package, so run
  # this under devtools::check() / after devtools::install().
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  run <- function() {
    efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 150,
              n_datasets = 6, criteria = "EKC", seed = 2024)
  }

  future::plan(future::sequential)
  one <- run()
  future::plan(future::multisession, workers = 2)
  two <- run()

  expect_equal(one$hit_rate, two$hit_rate)
  expect_equal(one$recovery, two$recovery)
  expect_equal(one$convergence, two$convergence)
  expect_equal(one$replicates, two$replicates)
})

test_that("simulation mode recovers a clean population", {
  sim <- efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 300,
                   n_datasets = 20, criteria = c("EKC", "MAP"), seed = 42)
  expect_s3_class(sim, "efa_power")
  expect_identical(sim$settings$mode, "simulation")
  expect_identical(sim$settings$rotation, "oblimin")
  expect_equal(sim$k_true, 3L)
  expect_named(sim, c("hit_rate", "hits", "recovery", "convergence", "replicates",
                      "k_true", "model_error", "settings"))

  # EKC and the loadings recovery should both be near-perfect on a clean 3-factor model.
  expect_gte(sim$hit_rate[["EKC_BvA2017"]], 0.9)
  expect_gte(sim$recovery$min_rate, 0.9)
  expect_equal(sim$recovery$threshold, 0.95)
  expect_equal(sim$convergence$convergence_rate, 1)
  expect_equal(sim$convergence$heywood_rate, 0)

  # Every fit completed, so no replicate carries a failure reason.
  expect_true(all(sim$replicates$fit_ok))
  expect_true(all(is.na(sim$replicates$fit_error)))

  # The hits data frame agrees with the hit_rate vector, over all replicates.
  expect_equal(sim$hits$hit_rate, unname(sim$hit_rate))
  expect_equal(sim$hits$hit_rate, sim$hits$hits / 20)
  expect_equal(sim$hits$mcse, sqrt(sim$hits$hit_rate * (1 - sim$hits$hit_rate) / 20))
})

test_that("a population without `mode` selects simulation mode", {
  sim <- efa_power(Lambda = sim_Lambda, Phi = sim_Phi, N = 120, n_datasets = 2,
                   criteria = "EKC", seed = 1)
  expect_identical(sim$settings$mode, "simulation")
})

test_that("simulation mode muffles the per-replicate fit messages", {
  expect_no_message(
    efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 120,
              n_datasets = 3, criteria = "EKC", seed = 1))
})

test_that("a criterion that fails on every replicate is reported with a hit-rate of 0", {
  # The criterion must still appear in the results -- with no valid replicates and a
  # hit-rate of 0, because a failure is a miss -- and be named in a classed warning.
  testthat::local_mocked_bindings(
    efa_map = function(...) cli::cli_abort("Mocked MAP failure.", class = "efa_mock_map"))
  expect_warning(
    sim <- efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 120,
                     n_datasets = 3, criteria = "MAP", seed = 1),
    class = "efa_power_criterion_failed")
  expect_s3_class(sim, "efa_power")
  expect_named(sim$hit_rate, "MAP")
  expect_equal(sim$hit_rate[["MAP"]], 0)
  expect_equal(sim$hits$n_valid, 0, ignore_attr = TRUE)
  expect_equal(sim$hits$hits, 0, ignore_attr = TRUE)
  expect_identical(sim$settings$criteria, "MAP")
  local_reproducible_output()
  expect_snapshot(print(sim), transform = scrub_num)
})

test_that("a criterion that fails on some replicates counts them as misses", {
  # MAP fails on the replicates whose first sample correlation is above its population value.
  real_map <- efa_map
  pop_r12 <- (sim_Lambda %*% sim_Phi %*% t(sim_Lambda))[1, 2]
  testthat::local_mocked_bindings(
    efa_map = function(x, ...) {
      if (x[1, 2] > pop_r12) cli::cli_abort("Mocked MAP failure.", class = "efa_mock_map")
      real_map(x, ...)
    })
  sim <- efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 150,
                   n_datasets = 8, criteria = "MAP", seed = 3)
  n_valid <- sim$hits$n_valid
  expect_true(all(n_valid > 0 & n_valid < 8))
  expect_equal(sim$hits$hit_rate, sim$hits$hits / 8)
  expect_true(all(sim$hits$hit_rate <= n_valid / 8))
})

test_that("a criterion that decides on every replicate raises no failure warning", {
  expect_no_warning(
    efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 150,
              n_datasets = 4, criteria = "EKC", seed = 3))
})

test_that("a recovery fit that fails everywhere is reported, not thrown away", {
  # The recovery fit fails on every replicate while the retention criteria are
  # unaffected: the completed retention results must still be returned, with a warning.
  # The mock carries a `rotation` formal, which efa_power() validates against.
  testthat::local_mocked_bindings(
    efa_fit = function(x, ..., rotation = "oblimin") {
      cli::cli_abort("Mocked recovery failure.", class = "efa_mock_fit_failure")
    })

  expect_warning(
    sim <- efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 150,
                     n_datasets = 3, criteria = "EKC", seed = 1),
    class = "efa_power_fit_failed")
  expect_s3_class(sim, "efa_power")
  expect_equal(sim$hits$n_valid, 3, ignore_attr = TRUE)
  expect_true(is.finite(sim$hit_rate[["EKC_BvA2017"]]))

  # No replicate is recovered, and every replicate keeps the reason its fit failed.
  expect_equal(sim$convergence$fit_rate, 0)
  expect_true(is.na(sim$convergence$convergence_rate))
  expect_equal(sim$recovery$min_rate, 0)
  expect_equal(sim$recovery$n_valid, 0L)
  expect_false(any(sim$replicates$fit_ok))
  expect_length(sim$replicates$fit_error, 3L)
  expect_true(all(grepl("Mocked recovery failure", sim$replicates$fit_error,
                        fixed = TRUE)))
})

test_that("a partly failing recovery fit keeps the replicates that did fit", {
  # A deterministic per-replicate property of the drawn data decides which fits fail,
  # so the run has both kinds of replicate under a fixed seed.
  real_fit <- efa_fit
  testthat::local_mocked_bindings(
    efa_fit = function(x, ..., rotation = "oblimin") {
      if (x[1, 1] > 0) {
        cli::cli_abort("Mocked recovery failure.", class = "efa_mock_fit_failure")
      }
      real_fit(x, ..., rotation = rotation)
    })

  expect_warning(
    sim <- efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 150,
                     n_datasets = 6, criteria = "EKC", seed = 11),
    class = "efa_power_fit_failed")
  failed <- !sim$replicates$fit_ok
  expect_true(any(failed))
  expect_false(all(failed))

  # Convergence is over the completed fits; recovery is over all replicates, with a
  # failed fit as not recovered.
  expect_equal(sim$convergence$n_fit_ok, sum(!failed))
  expect_equal(sim$convergence$fit_rate, mean(!failed))
  expect_equal(sim$recovery$n_valid, sum(!failed))
  expect_equal(sim$recovery$min_rate,
               sum(sim$replicates$rec_min >= 0.95, na.rm = TRUE) / 6)
  expect_lte(sim$recovery$min_rate, mean(!failed))
  expect_true(all(grepl("Mocked recovery failure",
                        sim$replicates$fit_error[failed], fixed = TRUE)))
  expect_true(all(is.na(sim$replicates$fit_error[!failed])))
})

test_that("an unusable recovery rotation aborts before any data is drawn", {
  # `rotation` is resolved against efa_fit()'s own choices, so a value it would reject
  # costs no simulation time.
  expect_error(
    efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 200,
              n_datasets = 500, criteria = "EKC", rotation = "nonesuch"),
    class = "efa_bad_choice")
  expect_error(
    efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 200,
              n_datasets = 500, criteria = "EKC",
              rotation = c("varimax", "promax")),
    class = "efa_bad_choice")
  # A valid rotation is matched case-insensitively and stored canonically, as elsewhere.
  vm <- efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 120,
                  n_datasets = 2, criteria = "EKC", rotation = "Varimax", seed = 1)
  expect_identical(vm$settings$rotation, "varimax")
  expect_true(all(is.na(vm$replicates$fit_error)))
})

test_that("simulation mode: R-only population turns recovery off but keeps the hit-rate", {
  R_pop <- efa_simulate(Lambda = sim_Lambda, Phi = sim_Phi, return_pop = TRUE)$population
  rr <- efa_power("simulation", R = R_pop, n_factors = 3, N = 200, n_datasets = 6,
                  criteria = "EKC", seed = 42)
  expect_null(rr$recovery)
  expect_equal(rr$k_true, 3L)
  expect_true(is.finite(rr$hit_rate[["EKC_BvA2017"]]))
  # Neither the exact-population note nor a recovery section applies to R.
  local_reproducible_output()
  out <- format(rr)
  expect_false(any(grepl("Model error|Structure recovery", out)))
})

test_that("model error follows each function's own documented default", {
  # The two defaults deliberately differ: efa_power() perturbs with TKL, whose minor
  # common factors degrade structure recovery as well as fit, while efa_simulate()
  # defaults to CB. The same target passed to both therefore builds different
  # populations unless `model_error` is set explicitly.
  pw <- efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 150,
                  n_datasets = 2, criteria = "EKC", target_rmsea = 0.05, seed = 1)
  expect_identical(pw$settings$model_error, "TKL")
  expect_identical(pw$model_error$method, "TKL")
  expect_equal(pw$model_error$target_rmsea, 0.05)
  expect_lt(abs(pw$model_error$rmsea - 0.05), 0.02)

  sim <- efa_simulate(Lambda = sim_Lambda, Phi = sim_Phi, target_rmsea = 0.05,
                      return_pop = TRUE, seed = 1)
  expect_identical(sim$model_error$method, "CB")

  # Without a target neither function perturbs the population, whatever the method.
  exact <- efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 120,
                     n_datasets = 2, criteria = "EKC", seed = 1)
  expect_null(exact$model_error)
})

test_that("simulation-mode classed conditions fire as expected", {
  R_pop <- efa_simulate(Lambda = sim_Lambda, Phi = sim_Phi, return_pop = TRUE)$population

  # The population must be given exactly one way.
  expect_error(efa_power("simulation", N = 200, n_datasets = 2),
               class = "efa_power_input")
  expect_error(efa_power("simulation", Lambda = sim_Lambda, R = R_pop, N = 200,
                         n_datasets = 2), class = "efa_power_input")
  # n_factors is required with a bare R, and must agree with Lambda otherwise.
  expect_error(efa_power("simulation", R = R_pop, N = 200, n_datasets = 2),
               class = "efa_power_missing_k")
  expect_error(efa_power("simulation", Lambda = sim_Lambda, n_factors = 2, N = 200,
                         n_datasets = 2), class = "efa_power_bad_k")
  # n_vars is read off the population: a conflicting value is an error on both
  # population routes rather than being silently replaced.
  expect_error(efa_power("simulation", Lambda = sim_Lambda, n_vars = 99, N = 200,
                         n_datasets = 2), class = "efa_power_bad_p")
  expect_error(efa_power("simulation", R = R_pop, n_factors = 3, n_vars = 99, N = 200,
                         n_datasets = 2), class = "efa_power_bad_p")
  # A matching n_vars is accepted and recorded.
  expect_identical(
    efa_power("simulation", Lambda = sim_Lambda, n_vars = nrow(sim_Lambda), N = 120,
              n_datasets = 2, criteria = "EKC", seed = 1)$settings$n_vars,
    nrow(sim_Lambda))
  # N is required in simulation mode (it is optional in RMSEA mode, where it is solved
  # for), both N and n_datasets must be positive whole numbers, and N must exceed the
  # number of variables before anything is drawn.
  expect_error(efa_power("simulation", Lambda = sim_Lambda), class = "efa_power_input")
  expect_error(efa_power("simulation", Lambda = sim_Lambda, N = -5, n_datasets = 2),
               class = "efa_power_input")
  expect_error(efa_power("simulation", Lambda = sim_Lambda, N = 200, n_datasets = 0),
               class = "efa_power_input")
  expect_error(efa_power("simulation", Lambda = sim_Lambda, N = nrow(sim_Lambda),
                         n_datasets = 500), class = "efa_power_input")
  # Out-of-range recovery threshold.
  expect_error(efa_power("simulation", Lambda = sim_Lambda, N = 200, n_datasets = 2,
                         recovery_threshold = 0), class = "efa_power_bad_threshold")
  # A visual criterion makes no numeric suggestion to score.
  expect_error(efa_power("simulation", Lambda = sim_Lambda, N = 200, n_datasets = 2,
                         criteria = "SCREE"))
  # A simulation object has no analytic power curve to plot.
  sim <- efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 120,
                   n_datasets = 3, criteria = "EKC", seed = 42)
  expect_error(plot(sim), class = "efa_power_no_plot")
})

test_that("simulation mode runs the internally-simulating criteria", {
  skip_if_not_slow()
  sim <- efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 250,
                   n_datasets = 10, criteria = c("EKC", "NEST", "PARALLEL"),
                   seed = 42)
  expect_true("EKC_BvA2017" %in% names(sim$hit_rate))
  expect_true(any(grepl("NEST", names(sim$hit_rate))))
  # the efa_retain() default: PCA eigenvalues
  expect_true("PARALLEL_PCA" %in% names(sim$hit_rate))
  expect_gte(sim$hit_rate[["EKC_BvA2017"]], 0.8)
})

test_that("simulation-mode printed output is stable", {
  local_reproducible_output()
  sim <- efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 200,
                   n_datasets = 10, criteria = c("EKC", "MAP"), seed = 42)
  # The rates move with the BLAS, so every decimal is masked; the snapshot holds the
  # layout, the section order, and the criterion labels. The rates are pinned by the
  # expect_gte() bounds in "simulation mode recovers a clean population".
  expect_snapshot(print(sim), transform = scrub_num)
})
