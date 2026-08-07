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

test_that("df is taken from p and k via the EFA df formula", {
  # p = 20, k = 5 -> ((20 - 5)^2 - (20 + 5)) / 2 = 100
  from_pk <- efa_power(p = 20, k = 5, N = 200)
  from_df <- efa_power(df = 100, N = 200)
  expect_equal(from_pk$settings$df, 100)
  expect_equal(from_pk$power, from_df$power)
})

test_that("the number of groups scales the noncentrality", {
  g1 <- efa_power(df = 100, N = 200, group = 1)
  g2 <- efa_power(df = 100, N = 200, group = 2)
  expect_equal(unname(g2$ncp), unname(g1$ncp) / 2)
})

test_that("the solved N is the total across groups, split by N_per_group", {
  # The 1 / group factor in the noncentrality makes the solved N a total, so it
  # grows linearly in `group` while the per-group requirement stays near-constant.
  g1 <- efa_power(df = 102, power = 0.80, group = 1)
  g2 <- efa_power(df = 102, power = 0.80, group = 2)
  expect_equal(g1$N, 130L)
  expect_equal(g2$N, 259L)
  expect_equal(g2$N_per_group, 129.5)
  expect_equal(g1$N_per_group, 130)
})

test_that("the solved total matches semTools::findRMSEAsamplesize per group", {
  skip_if_not_installed("semTools")

  # findRMSEAsamplesize() searches the same noncentrality but reports the per-group
  # figure (its last line divides the solved total by `group`), so its value is the
  # N_per_group counterpart of the total efa_power() returns.
  # Compared as named vectors so a failure names the group count that broke.
  groups <- c(1, 2, 3)
  ref <- vapply(groups, function(g) {
    semTools::findRMSEAsamplesize(rmsea0 = 0.05, rmseaA = 0.08, df = 102,
                                  power = 0.80, group = g)
  }, numeric(1))
  pw <- vapply(groups, function(g) {
    efa_power(df = 102, power = 0.80, group = g)$N_per_group
  }, numeric(1))
  tot <- vapply(groups, function(g) {
    as.numeric(efa_power(df = 102, power = 0.80, group = g)$N)
  }, numeric(1))
  names(ref) <- names(pw) <- names(tot) <- paste0("group", groups)

  expect_equal(pw, ref, tolerance = 1e-8)
  expect_equal(tot, ref * groups, tolerance = 1e-8, ignore_attr = "names")
})

test_that("output has the expected class and structure", {
  x <- efa_power(df = 100, N = 200)
  expect_s3_class(x, "efa_power")
  expect_named(x, c("power", "N", "N_per_group", "crit", "ncp", "solve_for",
                    "settings"))
  expect_named(x$settings, c("mode", "type", "eps0", "eps1", "df", "p", "k",
                             "alpha", "group", "power"))
})

test_that("classed conditions fire as expected", {
  # Equal null and alternative RMSEA: nothing to detect.
  expect_error(efa_power(df = 100, N = 200, eps0 = 0.05, eps1 = 0.05),
               class = "efa_power_unreachable")
  # Non-positive degrees of freedom, given directly or from an unidentified model.
  expect_error(efa_power(df = 0, N = 200), class = "efa_power_bad_df")
  expect_error(efa_power(p = 3, k = 3, N = 200), class = "efa_power_bad_df")
  # No way to obtain df.
  expect_error(efa_power(N = 200), class = "efa_power_missing_df")
  # Both N and power supplied: nothing to solve for.
  expect_error(efa_power(df = 100, N = 200, power = 0.80),
               class = "efa_power_overdetermined")
  # Out-of-range alpha and target power.
  expect_error(efa_power(df = 100, N = 200, alpha = 0),
               class = "efa_power_bad_alpha")
  expect_error(efa_power(df = 100, power = 1.2), class = "efa_power_bad_power")

  # eps0/eps1 on the wrong side for the chosen test: a message, and the value is
  # still returned (the branch follows `type`, not the ordering of the epsilons).
  expect_message(
    val <- efa_power(df = 100, N = 200, type = "close", eps0 = 0.08, eps1 = 0.05),
    class = "efa_power_wrong_side")
  expect_s3_class(val, "efa_power")
})

test_that("an unreachable target aborts cleanly without warning noise", {
  # Wrong-side epsilons make close-fit power fall short at every N, so the search
  # runs up to the cap and aborts (the wrong-side message is muffled here).
  expect_error(
    suppressMessages(
      efa_power(df = 100, power = 0.80, type = "close", eps0 = 0.08, eps1 = 0.05)),
    class = "efa_power_unreached")

  # The search evaluates the noncentral chi-square routines at extreme
  # noncentrality; those base-R precision/convergence notes must not leak.
  expect_no_warning(suppressMessages(tryCatch(
    efa_power(df = 100, power = 0.80, type = "close", eps0 = 0.08, eps1 = 0.05),
    error = function(e) NULL)))
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
  # Not-close fit with more than one group still plots.
  expect_s3_class(plot(efa_power(df = 100, N = 200, type = "notclose", group = 2)),
                  "ggplot")

  # Wrong-side epsilons: the automatic N axis cannot solve for 0.99 power, so the
  # curve span falls back to a multiple of the object's N rather than erroring out.
  wrong_side <- suppressMessages(
    efa_power(df = 100, N = 200, type = "close", eps0 = 0.08, eps1 = 0.05))
  expect_s3_class(suppressMessages(plot(wrong_side)), "ggplot")
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
  # stream to its index, so the result must be identical regardless of the number of
  # workers. The multisession workers are fresh R processes that load the installed
  # package, so run this under devtools::check() / after devtools::install().
  # (multicore is unavailable on Windows.)
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
  expect_equal(sim$k_true, 3L)
  expect_named(sim, c("hit_rate", "hits", "recovery", "convergence", "replicates",
                      "k_true", "model_error", "settings"))

  # EKC and the loadings recovery should both be near-perfect on a clean 3-factor model.
  expect_gte(sim$hit_rate[["EKC_BvA2017"]], 0.9)
  expect_gte(sim$recovery$min_rate, 0.9)
  expect_equal(sim$recovery$threshold, 0.95)
  expect_equal(sim$convergence$convergence_rate, 1)
  expect_equal(sim$convergence$heywood_rate, 0)

  # hits data frame is consistent with the hit_rate vector.
  expect_equal(sim$hits$hit_rate, unname(sim$hit_rate))
  expect_true(all(sim$hits$hits <= sim$hits$n_valid))
})

test_that("simulation mode muffles the per-replicate fit messages", {
  expect_no_message(
    efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 120,
              n_datasets = 3, criteria = "EKC", seed = 1))
})

test_that("a criterion that fails on every replicate is reported, not dropped", {
  # N < p makes each replicate's correlation matrix singular, so MAP errors on every
  # replicate. The criterion must still appear in the results -- with no valid replicates
  # and an NA hit-rate -- and be named in a classed warning, rather than vanishing while
  # the settings still record the request.
  expect_warning(
    sim <- efa_power("simulation", Lambda = sim_Lambda, N = 6, n_datasets = 3,
                     criteria = "MAP", seed = 1),
    class = "efa_power_criterion_failed")
  expect_s3_class(sim, "efa_power")
  expect_named(sim$hit_rate, "MAP")
  expect_true(is.na(sim$hit_rate[["MAP"]]))
  expect_equal(sim$hits$n_valid, 0L, ignore_attr = TRUE)
  expect_equal(sim$hits$hits, 0, ignore_attr = TRUE)
  expect_identical(sim$settings$criteria, "MAP")
  # The formatter reports the missing rate rather than dropping the line or erroring.
  local_reproducible_output()
  expect_snapshot(print(sim))
})

test_that("a criterion that decides on every replicate raises no failure warning", {
  expect_no_warning(
    efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 150,
              n_datasets = 4, criteria = "EKC", seed = 3))
})

test_that("simulation mode: R-only population turns recovery off but keeps the hit-rate", {
  R_pop <- efa_simulate(Lambda = sim_Lambda, Phi = sim_Phi, return_pop = TRUE)$population
  rr <- efa_power("simulation", R = R_pop, k = 3, N = 200, n_datasets = 6,
                  criteria = "EKC", seed = 42)
  expect_null(rr$recovery)
  expect_false(rr$settings$has_recovery)
  expect_equal(rr$k_true, 3L)
  expect_true(is.finite(rr$hit_rate[["EKC_BvA2017"]]))
})

test_that("simulation-mode classed conditions fire as expected", {
  R_pop <- efa_simulate(Lambda = sim_Lambda, Phi = sim_Phi, return_pop = TRUE)$population

  # The population must be given exactly one way.
  expect_error(efa_power("simulation", N = 200, n_datasets = 2),
               class = "efa_power_input")
  expect_error(efa_power("simulation", Lambda = sim_Lambda, R = R_pop, N = 200,
                         n_datasets = 2), class = "efa_power_input")
  # k is required with a bare R, and must agree with Lambda otherwise.
  expect_error(efa_power("simulation", R = R_pop, N = 200, n_datasets = 2),
               class = "efa_power_missing_k")
  expect_error(efa_power("simulation", Lambda = sim_Lambda, k = 2, N = 200,
                         n_datasets = 2), class = "efa_power_bad_k")
  # N is required in simulation mode (it is optional in RMSEA mode, where it is solved
  # for), and both N and n_datasets must be positive whole numbers.
  expect_error(efa_power("simulation", Lambda = sim_Lambda), class = "efa_power_input")
  expect_error(efa_power("simulation", Lambda = sim_Lambda, N = -5, n_datasets = 2),
               class = "efa_power_input")
  expect_error(efa_power("simulation", Lambda = sim_Lambda, N = 200, n_datasets = 0),
               class = "efa_power_input")
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
  expect_true(any(grepl("PARALLEL", names(sim$hit_rate))))
  expect_gte(sim$hit_rate[["EKC_BvA2017"]], 0.8)
})

test_that("simulation-mode printed output is stable", {
  local_reproducible_output()
  sim <- efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 200,
                   n_datasets = 10, criteria = c("EKC", "MAP"), seed = 42)
  # Every decimal here -- the retention hit rates and the median minimum congruence alike --
  # is averaged over ten simulated fits and moves with the BLAS, so all of them are masked.
  # What this snapshot holds fixed is the layout, the section order, the criterion labels and
  # the integer replicate counts. The rates themselves are pinned numerically instead, by the
  # expect_gte() bounds at lines 268 and 361.
  expect_snapshot(print(sim), transform = scrub_num)
})
