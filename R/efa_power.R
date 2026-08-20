#' Power analysis for exploratory factor analysis
#'
#' @description
#' Analyses power for exploratory factor analysis, in one of two modes chosen with
#' `mode`.
#'
#' `mode = "rmsea"` (the default) gives the analytic power of the root mean square
#' error of approximation (RMSEA) tests of close and not-close fit (MacCallum,
#' Browne, & Sugawara, 1996). Give a sample size to get the power of the test, or
#' give a target power to get the sample size needed to reach it.
#'
#' `mode = "simulation"` runs a Monte-Carlo study: it draws `n_datasets` samples from
#' a known population (via [efa_simulate()]), analyses each one, and reports how well
#' the analysis recovers that population. See *Details* for what is reported.
#'
#' Here the number of variables is `p` and the number of factors is `k` (elsewhere in
#' the package: `n_vars` and `n_factors`).
#'
#' @details
#' # RMSEA mode
#'
#' Power rises with a larger sample, a larger model (more degrees of freedom), and a
#' bigger gap between the null and alternative RMSEA (MacCallum, Browne, & Sugawara,
#' 1996).
#'
#' Two tests are supported, chosen with `type` (never by the order of `eps0` and
#' `eps1`):
#' \describe{
#'   \item{`"close"`}{Tests *close fit* (MacCallum et al., 1996). The null
#'     hypothesis is that the fit is close (RMSEA \eqn{\le} `eps0`; conventionally
#'     0.05). Power is the chance of detecting a worse alternative (`eps1`;
#'     conventionally 0.08, so `eps0 < eps1`), in the upper tail.}
#'   \item{`"notclose"`}{Tests *not-close fit*. The null hypothesis is that the fit
#'     is not close (RMSEA \eqn{\ge} `eps0`). Power is the chance of detecting a
#'     better alternative (`eps1`; conventionally 0.01, so `eps0 > eps1`), in the
#'     lower tail.}
#' }
#' When `eps0` and `eps1` are in the wrong order for the chosen `type`, a message is
#' shown but the requested test still runs. Equal `eps0` and `eps1` leave nothing to
#' detect and are an error.
#'
#' Power always increases with `N`, so the required sample size (the smallest `N`
#' reaching `power`) is found by bisection. `N` is the **total** sample size across
#' groups: with `group > 1` the power calculation divides by `group` (the
#' `1 / group` factor), so spreading a fixed total over more groups gives less
#' power. The matching per-group sample size, `N / group`, is returned as
#' `N_per_group`.
#'
#' The `1 / group` factor makes all `group` groups the same size, so a required
#' total is rounded up to the next multiple of `group`. A solved `N_per_group` is thus
#' a whole number of persons, and the reported `power` is the power at a total that a
#' study can collect. With `group = 2` and `df = 102`, for example, the required total
#' is 260, or 130 per group. Bisection on the total alone gives 259, which asks for
#' 129.5 persons in each group.
#'
#' # Simulation mode
#'
#' The population is passed to [efa_simulate()], which draws `n_datasets` samples of
#' size `N` from it. The population's true number of factors, `k_true`, is
#' `ncol(Lambda)` for a factor-model population, or `k` for a bare `R`. By default
#' the population fits the factor model exactly, which overstates how well the
#' criteria and the fit recover its structure; setting a misfit target makes the
#' population more realistic (MacCallum, 2003).
#'
#' Each replicate is analysed three ways:
#' \describe{
#'   \item{**Hit-rate**}{The share of replicates where a criterion's suggested factor
#'     count (from `criteria`) matches `k_true`. A replicate where the criterion
#'     errored or gave no answer is left out of this count -- it does not count as a
#'     miss.}
#'   \item{**Structure recovery** (factor-model populations only)}{The `k_true`-factor
#'     model is fitted with [efa_fit()], its loadings are matched to the population
#'     loadings, and the matched-factor Tucker congruences (Lorenzo-Seva & ten Berge,
#'     2006) are compared with `recovery_threshold`. A replicate succeeds when its
#'     smallest (`min`) or average (`mean`) matched congruence reaches the
#'     threshold.}
#'   \item{**Convergence**}{Among the replicates whose fit completed, the share that
#'     converged and the share that produced a Heywood case.}
#' }
#' A replicate whose fit fails completely is not counted in any of the three
#' measures above. If any fit fails, a warning reports how many failed and the cause
#' of the first failure.
#'
#' Replicates are analysed in parallel with \pkg{future.apply}; choose a parallel
#' plan with [future::plan()]. Each replicate uses its own reproducible
#' random-number stream, so with a fixed `seed` the result does not depend on the
#' number of workers, and the caller's random-number state is left unchanged.
#'
#' @param mode character. The kind of power analysis: `"rmsea"` (the default; analytic
#'   RMSEA power) or `"simulation"` (Monte-Carlo hit-rate and structure recovery).
#'   `type`, `eps0`, `eps1`, `df`, `alpha`, `power`, and `group` apply to RMSEA mode only;
#'   the arguments marked *Simulation mode* below apply to the other; `N`, `p`, and `k` are
#'   used in both.
#' @param type character. The RMSEA test: `"close"` (test of close fit) or
#'   `"notclose"` (test of not-close fit). See *Details*.
#' @param eps0 numeric. The null-hypothesis RMSEA. Default is `0.05`.
#' @param eps1 numeric. The alternative-hypothesis RMSEA (the true RMSEA power is
#'   evaluated at). Default is `0.08` for `type = "close"` and `0.01` for
#'   `type = "notclose"`.
#' @param N numeric. In `"rmsea"` mode, the total sample size across groups (the
#'   plain sample size when `group` is `1`): give `N` to compute power, or leave it
#'   `NULL` to solve for the required `N` at a target `power`. In `"simulation"`
#'   mode `N` is required and is the size of each drawn sample, with no sample size
#'   solved for and no `group` division.
#' @param p numeric. The number of observed variables. In `"rmsea"` mode, used with
#'   `k` to derive `df` when `df` is not given directly. In `"simulation"` mode it is
#'   read off the population, so leave it unset (or matching `nrow(Lambda)` / `nrow(R)`).
#' @param k numeric. The number of factors. In `"rmsea"` mode, used with `p` to
#'   derive `df` when `df` is not given directly. In `"simulation"` mode it is the
#'   true number of factors: it is required with an `R` population and must be left
#'   unset (or match `ncol(Lambda)`) with a factor-model population.
#' @param df numeric. The model degrees of freedom. Either supply `df` directly or
#'   supply both `p` and `k`, from which `df = ((p - k)^2 - (p + k)) / 2`. Must be
#'   positive.
#' @param alpha numeric. The significance level. Default is `0.05`.
#' @param power numeric. The target power. Give `power` (or leave both `power` and
#'   `N` `NULL`, defaulting to `0.80`) to solve for the required `N`; leave it
#'   `NULL` while giving `N` to compute power. Exactly one of `N` and `power` is
#'   solved for.
#' @param group numeric. The number of groups. Default is `1`. `N` is the total
#'   across all `group` groups, not the size of each one, and a solved `N` is a
#'   multiple of `group`. See *Details*.
#' @param Lambda matrix. Simulation mode. A `p` by `k_true` population loading matrix.
#'   Supply this (optionally with `Phi`/`Psi`) to build a factor-model population;
#'   structure recovery is available only with this form. Passed to [efa_simulate()].
#' @param Phi matrix. Simulation mode. The `k_true` by `k_true` population factor
#'   intercorrelations. Only used with `Lambda`; defaults to orthogonal factors. When
#'   `rotation` is unset, an oblique `Phi` selects a `"promax"` recovery fit and an
#'   orthogonal one a `"varimax"` fit.
#' @param Psi numeric or matrix. Simulation mode. The population unique variances (a
#'   length-`p` vector or a `p` by `p` matrix). Only used with `Lambda`. Passed to
#'   [efa_simulate()].
#' @param R matrix. Simulation mode. A `p` by `p` population correlation matrix to draw
#'   from directly, instead of a factor model. Structure recovery is not available for
#'   this form (there are no population loadings to recover), and `k` is required.
#' @param n_datasets numeric. Simulation mode. The number of samples to draw and
#'   analyse. Default is `500`.
#' @param criteria character. Simulation mode. The factor-retention criteria to
#'   evaluate the hit-rate for, any of `"CD"`, `"EKC"`, `"HULL"`, `"KGC"`, `"MAP"`,
#'   `"NEST"`, `"PARALLEL"`, and `"SMT"` (see [efa_retain()]). Default is
#'   `c("EKC", "MAP")`. Criteria that simulate internally (`"CD"`, `"HULL"`,
#'   `"NEST"`, `"PARALLEL"`) make each run substantially slower.
#' @param estimator character. Simulation mode. The estimator (`"PAF"`, `"ML"`,
#'   or `"ULS"`) used for the recovery fit and the retention criteria. Default
#'   is `"PAF"`.
#' @param rotation character. Simulation mode. The rotation for the recovery fit,
#'   passed to [efa_fit()]. Default is `NULL`, which matches the population: `"varimax"`
#'   for orthogonal factors and `"promax"` for oblique ones (a single factor is left
#'   unrotated). Recovery aligns the fitted loadings to the population pattern by
#'   permutation and sign only. A rotation that does not seek that structure -- for
#'   example `"none"` with more than one factor -- will understate recovery, so keep
#'   the default (or another structure-seeking rotation) for a meaningful recovery
#'   rate.
#' @param recovery_threshold numeric. Simulation mode. The matched-factor Tucker
#'   congruence a replicate must reach to count as recovered. Default is `0.95`:
#'   Lorenzo-Seva and ten Berge (2006) treat congruence at or above this level as
#'   indicating the same factor.
#' @param model_error character. Simulation mode. The [efa_simulate()] method that
#'   perturbs the population with model error: `"TKL"` (Tucker-Koopman-Linn, the
#'   default here), `"CB"` (Cudeck-Browne), `"WB"` (Wu-Browne), or `"none"` for an
#'   exact population. It only takes effect when a target is supplied (`target_rmsea`
#'   and/or `target_cfi`), and only for a factor-model population; without a target
#'   the population stays exact whatever the method. `"TKL"` adds minor common
#'   factors, giving a realistically imperfect population but lowering both the
#'   hit-rate and structure recovery; `"CB"` and `"WB"` target the RMSEA only. `"CB"`
#'   keeps the population loadings as the exact minimizer of the perturbed
#'   population, so recovery stays close to perfect; `"WB"`'s loadings are not the
#'   minimizer, so they carry no such guarantee. Note that [efa_simulate()] itself
#'   defaults to `"CB"`: the same `target_rmsea` passed to both functions gives an
#'   easier population there, unless `model_error` is also set explicitly here.
#' @param target_rmsea numeric. Simulation mode. The population RMSEA the model should
#'   have relative to the perturbed population, activating model error. Default is
#'   `NULL`. Passed to [efa_simulate()].
#' @param target_cfi numeric. Simulation mode. The population CFI target (only with
#'   `model_error = "TKL"`). Default is `NULL`. Passed to [efa_simulate()].
#' @param seed numeric. Simulation mode. Optional seed making the draws and analysis
#'   reproducible and worker-count independent; the caller's random-number stream is
#'   restored afterwards. Default is `NULL`.
#'
#' @returns An object of class `efa_power`. For `mode = "rmsea"`, a list containing:
#' \item{power}{The power of the test at `N` (the achieved power, which for a
#'   solved sample size is at least the target).}
#' \item{N}{The total sample size across groups: the supplied `N`, or the solved
#'   required sample size (a multiple of `group`).}
#' \item{N_per_group}{The per-group sample size `N / group`, equal to `N` when
#'   `group` is `1`. A whole number for a solved `N`; for a supplied `N` that is not
#'   a multiple of `group` it is the fraction that the noncentrality uses.}
#' \item{crit}{The critical chi-square value the fit statistic is compared against.}
#' \item{ncp}{The noncentrality parameters under the null (`H0`, from `eps0`) and
#'   the alternative (`H1`, from `eps1`).}
#' \item{solve_for}{`"power"` or `"N"`, recording which quantity was solved for.}
#' \item{settings}{A list of the inputs: `mode`, `type`, `eps0`, `eps1`, `df`,
#'   `p`, `k`, `alpha`, `group`, and the target `power` (the value solved to when
#'   `solve_for` is `"N"`, otherwise `NULL`).}
#'
#' For `mode = "simulation"`, a list containing:
#' \item{hit_rate}{A named numeric vector of the retention hit-rate per criterion (and,
#'   where a criterion has several variants, per variant); `NA` for a criterion that
#'   returned no suggestion on any replicate.}
#' \item{hits}{A data frame with one row per criterion (`criterion`) giving the
#'   number of replicates it returned a definite suggestion on (`n_valid`), the
#'   number of those that matched `k_true` (`hits`), and the `hit_rate`
#'   (`hits / n_valid`).}
#' \item{recovery}{For a factor-model population, a list with the structure-recovery
#'   rates (`min_rate`, `mean_rate`), the `threshold`, and the number of usable fits
#'   (`n_valid`); `NULL` for an `R` population. Rates are over every replicate whose fit
#'   returned loadings, including non-converged or Heywood solutions (their rates are
#'   reported separately in `convergence`).}
#' \item{convergence}{A list with the number of datasets (`n_datasets`), the number
#'   of fits that completed (`n_fit_ok`), how many of those converged
#'   (`n_converged`) and how many produced a Heywood case (`n_heywood`), and the
#'   corresponding rates: `fit_rate` (fits completed, over all datasets) and
#'   `convergence_rate` / `heywood_rate` (converged / Heywood, over the completed
#'   fits).}
#' \item{replicates}{The raw per-replicate values: the suggested factor counts
#'   (`n_hat`), the matched congruences (`rec_min`, `rec_mean`), the `converged`,
#'   `heywood`, and `fit_ok` flags, and `fit_error`, the message of the fit that did
#'   not complete (`NA` where it did).}
#' \item{k_true}{The true number of factors.}
#' \item{model_error}{The [efa_simulate()] model-error record, or `NULL`.}
#' \item{settings}{A list of the simulation inputs.}
#'
#' @references
#' MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power analysis and
#' determination of sample size for covariance structure modeling. *Psychological
#' Methods, 1*(2), 130-149. \doi{10.1037/1082-989X.1.2.130}
#'
#' MacCallum, R. C. (2003). 2001 Presidential Address: Working with imperfect models.
#' *Multivariate Behavioral Research, 38*(1), 113-139. \doi{10.1207/S15327906MBR3801_5}
#'
#' Lorenzo-Seva, U., & ten Berge, J. M. F. (2006). Tucker's congruence coefficient as a
#' meaningful index of factor similarity. *Methodology, 2*(2), 57-64.
#' \doi{10.1027/1614-2241.2.2.57}
#'
#' @seealso [efa_simulate()] draws the replicate datasets used in simulation mode.
#'   [efa_retain()] implements the retention criteria whose hit-rates simulation
#'   mode reports.
#'
#' @family power analysis
#'
#' @export
#'
#' @examples
#' # Power of the test of close fit at N = 200 for a 100-df model
#' efa_power(df = 100, N = 200)
#'
#' # Deriving df from the model dimensions instead of giving it directly
#' efa_power(p = 20, k = 3, N = 200)
#'
#' # Required total sample size for 80% power
#' efa_power(df = 100, power = 0.80)
#'
#' # Test of not-close fit
#' efa_power(df = 100, N = 200, type = "notclose")
#'
#' \donttest{
#' # Simulation mode: retention hit-rate and structure recovery for a known
#' # three-factor population at N = 300
#' efa_power("simulation", Lambda = population_models$loadings$baseline,
#'           Phi = population_models$phis_3$moderate, N = 300,
#'           n_datasets = 50, criteria = c("EKC", "MAP"), seed = 42)
#' }
#'
efa_power <- function(mode = c("rmsea", "simulation"),
                      type = c("close", "notclose"),
                      eps0 = NULL, eps1 = NULL, N = NULL, p = NULL, k = NULL,
                      df = NULL, alpha = 0.05, power = NULL, group = 1,
                      Lambda = NULL, Phi = NULL, Psi = NULL, R = NULL,
                      n_datasets = 500, criteria = c("EKC", "MAP"),
                      estimator = "PAF", rotation = NULL, recovery_threshold = 0.95,
                      model_error = c("TKL", "CB", "WB", "none"),
                      target_rmsea = NULL, target_cfi = NULL, seed = NULL) {

  mode <- .match_arg_ci(mode)

  # Simulation mode is a self-contained Monte-Carlo path with its own inputs (the
  # population, the retention criteria, the recovery fit); the analytic RMSEA
  # arguments do not apply to it, so it branches off before they are validated.
  if (mode == "simulation") {
    return(.efa_power_simulation(
      N = N, Lambda = Lambda, Phi = Phi, Psi = Psi, R = R, p = p, k = k,
      n_datasets = n_datasets, criteria = criteria, estimator = estimator,
      rotation = rotation, recovery_threshold = recovery_threshold,
      model_error = model_error, target_rmsea = target_rmsea,
      target_cfi = target_cfi, seed = seed))
  }

  type <- .match_arg_ci(type)

  # RMSEA defaults: a common null of .05, tested against a worse (.08) alternative
  # for close fit or a better (.01) alternative for not-close fit.
  if (is.null(eps0)) eps0 <- 0.05
  if (is.null(eps1)) eps1 <- if (type == "close") 0.08 else 0.01
  checkmate::assert_number(eps0, lower = 0, finite = TRUE)
  checkmate::assert_number(eps1, lower = 0, finite = TRUE)
  checkmate::assert_count(group, positive = TRUE)
  if (!checkmate::test_number(alpha) || alpha <= 0 || alpha >= 1) {
    cli::cli_abort("{.arg alpha} must be a single number strictly between 0 and 1.",
                   class = "efa_power_bad_alpha")
  }

  # Degrees of freedom: taken from `df` directly, or derived from the model
  # dimensions via the shared EFA df formula. A non-positive df has no test.
  # `p` / `k` are validated whenever supplied -- even alongside an explicit `df` --
  # so a malformed value never lands in the returned settings.
  if (!is.null(p)) checkmate::assert_count(p, positive = TRUE)
  if (!is.null(k)) checkmate::assert_count(k, positive = TRUE)
  if (is.null(df)) {
    if (is.null(p) || is.null(k)) {
      cli::cli_abort(
        c("The model degrees of freedom are missing.",
          "i" = "Supply {.arg df} directly, or both {.arg p} (variables) and {.arg k} (factors)."),
        class = "efa_power_missing_df")
    }
    df <- .efa_df(p, k)
  } else {
    checkmate::assert_number(df, finite = TRUE)
  }
  if (df <= 0) {
    cli::cli_abort(
      c("The model has {df} degree{?s} of freedom; a positive value is required.",
        "i" = "Use more variables or fewer factors."),
      class = "efa_power_bad_df")
  }

  # Equal null and alternative RMSEA cannot be told apart, so no sample size (and
  # no power beyond alpha) exists.
  if (eps0 == eps1) {
    cli::cli_abort(
      c("{.arg eps0} and {.arg eps1} are both {eps0}, so the test has no power to detect a difference.",
        "i" = "Give them different RMSEA values (for example {.val 0.05} vs {.val 0.08} for a close-fit test)."),
      class = "efa_power_unreachable")
  }

  # The branch is chosen by `type`, not by the ordering of eps0/eps1; a reversed
  # pair is reported but still computed as requested.
  wrong_side <- (type == "close" && eps0 > eps1) ||
    (type == "notclose" && eps0 < eps1)
  if (wrong_side) {
    expected <- if (type == "close") {
      "A close-fit test expects the null RMSEA below the alternative ({.arg eps0} < {.arg eps1})."
    } else {
      "A not-close-fit test expects the null RMSEA above the alternative ({.arg eps0} > {.arg eps1})."
    }
    cli::cli_inform(
      c("!" = "For a {.val {type}}-fit test, {.arg eps0} ({eps0}) and {.arg eps1} ({eps1}) look reversed.",
        "i" = expected),
      class = "efa_power_wrong_side")
  }

  # Solve for whichever of N / power is left NULL; supplying both is over-determined.
  if (!is.null(N) && !is.null(power)) {
    cli::cli_abort(
      c("Both {.arg N} and {.arg power} were supplied, so there is nothing to solve for.",
        "i" = "Leave {.arg N} as {.code NULL} to find the required sample size, or {.arg power} as {.code NULL} to find the power at {.arg N}."),
      class = "efa_power_overdetermined")
  }

  if (is.null(N)) {
    solve_for <- "N"
    target <- if (is.null(power)) 0.80 else power
    if (!checkmate::test_number(target) || target <= 0 || target >= 1) {
      cli::cli_abort("{.arg power} must be a single number strictly between 0 and 1.",
                     class = "efa_power_bad_power")
    }
    N <- .efa_power_solve_N(df, eps0, eps1, alpha, group, type, target)
    settings_power <- target
  } else {
    solve_for <- "power"
    checkmate::assert_count(N, positive = TRUE)
    N <- as.integer(N)
    settings_power <- NULL
  }

  res <- .efa_power_rmsea(N, df, eps0, eps1, alpha, group, type)

  structure(
    list(
      power = res$power,
      N = N,
      N_per_group = N / group,
      crit = res$crit,
      ncp = c(H0 = res$ncp0, H1 = res$ncp1),
      solve_for = solve_for,
      settings = list(mode = mode, type = type, eps0 = eps0, eps1 = eps1, df = df,
                      p = p, k = k, alpha = alpha, group = group,
                      power = settings_power)
    ),
    class = "efa_power"
  )
}

# Power of the RMSEA test of (not-)close fit at a single sample size, from the
# noncentral chi-square framework of MacCallum, Browne & Sugawara (1996). The
# noncentrality parameter is (N - 1) * df * eps^2 / group. For close fit the
# rejection region is the upper tail (critical value = upper-alpha point of the
# null); for not-close fit it is the lower tail. Returns the power, the critical
# value, and both noncentrality parameters.
.efa_power_rmsea <- function(N, df, eps0, eps1, alpha, group, type) {
  ncp0 <- (N - 1) * df * eps0^2 / group
  ncp1 <- (N - 1) * df * eps1^2 / group
  # At very large noncentrality (huge N) qchisq()/pchisq() emit precision and
  # non-convergence notes while still returning valid quantiles/probabilities. Those
  # notes carry no information for the caller, so suppress them and keep power
  # reporting to the function's own (classed) conditions.
  suppressWarnings(
    if (type == "close") {
      crit <- stats::qchisq(alpha, df, ncp = ncp0, lower.tail = FALSE)
      pow <- stats::pchisq(crit, df, ncp = ncp1, lower.tail = FALSE)
    } else {
      crit <- stats::qchisq(1 - alpha, df, ncp = ncp0, lower.tail = FALSE)
      pow <- 1 - stats::pchisq(crit, df, ncp = ncp1, lower.tail = FALSE)
    }
  )
  list(power = pow, crit = crit, ncp0 = ncp0, ncp1 = ncp1)
}

# Smallest total sample size reaching `target` power that a design with `group` equal
# groups can have. Power rises monotonically with N (from about `alpha` toward 1), so
# double an upper bound until it clears the target, then bisect. `cap` bounds the search:
# an alternative RMSEA too close to the null (or a target power too high) can otherwise
# never be reached.
.efa_power_solve_N <- function(df, eps0, eps1, alpha, group, type, target,
                               cap = 1e7) {
  powfun <- function(N) .efa_power_rmsea(N, df, eps0, eps1, alpha, group, type)$power

  # The noncentrality divides the total by `group`, so the design has `group` groups of
  # equal size and only a total that is a multiple of `group` is possible. Power increases
  # with the total, so the smallest possible total is the bisected one rounded up to the
  # next multiple. Without this step a total of, say, 259 over two groups asks for 129.5
  # persons in each group, which is not a sample size.
  design_N <- function(N) as.integer(ceiling(N / group) * group)

  # Rounding up happens after the search, so the search itself must stop at the largest
  # collectable total within `cap`; otherwise a returned total could pass the bound that
  # the abort message below states (10000002 for `group` = 3 at the default cap).
  cap <- floor(cap / group) * group

  lo <- 1
  if (powfun(lo) >= target) return(design_N(1))
  # Double the upper bound, clamped to `cap`, until it clears the target. Probing
  # `cap` itself before giving up keeps the reachable ceiling equal to the documented
  # bound (a plain `hi * 2` would abort once past the last power of two below `cap`).
  hi <- 2
  while (powfun(hi) < target) {
    if (hi >= cap) {
      cli::cli_abort(
        c("No sample size up to {.val {cap}} reaches a power of {target}.",
          "i" = "The alternative RMSEA may be too close to the null, or the target power too high."),
        class = "efa_power_unreached")
    }
    hi <- min(hi * 2, cap)
  }
  while (hi - lo > 1) {
    mid <- floor((lo + hi) / 2)
    if (powfun(mid) >= target) hi <- mid else lo <- mid
  }
  design_N(hi)
}


# Simulation-mode power: draw `n_datasets` samples of size `N` from a known
# population, and over the replicates report how often each factor-retention
# criterion recovers the true factor count (hit-rate), how often the fitted
# loadings recover the population structure (Tucker congruence >= threshold), and
# the convergence/Heywood rate of the k_true-factor fit. Parallelised over
# replicates with a per-replicate reproducible RNG stream.
.efa_power_simulation <- function(N, Lambda, Phi, Psi, R, p, k, n_datasets,
                                  criteria, estimator, rotation, recovery_threshold,
                                  model_error, target_rmsea, target_cfi, seed) {

  model_error <- .match_arg_ci(model_error, c("TKL", "CB", "WB", "none"))
  estimator <- .match_arg_ci(estimator, c("PAF", "ML", "ULS"))
  # Only criteria that make a numeric suggestion can score a hit-rate; the visual
  # scree plot is excluded.
  valid_ids <- names(.retention_registry)[
    !vapply(.retention_registry, function(e) isTRUE(e$visual), logical(1))]
  # `criteria` is always passed on from efa_power(), so this function's own formal carries no
  # default for .match_arg_ci() to fall back on when the caller supplied NULL. Read it off
  # efa_power() instead of restating it, so the documented default stays in one place.
  if (is.null(criteria)) criteria <- eval(formals(efa_power)$criteria)
  criteria <- .match_arg_ci(criteria, valid_ids, several.ok = TRUE)

  # `N` is optional in RMSEA mode (it is solved for at a target power), which makes
  # omitting it the most likely simulation-mode mistake; there is nothing to solve for
  # here, so report it as a requirement rather than as a bare type assertion.
  if (is.null(N)) {
    cli::cli_abort(
      c("{.arg N} is required in simulation mode.",
        "x" = "It is the size of each drawn sample; no sample size is solved for here.",
        "i" = "Set {.arg N}, e.g. {.code N = 300}, or use {.code mode = \"rmsea\"} to solve for the required sample size."),
      class = "efa_power_input")
  }
  if (!checkmate::test_count(N, positive = TRUE)) {
    cli::cli_abort(
      c("{.arg N} must be a single positive whole number.",
        "x" = "It is the number of cases in each drawn sample."),
      class = "efa_power_input")
  }
  if (!checkmate::test_count(n_datasets, positive = TRUE)) {
    cli::cli_abort(
      c("{.arg n_datasets} must be a single positive whole number.",
        "x" = "It is the number of samples drawn and analysed."),
      class = "efa_power_input")
  }
  if (!checkmate::test_number(recovery_threshold) ||
      recovery_threshold <= 0 || recovery_threshold > 1) {
    cli::cli_abort("{.arg recovery_threshold} must be a single number in (0, 1].",
                   class = "efa_power_bad_threshold")
  }

  # The population is given exactly one way; efa_simulate() owns the deep checks
  # (symmetry, positive-semidefiniteness, Heywood, dimensions), raised as classed
  # efa_simulate_* conditions.
  have_lambda <- !is.null(Lambda)
  have_R <- !is.null(R)
  if (have_lambda == have_R) {
    cli::cli_abort(
      c("Specify the simulation population in exactly one way.",
        "x" = "Provide either {.arg R}, or {.arg Lambda} (with optional {.arg Phi}/{.arg Psi})."),
      class = "efa_power_input")
  }

  # k_true (the true factor count) comes from Lambda for a factor-model population --
  # a conflicting `k` is an error -- and must be given for a bare correlation matrix.
  # `p` is likewise read off the population; a conflicting one is an error rather than
  # silently replaced, so the two dimension arguments are policed the same way.
  if (!is.null(p)) checkmate::assert_count(p, positive = TRUE)
  if (have_lambda) {
    Lambda <- as.matrix(Lambda)
    k_true <- ncol(Lambda)
    if (!is.null(k) && k != k_true) {
      cli::cli_abort(
        c("{.arg k} ({k}) does not match the number of columns of {.arg Lambda} ({k_true}).",
          "i" = "With a factor-model population the true factor count is {.code ncol(Lambda)}; leave {.arg k} unset."),
        class = "efa_power_bad_k")
    }
    pop_arg <- "Lambda"
    p_true <- nrow(Lambda)
  } else {
    if (is.null(k)) {
      cli::cli_abort(
        c("{.arg k} (the true number of factors) is required with an {.arg R} population.",
          "i" = "A bare correlation matrix carries no factor count."),
        class = "efa_power_missing_k")
    }
    checkmate::assert_count(k, positive = TRUE)
    k_true <- as.integer(k)
    pop_arg <- "R"
    p_true <- nrow(as.matrix(R))
  }
  if (!is.null(p) && p != p_true) {
    cli::cli_abort(
      c("{.arg p} ({p}) does not match the number of rows of {.arg {pop_arg}} ({p_true}).",
        "i" = "In simulation mode the number of variables comes from the population; leave {.arg p} unset."),
      class = "efa_power_bad_p")
  }
  p <- p_true

  # Structure recovery aligns the fitted loadings against the population loadings, so
  # it needs a factor-model population; a bare R carries no loadings to recover.
  has_recovery <- have_lambda

  # Recovery-fit rotation: resolved against efa_fit()'s own choices before any data is
  # drawn, so a value it would reject costs milliseconds rather than a full run whose
  # recovery, convergence, and Heywood columns then come back empty. The choices are
  # read off efa_fit()'s formal rather than restated, so the two cannot drift apart.
  # (A test that replaces efa_fit() with a stub replaces its formals too, so such a
  # stub must either carry a `rotation` formal or leave `rotation` unset here.)
  if (!is.null(rotation)) {
    rotation <- .match_arg_ci(rotation, eval(formals(efa_fit)$rotation))
  }

  # Match the population when unset -- varimax for orthogonal factors, promax for
  # oblique ones. Both are deterministic, so the run stays reproducible; a single
  # factor has no rotation.
  if (is.null(rotation)) {
    oblique_pop <- have_lambda && !is.null(Phi) &&
      !isTRUE(all.equal(unname(as.matrix(Phi)), diag(k_true)))
    rotation <- if (oblique_pop) "promax" else "varimax"
  }
  if (k_true == 1L) rotation <- "none"

  # One seed makes the whole call reproducible and independent of the worker count,
  # and the caller's random-number stream is left untouched (mirrors efa_simulate()).
  if (!is.null(seed)) {
    checkmate::assert_int(seed)
    .set_local_seed(seed)
  }

  # Draw the replicate datasets via efa_simulate(). seed = NULL inherits (and advances)
  # the umbrella stream above, decorrelating the data draw from the analysis draw while
  # staying deterministic; efa_simulate()'s own future.seed = TRUE keeps the datasets
  # worker-count independent.
  sim <- efa_simulate(N = N, Lambda = Lambda, Phi = Phi, Psi = Psi, R = R,
                      model_error = model_error, target_rmsea = target_rmsea,
                      target_cfi = target_cfi, n_datasets = n_datasets, seed = NULL)
  datasets <- sim$data
  if (!is.list(datasets)) datasets <- list(datasets)

  # Shared retention control list: the N_FACTORS() defaults with the sample size `N`
  # and `estimator` threaded in (`gof` follows `estimator` inside
  # .n_factors_ctl(): PAF supports only the CAF), both fixed across replicates.
  ctl <- .n_factors_ctl(N = N, estimator = estimator)

  # Analyse every replicate; future.seed = TRUE binds each to its own reproducible
  # L'Ecuyer stream, so criteria that simulate internally (PARALLEL, NEST, CD) are
  # reproducible and independent of the number of workers.
  # Recovery aligns against the population loadings, so pass `Lambda` (already NULL for
  # an R population, which disables recovery inside .efa_power_analyze_one).
  per_rep <- future.apply::future_lapply(
    datasets, .efa_power_analyze_one,
    Lambda = Lambda, k_true = k_true, criteria = criteria, ctl = ctl,
    estimator = estimator, rotation = rotation, future.seed = TRUE)

  # --- Aggregate over replicates ---
  # Hit-rate: union the criterion/variant keys (a criterion can fail on a replicate,
  # or emit several variants), then score each against k_true over the replicates
  # where it produced a suggestion.
  # unlist() drops the NULLs a failed criterion contributes, and returns NULL when every
  # criterion failed on every replicate, so coerce to a character vector before matching.
  observed <- as.character(unique(unlist(lapply(per_rep, function(r) names(r$n_hat)))))
  # A criterion that failed on *every* replicate contributes no key at all, so without
  # this it would silently vanish from the results while `settings$criteria` still records
  # the request -- worst for the internally-simulating criteria on a hard population,
  # where a long run would quietly drop the criterion the user was paying for. Seed the
  # missing ids so each gets a row with n_valid = 0 (and, via the guard below, an NA
  # hit-rate), and say so. A key is a criterion's own id or `id_variant`
  # (see .retention_key()).
  failed <- criteria[!vapply(criteria, function(id) {
    any(observed == id | startsWith(observed, paste0(id, "_")))
  }, logical(1))]
  if (length(failed)) {
    cli::cli_warn(
      c("{cli::qty(failed)}Criteri{?on/a} {.val {failed}} produced no suggestion on any of the {n_datasets} replicate{?s}.",
        "x" = "{cli::qty(failed)}{?It is/They are} reported with a missing hit-rate over zero valid replicates.",
        "i" = "Such a criterion errored or was undecided throughout; a larger {.arg N}, a different {.arg estimator}, or a less demanding population may let it decide."),
      class = "efa_power_criterion_failed")
  }
  # `criteria` always holds at least one id, so seeding the failed ones leaves `keys`
  # non-empty even when no criterion decided anywhere: that case now runs the ordinary
  # path and comes out as n_valid = 0 with an NA hit-rate.
  keys <- c(observed, failed)
  # Extract this replicate's suggestion for every key; a replicate where a criterion
  # failed (its key absent) contributes NA there, and one where all of them failed has a
  # NULL `n_hat` whose indexing gives the wrong length.
  pull <- function(r) {
    v <- r$n_hat[keys]
    if (length(v) != length(keys)) v <- rep(NA_real_, length(keys))
    v
  }
  hit_mat <- matrix(
    vapply(per_rep, pull, numeric(length(keys))),
    nrow = length(keys), dimnames = list(keys, NULL))
  hit_n_valid <- rowSums(!is.na(hit_mat))
  hit_hits <- rowSums(hit_mat == k_true, na.rm = TRUE)
  hit_rate <- ifelse(hit_n_valid > 0, hit_hits / hit_n_valid, NA_real_)
  names(hit_rate) <- keys
  hits <- data.frame(criterion = keys, n_valid = hit_n_valid, hits = hit_hits,
                     hit_rate = hit_rate, row.names = NULL)

  rec_min <- vapply(per_rep, function(r) r$rec_min, numeric(1))
  rec_mean <- vapply(per_rep, function(r) r$rec_mean, numeric(1))
  converged <- vapply(per_rep, function(r) r$converged, logical(1))
  heywood <- vapply(per_rep, function(r) r$heywood, logical(1))
  fit_ok <- vapply(per_rep, function(r) r$fit_ok, logical(1))
  fit_error <- vapply(per_rep, function(r) r$fit_error, character(1))

  recovery <- if (has_recovery) {
    n_rec <- sum(!is.na(rec_min))
    list(
      min_rate = if (n_rec > 0) mean(rec_min >= recovery_threshold, na.rm = TRUE) else NA_real_,
      mean_rate = if (n_rec > 0) mean(rec_mean >= recovery_threshold, na.rm = TRUE) else NA_real_,
      threshold = recovery_threshold, n_valid = n_rec)
  }

  # `fit_rate` is over all replicates; the convergence and Heywood rates are conditional
  # on a fit having completed (denominator `n_fit_ok`), so a replicate whose fit errored
  # is not silently counted as a converged, non-Heywood solution.
  n_fit_ok <- sum(fit_ok)
  convergence <- list(
    n_datasets = as.integer(n_datasets), n_fit_ok = n_fit_ok,
    n_converged = sum(converged), n_heywood = sum(heywood),
    fit_rate = mean(fit_ok),
    convergence_rate = if (n_fit_ok > 0) sum(converged) / n_fit_ok else NA_real_,
    heywood_rate = if (n_fit_ok > 0) sum(heywood) / n_fit_ok else NA_real_)

  # Mirrors the criterion-failure warning above. A systematic failure -- a rotation or
  # estimator that cannot fit this population, a population the fit cannot handle -- is
  # otherwise indistinguishable from per-replicate numerical trouble, and when every fit
  # fails the whole recovery/convergence half of the report is NA with no stated cause.
  # The replicate's own error message is carried through so the user does not have to
  # reproduce the failure to see it. The completed retention results are still returned.
  n_fit_failed <- as.integer(n_datasets) - n_fit_ok
  if (n_fit_failed > 0L) {
    first_failed <- which(!fit_ok)[1L]
    first_cause <- fit_error[first_failed]
    # An R population has no loadings to recover, so the same fit is only the source of
    # the convergence and Heywood rates there; name it for what it did on this run.
    what_fit <- if (has_recovery) "recovery fit" else "replicate fit"
    affected <- if (has_recovery) {
      "Structure recovery, the convergence rate, and the Heywood rate"
    } else {
      "The convergence rate and the Heywood rate"
    }
    cli::cli_warn(
      c("The {what_fit} failed on {n_fit_failed} of the {n_datasets} replicate{?s}.",
        "x" = if (n_fit_ok == 0L) {
          "{affected} are missing."
        } else {
          "{affected} are over the {n_fit_ok} fit{?s} that completed."
        },
        "i" = "Replicate {first_failed} failed with: {first_cause}",
        "i" = "Check {.arg rotation} ({.val {rotation}}) and {.arg estimator} ({.val {estimator}}) against what {.fn efa_fit} accepts."),
      class = "efa_power_fit_failed")
  }

  settings <- list(mode = "simulation", N = as.integer(N),
                   n_datasets = as.integer(n_datasets), p = p, k = k_true,
                   criteria = criteria, estimator = estimator, rotation = rotation,
                   recovery_threshold = recovery_threshold,
                   has_recovery = has_recovery, model_error = model_error,
                   target_rmsea = target_rmsea, target_cfi = target_cfi, seed = seed)

  structure(
    list(
      hit_rate = hit_rate, hits = hits, recovery = recovery,
      convergence = convergence,
      replicates = list(n_hat = t(hit_mat), rec_min = rec_min, rec_mean = rec_mean,
                        converged = converged, heywood = heywood, fit_ok = fit_ok,
                        fit_error = fit_error),
      k_true = k_true, model_error = sim$model_error, settings = settings),
    class = "efa_power")
}

# Analyse a single simulated dataset: (i) run each requested retention criterion
# straight off the registry and record its suggested factor count, (ii) fit the
# k_true-factor EFA once (its convergence code and Heywood flag), and (iii) align the
# fitted loadings to the population loadings and read the matched-factor Tucker
# congruences (recovery). Returns a per-replicate record for aggregation.
.efa_power_analyze_one <- function(dat, Lambda, k_true, criteria, ctl, estimator,
                                   rotation) {

  # Retention: the registry funs take raw data (needs_raw criteria) or a correlation
  # matrix; a criterion that errors on this replicate contributes NA.
  Rmat <- stats::cor(dat)
  n_hat <- unlist(lapply(criteria, function(id) {
    entry <- .retention_registry[[id]]
    xarg <- if (isTRUE(entry$needs_raw)) dat else Rmat
    # Muffle the criteria's routine progress/informational output (e.g. "computing
    # correlations from the raw data"). A criterion that errors contributes nothing on
    # this replicate (NULL is dropped by unlist below), lowering its `n_valid` -- rather
    # than a bare-id NA that would not match its variant-named success key.
    out <- suppressMessages(suppressWarnings(
      tryCatch(entry$fun(xarg, ctl), error = function(e) NULL)))
    if (is.null(out)) return(NULL)
    # Name the suggestion(s) exactly as N_FACTORS() does.
    .retention_key(id, out$n_factors)
  }))

  # k_true-factor fit (once): the source of the convergence/Heywood rate and the
  # loadings recovery aligns against. EFA's own Heywood/non-convergence warnings are
  # muffled here; the flags are read off the returned object instead. A failure keeps
  # its cause, so the aggregation can report why the fit did not complete instead of
  # leaving the user with a column of NAs.
  fit <- suppressMessages(suppressWarnings(tryCatch(
    efa_fit(dat, n_factors = k_true, estimator = estimator, rotation = rotation),
    error = function(e) e)))
  fit_error <- NA_character_
  if (inherits(fit, "error")) {
    fit_error <- conditionMessage(fit)
    fit <- NULL
  }
  fit_ok <- !is.null(fit)
  converged <- isTRUE(fit_ok && fit$convergence == 0)
  heywood <- isTRUE(fit_ok && length(fit$heywood) > 0)

  # Structure recovery: permute/sign-match the fitted loadings to the population
  # loadings, then take the matched-factor congruences (the diagonal). A degenerate
  # (near-zero) factor makes the congruence undefined -- caught as NA.
  rec_min <- NA_real_
  rec_mean <- NA_real_
  if (!is.null(Lambda) && fit_ok) {
    loadings <- if (rotation == "none") {
      .change_class(fit$unrot_loadings, "matrix")
    } else {
      .change_class(fit$rot_loadings, "matrix")
    }
    matched <- tryCatch({
      aligned <- .align_solution(Lambda, loadings)
      diag(.tucker_congruence(Lambda, aligned$loadings))
    },
    efa_zero_column = function(e) NULL,
    efa_undefined_congruence = function(e) NULL,
    error = function(e) NULL)
    if (!is.null(matched)) {
      rec_min <- min(matched)
      rec_mean <- mean(matched)
    }
  }

  list(n_hat = n_hat, converged = converged, heywood = heywood, fit_ok = fit_ok,
       fit_error = fit_error, rec_min = rec_min, rec_mean = rec_mean)
}

#' Print and format an efa_power object
#'
#' `print()` turns an [efa_power()] result into a short report, and `format()`
#' builds the same report as a character vector (`print()` is
#' `cat(format(x), sep = "\n")`).
#'
#' For an RMSEA-mode object, the report has a header naming the test, the null and
#' alternative hypotheses with the significance level and degrees of freedom, the
#' headline result (the power at the sample size, or the required sample size for
#' the target power), and the critical value and noncentrality parameters.
#'
#' For a simulation-mode object, the report instead has the population and design,
#' the retention hit-rate per criterion, the structure-recovery rate, and the
#' convergence and Heywood-case rate.
#'
#' The lines follow the active console theme, so they print as plain text when
#' colours are disabled -- for example when captured into a file, or stripped with
#' [cli::ansi_strip()].
#'
#' @param x An object of class `efa_power` (output from [efa_power()]).
#' @param digits Integer. The number of decimal places the reported values are
#'   rounded to. Default is 3.
#' @param ... Not used; for consistency with the generic.
#'
#' @returns `print()` returns its argument `x` invisibly. `format()` returns a
#'   character vector with the report lines.
#'
#' @family power analysis
#'
#' @export
#'
#' @method print efa_power
#'
#' @examples
#' pw <- efa_power(df = 100, N = 200)
#' pw
#'
#' # format() returns the same lines as a character vector:
#' writeLines(format(pw))
#'
print.efa_power <- function(x, digits = 3, ...) {
  cat(format(x, digits = digits, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.efa_power
#' @export
#' @method format efa_power
format.efa_power <- function(x, digits = 3, ...) {

  # Simulation-mode results carry a different payload (hit-rate, recovery,
  # convergence) and are rendered by their own formatter.
  if (identical(x$settings$mode, "simulation")) {
    return(.format_efa_power_simulation(x, digits = digits))
  }

  s <- x$settings
  test_lbl <- if (s$type == "close") "Test of close fit" else "Test of not-close fit"
  # <= bounds the null for close fit, >= for not-close fit.
  cmp <- if (s$type == "close") "\u2264" else "\u2265"
  e0 <- .efa_num(s$eps0, digits = digits, pad = FALSE)
  e1 <- .efa_num(s$eps1, digits = digits, pad = FALSE)
  a <- .efa_num(s$alpha, digits = digits, pad = FALSE)
  pw <- .efa_num(x$power, digits = digits, pad = FALSE)
  crit <- .efa_num(x$crit, digits = digits, pad = FALSE)
  n0 <- .efa_num(x$ncp[["H0"]], digits = digits, pad = FALSE)
  n1 <- .efa_num(x$ncp[["H1"]], digits = digits, pad = FALSE)

  cli::cli_format_method({
    cli::cli_text("")
    cli::cli_rule(left = "{.strong RMSEA power analysis}")
    cli::cli_text("")

    cli::cli_text("{test_lbl}: H0 RMSEA {cmp} {e0} vs. H1 RMSEA = {e1}.")
    if (s$group > 1) {
      cli::cli_text("alpha = {a} \u00b7 df = {s$df} \u00b7 groups = {s$group}")
    } else {
      cli::cli_text("alpha = {a} \u00b7 df = {s$df}")
    }
    cli::cli_text("")

    # With more than one group `N` is the total across them, which is the quantity a
    # multiple-group study is most likely to misread as per-group; name it and give the
    # per-group figure outright.
    # A sample size, not a coefficient: format it as a count rather than through
    # .efa_num(), whose fixed `digits` would render 100 as "100.000".
    # A solved `N` is a multiple of `group`, so the per-group figure is a whole number.
    # A supplied `N` is the caller's own total and is never changed: when it does not
    # divide by `group`, the fraction the noncentrality uses is shown as it is.
    # An object serialized before `N_per_group` existed carries `group` but not the field;
    # its `N` was per-group back then, so no note can be derived from it -- omit the note
    # rather than erroring in round(NULL) or dividing a per-group N by `group` again.
    n_note <- if (s$group > 1 && !is.null(x$N_per_group)) {
      paste0(" (total; ", format(round(x$N_per_group, 1)), " per group)")
    } else {
      ""
    }
    if (x$solve_for == "power") {
      cli::cli_text("{.strong Power = {pw}} at N = {x$N}{n_note}.")
    } else {
      tgt <- .efa_num(s$power, digits = digits, pad = FALSE)
      cli::cli_text("{.strong Required N = {x$N}}{n_note} for a power of {tgt} (achieved {pw}).")
    }
    cli::cli_text(
      "Critical value \u03c7\u00b2({s$df}) = {crit} \u00b7 noncentrality H0 = {n0}, H1 = {n1}.")
  })
}

# Report lines for a simulation-mode efa_power object: the population and design, the
# retention hit-rate per criterion, the structure-recovery rate, and the
# convergence/Heywood rate. Numbers are pre-formatted into plain strings so cli's
# `{...}` interpolation does not mistake `.efa_num(...)` for an inline style.
.format_efa_power_simulation <- function(x, digits = 3) {

  s <- x$settings
  conv <- x$convergence

  hit_lines <- vapply(seq_len(nrow(x$hits)), function(i) {
    h <- x$hits[i, ]
    paste0(h$criterion, ": ", .efa_num(h$hit_rate, digits = digits, pad = FALSE),
           " (n = ", h$n_valid, ")")
  }, character(1))

  conv_rate <- .efa_num(conv$convergence_rate, digits = digits, pad = FALSE)
  heywood_rate <- .efa_num(conv$heywood_rate, digits = digits, pad = FALSE)
  fit_rate <- .efa_num(conv$fit_rate, digits = digits, pad = FALSE)

  cli::cli_format_method({
    cli::cli_text("")
    cli::cli_rule(left = "{.strong EFA power simulation}")
    cli::cli_text("")

    cli::cli_text(
      "{s$p} variable{?s} \u00b7 {x$k_true} factor{?s} \u00b7 N = {s$N} \u00b7 {s$n_datasets} dataset{?s}")
    cli::cli_text("Estimation: {s$estimator} \u00b7 rotation: {s$rotation}")
    if (!is.null(x$model_error)) {
      me <- x$model_error
      me_rmsea <- .efa_num(me$rmsea, digits = digits, pad = FALSE)
      me_cfi <- .efa_num(me$cfi, digits = digits, pad = FALSE)
      cli::cli_text(
        "Model error ({me$method}): RMSEA = {me_rmsea} \u00b7 CFI = {me_cfi}")
    } else {
      cli::cli_text(
        "Model error: none. The population is exact, so the hit-rate and recovery are optimistic; set {.arg target_rmsea} for realism.")
    }
    cli::cli_text("")

    cli::cli_text("{.strong Retention hit-rate} P(k-hat = {x$k_true})")
    # Every requested criterion now gets a row, so a fresh object always has lines here;
    # an object serialized before that carries no row for a criterion that failed on every
    # replicate, and prints this instead of an empty section (as for `N_per_group` above).
    if (length(hit_lines) > 0) {
      cli::cli_ul(hit_lines)
    } else {
      cli::cli_text("No criterion produced a suggestion on any replicate.")
    }
    cli::cli_text("")

    if (!is.null(x$recovery)) {
      r <- x$recovery
      thr <- .efa_num(r$threshold, digits = digits, pad = FALSE)
      min_rate <- .efa_num(r$min_rate, digits = digits, pad = FALSE)
      mean_rate <- .efa_num(r$mean_rate, digits = digits, pad = FALSE)
      # min_rate/mean_rate are the proportions of replicates clearing the threshold,
      # not congruences. The median raw congruence is added alongside them so a rate
      # of .000 still says how far off the recovered structure was.
      med_min <- .efa_num(stats::median(x$replicates$rec_min, na.rm = TRUE),
                          digits = digits, pad = FALSE)
      cli::cli_text("{.strong Structure recovery} (Tucker congruence \u2265 {thr})")
      cli::cli_ul(c(
        paste0("recovery rate (min congruence): ", min_rate, " (n = ", r$n_valid, ")"),
        paste0("recovery rate (mean congruence): ", mean_rate, " (n = ", r$n_valid, ")"),
        paste0("median min congruence: ", med_min)))
    } else {
      cli::cli_text(
        "{.strong Structure recovery}: not available (needs a factor-model population).")
    }
    cli::cli_text("")

    cli::cli_text("{.strong Convergence}")
    cli::cli_ul(c(
      paste0("fits completed: ", fit_rate, " (", conv$n_fit_ok, "/", conv$n_datasets, ")"),
      paste0("converged (of completed): ", conv_rate),
      paste0("Heywood cases (of completed): ", heywood_rate)))
  })
}
