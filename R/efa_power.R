#' Power analysis for exploratory factor analysis
#'
#' @description
#' Analyses power for exploratory factor analysis in one of two modes, chosen with `mode`.
#'
#' `mode = "rmsea"` (the default) gives the analytic power of the root mean square error of
#' approximation (RMSEA) tests of close and not-close fit (MacCallum, Browne, & Sugawara,
#' 1996). Give a sample size to get the power of the test, or give a target power to get the
#' necessary sample size.
#'
#' `mode = "simulation"` runs a Monte-Carlo study: it draws `n_datasets` samples from a known
#' population (with [efa_simulate()]), analyses each sample, and reports how frequently the
#' analysis finds the population structure. A call with `Lambda` or `R` and no `mode` uses
#' this mode.
#'
#' @details
#' # RMSEA mode
#'
#' Power increases with the sample size, the degrees of freedom, and the difference between
#' the null and the alternative RMSEA (MacCallum et al., 1996). `type` selects the test:
#' \describe{
#'   \item{`"close"`}{The null hypothesis is close fit (RMSEA \eqn{\le} `eps0`, conventionally
#'     0.05). Power is the probability to reject it if the RMSEA is `eps1` (conventionally
#'     0.08). `eps0` must be smaller than `eps1`.}
#'   \item{`"notclose"`}{The null hypothesis is not-close fit (RMSEA \eqn{\ge} `eps0`). Power is
#'     the probability to reject it if the RMSEA is `eps1` (conventionally 0.01). `eps0` must be
#'     larger than `eps1`.}
#' }
#' The required sample size is the smallest `N` that reaches `power`. The result applies to the
#' RMSEA test only: it does not make sure that the estimates are precise. With many degrees of
#' freedom, the required `N` can be smaller than the number of variables. The sample correlation
#' matrix is then singular (MacCallum et al., 1996, p. 144). If `n_vars` is known, such a sample
#' size gives a warning.
#'
#' With `group > 1`, `N` is the total of `group` groups of equal size, and `df` is that of the
#' whole multiple-group model. The noncentrality parameter is
#' \eqn{(N - G)\, df\, \varepsilon^2 / G} for \eqn{G} groups (Steiger, 1998). Thus, at a given
#' `df`, power depends only on the group size `N / group`, and a required `N` is `group` times
#' the smallest group size that reaches `power`.
#'
#' # Simulation mode
#'
#' [efa_simulate()] draws `n_datasets` samples of size `N` from the population. The true number
#' of factors, `k_true`, is `ncol(Lambda)`, or `n_factors` for a correlation matrix `R`. Without
#' `target_rmsea` or `target_cfi`, a factor-model population fits the model exactly, which gives
#' optimistic results. A misfit target makes the population more realistic (MacCallum, 2003).
#'
#' The hit-rate and the recovery rates are proportions of all `n_datasets` replicates:
#' \describe{
#'   \item{**Hit-rate**}{A criterion in `criteria` suggests `k_true` factors. A replicate where
#'     the criterion fails or gives no suggestion counts as a miss.}
#'   \item{**Structure recovery** (factor-model populations only)}{The `k_true`-factor model is
#'     fitted with [efa_fit()], and its factors are matched to the population loadings. A
#'     replicate is recovered if the smallest (`min`) or the mean (`mean`) Tucker congruence of
#'     the matched factors reaches `recovery_threshold` (Lorenzo-Seva & ten Berge, 2006). A
#'     failed fit counts as not recovered.}
#' }
#' The Monte Carlo standard error of a rate \eqn{p} is \eqn{\sqrt{p(1 - p) / n}} for
#' \eqn{n} replicates (Morris, White, & Crowther, 2019). The convergence and Heywood rates are
#' proportions of the fits that completed. If fits fail, a warning gives their number and the
#' first cause.
#'
#' The replicates are analysed in parallel with \pkg{future.apply}; set a parallel plan with
#' [future::plan()]. With a fixed `seed`, the result does not depend on the number of workers,
#' and your random-number stream does not change.
#'
#' @param mode character. `"rmsea"` (the default) or `"simulation"`. If `mode` is not given,
#'   `Lambda` or `R` selects `"simulation"`. `type`, `eps0`, `eps1`, `df`, `alpha`, `power`, and
#'   `group` apply only to RMSEA mode. The arguments marked *Simulation mode* apply only to
#'   simulation mode. An argument of the other mode is an error.
#' @param type character. The RMSEA test: `"close"` or `"notclose"`. See *Details*.
#' @param eps0 numeric. The null-hypothesis RMSEA. Default is `0.05`.
#' @param eps1 numeric. The alternative-hypothesis RMSEA, at which power is calculated. Default
#'   is `0.08` for `type = "close"` and `0.01` for `type = "notclose"`.
#' @param N numeric. The sample size. In RMSEA mode, the total of all groups: give `N` to get the
#'   power, or leave it `NULL` to get the required `N`. In simulation mode (required), the size
#'   of each sample.
#' @param n_vars numeric. The number of variables. In RMSEA mode, it gives `df` together with
#'   `n_factors`. In simulation mode, it comes from the population; if given, it must agree.
#' @param n_factors numeric. The number of factors. In RMSEA mode, it gives `df` together with
#'   `n_vars`. In simulation mode, the true number of factors: necessary with `R`; with
#'   `Lambda`, if given, it must agree with `ncol(Lambda)`.
#' @param df numeric. The model degrees of freedom. Give `df`, or `n_vars` and `n_factors`,
#'   which give `df = group * ((n_vars - n_factors)^2 - (n_vars + n_factors)) / 2` (a separate
#'   EFA in each group).
#' @param alpha numeric. The significance level. Default is `0.05`.
#' @param power numeric. The target power for the required `N`. Default is `NULL`, which uses
#'   `0.80` if `N` is also `NULL`. Give only one of `N` and `power`.
#' @param group numeric. The number of groups of equal size. Default is `1`. See *Details*.
#' @param Lambda matrix. Simulation mode. The population loadings (variables by factors) of a
#'   factor-model population. Structure recovery is available only with `Lambda`.
#' @param Phi matrix. Simulation mode. The population factor correlations. Only used with
#'   `Lambda`. Default is orthogonal factors.
#' @param Psi numeric or matrix. Simulation mode. The population unique variances. Only used with
#'   `Lambda`. See [efa_simulate()].
#' @param R matrix. Simulation mode. A population correlation matrix, instead of `Lambda`. It
#'   needs `n_factors`.
#' @param n_datasets numeric. Simulation mode. The number of samples. Default is `500`.
#' @param criteria character. Simulation mode. The retention criteria for the hit-rate: any of
#'   `"CD"`, `"EKC"`, `"HULL"`, `"KGC"`, `"MAP"`, `"NEST"`, `"PARALLEL"`, and `"SMT"` (see
#'   [efa_retain()]). Default is `c("EKC", "MAP")`. The criteria use the [efa_retain()]
#'   defaults and `estimator`. `"CD"` and `"HULL"` are much slower.
#' @param estimator character. Simulation mode. The estimator of the recovery fit and the
#'   retention criteria: `"PAF"` (the default), `"ML"`, or `"ULS"`.
#' @param rotation character. Simulation mode. The rotation of the recovery fit (see
#'   [efa_fit()]). Default is `"oblimin"`, as in [efa_fit()]. Recovery only reorders and
#'   reflects the factors. Thus, `"none"` with more than one factor understates recovery.
#' @param recovery_threshold numeric. Simulation mode. The Tucker congruence that a replicate
#'   must reach to be recovered. Default is `0.95`, at which Lorenzo-Seva and ten Berge (2006)
#'   consider two factors equal.
#' @param model_error character. Simulation mode. The [efa_simulate()] method for model error:
#'   `"TKL"` (the default; it adds minor common factors), `"CB"`, `"WB"`, or `"none"`. It applies
#'   only with `target_rmsea` or `target_cfi`, which need `Lambda`. Note that [efa_simulate()]
#'   has the default `"CB"`.
#' @param target_rmsea numeric. Simulation mode. The population RMSEA. Default is `NULL` (no
#'   model error). See [efa_simulate()].
#' @param target_cfi numeric. Simulation mode. The population CFI (only with `"TKL"`). Default
#'   is `NULL`. See [efa_simulate()].
#' @param seed numeric. Simulation mode. An optional seed for a reproducible result. Default is
#'   `NULL`.
#' @param p,k `r lifecycle::badge("deprecated")` Use `n_vars` and `n_factors`.
#'
#' @returns An object of class `efa_power`. For `mode = "rmsea"`, a list with:
#' \item{power}{The power at `N` (for a required `N`, at least the target power).}
#' \item{N}{The total sample size: the given `N`, or the required `N` (a multiple of `group`).}
#' \item{N_per_group}{The group size `N / group`.}
#' \item{crit}{The critical chi-square value.}
#' \item{ncp}{The noncentrality parameters under the null (`H0`) and the alternative (`H1`).}
#' \item{solve_for}{`"power"` or `"N"`: the calculated quantity.}
#' \item{settings}{The inputs: `mode`, `type`, `eps0`, `eps1`, `df`, `n_vars`, `n_factors`,
#'   `alpha`, `group`, and the target `power` (`NULL` if `solve_for` is `"power"`).}
#'
#' For `mode = "simulation"`, a list with:
#' \item{hit_rate}{The hit-rate per criterion, and per variant if a criterion has more than one.}
#' \item{hits}{A data frame with one row per criterion or variant: `criterion`, the number of
#'   replicates
#'   with a suggestion (`n_valid`), the number of hits (`hits`), `hit_rate`
#'   (`hits / n_datasets`), and its Monte Carlo standard error (`mcse`).}
#' \item{recovery}{For `Lambda`, a list with the recovery rates (`min_rate`, `mean_rate`), the
#'   `threshold`, and the number of replicates with a congruence (`n_valid`); `NULL` for `R`.
#'   Non-converged and Heywood solutions are included.}
#' \item{convergence}{A list with `n_datasets`, the number of completed fits (`n_fit_ok`), the
#'   number of those that converged (`n_converged`) or have a Heywood case (`n_heywood`), and
#'   the rates `fit_rate` (of all datasets), `convergence_rate`, and `heywood_rate` (of the
#'   completed fits).}
#' \item{replicates}{The per-replicate values: the suggested numbers of factors (`n_hat`), the
#'   matched congruences (`rec_min`, `rec_mean`), the flags `converged`, `heywood`, and
#'   `fit_ok`, and `fit_error`, the error message of a failed fit (`NA` otherwise).}
#' \item{k_true}{The true number of factors.}
#' \item{model_error}{The [efa_simulate()] model-error record, or `NULL`.}
#' \item{settings}{The simulation inputs.}
#'
#' @references
#' Lorenzo-Seva, U., & ten Berge, J. M. F. (2006). Tucker's congruence coefficient as a
#' meaningful index of factor similarity. *Methodology, 2*(2), 57-64.
#' \doi{10.1027/1614-2241.2.2.57}
#'
#' MacCallum, R. C. (2003). 2001 Presidential Address: Working with imperfect models.
#' *Multivariate Behavioral Research, 38*(1), 113-139. \doi{10.1207/S15327906MBR3801_5}
#'
#' MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power analysis and
#' determination of sample size for covariance structure modeling. *Psychological
#' Methods, 1*(2), 130-149. \doi{10.1037/1082-989X.1.2.130}
#'
#' Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation studies to evaluate
#' statistical methods. *Statistics in Medicine, 38*(11), 2074-2102. \doi{10.1002/sim.8086}
#'
#' Steiger, J. H. (1998). A note on multiple sample extensions of the RMSEA fit index.
#' *Structural Equation Modeling, 5*(4), 411-419. \doi{10.1080/10705519809540115}
#'
#' @seealso [efa_simulate()] draws the samples of simulation mode. [efa_retain()] implements
#'   the retention criteria.
#'
#' @family power analysis
#'
#' @export
#'
#' @examples
#' # Power of the test of close fit at N = 200 for a 100-df model
#' efa_power(df = 100, N = 200)
#'
#' # df from the number of variables and factors
#' efa_power(n_vars = 20, n_factors = 3, N = 200)
#'
#' # Required sample size for 80% power
#' efa_power(df = 100, power = 0.80)
#'
#' # Test of not-close fit
#' efa_power(df = 100, N = 200, type = "notclose")
#'
#' \donttest{
#' # Simulation mode: hit-rate and structure recovery for a known three-factor
#' # population at N = 300
#' efa_power("simulation", Lambda = population_models$loadings$baseline,
#'           Phi = population_models$phis_3$moderate, N = 300,
#'           n_datasets = 50, criteria = c("EKC", "MAP"), seed = 42)
#' }
#'
efa_power <- function(mode = c("rmsea", "simulation"),
                      type = c("close", "notclose"),
                      eps0 = NULL, eps1 = NULL, N = NULL, n_vars = NULL, n_factors = NULL,
                      df = NULL, alpha = 0.05, power = NULL, group = 1,
                      Lambda = NULL, Phi = NULL, Psi = NULL, R = NULL,
                      n_datasets = 500, criteria = c("EKC", "MAP"),
                      estimator = "PAF", rotation = "oblimin", recovery_threshold = 0.95,
                      model_error = c("TKL", "CB", "WB", "none"),
                      target_rmsea = NULL, target_cfi = NULL, seed = NULL,
                      p = lifecycle::deprecated(), k = lifecycle::deprecated()) {

  if (lifecycle::is_present(p)) {
    lifecycle::deprecate_warn("1.2.0", "efa_power(p)", "efa_power(n_vars)")
    n_vars <- p
  }
  if (lifecycle::is_present(k)) {
    lifecycle::deprecate_warn("1.2.0", "efa_power(k)", "efa_power(n_factors)")
    n_factors <- k
  }

  # A population selects simulation mode; an argument of the other mode is an error, so that it
  # is never silently ignored.
  if (missing(mode) && (!is.null(Lambda) || !is.null(R))) mode <- "simulation"
  mode <- .match_arg_ci(mode)
  rmsea_args <- c("type", "eps0", "eps1", "df", "alpha", "power", "group")
  sim_args <- c("Lambda", "Phi", "Psi", "R", "n_datasets", "criteria", "estimator", "rotation",
                "recovery_threshold", "model_error", "target_rmsea", "target_cfi", "seed")
  other <- if (mode == "rmsea") "simulation" else "rmsea"
  wrong <- intersect(names(match.call())[-1L], if (mode == "rmsea") sim_args else rmsea_args)
  # An explicit NULL is the same as an argument that is not given.
  wrong <- wrong[!vapply(mget(wrong), is.null, logical(1))]
  if (length(wrong)) {
    cli::cli_abort(
      c("{.arg {wrong}} appl{?ies/y} only to {.code mode = \"{other}\"}.",
        "i" = "{cli::qty(wrong)}Remove {?it/them}, or set {.arg mode}."),
      class = "efa_power_wrong_mode")
  }

  if (mode == "simulation") {
    return(.efa_power_simulation(
      N = N, Lambda = Lambda, Phi = Phi, Psi = Psi, R = R, n_vars = n_vars,
      n_factors = n_factors, n_datasets = n_datasets, criteria = criteria,
      estimator = estimator, rotation = rotation, recovery_threshold = recovery_threshold,
      model_error = model_error, target_rmsea = target_rmsea, target_cfi = target_cfi,
      seed = seed))
  }

  type <- .match_arg_ci(type)

  # Conventional defaults: a null of .05 against a worse (.08) or a better (.01) alternative.
  if (is.null(eps0)) eps0 <- 0.05
  if (is.null(eps1)) eps1 <- if (type == "close") 0.08 else 0.01
  checkmate::assert_number(eps0, lower = 0, finite = TRUE)
  checkmate::assert_number(eps1, lower = 0, finite = TRUE)
  checkmate::assert_count(group, positive = TRUE)
  if (!checkmate::test_number(alpha) || alpha <= 0 || alpha >= 1) {
    cli::cli_abort("{.arg alpha} must be a single number strictly between 0 and 1.",
                   class = "efa_power_bad_alpha")
  }

  if (!is.null(n_vars)) checkmate::assert_count(n_vars, positive = TRUE)
  if (!is.null(n_factors)) checkmate::assert_count(n_factors, positive = TRUE)
  if (is.null(df)) {
    if (is.null(n_vars) || is.null(n_factors)) {
      cli::cli_abort(
        c("The model degrees of freedom are missing.",
          "i" = "Give {.arg df}, or both {.arg n_vars} and {.arg n_factors}."),
        class = "efa_power_missing_df")
    }
    # A separate EFA in each group: the multiple-group model has `group` times the df.
    df <- group * .efa_df(n_vars, n_factors)
  } else {
    checkmate::assert_number(df, finite = TRUE)
  }
  if (df <= 0) {
    cli::cli_abort(
      c("The model has {df} degree{?s} of freedom; a positive value is necessary.",
        "i" = "Use more variables or fewer factors."),
      class = "efa_power_bad_df")
  }

  if (eps0 == eps1) {
    cli::cli_abort(
      c("{.arg eps0} and {.arg eps1} are both {eps0}, so the test has no power to detect a difference.",
        "i" = "Give them different RMSEA values (for example {.val 0.05} vs {.val 0.08} for a close-fit test)."),
      class = "efa_power_unreachable")
  }
  # The rejection region follows `type`, so the alternative must be on that side of the null.
  if (if (type == "close") eps0 > eps1 else eps0 < eps1) {
    side <- if (type == "close") "smaller" else "larger"
    cli::cli_abort(
      c("For a {.val {type}}-fit test, {.arg eps0} ({eps0}) must be {side} than {.arg eps1} ({eps1}).",
        "i" = "Swap the two values, or change {.arg type}."),
      class = "efa_power_wrong_side")
  }

  if (!is.null(N) && !is.null(power)) {
    cli::cli_abort(
      c("Both {.arg N} and {.arg power} were given, so there is nothing to calculate.",
        "i" = "Leave {.arg N} as {.code NULL} to get the required sample size, or {.arg power} as {.code NULL} to get the power at {.arg N}."),
      class = "efa_power_overdetermined")
  }

  if (is.null(N)) {
    solve_for <- "N"
    if (is.null(power)) power <- 0.80
    if (!checkmate::test_number(power) || power <= 0 || power >= 1) {
      cli::cli_abort("{.arg power} must be a single number strictly between 0 and 1.",
                     class = "efa_power_bad_power")
    }
    N <- .efa_power_solve_N(df, eps0, eps1, alpha, group, type, power)
  } else {
    solve_for <- "power"
    checkmate::assert_count(N, positive = TRUE)
    if (N < group) {
      cli::cli_abort("{.arg N} ({N}) must be at least the number of groups ({group}).",
                     class = "efa_power_bad_n")
    }
    N <- as.integer(N)
  }

  # The power theory needs more cases than variables in each group (MacCallum et al., 1996,
  # p. 144); below that, the sample correlation matrix is singular.
  if (!is.null(n_vars) && N / group <= n_vars) {
    n_lbl <- if (group > 1) "sample size per group" else "sample size"
    cli::cli_warn(
      c("The {n_lbl} ({N / group}) is not larger than the number of variables ({n_vars}).",
        "x" = "The sample correlation matrix is then singular, and the RMSEA power does not apply.",
        "i" = "Use {.code mode = \"simulation\"} to examine the sample size for the factor analysis."),
      class = "efa_power_small_n")
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
                      n_vars = n_vars, n_factors = n_factors, alpha = alpha, group = group,
                      power = power)
    ),
    class = "efa_power"
  )
}

# Power of the RMSEA test of (not-)close fit at one total sample size N (MacCallum, Browne, &
# Sugawara, 1996). The noncentrality parameter is (N - group) * df * eps^2 / group (Steiger,
# 1998). Close fit rejects in the upper tail of the null distribution, not-close fit in the
# lower tail. Returns the power, the critical value, and both noncentrality parameters.
.efa_power_rmsea <- function(N, df, eps0, eps1, alpha, group, type) {
  ncp0 <- (N - group) * df * eps0^2 / group
  ncp1 <- (N - group) * df * eps1^2 / group
  lower <- type != "close"
  crit <- .qchisq_nc(alpha, df, ncp0, lower)
  list(power = .pchisq_nc(crit, df, ncp1, lower), crit = crit, ncp0 = ncp0, ncp1 = ncp1)
}

# Noncentral chi-square quantile. ?pchisq warns that a noncentrality above about 1e5 can give
# inaccurate results (qchisq() then returns wrong quantiles). Above 1e5, the three-moment
# approximation a + b * chisq(f) of Pearson (1959) is used; its probabilities are then within
# 4e-7 of the exact values.
.qchisq_nc <- function(p, df, ncp, lower) {
  # Below 1e5, qchisq() can give a convergence note with a correct quantile.
  if (ncp <= 1e5) return(suppressWarnings(stats::qchisq(p, df, ncp = ncp, lower.tail = lower)))
  m <- .pearson_chisq(df, ncp)
  m$a + m$b * stats::qchisq(p, m$f, lower.tail = lower)
}

# Noncentral chi-square probability (see .qchisq_nc()). Far in a tail, pchisq() can give a
# precision note for a probability near 0.
.pchisq_nc <- function(q, df, ncp, lower) {
  if (ncp <= 1e5) return(suppressWarnings(stats::pchisq(q, df, ncp = ncp, lower.tail = lower)))
  m <- .pearson_chisq(df, ncp)
  stats::pchisq((q - m$a) / m$b, m$f, lower.tail = lower)
}

# Pearson, E. S. (1959), Biometrika, 46, 364: a + b * chisq(f) with the mean, variance, and
# third cumulant of the noncentral chi-square (df + ncp, 2 (df + 2 ncp), 8 (df + 3 ncp)).
.pearson_chisq <- function(df, ncp) {
  list(a = -ncp^2 / (df + 3 * ncp), b = (df + 3 * ncp) / (df + 2 * ncp),
       f = (df + 2 * ncp)^3 / (df + 3 * ncp)^2)
}

# Smallest total sample size that reaches `target` power with `group` groups of equal size.
# Power depends only on the group size m = N / group and increases with it, so double an upper
# bound for m until it reaches the target, then bisect. `cap` bounds the total, because an
# alternative near the null can need a very large sample.
.efa_power_solve_N <- function(df, eps0, eps1, alpha, group, type, target, cap = 1e7) {
  powfun <- function(m) .efa_power_rmsea(m, df, eps0, eps1, alpha, 1, type)$power
  cap_m <- floor(cap / group)

  lo <- 1
  if (powfun(lo) >= target) return(as.integer(group))
  hi <- 2
  while (powfun(hi) < target) {
    if (hi >= cap_m) {
      cli::cli_abort(
        c("No sample size up to {.val {as.integer(cap_m * group)}} reaches a power of {target}.",
          "i" = "The alternative RMSEA can be too near the null, or the target power too high."),
        class = "efa_power_unreached")
    }
    hi <- min(hi * 2, cap_m)
  }
  while (hi - lo > 1) {
    mid <- floor((lo + hi) / 2)
    if (powfun(mid) >= target) hi <- mid else lo <- mid
  }
  as.integer(hi * group)
}


# Simulation mode: draw `n_datasets` samples from a known population and report, over all
# replicates, the hit-rate of each retention criterion, the structure recovery, and the
# convergence of the k_true-factor fit.
.efa_power_simulation <- function(N, Lambda, Phi, Psi, R, n_vars, n_factors, n_datasets,
                                  criteria, estimator, rotation, recovery_threshold,
                                  model_error, target_rmsea, target_cfi, seed) {

  model_error <- .match_arg_ci(model_error, c("TKL", "CB", "WB", "none"))
  estimator <- .match_arg_ci(estimator, c("PAF", "ML", "ULS"))
  # Checked here, so that a rotation that efa_fit() refuses costs no simulation time.
  rotation <- .match_arg_ci(rotation, eval(formals(efa_fit)$rotation))
  # The visual scree test gives no number of factors to score.
  valid_ids <- names(.retention_registry)[
    !vapply(.retention_registry, function(e) isTRUE(e$visual), logical(1))]
  # NULL selects the default, which is on the formal of efa_power().
  if (is.null(criteria)) criteria <- eval(formals(efa_power)$criteria)
  criteria <- .match_arg_ci(criteria, valid_ids, several.ok = TRUE)

  if (is.null(N)) {
    cli::cli_abort(
      c("{.arg N} is necessary in simulation mode.",
        "x" = "It is the size of each sample; no sample size is calculated here.",
        "i" = "Set {.arg N}, e.g. {.code N = 300}, or use {.code mode = \"rmsea\"} to get the required sample size."),
      class = "efa_power_input")
  }
  if (!checkmate::test_count(N, positive = TRUE)) {
    cli::cli_abort(
      c("{.arg N} must be a single positive whole number.",
        "x" = "It is the number of cases in each sample."),
      class = "efa_power_input")
  }
  if (!checkmate::test_count(n_datasets, positive = TRUE)) {
    cli::cli_abort(
      c("{.arg n_datasets} must be a single positive whole number.",
        "x" = "It is the number of samples."),
      class = "efa_power_input")
  }
  if (!checkmate::test_number(recovery_threshold) ||
      recovery_threshold <= 0 || recovery_threshold > 1) {
    cli::cli_abort("{.arg recovery_threshold} must be a single number in (0, 1].",
                   class = "efa_power_bad_threshold")
  }

  # efa_simulate() does the detailed checks of the population.
  have_lambda <- !is.null(Lambda)
  if (have_lambda == !is.null(R)) {
    cli::cli_abort(
      c("Specify the simulation population in exactly one way.",
        "x" = "Give either {.arg R}, or {.arg Lambda} (with optional {.arg Phi}/{.arg Psi})."),
      class = "efa_power_input")
  }
  if (!is.null(n_vars)) checkmate::assert_count(n_vars, positive = TRUE)
  if (have_lambda) {
    Lambda <- as.matrix(Lambda)
    k_true <- ncol(Lambda)
    if (!is.null(n_factors) && n_factors != k_true) {
      cli::cli_abort(
        c("{.arg n_factors} ({n_factors}) does not agree with the number of columns of {.arg Lambda} ({k_true}).",
          "i" = "With {.arg Lambda}, the true number of factors is {.code ncol(Lambda)}; leave {.arg n_factors} unset."),
        class = "efa_power_bad_k")
    }
    pop_arg <- "Lambda"
    p_true <- nrow(Lambda)
  } else {
    if (is.null(n_factors)) {
      cli::cli_abort(
        c("{.arg n_factors} (the true number of factors) is necessary with {.arg R}.",
          "i" = "A correlation matrix does not give the number of factors."),
        class = "efa_power_missing_k")
    }
    checkmate::assert_count(n_factors, positive = TRUE)
    k_true <- as.integer(n_factors)
    pop_arg <- "R"
    p_true <- nrow(as.matrix(R))
  }
  if (!is.null(n_vars) && n_vars != p_true) {
    cli::cli_abort(
      c("{.arg n_vars} ({n_vars}) does not agree with the number of rows of {.arg {pop_arg}} ({p_true}).",
        "i" = "In simulation mode, the number of variables comes from the population; leave {.arg n_vars} unset."),
      class = "efa_power_bad_p")
  }
  if (N <= p_true) {
    cli::cli_abort(
      c("{.arg N} ({N}) must be larger than the number of variables ({p_true}).",
        "x" = "With fewer cases, the correlation matrix of a sample is singular, and no criterion or fit can run."),
      class = "efa_power_input")
  }
  if (k_true == 1L) rotation <- "none"

  if (!is.null(seed)) {
    checkmate::assert_int(seed)
    .set_local_seed(seed)
  }

  # seed = NULL continues the stream set above; efa_simulate() binds each sample to its own
  # stream, so the samples do not depend on the number of workers.
  sim <- efa_simulate(N = N, Lambda = Lambda, Phi = Phi, Psi = Psi, R = R,
                      model_error = model_error, target_rmsea = target_rmsea,
                      target_cfi = target_cfi, n_datasets = n_datasets, seed = NULL)
  datasets <- sim$data
  if (!is.list(datasets)) datasets <- list(datasets)

  ctl <- .n_factors_ctl(N = N, estimator = estimator)

  # future.seed = TRUE gives each replicate its own reproducible stream, which the criteria
  # that simulate (PARALLEL, NEST, CD, HULL) and the random starts of the rotation use.
  per_rep <- future.apply::future_lapply(
    datasets, .efa_power_analyze_one,
    Lambda = Lambda, k_true = k_true, criteria = criteria, ctl = ctl,
    estimator = estimator, rotation = rotation, future.seed = TRUE)

  # Hit-rate. A key is a criterion id or `id_variant` (see .retention_key()). A criterion that
  # failed on every replicate has no key, so its id is added to keep its row.
  observed <- as.character(unique(unlist(lapply(per_rep, function(r) names(r$n_hat)))))
  failed <- criteria[!vapply(criteria, function(id) {
    any(observed == id | startsWith(observed, paste0(id, "_")))
  }, logical(1))]
  if (length(failed)) {
    cli::cli_warn(
      c("{cli::qty(failed)}Criteri{?on/a} {.val {failed}} gave no suggestion on any of the {n_datasets} replicate{?s}.",
        "x" = "{cli::qty(failed)}{?Its/Their} hit-rate is 0.",
        "i" = "A larger {.arg N} or a less demanding population can help."),
      class = "efa_power_criterion_failed")
  }
  keys <- c(observed, failed)
  # A replicate where all criteria failed has a NULL `n_hat`, which gives the wrong length.
  pull <- function(r) {
    v <- r$n_hat[keys]
    if (length(v) != length(keys)) v <- rep(NA_real_, length(keys))
    v
  }
  hit_mat <- matrix(
    vapply(per_rep, pull, numeric(length(keys))),
    nrow = length(keys), dimnames = list(keys, NULL))
  hit_hits <- rowSums(hit_mat == k_true, na.rm = TRUE)
  # A failed or undecided replicate counts as a miss, as it would in a study.
  hit_rate <- hit_hits / n_datasets
  hits <- data.frame(criterion = keys, n_valid = rowSums(!is.na(hit_mat)), hits = hit_hits,
                     hit_rate = hit_rate, mcse = sqrt(hit_rate * (1 - hit_rate) / n_datasets),
                     row.names = NULL)

  rec_min <- vapply(per_rep, function(r) r$rec_min, numeric(1))
  rec_mean <- vapply(per_rep, function(r) r$rec_mean, numeric(1))
  converged <- vapply(per_rep, function(r) r$converged, logical(1))
  heywood <- vapply(per_rep, function(r) r$heywood, logical(1))
  fit_ok <- vapply(per_rep, function(r) r$fit_ok, logical(1))
  fit_error <- vapply(per_rep, function(r) r$fit_error, character(1))

  # A replicate without a congruence (failed fit) counts as not recovered.
  recovery <- if (have_lambda) {
    list(min_rate = sum(rec_min >= recovery_threshold, na.rm = TRUE) / n_datasets,
         mean_rate = sum(rec_mean >= recovery_threshold, na.rm = TRUE) / n_datasets,
         threshold = recovery_threshold, n_valid = sum(!is.na(rec_min)))
  }

  n_fit_ok <- sum(fit_ok)
  convergence <- list(
    n_datasets = as.integer(n_datasets), n_fit_ok = n_fit_ok,
    n_converged = sum(converged), n_heywood = sum(heywood),
    fit_rate = mean(fit_ok),
    convergence_rate = if (n_fit_ok > 0) sum(converged) / n_fit_ok else NA_real_,
    heywood_rate = if (n_fit_ok > 0) sum(heywood) / n_fit_ok else NA_real_)

  n_fit_failed <- as.integer(n_datasets) - n_fit_ok
  if (n_fit_failed > 0L) {
    first_failed <- which(!fit_ok)[1L]
    first_cause <- strsplit(fit_error[first_failed], "\n", fixed = TRUE)[[1L]][1L]
    cli::cli_warn(
      c("The {k_true}-factor fit failed on {n_fit_failed} of the {n_datasets} replicate{?s}.",
        "x" = if (have_lambda) {
          "{cli::qty(n_fit_failed)}{?This replicate counts/These replicates count} as not recovered."
        },
        "i" = "Replicate {first_failed} failed with: {first_cause}"),
      class = "efa_power_fit_failed")
  }

  settings <- list(mode = "simulation", N = as.integer(N),
                   n_datasets = as.integer(n_datasets), n_vars = p_true, n_factors = k_true,
                   criteria = criteria, estimator = estimator, rotation = rotation,
                   recovery_threshold = recovery_threshold, model_error = model_error,
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

# Analyse one simulated sample: the suggestion of each retention criterion, one k_true-factor
# fit (convergence, Heywood case), and the matched-factor Tucker congruences with the
# population loadings. A criterion or fit that fails gives no value.
.efa_power_analyze_one <- function(dat, Lambda, k_true, criteria, ctl, estimator,
                                   rotation) {

  Rmat <- stats::cor(dat)
  n_hat <- unlist(lapply(criteria, function(id) {
    entry <- .retention_registry[[id]]
    xarg <- if (isTRUE(entry$needs_raw)) dat else Rmat
    out <- suppressMessages(suppressWarnings(
      tryCatch(entry$fun(xarg, ctl), error = function(e) NULL)))
    if (is.null(out)) return(NULL)
    .retention_key(id, out$n_factors)
  }))

  # The flags come from the returned object, so the fit's own conditions are muffled.
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

  # Reorder and reflect the fitted factors to match the population; a degenerate factor has no
  # congruence.
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
#' `print()` shows an [efa_power()] result as a short report, and `format()` gives the same
#' report as a character vector (`print()` is `cat(format(x), sep = "\n")`).
#'
#' For an RMSEA-mode object, the report shows the test, the hypotheses, the significance level,
#' the degrees of freedom, and the power at the sample size or the required sample size.
#'
#' For a simulation-mode object, the report shows the population and design, the hit-rate per
#' criterion, the structure recovery, the convergence, and an upper bound for the Monte Carlo
#' standard error of the rates.
#'
#' The lines follow the active console theme, so they are plain text when colours are off,
#' for example in a file, or after [cli::ansi_strip()].
#'
#' @param x An object of class `efa_power` (output from [efa_power()]).
#' @param digits Integer. The number of decimal places. Default is 3.
#' @param ... Not used; for consistency with the generic.
#'
#' @returns `print()` returns `x` invisibly. `format()` returns a character vector with the
#'   report lines.
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

  if (identical(x$settings$mode, "simulation")) {
    return(.format_efa_power_simulation(x, digits = digits))
  }

  s <- x$settings
  num <- function(v) .efa_num(v, digits = digits, pad = FALSE)
  test_lbl <- if (s$type == "close") "Test of close fit" else "Test of not-close fit"
  # The null bounds the RMSEA from above for close fit, from below for not-close fit.
  cmp <- if (s$type == "close") "\u2264" else "\u2265"
  e0 <- num(s$eps0)
  e1 <- num(s$eps1)
  a <- num(s$alpha)
  pw <- num(x$power)
  # With groups, `N` is the total; the group size is shown as is (a fraction for a given N
  # that is not a multiple of `group`).
  n_note <- if (s$group > 1) {
    paste0(" (total; ", format(round(x$N_per_group, 1), scientific = FALSE), " per group)")
  } else {
    ""
  }

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

    if (x$solve_for == "power") {
      cli::cli_text("{.strong Power = {pw}} at N = {x$N}{n_note}.")
    } else {
      tgt <- num(s$power)
      cli::cli_text("{.strong Required N = {x$N}}{n_note} for a power of {tgt} (achieved {pw}).")
    }
  })
}

# Report lines for a simulation-mode efa_power object. Numbers are formatted to strings
# first, so that cli's `{...}` interpolation does not read `.efa_num(...)` as a style.
.format_efa_power_simulation <- function(x, digits = 3) {

  s <- x$settings
  conv <- x$convergence
  n <- conv$n_datasets
  num <- function(v) .efa_num(v, digits = digits, pad = FALSE)

  n_missing <- n - x$hits$n_valid
  hit_lines <- paste0(x$hits$criterion, ": ", num(x$hits$hit_rate),
                      ifelse(n_missing > 0, paste0(" (no suggestion in ", n_missing, ")"), ""))

  conv_ok <- conv$n_fit_ok == n && conv$n_converged == n && conv$n_heywood == 0
  conv_lines <- c(
    paste0("fits completed: ", num(conv$fit_rate), " (", conv$n_fit_ok, "/", n, ")"),
    paste0("converged (of completed): ", num(conv$convergence_rate),
           " (", conv$n_converged, "/", conv$n_fit_ok, ")"),
    paste0("Heywood cases (of completed): ", num(conv$heywood_rate),
           " (", conv$n_heywood, "/", conv$n_fit_ok, ")"))
  # 0.5 / sqrt(n) is an upper bound for the Monte Carlo standard error of a rate over n
  # replicates; the convergence and Heywood rates have the smaller denominator n_fit_ok.
  mcse_max <- num(0.5 / sqrt(if (conv$n_fit_ok > 0) conv$n_fit_ok else n))

  cli::cli_format_method({
    cli::cli_text("")
    cli::cli_rule(left = "{.strong EFA power simulation}")
    cli::cli_text("")

    cli::cli_text(
      "{s$n_vars} variable{?s} \u00b7 {x$k_true} factor{?s} \u00b7 N = {s$N} \u00b7 {n} dataset{?s}")
    cli::cli_text("Estimation: {s$estimator} \u00b7 rotation: {s$rotation}")
    # Model error applies only to a factor-model population, which is also the one with recovery.
    if (!is.null(x$model_error)) {
      me <- x$model_error
      me_rmsea <- num(me$rmsea)
      me_cfi <- num(me$cfi)
      cli::cli_text("Model error ({me$method}): RMSEA = {me_rmsea} \u00b7 CFI = {me_cfi}")
    } else if (!is.null(x$recovery)) {
      cli::cli_text("Model error: none (an exact population gives optimistic rates).")
    }
    cli::cli_text("")

    cli::cli_text("{.strong Retention hit-rate} P(k-hat = {x$k_true})")
    cli::cli_ul(hit_lines)
    cli::cli_text("")

    if (!is.null(x$recovery)) {
      r <- x$recovery
      thr <- num(r$threshold)
      min_rate <- num(r$min_rate)
      mean_rate <- num(r$mean_rate)
      # The median congruence tells how far off the structure is, also when no replicate passes.
      med_min <- num(stats::median(x$replicates$rec_min, na.rm = TRUE))
      cli::cli_text("{.strong Structure recovery} (Tucker congruence \u2265 {thr})")
      cli::cli_ul(c(
        paste0("recovery rate (min congruence): ", min_rate),
        paste0("recovery rate (mean congruence): ", mean_rate),
        paste0("median min congruence: ", med_min)))
      cli::cli_text("")
    }

    if (conv_ok) {
      cli::cli_text("{.strong Convergence}: all {n} fit{?s} converged, without Heywood cases.")
    } else {
      cli::cli_text("{.strong Convergence}")
      cli::cli_ul(conv_lines)
    }
    cli::cli_text("")
    cli::cli_text("Monte Carlo standard error of each rate: {mcse_max} or less.")
  })
}
