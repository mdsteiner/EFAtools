# Power analysis for exploratory factor analysis

Power analysis for exploratory factor analysis, in two modes selected by
`mode`.

`mode = "rmsea"` (the default) gives the analytic statistical power of
the root mean square error of approximation (RMSEA) tests of close and
not-close fit (MacCallum, Browne, & Sugawara, 1996): given the model
degrees of freedom and a sample size it returns the power to reject the
null hypothesis, and given a target power instead of a sample size it
returns the sample size required to reach it.

`mode = "simulation"` runs a Monte-Carlo study of a known population: it
draws `n_datasets` samples of size `N` (via
[`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md))
and, over the replicates, reports how often each factor-retention
criterion recovers the true number of factors (the *hit-rate*), how
often the fitted loadings recover the population structure (*structure
recovery*, by Tucker congruence), and the convergence and Heywood-case
rate of the fit.

The model dimensions follow the notation of MacCallum, Browne and
Sugawara (1996): `p` is the number of observed variables (the `n_vars`
of the rest of the package) and `k` the number of factors (`n_factors`
elsewhere).

## Usage

``` r
efa_power(
  mode = c("rmsea", "simulation"),
  type = c("close", "notclose"),
  eps0 = NULL,
  eps1 = NULL,
  N = NULL,
  p = NULL,
  k = NULL,
  df = NULL,
  alpha = 0.05,
  power = NULL,
  group = 1,
  Lambda = NULL,
  Phi = NULL,
  Psi = NULL,
  R = NULL,
  n_datasets = 500,
  criteria = c("EKC", "MAP"),
  estimator = "PAF",
  rotation = NULL,
  recovery_threshold = 0.95,
  model_error = c("TKL", "CB", "WB", "none"),
  target_rmsea = NULL,
  target_cfi = NULL,
  seed = NULL
)
```

## Arguments

- mode:

  character. The kind of power analysis: `"rmsea"` (the default;
  analytic RMSEA power) or `"simulation"` (Monte-Carlo hit-rate and
  structure recovery). `type`, `eps0`, `eps1`, `df`, `alpha`, `power`,
  and `group` apply to RMSEA mode only; the arguments marked *Simulation
  mode* below apply to the other; `N`, `p`, and `k` are used in both.

- type:

  character. The RMSEA test: `"close"` (test of close fit) or
  `"notclose"` (test of not-close fit). See *Details*.

- eps0:

  numeric. The null-hypothesis RMSEA. Default is `0.05`.

- eps1:

  numeric. The alternative-hypothesis RMSEA (the true RMSEA power is
  evaluated at). Default is `0.08` for `type = "close"` and `0.01` for
  `type = "notclose"`.

- N:

  numeric. In `"rmsea"` mode, the total sample size across groups (the
  plain sample size when `group` is `1`): give `N` to compute power, or
  leave it `NULL` to solve for the required `N` at a target `power`. In
  `"simulation"` mode `N` is required and is the size of each drawn
  sample, with no sample size solved for and no `group` division.

- p:

  numeric. The number of observed variables. In `"rmsea"` mode, used
  with `k` to derive `df` when `df` is not given directly. In
  `"simulation"` mode it is read off the population, so leave it unset
  (or matching `nrow(Lambda)` / `nrow(R)`).

- k:

  numeric. The number of factors. In `"rmsea"` mode, used with `p` to
  derive `df` when `df` is not given directly. In `"simulation"` mode it
  is the true number of factors: it is required with an `R` population
  and must be left unset (or match `ncol(Lambda)`) with a factor-model
  population.

- df:

  numeric. The model degrees of freedom. Either supply `df` directly or
  supply both `p` and `k`, from which `df = ((p - k)^2 - (p + k)) / 2`.
  Must be positive.

- alpha:

  numeric. The significance level. Default is `0.05`.

- power:

  numeric. The target power. Give `power` (or leave both `power` and `N`
  `NULL`, defaulting to `0.80`) to solve for the required `N`; leave it
  `NULL` while giving `N` to compute power. Exactly one of `N` and
  `power` is solved for.

- group:

  numeric. The number of groups. Default is `1`. `N` is the total across
  all `group` groups, not the size of each one, and a solved `N` is a
  multiple of `group`. See *Details*.

- Lambda:

  matrix. Simulation mode. A `p` by `k_true` population loading matrix.
  Supply this (optionally with `Phi`/`Psi`) to build a factor-model
  population; structure recovery is available only with this form.
  Passed to
  [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md).

- Phi:

  matrix. Simulation mode. The `k_true` by `k_true` population factor
  intercorrelations. Only used with `Lambda`; defaults to orthogonal
  factors. When `rotation` is unset, an oblique `Phi` selects a
  `"promax"` recovery fit and an orthogonal one a `"varimax"` fit.

- Psi:

  numeric or matrix. Simulation mode. The population unique variances (a
  length-`p` vector or a `p` by `p` matrix). Only used with `Lambda`.
  Passed to
  [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md).

- R:

  matrix. Simulation mode. A `p` by `p` population correlation matrix to
  draw from directly, instead of a factor model. Structure recovery is
  not available for this form (there are no population loadings to
  recover), and `k` is required.

- n_datasets:

  numeric. Simulation mode. The number of samples to draw and analyse.
  Default is `500`.

- criteria:

  character. Simulation mode. The factor-retention criteria to evaluate
  the hit-rate for, any of `"CD"`, `"EKC"`, `"HULL"`, `"KGC"`, `"MAP"`,
  `"NEST"`, `"PARALLEL"`, and `"SMT"` (see
  [`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)).
  Default is `c("EKC", "MAP")`. Criteria that simulate internally
  (`"CD"`, `"HULL"`, `"NEST"`, `"PARALLEL"`) make each run substantially
  slower.

- estimator:

  character. Simulation mode. The estimator (`"PAF"`, `"ML"`, or
  `"ULS"`) used for the recovery fit and the retention criteria. Default
  is `"PAF"`.

- rotation:

  character. Simulation mode. The rotation for the recovery fit, passed
  to
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md).
  Default is `NULL`, which matches the population: `"varimax"` for
  orthogonal factors and `"promax"` for oblique ones (a single factor is
  left unrotated). Recovery aligns the fitted loadings to the population
  pattern by permutation and sign only, so a rotation that does not seek
  that structure – for example `"none"` with more than one factor –
  understates recovery; keep the default (or another structure-seeking
  rotation) for a meaningful recovery rate.

- recovery_threshold:

  numeric. Simulation mode. The matched-factor Tucker congruence a
  replicate must reach to count as recovered. Default is `0.95` (the
  Lorenzo-Seva & ten Berge, 2006, "equal" band).

- model_error:

  character. Simulation mode. The
  [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md)
  method used to perturb the population with model error: `"TKL"`
  (Tucker-Koopman-Linn, the default), `"CB"` (Cudeck-Browne), `"WB"`
  (Wu-Browne), or `"none"` for an exact population. It takes effect only
  when a target is supplied (`target_rmsea` and/or `target_cfi`);
  without one the population is exact whatever the method, and only a
  factor-model population can be perturbed. `"TKL"` adds minor common
  factors, so it degrades both the retention hit-rate and structure
  recovery realistically; `"CB"` and `"WB"` target the RMSEA only, and
  `"CB"` keeps the fitted loadings the exact minimizer so structure
  recovery stays near-perfect. That is why the default here is `"TKL"`
  while
  [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md)
  defaults to `"CB"`: the same target passed to the two functions
  perturbs the population differently unless `model_error` is set
  explicitly, and the difference is material – a `"CB"` population can
  leave recovery near 1 where the same target under `"TKL"` does not.

- target_rmsea:

  numeric. Simulation mode. The population RMSEA the model should have
  relative to the perturbed population, activating model error. Default
  is `NULL`. Passed to
  [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md).

- target_cfi:

  numeric. Simulation mode. The population CFI target (only with
  `model_error = "TKL"`). Default is `NULL`. Passed to
  [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md).

- seed:

  numeric. Simulation mode. Optional seed making the draws and analysis
  reproducible and worker-count independent; the caller's random-number
  stream is restored afterwards. Default is `NULL`.

## Value

An object of class `efa_power`. For `mode = "rmsea"`, a list containing:

- power:

  The power of the test at `N` (the achieved power, which for a solved
  sample size is at least the target).

- N:

  The total sample size across groups: the supplied `N`, or the solved
  required sample size (a multiple of `group`).

- N_per_group:

  The per-group sample size `N / group`, equal to `N` when `group` is
  `1`. A whole number for a solved `N`; for a supplied `N` that is not a
  multiple of `group` it is the fraction that the noncentrality uses.

- crit:

  The critical chi-square value the fit statistic is compared against.

- ncp:

  The noncentrality parameters under the null (`H0`, from `eps0`) and
  the alternative (`H1`, from `eps1`).

- solve_for:

  `"power"` or `"N"`, recording which quantity was solved for.

- settings:

  A list of the inputs: `mode`, `type`, `eps0`, `eps1`, `df`, `p`, `k`,
  `alpha`, `group`, and the target `power` (the value solved to when
  `solve_for` is `"N"`, otherwise `NULL`).

For `mode = "simulation"`, a list containing:

- hit_rate:

  A named numeric vector of the retention hit-rate per criterion (and,
  where a criterion has several variants, per variant); `NA` for a
  criterion that returned no suggestion on any replicate.

- hits:

  A data frame with, per criterion, the number of replicates it returned
  a definite suggestion on (`n_valid`), the number of those that matched
  `k_true` (`hits`), and the `hit_rate` (`hits / n_valid`).

- recovery:

  For a factor-model population, a list with the structure-recovery
  rates (`min_rate`, `mean_rate`), the `threshold`, and the number of
  usable fits (`n_valid`); `NULL` for an `R` population. Rates are over
  every replicate whose fit returned loadings, including non-converged
  or Heywood solutions (their rates are reported separately in
  `convergence`).

- convergence:

  A list with the number of datasets, the number of fits that completed
  (`n_fit_ok`), and rates: `fit_rate` (fits completed, over all
  datasets) and `convergence_rate` / `heywood_rate` (converged /
  Heywood, over the completed fits).

- replicates:

  The raw per-replicate values: the suggested factor counts (`n_hat`),
  the matched congruences (`rec_min`, `rec_mean`), the `converged`,
  `heywood`, and `fit_ok` flags, and `fit_error`, the message of the fit
  that did not complete (`NA` where it did).

- k_true:

  The true number of factors.

- model_error:

  The
  [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md)
  model-error record, or `NULL`.

- settings:

  A list of the simulation inputs.

## RMSEA mode

The RMSEA fit statistic is referred to a noncentral chi-square
distribution. For a model with `df` degrees of freedom and population
RMSEA \\\varepsilon\\, the noncentrality parameter is \$\$\lambda = (N -
1)\\ df\\ \varepsilon^2 / G,\$\$ where \\N\\ is the sample size and
\\G\\ the number of groups (`group`).

Two tests are supported, selected by `type` (never by the ordering of
`eps0` and `eps1`):

- `"close"`:

  The test of *close fit* (MacCallum et al., 1996): the null hypothesis
  is that the fit is close (RMSEA \\\le\\ `eps0`, conventionally 0.05),
  and power is computed against a worse alternative (`eps1`,
  conventionally 0.08, so `eps0 < eps1`), in the upper tail.

- `"notclose"`:

  The test of *not-close fit*: the null hypothesis is that the fit is
  not close (RMSEA \\\ge\\ `eps0`), and power is computed against a
  better alternative (`eps1`, conventionally 0.01, so `eps0 > eps1`), in
  the lower tail.

When `eps0` and `eps1` are ordered the wrong way round for the chosen
`type` a message is emitted but the requested test is still computed.
Equal `eps0` and `eps1` leave nothing to detect and are an error.

Power increases monotonically with `N`, so the required sample size (the
smallest `N` reaching `power`) is found by bisection. `N` is the
**total** sample size across groups: with `group > 1` the noncentrality
above carries the `1 / group` factor, so a fixed total spread over more
groups buys less power. The corresponding per-group sample size
`N / group` is returned as `N_per_group`.

The `1 / group` factor makes all `group` groups the same size, so a
required total is rounded up to the next multiple of `group`. A solved
`N_per_group` is thus a whole number of persons, and the reported
`power` is the power at a total that a study can collect. With
`group = 2` and `df = 102`, for example, the required total is 260, or
130 per group. Bisection on the total alone gives 259, which asks for
129.5 persons in each group.

## Simulation mode

The population is passed to
[`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md),
which draws `n_datasets` samples of size `N`. Its true number of factors
`k_true` is `ncol(Lambda)` for a factor-model population, or `k` for a
bare `R`. By default the population fits the factor model exactly, which
overstates how well the criteria and the fit recover the structure; a
misfit target makes it a more realistic one (MacCallum, 2003).

Each replicate is analysed three ways. **Hit-rate**: every criterion in
`criteria` is run and its suggested number of factors compared with
`k_true`; the hit-rate is the proportion of agreements over the
replicates on which the criterion returned a definite factor count
(replicates where it errored or was undecided are excluded, not counted
as misses). **Structure recovery** (only for a factor-model population):
the `k_true`-factor model is fitted with
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md),
its rotated loadings are matched to the population loadings, and the
matched-factor Tucker congruences (Lorenzo-Seva & ten Berge, 2006) are
compared with `recovery_threshold`; a replicate counts as a success when
the smallest (`min`) or the average (`mean`) matched congruence reaches
it. **Convergence**: the same fit supplies the proportion of replicates
whose fit completed and, among those, the proportion that converged and
that produced a Heywood case. A replicate whose fit fails outright
contributes to none of the three; if any fit fails, a warning reports
how many did and the first replicate's cause.

Replicates are analysed in parallel with future.apply; a plan is
selected with
[`future::plan()`](https://future.futureverse.org/reference/plan.html).
Each replicate is bound to its own reproducible random-number stream, so
with a fixed `seed` the result is identical regardless of the number of
workers, and the caller's random-number state is left unchanged.

## References

MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power
analysis and determination of sample size for covariance structure
modeling. *Psychological Methods, 1*(2), 130-149.
[doi:10.1037/1082-989X.1.2.130](https://doi.org/10.1037/1082-989X.1.2.130)

MacCallum, R. C. (2003). 2001 Presidential Address: Working with
imperfect models. *Multivariate Behavioral Research, 38*(1), 113-139.
[doi:10.1207/S15327906MBR3801_5](https://doi.org/10.1207/S15327906MBR3801_5)

Lorenzo-Seva, U., & ten Berge, J. M. F. (2006). Tucker's congruence
coefficient as a meaningful index of factor similarity. *Methodology,
2*(2), 57-64.
[doi:10.1027/1614-2241.2.2.57](https://doi.org/10.1027/1614-2241.2.2.57)

## See also

[`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md),
which draws the replicate data sets in simulation mode, and
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
for the retention criteria whose hit rates that mode reports.

Other power analysis:
[`plot.efa_power()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_power.md),
[`print.efa_power()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_power.md)

## Examples

``` r
# Power of the test of close fit at N = 200 for a 100-df model
efa_power(df = 100, N = 200)
#> 
#> ── RMSEA power analysis ────────────────────────────────────────────────────────
#> 
#> Test of close fit: H0 RMSEA ≤ .050 vs. H1 RMSEA = .080.
#> alpha = .050 · df = 100
#> 
#> Power = .955 at N = 200.
#> Critical value χ²(100) = 183.967 · noncentrality H0 = 49.750, H1 = 127.360.

# Deriving df from the model dimensions instead of giving it directly
efa_power(p = 20, k = 3, N = 200)
#> 
#> ── RMSEA power analysis ────────────────────────────────────────────────────────
#> 
#> Test of close fit: H0 RMSEA ≤ .050 vs. H1 RMSEA = .080.
#> alpha = .050 · df = 133
#> 
#> Power = .986 at N = 200.
#> Critical value χ²(133) = 238.428 · noncentrality H0 = 66.168, H1 = 169.389.

# Required total sample size for 80% power
efa_power(df = 100, power = 0.80)
#> 
#> ── RMSEA power analysis ────────────────────────────────────────────────────────
#> 
#> Test of close fit: H0 RMSEA ≤ .050 vs. H1 RMSEA = .080.
#> alpha = .050 · df = 100
#> 
#> Required N = 132 for a power of .800 (achieved .802).
#> Critical value χ²(100) = 163.977 · noncentrality H0 = 32.750, H1 = 83.840.

# Test of not-close fit
efa_power(df = 100, N = 200, type = "notclose")
#> 
#> ── RMSEA power analysis ────────────────────────────────────────────────────────
#> 
#> Test of not-close fit: H0 RMSEA ≥ .050 vs. H1 RMSEA = .010.
#> alpha = .050 · df = 100
#> 
#> Power = .870 at N = 200.
#> Critical value χ²(100) = 118.375 · noncentrality H0 = 49.750, H1 = 1.990.

# \donttest{
# Simulation mode: retention hit-rate and structure recovery for a known
# three-factor population at N = 300
efa_power("simulation", Lambda = population_models$loadings$baseline,
          Phi = population_models$phis_3$moderate, N = 300,
          n_datasets = 50, criteria = c("EKC", "MAP"), seed = 42)
#> 
#> ── EFA power simulation ────────────────────────────────────────────────────────
#> 
#> 18 variables · 3 factors · N = 300 · 50 datasets
#> Estimation: PAF · rotation: promax
#> Model error: none. The population is exact, so the hit-rate and recovery are
#> optimistic; set `target_rmsea` for realism.
#> 
#> Retention hit-rate P(k-hat = 3)
#> • EKC_BvA2017: .980 (n = 50)
#> • MAP_TR2: 1.000 (n = 50)
#> • MAP_TR4: 1.000 (n = 50)
#> 
#> Structure recovery (Tucker congruence ≥ .950)
#> • recovery rate (min congruence): 1.000 (n = 50)
#> • recovery rate (mean congruence): 1.000 (n = 50)
#> • median min congruence: .987
#> 
#> Convergence
#> • fits completed: 1.000 (50/50)
#> • converged (of completed): 1.000
#> • Heywood cases (of completed): .000
# }
```
