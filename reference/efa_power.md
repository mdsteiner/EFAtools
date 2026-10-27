# Power analysis for exploratory factor analysis

Analyses power for exploratory factor analysis in one of two modes,
chosen with `mode`.

`mode = "rmsea"` (the default) gives the analytic power of the root mean
square error of approximation (RMSEA) tests of close and not-close fit
(MacCallum, Browne, & Sugawara, 1996). Give a sample size to get the
power of the test, or give a target power to get the necessary sample
size.

`mode = "simulation"` runs a Monte-Carlo study: it draws `n_datasets`
samples from a known population (with
[`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md)),
analyses each sample, and reports how frequently the analysis finds the
population structure. A call with `Lambda` or `R` and no `mode` uses
this mode.

## Usage

``` r
efa_power(
  mode = c("rmsea", "simulation"),
  type = c("close", "notclose"),
  eps0 = NULL,
  eps1 = NULL,
  N = NULL,
  n_vars = NULL,
  n_factors = NULL,
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
  rotation = "oblimin",
  recovery_threshold = 0.95,
  model_error = c("TKL", "CB", "WB", "none"),
  target_rmsea = NULL,
  target_cfi = NULL,
  seed = NULL,
  p = lifecycle::deprecated(),
  k = lifecycle::deprecated()
)
```

## Arguments

- mode:

  character. `"rmsea"` (the default) or `"simulation"`. If `mode` is not
  given, `Lambda` or `R` selects `"simulation"`. `type`, `eps0`, `eps1`,
  `df`, `alpha`, `power`, and `group` apply only to RMSEA mode. The
  arguments marked *Simulation mode* apply only to simulation mode. An
  argument of the other mode is an error.

- type:

  character. The RMSEA test: `"close"` or `"notclose"`. See *Details*.

- eps0:

  numeric. The null-hypothesis RMSEA. Default is `0.05`.

- eps1:

  numeric. The alternative-hypothesis RMSEA, at which power is
  calculated. Default is `0.08` for `type = "close"` and `0.01` for
  `type = "notclose"`.

- N:

  numeric. The sample size. In RMSEA mode, the total of all groups: give
  `N` to get the power, or leave it `NULL` to get the required `N`. In
  simulation mode (required), the size of each sample.

- n_vars:

  numeric. The number of variables. In RMSEA mode, it gives `df`
  together with `n_factors`. In simulation mode, it comes from the
  population; if given, it must agree.

- n_factors:

  numeric. The number of factors. In RMSEA mode, it gives `df` together
  with `n_vars`. In simulation mode, the true number of factors:
  necessary with `R`; with `Lambda`, if given, it must agree with
  `ncol(Lambda)`.

- df:

  numeric. The model degrees of freedom. Give `df`, or `n_vars` and
  `n_factors`, which give
  `df = group * ((n_vars - n_factors)^2 - (n_vars + n_factors)) / 2` (a
  separate EFA in each group).

- alpha:

  numeric. The significance level. Default is `0.05`.

- power:

  numeric. The target power for the required `N`. Default is `NULL`,
  which uses `0.80` if `N` is also `NULL`. Give only one of `N` and
  `power`.

- group:

  numeric. The number of groups of equal size. Default is `1`. See
  *Details*.

- Lambda:

  matrix. Simulation mode. The population loadings (variables by
  factors) of a factor-model population. Structure recovery is available
  only with `Lambda`.

- Phi:

  matrix. Simulation mode. The population factor correlations. Only used
  with `Lambda`. Default is orthogonal factors.

- Psi:

  numeric or matrix. Simulation mode. The population unique variances.
  Only used with `Lambda`. See
  [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md).

- R:

  matrix. Simulation mode. A population correlation matrix, instead of
  `Lambda`. It needs `n_factors`.

- n_datasets:

  numeric. Simulation mode. The number of samples. Default is `500`.

- criteria:

  character. Simulation mode. The retention criteria for the hit-rate:
  any of `"CD"`, `"EKC"`, `"HULL"`, `"KGC"`, `"MAP"`, `"NEST"`,
  `"PARALLEL"`, and `"SMT"` (see
  [`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)).
  Default is `c("EKC", "MAP")`. The criteria use the
  [`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
  defaults and `estimator`. `"CD"` and `"HULL"` are much slower.

- estimator:

  character. Simulation mode. The estimator of the recovery fit and the
  retention criteria: `"PAF"` (the default), `"ML"`, or `"ULS"`.

- rotation:

  character. Simulation mode. The rotation of the recovery fit (see
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)).
  Default is `"oblimin"`, as in
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md).
  Recovery only reorders and reflects the factors. Thus, `"none"` with
  more than one factor understates recovery.

- recovery_threshold:

  numeric. Simulation mode. The Tucker congruence that a replicate must
  reach to be recovered. Default is `0.95`, at which Lorenzo-Seva and
  ten Berge (2006) consider two factors equal.

- model_error:

  character. Simulation mode. The
  [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md)
  method for model error: `"TKL"` (the default; it adds minor common
  factors), `"CB"`, `"WB"`, or `"none"`. It applies only with
  `target_rmsea` or `target_cfi`, which need `Lambda`. Note that
  [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md)
  has the default `"CB"`.

- target_rmsea:

  numeric. Simulation mode. The population RMSEA. Default is `NULL` (no
  model error). See
  [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md).

- target_cfi:

  numeric. Simulation mode. The population CFI (only with `"TKL"`).
  Default is `NULL`. See
  [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md).

- seed:

  numeric. Simulation mode. An optional seed for a reproducible result.
  Default is `NULL`.

- p, k:

  **\[deprecated\]** Use `n_vars` and `n_factors`.

## Value

An object of class `efa_power`. For `mode = "rmsea"`, a list with:

- power:

  The power at `N` (for a required `N`, at least the target power).

- N:

  The total sample size: the given `N`, or the required `N` (a multiple
  of `group`).

- N_per_group:

  The group size `N / group`.

- crit:

  The critical chi-square value.

- ncp:

  The noncentrality parameters under the null (`H0`) and the alternative
  (`H1`).

- solve_for:

  `"power"` or `"N"`: the calculated quantity.

- settings:

  The inputs: `mode`, `type`, `eps0`, `eps1`, `df`, `n_vars`,
  `n_factors`, `alpha`, `group`, and the target `power` (`NULL` if
  `solve_for` is `"power"`).

For `mode = "simulation"`, a list with:

- hit_rate:

  The hit-rate per criterion, and per variant if a criterion has more
  than one.

- hits:

  A data frame with one row per criterion or variant: `criterion`, the
  number of replicates with a suggestion (`n_valid`), the number of hits
  (`hits`), `hit_rate` (`hits / n_datasets`), and its Monte Carlo
  standard error (`mcse`).

- recovery:

  For `Lambda`, a list with the recovery rates (`min_rate`,
  `mean_rate`), the `threshold`, and the number of replicates with a
  congruence (`n_valid`); `NULL` for `R`. Non-converged and Heywood
  solutions are included.

- convergence:

  A list with `n_datasets`, the number of completed fits (`n_fit_ok`),
  the number of those that converged (`n_converged`) or have a Heywood
  case (`n_heywood`), and the rates `fit_rate` (of all datasets),
  `convergence_rate`, and `heywood_rate` (of the completed fits).

- replicates:

  The per-replicate values: the suggested numbers of factors (`n_hat`),
  the matched congruences (`rec_min`, `rec_mean`), the flags
  `converged`, `heywood`, and `fit_ok`, and `fit_error`, the error
  message of a failed fit (`NA` otherwise).

- k_true:

  The true number of factors.

- model_error:

  The
  [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md)
  model-error record, or `NULL`.

- settings:

  The simulation inputs.

## RMSEA mode

Power increases with the sample size, the degrees of freedom, and the
difference between the null and the alternative RMSEA (MacCallum et al.,
1996). `type` selects the test:

- `"close"`:

  The null hypothesis is close fit (RMSEA \\\le\\ `eps0`, conventionally
  0.05). Power is the probability to reject it if the RMSEA is `eps1`
  (conventionally 0.08). `eps0` must be smaller than `eps1`.

- `"notclose"`:

  The null hypothesis is not-close fit (RMSEA \\\ge\\ `eps0`). Power is
  the probability to reject it if the RMSEA is `eps1` (conventionally
  0.01). `eps0` must be larger than `eps1`.

The required sample size is the smallest `N` that reaches `power`. The
result applies to the RMSEA test only: it does not make sure that the
estimates are precise. With many degrees of freedom, the required `N`
can be smaller than the number of variables. The sample correlation
matrix is then singular (MacCallum et al., 1996, p. 144). If `n_vars` is
known, such a sample size gives a warning.

With `group > 1`, `N` is the total of `group` groups of equal size, and
`df` is that of the whole multiple-group model. The noncentrality
parameter is \\(N - G)\\ df\\ \varepsilon^2 / G\\ for \\G\\ groups
(Steiger, 1998). Thus, at a given `df`, power depends only on the group
size `N / group`, and a required `N` is `group` times the smallest group
size that reaches `power`.

## Simulation mode

[`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md)
draws `n_datasets` samples of size `N` from the population. The true
number of factors, `k_true`, is `ncol(Lambda)`, or `n_factors` for a
correlation matrix `R`. Without `target_rmsea` or `target_cfi`, a
factor-model population fits the model exactly, which gives optimistic
results. A misfit target makes the population more realistic (MacCallum,
2003).

The hit-rate and the recovery rates are proportions of all `n_datasets`
replicates:

- **Hit-rate**:

  A criterion in `criteria` suggests `k_true` factors. A replicate where
  the criterion fails or gives no suggestion counts as a miss.

- **Structure recovery** (factor-model populations only):

  The `k_true`-factor model is fitted with
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md),
  and its factors are matched to the population loadings. A replicate is
  recovered if the smallest (`min`) or the mean (`mean`) Tucker
  congruence of the matched factors reaches `recovery_threshold`
  (Lorenzo-Seva & ten Berge, 2006). A failed fit counts as not
  recovered.

The Monte Carlo standard error of a rate \\p\\ is \\\sqrt{p(1 - p) /
n}\\ for \\n\\ replicates (Morris, White, & Crowther, 2019). The
convergence and Heywood rates are proportions of the fits that
completed. If fits fail, a warning gives their number and the first
cause.

The replicates are analysed in parallel with future.apply; set a
parallel plan with
[`future::plan()`](https://future.futureverse.org/reference/plan.html).
With a fixed `seed`, the result does not depend on the number of
workers, and your random-number stream does not change.

## References

Lorenzo-Seva, U., & ten Berge, J. M. F. (2006). Tucker's congruence
coefficient as a meaningful index of factor similarity. *Methodology,
2*(2), 57-64.
[doi:10.1027/1614-2241.2.2.57](https://doi.org/10.1027/1614-2241.2.2.57)

MacCallum, R. C. (2003). 2001 Presidential Address: Working with
imperfect models. *Multivariate Behavioral Research, 38*(1), 113-139.
[doi:10.1207/S15327906MBR3801_5](https://doi.org/10.1207/S15327906MBR3801_5)

MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power
analysis and determination of sample size for covariance structure
modeling. *Psychological Methods, 1*(2), 130-149.
[doi:10.1037/1082-989X.1.2.130](https://doi.org/10.1037/1082-989X.1.2.130)

Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
studies to evaluate statistical methods. *Statistics in Medicine,
38*(11), 2074-2102.
[doi:10.1002/sim.8086](https://doi.org/10.1002/sim.8086)

Steiger, J. H. (1998). A note on multiple sample extensions of the RMSEA
fit index. *Structural Equation Modeling, 5*(4), 411-419.
[doi:10.1080/10705519809540115](https://doi.org/10.1080/10705519809540115)

## See also

[`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md)
draws the samples of simulation mode.
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
implements the retention criteria.

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

# df from the number of variables and factors
efa_power(n_vars = 20, n_factors = 3, N = 200)
#> 
#> ── RMSEA power analysis ────────────────────────────────────────────────────────
#> 
#> Test of close fit: H0 RMSEA ≤ .050 vs. H1 RMSEA = .080.
#> alpha = .050 · df = 133
#> 
#> Power = .986 at N = 200.

# Required sample size for 80% power
efa_power(df = 100, power = 0.80)
#> 
#> ── RMSEA power analysis ────────────────────────────────────────────────────────
#> 
#> Test of close fit: H0 RMSEA ≤ .050 vs. H1 RMSEA = .080.
#> alpha = .050 · df = 100
#> 
#> Required N = 132 for a power of .800 (achieved .802).

# Test of not-close fit
efa_power(df = 100, N = 200, type = "notclose")
#> 
#> ── RMSEA power analysis ────────────────────────────────────────────────────────
#> 
#> Test of not-close fit: H0 RMSEA ≥ .050 vs. H1 RMSEA = .010.
#> alpha = .050 · df = 100
#> 
#> Power = .870 at N = 200.

# \donttest{
# Simulation mode: hit-rate and structure recovery for a known three-factor
# population at N = 300
efa_power("simulation", Lambda = population_models$loadings$baseline,
          Phi = population_models$phis_3$moderate, N = 300,
          n_datasets = 50, criteria = c("EKC", "MAP"), seed = 42)
#> 
#> ── EFA power simulation ────────────────────────────────────────────────────────
#> 
#> 18 variables · 3 factors · N = 300 · 50 datasets
#> Estimation: PAF · rotation: oblimin
#> Model error: none (an exact population gives optimistic rates).
#> 
#> Retention hit-rate P(k-hat = 3)
#> • EKC_BvA2017: .980
#> • MAP_TR2: 1.000
#> • MAP_TR4: 1.000
#> 
#> Structure recovery (Tucker congruence ≥ .950)
#> • recovery rate (min congruence): 1.000
#> • recovery rate (mean congruence): 1.000
#> • median min congruence: .987
#> 
#> Convergence: all 50 fits converged, without Heywood cases.
#> 
#> Monte Carlo standard error of each rate: .071 or less.
# }
```
