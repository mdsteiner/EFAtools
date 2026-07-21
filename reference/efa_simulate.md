# Simulate data from a common-factor population model

Draws data from a population correlation matrix, given either directly
or built from a factor model. The population correlation is either
supplied in `R`, or assembled from a loading matrix `Lambda`, the factor
intercorrelations `Phi`, and the unique variances `Psi` as \\R =
Lambda\\ Phi\\ Lambda' + Psi\\, standardized to a correlation matrix. By
default (`marginals = "normal"`) cases are drawn with normal marginals
via a matrix square root of the population correlation (a Cholesky
factor, or a symmetric eigen square root for a positive-semidefinite but
singular population). With `marginals = "empirical"` the cases instead
reproduce the population correlation while carrying the empirical
marginal distributions of a supplied data set. With `marginals = "VM"`
or `"IG"` the cases reproduce the population correlation while carrying
non-normal marginals with a prescribed `skewness` and `kurtosis`.
Setting `categories` additionally discretizes the drawn data into
ordered categories, optionally so that the population polychoric
correlation of the categorized data equals the target correlation.
Setting `target_rmsea` or `target_cfi` perturbs the population with
model error, so the factor model fits it only approximately (a more
realistic simulation target).

## Usage

``` r
efa_simulate(
  N = NULL,
  Lambda = NULL,
  Phi = NULL,
  Psi = NULL,
  R = NULL,
  model_error = c("CB", "TKL", "WB", "none"),
  target_rmsea = NULL,
  target_cfi = NULL,
  marginals = c("normal", "empirical", "VM", "IG"),
  marginal_data = NULL,
  n_factors = NULL,
  skewness = NULL,
  kurtosis = NULL,
  force_pd = FALSE,
  categories = NULL,
  match = NULL,
  missing = c("none", "MCAR", "MAR", "MNAR"),
  missing_prop = NULL,
  missing_strength = NULL,
  missing_predictor = NULL,
  n_datasets = 1L,
  seed = NULL,
  return_pop = FALSE
)
```

## Arguments

- N:

  numeric. Number of cases (rows) to draw per dataset. Required unless
  `return_pop = TRUE`.

- Lambda:

  matrix. A `p` by `m` matrix of factor loadings. Supply this
  (optionally with `Phi` and `Psi`) instead of `R` to build the
  population from a factor model.

- Phi:

  matrix. The `m` by `m` factor intercorrelation matrix. Only used with
  `Lambda`. Default is `NULL`, in which case the factors are orthogonal
  (an identity matrix).

- Psi:

  numeric vector or matrix. The unique variances: either a length-`p`
  vector or a `p` by `p` matrix (added as the residual covariance). Only
  used with `Lambda`. Default is `NULL`, in which case the unique
  variances that standardize the population to a correlation matrix are
  used.

- R:

  matrix. A `p` by `p` population correlation matrix to draw from
  directly. Supply this instead of `Lambda`/`Phi`/`Psi`.

- model_error:

  character. The method used to perturb the population so the factor
  model fits it imperfectly ("model error"): one of `"CB"`
  (Cudeck-Browne, the default), `"TKL"` (Tucker-Koopman-Linn), `"WB"`
  (Wu-Browne), or `"none"`. Model error is only applied when a target is
  supplied in `target_rmsea` or `target_cfi`; without one the population
  is exact, whatever `model_error`. Only used with a factor-model
  population (`Lambda`).

- target_rmsea:

  numeric. The population RMSEA the factor model should have relative to
  the perturbed population, a single number strictly in `(0, 1)`.
  Supplying it activates model error. Simulating from an exact
  population overstates recovery, so a realistic value (around `0.05`)
  is recommended for simulation studies (MacCallum, 2003). Default is
  `NULL` (an exact population; do not pass `0`). Required for `"CB"` and
  `"WB"`; optional for `"TKL"`.

- target_cfi:

  numeric. Only used with `model_error = "TKL"`: the population CFI to
  target, a single number strictly in `(0, 1)`, on its own or together
  with `target_rmsea` (TKL then trades the two off). Default is `NULL`
  (do not pass `1`). `"CB"` and `"WB"` target the RMSEA only.

- marginals:

  character. The marginal distribution of the drawn data: one of
  `"normal"` (the default), which draws normal marginals; `"empirical"`,
  which reproduces the population correlation while preserving the
  empirical marginals supplied in `marginal_data`; or `"VM"`
  (Vale-Maurelli) and `"IG"` (independent generator), which draw
  non-normal marginals with the target `skewness` and `kurtosis`.

- marginal_data:

  matrix or data frame. Only used with `marginals = "empirical"`, where
  it is required: a data set with one numeric column per variable (`p`
  columns), each with at least two distinct values, whose per-column
  distributions are resampled as the marginals of the drawn data. Its
  correlations are ignored. Default is `NULL`.

- n_factors:

  numeric. Only used with `marginals = "empirical"`: the number of
  factors the rank-matching reproduction fits. Default is `NULL`, in
  which case it is the number of columns of `Lambda` when the population
  is built from a factor model; it must be given when the population is
  supplied via `R`.

- skewness:

  numeric. Only used with `marginals = "VM"` or `"IG"`: the target
  marginal skewness, as a single value applied to every variable or a
  length-`p` vector. Default is `NULL` (0, a symmetric marginal). At
  least one of `skewness` or `kurtosis` must be given for these
  marginals.

- kurtosis:

  numeric. Only used with `marginals = "VM"` or `"IG"`: the target
  marginal excess kurtosis (0 for a normal marginal), as a single value
  applied to every variable or a length-`p` vector. Default is `NULL`
  (0).

- force_pd:

  logical. Used with `marginals = "VM"` and with Cudeck-Browne model
  error (`model_error = "CB"`). If the Vale-Maurelli intermediate
  correlation matrix (or, for CB, the perturbed population at a large
  `target_rmsea`) is not positive definite, `FALSE` (the default)
  rejects it with an error, while `TRUE` projects it to the nearest
  correlation matrix (with a warning). Has no effect for the `"TKL"` or
  `"WB"` methods.

- categories:

  numeric or list. Requests ordinal output by discretizing each variable
  into ordered categories. Either a count of equally probable categories
  (a single value applied to every variable or a length-`p` vector), or
  a length-`p` list of numeric vectors giving the marginal category
  proportions per variable (each strictly positive and summing to 1).
  Default is `NULL`, which returns the continuous data.

- match:

  character. Only used with `categories`: an assertion about how the
  categorization relates to the population correlation. With a normal
  latent, cutting at the normal-scale thresholds already leaves the
  population polychoric correlation of the categorized data equal to the
  target correlation, so both values compute the same thresholds and
  produce identical data whenever both are legal. `"thresholds"` (the
  default) also cuts the `"VM"` and `"IG"` draws, whose ordinal Pearson
  and polychoric correlations then both depart from the population;
  `"polychoric"` states that the polychoric match is required and
  therefore rejects non-normal marginals. Not available with
  `marginals = "empirical"`. The value is matched case-insensitively.
  Default is `NULL` (`"thresholds"` when `categories` is set).

- missing:

  character. An optional missing-data mechanism to impose on the drawn
  data: one of `"none"` (the default, complete data), `"MCAR"` (missing
  completely at random), `"MAR"` (missing at random, depending on
  another variable), or `"MNAR"` (missing not at random, depending on
  the variable's own value). Introduced values become `NA`. Every
  variable is holed, so under `"MAR"` each variable's predictor is
  itself subject to missingness: the mechanism is MAR given the
  *complete* data and is not ignorable for an analyst who sees only the
  observed data (see Details).

- missing_prop:

  numeric. Only used when `missing` is not `"none"`, where it is
  required: the target marginal proportion of missing values per
  variable, a single number strictly between 0 and 1. This is the
  expected rate; the realized rate of a given draw varies around it.

- missing_strength:

  numeric. Only used with `missing = "MAR"` or `"MNAR"`: the slope of
  the logistic missingness model, setting how strongly the missing
  probability depends on the predictor. Default is `NULL` (1, a moderate
  dependence); `0` removes the dependence (equivalent to MCAR at the
  same rate) and large magnitudes make missingness nearly deterministic.

- missing_predictor:

  integer or character. Only used with `missing = "MAR"`: which variable
  drives each variable's missingness, as one column index or variable
  name per variable. Each must reference another variable, not itself
  (so a single shared predictor is not allowed – it would predict its
  own missingness). Default is `NULL`, in which case each variable's
  missingness is driven by the next variable cyclically (so variable
  order matters; supply this explicitly when the order is arbitrary).

- n_datasets:

  numeric. The number of datasets to draw. Default is 1. With more than
  one, a list of datasets is returned.

- seed:

  numeric. Optional seed for reproducible draws. When supplied, the
  caller's random-number stream is saved and restored, so the call
  leaves the global RNG state unchanged. Default is `NULL` (no seeding).

- return_pop:

  logical. If `TRUE`, return only the population correlation matrix and
  draw no data. Default is `FALSE`.

## Value

An object of class `efa_simulated`: a list with elements `data` (the
simulated data – an `N` by `p` numeric matrix, an integer matrix of
category codes when `categories` is set, or a length-`n_datasets` list
of these when `n_datasets > 1`; `NULL` when `return_pop = TRUE`),
`population` (the `p` by `p` population correlation matrix drawn from,
model-error-perturbed when requested; with `force_pd = TRUE` and
`marginals = "VM"` it stays the target matrix, from which the realized
correlations of the draw can drift), `model_error` (`NULL`, or a list of
the method and the target and achieved RMSEA/CFI when model error was
applied), and `settings`. Printing the object shows a compact summary.

## Details

Provide the population either as a ready correlation matrix in `R`, or
through the model components `Lambda`, `Phi`, and `Psi`; the two ways
are mutually exclusive. When the model components are used, `Phi`
defaults to the identity matrix (orthogonal factors) and `Psi` defaults
to the unique variances that make the population a correlation matrix
(\\1 - \mathrm{diag}(Lambda\\ Phi\\ Lambda')\\); the assembled
covariance is standardized with
[`stats::cov2cor()`](https://rdrr.io/r/stats/cor.html) so a
non-standardized `Psi` still yields a correlation matrix. With the
default `Psi`, a factor model whose implied communalities exceed 1 (a
Heywood case) leaves no unique variance and is rejected; a `Psi` you
supply is instead only required to give positive variances and a
positive-semidefinite population.

With `marginals = "empirical"`, the iterative rank-matching algorithm of
Ruscio and Kaczetow (2008) reproduces the population correlation while
each variable takes the empirical marginal distribution of the matching
column of `marginal_data` (resampled with replacement). Only the
marginals of `marginal_data` are used; its own correlations are ignored,
and the drawn columns follow the population's variables, not those of
`marginal_data`.

With `marginals = "VM"` (Vale-Maurelli, 1983) or `"IG"` (the
independent-generator method; Foldnes & Olsson, 2016), the cases
reproduce the population correlation while carrying non-normal marginals
with the target `skewness` and (excess) `kurtosis`. The Vale-Maurelli
family does not span every valid non-normal distribution (Foldnes &
Grønneberg, 2015); `"IG"` covers distributions `"VM"` cannot. Both
accept `skewness` and `kurtosis` as a single value (used for every
variable) or one value per variable, defaulting the unset one to 0. Not
every (`skewness`, `kurtosis`) pair is attainable: every distribution
needs excess kurtosis of at least `skewness^2 - 2`, and the method
covers a smaller region still, so an unreachable request is rejected.
For `marginals = "VM"`, the intermediate correlation matrix used for the
draw can itself be non-positive-definite; it is rejected unless
`force_pd = TRUE`, which projects it to the nearest correlation matrix
(via
[`psych::cor.smooth()`](https://rdrr.io/pkg/psych/man/cor.smooth.html))
with a warning.

With `categories`, the drawn data are discretized into ordered
categories (an integer code `1` to `K`). `categories` gives either the
number of equally probable categories (one count for every variable, or
one per variable) or, as a list of proportion vectors, the marginal
category proportions per variable. The cut points are the thresholds
that reproduce the requested proportions (Olsson, 1979): the
standard-normal quantiles for `marginals = "normal"`, and for
`marginals = "VM"` those quantiles mapped through the same Fleishman
cubic the draw uses, so the requested proportions are reproduced on the
non-normal scale too. Under `marginals = "IG"` the thresholds stay on
the standard-normal scale while the data do not, so the achieved
proportions depart from the request systematically rather than by
sampling noise; the departure grows with the non-normality, and only the
*number* of categories is guaranteed. The same holds for a `"VM"`
variable whose Fleishman cubic is not increasing over its own thresholds
and the tails beyond them, which keeps the normal-scale ones and is
reported with a warning; this arises when a turning point of the cubic
sits at or near an outer threshold – under a strongly platykurtic
marginal, or under substantial skewness or kurtosis combined with a
small outer-category proportion. Because categorization attenuates
product-moment correlations, the categorized data's Pearson correlation
is smaller in magnitude than the population correlation; under
non-normal marginals its polychoric correlation departs from the
population as well. `match` changes none of this – it asserts an intent
rather than selecting a computation, as described under that argument.
Ordinal output is not available with `marginals = "empirical"`. Empty
categories left by a draw are reported with a warning, as they
destabilize the polychoric correlation and the factor analysis.

With `missing`, missing values are introduced into the drawn data under
a chosen mechanism (Rubin, 1976), each variable holed at a target
expected rate `missing_prop`. `"MCAR"` draws an independent mask, so
missingness is unrelated to the data. `"MAR"` and `"MNAR"` set each
case's missing probability by a logistic model of a standardized
predictor: another variable for `"MAR"` (chosen by `missing_predictor`)
or the variable's own value for `"MNAR"`, with slope `missing_strength`.
The mechanism acts on the drawn (latent) values, so when `categories`
also discretizes the data the missingness is keyed on the underlying
value, not the category code. For `"MAR"` the predictor is evaluated on
the complete drawn values, but every variable is holed at rate
`missing_prop`, so a variable's MAR predictor is itself missing for
roughly a `missing_prop` fraction of the cases whose missingness it
drove. The mechanism is therefore MAR conditional on the *complete*
data, and **not** ignorable for an analyst who sees only the observed
data: estimators that are consistent under ignorable MAR, such as
`cor_method = "fiml"` in
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
and the multiple imputation behind
[`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md),
keep a residual bias here that grows with `missing_prop` and
`missing_strength`. The returned matrix carries the `NA`s, which the
correlation estimators handle downstream.

With `model_error`, the population is perturbed away from the exact
factor structure so the `q`-factor model (`q = ncol(Lambda)`) fits it
only approximately, at a prescribed misfit; exact factor structures are
unrealistic and overstate recovery in simulation studies (MacCallum,
2003). The perturbation is applied once to the population, and the
achieved fit is computed with the same fit-index formulas
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
uses and returned in the `model_error` element. It is applied only when
a target is supplied (`target_rmsea` and/or `target_cfi`), needs a
factor-model population (`Lambda`) with residual degrees of freedom and
an exact factor structure (a diagonal `Psi`), and is orthogonal to the
marginal, ordinal, and missing-data options. Three methods are
available. `"CB"` (Cudeck & Browne, 1992) matches the target RMSEA to
numerical precision and keeps the `q`-factor model the exact minimizer
(the CFI follows as a derived quantity). `"TKL"` (Tucker, Koopman &
Linn, 1969) adds minor common factors tuned so the achieved RMSEA – and,
optionally, CFI – match the target(s); with a single target the match is
close, with both it is a compromise. `"WB"` (Wu & Browne, 2015) draws
the population from an inverse-Wishart distribution around the
model-implied correlation; its calibration applies to the *best-fitting*
model, so the reported misfit of the *generating* model is
systematically larger than the target – by roughly
\\\sqrt{(p(p-1)/2)/df}\\, about 1.4 times for 12 variables and 3
factors. Use `"CB"` when the reported RMSEA must equal the target.
`"CB"` and `"WB"` target the RMSEA only; `"TKL"` can target the RMSEA
and/or the CFI. The reported RMSEA/CFI is the misfit of the specified
generating model.

Replicated draws (`n_datasets > 1`) are generated in parallel across
replicates with future.apply; a parallel plan can be selected with
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
(the default plan runs sequentially). Each replicate is assigned its own
reproducible random-number stream, so with a fixed `seed` the output is
identical regardless of the number of workers.

## References

Cudeck, R., & Browne, M. W. (1992). Constructing a covariance matrix
that yields a specified minimizer and a specified minimum discrepancy
function value. Psychometrika, 57(3), 357-369.
[doi:10.1007/BF02295424](https://doi.org/10.1007/BF02295424)

Fleishman, A. I. (1978). A method for simulating non-normal
distributions. Psychometrika, 43(4), 521-532.
[doi:10.1007/BF02293811](https://doi.org/10.1007/BF02293811)

Foldnes, N., & Grønneberg, S. (2015). How general is the Vale-Maurelli
simulation approach? Psychometrika, 80(4), 1066-1083.
[doi:10.1007/s11336-014-9414-0](https://doi.org/10.1007/s11336-014-9414-0)

Foldnes, N., & Olsson, U. H. (2016). A simple simulation technique for
nonnormal data with prespecified skewness, kurtosis, and covariance
matrix. Multivariate Behavioral Research, 51(2-3), 207-219.
[doi:10.1080/00273171.2015.1133274](https://doi.org/10.1080/00273171.2015.1133274)

MacCallum, R. C. (2003). 2001 Presidential Address: Working with
imperfect models. Multivariate Behavioral Research, 38(1), 113-139.
[doi:10.1207/S15327906MBR3801_5](https://doi.org/10.1207/S15327906MBR3801_5)

Olsson, U. (1979). Maximum likelihood estimation of the polychoric
correlation coefficient. Psychometrika, 44(4), 443-460.
[doi:10.1007/BF02296207](https://doi.org/10.1007/BF02296207)

Olvera Astivia, O. L., & Zumbo, B. D. (2019). A note on the solution
multiplicity of the Vale-Maurelli intermediate correlation equation.
Journal of Educational and Behavioral Statistics, 44(2), 127-143.
[doi:10.3102/1076998618803381](https://doi.org/10.3102/1076998618803381)

Rubin, D. B. (1976). Inference and missing data. Biometrika, 63(3),
581-592.
[doi:10.1093/biomet/63.3.581](https://doi.org/10.1093/biomet/63.3.581)

Ruscio, J., & Kaczetow, W. (2008). Simulating multivariate nonnormal
data using an iterative algorithm. Multivariate Behavioral Research,
43(3), 355-381.
[doi:10.1080/00273170802285693](https://doi.org/10.1080/00273170802285693)

Tucker, L. R., Koopman, R. F., & Linn, R. L. (1969). Evaluation of
factor analytic research procedures by means of simulated correlation
matrices. Psychometrika, 34(4), 421-459.
[doi:10.1007/BF02290601](https://doi.org/10.1007/BF02290601)

Vale, C. D., & Maurelli, V. A. (1983). Simulating multivariate nonnormal
distributions. Psychometrika, 48(3), 465-471.
[doi:10.1007/BF02293687](https://doi.org/10.1007/BF02293687)

Wu, H., & Browne, M. W. (2015). Quantifying adventitious error in a
covariance structure as a random effect. Psychometrika, 80(3), 571-600.
[doi:10.1007/s11336-015-9451-3](https://doi.org/10.1007/s11336-015-9451-3)

## See also

Other data simulation:
[`print.efa_simulated()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_simulated.md)

## Examples

``` r
# Build a population from a shipped loading pattern and factor correlations
Lambda <- population_models$loadings$baseline
Phi <- population_models$phis_3$moderate

# Draw one normal dataset of 500 cases (the data live in $data)
sim <- efa_simulate(N = 500, Lambda = Lambda, Phi = Phi, seed = 42)
dim(sim$data)
#> [1] 500  18

# Return only the population correlation matrix
R_pop <- efa_simulate(Lambda = Lambda, Phi = Phi, return_pop = TRUE)$population

# Draw several datasets at once from a supplied correlation matrix
sims <- efa_simulate(N = 500, R = R_pop, n_datasets = 3, seed = 42)
length(sims$data)
#> [1] 3

# Reproduce the population correlation but with skewed, empirical marginals
# (here from a chi-squared source with one column per variable)
src <- matrix(rchisq(200 * nrow(Lambda), df = 3), ncol = nrow(Lambda))
dat_emp <- efa_simulate(N = 500, Lambda = Lambda, Phi = Phi,
                        marginals = "empirical", marginal_data = src, seed = 42)

# Draw skewed, leptokurtic data with the Vale-Maurelli method
dat_vm <- efa_simulate(N = 500, Lambda = Lambda, Phi = Phi, marginals = "VM",
                       skewness = 1.5, kurtosis = 4, seed = 42)

# Draw five-category ordinal data whose polychoric correlation matches R
dat_ord <- efa_simulate(N = 500, Lambda = Lambda, Phi = Phi,
                        categories = 5, match = "polychoric", seed = 42)

# Draw data with 15% missing at random, driven by a neighbouring item
dat_mar <- efa_simulate(N = 500, Lambda = Lambda, Phi = Phi, missing = "MAR",
                        missing_prop = 0.15, seed = 42)
colMeans(is.na(dat_mar$data))
#>    V1    V2    V3    V4    V5    V6    V7    V8    V9   V10   V11   V12   V13 
#> 0.154 0.146 0.184 0.156 0.152 0.170 0.174 0.138 0.132 0.138 0.128 0.134 0.136 
#>   V14   V15   V16   V17   V18 
#> 0.132 0.152 0.144 0.166 0.140 

# Add realistic model error: a population the model fits with RMSEA of about .05
# (Cudeck-Browne, the default method; the achieved fit is reported)
sim_me <- efa_simulate(N = 500, Lambda = Lambda, Phi = Phi,
                       target_rmsea = 0.05, seed = 42)
sim_me$model_error$rmsea
#> [1] 0.05
```
