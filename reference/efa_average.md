# Model averaging across different EFA estimators and types

Not all EFA procedures always arrive at the same solution. This function
allows you perform a number of EFAs from different estimators (e.g.,
Maximum Likelihood and Principal Axis Factoring), with different
implementations (e.g., the SPSS and psych implementations of Principal
Axis Factoring), and across different rotations of the same type (e.g.,
multiple oblique rotations, like promax and oblimin). `efa_average()`
will then run all these EFAs (using the
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
function) and provide a summary across the different solutions.

## Usage

``` r
efa_average(
  x,
  n_factors,
  N = NA,
  estimator = "PAF",
  rotation = "promax",
  type = "none",
  averaging = c("mean", "median"),
  trim = 0,
  salience_threshold = 0.3,
  max_iter = 10000,
  init_comm = c("smc", "mac", "unity"),
  criterion = c(0.001),
  criterion_type = c("sum", "max_individual"),
  abs_eigen = c(TRUE),
  varimax_type = c("svd", "kaiser"),
  normalize = TRUE,
  k_promax = 2:4,
  k_simplimax = ncol(x),
  p_type = c("norm", "unnorm"),
  precision = 1e-05,
  start_method = c("psych", "factanal"),
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra", "fiml"),
  show_progress = TRUE,
  seed = NULL,
  P_type = lifecycle::deprecated()
)
```

## Source

Grieder, S., & Steiner, M. D. (2022). Algorithmic jingle jungle: A
comparison of implementations of principal axis factoring and promax
rotation in R and SPSS. Behavior Research Methods, 54, 54–74. doi:
10.3758/s13428-021-01581-x

Hendrickson, A. E., & White, P. O. (1964). Promax: A quick method for
rotation to oblique simple structure. British Journal of Statistical
Psychology, 17 , 65–70. doi: 10.1111/j.2044-8317.1964.tb00244.x

Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. L. (2011). The Hull
Method for Selecting the Number of Common Factors, Multivariate
Behavioral Research, 46, 340-364, doi: 10.1080/00273171.2011.564527

Kaiser, H. F. (1958). The varimax criterion for analytic rotation in
factor analysis. Psychometrika, 23, 187–200. doi: 10.1007/BF02289233

## Arguments

- x:

  data.frame or matrix. Dataframe or matrix of raw data or matrix with
  correlations. If raw data is entered, the correlation matrix is found
  from the data.

- n_factors:

  numeric. Number of factors to extract.

- N:

  numeric. The number of observations. Needs only be specified if a
  correlation matrix is used. If input is a correlation matrix and `N` =
  NA (default), not all fit indices can be computed.

- estimator:

  character vector. Any combination of "PAF", "ML", and "ULS", to use
  principal axis factoring, maximum likelihood, or unweighted least
  squares, respectively, to fit the EFAs. "MINRES" is accepted as a
  synonym for "ULS" (the same estimator). Default is "PAF". "DWLS",
  which
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  does accept, is deliberately not offered here: it weights each
  residual correlation by the inverse of its asymptotic variance, which
  is only available from raw ordinal data analysed with
  `cor_method = "poly"` or `"tetra"`, whereas every EFA in the grid is
  fitted to the single correlation matrix computed once from `x`. Fit a
  DWLS solution with
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  directly.

- rotation:

  character vector. Either perform no rotation ("none"), any combination
  of orthogonal rotations ("varimax", "equamax", "quartimax", "geominT",
  "bentlerT", and "bifactorT"; using "orthogonal" runs all of these), or
  of oblique rotations ("promax", "oblimin", "quartimin", "simplimax",
  "bentlerQ", "geominQ", and "bifactorQ"; using "oblique" runs all of
  these). Rotation types (no rotation, orthogonal rotations, and oblique
  rotations) cannot be mixed. Default is "promax".

- type:

  character vector. Any combination of "none" (default), "EFAtools",
  "psych", and "SPSS" can be entered. "none" allows the specification of
  various combinations of the arguments controlling both factor
  extraction methods and the rotations. The others ("EFAtools", "psych",
  and "SPSS"), control the execution of the respective factor extraction
  method and rotation to be in line with how it is executed in this
  package (i.e., the respective default procedure), in the psych
  package, and in SPSS. A specific psych implementation exists for PAF,
  ML, varimax, and promax. The SPSS implementation exists for PAF,
  varimax, and promax. For details, see
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md).

- averaging:

  character. One of "mean" (default), and "median". Controls whether the
  different results should be averaged using the (trimmed) mean, or the
  median.

- trim:

  numeric. If averaging is set to "mean", this argument controls the
  trimming of extremes (for details see
  [`base::mean()`](https://rdrr.io/r/base/mean.html)). By default no
  trimming is done (i.e., trim = 0).

- salience_threshold:

  numeric. The threshold to use to classify a pattern coefficient or
  loading as salient (i.e., substantial enough to assign it to a
  factor). Default is 0.3. Indicator-to-factor correspondences will be
  inferred based on this threshold. Note that this may not be meaningful
  if rotation = "none" and n_factors \> 1 are used, as no simple
  structure is present there.

- max_iter:

  numeric. The maximum number of iterations to perform after which the
  iterative PAF procedure is halted with a warning. Default is 10,000.
  It is only evaluated for the "PAF" solutions run under `type` "none":
  a named `type` brings the iteration cap that defines it ("SPSS" 25,
  "psych" 50, and "EFAtools" 300), and "ML" and "ULS" do not iterate
  this way. Note that non-converged procedures are excluded from the
  averaging procedure.

- init_comm:

  character vector. Any combination of "smc", "mac", and "unity".
  Controls the methods to estimate the initial communalities in `PAF` if
  "none" is among the specified types. "smc" will use squared multiple
  correlations, "mac" will use maximum absolute correlations, "unity"
  will use 1s (for details see
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)).
  Default is `c("smc", "mac", "unity")`.

- criterion:

  numeric vector. The convergence criterion used for PAF if "none" is
  among the specified types. If the change in communalities from one
  iteration to the next is smaller than this criterion the solution is
  accepted and the procedure ends. Default is `0.001`.

- criterion_type:

  character vector. Any combination of "max_individual" and "sum". Type
  of convergence criterion used for PAF if "none" is among the specified
  types. "max_individual" selects the maximum change in any of the
  communalities from one iteration to the next and tests it against the
  specified criterion. "sum" takes the difference of the sum of all
  communalities in one iteration and the sum of all communalities in the
  next iteration and tests this against the criterion (for details see
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)).
  Default is `c("sum", "max_individual")`.

- abs_eigen:

  logical vector. Any combination of TRUE and FALSE. Which algorithm to
  use in the PAF iterations if "none" is among the specified types. If
  FALSE, the loadings are computed from the eigenvalues. This is also
  used by the [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html)
  function. If TRUE the loadings are computed with the absolute
  eigenvalues as done by SPSS (for details see
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)).
  Default is `TRUE`.

- varimax_type:

  character vector. Any combination of "svd" and "kaiser". The type of
  the varimax rotation performed if "none" is among the specified types
  and "varimax", "promax", "orthogonal", or "oblique" is among the
  specified rotations. "svd" uses singular value decomposition, as
  [`stats::varimax()`](https://rdrr.io/r/stats/varimax.html) does, and
  "kaiser" uses the varimax procedure performed in SPSS. This is the
  original procedure from Kaiser (1958), but with slight alterations in
  the varimax criterion (for details, see
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  and Grieder & Steiner, 2022). Default is `c("svd", "kaiser")`.

- normalize:

  logical vector. Any combination of TRUE and FALSE. `TRUE` performs a
  kaiser normalization before the specified rotation(s). Default is
  `TRUE`.

- k_promax:

  numeric vector. The power used for computing the target matrix P in
  the promax rotation if "none" is among the specified types and
  "promax" or "oblique" is among the specified rotations. Default is
  `2:4`.

- k_simplimax:

  numeric. The number of 'close to zero loadings' for the simplimax
  rotation if "simplimax" or "oblique" is among the specified rotations.
  Default is `ncol(x)`, where x is the entered data.

- p_type:

  character vector. Any combination of "norm" and "unnorm". This
  specifies how the target matrix P is computed in promax rotation if
  "none" is among the specified types and "promax" or "oblique" is among
  the specified rotations. "unnorm" will use the unnormalized target
  matrix as originally done in Hendrickson and White (1964). "norm" will
  use a normalized target matrix (for details see
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)).
  Default is `c("norm", "unnorm")`.

- precision:

  numeric vector. The tolerance for stopping in the rotation
  procedure(s). Default is 10^-5.

- start_method:

  character vector. Any combination of "psych" and "factanal". How to
  specify the starting values for the optimization procedure for ML.
  "psych" takes the starting values specified in
  [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html). "factanal"
  takes the starting values specified in the
  [`stats::factanal()`](https://rdrr.io/r/stats/factanal.html) function.
  Default is `c("psych", "factanal")`.

- use:

  character. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) if raw data is
  given as input. Default is "pairwise.complete.obs". It is ignored when
  `cor_method = "fiml"`, which handles the missingness itself, so every
  case contributes.

- cor_method:

  character. Correlation computed from raw data: `"pearson"`,
  `"spearman"`, or `"kendall"` (passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)), `"poly"` /
  `"tetra"` for polychoric / tetrachoric correlations of ordinal /
  binary data (a two-step estimator), or `"fiml"` for a two-stage
  full-information maximum-likelihood correlation from raw data with
  missing values. With `"fiml"` the saturated multivariate-normal mean
  and covariance are estimated by an EM algorithm assuming the data are
  missing at random and the standardized covariance is analysed,
  reproducing
  [`psych::corFiml()`](https://rdrr.io/pkg/psych/man/corFiml.html)
  followed by [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html) and
  `lavaan(missing = "two.stage")`, *not* `lavaan::efa(missing = "ml")`
  (see
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  and the details). Default is "pearson".

- show_progress:

  logical. Whether a progress bar should be shown in the console.
  Default is TRUE.

- seed:

  numeric or `NULL`. An optional seed making the averaging run
  reproducible and independent of the number of parallel workers (the
  grid runs with
  [future_lapply()](https://future.apply.futureverse.org/reference/future_lapply.html),
  for which a parallel plan can be set via
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)).
  It matters whenever a criterion-based rotation is included, since
  those draw random starts. When supplied, the caller's random-number
  stream is restored afterwards, leaving no side effect. Default is
  `NULL`.

- P_type:

  **\[superseded\]** Former name of `p_type`. Still accepted (silently)
  for backwards compatibility; please use `p_type`.

## Value

A list of class `c("efa_average", "EFA_AVERAGE")` containing the
components below. Throughout, `range` is the width `maximum - minimum`
of each cell across the factor solutions, not the interval
`[minimum, maximum]`, and `average` is the (trimmed) mean or the median,
following `averaging`.

- orig_R:

  Original correlation matrix.

- h2:

  A list with the average, standard deviation, minimum, maximum, and
  range of the final communality estimates across the factor solutions.

- loadings:

  A list with the average, standard deviation, minimum, maximum, and
  range of the final loadings across the factor solutions. If rotation
  was "none", the unrotated loadings, otherwise the rotated loadings
  (pattern coefficients).

- Phi:

  A list with the average, standard deviation, minimum, maximum, and
  range of the factor intercorrelations across factor solutions obtained
  with oblique rotations.

- ind_fac_corres:

  A matrix with each cell containing the proportion of the factor
  solutions in which the respective indicator-to-factor correspondence
  occurred, i.e., in which the loading exceeded the specified salience
  threshold. Note: Rowsums can exceed 1 due to cross-loadings.

- vars_accounted:

  A list with the average, standard deviation, minimum, maximum, and
  range of explained variances and sums of squared loadings across the
  factor solutions. Based on the unrotated loadings if rotation was
  "none" or only one factor was extracted, otherwise on the rotated
  loadings.

- fit_indices:

  A matrix containing the average, standard deviation, minimum, maximum,
  and range for all applicable fit indices across the respective factor
  solutions, and the degrees of freedom (df). If the estimator argument
  contains ML or ULS: Fit indices derived from the unrotated factor
  loadings: Chi Square (chisq), including significance level,
  Comparative Fit Index (CFI), Tucker-Lewis Index (TLI), Root Mean
  Square Error of Approximation (RMSEA), Akaike Information Criterion
  (AIC), Bayesian Information Criterion (BIC), Expected Cross-Validation
  Index (ECVI), and the common part accounted for (CAF) index as
  proposed by Lorenzo-Seva, Timmerman, & Kiers (2011). The
  residual-based Standardized Root Mean Square Residual (SRMR) and Root
  Mean Square Residual (RMSR) and the CAF are also computed for PAF; for
  PAF the remaining (chi-square-based) indices are not available (see
  details).

- implementations_grid:

  A matrix containing, for each performed EFA, the setting combination,
  if an error occurred (logical), the error message (character), an
  integer convergence code (0 = converged; for ML and ULS the same codes
  as [`stats::optim()`](https://rdrr.io/r/stats/optim.html)'s
  "L-BFGS-B", for PAF 1 if the maximum number of iterations was reached
  without meeting the convergence criterion and 0 otherwise), if heywood
  cases occurred (logical, see details for definition), if the solution
  was admissible (logical, see details for definition), and the fit
  indices.

- efa_list:

  A list containing the outputs of all performed EFAs. The names
  correspond to the rownames from the implementations_grid.

- settings:

  A list of the settings used.

If the supplied arguments admit only a single EFA, there is nothing to
average across: that one
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
object is returned instead, with a warning.

## Details

As a first step in this function, a grid is produced containing the
setting combinations for the to-be-performed EFAs. These settings are
then entered as arguments to the
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
function and the EFAs are run in a second step. After all EFAs are run,
the factor solutions are averaged and their variability determined in a
third step.

When raw data are supplied, the correlation matrix is computed once
before the grid is run and reused for every EFA in it. Under
`cor_method = "fiml"` this means the saturated multivariate-normal
moments are EM-estimated a single time (from the raw data with missing
values, assuming the data are missing at random) and the resulting
two-stage correlation is analysed by every solution in the grid; the EM
is not re-run per solution. Under `cor_method = "fiml"`, `use` does not
select cases (every case contributes to the EM). The averaged loadings
and communalities are then the two-stage FIML estimates, but the
averaged Chi-Square and the indices derived from it (CFI, TLI, RMSEA,
AIC, BIC, ECVI) are the ordinary ML/ULS discrepancy statistics on the EM
correlation, not the corrected two-stage (Satorra-Bentler) statistics
that a standalone
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
with `cor_method = "fiml"` reports; in particular the averaged AIC and
BIC are finite here rather than `NA`.

The grid containing the setting combinations is produced based on the
entries to the respective arguments. To this end, all possible
combinations resulting in unique EFA models are considered: combinations
that resolve to the same model are run only once. Two combinations are
the same model only if every setting the fit consumes agrees, and for
"PAF" that includes the iteration cap `max_iter`. Since a named `type`
brings its own cap, a `type` of `c("none", "SPSS")` whose specific
settings match the SPSS combination in every other respect still gives
two "PAF" models, unless `max_iter` is also set to the cap of the SPSS
implementation. We include here a list of arguments that are only
evaluated under specific conditions:

The arguments `init_comm`, `criterion`, `criterion_type`, `abs_eigen`,
and `max_iter` are only evaluated if "PAF" is included in `estimator`
and "none" is included in `type`.

The argument `varimax_type` is only evaluated if "varimax", "promax",
"oblique", or "orthogonal" is included in `rotation` and "none" is
included in `type`.

The argument `normalize` is only evaluated if `rotation` is not set to
"none" and "none" is included in `type`.

The argument `k_simplimax` is only evaluated if "simplimax" or "oblique"
is included in `rotation`.

The arguments `k_promax` and `p_type` are only evaluated if "promax" or
"oblique" is included in `rotation` and "none" is included in `type`.

The argument `start_method` is only evaluated if "ML" is included in
`estimator`.

To avoid a bias in the averaged factor solutions from problematic
solutions, these are excluded prior to averaging. A solution is deemed
problematic if at least one of the following is true: an error occurred,
the model did not converge, or there is at least one Heywood (improper)
case (a communality at or above 1, or, for `ML`/`ULS`, a uniqueness
pinned at the estimator's lower bound). Information on errors,
convergence, and Heywood cases are returned in the implementations_grid
and a summary of these is given when printing the output. In addition to
these, information on the admissibility of the factor solutions is also
included. A solution was deemed admissible if (1) no error occurred, (2)
the model converged, (3) no Heywood cases are present, and (4) there are
at least two salient loadings (i.e., loadings exceeding the specified
`salience_threshold`) for each factor. So, solutions failing one of the
first three of these criteria of admissibility are also deemed
problematic and therefore excluded from averaging. However, solutions
failing only the fourth criterion of admissibility are still included
for averaging. Finally, if all solutions are problematic (e.g., all
solutions contain Heywood cases), no averaging is performed and the
respective outputs are NA. In this case, the implementations_grid should
be inspected to see if there are any error messages, and the separate
EFA solutions that are also included in the output can be inspected as
well, for example, to see where Heywood cases occurred.

A core output of this function includes the average, minimum, and
maximum loadings derived from all non-problematic (see above) factor
solutions. Please note that these are not entire solutions, but the
matrices include the average, minimum, or maximum value for each cell
(i.e., each loading separately). This means that, for example, the
matrix with the minimum loadings will contain the minimum value in any
of the factor solutions for each specific loading, and therefore most
likely contains loadings from different factor solutions. The matrices
containing the minimum and maximum factor solutions can therefore not be
interpreted as whole factor solutions.

The averaged loading matrix is likewise a cell-wise summary rather than
a fitted solution: it is not itself the solution of any EFA, does not in
general reproduce the correlation matrix, and need not reproduce the
averaged communalities. The fit indices described below are
correspondingly the mean (or, under `averaging = "median"`, the median)
of the per-solution fit indices, not the fit of the averaged loadings,
so the averaged loadings and the reported fit do not describe one and
the same model.

The output also includes information on the average, minimum, maximum,
and variability of the fit indices across the non-problematic factor
solutions. It is important to note that not all fit indices are computed
for all fit methods: For ML and ULS, all fit indices can be computed,
while for PAF the chi-square-based indices (the chi-square statistic and
its significance, CFI, TLI, RMSEA, AIC, BIC, and ECVI) are NA. The
common part accounted for (CAF) index (Lorenzo-Seva, Timmerman, & Kiers,
2011) and the residual-based SRMR and RMSR are still computed for PAF.
As a consequence, if only "PAF" is included in the `estimator` argument,
averaging is performed for the CAF, SRMR, and RMSR, while the
chi-square-based indices are NA. If a combination of "PAF" and "ML"
and/or "ULS" are included in the `estimator` argument, the CAF, SRMR,
and RMSR are averaged across all non-problematic factor solutions, while
the chi-square-based indices are only averaged across the ML and ULS
solutions. The user should therefore keep in mind that the number of
EFAs across which the fit indices are averaged can diverge for the CAF,
SRMR, and RMSR compared to the chi-square-based indices.

Each reported fit index is summarised across the (non-problematic)
solutions in the same descriptive way: the average, standard deviation,
minimum, and maximum of the per-solution values. This includes the
chi-square significance level (`p_chi`), which is therefore the mean (or
median) of the per-solution p-values and is purely descriptive; it is
not the p-value of any pooled chi-square test.

## See also

Other factor analysis:
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md),
[`efa_group()`](https://mdsteiner.github.io/EFAtools/reference/efa_group.md),
[`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md),
[`plot.efa_group()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_group.md),
[`print.efa_group()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_group.md)

## Examples

``` r
# Averaging across one implementation each of PAF (EFAtools type), ULS, and
# ML with one implementation of promax (EFAtools type) (3 EFAs)
Aver_meth <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                         estimator = c("PAF", "ULS", "ML"), type = "EFAtools",
                         start_method = "psych")
#> ℹ Extracting data
#> ✔ Extracting data [9ms]
#> 
#> ℹ Reordering factors
#> ✔ Reordering factors [25ms]
#> 
#> ℹ Averaging data
#> ✔ Averaging data [21ms]
#> 

# \donttest{
# Averaging across different implementations of PAF and promax rotation (72 EFAs)
Aver_PAF <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500)
#> ℹ Extracting data
#> ✔ Extracting data [10ms]
#> 
#> ℹ Reordering factors
#> ✔ Reordering factors [18ms]
#> 
#> ℹ Averaging data
#> ✔ Averaging data [15ms]
#> 

# Use median instead of mean for averaging (72 EFAs)
Aver_PAF_md <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                           averaging = "median")
#> ℹ Extracting data
#> ✔ Extracting data [10ms]
#> 
#> ℹ Reordering factors
#> ✔ Reordering factors [17ms]
#> 
#> ℹ Averaging data
#> ✔ Averaging data [18ms]
#> 

# Averaging across different implementations of PAF and promax rotation,
# and across ULS and different versions of ML (108 EFAs)
Aver_meth_ext <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                             estimator = c("PAF", "ULS", "ML"))
#> ℹ Extracting data
#> ✔ Extracting data [12ms]
#> 
#> ℹ Reordering factors
#> ✔ Reordering factors [27ms]
#> 
#> ℹ Averaging data
#> ✔ Averaging data [15ms]
#> 

# Averaging across different oblique rotation methods, using one implementation
# of ML and one implementation of promax (EFAtools type) (7 EFAs)
Aver_rot <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                         estimator = "ML", rotation = "oblique", type = "EFAtools",
                         start_method = "psych")
#> ℹ Extracting data
#> ✔ Extracting data [10ms]
#> 
#> ℹ Reordering factors
#> ✔ Reordering factors [17ms]
#> 
#> ℹ Averaging data
#> ✔ Averaging data [20ms]
#> 
# }

# \donttest{
# Two-stage FIML correlations from raw data with missing values: the EM
# saturated moments are estimated once and the resulting correlation is
# averaged across the grid of EFAs.
x_miss <- GRiPS_raw
x_miss[cbind(1:20, 1)] <- NA
Aver_fiml <- efa_average(x_miss, n_factors = 1, estimator = c("PAF", "ML"),
                         cor_method = "fiml")
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
#> ℹ `n_factors` is 1 but `rotation` is not "none"; setting `rotation` to "none",
#>   as single-factor solutions cannot be rotated.
#> ℹ Extracting data
#> ✔ Extracting data [7ms]
#> 
#> ℹ Averaging data
#> ✔ Averaging data [14ms]
#> 
# }
```
