# Model averaging across different EFA methods and types

**\[superseded\]**

`EFA_AVERAGE()` has been superseded by
[`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md),
which is the recommended interface going forward. It remains available
and unchanged so existing code keeps working.

## Usage

``` r
EFA_AVERAGE(
  x,
  n_factors,
  N = NA,
  method = "PAF",
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
  P_type = c("norm", "unnorm"),
  precision = 1e-05,
  start_method = c("psych", "factanal"),
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra", "fiml"),
  show_progress = TRUE
)
```

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

- method:

  character vector. Any combination of `"PAF"`, `"ML"`, and `"ULS"`, the
  estimators to average across; passed to
  [`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md)
  as its `estimator` argument. Default is `"PAF"`.

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

- P_type:

  character vector. Any combination of `"norm"` and `"unnorm"`. How the
  promax target matrix P is computed if `"none"` is among the specified
  types and `"promax"` or `"oblique"` is among the specified rotations:
  `"unnorm"` uses the unnormalized target matrix of Hendrickson and
  White (1964), `"norm"` a normalized one. This frozen argument keeps
  its original name;
  [`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md)
  takes the same setting as `p_type`. Default is `c("norm", "unnorm")`.

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

## Value

The value of
[`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md),
normally a list of class `c("efa_average", "EFA_AVERAGE")`; see there
for the components.

## See also

[`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md)
