# Exploratory factor analysis (EFA)

**\[superseded\]**

`EFA()` has been superseded by
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md),
which is the recommended interface going forward.
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
keeps the primary choices (data, factors, estimator, rotation, standard
errors) as top-level arguments and collects the estimation and rotation
tuning knobs into two control objects,
[`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
and
[`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md).
`EFA()` remains available and unchanged – its full flat argument list
still works exactly as before – so existing code keeps running.

## Usage

``` r
EFA(
  x,
  n_factors,
  N = NA,
  method = c("PAF", "ML", "ULS", "MINRES", "DWLS"),
  rotation = c("none", "varimax", "equamax", "quartimax", "geominT", "bentlerT",
    "bifactorT", "promax", "oblimin", "quartimin", "simplimax", "bentlerQ", "geominQ",
    "bifactorQ"),
  se = c("none", "information", "sandwich", "np-boot"),
  type = c("EFAtools", "psych", "SPSS", "none"),
  max_iter = NA,
  init_comm = NA,
  criterion = NA,
  criterion_type = NA,
  abs_eigen = NA,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  varimax_type = NA,
  k = NA,
  normalize = TRUE,
  p_type = NA,
  precision = 1e-05,
  order_type = NA,
  start_method = "psych",
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra", "fiml"),
  b_boot = 1000,
  ci = 0.95,
  random_starts = 100,
  seed = NULL,
  P_type = lifecycle::deprecated(),
  randomStarts = lifecycle::deprecated(),
  ...
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
  NA (default), not all fit indices can be computed. When raw data with
  missing values are entered and `use` is `"complete.obs"` or
  `"na.or.complete"`, rows are deleted listwise, so `N` is taken as the
  number of complete cases.

- method:

  character. The estimator used to fit the EFA: "PAF" (principal axis
  factoring), "ML" (maximum likelihood), "ULS" (unweighted least
  squares; "MINRES" is an accepted alias returning identical results),
  or "DWLS" (diagonally weighted least squares, for ordinal data). See
  the *Estimators* section in Details for their properties and data
  requirements.

- rotation:

  character. Either perform no rotation ("none"; default), an orthogonal
  rotation ("varimax", "equamax", "quartimax", "geominT", "bentlerT", or
  "bifactorT"), or an oblique rotation ("promax", "oblimin",
  "quartimin", "simplimax", "bentlerQ", "geominQ", or "bifactorQ"). See
  the *Rotations* section in Details for their properties and known
  issues.

- se:

  character. Whether and how to compute standard errors (and matching
  confidence intervals): "none" (default, no standard errors),
  "information" (analytic standard errors from the expected Fisher
  information of the ML solution), "sandwich" (robust Godambe sandwich
  standard errors from raw data), or "np-boot" (non-parametric
  bootstrap). The methods differ in their assumptions, their data
  requirements, and which estimator, rotation, and `cor_method`
  combinations they support; see the *Standard errors* section in
  Details.

- type:

  character. If one of "EFAtools" (default), "psych", or "SPSS" is used,
  and the following arguments with default NA are left with NA, these
  implementations are executed according to the respective program
  ("psych" and "SPSS") or according to the best solution found in
  Grieder & Steiner (2022; "EFAtools"). Individual properties can be
  adapted using one of the three types and specifying some of the
  following arguments. If set to "none" additional arguments must be
  specified depending on the `method` and `rotation` used (see details).

- max_iter:

  numeric. The maximum number of iterations to perform after which the
  iterative PAF procedure is halted with a warning. If `type` is one of
  "EFAtools", "SPSS", or "psych", this is automatically specified if
  `max_iter` is left to be `NA`, but can be overridden by entering a
  number. Default is `NA`.

- init_comm:

  character. The method to estimate the initial communalities in `PAF`.
  "smc" will use squared multiple correlations, "mac" will use maximum
  absolute correlations, "unity" will use 1s (see details). Default is
  `NA`.

- criterion:

  numeric. The convergence criterion used for PAF. If the change in
  communalities from one iteration to the next is smaller than this
  criterion the solution is accepted and the procedure ends. Default is
  `NA`.

- criterion_type:

  character. Type of convergence criterion used for PAF.
  "max_individual" selects the maximum change in any of the
  communalities from one iteration to the next and tests it against the
  specified criterion. This is also used by SPSS. "sum" takes the
  difference of the sum of all communalities in one iteration and the
  sum of all communalities in the next iteration and tests this against
  the criterion. This procedure is used by the
  [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html) function.
  Default is `NA`.

- abs_eigen:

  logical. Which algorithm to use in the PAF iterations. If FALSE, the
  loadings are computed from the eigenvalues. This is also used by the
  [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html) function. If
  TRUE the loadings are computed with the absolute eigenvalues as done
  by SPSS. Default is `NA`.

- use:

  character. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) if raw data is
  given as input. Default is "pairwise.complete.obs".

- varimax_type:

  character. The type of the varimax rotation performed. If "svd",
  singular value decomposition is used, as
  [`stats::varimax()`](https://rdrr.io/r/stats/varimax.html) does. If
  "kaiser", the varimax procedure performed in SPSS is used, following
  the original procedure from Kaiser (1958) (see details). Default is
  `NA`.

- k:

  numeric. Either the power used for computing the target matrix P in
  the promax rotation or the number of 'close to zero loadings' for the
  simplimax rotation. If left to `NA` (default), the value for promax
  depends on the specified type. For simplimax, `nrow(L)`, where L is
  the matrix of unrotated loadings, is used by default.

- normalize:

  logical. If `TRUE`, a kaiser normalization is performed before the
  specified rotation. Default is `TRUE`.

- p_type:

  character. This specifies how the target matrix P is computed in
  promax rotation. If "unnorm" it will use the unnormalized target
  matrix as originally done in Hendrickson and White (1964). This is
  also used in the psych and stats packages. If "norm" it will use the
  normalized target matrix as used in SPSS. Default is `NA`.

- precision:

  numeric. The tolerance for stopping in the rotation procedure. Default
  is 10^-5 for all rotation methods.

- order_type:

  character. How to order the factors. "eigen" reorders the factors by
  descending explained variance; "ss_factors" reorders the factors by
  descending (unweighted) sum of squared factor loadings per factor.
  Default is `NA`.

- start_method:

  character. How to specify the starting values for the optimization
  procedure for ML. Default is "psych" which takes the starting values
  specified in [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html).
  "factanal" takes the starting values specified in the
  [`stats::factanal()`](https://rdrr.io/r/stats/factanal.html) function.

- cor_method:

  character. How the correlation is computed from raw data: `"pearson"`,
  `"spearman"`, or `"kendall"` (passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)); `"poly"` /
  `"tetra"` for polychoric / tetrachoric correlations of ordinal /
  binary data; or `"fiml"` for a two-stage full-information
  maximum-likelihood correlation from raw data with missing values. See
  the *Correlation methods* section in Details for their properties and
  the combinations they support. Default is "pearson".

- b_boot:

  numeric. The number of bootstrap samples to draw. Default is 1000.
  Under `cor_method = "fiml"` each bootstrap sample re-runs the EM
  moment estimation, so a smaller value may be advisable.

- ci:

  numeric. The confidence interval to create from the bootstrap samples.
  Must be between 0 and 1. Default is .95 for 95% CIs.

- random_starts:

  numeric. The number of random starts to use in the rotation to guard
  against local minima. Default is 100.

- seed:

  numeric. An optional seed for the random-number generator used by the
  non-parametric bootstrap (`se = "np-boot"`), i.e. for the case
  resampling, the rotation random starts, and the Procrustes random
  starts. Setting it makes the bootstrap reproducible and independent of
  the number of parallel workers (see Details); the caller's
  random-number stream is restored afterwards, so supplying a seed
  leaves no lasting effect on it. Default is `NULL`, which uses (and
  advances) the current state of the generator.

- P_type, randomStarts:

  **\[superseded\]** Former names of `p_type` and `random_starts`. Still
  accepted (silently) for backwards compatibility; please use the new
  names.

- ...:

  Additional arguments passed to the rotation procedure (e.g., `maxit`
  for the maximum number of iterations).

## Value

The value of
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md),
a list of class `c("efa", "EFA")`; see there for the components.

## See also

[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md),
[`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md),
[`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
