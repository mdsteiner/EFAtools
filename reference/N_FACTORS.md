# Various factor retention criteria

**\[superseded\]**

`N_FACTORS()` has been superseded by
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md),
which is the recommended interface going forward. It remains available
and unchanged so existing code keeps working.

## Usage

``` r
N_FACTORS(
  x,
  criteria = c("CD", "EKC", "HULL", "MAP", "NEST", "PARALLEL"),
  suitability = TRUE,
  N = NA,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
  n_factors_max = NA,
  N_pop = 10000,
  N_samples = 500,
  alpha = 0.3,
  max_iter_CD = 50,
  n_fac_theor = NA,
  method = c("ML", "PAF", "ULS"),
  gof = c("CAF", "CFI", "RMSEA"),
  eigen_type_HULL = c("SMC", "PCA", "EFA"),
  eigen_type_other = c("SMC"),
  n_factors = 1,
  n_datasets = 1000,
  percent = 95,
  decision_rule = c("means", "percentile", "crawford"),
  ekc_type = lifecycle::deprecated(),
  n_datasets_nest = 1000,
  alpha_nest = 0.05,
  show_progress = FALSE,
  ...
)
```

## Arguments

- x:

  data.frame or matrix. Raw data, or a correlation matrix. If `"CD"` is
  included as a criterion, x must be raw data.

- criteria:

  character. A vector with the factor retention methods to perform.
  Possible inputs are: `"CD"`, `"EKC"`, `"HULL"`, `"KGC"`, `"MAP"`,
  `"NEST"`, `"PARALLEL"`, `"SCREE"`, and `"SMT"` (see the details in
  [`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)).
  By default, a subset of often used, well-performing methods are
  performed.

- suitability:

  logical. Whether the data should be checked for suitability for factor
  analysis using Bartlett's test of sphericity and the
  Kaiser-Meyer-Olkin criterion (see details). Default is `TRUE`.

- N:

  numeric. The number of observations. Only needed if x is a correlation
  matrix.

- use:

  character. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) if raw data is
  given as input. Default is `"pairwise.complete.obs"`.

- cor_method:

  character. Correlation computed from raw data: `"pearson"`,
  `"spearman"`, or `"kendall"` (passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)), or `"poly"` /
  `"tetra"` for polychoric / tetrachoric correlations (a two-step
  estimator). `CD`, `PARALLEL`, `NEST`, `HULL`, and `SMT` do not support
  `"poly"` / `"tetra"` and are skipped automatically if you request them
  together. Default is `"pearson"`.

- n_factors_max:

  numeric. Passed to
  [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md).
  The maximum number of factors to test against. Larger numbers will
  increase the duration the procedure takes, but test more possible
  solutions. If left NA (default), the maximum number of factors for
  which the model is still over-identified (df \> 0) is used.

- N_pop:

  numeric. Passed to
  [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md).
  Size of finite populations of comparison data. Default is 10000.

- N_samples:

  numeric. Passed to
  [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md).
  Number of samples drawn from each population. Default is 500.

- alpha:

  numeric. Passed to
  [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md).
  The alpha level used to test the significance of the improvement added
  by an additional factor. Default is .30.

- max_iter_CD:

  numeric. Passed to
  [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md).
  The maximum number of iterations to perform after which the iterative
  PAF procedure is halted. Default is 50.

- n_fac_theor:

  numeric. Passed to
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md).
  Theoretical number of factors to retain. The Hull method uses one plus
  the larger of this number and the number of factors suggested by
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  as its upper bound.

- method:

  character. The estimator to use in the criteria that fit EFA models;
  passed to
  [`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
  as its `estimator` argument. One of `"ML"`, `"PAF"`, or `"ULS"`.

- gof:

  character. Passed to
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md).
  The goodness of fit index to use. Either `"CAF"`, `"CFI"`, or
  `"RMSEA"`, or any combination of them. With the `"PAF"` estimator,
  only the CAF can be used as goodness of fit index. For details on the
  CAF, see Lorenzo-Seva, Timmerman, and Kiers (2011).

- eigen_type_HULL:

  character. Passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  in
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md).
  What the eigenvalues in the parallel analysis are based on. One of
  `"SMC"`, `"PCA"`, or `"EFA"` – different ways of estimating how much
  variance each indicator shares with the others before the eigenvalues
  are computed. `"SMC"` (default) uses each indicator's squared multiple
  correlation with the others (its diagonal value in the correlation
  matrix). `"PCA"` leaves the diagonal at 1, so each indicator's total
  variance – not just the shared part – feeds into the eigenvalues.
  `"EFA"` uses the communalities from a fitted EFA solution instead.

- eigen_type_other:

  character. Passed to
  [`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
  [`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
  and
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md).
  The same as eigen_type_HULL, but multiple inputs are possible here
  (any combination of `"PCA"`, `"SMC"`, and `"EFA"`). Default is
  `"SMC"`.

- n_factors:

  numeric. Passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  (also within
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)),
  [`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
  and
  [`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md).
  Number of factors to extract if `"EFA"` is included in
  `eigen_type_HULL` or `eigen_type_other`. Default is 1.

- n_datasets:

  numeric. Passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  (also within
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)).
  The number of datasets to simulate. Default is 1000.

- percent:

  numeric. Passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  (also within
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)).
  The percentile to take from the simulated eigenvalues. Default is 95.

- decision_rule:

  character. Passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  (also within
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)).
  Which rule to use to determine the number of factors to retain.
  Default is `"means"`, which uses the average simulated eigenvalues.
  `"percentile"` uses the percentiles specified in `percent`.
  `"crawford"` uses the 95th percentile for the first factor and the
  mean afterwards (based on Crawford et al., 2010).

- ekc_type:

  **\[deprecated\]** Accepted and ignored. It used to select between two
  ways to compute the
  [`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md)
  reference values. The `"AM2019"` reference values do not depend on the
  observed eigenvalues. They therefore skip the empirical correction
  that defines the criterion, so they are no longer computed.

- n_datasets_nest:

  numeric. Passed to
  [`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md).
  The number of datasets to simulate. Default is 1000.

- alpha_nest:

  numeric. Passed to
  [`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md).
  The alpha level to use. The reference values are the eigenvalues at
  the (1 - alpha_nest) percentile. Default is .05.

- show_progress:

  logical. Whether a progress bar should be shown in the console.
  Default is FALSE.

- ...:

  Further arguments passed on to the
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  fits, including the estimation tuning knobs (`type`, `init_comm`,
  `criterion`, `criterion_type`, `abs_eigen`, `start_method`), which are
  repacked into an
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  object so that they tune the fits exactly as they always did. The
  estimator is selected with `method`; `max_iter` is taken by the
  `max_iter_CD` argument (R matches an abbreviated name against the
  arguments before `...`) and so does not reach the fits.

## Value

A list of class `c("efa_retain", "N_FACTORS")`, identical to the value
of
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md);
see there for the components.

## See also

[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
