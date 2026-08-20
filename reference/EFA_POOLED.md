# Exploratory factor analysis on multiple data imputations

**\[superseded\]**

`EFA_POOLED()` has been superseded by
[`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md),
which is the recommended interface going forward. It remains available
and unchanged so existing code keeps working.

## Usage

``` r
EFA_POOLED(
  data_list,
  p = 0.05,
  target_method = c("first_target", "consensus"),
  align_unrotated = c("signed_tucker_congruence", "none", "procrustes"),
  fit_pool_method = c("D2"),
  consensus_args = list(),
  procrustes_args = list(),
  rmsea_ci_level = 0.9,
  rmsr_upper = TRUE,
  ...
)
```

## Arguments

- data_list:

  A list of length \\m\\, where \\m\\ is the number of imputations. Each
  list element is a data frame or matrix of raw data, or a correlation
  matrix. See argument `x` in
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md).
  A `mids` object from mice must be converted first, with
  `mice::complete(x, "all")`.

- p:

  Numeric in \\(0, 1)\\. One minus the confidence level for the pooled
  confidence intervals, whichever `se` method produced them
  (`"information"`, `"np-boot"`, or `"sandwich"`). For example,
  `p = .05` gives 95% intervals.

- target_method:

  Character. How rotated solutions are aligned across imputations before
  pooling: `"first_target"` (the default) aligns every imputation to the
  first imputation's rotated solution, while `"consensus"` refines a
  centroid target by Generalized Procrustes Analysis, started from the
  medoid imputation so that the pooled rotated solution does not depend
  on the order of `data_list` (orthogonal rotations only). See *Aligning
  solutions across imputations* in Details.

- align_unrotated:

  Character. How unrotated loadings are aligned before pooling:
  `"signed_tucker_congruence"` (the default; sign/permutation via Tucker
  congruence, anchored on the medoid imputation and returned in the
  extraction's canonical gauge), `"procrustes"` (orthogonal Procrustes
  to the first imputation), or `"none"`. See *Aligning solutions across
  imputations* in Details.

- fit_pool_method:

  Character. Only `"D2"` is implemented for pooling chi-square-type fit.
  If no chi-square is available, only residual-based fit and descriptive
  quantities are returned. See *Pooling the model chi-square and fit
  indices* in Details.

- consensus_args:

  List of additional arguments controlling the GPA-consensus iteration
  when `target_method = "consensus"`. Recognised tuning parameters
  include the convergence tolerances `tol` and `loss_tol`, the iteration
  bounds `min_iter` and `max_iter`, the target-update damping `alpha`,
  the multi-start controls `multi_start` and `starts`, and `start`,
  which overrides the medoid imputation the iteration is otherwise
  started from.

- procrustes_args:

  List of
  [`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md)
  algorithm controls for fixed-target alignment, for example
  `oblique_maxit` or `oblique_random_starts`. The loadings `A`, the
  alignment `Target`, the `rotation` family, and the cross-product `S`
  are derived from the imputations and cannot be set here.

- rmsea_ci_level:

  Numeric. Confidence level for the RMSEA CI.

- rmsr_upper:

  **\[deprecated\]** Deprecated and ignored.
  [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)
  now always computes RMSR the same way, from the unique off-diagonal
  residuals; SRMR is reported alongside it. Supplying it to
  [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)
  signals a deprecation warning; the superseded `EFA_POOLED()` accepts
  it silently.

- ...:

  Additional arguments passed to
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  (e.g. `estimator`, `rotation`, `se`, `n_factors`, `N`). These select
  the estimator, rotation, standard-error method, and fit indices used
  for every imputation; see
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  for the available options, their properties, and which combinations
  are valid. Two of them shape the pooled object rather than a single
  fit: `seed` sets the random state once for the whole
  [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)
  call – every component bootstrap and every random-start rotation draws
  from it, so a seeded call is reproducible as a whole, and the caller's
  random stream is restored afterwards – and `b_boot` sets the number of
  bootstrap replicates drawn per imputation under `se = "np-boot"`,
  which is what the pooled within-imputation variances are estimated
  from and is recorded in `settings$b_boot`. The
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  and
  [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  objects are accepted through `...` as well, although they are not
  declared formals: pass them as `estimate_control =` /
  `rotate_control =` exactly as you would to
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md).

## Value

The value of
[`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md),
normally a list of class `c("efa_mi", "EFA_POOLED", "efa", "EFA")`; see
there for the components.

## See also

[`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)
