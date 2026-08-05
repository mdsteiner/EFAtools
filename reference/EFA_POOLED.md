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

  Numeric in \\(0, 1)\\. One minus the confidence level used for pooled
  Wald-type bootstrap/MI confidence intervals when bootstrap replicates
  are available. For example, `p = .05` gives 95% intervals.

- target_method:

  Character. How rotated solutions are aligned across imputations before
  pooling: `"first_target"` (the default) aligns every imputation to the
  first imputation's rotated solution, while `"consensus"` refines a
  centroid target by Generalized Procrustes Analysis (orthogonal
  rotations only). See *Aligning solutions across imputations* in
  Details.

- align_unrotated:

  Character. How unrotated loadings are aligned before pooling:
  `"signed_tucker_congruence"` (the default; sign/permutation via Tucker
  congruence, anchored on the medoid imputation and returned in the
  extraction's canonical gauge), `"procrustes"` (orthogonal Procrustes
  to the first imputation), or `"none"`. See *Aligning solutions across
  imputations* in Details.

- fit_pool_method:

  Character. Currently only `"D2"` is implemented for chi-square-type
  fit. If no chi-square is available, only residual-based fit and
  descriptive quantities are returned. See *Pooling the model chi-square
  and fit indices* in Details.

- consensus_args:

  List of additional arguments controlling the GPA-consensus iteration
  when `target_method = "consensus"`. Recognised tuning parameters
  include the convergence tolerances `tol` and `loss_tol`, the iteration
  bounds `min_iter` and `max_iter`, the target-update damping `alpha`,
  and the multi-start controls `multi_start` and `starts`.

- procrustes_args:

  List of additional arguments passed to
  [`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md)
  for fixed-target alignment.

- rmsea_ci_level:

  Numeric. Confidence level for the RMSEA CI.

- rmsr_upper:

  **\[deprecated\]** Accepted and ignored. It selected between computing
  RMSR from the unique off-diagonal residual correlations and from the
  full off-diagonal matrix. The two element sets hold each residual pair
  once and twice respectively, so their sums and counts double together
  and the mean square is the same number whenever the residual matrix is
  symmetric, which the pooled residuals always are. RMSR is therefore
  always the root mean square of the unique off-diagonal residuals;
  SRMR, which divides the same sum by the number of non-redundant
  elements, is reported alongside it. Supplying it to
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
  are valid.

## Value

The value of
[`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md),
normally a list of class `c("efa_mi", "EFA_POOLED", "efa", "EFA")`; see
there for the components.

## See also

[`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)
