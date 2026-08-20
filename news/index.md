# Changelog

## EFAtools 1.0.0.9000

### Comparing and Averaging Solutions

- [`efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/efa_compare.md)
  display settings (as set via
  [`print()`](https://rdrr.io/r/base/print.html),
  [`format()`](https://rdrr.io/r/base/format.html), or
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html)) are now
  saved and reused, so you can change them after the fact without
  rerunning the comparison.

- [`efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/efa_compare.md)’s
  `are_equal` is now `NA` when the compared integer parts differ; `0` is
  reserved for values whose integer parts agree but decimals do not.

- Fixed
  [`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md)’s
  progress bar, which could previously trigger spurious warnings.

- [`efa_group()`](https://mdsteiner.github.io/EFAtools/reference/efa_group.md)
  now declares `se` as a formal argument; a supplied value is dropped
  with a warning pointing to `b_boot` (previously it was silently
  matched to `seed` and raised a confusing error).

- [`efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/efa_compare.md)’s
  elementwise differences are now signed (previously unsigned).

- `efa_compare(reorder = "congruence")` now also aligns the sign of
  single-factor solutions, so two one-factor solutions differing only in
  sign are no longer reported as maximally different.

- Fixed
  [`efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/efa_compare.md)’s
  handling of missing loadings, which previously inflated `diff_corres`
  and `diff_corres_cross`; both are now `NA` when no row can be
  compared.

- [`efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/efa_compare.md)
  now rejects invalid input: empty pairs, infinite values, non-numeric
  content, or a negative `thresh`.

- [`efa_group()`](https://mdsteiner.github.io/EFAtools/reference/efa_group.md)
  now rejects `b_boot = 1`.

### Data Screening and Simulation

- Fixed
  [`efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/efa_screen.md)’s
  multivariate-outlier check, which could previously discard the best
  solution when a robust subset’s covariance matrix was singular.

- Fixed
  [`efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/efa_screen.md)’s
  robust covariance estimation, which previously fell back to classical
  Mahalanobis distances when variable scales differed greatly.

- [`efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/efa_screen.md)’s
  multicollinearity verdict is now based on the condition index alone
  (Belsley, Kuh & Welsch 1980 thresholds of 10/30 unchanged); the
  determinant is reported as a plain number, since a small determinant
  no longer by itself flags multicollinearity.

- [`efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/efa_screen.md)
  now withholds the Henze-Zirkler p-value (previously an uninformative
  exact 0 or 1) when its null distribution has no resolvable spread; the
  `hz` result then carries class `efa_screen_no_hz`.

- Mardia’s skewness and kurtosis in
  [`efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/efa_screen.md)
  now always apply small-sample/exact-moment corrections (previously
  only below 20 observations, or asymptotic); both statistics, their
  p-values, and the normality verdict can change. The
  `mardia$small_sample` flag is dropped.

- [`efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/efa_screen.md)
  now requires at least 3 variables (error `efa_screen_too_few_vars`).

- [`efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/efa_screen.md)’s
  `outlier_cutoff` must now lie between 0.5 and 0.9999 (default 0.975
  unchanged).

- [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md)
  gains `missing_vars` to restrict missingness to selected variables.

- With `marginals = "VM"` and `force_pd = TRUE`,
  [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md)’s
  returned `population` now reflects the projected correlations actually
  used to generate the data.

- [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md)
  now gives a clear error when a requested empirical-marginal category
  is too rare to reproduce at the sample size.

- [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md)
  and `efa_power(mode = "simulation")` now reject an `R` whose diagonal
  is not 1.

### Factor Retention

- The `type` argument of
  [`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md)
  and `ekc_type` of
  [`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
  are deprecated and ignored: the empirical Kaiser criterion is now
  always computed as in Braeken and van Assen (2017), so suggested
  factor counts can differ from earlier results, and
  [`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
  no longer reports `EKC_AM2019`.

- [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md)
  now gives a clear error when data contain constant variables.

- [`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md)’s
  RMSEA rule now stops when the lower confidence bound cannot be
  computed, instead of continuing to later models.

- [`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md)
  results can now be plotted.

- Factor-retention functions now reject the unsupported arguments
  `seed`, `se`, `b_boot`, and `ci` (use
  [`set.seed()`](https://rdrr.io/r/base/Random.html) for randomness).

- [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  now shares simulated datasets between `"PCA"` and `"SMC"` references
  when both are requested, halving simulation time; results for a given
  seed can differ from before.

- A failed simulation block in
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  is now redrawn on its own instead of restarting the whole batch.

- [`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md),
  [`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md),
  and
  [`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md)
  now require a sample size larger than the number of variables (error
  `efa_n_too_small`, as
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  already required); in
  [`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
  the affected criterion is listed under `not_run` instead of stopping
  the call.

- [`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md)
  now stops the grid and warns (`efa_map_truncated`, recording the
  reached point in `m_last`) when a residual variance hits zero, instead
  of continuing with stale, uncomputed values.

- [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md)
  now gives the error `efa_cd_degenerate_population` when `N_pop` is
  smaller than 2.

- Count arguments of the retention criteria (`n_factors_max`,
  `N_samples`, `max_iter` of
  [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md);
  `n_factors`, `N`, `n_vars`, `n_datasets` of
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md);
  `n_datasets` of
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md);
  `N`, `n_datasets` of
  [`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md))
  must now be at least 1.

- Factor counts in retention output no longer print in scientific
  notation (e.g. `3e+00`) under a negative `options(scipen)`.

### Input Validation and Correlation Handling

- [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  and
  [`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md)
  now reject `N = 0` (`N = NA` still means an unknown sample size).

- `efa_*` arguments with fixed choices are now case-insensitive and
  accept unambiguous abbreviations.

- Correlation matrices supplied as data frames are now recognised and
  analysed correctly.

- Passing `NULL` to a choice-valued argument now uses its documented
  default.

- Correlation-matrix smoothing now consistently returns a
  positive-definite matrix, including for nearly singular matrices.

- A square, non-symmetric matrix is now refused as a correlation matrix
  (error `efa_input_not_symmetric`) instead of being treated as raw
  data.

- A correlation matrix is now considered singular when its
  smallest-to-largest eigenvalue ratio falls below
  `n_vars * machine epsilon`, rather than only when
  [`solve()`](https://rdrr.io/r/base/solve.html) fails; some previously
  accepted, ill-conditioned matrices now give the error
  `efa_cor_singular`.

- [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  now gives the error `efa_n_factors_object` when a factor-retention
  object is supplied as `n_factors`.

### Missing Data and Multiple Imputation

- [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)’s
  `rmsr_upper` argument is deprecated and ignored (it never affected the
  reported RMSR);
  [`EFA_POOLED()`](https://mdsteiner.github.io/EFAtools/reference/EFA_POOLED.md)
  still accepts it silently.

- [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)
  now rejects `mids` objects from `mice`; convert them first with
  `mice::complete(x, "all")`.

- With `target_method = "consensus"`,
  [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)
  now starts the search for the common rotational target at the
  imputation closest to all others rather than the first; pooled results
  no longer depend on the order of `data_list` (they can differ from
  earlier versions). Pass `consensus_args = list(start = 1)` to recover
  the previous behaviour.

- [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)’s
  pooled RMSEA is now capped at 1, and its confidence bounds are
  withheld when they do not contain the point estimate.

- [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)’s
  pooled communalities are now named `communalities` (`MI$h2` is gone;
  use `MI$communalities`, `SE$communalities`, `CI$communalities`).

- [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)’s
  pooled `fit_indices` now hold the same elements in the same order
  across every standard-error route; the sandwich route gains
  `pool_method` (as `NA`).

- An
  [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)
  solution now carries an `mi_admissibility` component, reporting the
  Heywood flags of the individual imputations.

- A `cor_method = "fiml"` fit with analytic standard errors now
  withholds loading standard errors (warning `efa_se_unreliable`) when
  the rotational orientation is weakly determined, matching the other
  standard-error routes.

### Model Estimation

- The `"uls"` (minimum residual) estimator’s search now matches its
  analytic gradient and reported `Fm`, minimising the full
  reduced-correlation residual including the diagonal (previously the
  search used only the off-diagonal residual). Results move little for
  well-supported solutions; over-factored models can show larger
  changes.

- Squared multiple correlations that start the `"paf"`, `"ml"`, and
  `"uls"` estimators are now held within \[0, 1\], matching
  [`psych::smc()`](https://rdrr.io/pkg/psych/man/smc.html); results can
  change for indefinite or badly conditioned correlation matrices
  (e.g. unsmoothed bootstrap resamples).

- [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  gains `fiml_max_iter` and `fiml_tol` to govern the FIML two-stage EM
  algorithm.

- An
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  fit with `cor_method = "fiml"` now returns a `fiml` component
  (convergence state, iteration count, missing-data patterns, sample
  size), and warns (`efa_fiml_em_nonconvergence`) when the EM algorithm
  does not converge.

- A `cor_method = "fiml"` fit that cannot form the corrected two-stage
  chi-square now records `chi_scaled_type = "uncorrected.lrt"` and warns
  (`efa_fiml_uncorrected_chisq`).

- An RMSEA confidence interval that cannot be computed is now reported
  as `NA` instead of stopping the fit.

### Ordinal Correlations

- Perfectly ordered polychoric or tetrachoric pairs are now reported as
  `0.9999`/`-0.9999`, with a warning listing the affected pairs.

- Binary pairs with an empty response cell now use a margin-preserving
  0.5 continuity correction.

- Unavailable asymptotic variances for polychoric/tetrachoric pairs are
  now reported as `NA`; `DWLS` stops with a clear error, and affected
  robust standard errors are withheld.

- Rare response combinations in strongly correlated pairs are now
  handled more accurately, so such pairs less often block `DWLS`
  estimation or robust standard errors.

### Power Analysis

- A sample size solved by
  [`efa_power()`](https://mdsteiner.github.io/EFAtools/reference/efa_power.md)
  with `group > 1` is now a multiple of `group`, so `N_per_group` is a
  whole number (previously, e.g., a required total of 259 across two
  groups gave 129.5 per group).

- In simulation mode,
  [`efa_power()`](https://mdsteiner.github.io/EFAtools/reference/efa_power.md)
  now reports an `NA` hit rate when a requested factor-retention
  criterion never produces a suggestion.

- Simulation mode now gives clear errors for missing or invalid `N` and
  `n_datasets`.

- In simulation mode,
  [`efa_power()`](https://mdsteiner.github.io/EFAtools/reference/efa_power.md)
  now rejects a `p` that disagrees with the population model (previously
  replaced silently).

- In simulation mode,
  [`efa_power()`](https://mdsteiner.github.io/EFAtools/reference/efa_power.md)
  now records the failure reason in `replicates$fit_error` when fits
  that recover the model fail (previously `NA` with no explanation).

### Printed Output

- Truncated variable names in loading tables no longer collide: names
  that would collide are abbreviated and numbered so they stay
  distinguishable.

- [`print()`](https://rdrr.io/r/base/print.html)/[`format()`](https://rdrr.io/r/base/format.html)
  for Schmid-Leiman loading matrices now honour `max_name_length`,
  `name_style`, `sort_loadings`, and `max_factors_per_block` (previously
  accepted but ignored).

- Errors from argument checks now carry the condition class
  `efa_invalid_argument` and name the function called.

### Reliability and Factor Scores

- The whole-scale row of an
  [`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md)
  table for a correlated-factors solution is now labelled
  `factor = "total"`, `level = "total"` (previously `"g"`/`"general"`,
  which implied a general factor such solutions do not have). Solutions
  with an actual general factor, and single-factor solutions, are
  unchanged.

- [`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md)
  and
  [`OMEGA()`](https://mdsteiner.github.io/EFAtools/reference/OMEGA.md)
  now compute every omega total from the model-implied common variance,
  counting contributions from cross-loadings; subscale totals change for
  solutions with cross-loadings.

- With `variance = "correlation"`, the whole-scale omega total is now
  the model-implied common variance over observed total variance
  (previously total variance minus unique variances).

- With `variance = "sums_load"`, a subscale composite’s model-implied
  variance now also includes what it receives from other group factors,
  and this setting now applies to solutions without a general factor,
  including correlated-factors
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  solutions (previously silently overridden to `"correlation"`).

- A `Phi` supplied to
  [`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md)
  without a `pattern` is now treated as the group-factor correlation
  matrix of `s_load` and enters the coefficients.

- A `Phi` supplied together with a loading matrix of two or more factors
  is no longer dropped: the pair is now scored as a correlated-factors
  solution regardless of the matrix’s class (previously a hierarchy’s
  coefficients were returned).

- [`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md)
  now refuses a `Phi` supplied together with an
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  loading table (previously dropped silently), since such a table
  already states it is a hierarchy.

- [`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md)
  and
  [`OMEGA()`](https://mdsteiner.github.io/EFAtools/reference/OMEGA.md)
  now include a `lavaan` fit’s residual covariances in the model-implied
  composite variances; a freed residual covariance previously overstated
  omegas and understated standardized alpha.

- [`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md)
  and
  [`OMEGA()`](https://mdsteiner.github.io/EFAtools/reference/OMEGA.md)
  now reject a bifactor or second-order `lavaan` fit whose latent
  variables are correlated (previously scored as though uncorrelated).

- [`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md)
  now correctly scores an
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  loading table as a bifactor matrix, and no longer misreads an oblique
  solution’s pattern matrix as one (pass the
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  object itself for that solution).

- [`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md)
  now scores a `lavaan` fit whose variables each load on a single factor
  as the correlated-factors solution it is.

- [`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md)
  now returns omega total, standardized alpha, and the H index for a
  single-factor solution consistently across all input routes, including
  a one-factor
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  solution or a one-column loading matrix (previously refused).

- The one factor of a single-factor solution is now reported under its
  input name, or `"F1"` if unnamed (previously always `"g"`).

- A `cormat` supplied in a different variable order than the solution is
  now reordered to match it (previously gave wrong subscale
  coefficients).

- [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  and [`SL()`](https://mdsteiner.github.io/EFAtools/reference/SL.md) now
  reorder a supplied `Phi` to match the loading columns (a differently
  ordered named `Phi` previously gave a silently wrong solution) and
  error if `Phi` does not match.

- [`efa_scores()`](https://mdsteiner.github.io/EFAtools/reference/efa_scores.md)
  and
  [`FACTOR_SCORES()`](https://mdsteiner.github.io/EFAtools/reference/FACTOR_SCORES.md)
  now check and reorder a named `Phi` and scoring correlation matrix
  (`rho`) to match the model, instead of silently using a mismatched
  order.

- [`efa_scores()`](https://mdsteiner.github.io/EFAtools/reference/efa_scores.md)
  and
  [`FACTOR_SCORES()`](https://mdsteiner.github.io/EFAtools/reference/FACTOR_SCORES.md)
  now check a correlation matrix supplied in `x` against the model and
  align it to the model variables, as already done for raw data; a
  mismatched matrix is now an error (previously returned the fitted
  solution’s weights regardless of the matrix supplied).

- [`efa_scores()`](https://mdsteiner.github.io/EFAtools/reference/efa_scores.md)/[`FACTOR_SCORES()`](https://mdsteiner.github.io/EFAtools/reference/FACTOR_SCORES.md)
  with `method = "Bartlett"` or `"Anderson"` now stop on a solution with
  a communality at or above 1, instead of returning unusable weights
  with a warning.

- A `factor_map` (and
  [`OMEGA()`](https://mdsteiner.github.io/EFAtools/reference/OMEGA.md)’s
  `factor_corres`) must now hold only 0 and 1, and is checked against
  the loading matrix’s dimensions (previously any value was accepted and
  silently multiplied into the coefficients).

- [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  and [`SL()`](https://mdsteiner.github.io/EFAtools/reference/SL.md) now
  reject a solution with a single first-order factor.

- [`efa_scores()`](https://mdsteiner.github.io/EFAtools/reference/efa_scores.md)
  and
  [`FACTOR_SCORES()`](https://mdsteiner.github.io/EFAtools/reference/FACTOR_SCORES.md)
  now reject scoring data where a model variable is constant, infinite,
  or observed fewer than twice (previously produced `NaN` scores).

- [`efa_scores()`](https://mdsteiner.github.io/EFAtools/reference/efa_scores.md)
  now reports the count of scored cases in `settings$n_scored`.

- [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  and [`SL()`](https://mdsteiner.github.io/EFAtools/reference/SL.md) now
  reject the unsupported arguments `se`, `b_boot`, `ci`, and `seed`.

### Rotation

- Oblique rotations now refuse a nearly singular transformation matrix
  (smallest singular value below 0.0001) at every evaluation point,
  including the oblique Procrustes solver; results can differ for
  near-degenerate solutions (typically more factors than the data
  support). Well-conditioned solutions are unaffected.

- A gradient-projection rotation is now reported as converged only when
  the projected gradient meets tolerance (except simplimax, whose kinked
  objective still uses a stalled-progress criterion); convergence flags
  and diagnostics can differ from before for the same data.

- Rotation criterion parameters from
  [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  (`gam`, `delta`, `maxit`, simplimax `k`) are now validated before
  rotation runs, rejecting invalid values that previously reached the
  rotation engine silently or with an opaque error.

- [`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md)
  now refuses an `S` that is not `crossprod(A)` (previously a mismatched
  matrix silently minimised a different criterion), and refuses a badly
  conditioned `T_init`, not only a singular one.

- `varimax`/`promax` with `varimax_type = "svd"` and Kaiser
  normalisation no longer fail on a solution containing a
  zero-communality variable.

### Standard Errors

- [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  now rejects `b_boot` below 2 (error `efa_b_boot_too_small`) and `ci`
  of 0 or 1 (error `efa_ci_out_of_bounds`).

- `efa_fit(se = "sandwich")` now withholds standard errors and
  confidence intervals when a Heywood case occurs.

- Bootstrap output now reports the number of usable replicates when
  fewer than requested are available (for both single fits and pooled
  [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)
  fits).

- A bootstrap replicate whose rotation cannot be aligned now warns under
  the classed condition `efa_boot_rotation_failed`.

- Analytic standard-error output now includes factor names in `SE$Phi`
  and variable-factor labels in `vcov_unrot_loadings`.

## EFAtools 1.0.0

CRAN release: 2026-07-23

### New Interface

- The recommended interface now consists of lowercase `efa_*` functions,
  with estimation and rotation settings defined through
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  and
  [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md).

- The existing uppercase functions remain available without warnings, so
  current scripts continue to work, but new features will generally be
  added only to the `efa_*` interface.

- Estimators, rotations, retention criteria, scoring methods, and
  presets can now be specified without matching capitalization exactly,
  while results continue to store and display the standard spelling.

### Changes When Using the `efa_*` Interface

- The `efa_*` functions now report an error for unused, misspelled, or
  misplaced arguments instead of silently ignoring them.

- Estimation and rotation options such as convergence settings,
  iteration limits, and rotation parameters should now be supplied
  through
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  or
  [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md).

- Factor-retention functions always fit unrotated models and therefore
  reject rotation settings that would not affect the retention result.

- The estimator argument is now named `estimator` rather than `method`,
  and rotation arguments use the names `p_type` and `random_starts`
  rather than `P_type` and `randomStarts`.

- [`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
  requires the full argument name `max_iter_CD`, because the shorter
  form `max_iter` is no longer matched automatically.

- The second-order estimation preset for
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  is now specified with `estimate_control(type = ...)`.

- Control objects are validated when they are created, allowing invalid
  settings to be detected before a model is fitted.

### Reproducibility

- A supplied `seed` now controls the complete analysis, including random
  rotation starts, bootstrap samples, group-specific fits, and model
  averaging where applicable.

- [`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md)
  now has a `seed` argument and produces the same result regardless of
  the number of parallel workers.

- Seeded results from
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  and
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)
  no longer depend on the parallel-processing setup.

- Because random-number handling and several calculations were
  corrected, the same seed may produce different results than in earlier
  versions for
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md),
  [`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md),
  [`efa_group()`](https://mdsteiner.github.io/EFAtools/reference/efa_group.md),
  and some
  [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md)
  settings.

### Standard Errors and Fit Statistics

- `efa_fit(se = "information")` now returns `NA` with an
  `efa_se_unreliable` warning when standard errors for unrotated
  loadings cannot be estimated reliably, such as for Heywood cases or
  poorly distinguished factors.

- Information-based standard errors are now calculated for correlation
  matrices rather than covariance matrices, improving their accuracy and
  agreement with sandwich and bootstrap estimates.

- Information-based standard errors are now supported only with Pearson
  correlations and FIML; use `se = "sandwich"` for polychoric,
  tetrachoric, or rank correlations.

- FIML-based analytic standard errors now use the same small-sample
  scaling as the package’s other standard errors and fit statistics.

- [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)
  no longer reports AIC, BIC, or ECVI when these measures are not
  meaningful for the underlying fitted models.

### Multiple Imputation and Group Comparisons

- The pooled unrotated solution from
  [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)
  no longer depends on the order of the imputed datasets.

- Pooled unrotated loadings from
  [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)
  now use the same orientation conventions as a regular
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  result, making the two easier to compare.

- These changes can affect unrotated loadings, their standard errors,
  and statistics derived from them, but rotated solutions and their fit
  indices remain unchanged.

- Consensus results from
  [`efa_group()`](https://mdsteiner.github.io/EFAtools/reference/efa_group.md)
  no longer depend on the order in which groups are supplied, so
  congruences, loading differences, salience classifications, and
  invariance conclusions are now stable across group orderings.

### Factor Retention and Rotation

- Fixed a regression in the revised MAP criterion, which could
  previously suggest an incorrect number of factors.

- [`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md)
  now rejects polychoric and tetrachoric correlations because its tests
  and fit measures are not valid for these inputs.

- Rotated factors are now labelled `F1`, `F2`, and so forth according to
  the variance they explain, rather than inheriting labels that could
  change with the random seed.

- Rotation diagnostics now distinguish between the total number of
  starting solutions and the number that were fully optimized; code
  using `n_starts` should use `n_starts_total` or `n_optimized` instead.

- Criterion-based rotations now respect the specified `maxit` limit
  during every stage of the optimization.

### Model Averaging

- [`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md)
  no longer prints `NaN%` when no models converge and now clearly states
  which fitted models are included in each reported rate.

- Named presets such as `"SPSS"`, `"psych"`, and `"EFAtools"` now use
  their intended iteration limits, so models that exceed those limits
  are marked as non-converged and excluded from the average.

### Data Simulation

- Cudeck–Browne model-error simulations are now reproducible across
  computing platforms and numerical libraries, although a given seed
  will produce different populations than in earlier versions.

- With `marginals = "VM"`, categorized data now reproduce the requested
  category proportions more accurately, with a warning when the
  requested thresholds cannot be transformed safely.

- The documentation now clarifies that `match = "thresholds"` and
  `match = "polychoric"` produce the same data for normal marginals, and
  the selected value is now stored in the result settings.

- The documentation now explains the missing-data mechanism produced by
  `missing = "MAR"` and the residual bias it may create for analyses
  based only on the observed data.

- The documentation now explains that `model_error = "WB"` generally
  produces a larger RMSEA for the generating model than the requested
  target.

### Input Validation and Correlation Handling

- [`efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/efa_screen.md)
  now handles collinear or nearly collinear variables without failing
  during its outlier checks.

- Polychoric covariance calculations now handle sparse or nearly
  perfectly associated response tables more consistently, allowing DWLS
  analyses to continue while warning about problematic variable pairs.

- Covariance matrices supplied where raw data or a correlation matrix is
  expected are now rejected with a message suggesting
  [`cov2cor()`](https://rdrr.io/r/stats/cor.html).

- Polychoric and tetrachoric correlations now require ordered factors or
  numeric response codes, preventing character or unordered-factor
  levels from being interpreted in alphabetical order.

### Power Analysis and Factor Scores

- [`efa_power()`](https://mdsteiner.github.io/EFAtools/reference/efa_power.md)
  now clearly documents that `N` is the total sample size across groups
  and returns the corresponding sample size per group in `N_per_group`.

- Simulation-based power summaries now distinguish structure-recovery
  rates from the underlying Tucker congruence values.

- The documentation for
  [`FACTOR_SCORES()`](https://mdsteiner.github.io/EFAtools/reference/FACTOR_SCORES.md)
  now correctly describes `R2` as squared factor-score determinacy.

## EFAtools 0.8.0

CRAN release: 2026-07-07

### New Functions

- [`efa_group()`](https://mdsteiner.github.io/EFAtools/reference/efa_group.md)
  performs multi-group EFA.
- [`efa_power()`](https://mdsteiner.github.io/EFAtools/reference/efa_power.md)
  to conduct power analysis (analytic power based on RMSEA and
  simulation based power analysis).
- [`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md)
  to calculate various reliability indices.
- [`efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/efa_screen.md)
  to screen data for multivariate normality and suitability for factor
  analysis.
- [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md)
  to simulate data with various distributions and missing data
  mechanisms.

### Changes to Functions

- [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md):
  - gains `cor_method = "fiml"` for raw data with missing values.
  - gains `se = "sandwich"` and `se = "information"`.
  - renamed the top-level fields `boot.SE`, `boot.CI`, and `boot.arrays`
    on `EFA` objects to `SE`, `CI`, and `replicates`, as the old names
    are no longer accurate given the additional SE implementations.
  - All criterion-based rotations are now computed by a built-in
    gradient-projection rotation engine instead of the `GPArotation`
    package. The rotated solutions are numerically equivalent but
    implemented in C++ for speed to allow high numbers of random starts.
  - The default number of random starts for the rotation in
    [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md)
    (`randomStarts`) has been raised from 10 to 100, making local minima
    less likely for the rotation criteria that are prone to them.
  - additionally reports the Tucker-Lewis index (TLI, also called the
    non-normed fit index), the expected cross-validation index (ECVI),
    and the standardized root mean square residual (SRMR)
  - gains a `seed` argument and now fits the non-parametric bootstrap
    replicates (`se = "np-boot"`) in parallel across replicates via the
    `future` framework.
  - now has a [`summary()`](https://rdrr.io/r/base/summary.html) method
    that prints the full, detailed solution (loadings with communalities
    and uniquenesses, explained variances, fit indices, and residuals);
    [`print()`](https://rdrr.io/r/base/print.html) now gives a more
    compact overview.
- [`EFA_AVERAGE()`](https://mdsteiner.github.io/EFAtools/reference/EFA_AVERAGE-superseded.md)
  gains `cor_method = "fiml"` (passed to
  [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md)).
- [`EFA_POOLED()`](https://mdsteiner.github.io/EFAtools/reference/EFA_POOLED.md):
  - now dispatches its multiple-imputation standard-error pooling
    automatically by the standard-error method of its component fits:
    Rubin’s rules for the information method (with Wald confidence
    intervals), the two-stage pooled-inputs (MI2S) approach for the
    sandwich method, and the existing bootstrap pooling for the
    non-parametric bootstrap method.
  - Renamed the top-level field `boot.MI` to `MI`.
  - now defaults to `target_method = "first_target"`, which aligns every
    imputation to the first imputation’s rotated solution with a single
    Procrustes rotation. For orthogonal rotations this reaches the same
    pooled estimate as `"consensus"` with substantially less work.
    `"consensus"` (Generalized Procrustes Analysis of the
    imputation-specific rotated loadings) is still available as an
    opt-in for orthogonal rotations, and now raises a classed condition
    (`efa_consensus_oblique_unsupported`) when combined with an oblique
    rotation, because the iterated-oblique-Procrustes centroid can drift
    to degenerate targets.
- `cor_method` now accepts `"poly"` and `"tetra"` to compute polychoric
  and tetrachoric correlations from raw ordinal (respectively binary)
  data, using a two-step estimator with no empty-cell continuity
  correction. Supported by
  [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md)
  (including its non-parametric bootstrap),
  [`EFA_AVERAGE()`](https://mdsteiner.github.io/EFAtools/reference/EFA_AVERAGE-superseded.md),
  the suitability tests
  [`KMO()`](https://mdsteiner.github.io/EFAtools/reference/KMO.md) and
  [`BARTLETT()`](https://mdsteiner.github.io/EFAtools/reference/BARTLETT.md),
  and the retention criteria
  [`EKC()`](https://mdsteiner.github.io/EFAtools/reference/EKC.md),
  [`KGC()`](https://mdsteiner.github.io/EFAtools/reference/KGC.md),
  [`MAP()`](https://mdsteiner.github.io/EFAtools/reference/MAP.md),
  [`SCREE()`](https://mdsteiner.github.io/EFAtools/reference/SCREE.md),
  [`SMT()`](https://mdsteiner.github.io/EFAtools/reference/SMT.md), and
  [`N_FACTORS()`](https://mdsteiner.github.io/EFAtools/reference/N_FACTORS.md).
  The criteria that compare the data against simulated continuous
  reference data
  ([`CD()`](https://mdsteiner.github.io/EFAtools/reference/CD.md),
  [`PARALLEL()`](https://mdsteiner.github.io/EFAtools/reference/PARALLEL.md),
  [`NEST()`](https://mdsteiner.github.io/EFAtools/reference/NEST.md),
  and
  [`HULL()`](https://mdsteiner.github.io/EFAtools/reference/HULL.md)) do
  not support `"poly"` / `"tetra"` and signal an error.
- The factor-retention functions
  [`CD()`](https://mdsteiner.github.io/EFAtools/reference/CD.md),
  [`EKC()`](https://mdsteiner.github.io/EFAtools/reference/EKC.md),
  [`HULL()`](https://mdsteiner.github.io/EFAtools/reference/HULL.md),
  [`KGC()`](https://mdsteiner.github.io/EFAtools/reference/KGC.md),
  [`MAP()`](https://mdsteiner.github.io/EFAtools/reference/MAP.md),
  [`NEST()`](https://mdsteiner.github.io/EFAtools/reference/NEST.md),
  [`PARALLEL()`](https://mdsteiner.github.io/EFAtools/reference/PARALLEL.md),
  [`SCREE()`](https://mdsteiner.github.io/EFAtools/reference/SCREE.md),
  and [`SMT()`](https://mdsteiner.github.io/EFAtools/reference/SMT.md)
  now return objects of a common `efa_retention` class with shared
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.
- [`COMPARE()`](https://mdsteiner.github.io/EFAtools/reference/COMPARE.md)
  objects now have a
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method that
  returns a `ggplot` object. The plot is no longer drawn as a side
  effect of [`print()`](https://rdrr.io/r/base/print.html).
- Console output (the `print`, `summary`, and `format` methods) is now
  produced with the `cli` package, and the messages, warnings, and
  errors emitted across the package are now S3-classed conditions, which
  makes them easier to handle programmatically.
- [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md),
  [`SL()`](https://mdsteiner.github.io/EFAtools/reference/SL.md), and
  [`EFA_AVERAGE()`](https://mdsteiner.github.io/EFAtools/reference/EFA_AVERAGE-superseded.md)
  now accept `method = "MINRES"` as a synonym for `method = "ULS"`.
  Minimum residual and unweighted least squares are two names for the
  same estimator and return identical results.

### Bug Fixes

- [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md) and
  [`EFA_POOLED()`](https://mdsteiner.github.io/EFAtools/reference/EFA_POOLED.md):
  The comparative fit index (CFI) now floors the model and baseline
  noncentralities at zero before taking their ratio (Bentler, 1990), so
  it is no longer deflated toward zero when the baseline (independence)
  model fits comparatively well. The value is unchanged for well-fitting
  models and now matches `lavaan`. CFI can change for solutions in which
  the baseline model does not misfit much.
- [`EFA_POOLED()`](https://mdsteiner.github.io/EFAtools/reference/EFA_POOLED.md):
  Corrected and extended the pooled chi-square-based model-fit indices.
  The D2 average relative increase in variance is now centred on the
  mean of the square-root statistics (Li, Meng, Raghunathan & Rubin,
  1991), removing a one-sided inflation of the pooled chi-square, RMSEA,
  and CFI. Bootstrap/MI confidence intervals for the pooled fit indices
  are now the Rubin-Wald multiple-imputation summaries; a miscalibrated
  bootstrap-percentile interval (obtained by re-running the pooling
  algorithm over matched replicate indices) was removed, as was a
  mislabeled pooled `Fm`.
- [`EFA_AVERAGE()`](https://mdsteiner.github.io/EFAtools/reference/EFA_AVERAGE-superseded.md):
  When every averaged solution fails (all runs error, fail to converge,
  or are Heywood cases), the function now returns an empty (`NA`)
  averaged result instead of erroring or averaging an empty set.
- [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md) with
  `method = "PAF"`: the reported number of iterations (`iter`) is now
  the number of iterations actually performed; it was previously one too
  high. The loadings, communalities, and convergence status are
  unchanged.
- `ULS` extraction: the minimised objective is now the sum of squared
  off-diagonal residuals, consistent with its analytic gradient and the
  reported minimum (`Fm`). The diagonal residuals were previously
  included in the objective but not in its gradient. The fitted loadings
  are unchanged to within optimiser tolerance.
- [`NEST()`](https://mdsteiner.github.io/EFAtools/reference/NEST.md) and
  [`PARALLEL()`](https://mdsteiner.github.io/EFAtools/reference/PARALLEL.md):
  A failed eigendecomposition of a degenerate simulated matrix now
  raises a clear error instead of resulting in undefined behaviour.
- The chi-square model-fit statistic is now the Bartlett-corrected
  maximum likelihood discrepancy evaluated at the model-implied
  correlation matrix, for both `ML` and `ULS` extraction. For
  `method = "ML"` this matches
  [`stats::factanal()`](https://rdrr.io/r/stats/factanal.html) and
  [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html); the
  small-sample Bartlett correction was previously omitted. For
  `method = "ULS"` it is now a proper chi-square-distributed statistic
  matching `psych::fa(fm = "minres")`; previously the least-squares
  residual sum of squares was multiplied by `N - 1` and read as if it
  were Wishart-distributed, which produced an invalid p-value. The
  independence (baseline) model used for the CFI is rescaled onto the
  same discrepancy scale. Consequently the p-value, CFI, RMSEA (and its
  confidence interval), AIC, and BIC change for ML and ULS solutions,
  and the number of factors suggested by
  [`SMT()`](https://mdsteiner.github.io/EFAtools/reference/SMT.md) and
  [`HULL()`](https://mdsteiner.github.io/EFAtools/reference/HULL.md) may
  change for these methods.
- [`PARALLEL()`](https://mdsteiner.github.io/EFAtools/reference/PARALLEL.md):
  The percentile reference eigenvalues are now computed with
  [`stats::quantile()`](https://rdrr.io/r/stats/quantile.html),
  correcting a slight off-by-one in the previous indexing.
- [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md): For
  oblique rotations, the factor intercorrelations (`Phi`), the structure
  matrix, and the explained variances are now reflected and reordered
  consistently with the rotated loadings. Previously, when a factor was
  reflected to a positive orientation the factor intercorrelations were
  not sign-adjusted (so the structure matrix and reported correlations
  did not match the loadings), and with `order_type = "ss_factors"` the
  factor intercorrelations were not reordered at all.
- [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md): The
  returned rotation matrix (`rotmat`) is now reflected and reordered
  consistently with the rotated loadings, for both orthogonal and
  oblique rotations, so that the rotated loadings can be reconstructed
  from the unrotated loadings and `rotmat`. Previously the sign
  reflection was not applied to it and its factors were left in a
  different order, so this reconstruction did not hold whenever a factor
  was reflected or reordered.
- [`HULL()`](https://mdsteiner.github.io/EFAtools/reference/HULL.md):
  The convex-hull elimination now tests every triplet of adjacent
  solutions, including the one formed by the last interior solution.
  Previously this final triplet was skipped, so a solution lying below
  the line connecting its neighbours could remain on the hull. This can
  change the number of factors suggested by
  [`HULL()`](https://mdsteiner.github.io/EFAtools/reference/HULL.md)
  (and hence by
  [`N_FACTORS()`](https://mdsteiner.github.io/EFAtools/reference/N_FACTORS.md)).
- [`NEST()`](https://mdsteiner.github.io/EFAtools/reference/NEST.md):
  When the test accepts every eigenvalue it examines (no empirical
  eigenvalue falls at or below its reference within the search range),
  it now retains that last accepted model rather than the model with one
  fewer factor. The search range is also bounded so that the reference
  model fitted at each step stays over-identified.
- [`PARALLEL()`](https://mdsteiner.github.io/EFAtools/reference/PARALLEL.md):
  When every real eigenvalue exceeds its reference, so the decision rule
  finds no crossing, the suggested number of factors is now all retained
  components, reported with a warning, instead of a silent `NA`
  (matching the convention used by
  [`EKC()`](https://mdsteiner.github.io/EFAtools/reference/EKC.md)).
- `EFA(se = "np-boot")`: The non-parametric bootstrap no longer repeats
  per-replicate warnings (about arguments pinned alongside `type`, and
  about the iterative fit reaching its maximum number of iterations)
  once per replicate. They are now suppressed during resampling, and
  non-convergence across replicates is reported once as a summary.
- [`COMPARE()`](https://mdsteiner.github.io/EFAtools/reference/COMPARE.md):
  With `reorder = "congruence"`, the columns of `y` are now matched to
  those of `x` by an optimal one-to-one assignment that maximizes the
  total Tucker congruence, rather than by an independent per-factor best
  match. The greedy match could assign two factors of `x` to the same
  column of `y`, duplicating that column and dropping another, which
  corrupted the reported differences, factor correspondences, and root
  mean squared distance. The result is unchanged whenever the greedy
  match was already one-to-one.

## EFAtools 0.7.1

CRAN release: 2026-05-08

### Changes to Functions

- [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md):
  Added `randomStarts` argument passed to GPArotation functions, as
  suggested by Coen Bernaards.
- [`FACTOR_SCORES()`](https://mdsteiner.github.io/EFAtools/reference/FACTOR_SCORES.md):
  Added `rho` argument, thanks to Andreas Soteriades.
- [`EFA_POOLED()`](https://mdsteiner.github.io/EFAtools/reference/EFA_POOLED.md):
  Fixed issue that could lead to averaged Phi not being symmetric.

## EFAtools 0.7.0

CRAN release: 2026-04-30

### New Functions

- [`MAP()`](https://mdsteiner.github.io/EFAtools/reference/MAP.md)
  computes the Velicer MAP criterion (both TR2 and TR4).
- [`PROCRUSTES()`](https://mdsteiner.github.io/EFAtools/reference/PROCRUSTES.md)
  to perform orthogonal and oblique Procrustes / target rotation.
- `CONSENSUS_PROCRUSTES()` to perform Procrustes on a list of targets to
  obtain a common target.
- `residuals.EFA()` to extract and print residuals and, if computed,
  standardized residuals.
- [`EFA_POOLED()`](https://mdsteiner.github.io/EFAtools/reference/EFA_POOLED.md)
  to run EFA on a list of multiple imputated datasets and pool the
  results. Thanks to Andreas Soteriades for the suggestion and first
  implementation.
- `print.EFA_POOLED()`, print method adapted from `print.EFA()`

### Changes to Functions

- [`N_FACTORS()`](https://mdsteiner.github.io/EFAtools/reference/N_FACTORS.md)
  can also compute the MAP criterion.
- [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md):
  - Now returns and prints residuals and, if SEs are computed,
    standardized residuals.
  - Calculates and prints RMSR.
  - Can calculate bootstrap standard errors and CIs of parameters and
    fit indices.
- `print.EFA()` now prints more information.

### Bug Fixes

- Fixed a bug in
  [`OMEGA()`](https://mdsteiner.github.io/EFAtools/reference/OMEGA.md)
  that led to incorrect omega, H, and ECV values for `lavaan` bifactor
  models. Tnanks to Christopher D. King for bug report and suggested
  fix.
- Small fix in the documentation of
  [`EFA_AVERAGE()`](https://mdsteiner.github.io/EFAtools/reference/EFA_AVERAGE-superseded.md).
- Fixed incorrect calculation of RMSEA in
  [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md).

## EFAtools 0.6.1

CRAN release: 2025-07-30

### Changes to Functions

- [`EKC()`](https://mdsteiner.github.io/EFAtools/reference/EKC.md): Now
  correctly returns the number of factors based on the first time the
  eigenvalues drop below the references values.

## EFAtools 0.6.0

CRAN release: 2025-06-19

### New Functions

- Added
  [`NEST()`](https://mdsteiner.github.io/EFAtools/reference/NEST.md) to
  perform the Next Eigenvalue Sufficiency Test (Achim, 2017).

### Changes to Functions

- [`EKC()`](https://mdsteiner.github.io/EFAtools/reference/EKC.md): The
  implementation based on Auerswald and Moshagen (2019) used in previous
  versions differed from the original implementation by Braeken and van
  Assen (2017). Now both versions are implemented and can be selected
  with the new `type` argument. Thanks to Luis Eduardo Garrido for
  pointing this out and to Johan Braeken for sharing sample code, based
  on which the original version is now implemented.
- [`N_FACTORS()`](https://mdsteiner.github.io/EFAtools/reference/N_FACTORS.md):
  - Updated default settings to only use often used and well performing
    factor retention methods (others, like the Kaiser Guttman criterion
    can still be used).
  - Added NEST as additional factor retention method.
  - New arguments: `ekc_type`, `alpha_nest`, and `n_datasets_nest` to
    account for the changes in
    [`EKC()`](https://mdsteiner.github.io/EFAtools/reference/EKC.md) and
    to control
    [`NEST()`](https://mdsteiner.github.io/EFAtools/reference/NEST.md).

## EFAtools 0.5.0

CRAN release: 2025-05-23

### Changes to Functions

- `print.EFA()`:
  - Now returns explained variance from rotated, rather than unrotated
    solution, if a rotation was performed.
  - Now prints communalities and uniquenesses in loading/pattern matrix
- [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md):
  Calculate and return model implied correlation matrix.

## EFAtools 0.4.6

CRAN release: 2025-03-21

### General

- Small change in test of `.gof()`: Changed some tests to take care of
  the ATLAS issue when using R-devel on x86_64 Fedora 34 Linux with
  alternative BLAS/LAPACK.

## EFAtools 0.4.5

CRAN release: 2024-12-22

### Bug Fixes

- Updated
  [`OMEGA()`](https://mdsteiner.github.io/EFAtools/reference/OMEGA.md)
  to accommodate changes in the upcoming version of
  [`psych::schmid()`](https://rdrr.io/pkg/psych/man/schmid.html)

## EFAtools 0.4.4

CRAN release: 2023-01-06

### Changes to Functions

- `print.EFA()`: Added arguments `cutoff`, `digits` and
  `max_name_length` that are passed to `print.LOADINGS()`.
- `print.LOADINGS()`: New Argument `max_name_length` to control the
  maximum length of the displayed variable names (names longer than this
  will be cut on the right side). Previously, this was fixed to 10
  (which is now the default).

### Misc

- Updated a test of a helper function (`.gof()`) that threw an error
  when using R-devel on x86_64 Fedora 36 Linux with alternative
  BLAS/LAPACK.
- added `dontrun` to examples of
  [`EFA_AVERAGE()`](https://mdsteiner.github.io/EFAtools/reference/EFA_AVERAGE-superseded.md)
  and its print and plot methods as these were causing issues on
  R-oldrel which were not directly related to EFAtools and thus could
  not be fixed from within the package.

## EFAtools 0.4.3

CRAN release: 2022-10-02

### Changes to Functions

- `.gof()`: Changed the helper function to take care of the MKL issue
  when using R-devel on x86_64 Fedora 34 Linux with alternative
  BLAS/LAPACK.

## EFAtools 0.4.2

CRAN release: 2022-09-27

### Changes to Functions

- `.is_cormat()`: Changed the helper function to better detect wheter a
  matrix is a correlation matrix.
- [`PARALLEL()`](https://mdsteiner.github.io/EFAtools/reference/PARALLEL.md):
  Added a check, testing whether N \> n_vars and throw an error if this
  is not the case.

### Bug Fixes

- Fixed some tests due to upcoming changes in the psych package which
  EFAtools depends on.

## EFAtools 0.4.1

CRAN release: 2022-04-24

### Bug Fixes

- Minor fixes in tests to solve problems on macOS m1.

## EFAtools 0.4.0

CRAN release: 2022-03-21

### Changes to Functions

- [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md):
  Changed error to warning when model is underidentified. This allows
  the Schmid-Leiman transformation to be performed on a two-factor
  solution.
- [`OMEGA()`](https://mdsteiner.github.io/EFAtools/reference/OMEGA.md):
  Added calculation of additional indices of interpretive relevance (H
  index, explained common variance \[ECV\], and percent of
  uncontaminated correlations \[PUC\]). This is optional and can be
  avoided by setting `add_ind = FALSE`.

### Bug Fixes

- [`CD()`](https://mdsteiner.github.io/EFAtools/reference/CD.md): Added
  [`na.omit()`](https://rdrr.io/r/stats/na.fail.html) to remove missing
  values from raw data to avoid an error in the comparison-data
  procedure.

## EFAtools 0.3.1

CRAN release: 2021-03-27

### General

- When testing for whether a matrix is singular and thus smoothing
  should be done, test against .Machine\$double.eps^.6 instead of 0, as
  suggested by Florian Scharf.

### Changes to Functions

- [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md):
  - Added warnings if `type = "SPSS"` was used with `method = "ML"` or
    `method = "ULS"`, or with a rotation other than `none`, `varimax` or
    `promax`.
  - Avoided smoothing of non-positive definite correlation matrices if
    `type = "SPSS"` is used.
  - Use Moore-Penrose Pseudo Inverse in computation of SMCs if
    `type = "psych"` is used, by calling
    [`psych::smc()`](https://rdrr.io/pkg/psych/man/smc.html).
  - Use `varimax_type = "kaiser"` if `type = "EFAtools"` is used with
    `varimax` or `promax`.

### Bug Fixes

- [`EFA_AVERAGE()`](https://mdsteiner.github.io/EFAtools/reference/EFA_AVERAGE-superseded.md):
  - Added `future.seed = TRUE` to call to
    [`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html)
    to prevent warnings.
  - Fixed test for Heywood cases from testing whether a communality or
    loading is greater than .998, to only test whether communalities
    exceed 1 + .Machine\$double.eps
- `print.EFA()`: Fixed test for Heywood cases from testing whether a
  communality or loading is greater than .998, to only test whether
  communalities exceed 1 + .Machine\$double.eps
- [`OMEGA()`](https://mdsteiner.github.io/EFAtools/reference/OMEGA.md):
  Small bugfix when `lavaan` second-order model is given as input

## EFAtools 0.3.0

CRAN release: 2020-11-04

### General

- Added examples for
  [`EFA_AVERAGE()`](https://mdsteiner.github.io/EFAtools/reference/EFA_AVERAGE-superseded.md)
  to readme and the EFAtools vignette
- Updated examples in readme and vignettes according to the updated
  `OMEGA` function

### New Functions

- Added function
  [`EFA_AVERAGE()`](https://mdsteiner.github.io/EFAtools/reference/EFA_AVERAGE-superseded.md)
  and respective print and plot methods, to allow running many EFAs
  across different implementations to obtain an average solution and
  test the stability of the results.

### Changes to Functions

- [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md):
  Defaults that were previously set to `NULL` are now mostly set to
  `NA`. This was necessary for
  [`EFA_AVERAGE()`](https://mdsteiner.github.io/EFAtools/reference/EFA_AVERAGE-superseded.md)
  to work correctly.
- [`PARALLEL()`](https://mdsteiner.github.io/EFAtools/reference/PARALLEL.md):
  Rewrote the generation of random data based eigenvalues to be more
  stable when SMCs are used.
- [`OMEGA()`](https://mdsteiner.github.io/EFAtools/reference/OMEGA.md):
  Changed expected input for argument `factor_corres` from vector to
  matrix. Can now be a logical matrix or a numeric matrix with 0’s and
  1’s of the same dimensions as the matrix of group factor loadings.
  This is more flexible and allows for cross-loadings.

## EFAtools 0.2.0

CRAN release: 2020-09-17

### General

- Created new vignette *Replicate_SPSS_psych* to show replication of
  original `psych` and `SPSS` EFA solutions with `EFAtools`.

### New Functions

- Added function
  [`FACTOR_SCORES()`](https://mdsteiner.github.io/EFAtools/reference/FACTOR_SCORES.md)
  to calculate factor scores from a solution from
  [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md). This
  is just a wrapper for the
  [`psych::factor.scores`](https://rdrr.io/pkg/psych/man/factor.scores.html)
  function.
- Added function
  [`SCREE()`](https://mdsteiner.github.io/EFAtools/reference/SCREE.md)
  that does a scree plot. Also added respective print and plot methods.

### Changes to Functions

- [`CD()`](https://mdsteiner.github.io/EFAtools/reference/CD.md): Added
  check for whether entered data is a tibble, and if so, convert to
  vanilla data.frame to avoid breaking the procedure.
- [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md):
  - Updated the EFAtools type in PAF and Promax.
  - Added p value for chi square value in output (calculated for ML and
    ULS fitting methods).
  - Updated the SPSS varimax implementation to fit SPSS results more
    closely.
  - Created an argument “varimax_type” that is set according to the
    specified type, but that can also be specified individually. With
    type R psych and EFAtools, the stats::varimax is called by default
    (`varimax_type = "svd"`), with type SPSS, the reproduced SPSS
    varimax implementation is used (`varimax_type = "kaiser"`).
  - Renamed the `kaiser` argument (controls if a Kaiser normalization is
    done or not) into `normalize` to avoid confusion with the
    `varimax_type` argument specifications.
- `ML()`: Changed default start method to “psych”.
- [`N_FACTORS()`](https://mdsteiner.github.io/EFAtools/reference/N_FACTORS.md):
  - Added option to do a scree plot if “SCREE” is included in the
    `criteria` argument.
  - Added a progress bar.
- [`OMEGA()`](https://mdsteiner.github.io/EFAtools/reference/OMEGA.md):
  Now also works with a lavaan second-order solution as input. In this
  case, it does a Schmid-Leiman transformation based on the first- and
  second-order loadings first and computes omegas based on this
  Schmid-Leiman solution.
- [`SL()`](https://mdsteiner.github.io/EFAtools/reference/SL.md): Now
  also works with a lavaan second-order solution as input (first- and
  second-order loadings taken directly from lavaan output).

### Bug Fixes

- `.get_compare_matrix()`: Fixed a bug that occurred when names of data
  were longer than n_char
- [`COMPARE()`](https://mdsteiner.github.io/EFAtools/reference/COMPARE.md):
  Fixed a bug that occurred when using `reorder = "names"`.
- [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md):
  RMSEA is now set to 1 if it is \> 1.
- [`HULL()`](https://mdsteiner.github.io/EFAtools/reference/HULL.md):
  Fixed a bug that occurred when no factors are located on the HULL
- [`KMO()`](https://mdsteiner.github.io/EFAtools/reference/KMO.md):
  Fixed a bug that the inverse of the correlation matrix was not taken
  anew after smoothing was necessary.
- [`PARALLEL()`](https://mdsteiner.github.io/EFAtools/reference/PARALLEL.md):
  - Fixed a bug that occurred when using `decision_rule = "percentile"`
  - Relocated error messages that were not evaluated if no data were
    entered (and should be)
- `print.COMPARE()`: Fixed a bug that occurred when using
  `print_diff = FALSE` in
  [`COMPARE()`](https://mdsteiner.github.io/EFAtools/reference/COMPARE.md).
- `print.KMO()`: Fixed a bug that printed two statements instead of one,
  when the KMO value was \< .6.

### Minor Changes

- [`OMEGA()`](https://mdsteiner.github.io/EFAtools/reference/OMEGA.md)
  and [`SL()`](https://mdsteiner.github.io/EFAtools/reference/SL.md):
  Added an error message if the entered term in `g_name` is invalid
  (i.e., it cannnot be found among the factor names of the entered
  lavaan solution).

## EFAtools 0.1.1

CRAN release: 2020-07-13

### Minor Changes

- Added an error message in
  [`PARALLEL()`](https://mdsteiner.github.io/EFAtools/reference/PARALLEL.md)
  if no solution has been found after 25 tries.

### Bug Fixes

- Updated different tests

- Deleted no longer used packages from Imports and Suggests in
  DESCRIPTION

- [`PARALLEL()`](https://mdsteiner.github.io/EFAtools/reference/PARALLEL.md):
  fixed a bug in indexing if method `"EFA"` was used.

## EFAtools 0.1.0

CRAN release: 2020-07-07

- Added a `NEWS.md` file to track changes to the package.
- Initial CRAN submission
