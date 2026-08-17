# EFAtools 1.0.0.9000

## Comparing and Averaging Solutions

* `efa_compare()` display settings can now be changed through the corresponding `print()`, `format()`, and `plot()` methods without rerunning the comparison. Otherwise, the saved settings are used.

* In `efa_compare()`, `are_equal` is now `NA` (printed as "none") when the integer parts differ. A value of `0` is reserved for values whose integer parts agree but decimal places do not.

* The `efa_compare()` report now separates summary statistics from elementwise differences. The minimum number of supplied decimal places is shown only for rounded inputs.

* The Model Fit section of `efa_average()` now shows how many solutions contributed to each index. Counts may differ because PAF contributes residual-based, but not chi-square-based, fit indices.

* `efa_average()` now clarifies that fit indices summarise the individual fitted solutions, not the averaged loading matrix. Range tables are labelled `max - min`.

* `efa_average()` now points out that `trim` is ignored when `averaging = "median"`, instead of recording a trim that never affected a value.

* The `efa_average()` progress bar now uses the correct number of steps and reaches 100%, avoiding spurious warnings.

* Post-fitting stages in `efa_average()` now use progress output that avoids stray blank lines in rendered documents and check logs.

* The `efa_group()` report now notes when all factors are rated "equal" despite substantial loading differences. Tucker's congruence does not detect proportional differences in loading magnitude.

* Wide `efa_group()` output is now wrapped or split into blocks to fit the console.

* The elementwise differences of `efa_compare()` are signed, and the table is now headed by the direction of the subtraction, named with `x_labels`. The line counting differing indicator-to-factor correspondences wraps to the console width.

* `efa_average()` now documents that every solution in the grid is fitted with the eigenvalue-based factor ordering, including under a named `type`, and what the settings of the single `efa_fit()` object record when the grid collapses to one combination.

## Data Screening

* The multivariate-outlier check in `efa_screen()` now correctly handles robust subsets with singular covariance matrices instead of discarding the best solution.

* Robust covariance estimation in `efa_screen()` no longer depends on the variables' measurement units, preventing unnecessary fallback to classical Mahalanobis distances when scales differ greatly.

## Data Simulation

* `efa_simulate()` gains `missing_vars` to restrict missingness to selected variables. With `missing = "MAR"`, a `missing_predictor` outside this set remains complete. If omitted, all variables can contain missing values as before.

* With `marginals = "VM"` and `force_pd = TRUE`, the returned `population` now reflects the projected correlations used to generate the data, and a warning reports their deviation from the target. `return_pop = TRUE` still returns the target because no data are generated.

* With `model_error = "TKL"`, `efa_simulate()` now warns when the requested RMSEA and CFI cannot both be achieved and reports the attained values.

* With empirical marginals, `efa_simulate()` now gives a clear error when a category is too rare to be reproduced at the requested sample size.

## Factor Retention

* `efa_cd()` now reports a clear error when the data contain constant variables.

* `efa_smt()` now warns when a recommended model has a Heywood case or has not converged. The recommendation is still returned but should be treated cautiously.

* The RMSEA rule in `efa_smt()` now stops when the lower bound of the RMSEA confidence interval cannot be computed instead of continuing to later models.

* `efa_nest()` results can now be plotted, showing the empirical and reference eigenvalues and the retained number of factors.

* `efa_retain()` now prints a one-line overview of the number and range of numeric recommendations, the number of contributing criteria, and the most common factor count.

* Factor-retention functions now reject the unsupported arguments `seed`, `se`, `b_boot`, and `ci`. Use `set.seed()` to control randomness.

* `efa_parallel()` now reads the `"PCA"` and `"SMC"` reference eigenvalues off the same simulated datasets when both are requested, halving the simulation and pairing the two reference series dataset by dataset. A draw that the `"SMC"` series cannot use is discarded for `"PCA"` as well.

* A failed block of the `efa_parallel()` simulation is now redrawn on its own instead of restarting the whole batch, and the abort that follows repeated failures reports the failure that caused it.

## Input Validation and Correlation Handling

* `efa_fit()` and `efa_average()` now reject `N = 0`. `N = NA` remains the way to say that the sample size is unknown.

* `efa_*` arguments with fixed choices are now case-insensitive. Those accepting one choice also allow unambiguous abbreviations, invalid choices produce clearer errors, and results use the standard spelling.

* Correlation matrices supplied as data frames are now recognised and analysed correctly.

* Passing `NULL` to a choice-valued argument now uses its documented default.

* Errors when computing correlations from raw data now identify non-numeric, infinite, and constant columns, and recommend `cor_method = "poly"` for ordinal factors or character variables.

* Correlation-matrix smoothing now consistently returns a positive-definite matrix, including for nearly singular matrices.

## Missing Data and Multiple Imputation

* The `rmsr_upper` argument of `efa_mi()` is deprecated and ignored because it never changed the reported RMSR. `EFA_POOLED()` continues to accept it without a warning.

* Printed `efa_mi()` output now identifies CFI and TLI as averages across imputations rather than indices calculated from separately pooled model and baseline statistics.

* The `efa_mi()` documentation now describes the `mi_diagnostics`, `fits`, `alignment`, and `settings` components, including how to calculate a `lavaan.mi`-style reference CFI from the pooled chi-square statistics.

* `efa_mi()` now rejects `mids` objects from `mice` with guidance to convert them first using `mice::complete(x, "all")`.

* The pooled header and the pooling settings printed by `efa_mi()` and its `summary()` now wrap to the console width instead of running on to a fixed length. Individual `setting = 'value'` entries are still never split across a line.

## Ordinal Correlations

* Perfectly ordered polychoric or tetrachoric pairs are now reported as `0.9999` or `-0.9999`, with a warning listing the affected pairs.

* Binary pairs with one empty response cell now use a margin-preserving 0.5 continuity correction, avoiding boundary estimates and making a positive-definite correlation matrix more likely.

* For both kinds of pair, unavailable asymptotic variances are reported as `NA`; `DWLS` now stops with a clear error, and affected robust standard errors are withheld.

* Rare response combinations in strongly correlated pairs are now handled more accurately, preventing such pairs from unnecessarily blocking `DWLS` estimation or robust standard errors.

## Power Analysis

* In simulation mode, `efa_power()` now warns when a requested factor-retention criterion never produces a suggestion and reports an `NA` hit rate over zero valid replicates.

* Simulation mode now gives clear errors for missing or invalid `N` and `n_datasets`.

## Printed Output

* Truncating long variable names no longer gives two items the same row label in a loading table. Names that would collide are abbreviated instead, and numbered if abbreviating is not enough, so items sharing a long prefix stay distinguishable.

* `print()` and `format()` for Schmid-Leiman loading matrices now honour `max_name_length`, `name_style`, `sort_loadings`, and `max_factors_per_block`, which were previously accepted and ignored. Row sorting compares the group factors only, since the general factor is the largest loading of nearly every item.

* A loading table split into stacked column blocks now shows a blank line between them, so a block's header no longer sits directly under the last row of the block before it.

* The model header of `efa_fit()` and the second-order header of `efa_schmid_leiman()` now wrap to the console width. As in every other report line of this kind, an individual `setting = 'value'` entry is never split across a line.

## Reliability and Factor Scores

* `efa_reliability()` and `OMEGA()` now compute every omega total from the model-implied common variance of the composite, counting what it receives from every factor its variables load on, so subscale totals change for any solution with cross-loadings.

* With `variance = "correlation"`, the whole-scale omega total of `efa_reliability()` and `OMEGA()` is now the model-implied common variance of all variables over their observed total variance, rather than that total variance less the unique variances.

* With `variance = "sums_load"`, a subscale composite's model-implied variance now also includes what the composite receives from the other group factors.

* `variance = "sums_load"` now applies to a solution without a general factor, including a correlated-factors `efa_fit()` one, instead of being overridden to `"correlation"` with a warning.

* A `Phi` supplied to `efa_reliability()` without a `pattern` is now the correlation matrix of the group factors of `s_load` and enters the coefficients, so manually supplied components can describe a correlated-factors solution.

* `efa_reliability()` and `OMEGA()` now warn when the variables of a solution are not all keyed in the same direction.

* `efa_reliability()` now warns when supplied uniquenesses do not complete the loadings to unit variance.

* The out-of-range warning of `efa_reliability()` and `OMEGA()` now also covers `variance = "sums_load"`.

* `efa_reliability()` now refuses a `Phi` supplied with a matrix it reads as a hierarchy -- a bifactor loading matrix, or the loading table of an `efa_schmid_leiman()` solution -- instead of dropping it and returning that hierarchy's coefficients. Such a matrix has its general factor in the first column and uncorrelated group factors, so there is nothing for `Phi` to describe, and the same combination supplied through the components was already an error. The pattern matrix of an oblique solution, which carries the loading class `efa_fit()` gives it, still reads `Phi` as its factor intercorrelations.

* `efa_reliability()` and `OMEGA()` now count the residual covariances of a `lavaan` fit in the model-implied composite variances. A fit with a freed covariance between two variables previously divided by the residual variances alone, which overstated the omegas and understated the standardized alpha of every composite holding such a pair. Fits whose residual covariance matrix is diagonal are unchanged.

* `efa_reliability()` and `OMEGA()` now reject a bifactor or second-order `lavaan` fit whose latent variables are correlated instead of scoring it as though they were not.

* `efa_reliability()` now scores the loading table of an `efa_schmid_leiman()` solution as the bifactor matrix it is, and no longer reads the pattern matrix of an oblique solution as a bifactor one; pass the `efa_fit()` object itself to score that solution.

* `efa_reliability()` now scores a `lavaan` fit whose variables each load on a single factor as the correlated-factors solution it is, instead of asking for a general factor it does not have.

* `efa_reliability()` now returns omega total, standardized alpha, and the H index for a single-factor solution on every input route, instead of a different subset of them depending on how the solution was supplied.

* `efa_reliability()` now scores a one-factor `efa_fit()` solution, which has nothing to rotate, instead of refusing it for want of an oblique rotation.

* `efa_reliability()` now scores a one-column loading matrix as the single factor it is, instead of refusing it and asking for a hierarchy or for factor intercorrelations that a solution with one factor does not have. Loading matrices with two or more factors are still refused.

* `efa_reliability()` now rejects a `fac_names` that does not have one name per factor, instead of failing with an opaque error from `rownames()` or, for a single-factor solution, silently falling back to the default label.

* The one factor of a single-factor solution is now reported in `efa_reliability()` under the name its input gives it, or as `"F1"` where the input names it none, rather than always as `"g"`.

* `efa_reliability()` now scores the pattern matrix of an oblique solution supplied together with its factor intercorrelations in `Phi` as a correlated-factors solution, and still rejects a matrix of two or more factors when `Phi` is absent.

* A `cormat` supplied to `efa_reliability()` in another variable order than the solution is now reordered to it instead of yielding wrong subscale coefficients.

* `efa_schmid_leiman()` and `SL()` no longer silence the second-order factor analysis, so its Heywood and non-convergence warnings are now shown.

* `efa_schmid_leiman()` and `SL()` now reject a solution with a single first-order factor, which has nothing to orthogonalise, instead of failing inside the second-order factor analysis with an error about the intercorrelation matrix.

* `efa_reliability()` now warns when a `factor_map` assigns items to a group factor on which they load weakly while loading on another factor.

* `efa_reliability()` now points out that `fac_names` is not used for a `lavaan` fit, whose factor labels come from the model syntax, and that `group_names` is only used for a `lavaan` multiple-group fit, instead of dropping either without a word.

* `efa_scores(method = "Anderson")` and `FACTOR_SCORES(method = "Anderson")` now warn when the fitted factors are correlated, because Anderson-Rubin scores are uncorrelated even when the fitted factors are not.

* `efa_scores()` and `FACTOR_SCORES()` now reject scoring data in which a model variable is constant, infinite, or observed fewer than twice, naming it, instead of standardising it into `NaN` scores for every case.

* `efa_scores()` and `FACTOR_SCORES()` now warn when cases cannot be scored because they are missing a model variable. `efa_scores()` additionally reports the number of scored cases in `settings$n_scored`, and its print method now says how many of the supplied observations were scored and where the scores are stored.

* `efa_scores()` and `FACTOR_SCORES()` now warn when `Phi` is supplied together with an `efa_fit()` solution, whose own factor intercorrelations are used instead.

* `efa_schmid_leiman()` and `SL()` now reject the unsupported `se`, `b_boot`, `ci`, and `seed` arguments.

## Standard Errors

* `efa_fit(se = "sandwich")` now withholds standard errors and confidence intervals when a Heywood case occurs.

* Bootstrap output now reports the number of usable replicates whenever fewer than requested are available, for both single fits and pooled `efa_mi()` fits.

* Analytic standard-error outputs now include factor names in `SE$Phi` and variable-factor labels in `vcov_unrot_loadings`.

# EFAtools 1.0.0

## New Interface

* The recommended interface now consists of lowercase `efa_*` functions, with estimation and rotation settings defined through `estimate_control()` and `rotate_control()`.

* The existing uppercase functions remain available without warnings, so current scripts continue to work, but new features will generally be added only to the `efa_*` interface.

* Estimators, rotations, retention criteria, scoring methods, and presets can now be specified without matching capitalization exactly, while results continue to store and display the standard spelling.

## Changes When Using the `efa_*` Interface

* The `efa_*` functions now report an error for unused, misspelled, or misplaced arguments instead of silently ignoring them.

* Estimation and rotation options such as convergence settings, iteration limits, and rotation parameters should now be supplied through `estimate_control()` or `rotate_control()`.

* Factor-retention functions always fit unrotated models and therefore reject rotation settings that would not affect the retention result.

* The estimator argument is now named `estimator` rather than `method`, and rotation arguments use the names `p_type` and `random_starts` rather than `P_type` and `randomStarts`.

* `efa_retain()` requires the full argument name `max_iter_CD`, because the shorter form `max_iter` is no longer matched automatically.

* The second-order estimation preset for `efa_schmid_leiman()` is now specified with `estimate_control(type = ...)`.

* Control objects are validated when they are created, allowing invalid settings to be detected before a model is fitted.

## Reproducibility

* A supplied `seed` now controls the complete analysis, including random rotation starts, bootstrap samples, group-specific fits, and model averaging where applicable.

* `efa_average()` now has a `seed` argument and produces the same result regardless of the number of parallel workers.

* Seeded results from `efa_parallel()` and `efa_hull()` no longer depend on the parallel-processing setup.

* Because random-number handling and several calculations were corrected, the same seed may produce different results than in earlier versions for `efa_parallel()`, `efa_hull()`, `efa_average()`, `efa_group()`, and some `efa_simulate()` settings.

## Standard Errors and Fit Statistics

* `efa_fit(se = "information")` now returns `NA` with an `efa_se_unreliable` warning when standard errors for unrotated loadings cannot be estimated reliably, such as for Heywood cases or poorly distinguished factors.

* Information-based standard errors are now calculated for correlation matrices rather than covariance matrices, improving their accuracy and agreement with sandwich and bootstrap estimates.

* Information-based standard errors are now supported only with Pearson correlations and FIML; use `se = "sandwich"` for polychoric, tetrachoric, or rank correlations.

* FIML-based analytic standard errors now use the same small-sample scaling as the package’s other standard errors and fit statistics.

* `efa_mi()` no longer reports AIC, BIC, or ECVI when these measures are not meaningful for the underlying fitted models.

## Multiple Imputation and Group Comparisons

* The pooled unrotated solution from `efa_mi()` no longer depends on the order of the imputed datasets.

* Pooled unrotated loadings from `efa_mi()` now use the same orientation conventions as a regular `efa_fit()` result, making the two easier to compare.

* These changes can affect unrotated loadings, their standard errors, and statistics derived from them, but rotated solutions and their fit indices remain unchanged.

* Consensus results from `efa_group()` no longer depend on the order in which groups are supplied, so congruences, loading differences, salience classifications, and invariance conclusions are now stable across group orderings.

## Factor Retention and Rotation

* Fixed a regression in the revised MAP criterion, which could previously suggest an incorrect number of factors.

* `efa_smt()` now rejects polychoric and tetrachoric correlations because its tests and fit measures are not valid for these inputs.

* Rotated factors are now labelled `F1`, `F2`, and so forth according to the variance they explain, rather than inheriting labels that could change with the random seed.

* Rotation diagnostics now distinguish between the total number of starting solutions and the number that were fully optimized; code using `n_starts` should use `n_starts_total` or `n_optimized` instead.

* Criterion-based rotations now respect the specified `maxit` limit during every stage of the optimization.

## Model Averaging

* `efa_average()` no longer prints `NaN%` when no models converge and now clearly states which fitted models are included in each reported rate.

* Named presets such as `"SPSS"`, `"psych"`, and `"EFAtools"` now use their intended iteration limits, so models that exceed those limits are marked as non-converged and excluded from the average.

## Data Simulation

* Cudeck–Browne model-error simulations are now reproducible across computing platforms and numerical libraries, although a given seed will produce different populations than in earlier versions.

* With `marginals = "VM"`, categorized data now reproduce the requested category proportions more accurately, with a warning when the requested thresholds cannot be transformed safely.

* The documentation now clarifies that `match = "thresholds"` and `match = "polychoric"` produce the same data for normal marginals, and the selected value is now stored in the result settings.

* The documentation now explains the missing-data mechanism produced by `missing = "MAR"` and the residual bias it may create for analyses based only on the observed data.

* The documentation now explains that `model_error = "WB"` generally produces a larger RMSEA for the generating model than the requested target.

## Input Validation and Correlation Handling

* `efa_screen()` now handles collinear or nearly collinear variables without failing during its outlier checks.

* Polychoric covariance calculations now handle sparse or nearly perfectly associated response tables more consistently, allowing DWLS analyses to continue while warning about problematic variable pairs.

* Covariance matrices supplied where raw data or a correlation matrix is expected are now rejected with a message suggesting `cov2cor()`.

* Polychoric and tetrachoric correlations now require ordered factors or numeric response codes, preventing character or unordered-factor levels from being interpreted in alphabetical order.

## Power Analysis and Factor Scores

* `efa_power()` now clearly documents that `N` is the total sample size across groups and returns the corresponding sample size per group in `N_per_group`.

* Simulation-based power summaries now distinguish structure-recovery rates from the underlying Tucker congruence values.

* The documentation for `FACTOR_SCORES()` now correctly describes `R2` as squared factor-score determinacy.


# EFAtools 0.8.0

## New Functions

* `efa_group()` performs multi-group EFA.
* `efa_power()` to conduct power analysis (analytic power based on RMSEA and simulation based power analysis).
* `efa_reliability()` to calculate various reliability indices.
* `efa_screen()` to screen data for multivariate normality and suitability for factor analysis.
* `efa_simulate()` to simulate data with various distributions and missing data mechanisms.


## Changes to Functions

* `EFA()`: 
  * gains `cor_method = "fiml"` for raw data with missing values.
  * gains `se = "sandwich"` and `se = "information"`.
  * renamed the top-level fields `boot.SE`, `boot.CI`, and `boot.arrays` on
    `EFA` objects to `SE`, `CI`, and `replicates`, as the old names are no longer 
    accurate given the additional SE implementations.
  * All criterion-based rotations are now computed by a built-in gradient-projection
    rotation engine instead of the `GPArotation` package. The rotated solutions
    are numerically equivalent but implemented in C++ for speed to allow high 
    numbers of random starts.
  * The default number of random starts for the rotation in `EFA()` (`randomStarts`)
    has been raised from 10 to 100, making local minima less likely for the 
    rotation criteria that are prone to them.
  * additionally reports the Tucker-Lewis index (TLI, also called the
    non-normed fit index), the expected cross-validation index (ECVI), and the
    standardized root mean square residual (SRMR)
  * gains a `seed` argument and now fits the non-parametric bootstrap replicates
    (`se = "np-boot"`) in parallel across replicates via the `future` framework.
  * now has a `summary()` method that prints the full, detailed solution 
    (loadings with communalities and uniquenesses, explained variances, fit
    indices, and residuals); `print()` now gives a more compact overview.
* `EFA_AVERAGE()` gains `cor_method = "fiml"` (passed to `EFA()`).
* `EFA_POOLED()`:
  * now dispatches its multiple-imputation standard-error pooling automatically
    by the standard-error method of its component fits: Rubin's rules for the
    information method (with Wald confidence intervals), the two-stage pooled-inputs
    (MI2S) approach for the sandwich method, and the existing bootstrap pooling
    for the non-parametric bootstrap method. 
  * Renamed the top-level field `boot.MI` to `MI`. 
  * now defaults to `target_method = "first_target"`, which aligns every imputation
    to the first imputation's rotated solution with a single Procrustes rotation.
    For orthogonal rotations this reaches the same pooled estimate as `"consensus"`
    with substantially less work. `"consensus"` (Generalized Procrustes Analysis of
    the imputation-specific rotated loadings) is still available as an opt-in for
    orthogonal rotations, and now raises a classed condition
    (`efa_consensus_oblique_unsupported`) when combined with an oblique rotation,
    because the iterated-oblique-Procrustes centroid can drift to degenerate targets.
* `cor_method` now accepts `"poly"` and `"tetra"` to compute polychoric and tetrachoric
  correlations from raw ordinal (respectively binary) data, using a two-step estimator with
  no empty-cell continuity correction. Supported by `EFA()` (including its non-parametric
  bootstrap), `EFA_AVERAGE()`, the suitability tests `KMO()` and `BARTLETT()`, and the
  retention criteria `EKC()`, `KGC()`, `MAP()`, `SCREE()`, `SMT()`, and `N_FACTORS()`. The
  criteria that compare the data against simulated continuous reference data (`CD()`,
  `PARALLEL()`, `NEST()`, and `HULL()`) do not support `"poly"` / `"tetra"` and signal an
  error.
* The factor-retention functions `CD()`, `EKC()`, `HULL()`, `KGC()`, `MAP()`, `NEST()`,
  `PARALLEL()`, `SCREE()`, and `SMT()` now return objects of a common `efa_retention` class
  with shared `print()` and `plot()` methods.
* `COMPARE()` objects now have a `plot()` method that returns a `ggplot` object. The plot is
  no longer drawn as a side effect of `print()`.
* Console output (the `print`, `summary`, and `format` methods) is now produced with the `cli`
  package, and the messages, warnings, and errors emitted across the package are now S3-classed
  conditions, which makes them easier to handle programmatically.
* `EFA()`, `SL()`, and `EFA_AVERAGE()` now accept `method = "MINRES"` as a synonym for
  `method = "ULS"`. Minimum residual and unweighted least squares are two names for the same
  estimator and return identical results.

## Bug Fixes

* `EFA()` and `EFA_POOLED()`: The comparative fit index (CFI) now floors the model and
  baseline noncentralities at zero before taking their ratio (Bentler, 1990), so it is no
  longer deflated toward zero when the baseline (independence) model fits comparatively well.
  The value is unchanged for well-fitting models and now matches `lavaan`. CFI can change for
  solutions in which the baseline model does not misfit much.
* `EFA_POOLED()`: Corrected and extended the pooled chi-square-based model-fit indices. The D2
  average relative increase in variance is now centred on the mean of the square-root
  statistics (Li, Meng, Raghunathan & Rubin, 1991), removing a one-sided inflation of the
  pooled chi-square, RMSEA, and CFI. Bootstrap/MI confidence intervals for
  the pooled fit indices are now the Rubin-Wald multiple-imputation summaries; a miscalibrated
  bootstrap-percentile interval (obtained by re-running the pooling algorithm over matched
  replicate indices) was removed, as was a mislabeled pooled `Fm`.
* `EFA_AVERAGE()`: When every averaged solution fails (all runs error, fail to converge, or
  are Heywood cases), the function now returns an empty (`NA`) averaged result instead of
  erroring or averaging an empty set.
* `EFA()` with `method = "PAF"`: the reported number of iterations (`iter`) is now the
  number of iterations actually performed; it was previously one too high. The loadings,
  communalities, and convergence status are unchanged.
* `ULS` extraction: the minimised objective is now the sum of squared off-diagonal residuals,
  consistent with its analytic gradient and the reported minimum (`Fm`). The diagonal
  residuals were previously included in the objective but not in its gradient. The fitted
  loadings are unchanged to within optimiser tolerance.
* `NEST()` and `PARALLEL()`: A failed eigendecomposition of a degenerate simulated matrix
  now raises a clear error instead of resulting in undefined behaviour.
* The chi-square model-fit statistic is now the Bartlett-corrected maximum likelihood
  discrepancy evaluated at the model-implied correlation matrix, for both `ML` and `ULS`
  extraction. For `method = "ML"` this matches `stats::factanal()` and `psych::fa()`; the
  small-sample Bartlett correction was previously omitted. For `method = "ULS"` it is now
  a proper chi-square-distributed statistic matching `psych::fa(fm = "minres")`;
  previously the least-squares residual sum of squares was multiplied by `N - 1` and read
  as if it were Wishart-distributed, which produced an invalid p-value. The independence
  (baseline) model used for the CFI is rescaled onto the same discrepancy scale.
  Consequently the p-value, CFI, RMSEA (and its confidence interval), AIC, and BIC change
  for ML and ULS solutions, and the number of factors suggested by `SMT()` and `HULL()`
  may change for these methods.
* `PARALLEL()`: The percentile reference eigenvalues are now computed with
  `stats::quantile()`, correcting a slight off-by-one in the previous indexing.
* `EFA()`: For oblique rotations, the factor intercorrelations (`Phi`), the structure
  matrix, and the explained variances are now reflected and reordered consistently with
  the rotated loadings. Previously, when a factor was reflected to a positive orientation
  the factor intercorrelations were not sign-adjusted (so the structure matrix and
  reported correlations did not match the loadings), and with `order_type = "ss_factors"`
  the factor intercorrelations were not reordered at all.
* `EFA()`: The returned rotation matrix (`rotmat`) is now reflected and reordered
  consistently with the rotated loadings, for both orthogonal and oblique rotations, so
  that the rotated loadings can be reconstructed from the unrotated loadings and `rotmat`.
  Previously the sign reflection was not applied to it and its factors were left in a
  different order, so this reconstruction did not hold whenever a factor was reflected or
  reordered.
* `HULL()`: The convex-hull elimination now tests every triplet of adjacent solutions,
  including the one formed by the last interior solution. Previously this final triplet
  was skipped, so a solution lying below the line connecting its neighbours could remain
  on the hull. This can change the number of factors suggested by `HULL()` (and hence by
  `N_FACTORS()`).
* `NEST()`: When the test accepts every eigenvalue it examines (no empirical eigenvalue
  falls at or below its reference within the search range), it now retains that last
  accepted model rather than the model with one fewer factor. The search range is also
  bounded so that the reference model fitted at each step stays over-identified.
* `PARALLEL()`: When every real eigenvalue exceeds its reference, so the decision rule finds
  no crossing, the suggested number of factors is now all retained components, reported with
  a warning, instead of a silent `NA` (matching the convention used by `EKC()`).
* `EFA(se = "np-boot")`: The non-parametric bootstrap no longer repeats per-replicate
  warnings (about arguments pinned alongside `type`, and about the iterative fit reaching its
  maximum number of iterations) once per replicate. They are now suppressed during
  resampling, and non-convergence across replicates is reported once as a summary.
* `COMPARE()`: With `reorder = "congruence"`, the columns of `y` are now matched to those of
  `x` by an optimal one-to-one assignment that maximizes the total Tucker congruence, rather
  than by an independent per-factor best match. The greedy match could assign two factors of
  `x` to the same column of `y`, duplicating that column and dropping another, which corrupted
  the reported differences, factor correspondences, and root mean squared distance. The result
  is unchanged whenever the greedy match was already one-to-one.


# EFAtools 0.7.1

## Changes to Functions

* `EFA()`: Added `randomStarts` argument passed to GPArotation functions, as suggested by Coen Bernaards.
* `FACTOR_SCORES()`: Added `rho` argument, thanks to Andreas Soteriades.
* `EFA_POOLED()`: Fixed issue that could lead to averaged Phi not being symmetric.

# EFAtools 0.7.0

## New Functions

* `MAP()` computes the Velicer MAP criterion (both TR2 and TR4).
* `PROCRUSTES()` to perform orthogonal and oblique Procrustes / target rotation.
* `CONSENSUS_PROCRUSTES()` to perform Procrustes on a list of targets to obtain a common target.
* `residuals.EFA()` to extract and print residuals and, if computed, standardized residuals.
* `EFA_POOLED()` to run EFA on a list of multiple imputated datasets and pool the results. Thanks to Andreas Soteriades for the suggestion and first implementation.
* `print.EFA_POOLED()`, print method adapted from `print.EFA()`

## Changes to Functions

* `N_FACTORS()` can also compute the MAP criterion.
* `EFA()`:
  * Now returns and prints residuals and, if SEs are computed, standardized residuals.
  * Calculates and prints RMSR.
  * Can calculate bootstrap standard errors and CIs of parameters and fit indices.
* `print.EFA()` now prints more information.

## Bug Fixes

* Fixed a bug in `OMEGA()` that led to incorrect omega, H, and ECV values for `lavaan` bifactor models. Tnanks to Christopher D. King for bug report and suggested fix.
* Small fix in the documentation of `EFA_AVERAGE()`.
* Fixed incorrect calculation of RMSEA in `EFA()`.


# EFAtools 0.6.1

## Changes to Functions

* `EKC()`: Now correctly returns the number of factors based on the first time the eigenvalues drop below the references values.

# EFAtools 0.6.0

## New Functions

* Added `NEST()` to perform the Next Eigenvalue Sufficiency Test (Achim, 2017).

## Changes to Functions

* `EKC()`: The implementation based on Auerswald and Moshagen (2019) used in previous versions differed from the original implementation by Braeken and van Assen (2017). Now both versions are implemented and can be selected with the new `type` argument. Thanks to Luis Eduardo Garrido for pointing this out and to Johan Braeken for sharing sample code, based on which the original version is now implemented.
* `N_FACTORS()`: 
  * Updated default settings to only use often used and well performing factor retention methods (others, like the Kaiser Guttman criterion can still be used).
  * Added NEST as additional factor retention method.
  * New arguments: `ekc_type`, `alpha_nest`, and `n_datasets_nest` to account for the changes in `EKC()` and to control `NEST()`.


# EFAtools 0.5.0

## Changes to Functions

* `print.EFA()`: 
  * Now returns explained variance from rotated, rather than unrotated solution, if a rotation was performed.
  * Now prints communalities and uniquenesses in loading/pattern matrix
* `EFA()`: Calculate and return model implied correlation matrix.

# EFAtools 0.4.6

## General

* Small change in test of `.gof()`: Changed some tests to take care of the ATLAS issue when using R-devel on x86_64 Fedora 34 Linux with alternative BLAS/LAPACK.

# EFAtools 0.4.5

## Bug Fixes

* Updated `OMEGA()` to accommodate changes in the upcoming version of `psych::schmid()`

# EFAtools 0.4.4

## Changes to Functions

* `print.EFA()`: Added arguments `cutoff`, `digits` and `max_name_length` that are passed to `print.LOADINGS()`.
* `print.LOADINGS()`: New Argument `max_name_length` to control the maximum length of the displayed variable names (names longer than this will be cut on the right side). Previously, this was fixed to 10 (which is now the default).

## Misc

* Updated a test of a helper function (`.gof()`) that threw an error when using R-devel on x86_64 Fedora 36 Linux with alternative BLAS/LAPACK.
* added `dontrun` to examples of `EFA_AVERAGE()` and its print and plot methods as these were causing issues on R-oldrel which were not directly related to EFAtools and thus could not be fixed from within the package.

# EFAtools 0.4.3

## Changes to Functions

* `.gof()`: Changed the helper function to take care of the MKL issue when using R-devel on x86_64 Fedora 34 Linux with alternative BLAS/LAPACK.


# EFAtools 0.4.2

## Changes to Functions

* `.is_cormat()`: Changed the helper function to better detect wheter a matrix is a correlation matrix.
* `PARALLEL()`: Added a check, testing whether N > n_vars and throw an error if this is not the case.

## Bug Fixes

* Fixed some tests due to upcoming changes in the psych package which EFAtools depends on.


# EFAtools 0.4.1

## Bug Fixes

* Minor fixes in tests to solve problems on macOS m1.


# EFAtools 0.4.0

## Changes to Functions

* `EFA()`: Changed error to warning when model is underidentified. This allows the Schmid-Leiman transformation to be performed on a two-factor solution.
* `OMEGA()`: Added calculation of additional indices of interpretive relevance (H index, explained common variance [ECV], and percent of uncontaminated correlations [PUC]). This is optional and can be avoided by setting `add_ind = FALSE`.

## Bug Fixes

* `CD()`: Added `na.omit()` to remove missing values from raw data to avoid an error in the comparison-data procedure.


# EFAtools 0.3.1

## General
* When testing for whether a matrix is singular and thus smoothing should be done, test against .Machine$double.eps^.6 instead of 0, as suggested by Florian Scharf. 

## Changes to Functions

* `EFA()`: 
    * Added warnings if `type = "SPSS"` was used with `method = "ML"` or `method = "ULS"`, or with a rotation other than `none`, `varimax` or `promax`.
    * Avoided smoothing of non-positive definite correlation matrices if `type = "SPSS"` is used.
    * Use Moore-Penrose Pseudo Inverse in computation of SMCs if `type = "psych"` is used, by calling `psych::smc()`.
    * Use `varimax_type = "kaiser"` if `type = "EFAtools"` is used with `varimax` or `promax`.

## Bug Fixes
* `EFA_AVERAGE()`:
    * Added `future.seed = TRUE` to call to `future.apply::future_lapply()` to prevent warnings.
    * Fixed test for Heywood cases from testing whether a communality or loading is greater than .998, to only test whether communalities exceed 1 + .Machine$double.eps
* `print.EFA()`: Fixed test for Heywood cases from testing whether a communality or loading is greater than .998, to only test whether communalities exceed 1 + .Machine$double.eps
* `OMEGA()`: Small bugfix when `lavaan` second-order model is given as input


# EFAtools 0.3.0

## General
* Added examples for `EFA_AVERAGE()` to readme and the EFAtools vignette
* Updated examples in readme and vignettes according to the updated `OMEGA` function

## New Functions

* Added function `EFA_AVERAGE()` and respective print and plot methods, to allow running many EFAs across different implementations to obtain an average solution and test the stability of the results.

## Changes to Functions

* `EFA()`: Defaults that were previously set to `NULL` are now mostly set to `NA`. This was necessary for `EFA_AVERAGE()` to work correctly.
* `PARALLEL()`: Rewrote the generation of random data based eigenvalues to be more stable when SMCs are used.
* `OMEGA()`: Changed expected input for argument `factor_corres` from vector to matrix. Can now be a logical matrix or a numeric matrix with 0's and 1's of the same dimensions as the matrix of group factor loadings. This is more flexible and allows for cross-loadings.


# EFAtools 0.2.0

## General

* Created new vignette *Replicate_SPSS_psych* to show replication of original `psych` and `SPSS` EFA solutions with `EFAtools`.

## New Functions

* Added function `FACTOR_SCORES()` to calculate factor scores from a solution from `EFA()`. This is just a wrapper for the `psych::factor.scores` function.
* Added function `SCREE()` that does a scree plot. Also added respective print and plot
methods.

## Changes to Functions

* `CD()`: Added check for whether entered data is a tibble, and if so, convert to vanilla data.frame to avoid breaking the procedure.
* `EFA()`: 
    * Updated the EFAtools type in PAF and Promax.
    * Added p value for chi square value in output (calculated for ML and ULS fitting methods).
    * Updated the SPSS varimax implementation to fit SPSS results more closely.
    * Created an argument "varimax_type" that is set according to the specified type, but that can also be specified individually. With type R psych and EFAtools, the stats::varimax is called by default (`varimax_type = "svd"`), with type SPSS, the reproduced SPSS varimax implementation is used (`varimax_type = "kaiser"`).
    * Renamed the `kaiser` argument (controls if a Kaiser normalization is done or not) into `normalize` to avoid confusion with the `varimax_type` argument specifications.
* `ML()`: Changed default start method to "psych".
* `N_FACTORS()`:
    * Added option to do a scree plot if "SCREE" is included in the `criteria` argument.
    * Added a progress bar.
* `OMEGA()`: Now also works with a lavaan second-order solution as input. In this case, it does a Schmid-Leiman transformation based on the first- and second-order loadings first and computes omegas based on this Schmid-Leiman solution.
* `SL()`: Now also works with a lavaan second-order solution as input (first- and second-order loadings taken directly from lavaan output).

## Bug Fixes

* `.get_compare_matrix()`: Fixed a bug that occurred when names of data were longer than n_char
* `COMPARE()`: Fixed a bug that occurred when using `reorder = "names"`.
* `EFA()`: RMSEA is now set to 1 if it is > 1.
* `HULL()`: Fixed a bug that occurred when no factors are located on the HULL
* `KMO()`: Fixed a bug that the inverse of the correlation matrix was not taken anew after smoothing was necessary.
* `PARALLEL()`:
    * Fixed a bug that occurred when using `decision_rule = "percentile"`
    * Relocated error messages that were not evaluated if no data were entered (and should be)
* `print.COMPARE()`: Fixed a bug that occurred when using `print_diff = FALSE` in `COMPARE()`.
* `print.KMO()`: Fixed a bug that printed two statements instead of one, when the KMO value was < .6.

## Minor Changes
* `OMEGA()` and `SL()`: Added an error message if the entered term in `g_name` is invalid (i.e., it cannnot be found among the factor names of the entered lavaan solution).


# EFAtools 0.1.1

## Minor Changes

* Added an error message in `PARALLEL()` if no solution has been found after 25 tries.

## Bug Fixes

* Updated different tests

* Deleted no longer used packages from Imports and Suggests in DESCRIPTION

* `PARALLEL()`: fixed a bug in indexing if method `"EFA"` was used.


# EFAtools 0.1.0

* Added a `NEWS.md` file to track changes to the package.
* Initial CRAN submission
