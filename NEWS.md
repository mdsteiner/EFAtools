# EFAtools 0.7.1.9000

## Changes to Functions

* `EFA_POOLED()` now dispatches its multiple-imputation standard-error pooling
  automatically by the standard-error method of its component fits: Rubin's
  rules for the information method (with Wald confidence intervals and
  Barnard-Rubin small-sample degrees of freedom), the two-stage pooled-inputs
  (MI2S) approach (Chung & Cai, 2019; Sriutaisuk, Liu, Chung, Kim, & Gu, 2025)
  for the sandwich method (a single fit on a Rubin-pooled correlation matrix and
  asymptotic covariance), and the existing bootstrap pooling for the
  non-parametric bootstrap method. Mixed standard-error methods across component
  fits now raise a classed `efa_pooled_mixed_se` condition rather than silently
  producing uninterpretable results, and a failure to produce pooled standard
  errors on the information or bootstrap route is signalled with a classed
  `efa_pooled_se_unavailable` warning (falling back to point-estimate-only
  pooling) rather than being silently dropped.

* `EFA_POOLED()` routes component fits with `se = "sandwich"` through the
  two-stage pooled-inputs (MI2S) approach (Chung & Cai, 2019; Sriutaisuk,
  Liu, Chung, Kim, & Gu, 2025): the correlation matrix and the asymptotic
  covariance of its off-diagonal entries are Rubin-pooled across imputations,
  and a single `EFA` fit on the pooled inputs produces native scaled-shifted
  chi-square test statistics and sandwich standard errors that reflect the
  multiple-imputation uncertainty. The single pooled fit is exposed as
  `mi_fit` on the returned object alongside the per-imputation `fits` list.
  Because there is a single fit there is one rotational gauge, so the route
  bypasses per-imputation rotation alignment (`target_method` and
  `align_unrotated` do not apply). The pooled asymptotic covariance is not
  smoothed: an indefinite pooled covariance at small m aborts with the
  classed condition `efa_pooled_mi2s_acov_not_psd`, and fewer than 20
  imputations raises a classed warning. The polychoric/tetrachoric case is the
  primary target; the continuous-Pearson case uses the same recipe but is less
  benchmarked under multiple imputation.

* `EFA_POOLED()` now pools the unrotated-loading and uniqueness standard errors
  produced by component fits under `se = "information"`, using Rubin's rules
  with Barnard-Rubin small-sample degrees of freedom and Wald confidence
  intervals. Under `align_unrotated = "procrustes"`, the per-imputation
  orthogonal Procrustes transform mixes loading columns, so the full
  unrotated covariance persisted on each fit (`vcov_unrot_loadings`) is
  propagated through the column-major Kronecker identity (treating the
  per-imputation Procrustes transform as fixed) before Rubin pooling; see
  `?EFA_POOLED` for the formula. Component fits missing the
  full covariance abort with the classed condition `efa_pooled_no_vcov`, and
  fits carrying an NA-filled covariance (the upstream signal of a Heywood
  case or singular bordered information matrix) abort with
  `efa_pooled_unreliable_vcov`.

* `EFA_POOLED()` now also pools the rotated-loading, communality, factor-
  correlation (`Phi`), and structure-coefficient standard errors from component
  fits under `se = "information"`. Orthogonal rotations use the closed-form
  orthogonal Procrustes alignment to the multiple-imputation rotated target and
  propagate the unrotated parameter covariance through the Kronecker identity;
  communalities are rotation-invariant and pool without alignment. Oblique
  rotations use a signed-permutation alignment of the per-imputation rotated
  standard errors as a documented approximation (flagged in
  `MI$<param>$method`); researchers requiring rigorous uncertainty
  quantification for oblique solutions should cross-check with `se = "np-boot"`.

* `EFA()` objects now carry through the full unrotated loading covariance matrix
  (`vcov_unrot_loadings`) for the analytic SE methods (`se = "information"` and
  `se = "sandwich"`, in both unrotated and rotated fits), and the asymptotic
  covariance of the off-diagonal correlations (`Gamma`) for `se = "sandwich"`.

* Renamed the top-level fields `boot.SE`, `boot.CI`, and `boot.arrays` on
  `EFA` objects to `SE`, `CI`, and `replicates`, and the additional
  `EFA_POOLED` field `boot.MI` to `MI`. The previous names were artefacts of
  internal list flattening and were misleading because `SE` and `CI` are
  populated under the `"information"` and `"sandwich"` SE methods as well, not
  only the nonparametric bootstrap. Code that read these fields directly must
  update its slot accessors.

* `cor_method` now accepts `"poly"` and `"tetra"` to compute polychoric and tetrachoric
  correlations from raw ordinal (respectively binary) data, using a two-step estimator with
  no empty-cell continuity correction. Supported by `EFA()` (including its non-parametric
  bootstrap), `EFA_AVERAGE()`, the suitability tests `KMO()` and `BARTLETT()`, and the
  retention criteria `EKC()`, `KGC()`, `MAP()`, `SCREE()`, `SMT()`, and `N_FACTORS()`. The
  criteria that compare the data against simulated continuous reference data (`CD()`,
  `PARALLEL()`, `NEST()`, and `HULL()`) do not support `"poly"` / `"tetra"` and signal an
  error.

* The `quartimax` and `equamax` orthogonal rotations in `EFA()` are now computed by a
  built-in gradient-projection rotation engine instead of the `GPArotation` package. The
  rotated solutions are numerically equivalent (the engine reaches the same minimum of the
  Crawford-Ferguson criterion); the remaining rotation criteria continue to use
  `GPArotation`.
* The `oblimin` and `quartimin` oblique rotations in `EFA()` are now computed by the same
  built-in gradient-projection rotation engine instead of the `GPArotation` package. The
  rotated solutions are numerically equivalent (the engine reaches the same minimum of the
  oblimin criterion); the remaining oblique rotation criteria continue to use
  `GPArotation`.
* The `geominT` and `geominQ` geomin rotations in `EFA()` are now computed by the same
  built-in gradient-projection rotation engine instead of the `GPArotation` package. Because
  the geomin criterion is prone to local minima, the engine reaches an equal-or-better minimum
  of the geomin criterion rather than reproducing `GPArotation`'s particular solution exactly;
  the remaining rotation criteria continue to use `GPArotation`.
* The `bentlerT` and `bentlerQ` Bentler rotations in `EFA()` are now computed by the same
  built-in gradient-projection rotation engine instead of the `GPArotation` package. The engine
  reaches an equal-or-better minimum of the Bentler invariant pattern simplicity criterion; the
  remaining rotation criteria continue to use `GPArotation`.
* The `bifactorT` and `bifactorQ` bifactor rotations in `EFA()` are now computed by the same
  built-in gradient-projection rotation engine instead of the `GPArotation` package. Because the
  bifactor criterion is prone to local minima, the engine reaches an equal-or-better minimum of the
  Jennrich-Bentler bifactor criterion rather than reproducing `GPArotation`'s particular solution
  exactly.
* The `simplimax` oblique rotation in `EFA()` is now computed by the same built-in
  gradient-projection rotation engine instead of the `GPArotation` package. Because the simplimax
  criterion is prone to local minima, the engine reaches an equal-or-better minimum of the
  simplimax criterion rather than reproducing `GPArotation`'s particular solution exactly. With
  this, every analytic rotation in `EFA()` is computed by the built-in engine, and `GPArotation`
  is no longer used to perform rotations.
* The default number of random starts for the rotation in `EFA()` (`randomStarts`) has been
  raised from 10 to 100, making local minima less likely for the rotation criteria that are
  prone to them. The rotation engine screens the random starts cheaply and fully optimizes only
  the most promising ones, so the higher default barely changes the runtime of a single `EFA()`
  for most criteria (a few seconds at most, for `simplimax`, which fully optimizes every start).
* `summary()` now reports the number of distinct local optima the rotation found across its
  random starts, and `EFA()` records the rotation diagnostics (the number of random starts, how
  many converged, the distinct-optima count, and the spread and best value of the attained
  criterion) in `settings$rotation_diagnostics`.
* The built-in gradient-projection rotation engine now uses a Barzilai-Borwein step size with a
  non-monotone line search. The analytic rotations in `EFA()`, and the `PROCRUSTES()` target
  rotation, reach the same solution in substantially fewer iterations; the speed-up is most
  noticeable on larger problems and in the resampling workflows that rotate many solutions
  (`EFA()` bootstrap standard errors and `EFA_AVERAGE()`).
* The oblique geomin rotation (`geominQ`) now triages more of its random starts by default, so it
  reliably reaches the global geomin optimum on problems where the previous setting occasionally
  settled at a nearby local optimum.
* `EFA()` now additionally reports the Tucker-Lewis index (TLI, also called the
  non-normed fit index), the expected cross-validation index (ECVI), and the
  standardized root mean square residual (SRMR) among its `fit_indices` for `ML` and
  `ULS` estimation (SRMR is also reported for `PAF`). The SRMR matches `lavaan`; the TLI
  and ECVI are based on the Bartlett-corrected chi square. The new indices are shown in
  `print()` and `summary()`.
* `EFA()` now detects Heywood (improper) cases in the fitted solution, records the affected
  variables in a new `heywood` element of the returned object, and emits a classed warning.
  Detection is consistent across `PAF`, `ML`, and `ULS`: a variable is flagged when its
  communality reaches or exceeds 1 (which can happen under `PAF`), or when its uniqueness is
  pinned at the estimator's lower bound (the boundary/improper case under `ML` and `ULS`).
* `EFA()` objects now have a `summary()` method that prints the full, detailed solution
  (loadings with communalities and uniquenesses, explained variances, fit indices, and
  residuals); `print()` now gives a more compact overview.
* The factor-retention functions `CD()`, `EKC()`, `HULL()`, `KGC()`, `MAP()`, `NEST()`,
  `PARALLEL()`, `SCREE()`, and `SMT()` now return objects of a common `efa_retention` class
  with shared `print()` and `plot()` methods.
* `COMPARE()` objects now have a `plot()` method that returns a `ggplot` object. The plot is
  no longer drawn as a side effect of `print()`.
* Console output (the `print`, `summary`, and `format` methods) is now produced with the `cli`
  package, and the messages, warnings, and errors emitted across the package are now S3-classed
  conditions, which makes them easier to handle programmatically.
* `lavaan` moved from Imports to Suggests. The `lavaan`-input paths of `OMEGA()` and `SL()` now
  require the `lavaan` package to be installed and raise a clear error if it is missing.
* `GPArotation` moved from Imports to Suggests. All rotations in `EFA()` are computed by the
  package's built-in rotation engines, so `GPArotation` is no longer needed to run the package;
  it is now used only as an optional reference in the test suite.
* `EFA()`, `SL()`, and `EFA_AVERAGE()` now accept `method = "MINRES"` as a synonym for
  `method = "ULS"`. Minimum residual and unweighted least squares are two names for the same
  estimator and return identical results.
* `EFA()` gains a `seed` argument and now fits the non-parametric bootstrap replicates
  (`se = "np-boot"`) in parallel across replicates via the `future` framework. By default the
  replicates are still fitted sequentially; register a plan with `future::plan()` (for example
  `future::plan(future::multisession, workers = 2)`) to fit them in parallel. With a fixed
  `seed`, the bootstrap is reproducible and returns the same result regardless of the number of
  workers. Because each worker runs its own linear algebra, keep the number of workers small if
  your `BLAS` is multi-threaded, to avoid over-subscribing the available cores.
* The bootstrapped standard errors and confidence intervals for rotated loadings from `EFA(se =
  "np-boot")` can differ slightly from earlier versions. The oblique re-rotation of the bootstrap
  replicates is now carried out in a single compiled pass (with its warm start computed in C++),
  and the order in which random numbers are drawn has changed. The results are unchanged in
  distribution and are now reproducible and independent of the number of workers.
* `EFA_POOLED()` now defaults to `target_method = "first_target"`, which aligns every imputation
  to the first imputation's rotated solution with a single Procrustes rotation. For orthogonal
  rotations this reaches the same pooled estimate as `"consensus"` with substantially less work.
  `"consensus"` (Generalized Procrustes Analysis of the imputation-specific rotated loadings,
  in the style of van Ginkel & Kroonenberg, 2014, and Lorenzo-Seva & Van Ginkel, 2016) is still
  available as an opt-in for orthogonal rotations, and now raises a classed condition
  (`efa_consensus_oblique_unsupported`) when combined with an oblique rotation, because the
  iterated-oblique-Procrustes centroid can drift to degenerate targets that the literature
  sidesteps by adding a Promin step that this engine does not implement.
* The internal GPA-consensus engine that used to be exported as `CONSENSUS_PROCRUSTES()` is no
  longer part of the public API. It is reached only via `EFA_POOLED(target_method = "consensus")`.

## Bug Fixes

* `EFA()` and `EFA_POOLED()`: The comparative fit index (CFI) now floors the model and
  baseline noncentralities at zero before taking their ratio (Bentler, 1990), so it is no
  longer deflated toward zero when the baseline (independence) model fits comparatively well.
  The value is unchanged for well-fitting models and now matches `lavaan`. CFI can change for
  solutions in which the baseline model does not misfit much.
* `EFA_POOLED()`: Corrected and extended the pooled chi-square-based model-fit indices. The D2
  average relative increase in variance is now centred on the mean of the square-root
  statistics (Li, Meng, Raghunathan & Rubin, 1991), removing a one-sided inflation of the
  pooled chi-square, RMSEA, and CFI; the pooled set now also reports the Tucker-Lewis index
  (TLI) and the expected cross-validation index (ECVI). Bootstrap/MI confidence intervals for
  the pooled fit indices are now the Rubin-Wald multiple-imputation summaries; a miscalibrated
  bootstrap-percentile interval (obtained by re-running the pooling algorithm over matched
  replicate indices) was removed, as was a mislabeled pooled `Fm`.
* `EFA_AVERAGE()`: When every averaged solution fails (all runs error, fail to converge, or
  are Heywood cases), the function now returns an empty (`NA`) averaged result instead of
  erroring or averaging an empty set.
* Factor extraction (`PAF`, `ML`, and `ULS`) now raises a clear error when the requested
  number of factors is not smaller than the number of variables, instead of reading past
  the available eigenvalues (undefined behaviour in builds without bounds checking).
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
  `stats::quantile()` (matching `psych::fa.parallel`), correcting a slight off-by-one in
  the previous indexing. In addition, the simulated datasets are now split across parallel
  workers with an exact integer partition, fixing an invalid negative chunk that could
  occur when the number of simulated datasets was close to the number of workers.
* Non-positive-definite correlation matrices are now reported with a single classed warning
  (`efa_cor_smoothed`) when they are smoothed, consistently across `EFA()`, `KMO()`,
  `BARTLETT()`, and the factor-retention functions.

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
* `EFA()` and `EFA_AVERAGE()`: The orthogonal `"quartimax"` rotation now optimizes the
  quartimax criterion. Previously it optimized the Bentler invariant pattern simplicity
  criterion, so the quartimax-rotated loadings change. `EFA(rotation = "quartimax")` now
  agrees with `psych::fa(rotate = "quartimax")` and `GPArotation::quartimax()`.
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
