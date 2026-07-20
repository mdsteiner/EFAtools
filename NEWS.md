# EFAtools 0.8.0.9000

## New Interface

* The public interface of the package is now the lowercase `efa_*` functions,
configured through `estimate_control()` and `rotate_control()`. `efa_fit()` fits
the factor models, `efa_retain()` and the individual retention criteria determine
the number of factors, and the other steps of an analysis are covered by
`efa_screen()`, `efa_average()`, `efa_mi()`, `efa_group()`,
`efa_schmid_leiman()`, `efa_procrustes()`, `efa_compare()`, `efa_reliability()`,
`efa_scores()`, `efa_simulate()`, and `efa_power()`.

* The uppercase function names are superseded by their `efa_*` equivalents. They
remain exported, keep their arguments, and emit no warning, so existing code
needs no changes. New arguments and features are added to the `efa_*` functions
only, which are the recommended interface for new code.

* The `efa_*` functions and the control constructors match the mixed-case choice
values case-insensitively: the estimator, the rotation, the factor-retention
criteria, the eigenvalue and goodness-of-fit types, the factor-score method, and
the estimation and rotation presets (`type`). For example,
`efa_fit(x, n_factors = 3, estimator = "ml", rotation = "geomint")` selects the
`"ML"` estimator and the `"geominT"` rotation. The canonical spelling is always
stored in the returned settings and used in the printed output.

## Breaking Changes

Code written against the superseded uppercase functions keeps working. The
following concern the `efa_*` interface.

* Arguments that do not belong to a fitting function are now an error rather than
silently ignored. `efa_fit()`, `efa_retain()`, `efa_schmid_leiman()`, and the
retention criteria that fit a model (`efa_hull()`, `efa_kgc()`, `efa_nest()`,
`efa_parallel()`, and `efa_scree()`) reject a tuning knob passed to them directly
(such as `type`, `max_iter`, `criterion`, or `p_type`), as these now live in
`estimate_control()` and `rotate_control()`. `efa_fit()` also rejects a `...`
argument the selected rotation's engine cannot consume — accepted are `maxit`
for the GPArotation-style rotations plus the selected criterion's parameter
(`gam` for oblimin, `delta` for geomin), and none at all for varimax, promax,
and `rotation = "none"`, where no engine runs. `rotate_control()` likewise
validates its extra engine arguments when the control is created. A misspelling
such as `gamma = 0.5` (for oblimin's `gam`) is therefore an immediate error
instead of a fit that silently runs the engine defaults. `EFA()` keeps silently
dropping an extra that its selected rotation cannot consume, so code written
against it is unaffected.

* `efa_retain()` takes the maximum number of comparison-data iterations as
`max_iter_CD`, and this name has to be given in full. The argument sits behind
`...`, where R no longer matches it from the abbreviation `max_iter`. Such a call
is rejected as a misplaced estimation knob, rather than silently capping the
comparison-data iterations as it does in `N_FACTORS()`.

* `estimate_control()` validates its settings when the control object is created,
so a convergence `criterion` of 1 or larger is now rejected there (a classed
`efa_control_input` error) instead of inside the fit.

* The `efa_*` functions select the estimator with `estimator` (e.g.
`efa_fit(x, n_factors = 3, estimator = "ML")`). Passing the former name `method`
to one of them is an error rather than a silently ignored argument; the superseded
uppercase functions keep their `method` argument, and the returned `settings` keep
a `method` entry duplicating `estimator`, so existing code keeps working.

* The rotation arguments `P_type` and `randomStarts` are now named `p_type` and
`random_starts`. `EFA()` still accepts the former names silently, and
`EFA_AVERAGE()` keeps `P_type` as its argument name, so existing code keeps
working.

* `efa_schmid_leiman()` has no `type` argument; the estimation preset for the
second-order fit is set with `estimate_control(type = ...)`. `SL()` keeps `type`.

## Reproducibility

* `seed` now governs the whole fit, not only the bootstrap. In `efa_fit()` it covers the
rotation's random starts on the point estimate as well as the case resampling and the
replicate fits, so `efa_fit(..., rotation = "simplimax", seed = 1)` is now reproducible
at any `se`. Previously `seed` was silently a no-op at the default `se = "none"`, and a
criterion-based rotation could reach a different optimum on every call. Seeded
`se = "np-boot"` results are unchanged: there the point fit already drew from the seeded
stream, so only the previously unreachable case is affected. `EFA()` forwards `seed` to
`efa_fit()` and inherits the wider coverage. `efa_group()` likewise applies `seed` to the
per-group fits whether or not a congruence bootstrap is requested, where it previously
applied it only with `b_boot > 0`.

* `efa_average()` gains a `seed` argument with the same contract: it makes the averaging
grid reproducible and independent of the number of parallel workers, and the caller's
random-number stream is restored afterwards.

* `efa_parallel()` results no longer depend on the number of parallel workers. The
simulated datasets were split into as many chunks as there were workers, and each chunk
drew its own random-number stream, so the reference eigenvalues from a single
`set.seed()` differed between a sequential and a multisession plan. The number of chunks
is now fixed, which makes a seeded run reproducible on any plan but changes the reference
eigenvalues a given seed produces. `efa_retain()` and `efa_hull()` inherit this.

* `efa_hull()` now runs its fits with a managed random-number stream, so a hull analysis
with a criterion-based rotation (such as `oblimin`) is reproducible from `set.seed()`
under a parallel plan. This changes the results a given seed produces for those
rotations.

## Bug Fixes

* `efa_fit(se = "information")` no longer reports enormous standard errors for the
unrotated loadings of a weakly determined solution. Two safeguards that were meant to
catch this never fired. The first withholds the analytic covariance at a Heywood case, but
tested for a uniqueness at or below zero — which the ML, ULS and DWLS fitters cannot
produce, because they constrain the uniquenesses to a small positive floor. An improper
solution is pinned *at* that floor rather than below zero, so it passed the check and was
reported with a Wald standard error taken at a boundary where the interval is not valid.
It is now recognised at the boundary the fitters can actually reach, the same one a
Heywood case is flagged by elsewhere. The second safeguard is new: when two of a solution's
canonical variances nearly coincide, the orientation of the unrotated loadings within that
plane is arbitrary and the constraint that identifies them stops pinning it down. The
parameter covariance stays finite and positive semidefinite, so nothing downstream flagged
it, yet the standard errors diverge along that direction — reaching 80 for a loading
bounded in [-1, 1]. Both cases now return `NA` for the affected standard errors and raise
the classed `efa_se_unreliable` warning. On a simulated design with a weak third factor at
N = 400, the 99th percentile of the mean unrotated loading standard error falls from 11.3
to 0.72, while the median is unchanged at 0.07.

Only the unrotated loadings are affected by the second case. The rotated loadings, the
factor correlations, the structure coefficients and the communalities do not depend on the
orientation and keep their standard errors, as does `$vcov_unrot_loadings`, so `efa_mi()`
continues to pool them. A Heywood case, where the Wald approximation fails altogether,
still withholds every analytic standard error as before. Where the unrotated loadings
themselves are the quantity of interest, `se = "np-boot"` remains available.

* `efa_mi()` no longer returns a different unrotated solution depending on the order of
`data_list`. The unrotated loadings were matched up to factor reordering and sign against
whichever imputation came first, and where that matching is ambiguous — a weakly
determined factor, or two factors of near-equal strength — a different first imputation
produced a different pooled unrotated solution, by up to 0.84 in a single loading. The
match is now anchored on the imputation closest to all the others, which is a property of
the set and so does not move when the list is reordered. Because the anchor has moved,
this changes the pooled unrotated solution for everyone, not only for those who reorder
their imputations, and with it everything computed from it: with `rotation = "none"` that
includes the communalities, the model-implied correlation matrix, the residuals, RMSR,
SRMR and CAF. On a weakly determined three-factor fit a communality moved by up to 0.14
and RMSR from 0.150 to 0.094. Chi-square, CFI, TLI and RMSEA are pooled from the
per-imputation fit indices and are unaffected. A whole factor column may also come back
with the opposite sign; a column sign is arbitrary in an unrotated solution and carries no
substantive meaning.

* The unrotated loadings pooled by `efa_mi()` are now returned in the same gauge the
extraction itself uses, so they can be compared with an `efa_fit()` solution
element-by-element. Averaging several aligned solutions does not preserve the constraint
that identifies the unrotated loadings — diagonal `L'L` for a principal-axis extraction,
diagonal `L' Psi^-1 L` for maximum likelihood — which left the pooled matrix in a gauge no
single fit uses. How far out of gauge the average drifts depends on how well separated the
factors are: on well-determined three-factor solutions the constraint was violated by well
under 1%, rising to 28% on a weakly determined one. This second correction, unlike the
change of anchor above, is an orthogonal rotation shared by all imputations, so on its own
it leaves the communalities, the total variance accounted for, the model-implied
correlation matrix, the residuals, RMSR, SRMR and every fit index exactly unchanged. What
it does move is the pooled *unrotated* loadings, their standard errors, confidence
intervals and MI diagnostics (RIV, FMI), and the unrotated per-factor variance table. The
standard errors move considerably more than the loadings do — typically by under 2%, but
by up to a factor of three where a factor is weakly determined.

* Neither change affects a rotated solution. If `efa_mi()` is called with a rotation, the
rotated loadings, `Phi`, `Structure`, the communalities, the rotated variance table and
every fit index are numerically identical to before; only the unrotated block that
`efa_fit()` also reports has changed.

* `se = "information"` now returns correlation-structure standard errors. The expected
information was assembled for a covariance structure — treating the uniquenesses as free
parameters and so attributing sampling variability to a diagonal that is fixed at 1 — which
made the loading standard errors conservative by up to a factor of two (most for the highest
loadings), with a nominal 95% Wald interval covering at about 99.8%. Real loadings and
cross-loadings could therefore be judged non-significant. The information is now built from
the normal-theory asymptotic covariance of the correlations at the model-implied `Sigma`
(Cudeck, 1989; Olkin & Siotani, 1976), which is Monte-Carlo calibrated and agrees with the
`"sandwich"` and `"np-boot"` standard errors on normal data. Loading standard errors,
confidence intervals, and the stored `vcov_unrot_loadings` (from which `efa_mi()` builds its
within-imputation variances) all change; `Phi`, the uniquenesses, and the communalities are
essentially unaffected.

* `se = "information"` is now rejected for `cor_method` other than `"pearson"` and `"fiml"`.
The formula presumes normal-theory sampling behaviour of the analysed matrix, which neither a
polychoric correlation (it carries first-stage threshold estimation error) nor a rank
correlation satisfies; both previously returned standard errors without comment. Use
`se = "sandwich"`, which is built for exactly those inputs.

* The corrected two-stage standard errors under `cor_method = "fiml"` now use the same
`N - 1` small-sample scaling as the other analytic standard errors and as the fit statistics.
They were previously scaled by `N`, so they differed from the complete-data Pearson path by a
systematic factor of `sqrt(N / (N - 1))`. Reported standard errors change by under 0.2%.

* Fixed a regression in `MAP()`: the revised (TR4) criterion was computed from
the element-wise fourth powers of the partial correlations instead of the trace
of the fourth matrix power, which could change the suggested number of factors.

* `efa_average()` no longer reports `NaN%` Heywood-case and admissibility rates
when no solution converged (or every solution errored). The printed summary now
names the denominator each rate is computed over — Heywood cases and
admissibility are assessed only for the solutions that converged — and omits a
rate when no solution reached that stage.

* `efa_simulate()` now reproduces its Cudeck-Browne model error across platforms. With
`model_error = "CB"` and a `target_rmsea`, the population perturbation drawn for a given
`seed` could differ by tens of percent between BLAS implementations, and between two
identical calls on a threaded BLAS. The Jacobian defining the perturbation is now
differentiated in closed form instead of by a finite difference, which removes the
rounding noise responsible (and is faster). Populations drawn from a given seed therefore
differ from those of earlier versions; the achieved RMSEA still equals the target.

* `efa_screen()` no longer aborts with a `system is computationally singular` error
from `solve()` when the outlier diagnostics fall back to the classical covariance
and the variables are collinear or nearly so. Whether the abort happened depended on
the platform's BLAS. The fallback now checks that the complete-case covariance is
invertible on the same terms as the Mahalanobis step that consumes it, so such data
reliably reach the documented `efa_screen_no_outliers` note instead.

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
