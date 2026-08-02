#' Exploratory factor analysis (EFA)
#'
#' This function does an EFA with either `PAF`, `ML`, `ULS`/`MINRES`,
#' or `DWLS` with or without subsequent rotation.
#' All arguments with default value `NA` can be left to default if `type`
#' is set to one of "EFAtools", "SPSS", or "psych". The respective specifications are
#' then handled according to the specified type (see details).
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#' correlations. If raw data is entered, the correlation matrix is found from the
#' data.
#' @param n_factors numeric. Number of factors to extract. Must be at least 1 and
#' smaller than the number of variables (the common factor model is not identified
#' otherwise). Use [efa_retain()] to decide on a value.
#' @param N numeric. The number of observations. Needs only be specified if a
#' correlation matrix is used. If input is a correlation matrix and `N` = NA
#' (default), not all fit indices can be computed. When raw data with missing
#' values are entered and `use` is `"complete.obs"` or `"na.or.complete"`, rows
#' are deleted listwise, so `N` is taken as the number of complete cases. The same
#' applies, whatever `use` asks for, whenever an asymptotic covariance is required
#' (the "DWLS" estimator, or `se = "sandwich"`); with `cor_method = "fiml"`, `N` is
#' instead the number of cases carrying at least one observed value.
#' @param estimator character. The estimator used to fit the EFA: "PAF" (principal axis
#' factoring), "ML" (maximum likelihood), "ULS" (unweighted least squares; "MINRES" is an
#' accepted alias returning identical results), or "DWLS" (diagonally weighted least
#' squares, for ordinal data). See the *Estimators* section in Details for their
#' properties and data requirements. The value is matched case-insensitively.
#' @param rotation character. Either perform no rotation ("none"; default),
#' an orthogonal rotation ("varimax", "equamax", "quartimax", "geominT",
#' "bentlerT", or "bifactorT"), or an oblique rotation ("promax", "oblimin",
#' "quartimin", "simplimax", "bentlerQ", "geominQ", or "bifactorQ"). See the
#' *Rotations* section in Details for their properties and known issues.
#' @param se character. Whether and how to compute standard errors (and matching
#'  confidence intervals): "none" (default, no standard errors), "information" (analytic
#'  standard errors from the expected Fisher information of the ML solution), "sandwich"
#'  (robust Godambe sandwich standard errors from raw data), or "np-boot" (non-parametric
#'  bootstrap). The methods differ in their assumptions, their data requirements, and which
#'  estimator, rotation, and `cor_method` combinations they support; see the *Standard
#'  errors* section in Details.
#' @param use character. Passed to [stats::cor()] if raw data
#' is given as input. Default is "pairwise.complete.obs". It is ignored when
#' `cor_method = "fiml"` (which handles the missingness itself, so every case
#' contributes), and it is overridden to listwise deletion whenever an asymptotic
#' covariance is required (the "DWLS" estimator, or `se = "sandwich"`), because the
#' covariance must describe the same cases as the correlation matrix.
#' @param cor_method character. How the correlation is computed from raw data:
#'   `"pearson"`, `"spearman"`, or `"kendall"` (passed to [stats::cor()]); `"poly"` /
#'   `"tetra"` for polychoric / tetrachoric correlations of ordinal / binary data; or
#'   `"fiml"` for a two-stage full-information maximum-likelihood correlation from raw data
#'   with missing values. See the *Correlation methods* section in Details for their
#'   properties and the combinations they support. Default is "pearson".
#' @param estimate_control a control object from [estimate_control()] bundling the
#'  estimation tuning knobs: the `type` preset; the principal-axis-factoring iteration
#'  settings `init_comm`, `criterion`, `criterion_type`, `max_iter`, and `abs_eigen`; and
#'  the maximum-likelihood `start_method`. Defaults to `estimate_control()`, which resolves
#'  every preset-driven knob from the `"EFAtools"` type. See [estimate_control()] for the
#'  individual knobs and the *Using the type presets* section in Details for how the preset
#'  fills them.
#' @param rotate_control a control object from [rotate_control()] bundling the rotation
#'  tuning knobs: the `type` preset; Kaiser `normalize`; the convergence `precision`; the
#'  factor `order_type`; the varimax/promax settings `varimax_type` and `p_type`; the
#'  simplimax/promax `k`; and `random_starts`. Defaults to `rotate_control()`, which resolves
#'  every preset-driven knob from the `"EFAtools"` type. The estimation and rotation presets
#'  are independent, so `estimate_control()` and `rotate_control()` may carry different
#'  `type`s. See [rotate_control()] for the individual knobs.
#' @param b_boot numeric. The number of bootstrap samples to draw. Default is 1000.
#'  Must be at least 2, the smallest number from which a standard error is defined.
#'  Under `cor_method = "fiml"` each bootstrap sample re-runs the EM moment
#'  estimation, so a smaller value may be advisable.
#' @param ci numeric. The confidence interval to create from the bootstrap samples.
#'  Must be between 0 and 1. Default is .95 for 95% CIs.
#' @param seed numeric. An optional seed for the random-number generator, governing every
#'  stochastic part of the fit: the rotation's random starts on the point estimate (the
#'  criterion-based rotations draw `random_starts` random starts; see *Rotations*) and,
#'  under `se = "np-boot"`, the case resampling, the replicate rotations, and the
#'  Procrustes random starts. Setting it makes the fit reproducible and the bootstrap
#'  additionally independent of the number of parallel workers (see Details); the
#'  caller's random-number stream is restored afterwards, so supplying a seed leaves no
#'  lasting effect on it. Default is `NULL`, which uses (and advances) the current state
#'  of the generator.
#' @param ... Additional arguments forwarded to the rotation engine. Only the arguments the
#'  selected `rotation` consumes are accepted: `maxit` (the maximum number of engine
#'  iterations) for the GPArotation-style rotations, plus the selected criterion's parameter
#'  (`gam` for "oblimin", where `gam = 0` is the recommended default and larger values can
#'  drive the solution toward factor collapse; `delta` for "geominT" and "geominQ").
#'  Anything else -- a misspelled name, another criterion's parameter, or any extra with
#'  "varimax", "promax", or "none", which consume no extras -- is an error rather than a
#'  setting that is silently ignored.
#'  The accepted arguments are merged with, and take precedence over, the extra arguments
#'  stored in `rotate_control()`. An estimation or rotation tuning knob (such as `type`,
#'  `max_iter`, or `k`) is likewise *not* accepted here: it belongs to [estimate_control()]
#'  or [rotate_control()].
#'
#' @details There are two main ways to use this function. The easiest way is to
#' use it with a specified `type` (see above), which sets most of the other
#' arguments accordingly. Another way is to use it more flexibly by explicitly
#' specifying all arguments used and set `type` to "none" (see examples).
#' A mix of the two can also be done by specifying a `type` as well as
#' additional arguments. However, this will throw warnings to avoid unintentional
#' deviations from the implementations according to the specified `type`.
#'
#' ## Estimators
#'
#' The estimator is chosen with `estimator`.
#'
#' - **PAF** (principal axis factoring) iteratively estimates the communalities and makes
#'   no distributional assumptions, which makes it robust and a good general-purpose
#'   default. Because it minimises no likelihood or weighted discrepancy it provides no
#'   model chi-square, and hence no chi-square-based fit indices (see *Fit indices*). The
#'   PAF iteration is governed by `init_comm`, `criterion`, `criterion_type`, `max_iter`,
#'   and `abs_eigen` (set by `type`; see *Using the type presets*).
#' - **ML** (maximum likelihood) maximises the normal-theory likelihood. It yields the
#'   full set of fit indices and is the only estimator with analytic expected-information
#'   standard errors (`se = "information"`), but it assumes multivariate normality and is
#'   the most prone to Heywood (improper) cases. Its starting values are set by
#'   `start_method`.
#' - **ULS** (unweighted least squares) minimises the sum of squared correlation
#'   residuals. "MINRES" (minimum residual) is the same estimator under a different name
#'   and returns identical results. It makes no normality assumption, is robust to mild
#'   non-normality, and yields the full set of fit indices.
#' - **DWLS** (diagonally weighted least squares) is the recommended estimator for ordinal
#'   data. It weights each off-diagonal correlation residual by the inverse asymptotic
#'   variance of the corresponding polychoric correlation (Muthén, 1984), reproducing the
#'   loadings of a diagonally weighted least squares fit (e.g.
#'   `lavaan::efa(..., estimator = "DWLS")`). It therefore requires raw ordinal data with
#'   `cor_method = "poly"` or `"tetra"` and has no fallback for a supplied correlation
#'   matrix or a continuous `cor_method`. Because the weighting follows the polychoric
#'   asymptotic covariance, the matrix and the weights are estimated on the
#'   listwise-complete cases. Its fit-index behaviour is described under *Fit indices*.
#'
#' ## Correlation methods
#'
#' When raw data are supplied, `cor_method` selects how the correlation matrix is computed
#' (it is ignored when a correlation matrix is entered directly).
#'
#' - **"pearson"** (default), **"spearman"**, and **"kendall"** are passed to
#'   [stats::cor()] for continuous or rank data.
#' - **"poly"** / **"tetra"** compute polychoric / tetrachoric correlations for ordinal /
#'   binary data, assuming an underlying bivariate-normal latent variable. They use a
#'   two-step estimator, matching `polycor::polychor()`. Every two-by-two response table
#'   with an empty cell is reproduced exactly by a correlation of 1 (or -1), so such a pair
#'   would otherwise be estimated at that boundary whatever the underlying correlation; a
#'   continuity correction of 0.5 is therefore added to the empty cell of a binary pair,
#'   preserving the table margins, as `lavaan` and `psych` do by default (Savalei, 2011).
#'   Larger tables get no correction: a response table showing a perfect ordering is instead
#'   reported at the boundary value, with a warning naming the pairs. Neither kind of pair
#'   has an asymptotic variance, so both are reported as `NA` and the `"DWLS"` estimator
#'   refuses data containing them. The polychoric asymptotic covariance that
#'   underlies both the DWLS weights and the scaled (sandwich) statistic relies on
#'   large-sample theory that degrades for empty or near-empty response-category
#'   combinations; with very sparse cells the resulting weights and standard errors can be
#'   unreliable (a warning is issued when empty cells are present), so interpret them with
#'   caution and consider collapsing rare categories. The response-category probabilities
#'   are integrated by Gauss-Legendre quadrature, and a strongly correlated pair -- at a
#'   lower correlation when its table has an empty response-category combination -- is
#'   re-integrated with a finer rule, because the base rule under-resolves the narrow
#'   conditional transition of such a pair. The resulting correlations agree with an exact
#'   bivariate-normal integrator to about `1e-5`, so no quadrature setting is exposed as an
#'   argument.
#' - **"fiml"** estimates a two-stage full-information maximum-likelihood correlation. The
#'   saturated multivariate-normal mean and covariance are estimated from raw data with
#'   missing values by an EM algorithm assuming the data are missing at random (Yuan,
#'   Marshall, & Bentler, 2002; Little & Rubin, 2002), and the standardized covariance is
#'   then analysed. This reproduces `psych::corFiml()` followed by `psych::fa()` and
#'   `lavaan(missing = "two.stage")`, *not* `lavaan::efa(missing = "ml")`, so the point
#'   estimates are not expected to match the latter. The model fit indices are corrected
#'   two-stage statistics (see *Fit indices*). `"fiml"` uses every case and handles the
#'   missingness itself, so `use` is ignored; it supplies a continuous (Pearson-type)
#'   correlation only and is therefore not compatible with `estimator = "DWLS"`. Standard
#'   errors are available analytically for `estimator = "ML"` or `"ULS"` and, for any estimator,
#'   by the non-parametric bootstrap (see *Standard errors*). For multiply imputed data,
#'   [efa_mi()] is the alternative route to handling missingness. Both routes assume
#'   the values are missing at random (MAR), and which one to prefer is largely
#'   practical: FIML is a single, efficient fit and is the simpler default when the
#'   analysis model is the whole story, whereas multiple imputation is more flexible when
#'   the imputation model should draw on auxiliary variables not in the factor model, or
#'   when the same imputations feed several downstream analyses.
#'
#' ## Rotations
#'
#' A rotation transforms the unrotated loadings toward a simpler, more interpretable
#' pattern; all rotations are performed by rotation engines built into the package.
#' Orthogonal rotations keep the factors uncorrelated, whereas oblique rotations let them
#' correlate (returning a pattern matrix, a structure matrix, and the factor
#' intercorrelations `Phi`) and are usually more realistic for psychological constructs.
#'
#' Orthogonal rotations:
#' - **varimax** maximises the variance of the squared loadings within each factor (column
#'   simplicity). It is the most widely used orthogonal rotation and spreads variance
#'   across factors rather than concentrating it in a general factor.
#' - **quartimax** simplifies the variables (rows) so that each loads mainly on one factor;
#'   it tends to produce a strong general factor.
#' - **equamax** is a Crawford-Ferguson compromise between varimax (column) and quartimax
#'   (row) simplicity.
#' - **geominT** uses a geometric-mean criterion that rewards a sparse pattern and
#'   tolerates variables with cross-loadings; a smaller offset `delta` gives a sparser
#'   solution but sharper local minima.
#' - **bentlerT** uses Bentler's invariant pattern simplicity criterion.
#' - **bifactorT** is the Jennrich-Bentler orthogonal bifactor criterion: a general factor
#'   plus group factors (bifactor simple structure).
#'
#' Oblique rotations:
#' - **promax** is a fast two-step rotation: a varimax solution is raised to a power
#'   (controlled by `k` and `p_type`) to form a target that is then fitted obliquely. It is
#'   the common, inexpensive oblique default.
#' - **oblimin** is a flexible oblique family controlled by `gam` (default 0); a good
#'   general-purpose criterion. `gam = 0` (quartimin) is the recommended setting: larger
#'   values increasingly reward correlated factors and can drive the solution toward
#'   factor collapse, so inspect `Phi` before interpreting a fit with `gam > 0`.
#' - **quartimin** is oblimin pinned at `gam = 0`; a robust default oblique criterion.
#' - **simplimax** drives the `k` smallest loadings toward zero. Its criterion is only
#'   piecewise smooth, so it is the most prone to local minima and relies on several
#'   `random_starts`.
#' - **bentlerQ** is the oblique Bentler invariant pattern simplicity criterion.
#' - **geominQ** is the oblique geomin criterion; it handles complex (cross-loading)
#'   structure well but is multimodal, so it benefits from more `random_starts` (and uses a
#'   more thorough multi-start search internally).
#' - **bifactorQ** is the oblique (correlated) Jennrich-Bentler bifactor criterion.
#'
#' The criterion-based rotations (all except varimax and promax) are fitted by gradient
#' projection with `random_starts` random starts to guard against local minima; the
#' complexity criteria (simplimax and geominQ in particular) are the most multimodal. The
#' starts are drawn from the random-number generator, so different starts can reach
#' genuinely different optima and such a fit is reproducible only when the generator is
#' controlled: pass `seed`, or call [set.seed()] beforehand. The
#' `type` argument changes the varimax and promax settings (see *Using the type presets*)
#' and, for every rotation, the factor `order_type`. A single factor cannot be rotated.
#'
#' ## Standard errors
#'
#' `se` selects whether and how standard errors (and matching confidence intervals) are
#' computed. Which quantities they cover depends on the method. The analytic methods
#' (`"information"` and `"sandwich"`) cover the unrotated loadings, the uniquenesses and
#' the communalities and, when a rotation is applied, the rotated loadings and -- for
#' oblique rotations -- the factor correlations and the structure coefficients. The
#' bootstrap (`"np-boot"`) covers the unrotated loadings, the residuals, and the fit
#' indices and, when a rotation is applied, the rotated loadings and -- for oblique
#' rotations -- the factor correlations and the structure coefficients; it reports no
#' uniqueness or communality standard errors (see the `SE` and `CI` entries in Value).
#'
#' - **"none"** (default) computes no standard errors.
#' - **"information"** returns analytic standard errors from the expected (Fisher)
#'   information matrix of the maximum-likelihood solution, and therefore requires
#'   `estimator = "ML"` and `cor_method = "pearson"` (or `"fiml"`, see below). The rotated
#'   standard errors are obtained by propagating the unrotated-loading covariance through
#'   the rotation by the delta method (Jennrich, 1973); because rotated quantities are
#'   identification-invariant they are directly comparable across programs. Unlike the
#'   bootstrap it also works from a correlation matrix as long as `N` is supplied.
#'
#'   [efa_fit()] analyses a *correlation* structure -- the diagonal of the analysed matrix
#'   is fixed at 1 and the uniquenesses \eqn{\psi_i = 1 - \sum_j \lambda_{ij}^2} are a
#'   function of the loadings rather than free parameters -- so the information is the
#'   correlation-structure one, \eqn{\Delta' \Gamma^{-1} \Delta} over the off-diagonal
#'   correlations, with \eqn{\Gamma} the normal-theory asymptotic covariance of the sample
#'   correlations (Olkin & Siotani, 1976) at the model-implied \eqn{\Sigma}
#'   (Cudeck, 1989). It is scaled by \eqn{1 / (N - 1)} and
#'   the confidence intervals are Wald intervals (estimate \eqn{\pm} z * SE). These
#'   standard errors assume multivariate normality and a correctly specified model; under
#'   heavy-tailed data or model misfit they can understate the sampling variability, where
#'   `"sandwich"` or `"np-boot"` are more robust.
#'
#'   The *rotated* quantities, the uniquenesses, and the communalities are
#'   identification-invariant and so are comparable across programs. The **unrotated**
#'   loading standard errors are not: they are reported in the orientation the solution
#'   itself is identified in (for ML, \eqn{\Lambda' \Psi^{-1} \Lambda} diagonal), and a
#'   program using a different identification will report different unrotated loading
#'   standard errors for the same fit.
#'
#'   That orientation can also fail to be determined. When two of the canonical variances
#'   \eqn{\mathrm{diag}(\Lambda' \Psi^{-1} \Lambda)} nearly coincide -- which happens with
#'   a weakly determined factor, or two factors of near-equal strength -- the loadings can
#'   be rotated within that plane without leaving the identifying constraint, so the
#'   unrotated loadings have no well-defined sampling distribution and their standard
#'   errors diverge. `efa_fit()` detects this and returns `NA` for the unrotated loading
#'   standard errors with an `efa_se_unreliable` warning, rather than reporting a number
#'   that looks like a standard error but is an artefact of the orientation. Everything
#'   that does not depend on the orientation -- the rotated loadings, `Phi`, the structure
#'   coefficients, the uniquenesses and the communalities -- is unaffected and still
#'   reported. Use those, or `se = "np-boot"`, when the unrotated loadings themselves are
#'   the quantity of interest. The detection covers the Pearson and polychoric paths; the
#'   two-stage `cor_method = "fiml"` sandwich carries no such diagnostic, so a weakly
#'   determined orientation is not flagged there. A Heywood case (a uniqueness at its lower
#'   boundary) is separate: the Wald approximation fails there for every parameter, so
#'   neither analytic method reports a standard error at all -- the whole `SE`/`CI` block is
#'   `NA` with an `efa_se_unreliable` warning. The `"sandwich"` scaled chi-square is not a
#'   Wald quantity and is still reported, so the fit indices are unaffected.
#' - **"sandwich"** returns robust (Godambe sandwich) standard errors from raw data,
#'   combining the estimator weight with an asymptotic-distribution-free covariance of the
#'   correlations, so it stays valid under non-normality and weight misspecification
#'   (Browne, 1984; Satorra & Bentler, 1994). It is available either for ordinal data with
#'   `cor_method = "poly"` or `"tetra"` and `estimator` one of `"ML"`, `"ULS"`, or `"DWLS"`
#'   (the meat is the polychoric / tetrachoric asymptotic covariance), or for continuous
#'   data with `cor_method = "pearson"` and `estimator = "ML"` or `"ULS"` (the meat is the
#'   fourth-moment ADF covariance of the sample correlations, the basis of the MLM / MLR
#'   robust statistics). It reports the same coverage as `"information"`, propagated
#'   by the same delta method, and additionally fills the model fit's chi-square block with
#'   a scaled (Satorra-Bentler / scaled-and-shifted) chi-square (see *Fit indices*). The
#'   statistic reported as `chi` is always the scaled-and-shifted one (the WLSMV default,
#'   flagged by `chi_scaled_type`), for the continuous Pearson path as well as the ordinal
#'   one; the mean-adjusted statistic that \pkg{lavaan}'s `MLM` reports for continuous data
#'   is returned alongside it as `chi_mean_adjusted`.
#'   Because the asymptotic covariance must describe the same cases as the correlation
#'   matrix, the sandwich (like `estimator = "DWLS"`) is computed on the listwise-complete
#'   cases; on data with missing values the reported `N`, the correlation matrix, and the
#'   point estimate therefore reflect the complete cases regardless of `use`.
#' - **"np-boot"** draws a non-parametric (case-resampling) bootstrap and needs raw data.
#'   A correlation matrix carries no cases to resample; alone among the unsupported
#'   combinations this one does not error but warns (condition class `"efa_boot_cormat"`)
#'   and downgrades `se` to `"none"`, so the fit is returned without an `SE` slot.
#'   It is the most general method -- available for any `estimator`, `rotation`, and
#'   `cor_method` -- and the most robust to non-normality and misfit, at the cost of speed;
#'   its intervals are bootstrap percentile intervals. The replicate fits are run across
#'   replicates with the `future` framework. By default they run sequentially; to run them
#'   in parallel, register a plan with [future::plan()] (e.g.
#'   `future::plan(future::multisession, workers = 2)`; see examples). With a fixed `seed`
#'   the bootstrap is reproducible and yields the same result regardless of the number of
#'   workers. Under `cor_method = "fiml"` each resample also
#'   re-runs the EM moment estimation and is therefore slow, so a smaller `b_boot` may be
#'   advisable.
#'
#'   The percentile intervals are centred on the point estimate for the loadings, the factor
#'   correlations, the structure coefficients and the residuals, but **not** for the indices
#'   derived from the chi-square (`RMSEA`, `AIC`, `BIC` and `ECVI`). A resample is drawn from
#'   the sample, which already carries the model's misfit, and is then refitted and evaluated
#'   against itself, so the replicate chi-square is on average about `chi + df` rather than
#'   `chi` and the whole interval rides upward with it -- often far enough that the point
#'   estimate falls below its own lower bound. That is the shift, not a miscomputed interval:
#'   read those intervals as a spread rather than as a range for the point estimate.
#'   Correcting the location needs resampling from a population transformed to fit the model
#'   (Bollen & Stine, 1992), which is not what this bootstrap does. `CFI` and `TLI` are not
#'   affected, being ratios in which the baseline chi-square shifts along with the model one.
#'
#' The analytic methods (`"information"` and `"sandwich"`) are not available with the
#' `"promax"` or `"simplimax"` rotations, which have no usable analytic rotation Jacobian;
#' use `"np-boot"` there. Under `cor_method = "fiml"`, `"information"` and `"sandwich"`
#' instead return, for `estimator = "ML"` or `"ULS"`, the corrected two-stage (Yuan & Bentler,
#' 2000; Savalei & Bentler, 2009) sandwich standard errors, built on the saturated FIML
#' asymptotic covariance with the estimator's own Stage-2 weight: the model is fitted to
#' the EM-estimated correlation, so the naive Stage-2 standard errors (treating that
#' correlation as complete data) are inconsistent under missingness and are not reported
#' (`estimator = "PAF"` carries no Stage-2 weight, so use `se = "np-boot"` there).
#'
#' ## Fit indices
#'
#' For ML and ULS, [efa_fit()] returns the model chi-square (with its p-value and degrees of
#' freedom), the Comparative Fit Index (CFI; Bentler, 1990), the Tucker-Lewis Index (TLI,
#' also called the non-normed fit index; Tucker & Lewis, 1973), the Root Mean Square Error
#' of Approximation (RMSEA) with its 90% confidence interval (Browne & Cudeck, 1992), the
#' Akaike and Bayesian Information Criteria (AIC, BIC), the Expected Cross-Validation Index
#' (ECVI; Browne & Cudeck, 1989), the Root Mean Squared Residual (RMSR), the Standardized
#' Root Mean Squared Residual (SRMR; Bentler, 1995), and the common-part-accounted-for
#' (CAF) index (Lorenzo-Seva, Timmerman, & Kiers, 2011). The print and summary methods show
#' SRMR, not RMSR, because the two residual summaries differ only by the fixed scaling
#' \eqn{\sqrt{(p - 1) / (p + 1)}} for a fixed number of variables; RMSR remains in the
#' returned object. The model chi-square is the
#' Bartlett-corrected discrepancy (matching [stats::factanal()] for ML); the AIC, BIC, and
#' ECVI are the minimum-fit-function (chi-square-based) forms (\eqn{\chi^2 - 2\,df} and
#' \eqn{\chi^2 - \log(N)\,df} for AIC and BIC, as in [psych::fa()]) and can therefore be
#' negative. The RMSEA, CFI, and TLI place the model and baseline
#' noncentrality on the uncorrected \eqn{N - 1} discrepancy scale on which these
#' approximate-fit indices are defined, so the Bartlett small-sample correction enters only
#' the chi-square test, not the approximate-fit indices.
#'
#' Which indices are reported depends on the estimator:
#' - **ML and ULS** compute the full set above.
#' - **PAF** returns only the descriptive residual indices (RMSR, SRMR, CAF) and df; the
#'   printed model-fit block shows CAF and SRMR. The chi-square-based indices are `NA`,
#'   because PAF minimises no discrepancy.
#' - **DWLS** by default returns only RMSR, SRMR, CAF, and df, because the ordinary
#'   maximum-likelihood discrepancy is not its fit function. When `se = "sandwich"`, a
#'   scaled (Satorra & Bentler, 1994; Asparouhov & Muthén, 2010) chi-square and the CFI,
#'   TLI, and RMSEA derived from it are reported (AIC and BIC remain `NA`). That scaled
#'   statistic is a two-stage correction applied to the polychoric-correlation residuals
#'   (Browne, 1984), so it is not identical to the full WLSMV test of \pkg{lavaan} or
#'   Mplus, which also projects the thresholds.
#' - **`cor_method = "fiml"`** (with ML or ULS) reports Satorra-Bentler-corrected two-stage
#'   statistics (Yuan, Marshall, & Bentler, 2002): the normal-theory discrepancy on the
#'   EM-estimated correlation, rescaled by the saturated FIML asymptotic covariance,
#'   because the plain two-stage likelihood-ratio statistic is not asymptotically
#'   \eqn{\chi^2(df)}. The CFI, TLI, and RMSEA follow from the scaled statistics; AIC, BIC,
#'   and ECVI are left `NA`, as for any scaled (moment-adjusted) chi-square.
#'
#' Whenever the chi-square is a scaled one (`se = "sandwich"`, or any `cor_method = "fiml"`
#' fit), the AIC, BIC, and ECVI are `NA` and the returned `fit_indices` additionally carry
#' the scaled-statistic components (see the `fit_indices` entry in Value). Note that
#' Lorenzo-Seva, Timmerman, and Kiers (2011) introduce the CAF as ranging between 0 and 1,
#' with values close to 1 indicating close fit; this does not match the formula they apply,
#' \eqn{1 - KMO(residuals)}, which only works if the diagonal of the residual
#' matrix is set to 1s and then approximates 0.5 with close fit.
#'
#' ## Available combinations
#'
#' Not every estimator, rotation, standard-error, and correlation method can be combined:
#'
#' - **Estimator and correlation method.** `estimator = "DWLS"` requires ordinal data with
#'   `cor_method = "poly"` or `"tetra"`. `cor_method = "fiml"` works with PAF, ML, and ULS
#'   (not DWLS) and needs raw data with missing values.
#' - **Standard errors.** `se = "information"` requires `estimator = "ML"` and
#'   `cor_method = "pearson"` or `"fiml"`, and can be computed from a correlation matrix
#'   when `N` is supplied. `se = "sandwich"` requires raw data,
#'   with either a polychoric/tetrachoric `cor_method` (ML, ULS, or DWLS) or a Pearson
#'   `cor_method` (ML or ULS); it is not available for PAF. Under `cor_method = "fiml"`,
#'   `"information"` and `"sandwich"` are available for ML and ULS only and both return the
#'   corrected two-stage sandwich. `se = "np-boot"` requires raw data and works with any
#'   estimator, rotation, and correlation method. Neither `"information"` nor `"sandwich"`
#'   is available with the `"promax"` or `"simplimax"` rotations.
#' - **Fit indices.** The chi-square-based indices are available for ML and ULS (and, as
#'   scaled statistics, for `cor_method = "fiml"` and for DWLS with `se = "sandwich"`); PAF
#'   and DWLS otherwise report only the descriptive residual indices.
#'
#' ## Using the type presets
#'
#' The `type` argument is evaluated for PAF and for all rotations (mainly
#' important for the varimax and promax rotations). The type-specific settings
#' for these functions are detailed below.
#'
#' For PAF, the values of `init_comm`, `criterion`, `criterion_type`,
#' `max_iter`, and `abs_eigen` depend on the `type` argument.
#'
#' `type = "EFAtools"` will use the following argument specification:
#' `init_comm = "smc", criterion = .001, criterion_type = "sum",
#' max_iter = 300, abs_eigen = TRUE`.
#'
#' `type = "psych"` will use the following argument specification:
#' `init_comm = "smc", criterion = .001, criterion_type = "sum",
#' max_iter = 50, abs_eigen = FALSE`.
#'
#' `type = "SPSS"` will use the following argument specification:
#' `init_comm = "smc", criterion = .001, criterion_type = "max_individual",
#' max_iter = 25, abs_eigen = TRUE`.
#'
#' If SMCs fail, SPSS takes "mac". However, as SPSS takes absolute eigenvalues,
#' this is hardly ever the case. Psych, on the other hand, takes "unity" if SMCs
#' fail, but uses the Moore-Penrose Psudo Inverse of a matrix, thus, taking "unity"
#' is only necessary if negative eigenvalues occur afterwards in the iterative
#' PAF procedure. The EFAtools type setting combination was the best in terms of accuracy
#' and number of Heywood cases compared to all the
#' other setting combinations tested in simulation studies in Grieder & Steiner
#' (2022), which is why this type is used as a default here.
#'
#' For varimax, the values of `varimax_type` and `order_type` depend on
#' the `type` argument.
#'
#' `type = "EFAtools"` will use the following argument specification:
#' `varimax_type = "kaiser", order_type = "eigen"`.
#'
#' `type = "psych"` will use the following argument specification:
#' `varimax_type = "svd", order_type = "eigen"`.
#'
#' `type = "SPSS"` will use the following argument specification:
#' `varimax_type = "kaiser", order_type = "ss_factors"`.
#'
#' For promax, the values of `p_type`,
#' `order_type`, and `k` depend on the `type` argument.
#'
#' `type = "EFAtools"` will use the following argument specification:
#' `p_type = "norm", order_type = "eigen", k = 4`.
#'
#' `type = "psych"` will use the following argument specification:
#' `p_type = "unnorm", order_type = "eigen", k = 4`.
#'
#' `type = "SPSS"` will use the following argument specification:
#' `p_type = "norm", order_type = "ss_factors", k = 4`.
#'
#' The `p_type` argument can take two values, "unnorm" and "norm". It controls
#' which formula is used to compute the target matrix P in the promax rotation.
#' "unnorm" uses the formula from Hendrickson and White (1964), specifically:
#' `P = abs(A^(k + 1)) / A`,
#' where A is the unnormalized matrix containing varimax rotated loadings.
#' "norm" uses the normalized varimax rotated loadings. Specifically it used the
#' following formula, which can be found in the SPSS 23 and SPSS 27 Algorithms manuals:
#' `P = abs(A / sqrt(rowSums(A^2))) ^(k + 1) * (sqrt(rowSums(A^2)) / A)`.
#' As for PAF, the EFAtools type setting combination for promax was the best
#' compared to the other setting combinations tested in simulation studies in
#' Grieder & Steiner (2022). Note that all `type` presets keep the EFAtools default
#' Kaiser normalization (`normalize = TRUE`), whereas [psych::fa()] does not
#' normalize before its promax target rotation; set `normalize = FALSE` to
#' reproduce the [psych::fa()] promax result to within the varimax convergence
#' tolerance (the residual difference is the convergence noise of the underlying
#' varimax base at `precision = 1e-5`, not an algorithmic difference).
#'
#' The `varimax_type` argument can take two values, "svd", and "kaiser". "svd" uses
#' singular value decomposition, by calling [stats::varimax()]. "kaiser"
#' performs the varimax procedure as described in the SPSS Algorithms manual and by
#' Kaiser (1958). The varimax simplicity criterion monitored for convergence is
#' `sum(n*colSums(lambda ^ 4) - colSums(lambda ^ 2) ^ 2) / n ^ 2`, where n is the
#' number of indicators, and lambda is the Kaiser-normalized rotated loadings matrix.
#'
#' For all other rotations except varimax and promax, the `type` argument
#' only controls the `order_type` argument with the same values as stated
#' above for the varimax and promax rotations. Additional arguments can also be
#' specified and will be passed to the rotation procedure (e.g., maxit to change the
#' maximum number of iterations). As for promax, every preset keeps the EFAtools
#' default Kaiser normalization (`normalize = TRUE`), whereas [psych::fa()] and
#' \pkg{GPArotation} do not normalize before these criterion rotations; set
#' `normalize = FALSE` to reproduce them.
#'
#' The `type` tuning arguments have no effect on ULS and ML; `type` itself still
#' governs the checks applied to the correlation matrix. For ULS, no additional
#' arguments are needed. For ML, an additional argument
#' `start_method` is needed to determine the starting values for the
#' optimization procedure. Default for this argument is "psych" which takes
#' the starting values specified in [psych::fa()].
#'
#'
#' @return A list of class `c("efa", "EFA")` containing (a subset of) the following:
#'
#' \item{orig_R}{Original correlation matrix.}
#' \item{h2_init}{Initial communality estimates from PAF.}
#' \item{h2}{Final communality estimates from the unrotated solution.}
#' \item{orig_eigen}{Eigen values of the original correlation matrix.}
#' \item{init_eigen}{Initial eigenvalues, obtained from the correlation matrix
#'  with the initial communality estimates as diagonal in PAF.}
#' \item{final_eigen}{Eigenvalues obtained from the correlation matrix
#'  with the final communality estimates as diagonal.}
#' \item{iter}{For PAF, the number of iterations until convergence. For ML, ULS, and
#'  DWLS, the number of objective-function evaluations used by the optimiser (not the
#'  number of optimiser iterations).}
#' \item{convergence}{Integer convergence code (0 = converged), using the codes of
#'  [`stats::optim()`][stats::optim]. For ML and ULS it is the code from the bounded
#'  optimiser; for DWLS the fit runs a bounded warm start followed by an
#'  unconstrained polish, and the reported code is from the final polish. For PAF
#'  it is 1 if the maximum number of iterations was reached
#'  without meeting the convergence criterion and 0 otherwise. A non-zero code is
#'  also reported with a warning.}
#' \item{heywood}{A named integer vector indicating which variables have a
#'  Heywood (improper) case in the unrotated solution; empty if there are none.}
#' \item{unrot_loadings}{Loading matrix containing the final unrotated loadings.}
#' \item{vars_accounted}{Matrix of explained variances and sums of squared loadings. Based on the unrotated loadings.}
#' \item{fit_indices}{A named list of fit indices computed from the unrotated
#' loadings. For ML and ULS it holds the model Chi Square (with its p-value and
#' df), CFI, TLI, RMSEA with its 90% confidence interval, AIC, BIC, ECVI, RMSR,
#' SRMR, and CAF; for PAF and DWLS only RMSR, SRMR, CAF, and df are populated and
#' the Chi-Square-based indices are `NA` (for DWLS with `se = "sandwich"` the full
#' block is filled from a scaled Chi Square instead). Whenever the Chi Square is a
#' scaled statistic (`se = "sandwich"`, or any `cor_method = "fiml"` fit) AIC,
#' BIC, and ECVI are `NA` and the list additionally carries the scaling
#' components: `chi_scaling` (the multiplier a in the scaled-and-shifted statistic
#' \eqn{aT + b}, i.e. the reciprocal of \pkg{lavaan}'s `chisq.scaling.factor`),
#' `chi_shift` (b), `chi_unscaled` (the unscaled statistic T), and the alternative
#' `chi_mean_adjusted` and `chi_mean_var` statistics with their `df_mean_var`.
#' RMSR is retained for programmatic use and backward compatibility, although the
#' print and summary methods display SRMR. See the *Fit indices* section in
#' Details for how each index is defined, scaled, and referenced.}
#' \item{model_implied_R}{The model implied correlation
#' matrix.}
#' \item{residuals}{Residual correlations, i.e., orig_R - model_implied_R}
#' \item{standardized_residuals}{Residual correlations standardized by their
#'  bootstrap standard errors. Only returned, if `se = "np-boot"`.}
#' \item{rot_loadings}{Loading matrix containing the final rotated loadings
#' (pattern matrix).}
#' \item{Phi}{The factor intercorrelations (only for oblique rotations).}
#' \item{Structure}{The structure matrix (only for oblique rotations).}
#' \item{rotmat}{The rotation matrix. The rotated loadings are recovered from the
#' unrotated loadings as `unrot_loadings %*% rotmat` for orthogonal rotations and
#' for promax, and as `unrot_loadings %*% t(solve(rotmat))` for the other oblique
#' rotations.}
#' \item{vars_accounted_rot}{Matrix of explained variances and sums of squared
#' loadings. Based on rotated loadings and, for oblique rotations, the factor
#' intercorrelations.}
#' \item{settings}{A list of the settings used. For the criterion rotations fitted by
#' gradient projection it additionally carries `rotation_diagnostics`, a list summarising
#' the multi-start run: `n_starts_total` (the `random_starts` random starts plus the
#' rational start), `n_optimized` (how many of those starts were actually optimized --
#' fewer than `n_starts_total` whenever the solver screens the random starts and optimizes
#' only the most promising ones), `n_converged` (how many optimized starts reached the
#' convergence tolerance), `n_distinct_minima` (how many distinct local optima those
#' converged starts found; more than one means the criterion is multimodal on these data),
#' `criterion_spread` (the range of the criterion values they attained), and
#' `criterion_best` (the criterion value of the returned solution). When
#' `normalize = TRUE`, `criterion_best` and `criterion_spread` are evaluated on the
#' Kaiser-normalized loadings the criterion is optimized on, not on the returned
#' `rot_loadings`, so they are not directly comparable to a criterion recomputed from the
#' returned loadings.}
#' \item{SE}{A named list of standard error matrices. For `se = "np-boot"`: bootstrap standard deviations of the unrotated and (when a rotation is applied) rotated loadings, the residuals, and the fit indices, plus -- for oblique rotations -- the factor correlations (`Phi`) and the structure coefficients; it additionally carries `valid_replicates`, the number of bootstrap replicates that were fitted and aligned successfully and that every entry above is therefore based on (replicates that failed are excluded and warned about), and, when a rotation is applied, `valid_target_rotations`, the number of those replicates that could also be aligned to the rotated point estimate and that the rotated entries (`rot_loadings`, `Phi`, `Structure`) are based on. For `se = "information"`: Wald standard errors from the expected (Fisher) information matrix for the unrotated loadings, the uniquenesses, and the communalities and, when a rotation is applied, the rotated loadings (and, for oblique rotations, `Phi` and the structure coefficients). Because \eqn{h^2_i = 1 - \psi_i} exactly, the communality and uniqueness standard errors are identical. For `se = "sandwich"`: robust Godambe sandwich standard errors with the same coverage as `"information"`, robust to non-normality and weight misspecification. Only returned if `se` is not `"none"`.}
#' \item{CI}{A named list of confidence intervals of width `ci`. For `se = "np-boot"`: percentile intervals matching the components of `SE`. For `se = "information"` and `se = "sandwich"`: Wald intervals matching the components of `SE`. Only returned if `se` is not `"none"`.}
#' \item{replicates}{A named list of bootstrap replicate arrays for the aligned unrotated and (where applicable) rotated loadings, structure coefficients, factor correlations (`Phi`), residuals, and fit indices. The replicate is the last dimension of the loading, residual, `Phi`, and structure cubes, and the first dimension of the `fit_indices` matrix (whose columns are named after the fit indices). Replicates that failed are left `NA`. Populated only for `se = "np-boot"`; `NULL` for the analytic SE methods.}
#' \item{vcov_unrot_loadings}{The full unrotated loading covariance matrix the marginal `SE$unrot_loadings` were derived from: a `p * n_factors` by `p * n_factors` numeric matrix in column-major `vec(Lambda)` order, with rows and columns labelled `"<variable>_<factor>"` so that ordering can be read off the object. Populated for `se = "information"` (expected-information block) and `se = "sandwich"` (robust V_AA), even when a rotation is applied (the persisted block is always the unrotated one); NA-filled if the analytic covariance is unreliable (a Heywood case or a singular bordered information matrix); `NULL` for `se = "np-boot"` and `se = "none"`. A weakly determined rotational orientation is the one case where this matrix is populated while `SE$unrot_loadings` is `NA`: the covariance itself is finite and valid, and only its gauge-dependent marginals are not (see *Standard errors*).}
#' \item{Gamma}{The asymptotic covariance of the off-diagonal sample correlations -- the meat of the robust sandwich SEs -- on the variance scale (`Var(rho-hat)`; lavaan's correlation NACOV is `N * Gamma`). A `p (p - 1) / 2` by `p (p - 1) / 2` numeric matrix; rows and columns ordered by [utils::combn()] over the column pairs and labelled `"<var_i>-<var_j>"`. Populated for `se = "sandwich"` on the polychoric/tetrachoric and Pearson paths; `NULL` otherwise, including under `cor_method = "fiml"`, whose meat is the saturated FIML asymptotic covariance and is not returned.}
#'
#' @source Bollen, K. A., & Stine, R. A. (1992). Bootstrapping goodness-of-fit measures
#' in structural equation models. Sociological Methods & Research, 21, 205–229.
#' doi: 10.1177/0049124192021002004
#' @source Grieder, S., & Steiner, M. D. (2022). Algorithmic jingle jungle: A comparison
#' of implementations of principal axis factoring and promax rotation in R and SPSS.
#' Behavior Research Methods, 54, 54–74. doi: 10.3758/s13428-021-01581-x
#' @source Hendrickson, A. E., & White, P. O. (1964). Promax: A quick method for
#' rotation to oblique simple structure. British Journal of Statistical Psychology,
#' 17 , 65–70. doi: 10.1111/j.2044-8317.1964.tb00244.x
#' @source Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. L. (2011). The
#' Hull Method for Selecting the Number of Common Factors, Multivariate Behavioral
#' Research, 46, 340-364, doi: 10.1080/00273171.2011.564527
#' @source Kaiser, H. F. (1958). The varimax criterion for analytic rotation in
#' factor analysis. Psychometrika, 23, 187–200. doi: 10.1007/BF02289233
#' @source Lawley, D. N., & Maxwell, A. E. (1971). Factor analysis as a statistical
#' method (2nd ed.). Butterworths.
#' @source Cudeck, R. (1989). Analysis of correlation matrices using covariance
#' structure models. Psychological Bulletin, 105, 317–327.
#' doi: 10.1037/0033-2909.105.2.317
#' @source Olkin, I., & Siotani, M. (1976). Asymptotic distribution of functions of a
#' correlation matrix. In S. Ikeda (Ed.), Essays in probability and statistics
#' (pp. 235–251). Shinko Tsusho.
#' @source Jennrich, R. I. (1973). Standard errors for obliquely rotated factor
#' loadings. Psychometrika, 38, 593–604. doi: 10.1007/BF02291497
#' @source Zhang, G., & Preacher, K. J. (2015). Factor rotation and standard errors
#' in exploratory factor analysis. Journal of Educational and Behavioral Statistics,
#' 40, 579–603. doi: 10.3102/1076998615606098
#' @source Browne, M. W. (1984). Asymptotically distribution-free methods for the
#' analysis of covariance structures. British Journal of Mathematical and Statistical
#' Psychology, 37, 62–83. doi: 10.1111/j.2044-8317.1984.tb00789.x
#' @source Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics and
#' standard errors in covariance structure analysis. In A. von Eye & C. C. Clogg (Eds.),
#' Latent variables analysis (pp. 399–419). Sage.
#' @source Asparouhov, T., & Muthén, B. (2010). Simple second order chi-square
#' correction. Mplus Technical Appendix.
#' @source Muthén, B. (1984). A general structural equation model with dichotomous,
#' ordered categorical, and continuous latent variable indicators. Psychometrika, 49,
#' 115–132. doi: 10.1007/BF02294210
#' @source Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods for mean
#' and covariance structure analysis with nonnormal missing data. Sociological Methodology,
#' 30, 165–200. doi: 10.1111/0081-1750.00078
#' @source Yuan, K.-H., Marshall, L. L., & Bentler, P. M. (2002). A unified approach to
#' exploratory factor analysis with missing data, nonnormal data, and in the presence of
#' outliers. Psychometrika, 67, 95–121. doi: 10.1007/BF02294711
#' @source Savalei, V., & Bentler, P. M. (2009). A two-stage approach to missing data:
#' Theory and application to auxiliary variables. Structural Equation Modeling, 16, 477–497.
#' doi: 10.1080/10705510903008238
#' @source Little, R. J. A., & Rubin, D. B. (2002). Statistical analysis with missing data
#' (2nd ed.). Wiley.
#' @source Bartlett, M. S. (1951). The effect of standardization on
#' approximation in factor analysis. Biometrika, 38, 337–344.
#' @source Bentler, P. M. (1990). Comparative fit indexes in structural models.
#' Psychological Bulletin, 107, 238–246. doi: 10.1037/0033-2909.107.2.238
#' @source Tucker, L. R., & Lewis, C. (1973). A reliability coefficient for maximum
#' likelihood factor analysis. Psychometrika, 38, 1–10. doi: 10.1007/BF02291170
#' @source Browne, M. W., & Cudeck, R. (1989). Single sample cross-validation indices
#' for covariance structures. Multivariate Behavioral Research, 24, 445–455.
#' doi: 10.1207/s15327906mbr2404_4
#' @source Browne, M. W., & Cudeck, R. (1992). Alternative ways of assessing model fit.
#' Sociological Methods & Research, 21, 230–258. doi: 10.1177/0049124192021002005
#' @source Bentler, P. M. (1995). EQS structural equations program manual. Multivariate
#' Software.
#'
#' @family factor analysis
#'
#' @seealso [estimate_control()] and [rotate_control()] for the estimation and rotation
#'  tuning knobs. [efa_retain()] for choosing `n_factors`, and [efa_scores()],
#'  [efa_reliability()], [efa_schmid_leiman()], and [efa_compare()] for working with the
#'  fitted solution.
#'
#' @export
#'
#' @examples
#'
#' # Principal axis factoring with oblimin rotation
#' mod_oblimin <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                        rotation = "oblimin")
#' mod_oblimin
#' summary(mod_oblimin)
#'
#' # ML estimation with oblimin rotation
#' mod_oblimin <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                        estimator = "ML", rotation = "oblimin")
#' mod_oblimin
#' summary(mod_oblimin)
#'
#' # Tuning knobs are supplied through the control objects. Here the SPSS preset is
#' # used for the estimation and rotation, with the maximum PAF iterations raised.
#' mod_spss <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                     rotation = "promax",
#'                     estimate_control = estimate_control(type = "SPSS", max_iter = 500),
#'                     rotate_control = rotate_control(type = "SPSS"))
#' mod_spss
#'
#' # Analytic (expected-information) standard errors for the above
#' ML_info <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                    estimator = "ML", rotation = "oblimin", se = "information")
#' ML_info
#' summary(ML_info)
#'
#' \donttest{
#' # Robust (sandwich) standard errors and a scaled chi-square for ordinal raw data.
#' # These need a polychoric/tetrachoric correlation method and estimator ML, ULS, or DWLS.
#' DWLS_rob <- efa_fit(DOSPERT_raw, n_factors = 6, cor_method = "poly",
#'                     estimator = "DWLS", rotation = "oblimin", se = "sandwich")
#' DWLS_rob
#' summary(DWLS_rob)
#'
#' # The same robust SEs and scaled chi-square for continuous data: a Pearson
#' # correlation with estimator ML or ULS (the fourth-moment ADF covariance).
#' ML_rob <- efa_fit(GRiPS_raw, n_factors = 1, cor_method = "pearson",
#'                   estimator = "ML", rotation = "none", se = "sandwich")
#' ML_rob
#' summary(ML_rob)
#' }
#'
#' \donttest{
#' # Two-stage FIML correlations from raw data with missing values: the saturated
#' # multivariate-normal moments are EM-estimated (assuming the data are missing at
#' # random) and the standardized covariance is analysed.
#' x_miss <- GRiPS_raw
#' x_miss[cbind(1:20, 1)] <- NA
#' efa_fiml <- efa_fit(x_miss, n_factors = 1, estimator = "ML", cor_method = "fiml")
#' efa_fiml
#' }
#'
#' \dontrun{
#' # Bootstrap standard errors from raw data, reproducible via a fixed seed and run
#' # in parallel across replicates.
#' future::plan(future::multisession, workers = 2)
#' efa_boot <- efa_fit(GRiPS_raw, n_factors = 1, estimator = "PAF", rotation = "none",
#'                     se = "np-boot", b_boot = 1000, seed = 42)
#' future::plan(future::sequential)
#' }
#'
efa_fit <- function(x, n_factors, N = NA,
                    estimator = c("PAF", "ML", "ULS", "MINRES", "DWLS"),
                    rotation = c("none", "varimax", "equamax", "quartimax", "geominT",
                                 "bentlerT", "bifactorT", "promax", "oblimin",
                                 "quartimin", "simplimax", "bentlerQ", "geominQ",
                                 "bifactorQ"),
                    se = c("none", "information", "sandwich", "np-boot"),
                    cor_method = c("pearson", "spearman", "kendall", "poly", "tetra",
                                   "fiml"),
                    use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                            "everything", "na.or.complete"),
                    estimate_control = NULL,
                    rotate_control = NULL,
                    b_boot = 1000, ci = .95, seed = NULL, ...) {

  # The control arguments share their names with the constructor functions, so a literal
  # `estimate_control = estimate_control()` default would recurse into the promise. Default
  # them to NULL and build the preset controls here instead; in call position R resolves the
  # namesake to the (function) constructor, skipping the NULL argument binding.
  if (is.null(estimate_control)) estimate_control <- estimate_control()
  if (is.null(rotate_control)) rotate_control <- rotate_control()

  # Perform argument checks
  .reject_flat_knobs(...names())
  .assert_cor_input(x)

  # n_factors has no default: catch its omission here so the most common first-call
  # mistake gets a classed condition pointing at the function that answers it,
  # rather than R's bare "argument is missing" error.
  if (missing(n_factors)) {
    cli::cli_abort(
      c("{.arg n_factors} is required.",
        "i" = "Use {.fn efa_retain} to decide how many factors to extract."),
      class = "efa_missing_n_factors"
    )
  }

  estimator <- .match_arg_ci(estimator)
  # "MINRES" is a synonym for "ULS" (same estimator); resolve it once here so the
  # rest of efa_fit() and the reported settings use the single canonical name.
  if (estimator == "MINRES") estimator <- "ULS"
  rotation <- .match_arg_ci(rotation)
  se <- match.arg(se)
  np_boot <- se == "np-boot"
  use <- match.arg(use)
  cor_method <- match.arg(cor_method)

  # The dots only carry extras for the rotation engine, which reads them by exact name:
  # `maxit` for the GPArotation-style engines, plus the selected criterion's parameter
  # (`gam` for oblimin, `delta` for geomin). Any other dot -- a misspelled name (e.g.
  # `gamma` for `gam`), another criterion's parameter, or any extra with varimax/promax,
  # which take none -- would be silently dropped and the fit would quietly run the engine
  # defaults instead, so it is rejected here. With no rotation there is no engine, so
  # nothing can consume the dots at all.
  if (...length() > 0L) {
    allowed <- .rotation_dot_extras(rotation)
    nms <- ...names()
    if (is.null(nms)) nms <- rep("", ...length())
    bad <- unique(setdiff(nms[nzchar(nms)], allowed))
    if (length(bad) > 0L || any(!nzchar(nms))) {
      if (rotation == "none") {
        msg <- if (length(bad) > 0L) {
          "{.arg {bad}} {?is/are} not used: with {.code rotation = \"none\"} no rotation engine
           runs, so {.arg ...} must be empty."
        } else {
          "{.arg ...} must be empty with {.code rotation = \"none\"}: no rotation engine runs,
           so the extra arguments are not used."
        }
        info <- c("i" = "Check for a misspelled argument name (for example {.arg estimate_control}
                         or {.arg rotate_control}).")
      } else {
        msg <- if (length(bad) > 0L) {
          "{.arg {bad}} {?is/are} not consumed by the {.val {rotation}} rotation."
        } else {
          "The arguments in {.arg ...} must be named; the {.val {rotation}} rotation engine
           reads its extras by exact name."
        }
        info <- c(
          if (length(allowed) > 0L) {
            c("i" = "The {.val {rotation}} rotation accepts {.arg {allowed}} through {.arg ...}.")
          } else {
            c("i" = "The {.val {rotation}} rotation takes no extra arguments through {.arg ...}.")
          },
          "i" = "Check for a misspelled engine argument; the rotation tuning knobs live in
                 {.fn rotate_control}."
        )
      }
      cli::cli_abort(c(msg, info), class = "efa_unused_dots")
    }
  }

  .assert_estimate_control(estimate_control)
  if (!inherits(rotate_control, "efa_rotate_control")) {
    cli::cli_abort(
      "{.arg rotate_control} must be a control object from {.fn rotate_control}.",
      class = "efa_control_input"
    )
  }

  # Unbundle the estimation and rotation tuning knobs from the two control objects (already
  # validated by their constructors) into the local variables the fitting pipeline expects.
  # Each control carries its own `type` preset, so estimation and rotation may resolve from
  # different presets. The extra rotation arguments the rotate control stored are merged with
  # `...` (with `...` winning on a name clash) and forwarded to the rotation engine.
  est_type <- estimate_control$type
  init_comm <- estimate_control$init_comm
  criterion <- estimate_control$criterion
  criterion_type <- estimate_control$criterion_type
  max_iter <- estimate_control$max_iter
  abs_eigen <- estimate_control$abs_eigen
  start_method <- estimate_control$start_method

  rot_type <- rotate_control$type
  normalize <- rotate_control$normalize
  precision <- rotate_control$precision
  order_type <- rotate_control$order_type
  varimax_type <- rotate_control$varimax_type
  p_type <- rotate_control$p_type
  k <- rotate_control$k
  random_starts <- rotate_control$random_starts
  rot_extra_args <- modifyList(rotate_control$extra_args, list(...))

  if (is.na(start_method) && estimator == "ML") {
    cli::cli_abort(
      c("{.arg start_method} must be set when {.code estimator = \"ML\"}.",
        "i" = "Set {.arg start_method} to {.val psych} or {.val factanal}."),
      class = "efa_ml_start_missing"
    )
  }

  # Detect a correlation-matrix input once, up front: the FIML guard below and the
  # bootstrap/DWLS guards further down all branch on it, and resolving it before R is
  # validated/smoothed keeps those warnings from being pre-empted by a singular or
  # non-positive-definite matrix.
  is_cormat <- .is_cormat(x)

  # Two-stage / full-information ML correlations are EM-estimated from raw data with missing
  # values (Yuan, Marshall, & Bentler, 2002). Reject the input and option combinations that path
  # cannot honour here, before any computation, so they fail with a dedicated message rather than
  # a downstream one. The analytic standard errors are handled by the corrected two-stage sandwich
  # (.se_fiml()) for ML and ULS; DWLS, a correlation-matrix input, and PAF with analytic standard
  # errors are unsupported.
  if (cor_method == "fiml") {
    if (is_cormat) {
      cli::cli_abort(
        c("{.code cor_method = \"fiml\"} needs raw data, not a correlation matrix.",
          "x" = "A correlation matrix carries no cases to estimate the FIML moments from.",
          "i" = "Supply raw data (with missing values), or choose another {.arg cor_method}."),
        class = "efa_fiml_needs_raw"
      )
    }
    if (estimator == "DWLS") {
      cli::cli_abort(
        c("{.code estimator = \"DWLS\"} is not compatible with {.code cor_method = \"fiml\"}.",
          "x" = "DWLS needs a polychoric asymptotic covariance, which the continuous FIML correlation does not provide.",
          "i" = "Use {.code estimator = \"ML\"}, {.val ULS}, or {.val PAF}, or {.code cor_method = \"poly\"}/{.val tetra} for DWLS."),
        class = "efa_fiml_unsupported_method"
      )
    }
    # The corrected two-stage sandwich (.se_fiml()) reuses the Stage-2 ML/ULS weight; PAF minimises
    # no discrepancy and so carries no weight for it, exactly as the polychoric/continuous sandwich
    # rejects PAF. Reject the analytic-SE request here (the bootstrap stays available for PAF).
    if (se %in% c("information", "sandwich") && !(estimator %in% c("ML", "ULS"))) {
      cli::cli_abort(
        c("Analytic standard errors under {.code cor_method = \"fiml\"} require {.code estimator = \"ML\"} or {.val ULS}.",
          "x" = "You requested {.code estimator = {.val {estimator}}}.",
          "i" = "{.val PAF} minimises no discrepancy, so it has no weight for the two-stage sandwich; use {.code se = \"np-boot\"} for {.val {estimator}}."),
        class = "efa_se_unsupported"
      )
    }
    if (use != "pairwise.complete.obs") {
      cli::cli_warn(
        c("{.arg use} is ignored when {.code cor_method = \"fiml\"}.",
          "i" = "FIML uses every case and handles the missingness itself, so {.code use = {.val {use}}} has no effect."),
        class = "efa_fiml_use_ignored"
      )
    }
  }

  # Analytic standard errors cover only a subset of estimators and rotations. Reject the
  # unsupported combinations here, before any computation, with a clear pointer to the
  # bootstrap. Information-matrix SEs are derived from the ML discrepancy and so require
  # estimator = "ML"; promax (a two-step target rotation) and simplimax (a non-smooth, piecewise
  # criterion) have no usable analytic rotation Jacobian; PAF has no discrepancy-based information
  # from which a sandwich could be built.
  if (se == "information" && estimator != "ML" && cor_method != "fiml") {
    cli::cli_abort(
      c("{.code se = \"information\"} is only available for {.code estimator = \"ML\"}.",
        "x" = "You requested {.code estimator = {.val {estimator}}}.",
        "i" = "Use {.code se = \"np-boot\"} for {.val {estimator}}."),
      class = "efa_se_unsupported"
    )
  }
  # The expected information is built from the normal-theory asymptotic covariance of PEARSON
  # correlations, so it presumes normal-theory sampling behaviour of the analysed matrix. A
  # polychoric/tetrachoric correlation additionally carries first-stage threshold estimation error
  # and a rank correlation is not a Pearson moment at all; neither is accounted for here. The
  # sandwich covers both, so point at it rather than at the bootstrap.
  if (se == "information" && !(cor_method %in% c("pearson", "fiml"))) {
    cli::cli_abort(
      c("{.code se = \"information\"} requires {.code cor_method = \"pearson\"} or {.val fiml}.",
        "x" = "You requested {.code cor_method = {.val {cor_method}}}.",
        "i" = "The expected information assumes normal-theory sampling behaviour of the analysed correlations.",
        # From a bare correlation matrix the sandwich is unavailable too (it needs the raw data its
        # meat is estimated from), so pointing there would be a dead end.
        "i" = if (is_cormat) {
          "No analytic standard error is available for a {.val {cor_method}} correlation matrix; supply the raw data and use {.code se = \"sandwich\"}."
        } else {
          "Use {.code se = \"sandwich\"} for {.val {cor_method}}."
        }),
      class = "efa_se_unsupported"
    )
  }
  if (se == "sandwich" && estimator == "PAF" && cor_method != "fiml") {
    cli::cli_abort(
      c("{.code se = \"sandwich\"} is not available for {.code estimator = \"PAF\"}.",
        "i" = "Use {.code se = \"np-boot\"} for {.val PAF}."),
      class = "efa_se_unsupported"
    )
  }
  if (se %in% c("information", "sandwich") && rotation %in% c("promax", "simplimax")) {
    cli::cli_abort(
      c("{.code se = {.val {se}}} is not available with {.code rotation = {.val {rotation}}}.",
        "i" = "Use {.code se = \"np-boot\"} for {rotation}-rotated solutions."),
      class = "efa_se_unsupported"
    )
  }
  # Sandwich (robust) SEs combine the estimator weight with an asymptotic-distribution-free
  # covariance of the correlations (the robust meat), so they require raw data (enforced once N
  # is resolved below) to estimate that covariance. Two paths supply it: the polychoric/
  # tetrachoric asymptotic covariance for ordinal data (any estimator), or the fourth-moment
  # (Browne, 1984) covariance for continuous data with cor_method = "pearson" and estimator ML
  # or ULS. Spearman/Kendall correlations and continuous DWLS have no such covariance.
  if (se == "sandwich" && cor_method != "fiml" && !.is_poly_cor(cor_method) &&
      !(cor_method == "pearson" && estimator %in% c("ML", "ULS"))) {
    cli::cli_abort(
      c("{.code se = \"sandwich\"} is not available for this correlation/estimator combination.",
        "x" = "You requested {.code cor_method = {.val {cor_method}}} with {.code estimator = {.val {estimator}}}.",
        "i" = "Use {.code cor_method = \"poly\"}/{.val tetra} (any estimator), or {.code cor_method = \"pearson\"} with {.code estimator = \"ML\"} or {.val ULS}; otherwise use {.code se = \"np-boot\"}."),
      class = "efa_se_unsupported"
    )
  }

  checkmate::assert_count(n_factors)
  checkmate::assert_count(N, na.ok = TRUE)
  checkmate::assert_count(max_iter, na.ok = TRUE)
  checkmate::assert_choice(init_comm, c("smc", "mac", "unity", NA))
  checkmate::assert_number(criterion, lower = 0, upper = 1, na.ok = TRUE)
  if (!is.na(criterion) && criterion >= 1) {
    cli::cli_abort(
      c("{.arg criterion} must be smaller than 1.",
        "x" = "You supplied {.arg criterion} = {criterion}.",
        "i" = "Use a small positive convergence tolerance such as {.val 0.001}."),
      class = "efa_criterion_too_large"
    )
  }
  checkmate::assert_choice(criterion_type, c("max_individual", "sum", NA))
  checkmate::assert_flag(abs_eigen, na.ok = TRUE)
  checkmate::assert_number(k, na.ok = TRUE)
  checkmate::assert_choice(varimax_type, c("svd", "kaiser", NA))
  checkmate::assert_flag(normalize, na.ok = TRUE)
  checkmate::assert_choice(p_type, c("unnorm", "norm", NA))
  checkmate::assert_number(precision, lower = 0, upper = 1)
  checkmate::assert_choice(order_type, c("eigen", "ss_factors", NA))
  checkmate::assert_integerish(b_boot, len = 1, any.missing = FALSE)
  # A bootstrap standard error is the dispersion across replicates, so two is the smallest number
  # from which one is defined at all: at b_boot = 1 every SE is the sd() of a single value and comes
  # back NA, and the percentile interval collapses onto that replicate.
  if (b_boot < 2) {
    cli::cli_abort(
      c("{.arg b_boot} must be at least 2.",
        "x" = "You supplied {.arg b_boot} = {b_boot}.",
        "i" = "A bootstrap standard error is the spread across replicates and is undefined below
               two of them. The default is {.val {1000}}."),
      class = "efa_b_boot_too_small"
    )
  }
  checkmate::assert_number(ci, lower = 0, upper = 1)
  checkmate::assert_int(seed, null.ok = TRUE)

  # The common-factor model needs at least one factor and fewer factors than
  # variables: with n_factors >= n_variables it is not identified and the
  # eigenvalue-based extraction in the ML, ULS, and PAF fitters reads past the
  # available eigenvalues (undefined behaviour in an unchecked build), and with
  # n_factors = 0 there is no model to fit. Both bounds are checked here so the
  # C++ kernels' own bound checks are never the user-facing error.
  n_vars <- ncol(x)
  # With a single variable no number of factors is admissible, so only state the
  # range when one exists (it would otherwise read "1 to 0").
  range_hint <- if (n_vars > 1) " The admissible range is 1 to {n_vars - 1}." else ""
  if (n_factors < 1) {
    cli::cli_abort(
      c("{.arg n_factors} must be at least 1.",
        "x" = paste0("You requested {n_factors} factor{?s}.", range_hint),
        "i" = "A factor retention criterion that returns 0 means no factor is worth
               extracting, so there is no exploratory factor model to fit."),
      class = "efa_too_few_factors"
    )
  }
  if (n_factors >= n_vars) {
    cli::cli_abort(
      c("{.arg n_factors} must be smaller than the number of variables.",
        "x" = paste0("You requested {n_factors} factor{?s} for {n_vars} variable{?s}.",
                     range_hint),
        "i" = "Extract fewer factors."),
      class = "efa_too_many_factors"
    )
  }

  # DWLS weights each polychoric correlation residual by the inverse of its asymptotic
  # variance, so it needs a polychoric/tetrachoric asymptotic covariance. That is only
  # available from raw ordinal data with cor_method = "poly" or "tetra"; there is no
  # fallback to unit weights. Resolve this before any computation so an unsupported
  # request fails with a single clear error rather than downstream.
  if (estimator == "DWLS" && !(.is_poly_cor(cor_method) && !is_cormat)) {
    cli::cli_abort(
      c("{.code estimator = \"DWLS\"} requires a polychoric asymptotic covariance.",
        "x" = if (is_cormat) {
          "You supplied a correlation matrix, so no asymptotic covariance can be estimated."
        } else {
          "{.code cor_method = {.val {cor_method}}} does not produce one."
        },
        "i" = "Supply raw ordinal data with {.code cor_method = \"poly\"} or {.code \"tetra\"}."),
      class = "efa_dwls_no_acov"
    )
  }

  # A correlation matrix carries no cases to resample, so the bootstrap is impossible.
  # Analytic SEs (information/sandwich) need only R, the loadings, the uniquenesses and N,
  # so they remain available from a correlation matrix as long as N is supplied (checked
  # below, once N is resolved).
  if (is_cormat && isTRUE(np_boot)) {

    cli::cli_warn(
      c("Cannot compute bootstrap standard errors from correlation matrix.",
      "x" = "You've supplied {.var se} = {.val {se}}, but {.var x} is a correlation matrix.",
      "i" = "Setting {.var se} to {.val none}. Rerun with raw data to calculate bootstrap SEs."),
      class = "efa_boot_cormat"
    )

    np_boot <- FALSE
    se <- "none"

  }

  # A fixed seed makes the whole fit reproducible: every stochastic step downstream of
  # this point draws from the state it sets. That covers the rotation's random starts on
  # the point estimate -- the criterion rotations draw `random_starts` random orthogonal
  # starts from the ambient stream, so without a seed a simplimax or geominQ fit is not
  # reproducible run to run -- and, when `se = "np-boot"`, the case resampling, the
  # replicate rotations, and the Procrustes random starts as well. The bootstrap is also
  # independent of the number of parallel workers: the case resampling advances the
  # global RNG by a b_boot-dependent amount, and the parallel replicate fit then adds a
  # fixed, worker-count-independent step when future.seed = TRUE derives a per-replicate
  # L'Ecuyer stream. Both advances are deterministic given the seed, so the downstream
  # draws -- and the result -- are identical at any number of workers. The caller's RNG
  # stream is saved and restored afterwards -- or, if none existed yet, the state
  # set.seed() creates is removed again -- so efa_fit() leaves no side effect on it.
  .set_local_seed(seed)

  # Detect or compute the correlation matrix, check it, and smooth it if needed
  prep <- .prepare_cor_input(x, N = N, use = use, cor_method = cor_method,
                             N_policy = "optional",
                             # Sandwich SEs need the full asymptotic covariance of the
                             # correlations (the robust meat); for DWLS its diagonal also
                             # supplies the per-element weights, so "full" subsumes the DWLS
                             # "diag" request. dwls = TRUE builds those weights; the ML / ULS
                             # sandwich path uses the meat alone. A correlation matrix carries no
                             # raw data to estimate the covariance from, so do not request one
                             # there (the sandwich is rejected just below); requesting it would
                             # only draw a spurious "acov ignored" warning before that abort.
                             acov = if (se == "sandwich" && !is_cormat && cor_method != "fiml") "full"
                                    else if (estimator == "DWLS") "diag" else "none",
                             dwls = estimator == "DWLS",
                             # The psych preset tolerates a singular correlation matrix only
                             # for PAF, whose starting communalities fall back to a
                             # pseudo-inverse. No other estimator has that fallback -- ML and
                             # ULS solve R for their starting values and the ML discrepancy
                             # needs log|R|, and DWLS would silently drop to a flat start --
                             # so the check stays on for them under every preset.
                             check_singular = !(est_type == "psych" && estimator == "PAF"),
                             posdef_abort = est_type == "SPSS")
  R <- prep$R
  N <- prep$N
  # DWLS weight matrix (1 / diag(Gamma)); NULL for the other estimators.
  weights <- prep$weights
  # Full asymptotic covariance of the off-diagonal correlations (the sandwich meat); NULL unless
  # se = "sandwich" on a raw-data path (polychoric/tetrachoric, or continuous Pearson).
  Gamma <- prep$Gamma
  # Two-stage / full-information ML moments for the fit-index likelihood-ratio chi-square; the
  # saturated mean/covariance and log-likelihood EM-estimated in .prepare_cor_input() are
  # carried alongside the raw data so .gof() can evaluate the model and baseline FIML log-
  # likelihoods. Reuse prep$fiml (no second EM run); NULL for every other cor_method.
  fiml_pt <- if (cor_method == "fiml") {
    list(data = x, mu = prep$fiml$mu, sigma = prep$fiml$sigma, logl = prep$fiml$logl)
  }
  # On an analytic-SE fit the saturated FIML correlation covariance is needed twice from the same
  # moments -- by the scaled chi-square (.fiml_scaled_test, via .gof()) and by the SE sandwich
  # (.se_fiml_core) -- so build it once here and cache it on the fiml list; both read it and skip the
  # rebuild. Only the analytic-SE path duplicates it: with se = "none"/"np-boot" only .gof() forms it
  # (once), and the bootstrap replicate lists omit the field and recompute per replicate. NULL on a
  # degenerate covariance, where each consumer falls back to its own guarded recompute.
  if (!is.null(fiml_pt) && se %in% c("information", "sandwich")) {
    fiml_pt$acov_cor <- tryCatch(
      .fiml_saturated_acov(fiml_pt$data, fiml_pt$mu, fiml_pt$sigma)$cor,
      error = function(e) NULL)
  }

  # Analytic SEs scale the inverse information by 1 / (N - 1), so they need the sample
  # size. Raw data always supplies it (N = number of rows); a correlation matrix does not,
  # so require it explicitly there.
  if (se %in% c("information", "sandwich") && is.na(N)) {
    cli::cli_abort(
      c("{.code se = {.val {se}}} requires the sample size {.arg N}.",
        "i" = "Supply {.arg N} when {.arg x} is a correlation matrix."),
      class = "efa_se_no_n"
    )
  }

  # The sandwich meat (the asymptotic covariance of the correlations) can only be estimated from
  # raw data; a correlation matrix carries no such covariance, so reject that combination
  # explicitly (rather than failing later when the meat is NULL).
  if (se == "sandwich" && is.null(Gamma) && cor_method != "fiml") {
    cli::cli_abort(
      c("{.code se = \"sandwich\"} requires raw data to estimate the asymptotic covariance of the correlations.",
        "x" = "You supplied a correlation matrix.",
        "i" = "Provide raw data (ordinal with {.code cor_method = \"poly\"}/{.val tetra}, or continuous with {.code cor_method = \"pearson\"}), or use {.code se = \"np-boot\"}."),
      class = "efa_se_unsupported"
    )
  }

  if (!is_cormat && isTRUE(np_boot)) {

    m <- ncol(R)
    # Resample the cases the correlation matrix was actually built from. Under
    # listwise deletion that is the complete cases (N of them), not the first N
    # row positions, so the case bootstrap stays a faithful resample of the
    # estimator that produced R (Efron & Tibshirani, 1993). FIML uses every row
    # that carries at least one observed value (fully-missing rows are dropped, so
    # N counts only these) and `use` does not apply; resampling exactly those rows
    # -- not all nrow(x) positions -- keeps each replicate's sample size equal to N.
    rows <- if (cor_method == "fiml") {
      which(rowSums(!is.na(x)) > 0L)
    } else if (.is_listwise_use(use) || estimator == "DWLS") {
      which(stats::complete.cases(x))
    } else {
      seq_len(nrow(x))
    }

    # create bootstrap samples and from these, correlation matrices
    R_boot_array <- array(NA_real_, c(m, m, b_boot), dimnames = list(colnames(x),
                                                               colnames(x),
                                                               NULL))

    poly_cor <- .is_poly_cor(cor_method)
    tetra_cor <- cor_method == "tetra"
    fiml_cor <- cor_method == "fiml"
    dwls <- estimator == "DWLS"

    # DWLS reweights each replicate by the inverse of its own polychoric asymptotic
    # variances, so the per-element weights are recomputed alongside the matrix and
    # carried into the lean fit; NULL for the other estimators, which need no weights.
    # The weights are positional, so no dimnames are needed.
    W_boot_array <- if (dwls) array(NA_real_, c(m, m, b_boot)) else NULL

    # FIML carries each replicate's resampled data and its own EM moments (saturated mean,
    # covariance, and log-likelihood) into the lean fit, so the per-replicate fit indices use
    # the same likelihood-ratio chi-square as the point estimate. One list element per
    # replicate (NULL for a dropped resample); NULL for the other cor_methods.
    fiml_boot <- if (fiml_cor) vector("list", b_boot) else NULL

    for (boot_i in seq_len(b_boot)) {
      ind <- sample(rows, size = N, replace = TRUE)
      if (poly_cor) {
        # A resample can be degenerate -- a constant column, a pair with no
        # overlapping cases, a numerically uncomputable matrix, or (for DWLS) a pair
        # whose asymptotic variance is non-positive -- and make .polychoric() or the
        # weight construction fail. Any genuine bug would already have surfaced on the
        # point-estimate fit over the full data above, so a failure here is necessarily
        # resample-specific; fall back to an all-NA matrix so the replicate is dropped at
        # the fit stage, mirroring how stats::cor() returns NA for a degenerate Pearson
        # resample and how .boot_fun() drops unfittable replicates. The bootstrap is
        # parallelised at the fit (across replicates, via future). DWLS requests the diagonal
        # ACOV and builds the weights inside the same try so the matrix and weights share one
        # resample and a degenerate weight drops the replicate too.
        rep_i <- tryCatch(
          suppressWarnings({
            poly <- .polychoric(x[ind, , drop = FALSE], nearest_pd = FALSE,
                                binary_only = tetra_cor,
                                acov = if (dwls) "diag" else "none",
                                label_acov = FALSE)
            list(R = poly$R,
                 W = if (dwls) .poly_weight_matrix(poly$acov, m) else NULL)
          }),
          error = function(e) NULL)
        if (is.null(rep_i)) {
          R_boot_array[,, boot_i] <- matrix(NA_real_, m, m)
        } else {
          R_boot_array[,, boot_i] <- rep_i$R
          if (dwls) W_boot_array[,, boot_i] <- rep_i$W
        }
      } else if (fiml_cor) {
        # Two-stage FIML: re-estimate the EM moments on the resample and standardise to a
        # correlation, mirroring the point-estimate path. The full EM object (moments and the
        # saturated log-likelihood) and the resampled data are retained so the lean fit can
        # form this replicate's FIML likelihood-ratio chi-square. A degenerate resample (a
        # constant or collinear column, an EM breakdown) makes .fiml_em_moments() abort; fall
        # back to an all-NA matrix so the replicate is dropped at the fit stage, as on the poly
        # and Pearson paths. The EM recompute is serial here; the replicate fits parallelise
        # downstream.
        em_i <- tryCatch(suppressWarnings(.fiml_em_moments(x[ind, , drop = FALSE])),
                         error = function(e) NULL)
        if (is.null(em_i)) {
          R_boot_array[,, boot_i] <- matrix(NA_real_, m, m)
        } else {
          R_boot_array[,, boot_i] <- stats::cov2cor(em_i$sigma)
          fiml_boot[[boot_i]] <- list(data = x[ind, , drop = FALSE], mu = em_i$mu,
                                      sigma = em_i$sigma, logl = em_i$logl)
        }
      } else {
        R_boot_array[,, boot_i] <- stats::cor(x[ind, , drop = FALSE], use = use,
                                              method = cor_method)
      }
    }

  }

  core_args <- list(
    R = R, N = N, weights = weights, Gamma = Gamma,
    R_boot_array = if (isTRUE(np_boot)) R_boot_array else NULL,
    W_boot_array = if (isTRUE(np_boot)) W_boot_array else NULL,
    fiml = fiml_pt, fiml_boot = if (isTRUE(np_boot)) fiml_boot else NULL,
    np_boot = np_boot, b_boot = b_boot, estimator = estimator, rotation = rotation,
    type = est_type, rot_type = rot_type, n_factors = n_factors, se = se, ci = ci,
    use = use, cor_method = cor_method, max_iter = max_iter, init_comm = init_comm,
    criterion = criterion, criterion_type = criterion_type,
    abs_eigen = abs_eigen, start_method = start_method, normalize = normalize,
    precision = precision, order_type = order_type, varimax_type = varimax_type,
    P_type = p_type, k = k, randomStarts = random_starts)

  # Splice the criterion-specific rotation extras (from `rotate_control()` or `...`) into
  # `.efa_core()`'s dots, exactly where the old flat `...` used to land. `.efa_core()` still
  # declares the former flat tuning knobs (`type`, `max_iter`, `k`, ...) as named formals, so
  # a caller that forwards one of them through `...` -- e.g. the retention criteria, efa_group,
  # or EFA_POOLED passing a former flat argument -- would otherwise supply that formal twice
  # and abort with "matched by multiple actual arguments". Drop any extra whose name is already
  # an explicit `core_args` entry (those knobs now live in the control objects, so a bare copy
  # is silently ignored); genuine rotation extras such as `maxit` or `gam` are not core_args and
  # pass through untouched.
  extra <- rot_extra_args[setdiff(names(rot_extra_args), names(core_args))]
  do.call(.efa_core, c(core_args, extra))
}

# Fit the common-factor model from already-prepared inputs: a correlation matrix R, the
# sample size N, optional DWLS weights, the optional sandwich meat Gamma, and -- for the
# bootstrap -- pre-resampled correlation/weight arrays. Split out from efa_fit() so multiple-
# imputation pooling can drive the same estimate -> rotate -> standard-error pipeline on
# pooled inputs (the MI2S route in EFA_POOLED()) without re-entering efa_fit()'s raw-data
# preparation and input guards. `fiml` (point estimate) and `fiml_boot` (per-replicate)
# carry the two-stage EM moments + raw data for the FIML likelihood-ratio fit indices and
# default to NULL; EFA_POOLED()'s MI2S route never supplies them (it is gated to poly/tetra/
# pearson), so its fit indices stay on the standard discrepancy path.
.efa_core <- function(R, N, weights = NULL, Gamma = NULL,
                      R_boot_array = NULL, W_boot_array = NULL,
                      fiml = NULL, fiml_boot = NULL,
                      np_boot = FALSE, b_boot = 1000, estimator, rotation, type,
                      rot_type = type,
                      n_factors, se = "none", ci = .95,
                      use = "pairwise.complete.obs", cor_method = "pearson",
                      max_iter = NA, init_comm = NA, criterion = NA,
                      criterion_type = NA, abs_eigen = NA, start_method = "psych",
                      normalize = TRUE, precision = 1e-5, order_type = NA,
                      varimax_type = NA, P_type = NA, k = NA, randomStarts = 100,
                      ...) {

  # Check if model is identified

  # calculate degrees of freedom
  m <- ncol(R)
  df <- .efa_df(m, n_factors)

  if(df < 0){

    cli::cli_warn(
      c("The model is underidentified; no fit indices were computed.",
        "i" = "Use fewer factors or more indicators."),
      class = "efa_underidentified"
    )

  } else if (df == 0){

    cli::cli_warn(
      c("The model is just identified ({.code df = 0}).",
        "i" = "Consider fewer factors or more indicators."),
      class = "efa_just_identified"
    )

  }

  # run factor analysis with the respective estimator

  if (estimator %in% c("ML", "ULS", "DWLS")) {

    if (type == "SPSS") {

      cli::cli_warn(
        c("Only {.val PAF} is validated against the SPSS implementation.",
          "i" = "{.val {estimator}} results may differ from those returned by SPSS."),
        class = "efa_spss_method_untested"
      )

    }

    if (is.na(N)) {

      cli::cli_warn(
        c("{.arg N} is {.val NA}; not all fit indices could be computed.",
          "i" = "Provide {.arg N} or raw data to compute all fit indices."),
        class = "efa_fit_na_n"
      )

    }

  }

  fit_out <- .estimate_model(R, method = estimator, n_factors = n_factors, N = N,
                             type = type, max_iter = max_iter,
                             init_comm = init_comm, criterion = criterion,
                             criterion_type = criterion_type,
                             abs_eigen = abs_eigen, start_method = start_method,
                             weights = weights, fiml = fiml)

  # Surface Heywood cases from the point-estimate solution (the detector runs in
  # .finalize_fit for every fit). This fires once per efa_fit() call; efa_average(),
  # which fits one EFA per grid cell, suppresses these per-model warnings and
  # reports a single summary instead.
  if (length(fit_out$heywood) > 0) {
    # Fall back to the indices if the correlation matrix carried no variable names.
    heywood_vars <- names(fit_out$heywood) %||% as.character(fit_out$heywood)
    cli::cli_warn(
      c(paste("{cli::qty(heywood_vars)}Heywood case{?s} detected for {.val {heywood_vars}}:",
              "the solution is improper (a communality at or above 1, or a uniqueness",
              "fixed at the estimation boundary)."),
        "i" = "Interpret the affected loadings and uniquenesses with caution."),
      class = "efa_heywood"
    )
  }

  # Surface optimiser non-convergence from the point-estimate solution for the
  # iterative estimators (ML, ULS, DWLS), whose fitters return a non-zero
  # convergence code when the optimiser stops before meeting its tolerance. PAF
  # raises its own non-convergence warning from inside .PAF(). This fires once per
  # efa_fit() call; the bootstrap replicates suppress their per-fit warnings and are
  # tallied separately in .boot_fun().
  if (estimator %in% c("ML", "ULS", "DWLS") && isTRUE(fit_out$convergence != 0)) {
    cli::cli_warn(
      c("The {.val {estimator}} optimiser did not converge (convergence code {fit_out$convergence}).",
        "i" = paste("It stopped before meeting the convergence tolerance (typically the maximum",
                    "number of iterations was reached); the results may not be interpretable."),
        "i" = paste("Try extracting fewer factors, or check the correlation matrix for",
                    "near-collinear variables.")),
      class = "efa_nonconvergence"
    )
  }

  # Only the bootstrap produces replicate fits; the analytic SE methods leave this NULL.
  boot_fits <- NULL
  if (isTRUE(np_boot)) {

    boot_fits <- .boot_fun(R_boot_array, b_boot, .estimate_model,
                           # .estimate_model arguments:
                           method = estimator, n_factors = n_factors, N = N,
                           type = type, max_iter = max_iter,
                           init_comm = init_comm, criterion = criterion,
                           criterion_type = criterion_type,
                           abs_eigen = abs_eigen, start_method = start_method,
                           # Each replicate fits only the quantities the bootstrap
                           # aggregation consumes (see .finalize_fit()).
                           lean = TRUE,
                           # DWLS carries the per-replicate weight matrices; NULL otherwise.
                           weights_array = W_boot_array,
                           # FIML carries the per-replicate EM moments + resampled data for the
                           # likelihood-ratio fit indices; NULL for the other cor_methods.
                           fiml_list = fiml_boot)

  }

  # rotate factor analysis results
  if (rotation == "none") {

    output <- fit_out
    boot_rot <- "none"

  } else {

    # Only promax and varimax are validated against the SPSS implementation.
    if (rot_type == "SPSS" && !rotation %in% c("promax", "varimax")) {

      cli::cli_warn(
        c("Only the {.val promax} and {.val varimax} rotations are validated against the SPSS implementation.",
          "i" = "{.val {rotation}} results may differ from those returned by SPSS."),
        class = "efa_spss_rotation_untested"
      )

    }

    rot_out <- .rotate_model(fit_out, rotation = rotation, type = rot_type,
                             normalize = normalize, precision = precision,
                             order_type = order_type, varimax_type = varimax_type,
                             P_type = P_type, k = k, randomStarts = randomStarts,
                             ...)

    boot_rot <- .rotation_family(rotation)

  }

  if (rotation != "none"){

    if(estimator %in% c("ULS", "DWLS")){

      settings <- rot_out$settings
      output <- c(fit_out, within(rot_out, rm(settings)),
                  settings = list(settings))

    } else {

      settings <- c(fit_out$settings, rot_out$settings)
      output <- c(within(fit_out, rm(settings)), within(rot_out, rm(settings)),
                  settings = list(settings))

    }

  }

  # Add settings used to output
  settings_EFA <- list(
    estimator = estimator,
    # back-compat alias: the output settings keep the former field name alongside the
    # current one, as they do for the frozen P_type/randomStarts rotate keys
    method = estimator,
    rotation = rotation,
    type = type,
    n_factors = n_factors,
    N = N,
    use = use,
    cor_method = cor_method,
    se = se,
    b_boot = b_boot,
    ci = ci
  )

  if(estimator %in% c("ULS", "DWLS") & rotation == "none"){

    output <- c(output, settings = list(settings_EFA))

  } else {

    settings <- c(settings_EFA, output$settings)

    output <- c(within(output, rm(settings)),
                settings = list(settings))

  }

  # Persist the full unrotated loading covariance (populated below for the analytic SE
  # methods) and the asymptotic covariance of the off-diagonal correlations (populated
  # for `se = "sandwich"`). Both stay present-but-NULL elsewhere so downstream consumers
  # (multiple-imputation pooling) can probe them by name without an `is.null(names(...))`
  # dance.
  output["vcov_unrot_loadings"] <- list(NULL)
  output["Gamma"] <- list(Gamma)

  if (se != "none") {
    if (rotation == "none") {
      L_rot <- NULL
      rot_info <- NULL
    } else {
      L_rot <- rot_out$rot_loadings
      # Analytic rotated SEs re-solve the rotation when finite-differencing its Jacobian, so they
      # need the converged transformation, the factor correlations, the resolved Kaiser-
      # normalization flag, and the criterion's tuning arguments (the same `.gpf_crit` defaults the
      # rotation itself used). NULL for the bootstrap, which carries its own replicate rotations.
      rot_info <- if (se %in% c("information", "sandwich")) {
        list(rotation = rotation,
             rotmat = rot_out$rotmat,
             rot_loadings = rot_out$rot_loadings,
             Phi = rot_out$Phi,
             normalize = rot_out$settings$normalize,
             crit_args = list(gam = .gpf_crit(list(...), "gam", 0),
                              delta = .gpf_crit(list(...), "delta", 0.01)))
      } else {
        NULL
      }
    }
    boot_out <- .compute_se_ci(fit_out, L_rot, se_method = se,
                               boot_fits = boot_fits, boot_rot = boot_rot,
                               ci = ci, b = b_boot, N = N, rot_info = rot_info,
                               gamma = Gamma, method = estimator, fiml = fiml)
    # The sandwich also returns the robust scaled chi-square block; it is patched into the
    # fit indices below rather than carried in the SE schema, so strip it before merging.
    scaled_test <- boot_out$scaled_test
    boot_out$scaled_test <- NULL
    output$SE <- boot_out$SE
    output$CI <- boot_out$CI
    # Single-bracket list assignment preserves a present-but-NULL slot for the
    # analytic SE methods (information, sandwich) where there are no replicate
    # arrays; `output$replicates <- NULL` would remove the slot entirely.
    output["replicates"] <- list(boot_out$replicates)
    # The analytic SE paths (information, sandwich) carry through the full unrotated
    # loading covariance the marginal SEs were derived from; the bootstrap path leaves it
    # NULL. Same single-bracket pattern as `replicates` so the slot is present-but-NULL
    # on the bootstrap rather than absent.
    output["vcov_unrot_loadings"] <- list(boot_out$vcov_unrot_loadings)
    # Only the bootstrap yields a sampling SD for every residual; the analytic methods do
    # not, so standardise the residuals only when those SEs are available.
    if (!is.null(boot_out$SE$residuals)) {
      # The residual diagonal is fixed at 0 with SE 0, so 0/0 would yield NaN on the
      # diagonal; set it to 0 so the off-diagonal standardised residuals are usable.
      std_resid <- output$residuals / boot_out$SE$residuals
      diag(std_resid) <- 0
      output$standardized_residuals <- std_resid
    }
    # Fill the chi-square block of the fit indices with the robust scaled chi-square
    # (Satorra-Bentler / scaled-shifted): .gof() leaves it undefined for DWLS, and for
    # ML/ULS the unscaled discrepancy chi-square is not robust to the ordinal weighting.
    if (se == "sandwich" && !is.null(scaled_test)) {
      output$fit_indices <- .apply_scaled_test(output$fit_indices, scaled_test, N)
    }
  }

  class(output) <- c("efa", "EFA")

  return(output)

}



.boot_fun <- function(x, b, call_fun, ..., weights_array = NULL, fiml_list = NULL) {

  # The per-replicate fits are independent and estimation is RNG-free, so they are
  # run in parallel across replicates at the R/process level with future.apply. A
  # parallel processing plan can be selected with future::plan(); the default plan
  # runs sequentially. future.seed = TRUE assigns each replicate its own reproducible
  # L'Ecuyer-CMRG stream, so the bootstrap is reproducible and independent of the
  # number of workers (any RNG a fitter might draw is bound to the replicate index,
  # never to the worker).
  #
  # A replicate whose (possibly degenerate) resampled correlation matrix cannot be
  # fit returns NULL and is skipped later, rather than aborting the whole call.
  #
  # Per-replicate warnings are suppressed here: they repeat identical information b
  # times. The type/preset-override notice depends only on the (type, pinned
  # arguments) combination, which is the same for every replicate and already
  # surfaced once by the point-estimate fit; the iterative fitter's max-iteration
  # warning would otherwise fire once per non-converged replicate. Non-convergence
  # is instead tallied and reported once, after all replicates have been fitted.
  boot_list <- future.apply::future_lapply(seq_len(b), function(boot_i) {
    # DWLS passes the replicate's own weight matrix; the other estimators ignore the
    # NULL (.estimate_model()'s weights default).
    w_i <- if (is.null(weights_array)) NULL else weights_array[,, boot_i]
    # FIML passes the replicate's own EM moments + resampled data for its likelihood-ratio fit
    # indices; the other paths ignore the NULL (.estimate_model()'s fiml default). Each element
    # carries an N x p resample, so the list adds to what future.apply serialises per worker.
    f_i <- if (is.null(fiml_list)) NULL else fiml_list[[boot_i]]
    tryCatch(suppressWarnings(call_fun(x[,, boot_i], ..., weights = w_i, fiml = f_i)),
             error = function(e) NULL)
  }, future.seed = TRUE)

  n_nonconverged <- sum(vapply(boot_list,
                               function(fit_i) isTRUE(fit_i$convergence != 0),
                               logical(1)))

  if (n_nonconverged > 0L) {
    cli::cli_warn(
      c("{n_nonconverged} of {b} bootstrap replicate{?s} did not converge.",
        "i" = "Their bootstrap standard errors and confidence intervals may be unreliable."),
      class = "efa_boot_nonconvergence"
    )
  }

  boot_list

}
