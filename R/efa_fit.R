#' Exploratory factor analysis (EFA)
#'
#' This function does an EFA with `PAF`, `ML`, `ULS`/`MINRES`, or `DWLS` estimation,
#' by default followed by an oblimin rotation.
#' Estimation and rotation are controlled through the control objects built by
#' [estimate_control()] and [rotate_control()]; each accepts a `type`
#' ("EFAtools", "SPSS", "psych", or "none") that fills in its remaining settings.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#' correlations. If raw data is entered, the correlation matrix is found from the
#' data.
#' @param n_factors numeric. Number of factors to extract. Must be at least 1 and
#' smaller than the number of variables (the common factor model is not identified
#' otherwise). Use [efa_retain()] to decide on a value.
#' @param N numeric. The number of observations. Needs only be specified if a
#' correlation matrix is used. Without `N`, the chi-square-based fit indices are not
#' available, with a warning; the same applies to a positive `N` that is very small
#' relative to the number of variables. With raw data, `N` is found from the data:
#'  - the number of rows under the default pairwise deletion (`use`), with a warning
#'    if values are missing, because each correlation then rests on fewer cases, so the
#'    chi-square is too large and the standard errors are too small (Marsh, 1998);
#'  - the number of complete cases under listwise deletion, and whenever an asymptotic
#'    covariance is needed (`estimator = "DWLS"` or `se = "sandwich"`, except under
#'    `cor_method = "fiml"`);
#'  - the number of rows with at least one observed value under `cor_method = "fiml"`.
#' @param estimator character. The estimator used to fit the EFA: "PAF" (principal axis
#' factoring), "ML" (maximum likelihood), "ULS" (unweighted least squares; "MINRES" is an
#' accepted alias returning identical results), or "DWLS" (diagonally weighted least
#' squares, for ordinal data). See the *Estimators* section in Details for their
#' properties and data requirements. Lower-case versions (e.g., "paf") are also accepted.
#' @param rotation character. An oblique rotation ("oblimin", the default, "promax",
#' "quartimin", "simplimax", "bentlerQ", "geominQ", or "bifactorQ"), an orthogonal
#' rotation ("varimax", "equamax", "quartimax", "geominT", "bentlerT", or "bifactorT"),
#' or no rotation ("none"). A single factor cannot be rotated, so with `n_factors = 1`
#' the default is "none". See the *Rotations* section in Details for their properties
#' and known issues.
#' @param se character. Whether and how to compute standard errors (and matching
#'  confidence intervals): "none" (default), "information" (analytic standard errors
#'  from the expected Fisher information of the ML solution), "sandwich" (robust
#'  standard errors from raw data, which stay reliable under non-normality or a
#'  misspecified estimator weight), or "np-boot" (non-parametric bootstrap). The methods
#'  differ in their assumptions, their data requirements, and which estimator, rotation,
#'  and `cor_method` combinations they support; see the *Standard errors* section in
#'  Details.
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
#'  estimation control arguments: the `type` preset; the principal-axis-factoring iteration
#'  settings `init_comm`, `criterion`, `criterion_type`, `max_iter`, and `abs_eigen`; and
#'  the maximum-likelihood `start_method`. Defaults to `estimate_control()`, which uses
#'  the `"EFAtools"` type. See [estimate_control()].
#' @param rotate_control a control object from [rotate_control()] bundling the rotation
#'  control arguments: the `type` preset; Kaiser `normalize`; the convergence `precision`; the
#'  factor `order_type`; the varimax/promax settings `varimax_type` and `p_type`; the
#'  simplimax/promax `k`; `random_starts`; and the rotation-engine arguments `gam`
#'  (oblimin), `delta` (geomin), and `maxit`. Defaults to `rotate_control()`, which uses
#'  the `"EFAtools"` type. The estimation and rotation presets are independent. See
#'  [rotate_control()] for details.
#' @param b_boot numeric. The number of bootstrap samples to draw. Default is 1000.
#'  Must be at least 2, the smallest number from which a standard error is defined.
#'  Under `cor_method = "fiml"` each bootstrap sample re-runs the EM moment
#'  estimation, so a smaller value may be advisable.
#' @param ci numeric. The level of the confidence intervals: the percentile intervals from
#'  the bootstrap samples under `se = "np-boot"`, and the analytic Wald intervals under
#'  `se = "information"` and `se = "sandwich"`, the corrected two-stage intervals of
#'  `cor_method = "fiml"` included. Must be greater than 0 and smaller than 1. Default is
#'  .95 for 95% CIs.
#' @param seed numeric. An optional seed for the random-number generator.
#'
#' @details
#'
#' ## Estimators
#'
#' The estimator is chosen with `estimator`.
#'
#' - **PAF** (principal axis factoring) iteratively estimates the communalities and makes
#'   no distributional assumptions, which makes it robust and a good general-purpose
#'   default. A converged PAF solution is a stationary point of the ULS criterion
#'   (Jöreskog & van Thillo, 1971), so PAF and ULS give practically the same solution and
#'   the same fit indices (see *Fit indices*).
#' - **ML** (maximum likelihood) maximises the normal-theory likelihood. It is the only
#'   estimator with analytic expected-information standard errors
#'   (`se = "information"`), but it assumes multivariate normality and is the most prone
#'   to Heywood (improper) cases. Its starting values are set by `start_method`.
#' - **ULS** (unweighted least squares) minimises the sum of squared correlation
#'   residuals. "MINRES" (minimum residual) is the same estimator under a different name
#'   and returns identical results. It makes no normality assumption and is robust to mild
#'   non-normality.
#' - **DWLS** (diagonally weighted least squares) is the recommended estimator for ordinal
#'   data. It weights each off-diagonal correlation residual by the inverse asymptotic
#'   variance of the corresponding polychoric correlation (Muthén, du Toit, & Spisic, 1997),
#'   reproducing the loadings of a diagonally weighted least squares fit. It therefore
#'   requires raw ordinal data with `cor_method = "poly"` or `"tetra"`. Because the
#'   weighting follows the polychoric asymptotic covariance, the matrix and the weights are
#'   estimated on the listwise-complete cases. Its fit-index behaviour is described under
#'   *Fit indices*.
#'
#' ## Correlation methods
#'
#' When raw data are supplied, `cor_method` selects how the correlation matrix is computed
#' (it is ignored when a correlation matrix is entered directly).
#'
#' - **"pearson"** (default), **"spearman"**, and **"kendall"** are passed to
#'   [stats::cor()]. The factor model assumes a Pearson correlation, but a rank
#'   correlation is analysed on its own scale, not converted to a Pearson-equivalent
#'   value; Kendall's tau in particular gives the most attenuated loadings of the three.
#'   For ordinal items prefer `"poly"` / `"tetra"` below, which estimate the correlation
#'   of the underlying continuous variables.
#' - **"poly"** / **"tetra"** compute polychoric / tetrachoric correlations for ordinal /
#'   binary data, assuming an underlying bivariate-normal latent variable. They use a
#'   two-step estimator. The polychoric asymptotic covariance that underlies both the DWLS
#'   weights and the scaled (sandwich) statistic relies on large-sample theory that
#'   degrades for empty or near-empty response-category combinations; with very sparse
#'   cells the resulting weights and standard errors can be unreliable (a warning is issued
#'   when empty cells are present), so interpret them with caution and consider collapsing
#'   rare categories. Each of the `p(p - 1)/2` variable pairs is a separate numerical
#'   optimisation, so a polychoric matrix takes much longer to compute than a Pearson one,
#'   and the difference grows quadratically in the number of variables; with
#'   `se = "np-boot"` the whole matrix is re-estimated for every bootstrap replicate.
#' - **"fiml"** estimates a two-stage full-information maximum-likelihood correlation. The
#'   saturated multivariate-normal mean and covariance are estimated from raw data with
#'   missing values by an EM algorithm assuming the data are missing at random (Yuan,
#'   Marshall, & Bentler, 2002; Little & Rubin, 2002), and the standardized covariance is
#'   then analysed. The model fit indices are corrected two-stage statistics wherever the
#'   correction can be formed (see *Fit indices*). `"fiml"` uses every case and handles the
#'   missingness itself, so `use` is ignored; it supplies a continuous (Pearson-type)
#'   correlation only and is therefore not compatible with `estimator = "DWLS"`. Standard
#'   errors are available analytically for `estimator = "ML"` or `"ULS"` and, for any
#'   estimator, by the non-parametric bootstrap (see *Standard errors*). For multiply
#'   imputed data, [efa_mi()] is the alternative route to handling missingness. Both routes
#'   assume the values are missing at random (MAR), and which one to prefer is largely
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
#' The default, oblimin, is therefore oblique: if `Phi` comes back near zero the oblique
#' solution is essentially the orthogonal one anyway, whereas imposing orthogonality on
#' genuinely correlated factors distorts the pattern.
#'
#' For an oblique solution the **pattern** matrix (`rot_loadings`) holds each variable's
#' unique contribution from each factor, with the other factors partialled out; it is what
#' is normally interpreted and reported. The **structure** matrix
#' (`Structure = rot_loadings %*% Phi`) holds the plain variable-factor correlations, which
#' are inflated by the factor intercorrelations. The two coincide only when `Phi` is the
#' identity, which is why an orthogonal rotation returns `rot_loadings` alone.
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
#'   plus group factors (bifactor simple structure). It needs at least two group factors
#'   (`n_factors >= 3`): with two factors the criterion is identically zero, so no rotation
#'   is performed and the unrotated loadings are returned with a warning.
#'
#' Oblique rotations:
#' - **oblimin** (the default) is quartimin with a tunable argument, `gam` (default 0, i.e.
#'   quartimin itself); turning it up trades some row simplicity for more strongly
#'   correlated factors, and can drive the solution toward factor collapse, so inspect
#'   `Phi` before interpreting a fit with `gam > 0`.
#' - **promax** is a fast two-step rotation: a varimax solution is raised to a power
#'   (controlled by `k` and `p_type`) to form a target that is then fitted obliquely.
#' - **quartimin** simplifies the variables (rows), like quartimax, so each loads mainly
#'   on one factor, but factors are allowed to correlate.
#' - **simplimax** drives the `k` smallest loadings toward zero. Its criterion is only
#'   piecewise smooth and strongly multimodal, which makes it by far the most
#'   start-dependent rotation offered here: different random seeds can reach noticeably
#'   different solutions. Because every start is fully optimised (there is no screening
#'   stage), raising `random_starts` to several hundred costs proportionally more time but
#'   buys a better optimum.
#' - **bentlerQ** is the oblique Bentler invariant pattern simplicity criterion.
#' - **geominQ** is the oblique geomin criterion; it handles complex (cross-loading)
#'   structure well but is multimodal, so it benefits from more `random_starts` (and uses a
#'   more thorough multi-start search internally).
#' - **bifactorQ** is the oblique (correlated) Jennrich-Bentler bifactor criterion, with the
#'   same two-group-factor requirement as **bifactorT**.
#'
#' The criterion-based rotations (all except varimax and promax) are fitted by gradient
#' projection with `random_starts` random starts to guard against local minima; the
#' complexity criteria (simplimax and geominQ in particular) are the most multimodal. The
#' starts are drawn from the random-number generator, so different starts can reach
#' genuinely different optima and such a fit is reproducible only when the generator is
#' controlled: pass `seed`, or call [base::set.seed()] beforehand.
#'
#' ## Standard errors
#'
#' `se` selects whether and how standard errors (and matching confidence intervals) are
#' computed. The analytic methods (`"information"` and `"sandwich"`) cover the unrotated
#' loadings, the uniquenesses and the communalities and, when a rotation is applied, the
#' rotated loadings and -- for oblique rotations -- the factor correlations and the
#' structure coefficients. The bootstrap (`"np-boot"`) covers the unrotated loadings, the
#' residuals, where available RMSEA, CFI, and TLI, and, when a rotation is applied, the
#' rotated loadings and -- for oblique rotations -- the factor correlations and the
#' structure coefficients; it reports no uniqueness or communality standard errors.
#'
#' - **"none"** (default) computes no standard errors.
#' - **"information"** returns analytic standard errors from the expected (Fisher)
#'   information matrix of the maximum-likelihood solution, and therefore requires
#'   `estimator = "ML"` and `cor_method = "pearson"` (or `"fiml"`, see below). The rotated
#'   standard errors are obtained by propagating the unrotated-loading covariance through
#'   the rotation (Jennrich, 1973). Unlike the bootstrap it also works from a correlation
#'   matrix as long as `N` is supplied. These standard errors assume multivariate normality
#'   and a correctly specified model; under heavy-tailed data or model misfit they can
#'   understate the sampling variability, where `"sandwich"` or `"np-boot"` are more robust.
#'
#'   The rotated loadings, `Phi`, the structure coefficients, the uniquenesses, and the
#'   communalities do not depend on how the unrotated solution happens to be oriented, and
#'   so are comparable across programs. The **unrotated** loading standard errors are not:
#'   a program using a different orientation will report different unrotated loading
#'   standard errors for the same fit. That orientation can also fail to be well defined --
#'   for example when two factors are of near-equal strength -- in which case the
#'   unrotated loadings have no well-defined sampling distribution and their standard
#'   errors diverge. `efa_fit()` detects this and returns `NA` for the unrotated loading
#'   standard errors with a warning; the other quantities are unaffected. Use those, or
#'   `se = "np-boot"`, when the unrotated loadings themselves are the quantity of
#'   interest. This detection applies to every analytic route. A Heywood case (a
#'   uniqueness at its lower boundary) is more severe: there the whole `SE`/`CI` block
#'   comes back `NA` with a warning, because the standard-error approximation fails for
#'   every parameter. The scaled chi-square (from `"sandwich"`, or from
#'   `cor_method = "fiml"`) does not rely on that approximation and is still reported.
#' - **"sandwich"** returns robust ("sandwich") standard errors estimated from raw data,
#'   combining the estimator weight with a distribution-free covariance of the
#'   correlations, so it stays valid under non-normality and weight misspecification
#'   (Browne, 1984; Satorra & Bentler, 1994). It is available either for ordinal data with
#'   `cor_method = "poly"` or `"tetra"` and `estimator` one of `"ML"`, `"ULS"`, or
#'   `"DWLS"`, or for continuous data with `cor_method = "pearson"` and
#'   `estimator = "ML"` or `"ULS"`. It reports the same coverage as `"information"` and
#'   additionally fills the model fit's chi-square block with a scaled chi-square (see
#'   *Fit indices*); the statistic reported as `chi` is always the scaled-and-shifted one
#'   (flagged by `chi_scaled_type`), and a mean-adjusted alternative is returned alongside
#'   it as `chi_mean_adjusted`. Because the covariance must describe the same cases as the
#'   correlation matrix, the sandwich (like `estimator = "DWLS"`) is computed on the
#'   listwise-complete cases; on data with missing values the reported `N`, the
#'   correlation matrix, and the point estimate therefore reflect the complete cases
#'   regardless of `use`.
#' - **"np-boot"** draws a non-parametric (case-resampling) bootstrap and needs raw data.
#'   It is the most general method -- available for any `estimator`, `rotation`, and
#'   `cor_method` -- and the most robust to non-normality and misfit, at the cost of speed;
#'   its intervals are bootstrap percentile intervals. The replicates run with the `future`
#'   framework; by default they run sequentially, but registering a plan with
#'   [future::plan()] runs them in parallel instead, as the examples show. With a fixed
#'   `seed` the bootstrap is reproducible and yields the same result regardless of the
#'   number of workers. Under `cor_method = "fiml"` each resample also re-runs the EM
#'   moment estimation and is therefore slow, so a smaller `b_boot` may be advisable.
#'
#'   The intervals of RMSEA, CFI, and TLI come from a second set of fits. A resample of
#'   the raw cases adds sampling misfit to the misfit of the sample, which moves these
#'   intervals toward worse fit. Thus, `efa_fit()` first transforms the cases so that their
#'   misfit is equal to the estimated population misfit. Then it resamples the same rows
#'   (Yuan, Hayashi, & Yanagihara, 2007). If the RMSEA is 0, the transformed cases fit the
#'   model exactly. In simulations with normal data, the intervals of RMSEA and CFI had good
#'   coverage. Use them with caution if the CFI is below approximately .76 or the RMSEA is
#'   above approximately .10 (Zhang & Savalei, 2016). TLI is CFI rescaled by a constant, so
#'   its interval behaves in the same way. The transformation needs complete cases and
#'   Pearson correlations. Thus, these intervals need `cor_method = "pearson"`, no missing
#'   values in the analysed cases (complete data, or listwise deletion), and a chi-square
#'   with df > 0. Otherwise, CFI and TLI are printed without an interval. The other fit
#'   indices get no bootstrap interval: SRMR had poor coverage with either bootstrap, and
#'   Zhang and Savalei (2016) did not study the others. The second set of fits is not
#'   rotated, so it usually takes less time than the first.
#'
#' The analytic methods (`"information"` and `"sandwich"`) are not available with the
#' `"promax"` or `"simplimax"` rotations, which have no supported analytic route for the
#' rotated standard errors; use `"np-boot"` there. They are not available for PAF either;
#' use ULS, which gives the same solution at convergence, or `"np-boot"`. Under
#' `cor_method = "fiml"`, `"information"` and `"sandwich"` both return, for
#' `estimator = "ML"` or `"ULS"`, corrected two-stage sandwich standard errors (Yuan &
#' Bentler, 2000; Savalei & Bentler, 2009).
#'
#' ## Fit indices
#'
#' For ML, ULS, and PAF, [efa_fit()] returns the model chi-square (with its p-value and
#' degrees of freedom), the Comparative Fit Index (CFI; Bentler, 1990), the Tucker-Lewis
#' Index (TLI, also called the non-normed fit index; Tucker & Lewis, 1973), the Root Mean
#' Square Error of Approximation (RMSEA) with its 90% confidence interval (Browne & Cudeck,
#' 1992), the Akaike and Bayesian Information Criteria (AIC, BIC), the Expected
#' Cross-Validation Index (ECVI; Browne & Cudeck, 1989), the Root Mean Squared Residual
#' (RMSR), the Standardized Root Mean Squared Residual (SRMR; Bentler, 1995), and the
#' common-part-accounted-for (CAF) index (Lorenzo-Seva, Timmerman, & Kiers, 2011). They
#' come with the independence-baseline statistics `chi_null`, `df_null`, and `p_null`. On
#' the unscaled paths `chi_null` is Bartlett's test of sphericity; a scaled chi-square, and
#' the two-stage statistic of a `cor_method = "fiml"` fit, each carry their own baseline
#' instead. The baseline degrees of freedom `df_null` are \eqn{p(p - 1)/2} for \eqn{p}
#' variables.
#'
#' The print method shows the chi-square, CFI, TLI, RMSEA, and SRMR; `summary()` adds AIC,
#' BIC, and ECVI, which are only useful to compare models. With `se = "np-boot"`, CFI and TLI
#' also show their bootstrap interval (see *Standard errors*). RMSEA always shows its
#' analytic 90% interval; its bootstrap interval is in `CI$fit_indices`. RMSR is the root
#' mean square of the off-diagonal residuals; SRMR rescales it by a fixed factor that
#' depends only on the number of variables, and both are `NA` if a residual is unavailable. (\pkg{psych}'s
#' `rms` uses a different divisor and equals RMSR\eqn{/\sqrt{2}}.) CAF is returned for
#' [efa_hull()] but not printed: although Lorenzo-Seva et al. describe it as ranging from 0
#' to 1, a correctly specified model gives a CAF near 0.5, because the KMO of residuals
#' that hold only sampling noise is near 0.5. Read it as a relative measure. It is `NA`
#' where that KMO is undefined.
#'
#' The model chi-square is the Bartlett-corrected discrepancy (matching
#' [stats::factanal()] for ML). For ULS and PAF it is the same maximum-likelihood
#' discrepancy, evaluated at the fitted solution, as in [psych::fa()]; its chi-square
#' distribution is then an approximation. It is not the least-squares criterion
#' \pkg{lavaan} reports as its standard ULS test statistic, so the two are not comparable.
#' AIC and BIC are built on this chi-square and can therefore be negative; ECVI (built on
#' the same chi-square plus a non-negative penalty) cannot. Because of the Bartlett
#' correction, ECVI differs slightly from the uncorrected Browne-Cudeck form reported by
#' \pkg{lavaan} and Mplus (their AIC/BIC use an unrelated log-likelihood-based formula, so
#' they are not comparable). On the unscaled paths, CFI, TLI, and RMSEA are computed on the
#' uncorrected \eqn{(N - 1)} scale rather than on the Bartlett-corrected chi-square, so you
#' cannot recompute one from the other by hand there; on the scaled (sandwich) and FIML
#' paths, the reported chi and chi_null are exactly the pair the indices use.
#'
#' Two cases differ from the full set:
#' - **DWLS** by default returns only RMSR, SRMR, CAF, and df, because the ordinary
#'   maximum-likelihood discrepancy is not its fit function. When `se = "sandwich"`, a
#'   scaled chi-square and the CFI, TLI, and RMSEA derived from it are reported; that
#'   statistic is a two-stage correction applied to the polychoric correlation residuals
#'   (Browne, 1984), not identical to the full WLSMV test of \pkg{lavaan} or Mplus, which
#'   also corrects for the response thresholds.
#' - **`cor_method = "fiml"`** reports two-stage-corrected statistics (Yuan, Marshall, &
#'   Bentler, 2002). The correction itself can be degenerate -- typically with a small
#'   sample, a high proportion of missing values, or near-collinear variables -- in which
#'   case an uncorrected likelihood-ratio statistic is reported in its place with a
#'   warning, and the print methods label that line `uncorrected` rather than `scaled`;
#'   read its p-value and the CFI, TLI, and RMSEA derived from it as indicative only.
#'
#' Whenever the chi-square is a scaled one (`se = "sandwich"`, or a `cor_method = "fiml"`
#' fit), AIC, BIC, and ECVI are `NA`; see the `fit_indices` entry in Value for the
#' additional components then returned. Beyond that, the chi-square and everything derived
#' from it are `NA` whenever the statistic itself is undefined: when `N` is not supplied,
#' when the model is underidentified (a negative df), and when a positive `N` is too small
#' relative to the number of variables and factors for the small-sample correction to
#' remain valid. Each case raises its own warning. The residual summaries (RMSR, SRMR,
#' CAF) and the degrees of freedom are still returned there, but residual size does not
#' establish that a model is identified: below zero degrees of freedom a near-zero
#' residual is an artefact of over-parameterisation, not close fit. A just-identified model
#' (df = 0) has no test of fit, so the print method shows only SRMR and df.
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
#' \item{vars_accounted}{Matrix of explained variances and sums of squared loadings. Based
#' on the unrotated loadings. Its rows are `"SS loadings"` and `"Prop Tot Var"`, followed by
#' `"Cum Prop Tot Var"`, `"Prop Comm Var"`, and `"Cum Prop Comm Var"`; those last three are
#' omitted for a single-factor solution, where they would only repeat the two above them. So
#' code that indexes a row by name (for example `["Prop Comm Var", ]`) must allow for the
#' two-row form at `n_factors = 1`.}
#' \item{fit_indices}{A named list of fit indices computed from the unrotated
#' loadings (see *Fit indices* in Details): the model Chi Square (`chi`, with `p_chi` and
#' `df`), `CFI`, `TLI`, `RMSEA` with its 90% confidence interval (`RMSEA_LB`, `RMSEA_UB`),
#' `AIC`, `BIC`, `ECVI`, `RMSR`, `SRMR`, `CAF`, and the independence-model baseline
#' `chi_null`, `df_null`, and `p_null` that CFI and TLI are computed from. Whenever the Chi
#' Square is scaled, a few extra scaled-statistic fields are appended for advanced
#' diagnostics. The list also carries `Fm`, the estimator's own objective value at the
#' solution -- on its own scale, not the discrepancy the Chi Square is built from.}
#' \item{model_implied_R}{The model implied correlation
#' matrix.}
#' \item{residuals}{Residual correlations, i.e., orig_R - model_implied_R}
#' \item{standardized_residuals}{Residual correlations standardized by their
#'  bootstrap standard errors. Only returned, if `se = "np-boot"`.}
#' \item{rot_loadings}{Loading matrix containing the final rotated loadings. For an oblique
#' rotation this is the pattern matrix -- each variable's unique contribution from each
#' factor, with the other factors partialled out -- and is the matrix normally interpreted
#' (see *Rotations*).}
#' \item{Phi}{The factor intercorrelations (only for oblique rotations).}
#' \item{Structure}{The structure matrix `rot_loadings %*% Phi`, holding the plain
#' variable-factor correlations, which are inflated by the factor intercorrelations (only
#' for oblique rotations).}
#' \item{rotmat}{The rotation matrix. The rotated loadings are recovered from the
#' unrotated loadings as `unrot_loadings %*% rotmat` for orthogonal rotations and
#' for promax, and as `unrot_loadings %*% t(solve(rotmat))` for the other oblique
#' rotations.}
#' \item{vars_accounted_rot}{Matrix of explained variances and sums of squared
#' loadings. Based on rotated loadings and, for oblique rotations, the factor
#' intercorrelations. Same rows as `vars_accounted`; it is returned only when a rotation
#' was actually applied, which never happens for a single-factor solution, so it always
#' has the five-row form.}
#' \item{settings}{A list of the settings used, including `seed` (`NULL` when none was
#' supplied), `input_type` (`"raw"` or `"correlation"`, what `x` was), and
#' `cor_method_used` (the correlation method that actually ran, `NA_character_` for a
#' correlation-matrix input, which consumes none). `cor_method` keeps the requested value
#' whether or not it was used. For the criterion rotations fitted by gradient projection
#' it additionally carries `rotation_diagnostics`, a list summarising the multi-start run:
#'  - `converged`: whether the start whose solution is returned reached the convergence
#'    tolerance (reported separately from the counts below because the returned solution
#'    is the one with the lowest criterion value, which need not be a converged one).
#'  - `n_starts_total`: the `random_starts` random starts plus the rational start.
#'  - `n_optimized`: how many of those starts were actually optimized -- fewer than
#'    `n_starts_total` whenever the solver screens the random starts and optimizes only
#'    the most promising ones.
#'  - `n_converged`: how many optimized starts reached the convergence tolerance.
#'  - `n_distinct_minima`: how many distinct local optima those converged starts found;
#'    more than one means the criterion is multimodal on these data.
#'  - `criterion_spread`: the range of the criterion values they attained.
#'  - `criterion_best`: the criterion value of the returned solution.
#'
#'  When `normalize = TRUE`, `criterion_best` and `criterion_spread` are evaluated on
#'  the Kaiser-normalized loadings the criterion is optimized on, not on the returned
#'  `rot_loadings`, so they are not directly comparable to a criterion recomputed from
#'  the returned loadings.}
#' \item{fiml}{Diagnostics of the FIML correlation's EM estimation, present only for
#' `cor_method = "fiml"`: whether it `converged` before `fiml_max_iter`, how many
#' iterations it took (`iter`), the number of distinct missingness patterns
#' (`n_patterns`), and the number of cases used (`n`, the reported `N`). A `converged`
#' of `FALSE` means the analysed correlation is the EM's last iterate, not the true
#' FIML estimate; this is flagged in the printed output. Set
#' `fiml_max_iter` and `fiml_tol` through [estimate_control()].}
#' \item{SE}{A named list of standard error matrices, returned only when `se` is not
#' `"none"` (see *Standard errors* in Details for what each method assumes and when a
#' component comes back `NA`). For `se = "np-boot"`: bootstrap SDs of the loadings
#' (rotated too, if a rotation was applied), the residuals, and, where available, RMSEA,
#' CFI, and TLI (`fit_indices`, from the transformed cases), plus
#' -- for oblique rotations -- `Phi` and the structure coefficients; `valid_replicates`
#' (and, when rotated, `valid_target_rotations`) records how many bootstrap replicates
#' each of those is based on (for `fit_indices`, count the rows of
#' `replicates$fit_indices` without `NA`). For `se = "information"` and `se = "sandwich"`: Wald-type
#' SEs for the loadings, the uniquenesses, and the communalities -- identical to each
#' other, since communality is just `1 - uniqueness` -- plus `Phi` and the structure
#' coefficients for oblique rotations. `"sandwich"` stays valid under non-normality;
#' `"information"` assumes the model is correctly specified.}
#' \item{CI}{A named list of confidence intervals of width `ci`, matching the
#' components of `SE`: percentile intervals for `se = "np-boot"`, Wald intervals for
#' `se = "information"` and `se = "sandwich"`. Only returned if `se` is not `"none"`.
#' See *Standard errors* in Details for when the bootstrap gives `fit_indices`.}
#' \item{replicates}{A named list of raw bootstrap replicate arrays -- the aligned
#' loadings, `Phi`, structure coefficients, residuals, and RMSEA, CFI, and TLI behind `SE`
#' and `CI` -- with one replicate per array's last dimension (first dimension for
#' `fit_indices`). Failed replicates are left `NA`. Populated only for
#' `se = "np-boot"`; `NULL` otherwise.}
#' \item{vcov_unrot_loadings}{The full unrotated-loading covariance matrix behind
#' `SE$unrot_loadings`: a `p * n_factors` by `p * n_factors` matrix in column-major
#' `vec(Lambda)` order, with rows and columns labelled `"<variable>_<factor>"`.
#' Populated for `se = "information"` and `se = "sandwich"` -- always the unrotated
#' block, even when a rotation is applied -- and `NA`-filled if the analytic covariance
#' is unreliable (a Heywood case or a singular information matrix); `NULL` for
#' `se = "np-boot"` and `se = "none"`. It can be populated even when
#' `SE$unrot_loadings` is `NA`: a weakly determined rotational orientation invalidates
#' only the marginal SEs, not the underlying covariance.}
#'
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
#' Latent variables analysis: Applications for developmental research (pp. 399–419). Sage.
#' @source Asparouhov, T., & Muthén, B. (2010). Simple second order chi-square
#' correction. Mplus Technical Appendix.
#' @source Muthén, B., du Toit, S. H. C., & Spisic, D. (1997). Robust inference using
#' weighted least squares and quadratic estimating equations in latent variable modeling
#' with categorical and continuous outcomes. Unpublished manuscript.
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
#' @source Bartlett, M. S. (1951). The effect of standardization on a Chi-square
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
#' @source Jöreskog, K. G., & van Thillo, M. (1971). New rapid algorithms for factor
#' analysis by unweighted least squares, generalized least squares and maximum likelihood
#' (Research Memorandum RM-71-5). Educational Testing Service.
#' @source Marsh, H. W. (1998). Pairwise deletion for missing data in structural equation
#' models: Nonpositive definite matrices, parameter estimates, goodness of fit, and adjusted
#' sample sizes. Structural Equation Modeling, 5, 22–36. doi: 10.1080/10705519809540087
#' @source Yuan, K.-H., Hayashi, K., & Yanagihara, H. (2007). A class of population
#' covariance matrices in the bootstrap approach to covariance structure analysis.
#' Multivariate Behavioral Research, 42, 261–281. doi: 10.1080/00273170701360662
#' @source Zhang, X., & Savalei, V. (2016). Bootstrapping confidence intervals for fit
#' indexes in structural equation modeling. Structural Equation Modeling, 23, 392–408.
#' doi: 10.1080/10705511.2015.1118692
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
#' # Principal axis factoring with oblimin rotation (the defaults)
#' mod_paf <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500)
#' mod_paf
#' summary(mod_paf)
#'
#' # ML estimation with promax rotation
#' mod_ml <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                   estimator = "ML", rotation = "promax")
#' mod_ml
#'
#' # Tuning knobs are supplied through the control objects. Here the SPSS preset is
#' # used for the estimation and rotation, with the maximum PAF iterations raised.
#' mod_spss <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                     rotation = "promax",
#'                     estimate_control = estimate_control(type = "SPSS", max_iter = 500),
#'                     rotate_control = rotate_control(type = "SPSS"))
#' mod_spss
#'
#' # Analytic (expected-information) standard errors of an ML fit
#' ML_info <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                    estimator = "ML", se = "information")
#' ML_info
#' summary(ML_info)
#'
#' \donttest{
#' # Robust (sandwich) standard errors and a scaled chi-square for ordinal raw data.
#' # These need a polychoric/tetrachoric correlation method and estimator ML, ULS, or DWLS.
#' DWLS_rob <- efa_fit(DOSPERT_raw, n_factors = 6, cor_method = "poly",
#'                     estimator = "DWLS", se = "sandwich")
#' DWLS_rob
#' summary(DWLS_rob)
#'
#' # The same robust SEs and scaled chi-square for continuous data: a Pearson
#' # correlation with estimator ML or ULS (the fourth-moment ADF covariance).
#' ML_rob <- efa_fit(GRiPS_raw, n_factors = 1, cor_method = "pearson",
#'                   estimator = "ML", se = "sandwich")
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
#' # in parallel across replicates. future::plan() returns the plan it replaces, so
#' # on.exit() puts the session back as it was -- also if the fit fails.
#' efa_boot <- local({
#'   old_plan <- future::plan(future::multisession, workers = 2)
#'   on.exit(future::plan(old_plan), add = TRUE)
#'   efa_fit(GRiPS_raw, n_factors = 1, estimator = "PAF",
#'           se = "np-boot", b_boot = 1000, seed = 42)
#' })
#' }
#'
efa_fit <- function(x, n_factors, N = NA,
                    estimator = c("PAF", "ML", "ULS", "MINRES", "DWLS"),
                    rotation = c("oblimin", "none", "varimax", "equamax", "quartimax",
                                 "geominT", "bentlerT", "bifactorT", "promax", "quartimin",
                                 "simplimax", "bentlerQ", "geominQ", "bifactorQ"),
                    se = c("none", "information", "sandwich", "np-boot"),
                    cor_method = c("pearson", "spearman", "kendall", "poly", "tetra",
                                   "fiml"),
                    use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                            "everything", "na.or.complete"),
                    estimate_control = NULL,
                    rotate_control = NULL,
                    b_boot = 1000, ci = .95, seed = NULL) {

  # A literal `estimate_control = estimate_control()` default would recurse into its own
  # promise, so the controls default to NULL and are built here.
  if (is.null(estimate_control)) estimate_control <- estimate_control()
  if (is.null(rotate_control)) rotate_control <- rotate_control()

  .assert_cor_input(x)

  if (missing(n_factors)) {
    cli::cli_abort(
      c("{.arg n_factors} is required.",
        "i" = "Use {.fn efa_retain} to decide how many factors to extract."),
      class = "efa_missing_n_factors"
    )
  }
  # A retention object holds one suggested count per criterion or variant, so name the
  # component that holds them and ask the reader to choose.
  if (inherits(n_factors, c("efa_retain", "efa_retention"))) {
    n_suggested <- length(n_factors$n_factors)
    cli::cli_abort(
      c("{.arg n_factors} must be a single count, not an {.cls {class(n_factors)[[1]]}} object.",
        "i" = "Its {.code $n_factors} holds {n_suggested} suggested count{?s}; supply
               {?it/one of them}."),
      class = "efa_n_factors_object"
    )
  }

  estimator <- .match_arg_ci(estimator)
  # "MINRES" is an alias of "ULS"; the settings report the canonical name.
  if (estimator == "MINRES") estimator <- "ULS"
  # A single factor cannot be rotated, so the default rotation becomes "none" there without
  # the warning that an explicit rotation request gives.
  if (missing(rotation) && isTRUE(n_factors == 1)) rotation <- "none"
  rotation <- .match_arg_ci(rotation)
  se <- .match_arg_ci(se)
  use <- .match_arg_ci(use)
  cor_method <- .match_arg_ci(cor_method)

  .assert_estimate_control(estimate_control)
  if (!inherits(rotate_control, "efa_rotate_control")) {
    cli::cli_abort(
      "{.arg rotate_control} must be a control object from {.fn rotate_control}.",
      class = "efa_control_input"
    )
  }
  ec <- estimate_control
  rc <- rotate_control

  if (is.na(ec$start_method) && estimator == "ML") {
    cli::cli_abort(
      c("{.arg start_method} must be set when {.code estimator = \"ML\"}.",
        "i" = "Set {.arg start_method} to {.val psych} or {.val factanal}."),
      class = "efa_ml_start_missing"
    )
  }

  .assert_args({
    checkmate::assert_count(n_factors)
    # N = 0 is an impossible sample size, not an unknown one; NA means "unknown".
    checkmate::assert_count(N, positive = TRUE, na.ok = TRUE)
    checkmate::assert_integerish(b_boot, len = 1, any.missing = FALSE)
    checkmate::assert_number(ci)
    checkmate::assert_int(seed, null.ok = TRUE)
  })
  # A bootstrap standard error is the spread across replicates, so it needs two of them.
  if (b_boot < 2) {
    cli::cli_abort(
      c("{.arg b_boot} must be at least 2.",
        "x" = "You supplied {.arg b_boot} = {b_boot}.",
        "i" = "A bootstrap standard error is the spread across replicates and is undefined below
               two of them. The default is {.val {1000}}."),
      class = "efa_b_boot_too_small"
    )
  }
  # At ci = 1 every Wald bound is infinite, and at ci = 0 every interval has zero width.
  if (ci <= 0 || ci >= 1) {
    cli::cli_abort(
      c("{.arg ci} must be greater than 0 and smaller than 1.",
        "x" = "You supplied {.arg ci} = {ci}.",
        "i" = "At 1 every analytic bound is infinite, and at 0 every interval has zero width.
               The default is {.val {0.95}}."),
      class = "efa_ci_out_of_bounds"
    )
  }

  # The common factor model needs at least one factor and fewer factors than variables. Check
  # both here, so that the bound checks of the compiled fitters are never the user-facing error.
  n_vars <- ncol(x)
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

  is_cormat <- .is_cormat(x)
  .check_fit_support(estimator, rotation, se, cor_method, is_cormat)
  # Analytic standard errors scale with N, which a correlation matrix does not carry.
  if (se == "information" && is_cormat && is.na(N)) {
    cli::cli_abort(
      c("{.code se = \"information\"} requires the sample size {.arg N}.",
        "i" = "Supply {.arg N} when {.arg x} is a correlation matrix."),
      class = "efa_se_no_n"
    )
  }
  if (cor_method == "fiml" && use != "pairwise.complete.obs") {
    cli::cli_warn(
      c("{.arg use} is ignored when {.code cor_method = \"fiml\"}.",
        "i" = "FIML uses every case and handles the missing values itself, so
               {.code use = {.val {use}}} has no effect."),
      class = "efa_fiml_use_ignored"
    )
  }

  # A fixed seed makes the whole fit reproducible (the random starts of the rotation and the
  # bootstrap resamples) and leaves the caller's random-number stream as it was.
  .set_local_seed(seed)

  prep <- .prepare_cor_input(x, N = N, use = use, cor_method = cor_method,
                             N_policy = "optional", inform_from_data = FALSE,
                             # The sandwich needs the full asymptotic covariance of the
                             # correlations; DWLS needs its diagonal for the weights.
                             acov = if (se == "sandwich" && cor_method != "fiml") "full"
                                    else if (estimator == "DWLS") "diag" else "none",
                             dwls = estimator == "DWLS",
                             fiml_max_iter = ec$fiml_max_iter, fiml_tol = ec$fiml_tol,
                             # Only PAF has a fallback (a pseudo-inverse start) for a singular
                             # matrix, so only the psych preset with PAF skips the check.
                             check_singular = !(ec$type == "psych" && estimator == "PAF"),
                             posdef_abort = ec$type == "SPSS")
  N <- prep$N

  # Pairwise deletion estimates each correlation from its own subset of cases, but N counts every
  # row (as in psych and lavaan). The chi-square is then too large and the standard errors too
  # small (Marsh, 1998).
  if (!is_cormat && use == "pairwise.complete.obs" && cor_method != "fiml" &&
      anyNA(x) && N == nrow(x)) {
    n_pair <- crossprod(!is.na(x))
    n_pair <- unique(range(n_pair[upper.tri(n_pair)]))
    n_text <- paste(n_pair, collapse = " to ")
    cli::cli_warn(
      c("{.arg N} is the number of rows ({N}), but with pairwise deletion each correlation
         uses only {n_text} cases.",
        "!" = "The chi-square test therefore rejects the model too often, and the standard
               errors are too small.",
        "i" = if (.is_poly_cor(cor_method)) {
          "For ordinal data with missing values, use multiple imputation with {.fn efa_mi}."
        } else {
          "For data with missing values, use {.code cor_method = \"fiml\"}."
        }),
      class = "efa_pairwise_n"
    )
  }

  # The EM moments and the raw data feed the two-stage FIML fit indices.
  fiml_pt <- if (cor_method == "fiml") c(list(data = x), prep$fiml)
  # An analytic-SE fit needs the saturated FIML covariance of the correlations twice (for the
  # scaled chi-square and for the sandwich), so build it once here. NULL on a degenerate
  # covariance, where each consumer recomputes it with its own guard.
  if (!is.null(fiml_pt) && se %in% c("information", "sandwich")) {
    fiml_pt$acov_cor <- tryCatch(
      .fiml_saturated_acov(fiml_pt$data, fiml_pt$mu, fiml_pt$sigma)$cor,
      error = function(e) NULL)
  }

  boot_rep <- NULL
  boot_idx <- NULL
  if (se == "np-boot") {
    # Resample the cases R was computed from (Efron & Tibshirani, 1993): the complete cases
    # under listwise deletion and for DWLS, the rows with at least one observed value under
    # FIML, and every row otherwise.
    rows <- if (cor_method == "fiml") {
      which(rowSums(!is.na(x)) > 0L)
    } else if (.is_listwise_use(use) || estimator == "DWLS") {
      which(stats::complete.cases(x))
    } else {
      seq_len(nrow(x))
    }
    # Draw all resamples here, in a fixed order, so that a seeded bootstrap gives the same
    # result for any number of parallel workers. The workers compute the correlations.
    boot_idx <- lapply(seq_len(b_boot), function(i) sample(rows, size = N, replace = TRUE))
    boot_rep <- .boot_rep_maker(x, cor_method = cor_method, use = use,
                                dwls = estimator == "DWLS",
                                fiml_max_iter = ec$fiml_max_iter, fiml_tol = ec$fiml_tol)
  }

  core_args <- list(
    R = prep$R, N = N, weights = prep$weights, Gamma = prep$Gamma, fiml = fiml_pt,
    boot_rep = boot_rep, boot_idx = boot_idx, np_boot = se == "np-boot", b_boot = b_boot,
    estimator = estimator, rotation = rotation, type = ec$type, rot_type = rc$type,
    n_factors = n_factors, se = se, ci = ci, use = use, cor_method = cor_method,
    max_iter = ec$max_iter, init_comm = ec$init_comm, criterion = ec$criterion,
    criterion_type = ec$criterion_type, abs_eigen = ec$abs_eigen,
    start_method = ec$start_method, normalize = rc$normalize, precision = rc$precision,
    order_type = rc$order_type, varimax_type = rc$varimax_type, P_type = rc$p_type,
    k = rc$k, randomStarts = rc$random_starts)
  out <- do.call(.efa_core, c(core_args, rc$extra_args))

  # Bootstrap intervals of RMSEA, CFI, and TLI from transformed cases (.yhy_fit_boot()). The
  # transformation is linear, so it needs Pearson correlations.
  if (se == "np-boot" && cor_method == "pearson" && !isTRUE(.efa_state$no_fit_boot)) {
    out <- .yhy_fit_boot(out, x, rows, boot_idx, use = use, ci = ci,
                         method = estimator, n_factors = n_factors, N = N, type = ec$type,
                         max_iter = ec$max_iter, init_comm = ec$init_comm,
                         criterion = ec$criterion, criterion_type = ec$criterion_type,
                         abs_eigen = ec$abs_eigen, start_method = ec$start_method)
  }

  # Provenance that .efa_core() cannot see: what `x` was and how the run was seeded.
  # `cor_method` keeps the requested value; `cor_method_used` says what actually ran.
  out$settings <- append(
    out$settings,
    list(input_type = if (is_cormat) "correlation" else "raw",
         cor_method_used = if (is_cormat) NA_character_ else cor_method),
    after = match("cor_method", names(out$settings)))
  out$settings <- append(out$settings, list(seed = seed),
                         after = match("ci", names(out$settings)))
  # The asymptotic covariance of the correlations is large (about 11 MB at p = 50) and only
  # efa_mi() reads it, so it is kept only on request (see .efa_state).
  if (!isTRUE(.efa_state$keep_gamma)) out$Gamma <- NULL
  out
}

# Package-level state of efa_fit(): `keep_gamma` is set by efa_mi() while it fits the
# imputations, because its sandwich route pools the per-imputation `Gamma` that efa_fit()
# otherwise drops. `no_fit_boot` is set by efa_group(), which reads only the replicate
# loadings, to skip the bootstrap of the fit indices.
.efa_state <- new.env(parent = emptyenv())

# Reject the estimator, rotation, standard-error, and correlation-method combinations that
# efa_fit() cannot fit, before any computation. Expected-information SEs come from the ML
# discrepancy and the normal-theory covariance of Pearson (or FIML) correlations. Sandwich SEs
# need a covariance of the correlations estimated from raw data: the polychoric one, or the
# fourth-moment one of Pearson correlations (Browne, 1984), which goes with ML or ULS. Neither
# analytic method has a rotation Jacobian for promax (a two-step target rotation) or simplimax
# (a non-smooth criterion).
.check_fit_support <- function(estimator, rotation, se, cor_method, is_cormat,
                               call = rlang::caller_env()) {
  stop_fit <- function(msg, class = "efa_se_unsupported") {
    cli::cli_abort(msg, class = class, call = call, .envir = parent.frame())
  }
  analytic <- se %in% c("information", "sandwich")
  boot_hint <- "Use {.code se = \"np-boot\"} instead."

  if (cor_method == "fiml") {
    if (is_cormat) {
      stop_fit(c("{.code cor_method = \"fiml\"} needs raw data, not a correlation matrix.",
                 "i" = "Supply raw data with missing values, or choose another {.arg cor_method}."),
               class = "efa_fiml_needs_raw")
    }
    if (estimator == "DWLS") {
      stop_fit(c("{.code estimator = \"DWLS\"} cannot be used with {.code cor_method = \"fiml\"}.",
                 "i" = "DWLS needs polychoric correlations. With FIML, use {.val ML}, {.val ULS},
                        or {.val PAF}."),
               class = "efa_fiml_unsupported_method")
    }
  }

  if (analytic && rotation %in% c("promax", "simplimax")) {
    stop_fit(c("{.code se = {.val {se}}} is not available with {.code rotation = {.val {rotation}}}.",
               "i" = boot_hint))
  }
  # PAF converges to the ULS solution, so ULS gives its analytic standard errors.
  paf_hint <- "Use {.code estimator = \"ULS\"}, which gives the same solution as PAF at
               convergence, or {.code se = \"np-boot\"}."
  if (analytic && cor_method == "fiml" && !estimator %in% c("ML", "ULS")) {
    stop_fit(c("Analytic standard errors with {.code cor_method = \"fiml\"} require
                {.code estimator = \"ML\"} or {.val ULS}.",
               "i" = paf_hint))
  }
  if (se == "information" && cor_method != "fiml") {
    if (estimator != "ML") {
      stop_fit(c("{.code se = \"information\"} is only available for {.code estimator = \"ML\"}.",
                 "x" = "You requested {.code estimator = {.val {estimator}}}.",
                 "i" = boot_hint))
    }
    if (cor_method != "pearson") {
      stop_fit(c("{.code se = \"information\"} requires {.code cor_method = \"pearson\"} or {.val fiml}.",
                 "x" = "You requested {.code cor_method = {.val {cor_method}}}.",
                 "i" = "The expected information assumes the normal-theory sampling behaviour of
                        Pearson correlations.",
                 "i" = if (is_cormat) {
                   "Supply the raw data and use {.code se = \"sandwich\"}."
                 } else {
                   "Use {.code se = \"sandwich\"} for {.val {cor_method}}."
                 }))
    }
  }
  if (se == "sandwich" && cor_method != "fiml") {
    if (estimator == "PAF") {
      stop_fit(c("{.code se = \"sandwich\"} is not available for {.code estimator = \"PAF\"}.",
                 "i" = paf_hint))
    }
    if (!.is_poly_cor(cor_method) && !(cor_method == "pearson" && estimator %in% c("ML", "ULS"))) {
      stop_fit(c("{.code se = \"sandwich\"} is not available for this correlation and estimator.",
                 "x" = "You requested {.code cor_method = {.val {cor_method}}} with
                        {.code estimator = {.val {estimator}}}.",
                 "i" = "Use {.code cor_method = \"poly\"} or {.val tetra} (any estimator but PAF),
                        or {.code cor_method = \"pearson\"} with {.val ML} or {.val ULS}.
                        Otherwise, use {.code se = \"np-boot\"}."))
    }
  }

  if (estimator == "DWLS" && !(.is_poly_cor(cor_method) && !is_cormat)) {
    stop_fit(c("{.code estimator = \"DWLS\"} requires a polychoric asymptotic covariance.",
               "x" = if (is_cormat) {
                 "You supplied a correlation matrix, so no asymptotic covariance can be estimated."
               } else {
                 "{.code cor_method = {.val {cor_method}}} does not produce one."
               },
               "i" = "Supply raw ordinal data with {.code cor_method = \"poly\"} or {.code \"tetra\"}."),
             class = "efa_dwls_no_acov")
  }
  if (se == "sandwich" && cor_method != "fiml" && is_cormat) {
    stop_fit(c("{.code se = \"sandwich\"} needs raw data, not a correlation matrix.",
               "i" = "The sandwich estimates the covariance of the correlations from the cases.
                      Supply raw data, or use {.code se = \"information\"} with {.arg N}."))
  }
  if (se == "np-boot" && is_cormat) {
    stop_fit(c("{.code se = \"np-boot\"} needs raw data, not a correlation matrix.",
               "i" = "A correlation matrix has no cases to resample. Supply raw data."),
             class = "efa_boot_cormat")
  }

  invisible(NULL)
}

# Explanation for a non-zero optimiser convergence code, branched on what the code means:
# 1 is the fixed iteration limit, 51 and 52 are L-BFGS-B's line-search warning and error
# (which are not iteration limits and are not fixed by running longer).
.nonconvergence_cause <- function(convergence) {
  if (isTRUE(convergence == 1)) {
    paste("It reached the optimiser's fixed iteration limit before meeting the",
          "convergence tolerance; that limit is not user-adjustable.")
  } else if (isTRUE(convergence %in% c(51, 52))) {
    paste("The bounded optimiser's line search failed, which typically means the",
          "objective is flat or ill-conditioned near the point it stopped at.")
  } else {
    "It stopped before meeting the convergence tolerance."
  }
}

# Fit the common-factor model from already-prepared inputs: a correlation matrix R, the
# sample size N, optional DWLS weights, the optional sandwich meat Gamma, the optional FIML
# moments, and -- for the bootstrap -- `boot_rep`, the per-replicate input function from
# .boot_rep_maker(), with `boot_idx`, the list of resample indices it takes. Split out from
# efa_fit() so that efa_mi() can run the same estimate -> rotate -> standard-error pipeline on
# pooled inputs without the raw-data preparation.
.efa_core <- function(R, N, weights = NULL, Gamma = NULL, fiml = NULL, boot_rep = NULL,
                      boot_idx = NULL,
                      np_boot = FALSE, b_boot = 1000, estimator, rotation, type,
                      rot_type = type,
                      n_factors, se = "none", ci = .95,
                      use = "pairwise.complete.obs", cor_method = "pearson",
                      max_iter = NA, init_comm = NA, criterion = NA,
                      criterion_type = NA, abs_eigen = NA, start_method = "psych",
                      normalize = TRUE, precision = 1e-5, order_type = NA,
                      varimax_type = NA, P_type = NA, k = NA, randomStarts = 100,
                      ...) {

  m <- ncol(R)
  df <- .efa_df(m, n_factors)

  if(df < 0){

    chisq_block <- .fit_unavailable_text("chisq_block")
    residuals_kept <- .fit_unavailable_text("residuals_kept")
    not_identification <- .fit_unavailable_text("residuals_not_identification")
    cli::cli_warn(
      c("The model is underidentified ({.code df = {df}}); {chisq_block}.",
        "!" = "{residuals_kept}",
        "!" = "{not_identification}",
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

  if (estimator != "PAF" && type == "SPSS") {
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
  } else if (estimator != "DWLS" && is.null(fiml) && se != "sandwich" &&
             df >= 0 && .bartlett_mult(N, m, n_factors) <= 0) {
    # A very small N turns Bartlett's multiplier non-positive, so .gof() leaves the chi-square
    # block undefined. Warn here, not in .gof(), which runs once per bootstrap replicate. The
    # FIML statistic has no Bartlett correction, and the sandwich refills the block with a
    # scaled statistic, so neither loses it.
    chisq_block <- .fit_unavailable_text("chisq_block")
    mult_note <- .fit_unavailable_text("bartlett_mult", N = N, p = m, q = n_factors)
    residuals_kept <- .fit_unavailable_text("residuals_kept")
    cli::cli_warn(
      c("{.arg N} is too small for the small-sample correction; {chisq_block}.",
        "i" = "{mult_note}",
        "i" = "{residuals_kept}"),
      class = "efa_fit_indices_undefined"
    )
  }

  fit_out <- .estimate_model(R, method = estimator, n_factors = n_factors, N = N,
                             type = type, max_iter = max_iter,
                             init_comm = init_comm, criterion = criterion,
                             criterion_type = criterion_type,
                             abs_eigen = abs_eigen, start_method = start_method,
                             weights = weights, fiml = fiml)

  # Heywood cases and non-convergence of the point estimate; the bootstrap replicates suppress
  # their own warnings (PAF warns about non-convergence from inside .PAF()).
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

  if (estimator != "PAF" && isTRUE(fit_out$convergence != 0)) {
    cli::cli_warn(
      c("The {.val {estimator}} optimiser did not converge (convergence code {fit_out$convergence}).",
        "i" = paste(.nonconvergence_cause(fit_out$convergence),
                    "The results may not be interpretable."),
        "i" = paste("Try extracting fewer factors, or check the correlation matrix for",
                    "near-collinear variables.")),
      class = "efa_nonconvergence"
    )
  }

  # Only the bootstrap produces replicate fits; the analytic SE methods leave this NULL.
  boot_fits <- NULL
  if (isTRUE(np_boot)) {
    # Each replicate fits only the quantities the bootstrap aggregation reads (lean = TRUE).
    # Its fit indices are not read: efa_fit() bootstraps them separately (.yhy_fit_boot()).
    boot_fits <- .boot_fun(boot_rep, boot_idx, .estimate_model,
                           method = estimator, n_factors = n_factors, N = N, type = type,
                           max_iter = max_iter, init_comm = init_comm, criterion = criterion,
                           criterion_type = criterion_type, abs_eigen = abs_eigen,
                           start_method = start_method, lean = TRUE, gof = FALSE)
  }

  rot_out <- NULL
  boot_rot <- "none"
  if (rotation != "none") {
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
  # The single-factor fallback of .rotate_model() returns NULL slots (Phi, Structure,
  # vars_accounted_rot); leave them out, as for an unrotated fit.
  drop_slots <- function(l) l[names(l) != "settings" & !vapply(l, is.null, logical(1))]
  output <- c(drop_slots(fit_out), drop_slots(rot_out),
              list(settings = c(settings_EFA, fit_out$settings, rot_out$settings)))

  # Stage-1 diagnostics of the FIML correlation, kept on the object because callers that
  # suppress per-fit warnings would otherwise lose the EM's non-convergence.
  if (!is.null(fiml)) {
    output$fiml <- list(converged = fiml$converged, iter = fiml$iter,
                        n_patterns = fiml$n_patterns, n = fiml$n)
  }

  # The unrotated loading covariance (filled below for the analytic SEs) and the covariance of
  # the correlations (sandwich only) stay present-but-NULL elsewhere, for efa_mi(); efa_fit()
  # drops `Gamma` unless efa_mi() asks for it.
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
    # Single-bracket assignment keeps a present-but-NULL slot where the method has none.
    output["replicates"] <- list(boot_out$replicates)
    output["vcov_unrot_loadings"] <- list(boot_out$vcov_unrot_loadings)
    # Only the bootstrap gives an SE for every residual. The diagonal (residual 0, SE 0)
    # is set to 0 instead of NaN.
    if (!is.null(boot_out$SE$residuals)) {
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



# Build the per-replicate input of the case bootstrap: a function of the resample's row indices
# that computes its correlation matrix (plus the DWLS weights). It runs inside
# the parallel workers, because polychoric and FIML correlations take much longer than the fit to
# them. It closes over the raw data only, so a parallel plan sends the data once; the indices
# travel as the elements of the parallel map, so each worker gets only its own share. A
# degenerate resample (a constant column, a failed polychoric or EM estimation) gives NULL, and
# the replicate is dropped.
.boot_rep_maker <- function(x, cor_method, use, dwls, fiml_max_iter, fiml_tol) {
  force(x); force(cor_method); force(use); force(dwls); force(fiml_max_iter); force(fiml_tol)
  vars <- colnames(x)
  function(ind) {
    xi <- x[ind, , drop = FALSE]
    rep_i <- tryCatch(suppressWarnings(
      if (.is_poly_cor(cor_method)) {
        poly <- .polychoric(xi, nearest_pd = FALSE, binary_only = cor_method == "tetra",
                            acov = if (dwls) "diag" else "none", label_acov = FALSE)
        list(R = poly$R, W = if (dwls) .poly_weight_matrix(poly$acov, ncol(xi)))
      } else if (cor_method == "fiml") {
        # The replicate fit indices are not used, so the EM moments that feed the two-stage
        # statistic are not passed on.
        em <- .fiml_em_moments(xi, max_iter = fiml_max_iter, tol = fiml_tol)
        list(R = stats::cov2cor(em$sigma))
      } else {
        list(R = stats::cor(xi, use = use, method = cor_method))
      }), error = function(e) NULL)
    if (!is.null(rep_i)) dimnames(rep_i$R) <- list(vars, vars)
    rep_i
  }
}

# Fit the bootstrap replicates, one per element of `units` (the input `boot_rep` takes), in
# parallel across replicates with future.apply; the default plan runs them sequentially.
# future.seed = TRUE gives each replicate its own L'Ecuyer-CMRG stream, so the result does not
# depend on the number of workers. A replicate that cannot be built or fitted is NULL and is
# skipped later. Per-replicate warnings would repeat the same information for every replicate,
# so they are suppressed; non-convergence is counted and reported once. `future_seed = FALSE`
# is for fits that draw no random numbers: it leaves the caller's random-number stream as it was,
# which future.seed = TRUE advances.
.boot_fun <- function(boot_rep, units, call_fun, ..., future_seed = TRUE) {
  boot_list <- future.apply::future_lapply(units, function(unit) {
    rep_i <- boot_rep(unit)
    if (is.null(rep_i)) return(NULL)
    tryCatch(suppressWarnings(call_fun(rep_i$R, ..., weights = rep_i$W)),
             error = function(e) NULL)
  }, future.seed = future_seed)

  n_nonconverged <- sum(vapply(boot_list,
                               function(fit_i) isTRUE(fit_i$convergence != 0),
                               logical(1)))

  if (n_nonconverged > 0L) {
    cli::cli_warn(
      c("{n_nonconverged} of {length(units)} bootstrap replicate{?s} did not converge.",
        "i" = "Their bootstrap standard errors and confidence intervals may be unreliable."),
      class = "efa_boot_nonconvergence"
    )
  }

  boot_list
}

# Bootstrap intervals of RMSEA, CFI, and TLI (Zhang & Savalei, 2016). A case bootstrap of the raw
# data adds sampling misfit to the misfit the sample already has, so its intervals of these
# indices are shifted toward worse fit. Here the cases are first transformed to the estimated
# misfit (.yhy_transform()) and then resampled with the same row indices as the case bootstrap.
# TLI is CFI rescaled by a constant, so its interval behaves in the same way (Zhang & Savalei,
# 2016, p. 394). The other fit indices get no interval: SRMR had poor coverage under both
# bootstraps, and Zhang and Savalei did not study the others. The replicate fits draw no random
# numbers, so they leave the caller's random-number stream as it was. `...` goes to
# .estimate_model(). Returns `out` unchanged where the transformation is undefined: with missing
# values in the analysed cases (it is linear in the cases), without a chi-square, or with df <= 0.
.yhy_fit_boot <- function(out, x, rows, boot_idx, use, ci, ...) {
  fit <- out$fit_indices
  if (!isTRUE(fit$df > 0) || !is.finite(fit$chi)) return(out)
  x_rows <- data.matrix(x[rows, , drop = FALSE])
  if (anyNA(x_rows)) return(out)
  # The model-implied matrix of the chi-square (.gof()): LL' with a unit diagonal.
  R_hat <- tcrossprod(unclass(out$unrot_loadings))
  diag(R_hat) <- 1
  yhy <- .yhy_transform(scale(x_rows), R_hat, out$settings$N, fit$df)
  if (is.null(yhy)) return(out)

  # The resample indices are row numbers of `x`, so the transformed cases keep those rows.
  x_yhy <- matrix(NA_real_, nrow(x), ncol(x), dimnames = list(NULL, colnames(x)))
  x_yhy[rows, ] <- yhy$data
  boot_rep <- .boot_rep_maker(x_yhy, cor_method = "pearson", use = use, dwls = FALSE,
                              fiml_max_iter = NA, fiml_tol = NA)
  # One warning about non-converged replicates is enough: these fits enter the percentiles in
  # the same way as those of the case bootstrap, which has already warned.
  fits <- suppressWarnings(
    .boot_fun(boot_rep, boot_idx, .yhy_rep_fit, ..., future_seed = FALSE),
    classes = "efa_boot_nonconvergence")

  keep <- c("RMSEA", "CFI", "TLI")
  reps <- t(vapply(fits, function(fit_i) {
    if (is.null(fit_i)) rep(NA_real_, length(keep)) else fit_i$fit_indices
  }, numeric(length(keep))))
  colnames(reps) <- keep
  l_ci <- (1 - ci) / 2
  agg <- .array_se_ci(reps, c(l_ci, ci + l_ci), M = 2)
  # Each slot goes before `residuals`, in the same order as the other bootstrap components.
  put <- function(l, value) {
    append(l, list(fit_indices = value), after = match("residuals", names(l)) - 1L)
  }
  out$SE <- put(out$SE, agg$se)
  out$CI <- put(out$CI, agg$ci)
  out$replicates <- put(out$replicates, reps)
  out
}

# One replicate of .yhy_fit_boot(): RMSEA, CFI, and TLI of a lean fit, and its convergence code.
# A top-level function, so the parallel workers do not receive the caller's data with it.
.yhy_rep_fit <- function(R, ...) {
  fit <- .estimate_model(R, ..., lean = TRUE)
  list(fit_indices = unlist(fit$fit_indices[c("RMSEA", "CFI", "TLI")]),
       convergence = fit$convergence)
}

# Yuan, Hayashi, and Yanagihara (2007) transformation of the standardized cases `z`, whose
# correlation matrix R is the one the model was fitted to. The transformed cases have the
# correlation matrix
# R_a = a R + (1 - a) R_hat, with R_hat the model-implied matrix (unit diagonal), so R_a also has
# a unit diagonal. `a` makes the ML discrepancy F(R_a, R_hat) equal the estimated misfit
# tau = max(0, F - df / (N - 1)), the noncentrality estimate behind the reported RMSEA and CFI
# (Zhang & Savalei, 2016, eqs. 4 and 11). With d the eigenvalues of R_hat^-1 R,
# F(R_a, R_hat) = g(a) = sum(a (d - 1) - log(1 + a (d - 1))), which increases strictly from
# g(0) = 0 to g(1) = F (YHY, Theorem 2), so a is in [0, 1]; tau = 0 gives a = 0, the Bollen-Stine
# population (YHY, p. 275). The fitted solution also fits R_a best, because the ML and
# least-squares estimating equations are linear in R (YHY, Theorem 1 and p. 275); this holds for
# PAF too, which converges to the ULS solution. The square roots are symmetric, as in lavaan's
# bootstrapLavaan(type = "yuan"). Returns the transformed cases and `a`, or NULL when R is
# numerically singular or R_hat is not positive definite.
.yhy_transform <- function(z, R_hat, N, df) {
  # R is taken from the cases, so that the transformed cases have exactly R_a. It needs an
  # inverse square root. The rank tolerance p * eps is the one of .prepare_cor_input(), which
  # explains it; the psych preset of PAF lets a singular R through (smoothed).
  R <- crossprod(z) / (nrow(z) - 1)
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) / max(ev) < ncol(R) * .Machine$double.eps) return(NULL)
  U <- tryCatch(chol(R_hat), error = function(e) NULL)
  if (is.null(U)) return(NULL)
  # R_hat = U'U, so U^-T R U^-1 has the eigenvalues of R_hat^-1 R.
  E <- backsolve(U, t(backsolve(U, R, transpose = TRUE)), transpose = TRUE)
  d <- eigen(E, symmetric = TRUE, only.values = TRUE)$values
  g <- function(a) sum(a * (d - 1) - log1p(a * (d - 1)))
  F_hat <- g(1)
  if (!is.finite(F_hat)) return(NULL)
  tau <- max(0, F_hat - df / (N - 1))
  # g(0) - tau = -tau < 0 and g(1) - tau = df / (N - 1) > 0, so the root is bracketed.
  a <- if (tau > 0) {
    stats::uniroot(function(a) g(a) - tau, c(0, 1), tol = 1e-10)$root
  } else {
    0
  }
  R_a <- a * R + (1 - a) * R_hat
  list(data = z %*% .inv_mat_sqrt(R) %*% .mat_sqrt(R_a), a = a)
}
