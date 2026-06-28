#' Exploratory factor analysis (EFA)
#'
#' This function does an EFA with either `PAF`, `ML`, `ULS`/`MINRES`,
#' or `DWLS` with or without subsequent rotation.
#' All arguments with default value `NA` can be left to default if `type`
#' is set to one of "EFAtools", "SPSS", or "psych". The respective specifications are
#' then handled according to the specified type (see details). All rotations are
#' performed by rotation engines built into the package.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#' correlations. If raw data is entered, the correlation matrix is found from the
#' data.
#' @param n_factors numeric. Number of factors to extract.
#' @param N numeric. The number of observations. Needs only be specified if a
#' correlation matrix is used. If input is a correlation matrix and `N` = NA
#' (default), not all fit indices can be computed. When raw data with missing
#' values are entered and `use` is `"complete.obs"` or `"na.or.complete"`, rows
#' are deleted listwise, so `N` is taken as the number of complete cases.
#' @param method character. The estimator used to fit the EFA: "PAF" (principal axis
#' factoring), "ML" (maximum likelihood), "ULS" (unweighted least squares; "MINRES" is an
#' accepted alias returning identical results), or "DWLS" (diagonally weighted least
#' squares, for ordinal data). See the *Estimators* section in Details for their
#' properties and data requirements.
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
#' @param type character. If one of "EFAtools" (default), "psych", or "SPSS" is
#'  used, and the following arguments with default NA are left with
#'  NA, these implementations are executed according to the respective program
#'  ("psych" and "SPSS") or according to the best solution found in Grieder &
#'  Steiner (2022; "EFAtools"). Individual properties can be adapted using one of
#'  the three types and specifying some of the following arguments. If set to
#'  "none" additional arguments must be specified depending on the `method`
#'  and `rotation` used (see details).
#' @param max_iter numeric. The maximum number of iterations to perform after which
#' the iterative PAF procedure is halted with a warning. If `type` is one of
#' "EFAtools", "SPSS", or "psych", this is automatically specified if `max_iter` is
#' left to be `NA`, but can be overridden by entering a number. Default is
#' `NA`.
#' @param init_comm character. The method to estimate the initial communalities
#' in `PAF`. "smc" will use squared multiple correlations, "mac" will use
#' maximum absolute correlations, "unity" will use 1s (see details).
#' Default is `NA`.
#' @param criterion numeric. The convergence criterion used for PAF.
#' If the change in communalities from one iteration to the next is smaller than
#' this criterion the solution is accepted and the procedure ends.
#' Default is `NA`.
#' @param criterion_type character. Type of convergence criterion used for
#' PAF. "max_individual" selects the maximum change in any of the
#' communalities from one iteration to the next and tests it against the
#' specified criterion. This is also used by SPSS. "sum" takes the difference of
#' the sum of all communalities in one iteration and the sum of all communalities
#' in the next iteration and tests this against the criterion. This procedure is
#' used by the [psych::fa()] function. Default is `NA`.
#' @param abs_eigen logical. Which algorithm to use in the PAF
#' iterations. If FALSE, the loadings are computed from the eigenvalues. This is
#' also used by the [psych::fa()] function. If TRUE the
#' loadings are computed with the absolute eigenvalues as done by SPSS.
#' Default is `NA`.
#' @param use character. Passed to [stats::cor()] if raw data
#' is given as input. Default is "pairwise.complete.obs".
#' @param cor_method character. How the correlation is computed from raw data:
#'   `"pearson"`, `"spearman"`, or `"kendall"` (passed to [stats::cor()]); `"poly"` /
#'   `"tetra"` for polychoric / tetrachoric correlations of ordinal / binary data; or
#'   `"fiml"` for a two-stage full-information maximum-likelihood correlation from raw data
#'   with missing values. See the *Correlation methods* section in Details for their
#'   properties and the combinations they support. Default is "pearson".
#' @param k numeric. Either the power used for computing the target matrix P in
#' the promax rotation or the number of 'close to zero loadings' for the simplimax
#' rotation. If left to `NA` (default), the value for promax depends on the specified type.
#' For simplimax, `nrow(L)`, where L is the matrix of unrotated loadings,
#' is used by default.
#' @param normalize logical. If `TRUE`, a kaiser normalization is
#' performed before the specified rotation. Default is `TRUE`.
#' @param P_type character. This specifies how the target
#' matrix P is computed in promax rotation. If "unnorm" it will use the
#' unnormalized target matrix as originally done in Hendrickson and White (1964).
#' This is also used in the psych and stats packages. If "norm" it will use the
#' normalized target matrix as used in SPSS. Default is `NA`.
#' @param precision numeric. The tolerance for stopping in the rotation
#' procedure. Default is 10^-5 for all rotation methods.
#' @param varimax_type character. The type of the varimax rotation performed.
#' If "svd", singular value decomposition is used, as [stats::varimax()] does. If
#' "kaiser", the varimax procedure performed in SPSS is used, following the original
#' procedure from Kaiser (1958) (see details). Default is `NA`.
#' @param order_type character. How to order the factors. "eigen" reorders the
#' factors by descending explained variance, i.e. by their reported sums of squared
#' loadings ("SS loadings"): the column sums of squares for orthogonal solutions and
#' the factor-intercorrelation-weighted sums of squares for oblique solutions, so the
#' reported variances decrease monotonically (as in the psych package). "ss_factors"
#' reorders the factors by descending (unweighted) sum of squared factor loadings per
#' factor; for oblique solutions this can differ from "eigen", whereas for orthogonal
#' solutions the two coincide. Default is `NA`.
#' @param start_method character. How to specify the starting values for the
#' optimization procedure for ML. Default is "psych" which takes the
#' starting values specified in [psych::fa()]. "factanal" takes the
#' starting values specified in the [stats::factanal()] function.
#' Solutions are very similar.
#' @param b_boot numeric. The number of bootstrap samples to draw. Default is 1000.
#'  Under `cor_method = "fiml"` each bootstrap sample re-runs the EM moment
#'  estimation, so a smaller value may be advisable.
#' @param ci numeric. The confidence interval to create from the bootstrap samples.
#'  Must be between 0 and 1. Default ist .95 for 95% CIs.
#' @param randomStarts numeric. The number of random starts to use in the
#'  rotation. Some rotation criteria are prone to produce local minima, and
#'  several random starts are usually needed to locate the best solution. The
#'  rotation screens the random starts cheaply and fully optimises only the most
#'  promising ones, so a large value adds little cost for most criteria. The
#'  complexity criteria (simplimax and, to a lesser extent, geomin) are the most
#'  multimodal and may need a larger value on difficult data. Default is 100.
#' @param seed numeric. An optional seed for the random-number generator used by the
#'  non-parametric bootstrap (`se = "np-boot"`), i.e. for the case resampling, the
#'  rotation random starts, and the Procrustes random starts. Setting it makes the
#'  bootstrap reproducible and, importantly, independent of the number of parallel
#'  workers (see Details); the caller's random-number stream is restored afterwards,
#'  so supplying a seed leaves no lasting effect on it. Default is `NULL`, which uses
#'  (and advances) the current state of the generator.
#' @param ... Additional arguments passed to the rotation procedure (e.g., `maxit` for the maximum number of iterations).
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
#' The estimator is chosen with `method`.
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
#'   two-step estimator with no empty-cell continuity correction, matching
#'   `polycor::polychor()` and `lavaan`. The polychoric asymptotic covariance that
#'   underlies both the DWLS weights and the scaled (sandwich) statistic relies on
#'   large-sample theory that degrades for empty or near-empty response-category
#'   combinations; with very sparse cells the resulting weights and standard errors can be
#'   unreliable (a warning is issued when empty cells are present), so interpret them with
#'   caution and consider collapsing rare categories.
#' - **"fiml"** estimates a two-stage full-information maximum-likelihood correlation. The
#'   saturated multivariate-normal mean and covariance are estimated from raw data with
#'   missing values by an EM algorithm assuming the data are missing at random (Yuan,
#'   Marshall, & Bentler, 2002; Little & Rubin, 2002), and the standardized covariance is
#'   then analysed. This reproduces `psych::corFiml()` followed by `psych::fa()` and
#'   `lavaan(missing = "two.stage")`, *not* `lavaan::efa(missing = "ml")`, so the point
#'   estimates are not expected to match the latter. The model fit indices are corrected
#'   two-stage statistics (see *Fit indices*). `"fiml"` uses every case and handles the
#'   missingness itself, so `use` is ignored; it supplies a continuous (Pearson-type)
#'   correlation only and is therefore not compatible with `method = "DWLS"`. Standard
#'   errors are available analytically for `method = "ML"` or `"ULS"` and, for any method,
#'   by the non-parametric bootstrap (see *Standard errors*). For multiply imputed data,
#'   [EFA_POOLED()] is the alternative route to handling missingness.
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
#'   (controlled by `k` and `P_type`) to form a target that is then fitted obliquely. It is
#'   the common, inexpensive oblique default.
#' - **oblimin** is a flexible oblique family controlled by `gam` (default 0); a good
#'   general-purpose criterion.
#' - **quartimin** is oblimin pinned at `gam = 0`; a robust default oblique criterion.
#' - **simplimax** drives the `k` smallest loadings toward zero. Its criterion is only
#'   piecewise smooth, so it is the most prone to local minima and relies on several
#'   `randomStarts`.
#' - **bentlerQ** is the oblique Bentler invariant pattern simplicity criterion.
#' - **geominQ** is the oblique geomin criterion; it handles complex (cross-loading)
#'   structure well but is multimodal, so it benefits from more `randomStarts` (and uses a
#'   more thorough multi-start search internally).
#' - **bifactorQ** is the oblique (correlated) Jennrich-Bentler bifactor criterion.
#'
#' The criterion-based rotations (all except varimax and promax) are fitted by gradient
#' projection with `randomStarts` random starts to guard against local minima; the
#' complexity criteria (simplimax and geominQ in particular) are the most multimodal. The
#' `type` argument changes the varimax and promax settings (see *Using the type presets*)
#' and, for every rotation, the factor `order_type`. A single factor cannot be rotated.
#'
#' ## Standard errors
#'
#' `se` selects whether and how standard errors (and matching confidence intervals) are
#' computed. They cover the unrotated loadings and uniquenesses and, when a rotation is
#' applied, the rotated loadings, the communalities, and -- for oblique rotations -- the
#' factor correlations and the structure coefficients (see the `SE` and `CI` entries in
#' Value).
#'
#' - **"none"** (default) computes no standard errors.
#' - **"information"** returns analytic standard errors from the expected (Fisher)
#'   information matrix of the maximum-likelihood solution, and therefore requires
#'   `method = "ML"`. The rotated standard errors are obtained by propagating the
#'   unrotated-loading covariance through the rotation by the delta method (Jennrich,
#'   1973); because rotated quantities are identification-invariant they are directly
#'   comparable across programs. Unlike the bootstrap it also works from a correlation
#'   matrix as long as `N` is supplied. The covariance is the inverse expected information
#'   under the identification constraint that \eqn{\Lambda' \Psi^{-1} \Lambda} is diagonal,
#'   scaled by \eqn{1 / (N - 1)}; the confidence intervals are Wald intervals (estimate
#'   \eqn{\pm} z * SE). These standard errors assume multivariate normality and a correctly
#'   specified model; under heavy-tailed data or model misfit they can understate the
#'   sampling variability, where a bootstrap is more robust. The rotated
#'   structure-coefficient intervals are somewhat conservative for
#'   high-communality variables, where `"sandwich"` or `"np-boot"` give sharper intervals.
#' - **"sandwich"** returns robust (Godambe sandwich) standard errors from raw data,
#'   combining the estimator weight with an asymptotic-distribution-free covariance of the
#'   correlations, so it stays valid under non-normality and weight misspecification
#'   (Browne, 1984; Satorra & Bentler, 1994). It is available either for ordinal data with
#'   `cor_method = "poly"` or `"tetra"` and `method` one of `"ML"`, `"ULS"`, or `"DWLS"`
#'   (the meat is the polychoric / tetrachoric asymptotic covariance), or for continuous
#'   data with `cor_method = "pearson"` and `method = "ML"` or `"ULS"` (the meat is the
#'   fourth-moment ADF covariance of the sample correlations, the basis of the MLM / MLR
#'   robust statistics). It reports the same standard errors as `"information"`, propagated
#'   by the same delta method, and additionally fills the model fit's chi-square block with
#'   a scaled (Satorra-Bentler / scaled-and-shifted) chi-square (see *Fit indices*).
#'   Because the asymptotic covariance must describe the same cases as the correlation
#'   matrix, the sandwich (like `method = "DWLS"`) is computed on the listwise-complete
#'   cases; on data with missing values the reported `N`, the correlation matrix, and the
#'   point estimate therefore reflect the complete cases regardless of `use`.
#' - **"np-boot"** draws a non-parametric (case-resampling) bootstrap and needs raw data.
#'   It is the most general method -- available for any `method`, rotation, and
#'   `cor_method` -- and the most robust to non-normality and misfit, at the cost of speed;
#'   its intervals are bootstrap percentile intervals. The replicate fits are run across
#'   replicates with the `future` framework. By default they run sequentially; to run them
#'   in parallel, register a plan with [future::plan()] (e.g.
#'   `future::plan(future::multisession, workers = 2)`; see examples). With a fixed `seed`
#'   the bootstrap is reproducible and yields the same result regardless of the number of
#'   workers. Each worker runs its own (Armadillo) linear algebra, so if your `BLAS` is
#'   multi-threaded, limit the number of workers (or the BLAS threads) to avoid
#'   over-subscribing the available cores. Under `cor_method = "fiml"` each resample also
#'   re-runs the EM moment estimation and is therefore slow, so a smaller `b_boot` may be
#'   advisable.
#'
#' The analytic methods (`"information"` and `"sandwich"`) are not available with the
#' `"promax"` or `"simplimax"` rotations, which have no usable analytic rotation Jacobian;
#' use `"np-boot"` there. Under `cor_method = "fiml"`, `"information"` and `"sandwich"`
#' instead return, for `method = "ML"` or `"ULS"`, the corrected two-stage (Yuan & Bentler,
#' 2000; Savalei & Bentler, 2009) sandwich standard errors, built on the saturated FIML
#' asymptotic covariance with the estimator's own Stage-2 weight: the model is fitted to
#' the EM-estimated correlation, so the naive Stage-2 standard errors (treating that
#' correlation as complete data) are inconsistent under missingness and are never reported
#' (`method = "PAF"` carries no Stage-2 weight, so use `se = "np-boot"` there).
#'
#' See Lawley and Maxwell (1971) and Jennrich and Thayer (1973) for the information-matrix
#' standard errors; Jennrich (1973) and Zhang and Preacher (2015) for the rotated
#' quantities; Browne (1984), Satorra and Bentler (1994), and Asparouhov and Muthén (2010)
#' for the robust standard errors and scaled chi-square; and Yuan and Bentler (2000),
#' Savalei and Bentler (2009), and Yuan, Marshall, and Bentler (2002) for the two-stage
#' FIML standard errors and rescaled statistic.
#'
#' ## Fit indices
#'
#' For ML and ULS, [EFA()] reports the model chi-square (with its p-value and degrees of
#' freedom), the Comparative Fit Index (CFI; Bentler, 1990), the Tucker-Lewis Index (TLI,
#' also called the non-normed fit index; Tucker & Lewis, 1973), the Root Mean Square Error
#' of Approximation (RMSEA) with its 90% confidence interval (Browne & Cudeck, 1992), the
#' Akaike and Bayesian Information Criteria (AIC, BIC), the Expected Cross-Validation Index
#' (ECVI; Browne & Cudeck, 1989), the Root Mean Squared Residual (RMSR), the Standardized
#' Root Mean Squared Residual (SRMR; Bentler, 1995), and the common-part-accounted-for
#' (CAF) index (Lorenzo-Seva, Timmerman, & Kiers, 2011). The model chi-square is the
#' Bartlett-corrected discrepancy (matching [stats::factanal()] for ML); the AIC, BIC, and
#' ECVI are the minimum-fit-function (chi-square-based) forms (as in [psych::fa()]) and can
#' therefore be negative. The RMSEA, CFI, and TLI place the model and baseline
#' noncentrality on the uncorrected \eqn{N - 1} discrepancy scale on which these
#' approximate-fit indices are defined, so the Bartlett small-sample correction enters only
#' the chi-square test, not the approximate-fit indices.
#'
#' Which indices are reported depends on the estimator:
#' - **ML and ULS** report the full set above.
#' - **PAF** reports only the descriptive residual indices (RMSR, SRMR, CAF) and df; the
#'   chi-square-based indices are `NA`, because PAF minimises no discrepancy.
#' - **DWLS** by default reports only RMSR, SRMR, CAF, and df, because the ordinary
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
#' with values close to 1 indicating close fit; this does not match the formula they give
#' for it, \eqn{1 - KMO(residuals)}, which only works if the diagonal of the residual
#' matrix is set to 1s and then approximates 0.5 with close fit.
#'
#' ## Available combinations
#'
#' Not every estimator, rotation, standard-error, and correlation method can be combined:
#'
#' - **Estimator and correlation method.** `method = "DWLS"` requires ordinal data with
#'   `cor_method = "poly"` or `"tetra"`. `cor_method = "fiml"` works with PAF, ML, and ULS
#'   (not DWLS) and needs raw data with missing values.
#' - **Standard errors.** `se = "information"` requires `method = "ML"` and can be computed
#'   from a correlation matrix when `N` is supplied. `se = "sandwich"` requires raw data,
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
#' For promax, the values of `P_type`,
#' `order_type`, and `k` depend on the `type` argument.
#'
#' `type = "EFAtools"` will use the following argument specification:
#' `P_type = "norm", order_type = "eigen", k = 4`.
#'
#' `type = "psych"` will use the following argument specification:
#' `P_type = "unnorm", order_type = "eigen", k = 4`.
#'
#' `type = "SPSS"` will use the following argument specification:
#' `P_type = "norm", order_type = "ss_factors", k = 4`.
#'
#' The `P_type` argument can take two values, "unnorm" and "norm". It controls
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
#' reproduce the [psych::fa()] promax result exactly.
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
#' maximum number of iterations).
#'
#' The `type` argument has no effect on ULS and ML. For ULS, no additional
#' arguments are needed. For ML, an additional argument
#' `start_method` is needed to determine the starting values for the
#' optimization procedure. Default for this argument is "psych" which takes
#' the starting values specified in [psych::fa()].
#'
#'
#' @return A list of class EFA containing (a subset of) the following:
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
#' \item{convergence}{Integer convergence code (0 = converged). For ML, ULS, and
#'  DWLS this is the convergence code from the bounded optimiser (the same codes as
#'  [`stats::optim()`][stats::optim]'s "L-BFGS-B"); for PAF it is 1 if the maximum
#'  number of iterations was reached without meeting the convergence criterion and 0
#'  otherwise. A non-zero code is also reported with a warning.}
#' \item{heywood}{A named integer vector indicating which variables have a
#'  Heywood (improper) case in the unrotated solution; empty if there are none.}
#' \item{unrot_loadings}{Loading matrix containing the final unrotated loadings.}
#' \item{vars_accounted}{Matrix of explained variances and sums of squared loadings. Based on the unrotated loadings.}
#' \item{fit_indices}{For ML and ULS: Fit indices derived from the unrotated
#' factor loadings: Chi Square, including significance level, degrees of freedom
#' (df), Comparative Fit Index (CFI; Bentler, 1990), Tucker-Lewis Index (TLI, also
#' called the non-normed fit index; Tucker & Lewis, 1973), Root Mean Square Error of
#' Approximation (RMSEA), including its 90% confidence interval (Browne & Cudeck,
#' 1992), Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC),
#' Expected Cross-Validation Index (ECVI; Browne & Cudeck, 1989), Root Mean Squared
#' Residual (RMSR), Standardized Root Mean Squared Residual (SRMR; Bentler, 1995),
#' and the common part accounted for (CAF) index as proposed by Lorenzo-Seva,
#' Timmerman, & Kiers (2011). The model Chi Square is the Bartlett-corrected
#' discrepancy (matching [stats::factanal()] for ML), and the AIC, BIC, and ECVI are
#' derived from it. When `cor_method = "fiml"` the model and baseline Chi Square are
#' instead the Satorra-Bentler-corrected two-stage statistics (Yuan, Marshall, &
#' Bentler, 2002): the normal-theory discrepancy on the EM-estimated correlation,
#' rescaled by the saturated FIML asymptotic covariance, because the plain two-stage
#' likelihood-ratio statistic is not asymptotically `Chi^2(df)` under the two-stage
#' estimator. For ML and ULS the AIC, BIC, and ECVI are the minimum-fit-function
#' (Chi-Square-based) forms `Chi Square - 2 * df` and `Chi Square - log(N) * df` (as in
#' [psych::fa()]), not the likelihood-based criteria, and can therefore be negative; they
#' are left `NA` for any scaled (moment-adjusted) Chi Square, i.e. under
#' `cor_method = "fiml"` and with `se = "sandwich"`. The RMSEA, CFI, and
#' TLI instead place the model (and, for CFI and TLI, the baseline) noncentrality on the
#' uncorrected `N - 1` discrepancy scale (Browne & Cudeck, 1992; Bentler, 1990; Tucker &
#' Lewis, 1973) on which these approximation indices are defined, so the Bartlett
#' small-sample correction enters only the Chi Square test, not the approximate-fit
#' indices. For PAF
#' and DWLS, only CAF, RMSR, SRMR, and df are computed and the Chi-Square-derived indices
#' are returned as `NA`; for DWLS with `se = "sandwich"` the full block is instead filled
#' from a scaled Chi Square (see Details). Whenever the Chi Square is a scaled one
#' (`se = "sandwich"`, or any `cor_method = "fiml"` fit) the `fit_indices` also
#' carry the scaled-statistic components: `chi_scaling` (the multiplier a in the
#' scaled-and-shifted statistic a * T + b, i.e. the reciprocal of \pkg{lavaan}'s
#' `chisq.scaling.factor`), `chi_shift` (b), `chi_unscaled` (the unscaled statistic T),
#' and the alternative `chi_mean_adjusted` and `chi_mean_var` statistics with their
#' `df_mean_var`. Note that Lorenzo-Seva, Timmerman, & Kiers (2011)
#' introduce the CAF as ranging between 0 and 1, with values close to 1 indicating close fit.
#' This does not match the formula they give for it, `1 - KMO(residuals)`, which only works if the
#' diagonal of the residual matrix is set to 1s and then approximates 0.5 with close fit.}
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
#' \item{settings}{A list of the settings used.}
#' \item{SE}{A named list of standard error matrices. For `se = "np-boot"`: bootstrap standard deviations of the unrotated and (when a rotation is applied) rotated loadings, the residuals, and the fit indices, plus -- for oblique rotations -- the factor correlations (`Phi`) and the structure coefficients. For `se = "information"`: Wald standard errors from the expected (Fisher) information matrix for the unrotated loadings and the uniquenesses and, when a rotation is applied, the rotated loadings and the communalities (and, for oblique rotations, `Phi` and the structure coefficients). For `se = "sandwich"`: robust Godambe sandwich standard errors with the same coverage as `"information"`, robust to non-normality and weight misspecification. Only returned if `se` is not `"none"`.}
#' \item{CI}{A named list of confidence intervals of width `ci`. For `se = "np-boot"`: percentile intervals matching the components of `SE`. For `se = "information"` and `se = "sandwich"`: Wald intervals matching the components of `SE`. Only returned if `se` is not `"none"`.}
#' \item{replicates}{A named list of bootstrap replicate cubes for the aligned unrotated and (where applicable) rotated loadings, structure coefficients, factor correlations (`Phi`), residuals, and fit indices. Each cube's last dimension indexes the replicate. Populated only for `se = "np-boot"`; `NULL` for the analytic SE methods.}
#' \item{vcov_unrot_loadings}{The full unrotated loading covariance matrix the marginal `SE$unrot_loadings` were derived from: a `p * n_factors` by `p * n_factors` numeric matrix in column-major `vec(Lambda)` order. Populated for `se = "information"` (expected-information block) and `se = "sandwich"` (robust V_AA), even when a rotation is applied (the persisted block is always the unrotated one); NA-filled if the analytic covariance is unreliable (a Heywood case or a singular bordered information matrix); `NULL` for `se = "np-boot"` and `se = "none"`.}
#' \item{Gamma}{The asymptotic covariance of the off-diagonal sample correlations -- the meat of the robust sandwich SEs -- on the variance scale (`Var(rho-hat)`; lavaan's correlation NACOV is `N * Gamma`). A `p (p - 1) / 2` by `p (p - 1) / 2` numeric matrix; rows and columns ordered by [utils::combn()] over the column pairs and labelled `"<var_i>-<var_j>"`. Populated only for `se = "sandwich"`; `NULL` otherwise.}
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
#' @source Jennrich, R. I., & Thayer, D. T. (1973). A note on Lawley's formulas for
#' standard errors in maximum likelihood factor analysis. Psychometrika, 38, 571–580.
#' doi: 10.1007/BF02291495
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
#' @export
#'
#' @examples
#'
#' # Principal axis factoring with oblimin rotation
#' mod_oblimin <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                    rotation = "oblimin")
#' mod_oblimin
#' summary(mod_oblimin)
#'
#' # ML estimation with oblimin rotation
#' mod_oblimin <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                    method = "ML", rotation = "oblimin")
#' mod_oblimin
#' summary(mod_oblimin)
#'
#' # Analytic (expected-information) standard errors for the above
#' ML_info <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                method = "ML", rotation = "oblimin", se = "information")
#' ML_info
#' summary(ML_info)
#'
#' \donttest{
#' # Robust (sandwich) standard errors and a scaled chi-square for ordinal raw data.
#' # These need a polychoric/tetrachoric correlation method and method ML, ULS, or DWLS.
#' DWLS_rob <- EFA(DOSPERT_raw, n_factors = 6, cor_method = "poly",
#'                 method = "DWLS", rotation = "oblimin", se = "sandwich")
#' DWLS_rob
#' summary(DWLS_rob)
#'
#' # The same robust SEs and scaled chi-square for continuous data: a Pearson
#' # correlation with method ML or ULS (the fourth-moment ADF covariance).
#' ML_rob <- EFA(GRiPS_raw, n_factors = 1, cor_method = "pearson",
#'               method = "ML", rotation = "none", se = "sandwich")
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
#' EFA_fiml <- EFA(x_miss, n_factors = 1, method = "ML", cor_method = "fiml")
#' EFA_fiml
#' }
#'
#' \dontrun{
#' # Bootstrap standard errors from raw data, reproducible via a fixed seed and run
#' # in parallel across replicates.
#' future::plan(future::multisession, workers = 2)
#' EFA_boot <- EFA(GRiPS_raw, n_factors = 1, method = "PAF", rotation = "none",
#'                 se = "np-boot", b_boot = 1000, seed = 42)
#' future::plan(future::sequential)
#' }
#'
EFA <- function(x, n_factors, N = NA, method = c("PAF", "ML", "ULS", "MINRES", "DWLS"),
                rotation = c("none", "varimax", "equamax", "quartimax", "geominT",
                             "bentlerT", "bifactorT", "promax", "oblimin",
                             "quartimin", "simplimax", "bentlerQ", "geominQ",
                             "bifactorQ"),
                se = c("none", "information", "sandwich", "np-boot"),
                type = c("EFAtools", "psych", "SPSS", "none"), max_iter = NA,
                init_comm = NA, criterion = NA, criterion_type = NA,
                abs_eigen = NA, use = c("pairwise.complete.obs", "all.obs",
                                          "complete.obs", "everything",
                                          "na.or.complete"),
                varimax_type = NA,
                k = NA, normalize = TRUE, P_type = NA, precision = 1e-5,
                order_type = NA, start_method = "psych",
                cor_method = c("pearson", "spearman", "kendall", "poly", "tetra",
                               "fiml"),
                b_boot = 1000, ci = .95,
                randomStarts = 100, seed = NULL,
                ...) {

  # Perform argument checks
  .assert_cor_input(x)

  method <- match.arg(method)
  # "MINRES" is a synonym for "ULS" (same estimator); resolve it once here so the
  # rest of EFA() and the reported settings use the single canonical name.
  if (method == "MINRES") method <- "ULS"
  rotation <- match.arg(rotation)
  se <- match.arg(se)
  np_boot <- se == "np-boot"
  use <- match.arg(use)
  cor_method <- match.arg(cor_method)
  type <- match.arg(type)
  start_method <- checkmate::matchArg(start_method, c("psych", "factanal", NA))

  if (is.na(start_method) && method == "ML") {
    cli::cli_abort(
      c("{.arg start_method} must be set when {.code method = \"ML\"}.",
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
    if (method == "DWLS") {
      cli::cli_abort(
        c("{.code method = \"DWLS\"} is not compatible with {.code cor_method = \"fiml\"}.",
          "x" = "DWLS needs a polychoric asymptotic covariance, which the continuous FIML correlation does not provide.",
          "i" = "Use {.code method = \"ML\"}, {.val ULS}, or {.val PAF}, or {.code cor_method = \"poly\"}/{.val tetra} for DWLS."),
        class = "efa_fiml_unsupported_method"
      )
    }
    # The corrected two-stage sandwich (.se_fiml()) reuses the Stage-2 ML/ULS weight; PAF minimises
    # no discrepancy and so carries no weight for it, exactly as the polychoric/continuous sandwich
    # rejects PAF. Reject the analytic-SE request here (the bootstrap stays available for PAF).
    if (se %in% c("information", "sandwich") && !(method %in% c("ML", "ULS"))) {
      cli::cli_abort(
        c("Analytic standard errors under {.code cor_method = \"fiml\"} require {.code method = \"ML\"} or {.val ULS}.",
          "x" = "You requested {.code method = {.val {method}}}.",
          "i" = "{.val PAF} minimises no discrepancy, so it has no weight for the two-stage sandwich; use {.code se = \"np-boot\"} for {.val {method}}."),
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
  # method = "ML"; promax (a two-step target rotation) and simplimax (a non-smooth, piecewise
  # criterion) have no usable analytic rotation Jacobian; PAF has no discrepancy-based information
  # from which a sandwich could be built.
  if (se == "information" && method != "ML" && cor_method != "fiml") {
    cli::cli_abort(
      c("{.code se = \"information\"} is only available for {.code method = \"ML\"}.",
        "x" = "You requested {.code method = {.val {method}}}.",
        "i" = "Use {.code se = \"np-boot\"} for {.val {method}}."),
      class = "efa_se_unsupported"
    )
  }
  if (se == "sandwich" && method == "PAF" && cor_method != "fiml") {
    cli::cli_abort(
      c("{.code se = \"sandwich\"} is not available for {.code method = \"PAF\"}.",
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
  # (Browne, 1984) covariance for continuous data with cor_method = "pearson" and method ML or
  # ULS. Spearman/Kendall correlations and continuous DWLS have no such covariance.
  if (se == "sandwich" && cor_method != "fiml" && !.is_poly_cor(cor_method) &&
      !(cor_method == "pearson" && method %in% c("ML", "ULS"))) {
    cli::cli_abort(
      c("{.code se = \"sandwich\"} is not available for this correlation/estimator combination.",
        "x" = "You requested {.code cor_method = {.val {cor_method}}} with {.code method = {.val {method}}}.",
        "i" = "Use {.code cor_method = \"poly\"}/{.val tetra} (any estimator), or {.code cor_method = \"pearson\"} with {.code method = \"ML\"} or {.val ULS}; otherwise use {.code se = \"np-boot\"}."),
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
  checkmate::assert_choice(P_type, c("unnorm", "norm", NA))
  checkmate::assert_number(precision, lower = 0, upper = 1)
  checkmate::assert_choice(order_type, c("eigen", "ss_factors", NA))
  checkmate::assert_integerish(b_boot, len = 1, any.missing = FALSE)
  checkmate::assert_number(ci, lower = 0, upper = 1)
  checkmate::assert_int(seed, null.ok = TRUE)

  # The common-factor model requires fewer factors than variables; with
  # n_factors >= n_variables it is not identified and the eigenvalue-based
  # extraction in the ML, ULS, and PAF fitters reads past the available
  # eigenvalues (undefined behaviour in an unchecked build).
  n_vars <- ncol(x)
  if (n_factors >= n_vars) {
    cli::cli_abort(
      c("{.arg n_factors} must be smaller than the number of variables.",
        "x" = "You requested {n_factors} factor{?s} for {n_vars} variable{?s}.",
        "i" = "Extract fewer factors."),
      class = "efa_too_many_factors"
    )
  }

  # DWLS weights each polychoric correlation residual by the inverse of its asymptotic
  # variance, so it needs a polychoric/tetrachoric asymptotic covariance. That is only
  # available from raw ordinal data with cor_method = "poly" or "tetra"; there is no
  # fallback to unit weights. Resolve this before any computation so an unsupported
  # request fails with a single clear error rather than downstream.
  if (method == "DWLS" && !(.is_poly_cor(cor_method) && !is_cormat)) {
    cli::cli_abort(
      c("{.code method = \"DWLS\"} requires a polychoric asymptotic covariance.",
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
      "i" = "Setting {.var se} to {.val none}. Rerun with raw data to calculate bootstrap SEs.")
    )

    np_boot <- FALSE
    se <- "none"

  }

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
                                    else if (method == "DWLS") "diag" else "none",
                             dwls = method == "DWLS",
                             check_singular = type != "psych",
                             posdef_abort = type == "SPSS")
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

    # A fixed seed makes the whole bootstrap reproducible and independent of the
    # number of parallel workers: the case resampling, the rotation random starts,
    # and the Procrustes random starts all run from this state. The case resampling
    # advances the global RNG by a b_boot-dependent amount; the parallel replicate
    # fit then adds a fixed, worker-count-independent step when future.seed = TRUE
    # derives a per-replicate L'Ecuyer stream. Both advances are deterministic given
    # the seed, so the downstream draws -- and the result -- are identical at any
    # number of workers. The caller's RNG stream is saved and restored afterwards --
    # or, if none existed yet, the seed set.seed() creates is removed again -- so
    # EFA() leaves no side effect on it.
    if (!is.null(seed)) {
      seed_existed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
      saved_seed <- if (seed_existed) {
        get(".Random.seed", envir = globalenv(), inherits = FALSE)
      }
      on.exit({
        if (seed_existed) {
          assign(".Random.seed", saved_seed, envir = globalenv())
        } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
          rm(".Random.seed", envir = globalenv())
        }
      }, add = TRUE)
      set.seed(seed)
    }

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
    } else if (.is_listwise_use(use) || method == "DWLS") {
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
    dwls <- method == "DWLS"

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

  .efa_core(
    R = R, N = N, weights = weights, Gamma = Gamma,
    R_boot_array = if (isTRUE(np_boot)) R_boot_array else NULL,
    W_boot_array = if (isTRUE(np_boot)) W_boot_array else NULL,
    fiml = fiml_pt, fiml_boot = if (isTRUE(np_boot)) fiml_boot else NULL,
    np_boot = np_boot, b_boot = b_boot, method = method, rotation = rotation,
    type = type, n_factors = n_factors, se = se, ci = ci, use = use,
    cor_method = cor_method, max_iter = max_iter, init_comm = init_comm,
    criterion = criterion, criterion_type = criterion_type,
    abs_eigen = abs_eigen, start_method = start_method, normalize = normalize,
    precision = precision, order_type = order_type, varimax_type = varimax_type,
    P_type = P_type, k = k, randomStarts = randomStarts, ...)
}

# Fit the common-factor model from already-prepared inputs: a correlation matrix R, the
# sample size N, optional DWLS weights, the optional sandwich meat Gamma, and -- for the
# bootstrap -- pre-resampled correlation/weight arrays. Split out from EFA() so multiple-
# imputation pooling can drive the same estimate -> rotate -> standard-error pipeline on
# pooled inputs (the MI2S route in EFA_POOLED()) without re-entering EFA()'s raw-data
# preparation and input guards. `fiml` (point estimate) and `fiml_boot` (per-replicate)
# carry the two-stage EM moments + raw data for the FIML likelihood-ratio fit indices and
# default to NULL; EFA_POOLED()'s MI2S route never supplies them (it is gated to poly/tetra/
# pearson), so its fit indices stay on the standard discrepancy path.
.efa_core <- function(R, N, weights = NULL, Gamma = NULL,
                      R_boot_array = NULL, W_boot_array = NULL,
                      fiml = NULL, fiml_boot = NULL,
                      np_boot = FALSE, b_boot = 1000, method, rotation, type,
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

  # run factor analysis with respective fit method

  if (method %in% c("ML", "ULS", "DWLS")) {

    if (type == "SPSS") {

      cli::cli_warn(
        c("Only {.val PAF} is validated against the SPSS implementation.",
          "i" = "{.val {method}} results may differ from those returned by SPSS."),
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

  fit_out <- .estimate_model(R, method = method, n_factors = n_factors, N = N,
                             type = type, max_iter = max_iter,
                             init_comm = init_comm, criterion = criterion,
                             criterion_type = criterion_type,
                             abs_eigen = abs_eigen, start_method = start_method,
                             weights = weights, fiml = fiml)

  # Surface Heywood cases from the point-estimate solution (the detector runs in
  # .finalize_fit for every fit). This fires once per EFA() call; EFA_AVERAGE,
  # which fits one EFA per grid cell, suppresses these per-model warnings and
  # reports a single summary instead.
  if (length(fit_out$heywood) > 0) {
    heywood_vars <- names(fit_out$heywood)
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
  # EFA() call; the bootstrap replicates suppress their per-fit warnings and are
  # tallied separately in .boot_fun().
  if (method %in% c("ML", "ULS", "DWLS") && isTRUE(fit_out$convergence != 0)) {
    cli::cli_warn(
      c("The {.val {method}} optimiser did not converge (convergence code {fit_out$convergence}).",
        "i" = paste("It stopped before meeting the convergence tolerance (typically the maximum",
                    "number of iterations was reached); the results may not be interpretable.")),
      class = "efa_nonconvergence"
    )
  }

  # Only the bootstrap produces replicate fits; the analytic SE methods leave this NULL.
  boot_fits <- NULL
  if (isTRUE(np_boot)) {

    boot_fits <- .boot_fun(R_boot_array, b_boot, .estimate_model,
                           # .estimate_model arguments:
                           method = method, n_factors = n_factors, N = N,
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
    if (type == "SPSS" && !rotation %in% c("promax", "varimax")) {

      cli::cli_warn(
        c("Only the {.val promax} and {.val varimax} rotations are validated against the SPSS implementation.",
          "i" = "{.val {rotation}} results may differ from those returned by SPSS."),
        class = "efa_spss_rotation_untested"
      )

    }

    rot_out <- .rotate_model(fit_out, rotation = rotation, type = type,
                             normalize = normalize, precision = precision,
                             order_type = order_type, varimax_type = varimax_type,
                             P_type = P_type, k = k, randomStarts = randomStarts,
                             ...)

    boot_rot <- .rotation_family(rotation)

  }

  if (rotation != "none"){

    if(method %in% c("ULS", "DWLS")){

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
    method = method,
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

  if(method %in% c("ULS", "DWLS") & rotation == "none"){

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
                               gamma = Gamma, method = method, fiml = fiml)
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

  class(output) <- "EFA"

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

# Dispatch standard-error/confidence-interval computation on the requested method. The
# nonparametric bootstrap aggregates resampled fits; the analytic methods derive SEs from
# the fitted model itself. Every branch returns the same list(SE, CI, replicates) schema, so
# the print and summary methods are agnostic to how the SEs were produced. The sandwich branch
# additionally returns a `scaled_test` element (the robust scaled chi-square), which EFA()
# strips off and folds into the fit indices. `gamma` (the polychoric ACOV meat) and `method`
# are used only by the sandwich. `fiml` (the two-stage EM moments + raw data) is non-NULL only
# for cor_method = "fiml", where both analytic settings route to the corrected two-stage sandwich.
.compute_se_ci <- function(fit_out, L_rot, se_method, boot_fits = NULL,
                           boot_rot = "none", ci = .95, b = NULL, N = NULL,
                           rot_info = NULL, gamma = NULL, method = NULL,
                           fiml = NULL) {

  # cor_method = "fiml": the EM correlation is a two-stage estimate, so both analytic settings
  # ("information"/"sandwich") return the corrected two-stage sandwich SE built on the saturated
  # FIML covariance. The naive Stage-2 information SE (treating the EM correlation as complete data)
  # is inconsistent under missingness and is never shipped. The bootstrap path is unchanged.
  if (!is.null(fiml) && se_method %in% c("information", "sandwich")) {
    return(.se_fiml(fit_out, rot_info, N, ci, fiml, method))
  }

  switch(se_method,
    "np-boot" = .boot_se_ci(fit_out, L_rot, boot_fits, boot_rot, ci, b),
    # rot_info is non-NULL exactly for an analytic-SE fit under a real rotation; the unrotated
    # path (rotation = "none") leaves it NULL.
    "information" = if (is.null(rot_info)) {
      .se_information(fit_out, N, ci)
    } else {
      .se_information_rotated(fit_out, rot_info, N, ci)
    },
    "sandwich" = .se_sandwich_dispatch(fit_out, rot_info, N, ci, gamma, method),
    NULL
  )
}


# Leading `keep`-block of the inverse of the bordered matrix [A C'; C 0] -- the constrained
# (reflexive generalised) inverse of a gauge-singular information matrix A under the identification
# constraints whose Jacobian is `Cmat` (k(k-1)/2 rows, zero rows for a single factor). Returns NULL
# if the augmented system is singular. Shared by the expected-information and sandwich SEs.
.bordered_inverse_block <- function(A, Cmat, keep) {
  nc <- nrow(Cmat)
  Aug <- if (nc > 0L) {
    rbind(cbind(A, t(Cmat)), cbind(Cmat, matrix(0, nc, nc)))
  } else {
    A
  }
  inv <- tryCatch(solve(Aug), error = function(e) NULL)
  if (is.null(inv)) return(NULL)
  inv[seq_len(keep), seq_len(keep), drop = FALSE]
}

# Delta-method SEs of the communalities h2_i = rowSums(Lambda^2)_i from a loading covariance V
# (p*k x p*k over column-major vec(Lambda)). The gradient of h2_i is 2 Lambda[i, ], nonzero only in
# variable i's loading columns. Used by the sandwich path, where it supplies the coinciding
# uniqueness and communality SEs from the robust V_AA (the expected-information path instead reports
# the psi-block value, which the communalities share with the uniquenesses; see `.se_information_rotated`).
.communality_se <- function(L, V) {
  p <- nrow(L)
  k <- ncol(L)
  G_h <- matrix(0, p, p * k)
  for (i in seq_len(p)) G_h[i, (seq_len(k) - 1L) * p + i] <- 2 * L[i, ]
  sqrt(pmax(rowSums((G_h %*% V) * G_h), 0))
}

# Skeleton gradient of the rotational identification constraint for the factor pair (u, v): a
# p x k matrix with column u set to X[, v] and column v to X[, u]. The off-diagonal (u, v) entry
# of X' Lambda has exactly this gradient in vec(Lambda) -- with X = Lambda it is the gradient of
# off-diag(Lambda' Lambda), and with X = Psi^-1 Lambda the direct part of off-diag(Lambda' Psi^-1
# Lambda); callers add any free-psi column or psi(Lambda) chain-rule term. Shared by the expected-
# information and sandwich gauge constraints.
.gauge_grad <- function(X, u, v) {
  g <- matrix(0, nrow(X), ncol(X))
  g[, u] <- X[, v]
  g[, v] <- X[, u]
  g
}

# TRUE if M is a (numerically) positive-semidefinite covariance: finite, with a non-negative
# diagonal and no eigenvalue below a small tolerance scaled by its largest diagonal entry. Used to
# gate the analytic parameter covariances before their square roots are reported: checking only the
# diagonal would pass a covariance that is non-PSD off the diagonal (where `pmax(., 0)` would then
# silently floor a negative rotated variance to an understated SE rather than flag it).
.is_psd <- function(M) {
  if (!all(is.finite(M))) return(FALSE)
  d <- diag(M)
  if (any(d < 0)) return(FALSE)
  sym <- (M + t(M)) / 2
  ev <- tryCatch(min(eigen(sym, symmetric = TRUE, only.values = TRUE)$values),
                 error = function(e) NA_real_)
  is.finite(ev) && ev >= -1e-8 * max(d, 1)
}


# Analytic information-matrix SEs and Wald CIs for the unrotated ML loadings and
# uniquenesses. Fills the same SE/CI schema as the bootstrap; the rotated-loading,
# fit-index and residual slots have no closed-form information here and stay absent, and
# there are no replicate arrays.
.se_information <- function(fit_out, N, ci) {

  L <- fit_out$unrot_loadings
  # The ML solution reproduces the unit diagonal of the correlation matrix, so the
  # uniquenesses are 1 minus the communalities of the unrotated loadings.
  psi <- 1 - rowSums(L^2)

  se <- .se_information_ml(L, psi, N)
  SE_L <- se$loadings_se
  SE_psi <- se$uniquenesses_se
  pk <- length(L)
  # The full bordered vcov covers loadings AND uniquenesses (q = pk + p); the persisted slot
  # is the leading pk x pk loading block (matches the sandwich V_AA schema and what
  # `.se_information_rotated()` already extracts inline).
  V_unrot <- se$vcov[seq_len(pk), seq_len(pk), drop = FALSE]
  # `.se_information_ml()` returns an all-NA vcov together with NA marginal SEs whenever the
  # parameter covariance is unusable (a Heywood / singular early return, or a non-PSD covariance
  # caught by the PSD gate). Mirror that NA into the persisted slot defensively, so a finite vcov is
  # never shipped next to NA SEs (which would otherwise propagate sqrt(NaN) downstream).
  if (anyNA(SE_L)) V_unrot[] <- NA_real_

  if (anyNA(SE_L) || anyNA(SE_psi)) {
    cli::cli_warn(
      c("Analytic standard errors could not be computed for all parameters.",
        "i" = "This occurs at a Heywood case (a uniqueness at or below its zero boundary) or when the expected information matrix is singular."),
      class = "efa_se_unreliable"
    )
  }

  dimnames(SE_L) <- dimnames(L)
  names(SE_psi) <- rownames(L)

  # Wald intervals around the point estimates.
  z <- stats::qnorm(1 - (1 - ci) / 2)
  L_lower <- L - z * SE_L
  L_upper <- L + z * SE_L
  psi_lower <- psi - z * SE_psi
  psi_upper <- psi + z * SE_psi
  names(psi) <- names(psi_lower) <- names(psi_upper) <- rownames(L)

  list(
    SE = list(
      unrot_loadings = SE_L,
      uniquenesses = SE_psi
    ),
    CI = list(
      unrot_loadings = list(lower = L_lower, upper = L_upper),
      uniquenesses = list(lower = psi_lower, upper = psi_upper)
    ),
    replicates = NULL,
    vcov_unrot_loadings = V_unrot
  )
}


# Expected-information (Fisher) covariance of the unrotated ML loadings and uniquenesses,
# under the EFA identification constraint that Lambda' Psi^{-1} Lambda is diagonal.
#
# theta = (vec(Lambda) [p*k], psi [p]); the unit expected information
# A[i, j] = 1/2 tr(W dSig_i W dSig_j), with W = Sigma^{-1} and Sigma = Lambda Lambda' +
# diag(psi), is assembled in closed form from W, B = Sigma^{-1} Lambda and
# G = Lambda' Sigma^{-1} Lambda (Jennrich & Thayer, 1973; Lawley & Maxwell). The constraint
# removes the k(k-1)/2 rotational degrees of freedom and is imposed with a bordered
# information matrix; the parameter covariance is the leading block of its inverse, scaled
# by 1 / (N - 1) (matching the N - 1 convention of the fit statistics). R is accepted for
# forward compatibility (sandwich SEs) but is unused by the expected-information path.
#
# Validity: this is the normal-theory expected information, consistent only under
# multivariate normality and a correctly specified model (Yuan & Hayashi, 2006). In that
# regime it agrees with resampling SEs; under heavy-tailed data or model misfit it tends to
# UNDERESTIMATE the true sampling variability (analytic < bootstrap), so a robust sandwich
# or bootstrap SE is preferable there (Zhang, Preacher, & Jennrich, 2012; Zhang et al.,
# 2019). Cross-loading (nonmarker) variables also carry larger SEs than single-factor
# markers (Zhang & Preacher, 2015).
.se_information_ml <- function(L, psi, N, R = NULL) {

  p <- nrow(L)
  k <- ncol(L)
  pk <- p * k
  q <- pk + p

  na_out <- list(
    vcov = matrix(NA_real_, q, q),
    loadings_se = matrix(NA_real_, p, k),
    uniquenesses_se = rep(NA_real_, p)
  )

  # A Heywood case (a uniqueness at or below its zero boundary) makes Psi^{-1} and the
  # Lambda' Psi^{-1} Lambda identification ill-defined, so no analytic SE is available.
  if (anyNA(psi) || any(psi <= 0)) {
    return(na_out)
  }

  Sigma <- tcrossprod(L)
  diag(Sigma) <- diag(Sigma) + psi
  W <- tryCatch(solve(Sigma), error = function(e) NULL)
  if (is.null(W)) {
    return(na_out)
  }

  B <- W %*% L              # Sigma^{-1} Lambda          (p x k)
  G <- crossprod(L, B)      # Lambda' Sigma^{-1} Lambda  (k x k)

  lidx <- function(a, b) (b - 1L) * p + a   # column-major index of Lambda[a, b]

  # --- unit expected information A (q x q) -----------------------------------------------
  A <- matrix(0, q, q)

  # loading x loading: A[(a,b),(c,d)] = G[b,d] W[a,c] + B[a,d] B[c,b].
  # The first term is kronecker(G, W) under the column-major loading order.
  A[seq_len(pk), seq_len(pk)] <- kronecker(G, W)
  for (b in seq_len(k)) {
    rows <- lidx(seq_len(p), b)
    for (d in seq_len(k)) {
      cols <- lidx(seq_len(p), d)
      A[rows, cols] <- A[rows, cols] + outer(B[, d], B[, b])
    }
  }

  # loading x uniqueness: A[(a,b), psi_c] = W[a,c] B[c,b].
  for (b in seq_len(k)) {
    rows <- lidx(seq_len(p), b)
    A[rows, pk + seq_len(p)] <- W * matrix(B[, b], p, p, byrow = TRUE)
  }
  A[pk + seq_len(p), seq_len(pk)] <- t(A[seq_len(pk), pk + seq_len(p)])

  # uniqueness x uniqueness: A[psi_a, psi_c] = 1/2 W[a,c]^2.
  A[pk + seq_len(p), pk + seq_len(p)] <- 0.5 * W^2

  # --- identification constraint: off-diagonals of Lambda' Psi^{-1} Lambda = 0 -----------
  nc <- k * (k - 1L) / 2L
  Astar <- L / psi          # Psi^{-1} Lambda (row a divided by psi[a])
  Cmat <- matrix(0, nc, q)
  uv <- 0L
  if (k > 1L) {
    for (v in 2:k) {
      for (u in seq_len(v - 1L)) {
        uv <- uv + 1L
        Cmat[uv, seq_len(pk)] <- as.vector(.gauge_grad(Astar, u, v))
        Cmat[uv, pk + seq_len(p)] <- -Astar[, u] * Astar[, v]
      }
    }
  }

  # --- bordered information, inverse, scaling --------------------------------------------
  Vblock <- .bordered_inverse_block(A, Cmat, q)
  if (is.null(Vblock)) {
    return(na_out)
  }

  V <- Vblock / (N - 1)
  # The parameter covariance must be positive semidefinite. A near-degenerate orientation can
  # leave it non-PSD -- a tiny negative variance on the diagonal, or a negative eigenvalue with a
  # still-positive diagonal -- and no trustworthy standard error then exists, so fall back to the
  # all-NA return rather than report the square root of a negative variance.
  if (!.is_psd(V)) {
    return(na_out)
  }
  se <- sqrt(diag(V))

  list(
    vcov = V,
    loadings_se = matrix(se[seq_len(pk)], p, k),
    uniquenesses_se = se[pk + seq_len(p)]
  )
}


# Map a rotation name to the warm-start criterion family and tuning argument used by the compiled
# `.rotation_se_jacobian()`, mirroring the engine selection in `.orth_engines`/`.oblq_engines`
# (R/rotate_model.R). varimax shares the Crawford-Ferguson criterion at kappa = 1 / p, quartimax at
# kappa = 0, equamax at kappa = k / (2 p); quartimin is oblimin at gam = 0. The geomin offset
# `delta` and oblimin `gam` are taken from the resolved criterion arguments (the same `.gpf_crit`
# defaults the rotation itself used). promax and simplimax have no usable analytic Jacobian and are
# rejected before reaching here, so an unrecognised rotation returns NULL.
.rotation_se_method <- function(rotation, p, k, crit_args = list()) {
  gam <- if (is.null(crit_args$gam)) 0 else crit_args$gam
  delta <- if (is.null(crit_args$delta)) 0.01 else crit_args$delta
  switch(rotation,
    varimax   = list(method = "cf",       param = 1 / p,       oblique = FALSE),
    quartimax = list(method = "cf",       param = 0,           oblique = FALSE),
    equamax   = list(method = "cf",       param = k / (2 * p), oblique = FALSE),
    bentlerT  = list(method = "bentler",  param = 0,           oblique = FALSE),
    geominT   = list(method = "geomin",   param = delta,       oblique = FALSE),
    bifactorT = list(method = "bifactor", param = 0,           oblique = FALSE),
    oblimin   = list(method = "oblimin",  param = gam,         oblique = TRUE),
    quartimin = list(method = "oblimin",  param = 0,           oblique = TRUE),
    bentlerQ  = list(method = "bentler",  param = 0,           oblique = TRUE),
    geominQ   = list(method = "geomin",   param = delta,       oblique = TRUE),
    bifactorQ = list(method = "bifactor", param = 0,           oblique = TRUE),
    NULL)
}

# Wald lower/upper interval around an estimate from its standard error. `est` is unclassed so the
# returned bounds are plain numeric matrices/vectors (the rotated estimates carry the "LOADINGS"
# class).
.wald_ci <- function(est, se, z) {
  est <- unclass(est)
  list(lower = est - z * se, upper = est + z * se)
}


# Analytic standard errors and Wald CIs for an obliquely or orthogonally ROTATED ML solution.
# The expected-information covariance of the unrotated loadings (P6.1, `.se_information_ml`) is
# propagated through the rotation by the delta method. The rotation maps the unrotated loadings A
# to the rotated pattern L = A T(A)^{-1 T} (oblique) / A T(A) (orthogonal) and, for oblique
# rotations, the factor correlations Phi = T(A)' T(A); both depend on A alone (not the
# uniquenesses), so only the p*k loading block V of the parameter covariance is needed. The
# rotation Jacobians d vec(L) / d vec(A) and d vec(Phi) / d vec(A) are obtained by
# finite-differencing a warm-started re-rotation of A (`.rotation_se_jacobian()`), which re-solves
# the rotation from the converged transformation so the optimum tracks each perturbation smoothly.
# This is the delta-method standard error of Jennrich (1973), matching lavaan's rotation.se =
# "delta".
# Because rotated quantities are identification-invariant they are comparable across packages
# (unlike the unrotated loadings, whose orientation is identification-dependent).
#
#   Var(vec L)   = J_L V J_L'                                   -> rotated-loading SEs
#   Var(vec Phi) = J_Phi V J_Phi'                              -> factor-correlation SEs (oblique)
#   S = L Phi:  J_S = (Phi' (x) I_p) J_L + (I_k (x) L) J_Phi;  Var(vec S) = J_S V J_S'  (oblique)
#   h2_i = rowSums(A^2)_i = 1 - psi_i (rotation-invariant); SE(h2_i) = SE(psi_i), the uniqueness SE
#
# Fills the same SE/CI schema as the bootstrap (rot_loadings/Phi/Structure plus the unrotated
# loadings and uniquenesses); there are no replicate arrays. Falls back to NA rotated SEs (with a
# classed warning) at a Heywood case, a singular information matrix, or a rotation the warm start
# cannot reproduce (e.g. a non-converged transformation). References: Jennrich (1973); Zhang &
# Preacher (2015).
.se_information_rotated <- function(fit_out, rot_info, N, ci, se0 = NULL) {

  A <- unclass(fit_out$unrot_loadings)
  p <- nrow(A)
  k <- ncol(A)
  pk <- p * k
  psi <- 1 - rowSums(A^2)
  z <- stats::qnorm(1 - (1 - ci) / 2)

  # Unrotated pieces (always reported, mirroring the bootstrap's unrot_loadings) and the
  # loading covariance block propagated through the rotation. The covariance is the expected-
  # information one by default; the sandwich path supplies its own robust covariance in `se0`
  # (with the same $vcov/$loadings_se/$uniquenesses_se schema), so the rotation propagation
  # below is identical for both and is purely a function of the loading covariance V.
  if (is.null(se0)) se0 <- .se_information_ml(A, psi, N)
  V <- se0$vcov[seq_len(pk), seq_len(pk), drop = FALSE]
  SE_unrot <- se0$loadings_se
  SE_psi <- se0$uniquenesses_se
  dimnames(SE_unrot) <- dimnames(A)
  names(SE_psi) <- rownames(A)
  # `.se_information_ml()` NA's the unrotated SEs whenever the parameter covariance is non-finite
  # OR carries a negative variance on its diagonal (a near-degenerate orientation `solve()` still
  # inverts). Reuse that as the single reliability signal for the rotated quantities too: they are
  # propagated from the same covariance, so if it was too ill-conditioned for the unrotated SEs it
  # must not yield finite rotated SEs either.
  info_reliable <- !anyNA(SE_unrot)

  spec <- .rotation_se_method(rot_info$rotation, p, k, rot_info$crit_args)
  rotmat <- rot_info$rotmat
  L_rot <- unclass(rot_info$rot_loadings)
  Phi_pt <- rot_info$Phi
  # Orthogonal rotations and the single-factor fall-back carry no factor correlations, so they
  # report rotated loadings and communalities only.
  oblique <- !is.null(Phi_pt)

  # The communalities are the exact complement of the uniquenesses (h2_i = rowSums(A^2)_i =
  # 1 - psi_i = (L Phi L')_ii), so h2_i and psi_i are one estimand up to sign and must share a
  # standard error. SE_psi already carries it -- from the psi-block of the bordered parameter
  # covariance on the expected-information path, and from the loading gradient on the robust V_AA
  # on the sandwich path (where the uniqueness and communality gradients coincide). Report it for
  # both. (A loading-gradient delta 2A[i, ]' V 2A[i, ] on the expected-information covariance is a
  # second, less accurate route to the same variance: it overstates Var(h2_i) against the sampling
  # distribution, whereas the psi-block matches it.)
  h2 <- rowSums(A^2)
  names(h2) <- rownames(A)
  SE_h2 <- SE_psi
  names(SE_h2) <- rownames(A)

  na_mat <- function() {
    m <- matrix(NA_real_, p, k)
    dimnames(m) <- dimnames(L_rot)
    m
  }
  SE_rot <- na_mat()
  SE_Phi <- if (oblique) matrix(NA_real_, k, k) else NULL
  SE_S <- if (oblique) na_mat() else NULL
  S_pt <- if (oblique) `dimnames<-`(L_rot %*% Phi_pt, dimnames(L_rot)) else NULL

  can_rotate <- !is.null(spec) && is.matrix(rotmat) && !anyNA(rotmat) &&
    k >= 2L && !is.null(L_rot) && info_reliable

  if (can_rotate) {
    normalize <- isTRUE(rot_info$normalize)
    # The bifactor criterion exempts a fixed column as the general factor, but `.reflect_and_order`
    # may have moved the general factor out of the first column. It is the one column that loads on
    # every variable -- so the column whose smallest absolute loading is largest, the others being
    # group factors with near-zero out-group loadings. Identify it so the re-rotation optimises the
    # same criterion the point estimate did; the other criteria ignore general_col.
    general_col <- if (identical(spec$method, "bifactor")) {
      which.max(apply(abs(L_rot), 2, min)) - 1L
    } else {
      0L
    }
    # Forward difference of the warm-started re-rotation (the stencil lavaan's delta method uses);
    # the whole p*k finite-difference loop runs in compiled code (.rotation_se_jacobian()) to avoid
    # p*k round trips to R. This is a one-shot cost per fit: a few hundredths of a second for a
    # typical problem (p = 18, k = 3) and about 1 s for a large one (p = 40, k = 8) -- far below a
    # bootstrap, and an order of magnitude faster than the equivalent pure-R numeric differentiation.
    eps <- 1e-4
    jac <- .rotation_se_jacobian(A, rotmat, spec$method, spec$param, normalize,
                                 spec$oblique, eps, general_col)
    # The re-rotation at the unperturbed A must reproduce the reported rotated loadings; a gross
    # mismatch flags a criterion/parameter mismatch or a transformation outside the criterion's
    # basin, for which the Jacobian is not trustworthy. The tolerance absorbs the small drift
    # between the reported optimum and a re-solve to a tighter tolerance -- and, for varimax, between
    # the SVD/SPSS varimax algorithm and the Crawford-Ferguson re-rotation at a flat optimum -- while
    # still catching a genuine basin/criterion mismatch (which moves loadings an order of magnitude
    # more).
    if (isTRUE(jac$valid) && max(abs(jac$base_loadings - L_rot)) < 5e-2) {
      J_L <- jac$J_L
      SE_rot <- matrix(sqrt(pmax(rowSums((J_L %*% V) * J_L), 0)), p, k)
      dimnames(SE_rot) <- dimnames(L_rot)
      if (spec$oblique) {
        J_Phi <- jac$J_Phi
        SE_Phi <- matrix(sqrt(pmax(rowSums((J_Phi %*% V) * J_Phi), 0)), k, k)
        diag(SE_Phi) <- 0       # the unit diagonal of Phi is fixed, so it has no variance
        SE_Phi <- (SE_Phi + t(SE_Phi)) / 2
        # Structure S = L Phi: dvec(S) = (Phi' (x) I_p) dvec(L) + (I_k (x) L) dvec(Phi); the
        # rotation matrices are the re-rotated point estimates, consistent with the Jacobians.
        L0 <- jac$base_loadings
        Phi0 <- jac$base_Phi
        J_S <- kronecker(t(Phi0), diag(p)) %*% J_L + kronecker(diag(k), L0) %*% J_Phi
        SE_S <- matrix(sqrt(pmax(rowSums((J_S %*% V) * J_S), 0)), p, k)
        dimnames(SE_S) <- dimnames(L_rot)
      }
    }
  } else if (info_reliable && k < 2L) {
    # No rotation was actually performed (single factor): the rotated loadings equal the unrotated
    # ones, so their SEs are the unrotated loading SEs; there is no Phi or structure matrix.
    SE_rot <- SE_unrot
    dimnames(SE_rot) <- dimnames(L_rot)
  }

  SE <- list(unrot_loadings = SE_unrot, uniquenesses = SE_psi,
             rot_loadings = SE_rot, communalities = SE_h2)
  CI <- list(unrot_loadings = .wald_ci(A, SE_unrot, z),
             uniquenesses = .wald_ci(psi, SE_psi, z),
             rot_loadings = .wald_ci(L_rot, SE_rot, z),
             communalities = .wald_ci(h2, SE_h2, z))
  if (oblique) {
    SE$Phi <- SE_Phi
    SE$Structure <- SE_S
    CI$Phi <- .wald_ci(Phi_pt, SE_Phi, z)
    CI$Structure <- .wald_ci(S_pt, SE_S, z)
  }

  if (anyNA(SE_unrot) || anyNA(SE_psi) || anyNA(SE_rot) || anyNA(SE_h2) ||
      (oblique && (anyNA(SE_Phi) || anyNA(SE_S)))) {
    cli::cli_warn(
      c("Analytic standard errors could not be computed for all parameters.",
        "i" = "This occurs at a Heywood case (a uniqueness at or below its zero boundary), when the parameter covariance is singular, or when the rotation could not be reproduced for the standard-error Jacobian."),
      class = "efa_se_unreliable"
    )
  }

  # Surface the pk x pk unrotated loading vcov (the information block, or the sandwich V_AA
  # when `se0` was supplied) so it can be persisted on the EFA object alongside the rotated SEs.
  # NA-fill if the underlying unrotated SEs were unreliable (info_reliable signals both the
  # Heywood early-return and the finite-but-not-PSD covariance path) so the slot is consistent
  # with the marginal SEs and downstream consumers can fail closed on `anyNA()`.
  if (!info_reliable) V[] <- NA_real_
  list(SE = SE, CI = CI, replicates = NULL, vcov_unrot_loadings = V)
}


# Robust (sandwich) standard errors and a scaled chi-square for the ordinal/polychoric path.
# The factor model is fitted by minimising a weighted off-diagonal discrepancy
# F(theta) = (s - sigma(theta))' V (s - sigma(theta)) over the free loadings theta = vec(Lambda),
# where s = vech of the off-diagonal sample (polychoric/tetrachoric) correlations,
# sigma_ij(theta) = (Lambda Lambda')_ij, and V is the estimator's weight matrix (DWLS: diagonal
# inverse asymptotic variances; ULS: identity; ML: the normal-theory GLS weight). When V is not
# the inverse of the true asymptotic covariance of s (which it never is for ordinal data), the
# naive inverse-information SEs are biased, and the Godambe sandwich
#   Var(theta) = (1/(N-1)) A^- (Delta' V Gamma V Delta) A^-,   A = Delta' V Delta,
# with Gamma the (threshold-adjusted) asymptotic covariance of s, gives consistent SEs (Browne,
# 1984; Muthen, 1984; Satorra & Bentler, 1994). This mirrors lavaan's se = "robust.sem". The
# rotational over-parameterisation of A is handled by bordering it with the identification-
# constraint Jacobian, exactly as the expected-information path does. The bordered loading
# covariance is propagated through a rotation by the same Jacobian machinery as the information
# SEs (.se_information_rotated), so promax/simplimax stay bootstrap-only.
.se_sandwich_dispatch <- function(fit_out, rot_info, N, ci, Gamma, method) {

  core <- .se_sandwich_core(fit_out, N, Gamma, method)

  res <- if (is.null(rot_info)) {
    .se_sandwich_unrotated(fit_out, core, ci)
  } else {
    # Reuse the information-SE rotation propagation verbatim by supplying the robust loading
    # covariance in place of the expected-information one (identical $vcov schema). Because the
    # rotation Jacobian annihilates the rotational gauge directions, the rotated SEs are
    # identification-invariant and do not depend on which bordering constraint produced V_AA.
    se0 <- list(vcov = core$V_AA,
                loadings_se = core$loadings_se,
                uniquenesses_se = core$uniquenesses_se)
    .se_information_rotated(fit_out, rot_info, N, ci, se0 = se0)
  }

  res$scaled_test <- core$scaled_test
  res
}

# Identification-constraint Jacobian for the rotational gauge: the gradient, with respect to the
# column-major vec(Lambda), of the off-diagonal entries the estimator's canonical orientation sets
# to zero. nc = k(k-1)/2 rows, p*k columns; for k = 1 there is no gauge freedom and the matrix has
# zero rows. The constraint must match the orientation the loadings are reported in so the unrotated
# loading SEs are scaled in that gauge (the rotated SEs are gauge-invariant and the same for any
# transversal choice). With `psi = NULL` (eigen-based ULS/DWLS) it fixes off-diag(Lambda'Lambda);
# with `psi` supplied (ML, where psi = 1 - rowSums(Lambda^2)) it fixes off-diag(Lambda' Psi^-1
# Lambda) (Lawley & Maxwell), whose gradient carries a chain-rule term through psi(Lambda).
.se_sandwich_constraint <- function(L, psi = NULL) {
  p <- nrow(L)
  k <- ncol(L)
  nc <- k * (k - 1L) / 2L
  Cmat <- matrix(0, nc, p * k)
  if (k > 1L) {
    Astar <- if (is.null(psi)) NULL else L / psi   # Psi^-1 Lambda for the ML gauge
    uv <- 0L
    for (v in 2:k) {
      for (u in seq_len(v - 1L)) {
        uv <- uv + 1L
        grad <- if (is.null(psi)) {
          .gauge_grad(L, u, v)
        } else {
          # d off-diag(Lambda' Psi^-1 Lambda)_uv / d vec(Lambda): the direct part plus the term
          # from psi_a = 1 - rowSums(Lambda^2)_a depending on Lambda (d psi_a / d lambda_ac = -2
          # lambda_ac), which contributes 2 (Astar[, u] Astar[, v]) Lambda across all columns.
          .gauge_grad(Astar, u, v) + 2 * (Astar[, u] * Astar[, v]) * L
        }
        Cmat[uv, ] <- as.vector(grad)
      }
    }
  }
  Cmat
}

# Build the robust loading covariance V_AA (p*k x p*k), the unrotated loading/uniqueness SEs, and
# the scaled chi-square from the fitted loadings and the polychoric ACOV `Gamma`. `Gamma` enters on
# the variance scale Var(rho-hat) and is converted here to the unit asymptotic-variance scale
# (N * Var) used by the WLS/sandwich formulas (= lavaan's NACOV). Returns NA SEs (reliable = FALSE)
# at a singular bordered information matrix or an unusable covariance.
#
# The linear algebra is O(n^2 q) for the meat and O(n^3) for the chi-square trace (n = p(p-1)/2
# pairs, q = p*k), and is kept in R: it is a one-shot cost of a few hundredths of a second for a
# typical problem and about 0.3 s (p = 30) to 1.5 s (p = 40, n = 780 pairs) for a large one --
# negligible next to the polychoric estimate it builds on, and far below a bootstrap.
.se_sandwich_core <- function(fit_out, N, Gamma, method) {

  L <- unclass(fit_out$unrot_loadings)
  p <- nrow(L)
  k <- ncol(L)
  pk <- p * k

  na_core <- list(
    V_AA = matrix(NA_real_, pk, pk),
    loadings_se = matrix(NA_real_, p, k),
    uniquenesses_se = rep(NA_real_, p),
    scaled_test = NULL,
    reliable = FALSE
  )

  if (is.null(Gamma) || anyNA(Gamma) || nrow(Gamma) != p * (p - 1L) / 2L) {
    return(na_core)
  }

  # Variance scale (Var(rho-hat)) -> unit asymptotic-variance scale (N * Var); the 1/(N-1) on the
  # final covariance and the N on the chi-square then follow lavaan's robust.sem conventions.
  Gamma <- N * Gamma

  pairs <- utils::combn(p, 2L)
  pi <- pairs[1, ]
  pj <- pairs[2, ]
  n <- ncol(pairs)

  # Model Jacobian Delta = d sigma_offdiag / d vec(Lambda) (n x pk). For pair (i, j),
  # d sigma_ij / d Lambda[a, f] = (a == i) Lambda[j, f] + (a == j) Lambda[i, f]; rows in
  # utils::combn(p, 2) order, columns in column-major vec(Lambda) order.
  Delta <- matrix(0, n, pk)
  for (f in seq_len(k)) {
    Delta[cbind(seq_len(n), (f - 1L) * p + pi)] <- L[pj, f]
    Delta[cbind(seq_len(n), (f - 1L) * p + pj)] <- L[pi, f]
  }

  # Estimator weight V (unit scale). DWLS: diagonal inverse variances; ULS: identity; ML: the
  # normal-theory GLS weight 1/2 (Sigma^-1 (x) Sigma^-1) restricted to the off-diagonal pairs, at
  # the model-implied correlation matrix Sigma = Lambda Lambda' (unit diagonal).
  if (method == "DWLS") {
    vdiag <- 1 / diag(Gamma)
    Vmat <- NULL
  } else if (method == "ULS") {
    vdiag <- rep(1, n)
    Vmat <- NULL
  } else {
    Sigma <- tcrossprod(L)
    diag(Sigma) <- 1
    P <- tryCatch(solve(Sigma), error = function(e) NULL)
    if (is.null(P)) return(na_core)
    vdiag <- NULL
    Vmat <- matrix(0, n, n)
    for (s in seq_len(n)) {
      a <- pi[s]
      b <- pj[s]
      Vmat[, s] <- 0.5 * (P[pi, a] * P[pj, b] + P[pi, b] * P[pj, a])
    }
  }
  is_diag <- is.null(Vmat)

  # VD = V Delta; bread A = Delta' V Delta (singular by the k(k-1)/2 rotational gauge freedoms).
  VD <- if (is_diag) vdiag * Delta else Vmat %*% Delta
  A <- crossprod(Delta, VD)

  # Border A with the gauge-fixing constraint so the augmented system is invertible; the leading
  # pk-block of its inverse is the reflexive generalised inverse of A. The constraint fixes the
  # rotational orientation, so it must match the one the solution's loadings are reported in: the
  # unrotated loading SEs are scaled in that gauge (the rotated SEs and scaled chi-square are
  # gauge-invariant, so the choice does not affect them). Rather than assume the orientation from the
  # estimator label, detect it from the loadings: an eigen-based ULS/DWLS solution is
  # Lambda'Lambda-diagonal, an ML solution is Lambda'Psi^-1 Lambda-diagonal (Lawley & Maxwell). The
  # orientation the solution is in leaves a (scale-free) relative off-diagonal at the optimiser floor
  # while the other is O(0.01-1), so the smaller ratio identifies the gauge -- and detecting keeps it
  # tied to the actual solution, so a future estimator with either identification is handled without
  # special-casing. Two cases the loadings cannot resolve fall back to a fixed choice: the
  # Lambda'Psi^-1 Lambda gauge needs Psi^-1, so it is unavailable at a Heywood case (psi <= 0) or a
  # non-finite solution -- routed to Lambda'Lambda; and homogeneous uniquenesses (Psi proportional to
  # I) make BOTH orientations diagonal, so the gauge is undetermined by the loadings (yet the two
  # still give different SEs through the chain-rule term), broken by the estimator's identification.
  # A single factor has no rotational freedom (the constraint is empty either way).
  psi <- 1 - rowSums(L^2)
  use_ltpil_gauge <- if (k < 2L || anyNA(psi) || any(psi <= 0)) {
    FALSE
  } else {
    rel_off <- function(M) max(abs(M[upper.tri(M)])) / max(abs(diag(M)), .Machine$double.eps)
    r_ltl <- rel_off(crossprod(L))
    r_ltpil <- rel_off(crossprod(L, L / psi))
    # Both orientations near-diagonal (Psi proportional to I): the loadings cannot distinguish the
    # gauge, so use the estimator's known identification; otherwise the smaller ratio wins.
    if (max(r_ltl, r_ltpil) < 1e-4) method == "ML" else r_ltpil < r_ltl
  }
  Cmat <- if (use_ltpil_gauge) {
    .se_sandwich_constraint(L, psi = psi)
  } else {
    .se_sandwich_constraint(L)
  }
  Abread <- .bordered_inverse_block(A, Cmat, pk)
  if (is.null(Abread)) return(na_core)

  # Meat Delta' V Gamma V Delta = (V Delta)' Gamma (V Delta); robust covariance A^- meat A^- /(N-1).
  Gamma_theta <- crossprod(VD, Gamma %*% VD)
  V_AA <- (Abread %*% Gamma_theta %*% Abread) / (N - 1)

  # The covariance must be positive semidefinite: a non-finite entry or a negative eigenvalue (even
  # with a still-positive diagonal) would otherwise corrupt the marginal, uniqueness, and rotated
  # SEs, which read the full V_AA.
  reliable <- .is_psd(V_AA)
  loadings_se <- if (reliable) sqrt(diag(V_AA)) else rep(NA_real_, pk)

  # Uniqueness SE = communality SE (psi_i = 1 - rowSums(Lambda^2)_i), via the shared gradient.
  uniq_se <- if (!reliable) rep(NA_real_, p) else .communality_se(L, V_AA)

  scaled_test <- if (!reliable) NULL else {
    .scaled_chisq(fit_out, Gamma, pairs, VD, vdiag, Vmat, Abread, N)
  }

  list(
    V_AA = V_AA,
    loadings_se = matrix(loadings_se, p, k),
    uniquenesses_se = uniq_se,
    scaled_test = scaled_test,
    reliable = reliable
  )
}

# Scaled chi-square test statistics (Satorra & Bentler, 1994; Asparouhov & Muthen, 2010) for the
# weighted off-diagonal fit. T = N * (s - sigma)' V (s - sigma) is the (unscaled) fit statistic;
# U = V - V Delta A^- Delta' V is the residual projector; c1 = tr(U Gamma), c2 = tr((U Gamma)^2)
# its trace coefficients (Gamma on the unit scale). The projector spans only the off-diagonal
# correlation residuals, so this is the two-stage correlation-structure correction (Browne, 1984)
# and is not identical to the full WLSMV statistic of lavaan/Mplus, which also projects the
# thresholds. Returns the scaled-shifted, mean-adjusted, and mean-and-variance-adjusted statistics
# plus the scaled baseline statistic for the robust CFI/TLI/RMSEA. NULL when the model is
# just-identified (df <= 0) or the traces degenerate.
.scaled_chisq <- function(fit_out, Gamma, pairs, VD, vdiag, Vmat, Abread, N) {

  L <- unclass(fit_out$unrot_loadings)
  df <- fit_out$fit_indices$df
  if (is.null(df) || is.na(df) || df <= 0) return(NULL)

  n <- ncol(pairs)
  is_diag <- is.null(Vmat)
  apply_V <- function(x) if (is_diag) vdiag * x else as.vector(Vmat %*% x)

  # The sample correlations s and residuals e are read from the (possibly smoothed) analysis matrix,
  # while Gamma is the asymptotic covariance of the observed un-projected polychoric correlations --
  # the same convention the DWLS point estimate and lavaan use (the weights/ACOV describe the
  # observed correlations even when the matrix is projected to positive definiteness).
  R <- unclass(fit_out$orig_R)
  s <- R[t(pairs)]                       # off-diagonal sample correlations
  e <- (R - tcrossprod(L))[t(pairs)]     # residuals r_ij - (Lambda Lambda')_ij

  Tstat <- N * sum(e * apply_V(e))

  # Residual projector U = V - V Delta A^- Delta' V. The dense V Delta A^- Delta' V term makes U
  # dense anyway, so add the (diagonal or full) weight onto it in place instead of materialising a
  # dense V.
  U <- -(VD %*% Abread %*% t(VD))
  if (is_diag) diag(U) <- diag(U) + vdiag else U <- U + Vmat
  c1 <- sum(U * Gamma)
  UG <- U %*% Gamma
  c2 <- sum(UG * t(UG))
  if (!is.finite(c1) || !is.finite(c2) || c1 <= 0 || c2 <= 0) return(NULL)

  # Independence baseline (all off-diagonal correlations fixed to 0): residual is s, projector the
  # baseline weight V0. V0 is diagonal for every estimator -- it equals the (model-independent)
  # weight for DWLS/ULS, and 1/2 I for ML (the normal-theory weight at Sigma = I) -- so the baseline
  # traces use only its diagonal.
  v0 <- if (is_diag) vdiag else rep(0.5, n)
  Tbase <- N * sum(s * (v0 * s))
  df_base <- n
  c1_base <- sum(v0 * diag(Gamma))       # tr(diag(v0) Gamma)
  V0G <- v0 * Gamma                      # diag(v0) %*% Gamma (row scaling)
  c2_base <- sum(V0G * t(V0G))           # tr((diag(v0) Gamma)^2)

  scaled <- .scaled_variants(Tstat, df, c1, c2)
  chi_null <- if (is.finite(c1_base) && is.finite(c2_base) && c2_base > 0) {
    .scaled_variants(Tbase, df_base, c1_base, c2_base)$T_ss
  } else {
    NA_real_
  }

  list(
    m = nrow(L),
    df = df,
    chi = scaled$T_ss,                   # scaled-shifted (WLSMV default) -> the reported chi-square
    chi_scaling = scaled$a,
    chi_shift = scaled$b,
    chi_unscaled = Tstat,
    chi_mean_adjusted = scaled$T_mean,
    chi_mean_var = scaled$T_mv,
    df_mean_var = scaled$df_mv,
    chi_null = chi_null,
    df_null = df_base
  )
}

# The three scalings of a fit statistic T (Satorra & Bentler, 1994; Asparouhov & Muthen, 2010)
# from its trace coefficients c1 = tr(U Gamma), c2 = tr((U Gamma)^2) and degrees of freedom df:
# mean-adjusted (T * df/c1, reference df); scaled-shifted a*T + b with a = sqrt(df/c2),
# b = df - a*c1 (reference df, the WLSMV default); mean-and-variance-adjusted T * c1/c2 with the
# Satterthwaite df* = c1^2/c2.
.scaled_variants <- function(Tstat, df, c1, c2) {
  a <- sqrt(df / c2)
  list(
    T_mean = Tstat * df / c1,
    a = a,
    b = df - a * c1,
    T_ss = a * Tstat + (df - a * c1),
    df_mv = c1^2 / c2,
    T_mv = Tstat * c1 / c2
  )
}

# Unrotated robust SE/CI wrapper: fills the uniform list(SE, CI, replicates) schema from the
# core's loading/uniqueness SEs, mirroring .se_information() but with the sandwich covariance.
.se_sandwich_unrotated <- function(fit_out, core, ci) {

  L <- unclass(fit_out$unrot_loadings)
  SE_L <- core$loadings_se
  SE_psi <- core$uniquenesses_se
  dimnames(SE_L) <- dimnames(L)
  names(SE_psi) <- rownames(L)

  if (anyNA(SE_L) || anyNA(SE_psi)) {
    cli::cli_warn(
      c("Robust standard errors could not be computed for all parameters.",
        "i" = "This occurs when the bordered information matrix is singular or the asymptotic covariance is not usable."),
      class = "efa_se_unreliable"
    )
  }

  psi <- 1 - rowSums(L^2)
  names(psi) <- rownames(L)
  z <- stats::qnorm(1 - (1 - ci) / 2)

  # `.se_sandwich_core()` always returns a numeric V_AA, even when `reliable = FALSE` (it is
  # only the marginal `loadings_se` that gets NA-filled in the unreliable branch). NA-fill the
  # persisted covariance to match, so a finite-but-not-PSD V_AA does not silently ship next to
  # NA SEs and propagate sqrt(NaN) into downstream pooling.
  V_AA <- core$V_AA
  if (!isTRUE(core$reliable)) V_AA[] <- NA_real_

  list(
    SE = list(unrot_loadings = SE_L, uniquenesses = SE_psi),
    CI = list(unrot_loadings = .wald_ci(L, SE_L, z),
              uniquenesses = .wald_ci(psi, SE_psi, z)),
    replicates = NULL,
    vcov_unrot_loadings = V_AA
  )
}

# Fold the robust scaled chi-square into the fit indices, overwriting the chi-square-derived block
# (.gof() leaves it undefined for DWLS and reports the non-robust discrepancy for ML/ULS). The
# reported chi-square is the scaled-shifted (WLSMV-default) statistic; the mean-adjusted and
# mean-and-variance-adjusted statistics and the scaling/shift are added as extra fields. CFI/TLI/
# RMSEA come from the scaled model and baseline statistics via the shared .chi_fit_indices().
.apply_scaled_test <- function(fit_indices, st, N) {

  if (is.null(st)) return(fit_indices)

  # The scaled (sandwich) model and baseline statistics are already on a comparable scale, so they
  # serve directly as the CFI/TLI/RMSEA noncentrality inputs (chi_cfi / chi_null_cfi).
  idx <- .chi_fit_indices(st$chi, st$df, st$chi_null, st$df_null, N, st$m, ci = TRUE,
                          chi_cfi = st$chi, chi_null_cfi = st$chi_null)

  fit_indices$chi <- st$chi
  fit_indices$df <- st$df
  fit_indices$p_chi <- idx$p_chi
  fit_indices$CFI <- idx$CFI
  fit_indices$TLI <- idx$TLI
  fit_indices$RMSEA <- idx$RMSEA
  fit_indices$RMSEA_LB <- idx$RMSEA_LB
  fit_indices$RMSEA_UB <- idx$RMSEA_UB
  # AIC/BIC/ECVI are likelihood-ratio chi-square information criteria; they have no standard
  # interpretation when built on the moment-scaled (Satorra-Bentler) statistic, so leave them NA
  # rather than report a misleading model-comparison number.
  fit_indices$AIC <- NA_real_
  fit_indices$BIC <- NA_real_
  fit_indices$ECVI <- NA_real_
  fit_indices$chi_null <- st$chi_null
  fit_indices$df_null <- st$df_null
  fit_indices$p_null <- idx$p_null

  fit_indices$chi_scaled_type <- "scaled.shifted"
  fit_indices$chi_scaling <- st$chi_scaling
  fit_indices$chi_shift <- st$chi_shift
  fit_indices$chi_unscaled <- st$chi_unscaled
  fit_indices$chi_mean_adjusted <- st$chi_mean_adjusted
  fit_indices$chi_mean_var <- st$chi_mean_var
  fit_indices$df_mean_var <- st$df_mean_var

  fit_indices
}


.boot_se_ci <- function(fit_target, L_rot, boot_fit, boot_rot, ci, b) {

  l_ci <- (1 - ci) / 2
  ps <- c(l_ci, ci + l_ci)

  ### calculate stats for unrot loadings and gof measures
  L_unrot <- fit_target$unrot_loadings

  ncol_L <- ncol(L_unrot)
  nrow_L <- nrow(L_unrot)
  colnam_L <- colnames(L_unrot)
  rownam_L <- rownames(L_unrot)

  L_unrot_boot <- array(NA_real_, c(nrow_L, ncol_L, b),
                        dimnames = list(rownam_L, colnam_L,
                                        NULL))
  # The bootstrap aggregates only the numeric fit indices: a scaled chi-square fit (the
  # cor_method = "fiml" two-stage statistic, or any se = "sandwich" fit) adds the character
  # `chi_scaled_type` tag and the extra scaled-statistic components, which the quantile/SD
  # aggregation below cannot consume. Aligning each replicate to these names also lets a
  # replicate whose scaled statistic degenerated (and so lacks a component) contribute NA
  # there rather than shifting every column.
  gof_names <- names(fit_target$fit_indices)[
    vapply(fit_target$fit_indices, is.numeric, logical(1))]
  gof_boot <- matrix(NA_real_, ncol = length(gof_names), nrow = b)
  residuals_boot <- array(NA_real_, c(nrow_L, nrow_L, b),
                          dimnames = list(rownam_L, rownam_L,
                                          NULL))

  # Track replicates that could not be fit (NULL) or aligned, so they are
  # excluded from the bootstrap statistics rather than aborting the whole call.
  failed <- vapply(boot_fit, is.null, logical(1))

  for (boot_i in seq_len(b)) {

    if (failed[boot_i]) next

    # save aligned loading matrix
    aligned <- tryCatch(
      .align_solution(L_unrot, boot_fit[[boot_i]]$unrot_loadings),
      error = function(e) NULL
    )
    if (is.null(aligned)) {
      failed[boot_i] <- TRUE
      next
    }

    L_unrot_boot[,, boot_i] <- aligned$loadings
    fi_i <- boot_fit[[boot_i]]$fit_indices
    gof_boot[boot_i, ] <- vapply(gof_names, function(nm) {
      v <- fi_i[[nm]]
      if (is.null(v) || !is.numeric(v)) NA_real_ else as.numeric(v[1L])
    }, numeric(1))
    residuals_boot[,, boot_i] <- boot_fit[[boot_i]]$residuals

  }

  n_failed <- sum(failed)
  if (n_failed == b) {
    cli::cli_abort(
      c("All {b} bootstrap replicates failed; no bootstrap standard errors could be computed.",
        "i" = "The resampled correlation matrices may be degenerate; try more observations or fewer factors."),
      class = "efa_boot_all_failed"
    )
  }
  if (n_failed > 0) {
    cli::cli_warn(
      c("{n_failed} bootstrap replicate{?s} failed and {?was/were} excluded.",
        "i" = "Bootstrap standard errors and confidence intervals are based on {b - n_failed} replicate{?s}."),
      class = "efa_boot_replicate_failed"
    )
  }

  # se = sd of bootstrap replications (Zhang, 2014, Estimating Standard Errors
  # in Exploratory Factor Analysis)
  L_unrot_se_ci <- .array_se_ci(L_unrot_boot, ps)
  gof_se_ci <- .array_se_ci(gof_boot, ps, M = 2)
  residuals_se_ci <- .array_se_ci(residuals_boot, ps)
  names(gof_se_ci$se) <- gof_names
  names(gof_se_ci$ci$lower) <- gof_names
  names(gof_se_ci$ci$upper) <- gof_names

  if(boot_rot == "oblique") {

    colnam_L <- colnames(L_rot)

    L_rot_boot <- array(NA_real_, c(nrow_L, ncol_L, b),
                        dimnames = list(rownam_L, colnam_L,
                                        NULL))
    Phi_rot_boot <- array(NA_real_, c(ncol_L, ncol_L, b),
                          dimnames = list(colnam_L, colnam_L,
                                          NULL))
    Structure_boot <- array(NA_real_, c(nrow_L, ncol_L, b),
                        dimnames = list(rownam_L, colnam_L,
                                        NULL))

    failed_rot <- 0

    # Align every successfully fit replicate to the target in a single compiled
    # call over the loading cube, rather than one PROCRUSTES() round trip per
    # replicate. Replicates that failed to fit or could not be aligned to the
    # point estimate are excluded from the cube and stay NA in the output arrays.
    keep <- which(!failed)
    if (length(keep) > 0) {
      A_cube <- array(NA_real_, c(nrow_L, ncol_L, length(keep)))
      for (j in seq_along(keep)) {
        A_cube[, , j] <- boot_fit[[keep[j]]]$unrot_loadings
      }

      aligned <- .oblique_procrustes_batch(A_cube, L_rot, random_starts = 5)

      for (j in seq_along(keep)) {
        boot_i <- keep[j]
        # Exclude a replicate only when no valid alignment could be produced. The
        # best multi-start fit is kept even if it did not formally converge: it is
        # the lowest-objective alignment available and its loadings are well-defined.
        if (!isTRUE(aligned$valid[j])) {
          failed_rot <- failed_rot + 1
          next
        }
        L_j <- matrix(aligned$loadings[, , j], nrow_L, ncol_L)
        Phi_j <- matrix(aligned$Phi[, , j], ncol_L, ncol_L)
        L_rot_boot[, , boot_i] <- L_j
        Phi_rot_boot[, , boot_i] <- Phi_j
        Structure_boot[, , boot_i] <- L_j %*% Phi_j
      }
    }

    valid_rot <- b - n_failed - failed_rot
    if (failed_rot > 0) {
      cli::cli_warn(c("{failed_rot} target rotation{?s} in the bootstrap procedure could not be aligned.",
                    "i" = "Bootstrap SE and CI of rotated loadings, factor correlations and structure coefficients are based on {valid_rot} bootstrap sample{?s}."))
    }

    L_rot_se_ci <- .array_se_ci(L_rot_boot, ps)
    Phi_rot_se_ci <- .array_se_ci(Phi_rot_boot, ps)
    Structure_se_ci <- .array_se_ci(Structure_boot, ps)

    out <- list(
      SE = list(
        unrot_loadings = L_unrot_se_ci$se,
        rot_loadings = L_rot_se_ci$se,
        Phi = Phi_rot_se_ci$se,
        Structure = Structure_se_ci$se,
        fit_indices = gof_se_ci$se,
        residuals = residuals_se_ci$se,
        valid_target_rotations = valid_rot
      ),
      CI = list(
        unrot_loadings = L_unrot_se_ci$ci,
        rot_loadings = L_rot_se_ci$ci,
        Phi = Phi_rot_se_ci$ci,
        Structure = Structure_se_ci$ci,
        fit_indices = gof_se_ci$ci,
        residuals = residuals_se_ci$ci
      ),
      replicates = list(
        unrot_loadings = L_unrot_boot,
        rot_loadings = L_rot_boot,
        Phi = Phi_rot_boot,
        Structure = Structure_boot,
        fit_indices = gof_boot,
        residuals = residuals_boot
      )
    )

  } else if (boot_rot == "orthogonal") {
    colnam_L <- colnames(L_rot)

    L_rot_boot <- array(NA_real_, c(nrow_L, ncol_L, b),
                        dimnames = list(rownam_L, colnam_L,
                                        NULL))

    failed_rot <- 0

    for (boot_i in seq_len(b)) {

      if (failed[boot_i]) next

      # save target-rotated loading matrix
      aligned_i <- tryCatch(
        PROCRUSTES(boot_fit[[boot_i]]$unrot_loadings,
                   Target = L_rot, rotation = "orthogonal"),
        error = function(e) NULL
      )
      if (is.null(aligned_i)) {
        failed_rot <- failed_rot + 1
        next
      }
      L_rot_boot[,, boot_i] <- aligned_i$loadings
    }

    valid_rot <- b - n_failed - failed_rot
    if (failed_rot > 0) {
      cli::cli_warn(c("{failed_rot} target rotation{?s} in the bootstrap procedure could not be aligned.",
                    "i" = "Bootstrap SE and CI of rotated loadings are based on {valid_rot} bootstrap sample{?s}."))
    }

    L_rot_se_ci <- .array_se_ci(L_rot_boot, ps)

    out <- list(
      SE = list(
        unrot_loadings = L_unrot_se_ci$se,
        rot_loadings = L_rot_se_ci$se,
        fit_indices = gof_se_ci$se,
        residuals = residuals_se_ci$se,
        valid_target_rotations = valid_rot
      ),
      CI = list(
        unrot_loadings = L_unrot_se_ci$ci,
        rot_loadings = L_rot_se_ci$ci,
        fit_indices = gof_se_ci$ci,
        residuals = residuals_se_ci$ci
      ),
      replicates = list(
        unrot_loadings = L_unrot_boot,
        rot_loadings = L_rot_boot,
        fit_indices = gof_boot,
        residuals = residuals_boot
      )
    )


  } else {

    out <- list(
      SE = list(
        unrot_loadings = L_unrot_se_ci$se,
        fit_indices = gof_se_ci$se,
        residuals = residuals_se_ci$se
      ),
      CI = list(
        unrot_loadings = L_unrot_se_ci$ci,
        fit_indices = gof_se_ci$ci,
        residuals = residuals_se_ci$ci
      ),
      replicates = list(
        unrot_loadings = L_unrot_boot,
        fit_indices = gof_boot,
        residuals = residuals_boot
      )
    )

  }


  out

}

.array_se_ci <- function(x, probs, M = c(1, 2)) {
  se <- apply(x, M, stats::sd, na.rm = TRUE)
  ci <- lapply(probs, function(p) {
    apply(x, M, stats::quantile, probs = p, na.rm = TRUE)
  })
  names(ci) <- c("lower", "upper")

  list(
    se = se,
    ci = ci
  )
}
