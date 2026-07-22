# Exploratory factor analysis (EFA)

This function does an EFA with either `PAF`, `ML`, `ULS`/`MINRES`, or
`DWLS` with or without subsequent rotation. All arguments with default
value `NA` can be left to default if `type` is set to one of "EFAtools",
"SPSS", or "psych". The respective specifications are then handled
according to the specified type (see details).

## Usage

``` r
efa_fit(
  x,
  n_factors,
  N = NA,
  estimator = c("PAF", "ML", "ULS", "MINRES", "DWLS"),
  rotation = c("none", "varimax", "equamax", "quartimax", "geominT", "bentlerT",
    "bifactorT", "promax", "oblimin", "quartimin", "simplimax", "bentlerQ", "geominQ",
    "bifactorQ"),
  se = c("none", "information", "sandwich", "np-boot"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra", "fiml"),
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  estimate_control = NULL,
  rotate_control = NULL,
  b_boot = 1000,
  ci = 0.95,
  seed = NULL,
  ...
)
```

## Source

Bollen, K. A., & Stine, R. A. (1992). Bootstrapping goodness-of-fit
measures in structural equation models. Sociological Methods & Research,
21, 205–229. doi: 10.1177/0049124192021002004

Grieder, S., & Steiner, M. D. (2022). Algorithmic jingle jungle: A
comparison of implementations of principal axis factoring and promax
rotation in R and SPSS. Behavior Research Methods, 54, 54–74. doi:
10.3758/s13428-021-01581-x

Hendrickson, A. E., & White, P. O. (1964). Promax: A quick method for
rotation to oblique simple structure. British Journal of Statistical
Psychology, 17 , 65–70. doi: 10.1111/j.2044-8317.1964.tb00244.x

Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. L. (2011). The Hull
Method for Selecting the Number of Common Factors, Multivariate
Behavioral Research, 46, 340-364, doi: 10.1080/00273171.2011.564527

Kaiser, H. F. (1958). The varimax criterion for analytic rotation in
factor analysis. Psychometrika, 23, 187–200. doi: 10.1007/BF02289233

Lawley, D. N., & Maxwell, A. E. (1971). Factor analysis as a statistical
method (2nd ed.). Butterworths.

Cudeck, R. (1989). Analysis of correlation matrices using covariance
structure models. Psychological Bulletin, 105, 317–327. doi:
10.1037/0033-2909.105.2.317

Olkin, I., & Siotani, M. (1976). Asymptotic distribution of functions of
a correlation matrix. In S. Ikeda (Ed.), Essays in probability and
statistics (pp. 235–251). Shinko Tsusho.

Jennrich, R. I. (1973). Standard errors for obliquely rotated factor
loadings. Psychometrika, 38, 593–604. doi: 10.1007/BF02291497

Zhang, G., & Preacher, K. J. (2015). Factor rotation and standard errors
in exploratory factor analysis. Journal of Educational and Behavioral
Statistics, 40, 579–603. doi: 10.3102/1076998615606098

Browne, M. W. (1984). Asymptotically distribution-free methods for the
analysis of covariance structures. British Journal of Mathematical and
Statistical Psychology, 37, 62–83. doi:
10.1111/j.2044-8317.1984.tb00789.x

Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics and
standard errors in covariance structure analysis. In A. von Eye & C. C.
Clogg (Eds.), Latent variables analysis (pp. 399–419). Sage.

Asparouhov, T., & Muthén, B. (2010). Simple second order chi-square
correction. Mplus Technical Appendix.

Muthén, B. (1984). A general structural equation model with dichotomous,
ordered categorical, and continuous latent variable indicators.
Psychometrika, 49, 115–132. doi: 10.1007/BF02294210

Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods for
mean and covariance structure analysis with nonnormal missing data.
Sociological Methodology, 30, 165–200. doi: 10.1111/0081-1750.00078

Yuan, K.-H., Marshall, L. L., & Bentler, P. M. (2002). A unified
approach to exploratory factor analysis with missing data, nonnormal
data, and in the presence of outliers. Psychometrika, 67, 95–121. doi:
10.1007/BF02294711

Savalei, V., & Bentler, P. M. (2009). A two-stage approach to missing
data: Theory and application to auxiliary variables. Structural Equation
Modeling, 16, 477–497. doi: 10.1080/10705510903008238

Little, R. J. A., & Rubin, D. B. (2002). Statistical analysis with
missing data (2nd ed.). Wiley.

Bartlett, M. S. (1951). The effect of standardization on approximation
in factor analysis. Biometrika, 38, 337–344.

Bentler, P. M. (1990). Comparative fit indexes in structural models.
Psychological Bulletin, 107, 238–246. doi: 10.1037/0033-2909.107.2.238

Tucker, L. R., & Lewis, C. (1973). A reliability coefficient for maximum
likelihood factor analysis. Psychometrika, 38, 1–10. doi:
10.1007/BF02291170

Browne, M. W., & Cudeck, R. (1989). Single sample cross-validation
indices for covariance structures. Multivariate Behavioral Research, 24,
445–455. doi: 10.1207/s15327906mbr2404_4

Browne, M. W., & Cudeck, R. (1992). Alternative ways of assessing model
fit. Sociological Methods & Research, 21, 230–258. doi:
10.1177/0049124192021002005

Bentler, P. M. (1995). EQS structural equations program manual.
Multivariate Software.

## Arguments

- x:

  data.frame or matrix. Dataframe or matrix of raw data or matrix with
  correlations. If raw data is entered, the correlation matrix is found
  from the data.

- n_factors:

  numeric. Number of factors to extract. Must be at least 1 and smaller
  than the number of variables (the common factor model is not
  identified otherwise). Use
  [`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
  to decide on a value.

- N:

  numeric. The number of observations. Needs only be specified if a
  correlation matrix is used. If input is a correlation matrix and `N` =
  NA (default), not all fit indices can be computed. When raw data with
  missing values are entered and `use` is `"complete.obs"` or
  `"na.or.complete"`, rows are deleted listwise, so `N` is taken as the
  number of complete cases.

- estimator:

  character. The estimator used to fit the EFA: "PAF" (principal axis
  factoring), "ML" (maximum likelihood), "ULS" (unweighted least
  squares; "MINRES" is an accepted alias returning identical results),
  or "DWLS" (diagonally weighted least squares, for ordinal data). See
  the *Estimators* section in Details for their properties and data
  requirements. The value is matched case-insensitively.

- rotation:

  character. Either perform no rotation ("none"; default), an orthogonal
  rotation ("varimax", "equamax", "quartimax", "geominT", "bentlerT", or
  "bifactorT"), or an oblique rotation ("promax", "oblimin",
  "quartimin", "simplimax", "bentlerQ", "geominQ", or "bifactorQ"). See
  the *Rotations* section in Details for their properties and known
  issues.

- se:

  character. Whether and how to compute standard errors (and matching
  confidence intervals): "none" (default, no standard errors),
  "information" (analytic standard errors from the expected Fisher
  information of the ML solution), "sandwich" (robust Godambe sandwich
  standard errors from raw data), or "np-boot" (non-parametric
  bootstrap). The methods differ in their assumptions, their data
  requirements, and which estimator, rotation, and `cor_method`
  combinations they support; see the *Standard errors* section in
  Details.

- cor_method:

  character. How the correlation is computed from raw data: `"pearson"`,
  `"spearman"`, or `"kendall"` (passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)); `"poly"` /
  `"tetra"` for polychoric / tetrachoric correlations of ordinal /
  binary data; or `"fiml"` for a two-stage full-information
  maximum-likelihood correlation from raw data with missing values. See
  the *Correlation methods* section in Details for their properties and
  the combinations they support. Default is "pearson".

- use:

  character. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) if raw data is
  given as input. Default is "pairwise.complete.obs".

- estimate_control:

  a control object from
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  bundling the estimation tuning knobs: the `type` preset; the
  principal-axis-factoring iteration settings `init_comm`, `criterion`,
  `criterion_type`, `max_iter`, and `abs_eigen`; and the
  maximum-likelihood `start_method`. Defaults to
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md),
  which resolves every preset-driven knob from the `"EFAtools"` type.
  See
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  for the individual knobs and the *Using the type presets* section in
  Details for how the preset fills them.

- rotate_control:

  a control object from
  [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  bundling the rotation tuning knobs: the `type` preset; Kaiser
  `normalize`; the convergence `precision`; the factor `order_type`; the
  varimax/promax settings `varimax_type` and `p_type`; the
  simplimax/promax `k`; and `random_starts`. Defaults to
  [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md),
  which resolves every preset-driven knob from the `"EFAtools"` type.
  The estimation and rotation presets are independent, so
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  and
  [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  may carry different `type`s. See
  [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  for the individual knobs.

- b_boot:

  numeric. The number of bootstrap samples to draw. Default is 1000.
  Under `cor_method = "fiml"` each bootstrap sample re-runs the EM
  moment estimation, so a smaller value may be advisable.

- ci:

  numeric. The confidence interval to create from the bootstrap samples.
  Must be between 0 and 1. Default is .95 for 95% CIs.

- seed:

  numeric. An optional seed for the random-number generator, governing
  every stochastic part of the fit: the rotation's random starts on the
  point estimate (the criterion-based rotations draw `random_starts`
  random starts; see *Rotations*) and, under `se = "np-boot"`, the case
  resampling, the replicate rotations, and the Procrustes random starts.
  Setting it makes the fit reproducible and the bootstrap additionally
  independent of the number of parallel workers (see Details); the
  caller's random-number stream is restored afterwards, so supplying a
  seed leaves no lasting effect on it. Default is `NULL`, which uses
  (and advances) the current state of the generator.

- ...:

  Additional arguments forwarded to the rotation engine. Only the
  arguments the selected `rotation` consumes are accepted: `maxit` (the
  maximum number of engine iterations) for the GPArotation-style
  rotations, plus the selected criterion's parameter (`gam` for
  "oblimin", `delta` for "geominT" and "geominQ"). Anything else – a
  misspelled name, another criterion's parameter, or any extra with
  "varimax", "promax", or "none", which consume no extras – is an error
  rather than a setting that is silently ignored. The accepted arguments
  are merged with, and take precedence over, the extra arguments stored
  in
  [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md).
  An estimation or rotation tuning knob (such as `type`, `max_iter`, or
  `k`) is likewise *not* accepted here: it belongs to
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  or
  [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md).

## Value

A list of class `c("efa", "EFA")` containing (a subset of) the
following:

- orig_R:

  Original correlation matrix.

- h2_init:

  Initial communality estimates from PAF.

- h2:

  Final communality estimates from the unrotated solution.

- orig_eigen:

  Eigen values of the original correlation matrix.

- init_eigen:

  Initial eigenvalues, obtained from the correlation matrix with the
  initial communality estimates as diagonal in PAF.

- final_eigen:

  Eigenvalues obtained from the correlation matrix with the final
  communality estimates as diagonal.

- iter:

  For PAF, the number of iterations until convergence. For ML, ULS, and
  DWLS, the number of objective-function evaluations used by the
  optimiser (not the number of optimiser iterations).

- convergence:

  Integer convergence code (0 = converged), using the codes of
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html). For ML and ULS
  it is the code from the bounded optimiser; for DWLS the fit runs a
  bounded warm start followed by an unconstrained polish, and the
  reported code is from the final polish. For PAF it is 1 if the maximum
  number of iterations was reached without meeting the convergence
  criterion and 0 otherwise. A non-zero code is also reported with a
  warning.

- heywood:

  A named integer vector indicating which variables have a Heywood
  (improper) case in the unrotated solution; empty if there are none.

- unrot_loadings:

  Loading matrix containing the final unrotated loadings.

- vars_accounted:

  Matrix of explained variances and sums of squared loadings. Based on
  the unrotated loadings.

- fit_indices:

  A named list of fit indices computed from the unrotated loadings. For
  ML and ULS it holds the model Chi Square (with its p-value and df),
  CFI, TLI, RMSEA with its 90% confidence interval, AIC, BIC, ECVI,
  RMSR, SRMR, and CAF; for PAF and DWLS only RMSR, SRMR, CAF, and df are
  populated and the Chi-Square-based indices are `NA` (for DWLS with
  `se = "sandwich"` the full block is filled from a scaled Chi Square
  instead). Whenever the Chi Square is a scaled statistic
  (`se = "sandwich"`, or any `cor_method = "fiml"` fit) AIC, BIC, and
  ECVI are `NA` and the list additionally carries the scaling
  components: `chi_scaling` (the multiplier a in the scaled-and-shifted
  statistic \\aT + b\\, i.e. the reciprocal of lavaan's
  `chisq.scaling.factor`), `chi_shift` (b), `chi_unscaled` (the unscaled
  statistic T), and the alternative `chi_mean_adjusted` and
  `chi_mean_var` statistics with their `df_mean_var`. RMSR is retained
  for programmatic use and backward compatibility, although the print
  and summary methods display SRMR. See the *Fit indices* section in
  Details for how each index is defined, scaled, and referenced.

- model_implied_R:

  The model implied correlation matrix.

- residuals:

  Residual correlations, i.e., orig_R - model_implied_R

- standardized_residuals:

  Residual correlations standardized by their bootstrap standard errors.
  Only returned, if `se = "np-boot"`.

- rot_loadings:

  Loading matrix containing the final rotated loadings (pattern matrix).

- Phi:

  The factor intercorrelations (only for oblique rotations).

- Structure:

  The structure matrix (only for oblique rotations).

- rotmat:

  The rotation matrix. The rotated loadings are recovered from the
  unrotated loadings as `unrot_loadings %*% rotmat` for orthogonal
  rotations and for promax, and as `unrot_loadings %*% t(solve(rotmat))`
  for the other oblique rotations.

- vars_accounted_rot:

  Matrix of explained variances and sums of squared loadings. Based on
  rotated loadings and, for oblique rotations, the factor
  intercorrelations.

- settings:

  A list of the settings used. For the criterion rotations fitted by
  gradient projection it additionally carries `rotation_diagnostics`, a
  list summarising the multi-start run: `n_starts_total` (the
  `random_starts` random starts plus the rational start), `n_optimized`
  (how many of those starts were actually optimized – fewer than
  `n_starts_total` whenever the solver screens the random starts and
  optimizes only the most promising ones), `n_converged` (how many
  optimized starts reached the convergence tolerance),
  `n_distinct_minima` (how many distinct local optima those converged
  starts found; more than one means the criterion is multimodal on these
  data), `criterion_spread` (the range of the criterion values they
  attained), and `criterion_best` (the criterion value of the returned
  solution). When `normalize = TRUE`, `criterion_best` and
  `criterion_spread` are evaluated on the Kaiser-normalized loadings the
  criterion is optimized on, not on the returned `rot_loadings`, so they
  are not directly comparable to a criterion recomputed from the
  returned loadings.

- SE:

  A named list of standard error matrices. For `se = "np-boot"`:
  bootstrap standard deviations of the unrotated and (when a rotation is
  applied) rotated loadings, the residuals, and the fit indices, plus –
  for oblique rotations – the factor correlations (`Phi`) and the
  structure coefficients; it additionally carries `valid_replicates`,
  the number of bootstrap replicates that were fitted and aligned
  successfully and that every entry above is therefore based on
  (replicates that failed are excluded and warned about), and, when a
  rotation is applied, `valid_target_rotations`, the number of those
  replicates that could also be aligned to the rotated point estimate
  and that the rotated entries (`rot_loadings`, `Phi`, `Structure`) are
  based on. For `se = "information"`: Wald standard errors from the
  expected (Fisher) information matrix for the unrotated loadings, the
  uniquenesses, and the communalities and, when a rotation is applied,
  the rotated loadings (and, for oblique rotations, `Phi` and the
  structure coefficients). Because \\h^2_i = 1 - \psi_i\\ exactly, the
  communality and uniqueness standard errors are identical. For
  `se = "sandwich"`: robust Godambe sandwich standard errors with the
  same coverage as `"information"`, robust to non-normality and weight
  misspecification. Only returned if `se` is not `"none"`.

- CI:

  A named list of confidence intervals of width `ci`. For
  `se = "np-boot"`: percentile intervals matching the components of
  `SE`. For `se = "information"` and `se = "sandwich"`: Wald intervals
  matching the components of `SE`. Only returned if `se` is not
  `"none"`.

- replicates:

  A named list of bootstrap replicate arrays for the aligned unrotated
  and (where applicable) rotated loadings, structure coefficients,
  factor correlations (`Phi`), residuals, and fit indices. The replicate
  is the last dimension of the loading, residual, `Phi`, and structure
  cubes, and the first dimension of the `fit_indices` matrix (whose
  columns are named after the fit indices). Replicates that failed are
  left `NA`. Populated only for `se = "np-boot"`; `NULL` for the
  analytic SE methods.

- vcov_unrot_loadings:

  The full unrotated loading covariance matrix the marginal
  `SE$unrot_loadings` were derived from: a `p * n_factors` by
  `p * n_factors` numeric matrix in column-major `vec(Lambda)` order.
  Populated for `se = "information"` (expected-information block) and
  `se = "sandwich"` (robust V_AA), even when a rotation is applied (the
  persisted block is always the unrotated one); NA-filled if the
  analytic covariance is unreliable (a Heywood case or a singular
  bordered information matrix); `NULL` for `se = "np-boot"` and
  `se = "none"`. A weakly determined rotational orientation is the one
  case where this matrix is populated while `SE$unrot_loadings` is `NA`:
  the covariance itself is finite and valid, and only its
  gauge-dependent marginals are not (see *Standard errors*).

- Gamma:

  The asymptotic covariance of the off-diagonal sample correlations –
  the meat of the robust sandwich SEs – on the variance scale
  (`Var(rho-hat)`; lavaan's correlation NACOV is `N * Gamma`). A
  `p (p - 1) / 2` by `p (p - 1) / 2` numeric matrix; rows and columns
  ordered by [`utils::combn()`](https://rdrr.io/r/utils/combn.html) over
  the column pairs and labelled `"<var_i>-<var_j>"`. Populated for
  `se = "sandwich"` on the polychoric/tetrachoric and Pearson paths;
  `NULL` otherwise, including under `cor_method = "fiml"`, whose meat is
  the saturated FIML asymptotic covariance and is not returned.

## Details

There are two main ways to use this function. The easiest way is to use
it with a specified `type` (see above), which sets most of the other
arguments accordingly. Another way is to use it more flexibly by
explicitly specifying all arguments used and set `type` to "none" (see
examples). A mix of the two can also be done by specifying a `type` as
well as additional arguments. However, this will throw warnings to avoid
unintentional deviations from the implementations according to the
specified `type`.

### Estimators

The estimator is chosen with `estimator`.

- **PAF** (principal axis factoring) iteratively estimates the
  communalities and makes no distributional assumptions, which makes it
  robust and a good general-purpose default. Because it minimises no
  likelihood or weighted discrepancy it provides no model chi-square,
  and hence no chi-square-based fit indices (see *Fit indices*). The PAF
  iteration is governed by `init_comm`, `criterion`, `criterion_type`,
  `max_iter`, and `abs_eigen` (set by `type`; see *Using the type
  presets*).

- **ML** (maximum likelihood) maximises the normal-theory likelihood. It
  yields the full set of fit indices and is the only estimator with
  analytic expected-information standard errors (`se = "information"`),
  but it assumes multivariate normality and is the most prone to Heywood
  (improper) cases. Its starting values are set by `start_method`.

- **ULS** (unweighted least squares) minimises the sum of squared
  correlation residuals. "MINRES" (minimum residual) is the same
  estimator under a different name and returns identical results. It
  makes no normality assumption, is robust to mild non-normality, and
  yields the full set of fit indices.

- **DWLS** (diagonally weighted least squares) is the recommended
  estimator for ordinal data. It weights each off-diagonal correlation
  residual by the inverse asymptotic variance of the corresponding
  polychoric correlation (Muthén, 1984), reproducing the loadings of a
  diagonally weighted least squares fit (e.g.
  `lavaan::efa(..., estimator = "DWLS")`). It therefore requires raw
  ordinal data with `cor_method = "poly"` or `"tetra"` and has no
  fallback for a supplied correlation matrix or a continuous
  `cor_method`. Because the weighting follows the polychoric asymptotic
  covariance, the matrix and the weights are estimated on the
  listwise-complete cases. Its fit-index behaviour is described under
  *Fit indices*.

### Correlation methods

When raw data are supplied, `cor_method` selects how the correlation
matrix is computed (it is ignored when a correlation matrix is entered
directly).

- **"pearson"** (default), **"spearman"**, and **"kendall"** are passed
  to [`stats::cor()`](https://rdrr.io/r/stats/cor.html) for continuous
  or rank data.

- **"poly"** / **"tetra"** compute polychoric / tetrachoric correlations
  for ordinal / binary data, assuming an underlying bivariate-normal
  latent variable. They use a two-step estimator with no empty-cell
  continuity correction, matching
  [`polycor::polychor()`](https://rdrr.io/pkg/polycor/man/polychor.html)
  and `lavaan`. The polychoric asymptotic covariance that underlies both
  the DWLS weights and the scaled (sandwich) statistic relies on
  large-sample theory that degrades for empty or near-empty
  response-category combinations; with very sparse cells the resulting
  weights and standard errors can be unreliable (a warning is issued
  when empty cells are present), so interpret them with caution and
  consider collapsing rare categories.

- **"fiml"** estimates a two-stage full-information maximum-likelihood
  correlation. The saturated multivariate-normal mean and covariance are
  estimated from raw data with missing values by an EM algorithm
  assuming the data are missing at random (Yuan, Marshall, & Bentler,
  2002; Little & Rubin, 2002), and the standardized covariance is then
  analysed. This reproduces
  [`psych::corFiml()`](https://rdrr.io/pkg/psych/man/corFiml.html)
  followed by [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html) and
  `lavaan(missing = "two.stage")`, *not* `lavaan::efa(missing = "ml")`,
  so the point estimates are not expected to match the latter. The model
  fit indices are corrected two-stage statistics (see *Fit indices*).
  `"fiml"` uses every case and handles the missingness itself, so `use`
  is ignored; it supplies a continuous (Pearson-type) correlation only
  and is therefore not compatible with `estimator = "DWLS"`. Standard
  errors are available analytically for `estimator = "ML"` or `"ULS"`
  and, for any estimator, by the non-parametric bootstrap (see *Standard
  errors*). For multiply imputed data,
  [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)
  is the alternative route to handling missingness.

### Rotations

A rotation transforms the unrotated loadings toward a simpler, more
interpretable pattern; all rotations are performed by rotation engines
built into the package. Orthogonal rotations keep the factors
uncorrelated, whereas oblique rotations let them correlate (returning a
pattern matrix, a structure matrix, and the factor intercorrelations
`Phi`) and are usually more realistic for psychological constructs.

Orthogonal rotations:

- **varimax** maximises the variance of the squared loadings within each
  factor (column simplicity). It is the most widely used orthogonal
  rotation and spreads variance across factors rather than concentrating
  it in a general factor.

- **quartimax** simplifies the variables (rows) so that each loads
  mainly on one factor; it tends to produce a strong general factor.

- **equamax** is a Crawford-Ferguson compromise between varimax (column)
  and quartimax (row) simplicity.

- **geominT** uses a geometric-mean criterion that rewards a sparse
  pattern and tolerates variables with cross-loadings; a smaller offset
  `delta` gives a sparser solution but sharper local minima.

- **bentlerT** uses Bentler's invariant pattern simplicity criterion.

- **bifactorT** is the Jennrich-Bentler orthogonal bifactor criterion: a
  general factor plus group factors (bifactor simple structure).

Oblique rotations:

- **promax** is a fast two-step rotation: a varimax solution is raised
  to a power (controlled by `k` and `p_type`) to form a target that is
  then fitted obliquely. It is the common, inexpensive oblique default.

- **oblimin** is a flexible oblique family controlled by `gam` (default
  0); a good general-purpose criterion.

- **quartimin** is oblimin pinned at `gam = 0`; a robust default oblique
  criterion.

- **simplimax** drives the `k` smallest loadings toward zero. Its
  criterion is only piecewise smooth, so it is the most prone to local
  minima and relies on several `random_starts`.

- **bentlerQ** is the oblique Bentler invariant pattern simplicity
  criterion.

- **geominQ** is the oblique geomin criterion; it handles complex
  (cross-loading) structure well but is multimodal, so it benefits from
  more `random_starts` (and uses a more thorough multi-start search
  internally).

- **bifactorQ** is the oblique (correlated) Jennrich-Bentler bifactor
  criterion.

The criterion-based rotations (all except varimax and promax) are fitted
by gradient projection with `random_starts` random starts to guard
against local minima; the complexity criteria (simplimax and geominQ in
particular) are the most multimodal. The starts are drawn from the
random-number generator, so different starts can reach genuinely
different optima and such a fit is reproducible only when the generator
is controlled: pass `seed`, or call
[`set.seed()`](https://rdrr.io/r/base/Random.html) beforehand. The
`type` argument changes the varimax and promax settings (see *Using the
type presets*) and, for every rotation, the factor `order_type`. A
single factor cannot be rotated.

### Standard errors

`se` selects whether and how standard errors (and matching confidence
intervals) are computed. Which quantities they cover depends on the
method. The analytic methods (`"information"` and `"sandwich"`) cover
the unrotated loadings, the uniquenesses and the communalities and, when
a rotation is applied, the rotated loadings and – for oblique rotations
– the factor correlations and the structure coefficients. The bootstrap
(`"np-boot"`) covers the unrotated loadings, the residuals, and the fit
indices and, when a rotation is applied, the rotated loadings and – for
oblique rotations – the factor correlations and the structure
coefficients; it reports no uniqueness or communality standard errors
(see the `SE` and `CI` entries in Value).

- **"none"** (default) computes no standard errors.

- **"information"** returns analytic standard errors from the expected
  (Fisher) information matrix of the maximum-likelihood solution, and
  therefore requires `estimator = "ML"` and `cor_method = "pearson"` (or
  `"fiml"`, see below). The rotated standard errors are obtained by
  propagating the unrotated-loading covariance through the rotation by
  the delta method (Jennrich, 1973); because rotated quantities are
  identification-invariant they are directly comparable across programs.
  Unlike the bootstrap it also works from a correlation matrix as long
  as `N` is supplied.

  `efa_fit()` analyses a *correlation* structure – the diagonal of the
  analysed matrix is fixed at 1 and the uniquenesses \\\psi_i = 1 -
  \sum_j \lambda\_{ij}^2\\ are a function of the loadings rather than
  free parameters – so the information is the correlation-structure one,
  \\\Delta' \Gamma^{-1} \Delta\\ over the off-diagonal correlations,
  with \\\Gamma\\ the normal-theory asymptotic covariance of the sample
  correlations (Olkin & Siotani, 1976) at the model-implied \\\Sigma\\
  (Cudeck, 1989). It is scaled by \\1 / (N - 1)\\ and the confidence
  intervals are Wald intervals (estimate \\\pm\\ z \* SE). These
  standard errors assume multivariate normality and a correctly
  specified model; under heavy-tailed data or model misfit they can
  understate the sampling variability, where `"sandwich"` or `"np-boot"`
  are more robust.

  The *rotated* quantities, the uniquenesses, and the communalities are
  identification-invariant and so are comparable across programs. The
  **unrotated** loading standard errors are not: they are reported in
  the orientation the solution itself is identified in (for ML,
  \\\Lambda' \Psi^{-1} \Lambda\\ diagonal), and a program using a
  different identification will report different unrotated loading
  standard errors for the same fit.

  That orientation can also fail to be determined. When two of the
  canonical variances \\\mathrm{diag}(\Lambda' \Psi^{-1} \Lambda)\\
  nearly coincide – which happens with a weakly determined factor, or
  two factors of near-equal strength – the loadings can be rotated
  within that plane without leaving the identifying constraint, so the
  unrotated loadings have no well-defined sampling distribution and
  their standard errors diverge. `efa_fit()` detects this and returns
  `NA` for the unrotated loading standard errors with an
  `efa_se_unreliable` warning, rather than reporting a number that looks
  like a standard error but is an artefact of the orientation.
  Everything that does not depend on the orientation – the rotated
  loadings, `Phi`, the structure coefficients, the uniquenesses and the
  communalities – is unaffected and still reported. Use those, or
  `se = "np-boot"`, when the unrotated loadings themselves are the
  quantity of interest. The detection covers the Pearson and polychoric
  paths; the two-stage `cor_method = "fiml"` sandwich carries no such
  diagnostic, so a weakly determined orientation is not flagged there. A
  Heywood case (a uniqueness at its lower boundary) is separate: the
  Wald approximation fails there for every parameter, so no analytic
  standard error is reported at all.

- **"sandwich"** returns robust (Godambe sandwich) standard errors from
  raw data, combining the estimator weight with an
  asymptotic-distribution-free covariance of the correlations, so it
  stays valid under non-normality and weight misspecification (Browne,
  1984; Satorra & Bentler, 1994). It is available either for ordinal
  data with `cor_method = "poly"` or `"tetra"` and `estimator` one of
  `"ML"`, `"ULS"`, or `"DWLS"` (the meat is the polychoric / tetrachoric
  asymptotic covariance), or for continuous data with
  `cor_method = "pearson"` and `estimator = "ML"` or `"ULS"` (the meat
  is the fourth-moment ADF covariance of the sample correlations, the
  basis of the MLM / MLR robust statistics). It reports the same
  coverage as `"information"`, propagated by the same delta method, and
  additionally fills the model fit's chi-square block with a scaled
  (Satorra-Bentler / scaled-and-shifted) chi-square (see *Fit indices*).
  The statistic reported as `chi` is always the scaled-and-shifted one
  (the WLSMV default, flagged by `chi_scaled_type`), for the continuous
  Pearson path as well as the ordinal one; the mean-adjusted statistic
  that lavaan's `MLM` reports for continuous data is returned alongside
  it as `chi_mean_adjusted`. Because the asymptotic covariance must
  describe the same cases as the correlation matrix, the sandwich (like
  `estimator = "DWLS"`) is computed on the listwise-complete cases; on
  data with missing values the reported `N`, the correlation matrix, and
  the point estimate therefore reflect the complete cases regardless of
  `use`.

- **"np-boot"** draws a non-parametric (case-resampling) bootstrap and
  needs raw data. A correlation matrix carries no cases to resample;
  alone among the unsupported combinations this one does not error but
  warns (condition class `"efa_boot_cormat"`) and downgrades `se` to
  `"none"`, so the fit is returned without an `SE` slot. It is the most
  general method – available for any `estimator`, `rotation`, and
  `cor_method` – and the most robust to non-normality and misfit, at the
  cost of speed; its intervals are bootstrap percentile intervals. The
  replicate fits are run across replicates with the `future` framework.
  By default they run sequentially; to run them in parallel, register a
  plan with
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  (e.g. `future::plan(future::multisession, workers = 2)`; see
  examples). With a fixed `seed` the bootstrap is reproducible and
  yields the same result regardless of the number of workers. Under
  `cor_method = "fiml"` each resample also re-runs the EM moment
  estimation and is therefore slow, so a smaller `b_boot` may be
  advisable.

  The percentile intervals are centred on the point estimate for the
  loadings, the factor correlations, the structure coefficients and the
  residuals, but **not** for the indices derived from the chi-square
  (`RMSEA`, `AIC`, `BIC` and `ECVI`). A resample is drawn from the
  sample, which already carries the model's misfit, and is then refitted
  and evaluated against itself, so the replicate chi-square is on
  average about `chi + df` rather than `chi` and the whole interval
  rides upward with it – often far enough that the point estimate falls
  below its own lower bound. That is the shift, not a miscomputed
  interval: read those intervals as a spread rather than as a range for
  the point estimate. Correcting the location needs resampling from a
  population transformed to fit the model (Bollen & Stine, 1992), which
  is not what this bootstrap does. `CFI` and `TLI` are not affected,
  being ratios in which the baseline chi-square shifts along with the
  model one.

The analytic methods (`"information"` and `"sandwich"`) are not
available with the `"promax"` or `"simplimax"` rotations, which have no
usable analytic rotation Jacobian; use `"np-boot"` there. Under
`cor_method = "fiml"`, `"information"` and `"sandwich"` instead return,
for `estimator = "ML"` or `"ULS"`, the corrected two-stage (Yuan &
Bentler, 2000; Savalei & Bentler, 2009) sandwich standard errors, built
on the saturated FIML asymptotic covariance with the estimator's own
Stage-2 weight: the model is fitted to the EM-estimated correlation, so
the naive Stage-2 standard errors (treating that correlation as complete
data) are inconsistent under missingness and are not reported
(`estimator = "PAF"` carries no Stage-2 weight, so use `se = "np-boot"`
there).

### Fit indices

For ML and ULS, `efa_fit()` returns the model chi-square (with its
p-value and degrees of freedom), the Comparative Fit Index (CFI;
Bentler, 1990), the Tucker-Lewis Index (TLI, also called the non-normed
fit index; Tucker & Lewis, 1973), the Root Mean Square Error of
Approximation (RMSEA) with its 90% confidence interval (Browne & Cudeck,
1992), the Akaike and Bayesian Information Criteria (AIC, BIC), the
Expected Cross-Validation Index (ECVI; Browne & Cudeck, 1989), the Root
Mean Squared Residual (RMSR), the Standardized Root Mean Squared
Residual (SRMR; Bentler, 1995), and the common-part-accounted-for (CAF)
index (Lorenzo-Seva, Timmerman, & Kiers, 2011). The print and summary
methods show SRMR, not RMSR, because the two residual summaries differ
only by the fixed scaling \\\sqrt{(p - 1) / (p + 1)}\\ for a fixed
number of variables; RMSR remains in the returned object. The model
chi-square is the Bartlett-corrected discrepancy (matching
[`stats::factanal()`](https://rdrr.io/r/stats/factanal.html) for ML);
the AIC, BIC, and ECVI are the minimum-fit-function (chi-square-based)
forms (\\\chi^2 - 2\\df\\ and \\\chi^2 - \log(N)\\df\\ for AIC and BIC,
as in [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html)) and can
therefore be negative. The RMSEA, CFI, and TLI place the model and
baseline noncentrality on the uncorrected \\N - 1\\ discrepancy scale on
which these approximate-fit indices are defined, so the Bartlett
small-sample correction enters only the chi-square test, not the
approximate-fit indices.

Which indices are reported depends on the estimator:

- **ML and ULS** compute the full set above.

- **PAF** returns only the descriptive residual indices (RMSR, SRMR,
  CAF) and df; the printed model-fit block shows CAF and SRMR. The
  chi-square-based indices are `NA`, because PAF minimises no
  discrepancy.

- **DWLS** by default returns only RMSR, SRMR, CAF, and df, because the
  ordinary maximum-likelihood discrepancy is not its fit function. When
  `se = "sandwich"`, a scaled (Satorra & Bentler, 1994; Asparouhov &
  Muthén, 2010) chi-square and the CFI, TLI, and RMSEA derived from it
  are reported (AIC and BIC remain `NA`). That scaled statistic is a
  two-stage correction applied to the polychoric-correlation residuals
  (Browne, 1984), so it is not identical to the full WLSMV test of
  lavaan or Mplus, which also projects the thresholds.

- **`cor_method = "fiml"`** (with ML or ULS) reports
  Satorra-Bentler-corrected two-stage statistics (Yuan, Marshall, &
  Bentler, 2002): the normal-theory discrepancy on the EM-estimated
  correlation, rescaled by the saturated FIML asymptotic covariance,
  because the plain two-stage likelihood-ratio statistic is not
  asymptotically \\\chi^2(df)\\. The CFI, TLI, and RMSEA follow from the
  scaled statistics; AIC, BIC, and ECVI are left `NA`, as for any scaled
  (moment-adjusted) chi-square.

Whenever the chi-square is a scaled one (`se = "sandwich"`, or any
`cor_method = "fiml"` fit), the AIC, BIC, and ECVI are `NA` and the
returned `fit_indices` additionally carry the scaled-statistic
components (see the `fit_indices` entry in Value). Note that
Lorenzo-Seva, Timmerman, and Kiers (2011) introduce the CAF as ranging
between 0 and 1, with values close to 1 indicating close fit; this does
not match the formula they apply, \\1 - KMO(residuals)\\, which only
works if the diagonal of the residual matrix is set to 1s and then
approximates 0.5 with close fit.

### Available combinations

Not every estimator, rotation, standard-error, and correlation method
can be combined:

- **Estimator and correlation method.** `estimator = "DWLS"` requires
  ordinal data with `cor_method = "poly"` or `"tetra"`.
  `cor_method = "fiml"` works with PAF, ML, and ULS (not DWLS) and needs
  raw data with missing values.

- **Standard errors.** `se = "information"` requires `estimator = "ML"`
  and `cor_method = "pearson"` or `"fiml"`, and can be computed from a
  correlation matrix when `N` is supplied. `se = "sandwich"` requires
  raw data, with either a polychoric/tetrachoric `cor_method` (ML, ULS,
  or DWLS) or a Pearson `cor_method` (ML or ULS); it is not available
  for PAF. Under `cor_method = "fiml"`, `"information"` and `"sandwich"`
  are available for ML and ULS only and both return the corrected
  two-stage sandwich. `se = "np-boot"` requires raw data and works with
  any estimator, rotation, and correlation method. Neither
  `"information"` nor `"sandwich"` is available with the `"promax"` or
  `"simplimax"` rotations.

- **Fit indices.** The chi-square-based indices are available for ML and
  ULS (and, as scaled statistics, for `cor_method = "fiml"` and for DWLS
  with `se = "sandwich"`); PAF and DWLS otherwise report only the
  descriptive residual indices.

### Using the type presets

The `type` argument is evaluated for PAF and for all rotations (mainly
important for the varimax and promax rotations). The type-specific
settings for these functions are detailed below.

For PAF, the values of `init_comm`, `criterion`, `criterion_type`,
`max_iter`, and `abs_eigen` depend on the `type` argument.

`type = "EFAtools"` will use the following argument specification:
`init_comm = "smc", criterion = .001, criterion_type = "sum", max_iter = 300, abs_eigen = TRUE`.

`type = "psych"` will use the following argument specification:
`init_comm = "smc", criterion = .001, criterion_type = "sum", max_iter = 50, abs_eigen = FALSE`.

`type = "SPSS"` will use the following argument specification:
`init_comm = "smc", criterion = .001, criterion_type = "max_individual", max_iter = 25, abs_eigen = TRUE`.

If SMCs fail, SPSS takes "mac". However, as SPSS takes absolute
eigenvalues, this is hardly ever the case. Psych, on the other hand,
takes "unity" if SMCs fail, but uses the Moore-Penrose Psudo Inverse of
a matrix, thus, taking "unity" is only necessary if negative eigenvalues
occur afterwards in the iterative PAF procedure. The EFAtools type
setting combination was the best in terms of accuracy and number of
Heywood cases compared to all the other setting combinations tested in
simulation studies in Grieder & Steiner (2022), which is why this type
is used as a default here.

For varimax, the values of `varimax_type` and `order_type` depend on the
`type` argument.

`type = "EFAtools"` will use the following argument specification:
`varimax_type = "kaiser", order_type = "eigen"`.

`type = "psych"` will use the following argument specification:
`varimax_type = "svd", order_type = "eigen"`.

`type = "SPSS"` will use the following argument specification:
`varimax_type = "kaiser", order_type = "ss_factors"`.

For promax, the values of `p_type`, `order_type`, and `k` depend on the
`type` argument.

`type = "EFAtools"` will use the following argument specification:
`p_type = "norm", order_type = "eigen", k = 4`.

`type = "psych"` will use the following argument specification:
`p_type = "unnorm", order_type = "eigen", k = 4`.

`type = "SPSS"` will use the following argument specification:
`p_type = "norm", order_type = "ss_factors", k = 4`.

The `p_type` argument can take two values, "unnorm" and "norm". It
controls which formula is used to compute the target matrix P in the
promax rotation. "unnorm" uses the formula from Hendrickson and White
(1964), specifically: `P = abs(A^(k + 1)) / A`, where A is the
unnormalized matrix containing varimax rotated loadings. "norm" uses the
normalized varimax rotated loadings. Specifically it used the following
formula, which can be found in the SPSS 23 and SPSS 27 Algorithms
manuals:
`P = abs(A / sqrt(rowSums(A^2))) ^(k + 1) * (sqrt(rowSums(A^2)) / A)`.
As for PAF, the EFAtools type setting combination for promax was the
best compared to the other setting combinations tested in simulation
studies in Grieder & Steiner (2022). Note that all `type` presets keep
the EFAtools default Kaiser normalization (`normalize = TRUE`), whereas
[`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html) does not
normalize before its promax target rotation; set `normalize = FALSE` to
reproduce the [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html)
promax result to within the varimax convergence tolerance (the residual
difference is the convergence noise of the underlying varimax base at
`precision = 1e-5`, not an algorithmic difference).

The `varimax_type` argument can take two values, "svd", and "kaiser".
"svd" uses singular value decomposition, by calling
[`stats::varimax()`](https://rdrr.io/r/stats/varimax.html). "kaiser"
performs the varimax procedure as described in the SPSS Algorithms
manual and by Kaiser (1958). The varimax simplicity criterion monitored
for convergence is
`sum(n*colSums(lambda ^ 4) - colSums(lambda ^ 2) ^ 2) / n ^ 2`, where n
is the number of indicators, and lambda is the Kaiser-normalized rotated
loadings matrix.

For all other rotations except varimax and promax, the `type` argument
only controls the `order_type` argument with the same values as stated
above for the varimax and promax rotations. Additional arguments can
also be specified and will be passed to the rotation procedure (e.g.,
maxit to change the maximum number of iterations). As for promax, every
preset keeps the EFAtools default Kaiser normalization
(`normalize = TRUE`), whereas
[`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html) and GPArotation
do not normalize before these criterion rotations; set
`normalize = FALSE` to reproduce them.

The `type` tuning arguments have no effect on ULS and ML; `type` itself
still governs the checks applied to the correlation matrix. For ULS, no
additional arguments are needed. For ML, an additional argument
`start_method` is needed to determine the starting values for the
optimization procedure. Default for this argument is "psych" which takes
the starting values specified in
[`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html).

## See also

[`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
and
[`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
for the estimation and rotation tuning knobs.
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
for choosing `n_factors`, and
[`efa_scores()`](https://mdsteiner.github.io/EFAtools/reference/efa_scores.md),
[`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md),
[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md),
and
[`efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/efa_compare.md)
for working with the fitted solution.

Other factor analysis:
[`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md),
[`efa_group()`](https://mdsteiner.github.io/EFAtools/reference/efa_group.md),
[`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md),
[`plot.efa_group()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_group.md),
[`print.efa_group()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_group.md)

## Examples

``` r

# Principal axis factoring with oblimin rotation
mod_oblimin <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
                       rotation = "oblimin")
mod_oblimin
#> 
#> EFA performed with estimator = 'PAF' and rotation = 'oblimin'.
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>        F1     F2     F3    h2    u2
#> V1   -.035   .049   .597  .367  .633
#> V2    .010   .077   .470  .277  .723
#> V3    .069   .066   .442  .283  .717
#> V4    .110   .007   .536  .378  .622
#> V5    .164  -.005   .427  .293  .707
#> V6   -.058  -.032   .684  .399  .601
#> V7    .012   .525   .099  .357  .643
#> V8   -.005   .570   .039  .349  .651
#> V9    .047   .540   .008  .330  .670
#> V10  -.011   .659  -.058  .383  .617
#> V11   .025   .355   .232  .297  .703
#> V12   .031   .638   .001  .432  .568
#> V13   .606   .092  -.059  .394  .606
#> V14   .541  -.056   .088  .322  .678
#> V15   .554   .133  -.062  .363  .637
#> V16   .548  -.039   .091  .344  .656
#> V17   .654  -.027  -.023  .390  .610
#> V18   .549   .014   .052  .350  .650
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1     F2     F3
#> F1  1.000
#> F2   .591  1.000
#> F3   .623   .598  1.000
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F2     F3
#> SS loadings        2.218  2.083  2.005
#> Prop Tot Var        .123   .116   .111
#> Cum Prop Tot Var    .123   .239   .350
#> Prop Comm Var       .352   .330   .318
#> Cum Prop Comm Var   .352   .682  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> CAF: .50
#> SRMR: .02
#> df: 102
summary(mod_oblimin)
#> 
#> EFA performed with estimator = 'PAF' and rotation = 'oblimin'.
#> 
#> ── Model Diagnostics ───────────────────────────────────────────────────────────
#> 
#> Factors: 3
#> Variables: 18
#> N: 500
#> Rotation local optima: 1 distinct from 6 of 101 starts
#> Heywood cases: 0
#> Cross-loading items (|loading| >= .300): 0
#> Items without salient loading (|loading| >= .300): 0
#> Factors with fewer than 3 salient indicators: 0
#> Items with primary-loading gap < .200: 1
#> Largest |residual|: .069
#> Factor intercorrelations > .85: none
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>        F1     F2     F3    h2    u2
#> V1   -.035   .049   .597  .367  .633
#> V2    .010   .077   .470  .277  .723
#> V3    .069   .066   .442  .283  .717
#> V4    .110   .007   .536  .378  .622
#> V5    .164  -.005   .427  .293  .707
#> V6   -.058  -.032   .684  .399  .601
#> V7    .012   .525   .099  .357  .643
#> V8   -.005   .570   .039  .349  .651
#> V9    .047   .540   .008  .330  .670
#> V10  -.011   .659  -.058  .383  .617
#> V11   .025   .355   .232  .297  .703
#> V12   .031   .638   .001  .432  .568
#> V13   .606   .092  -.059  .394  .606
#> V14   .541  -.056   .088  .322  .678
#> V15   .554   .133  -.062  .363  .637
#> V16   .548  -.039   .091  .344  .656
#> V17   .654  -.027  -.023  .390  .610
#> V18   .549   .014   .052  .350  .650
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1     F2     F3
#> F1  1.000
#> F2   .591  1.000
#> F3   .623   .598  1.000
#> 
#> ── Structure Matrix ────────────────────────────────────────────────────────────
#> 
#>       F1    F2    F3
#> V1   .366  .385  .604
#> V2   .348  .364  .522
#> V3   .384  .371  .525
#> V4   .448  .392  .609
#> V5   .427  .347  .526
#> V6   .350  .343  .629
#> V7   .384  .591  .420
#> V8   .356  .590  .376
#> V9   .371  .573  .361
#> V10  .341  .617  .328
#> V11  .379  .508  .460
#> V12  .408  .657  .401
#> V13  .624  .415  .374
#> V14  .563  .317  .392
#> V15  .594  .423  .363
#> V16  .582  .340  .410
#> V17  .624  .346  .369
#> V18  .589  .369  .402
#> 
#> ── Simple Structure Diagnostics ────────────────────────────────────────────────
#> 
#> Items with primary-loading gap < .200:
#> • V11: F2 = .355, F3 = .232
#> 
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F2     F3
#> SS loadings        2.218  2.083  2.005
#> Prop Tot Var        .123   .116   .111
#> Cum Prop Tot Var    .123   .239   .350
#> Prop Comm Var       .352   .330   .318
#> Cum Prop Comm Var   .352   .682  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> CAF: .50
#> SRMR: .02
#> df: 102
#> 
#> ── Residual Diagnostics ────────────────────────────────────────────────────────
#> 
#> Residual cutoff: |r| > .100
#> Number of large residuals: 0
#> Largest absolute residual: .069
#> 
#> No absolute residuals > .100 occurred.
#> 
#> Inspect the residual matrix for details (e.g., with residuals()).

# ML estimation with oblimin rotation
mod_oblimin <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
                       estimator = "ML", rotation = "oblimin")
mod_oblimin
#> 
#> EFA performed with estimator = 'ML' and rotation = 'oblimin'.
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>        F1     F2     F3    h2    u2
#> V1   -.036   .043   .607  .373  .627
#> V2    .013   .087   .458  .274  .726
#> V3    .074   .074   .430  .280  .720
#> V4    .111   .007   .536  .379  .621
#> V5    .164   .005   .418  .290  .710
#> V6   -.055  -.036   .687  .402  .598
#> V7    .017   .524   .095  .355  .645
#> V8   -.003   .562   .044  .345  .655
#> V9    .044   .535   .017  .328  .672
#> V10  -.019   .661  -.051  .385  .615
#> V11   .030   .352   .230  .296  .704
#> V12   .034   .649  -.015  .437  .563
#> V13   .612   .095  -.068  .397  .603
#> V14   .540  -.053   .086  .320  .680
#> V15   .552   .137  -.065  .363  .637
#> V16   .550  -.039   .092  .345  .655
#> V17   .652  -.035  -.013  .390  .610
#> V18   .549   .012   .052  .349  .651
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1     F2     F3
#> F1  1.000
#> F2   .591  1.000
#> F3   .621   .596  1.000
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F2     F3
#> SS loadings        2.225  2.088  1.994
#> Prop Tot Var        .124   .116   .111
#> Cum Prop Tot Var    .124   .240   .350
#> Prop Comm Var       .353   .331   .316
#> Cum Prop Comm Var   .353   .684  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> χ²(102) = 123.75, p = .070
#> CFI: .99
#> TLI: .98
#> RMSEA [90% CI]: .02 [.00; .03]
#> AIC: -80.25
#> BIC: -510.14
#> ECVI: 0.52
#> CAF: .50
#> SRMR: .03
summary(mod_oblimin)
#> 
#> EFA performed with estimator = 'ML' and rotation = 'oblimin'.
#> 
#> ── Model Diagnostics ───────────────────────────────────────────────────────────
#> 
#> Factors: 3
#> Variables: 18
#> N: 500
#> Rotation local optima: 1 distinct from 6 of 101 starts
#> Heywood cases: 0
#> Cross-loading items (|loading| >= .300): 0
#> Items without salient loading (|loading| >= .300): 0
#> Factors with fewer than 3 salient indicators: 0
#> Items with primary-loading gap < .200: 1
#> Largest |residual|: .069
#> Factor intercorrelations > .85: none
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>        F1     F2     F3    h2    u2
#> V1   -.036   .043   .607  .373  .627
#> V2    .013   .087   .458  .274  .726
#> V3    .074   .074   .430  .280  .720
#> V4    .111   .007   .536  .379  .621
#> V5    .164   .005   .418  .290  .710
#> V6   -.055  -.036   .687  .402  .598
#> V7    .017   .524   .095  .355  .645
#> V8   -.003   .562   .044  .345  .655
#> V9    .044   .535   .017  .328  .672
#> V10  -.019   .661  -.051  .385  .615
#> V11   .030   .352   .230  .296  .704
#> V12   .034   .649  -.015  .437  .563
#> V13   .612   .095  -.068  .397  .603
#> V14   .540  -.053   .086  .320  .680
#> V15   .552   .137  -.065  .363  .637
#> V16   .550  -.039   .092  .345  .655
#> V17   .652  -.035  -.013  .390  .610
#> V18   .549   .012   .052  .349  .651
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1     F2     F3
#> F1  1.000
#> F2   .591  1.000
#> F3   .621   .596  1.000
#> 
#> ── Structure Matrix ────────────────────────────────────────────────────────────
#> 
#>       F1    F2    F3
#> V1   .366  .384  .610
#> V2   .349  .367  .518
#> V3   .385  .374  .520
#> V4   .448  .392  .609
#> V5   .427  .351  .523
#> V6   .350  .341  .631
#> V7   .385  .590  .418
#> V8   .357  .587  .377
#> V9   .371  .571  .363
#> V10  .339  .619  .331
#> V11  .381  .507  .459
#> V12  .409  .661  .393
#> V13  .626  .416  .369
#> V14  .562  .317  .390
#> V15  .593  .425  .360
#> V16  .584  .341  .410
#> V17  .623  .343  .371
#> V18  .589  .368  .401
#> 
#> ── Simple Structure Diagnostics ────────────────────────────────────────────────
#> 
#> Items with primary-loading gap < .200:
#> • V11: F2 = .352, F3 = .230
#> 
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F2     F3
#> SS loadings        2.225  2.088  1.994
#> Prop Tot Var        .124   .116   .111
#> Cum Prop Tot Var    .124   .240   .350
#> Prop Comm Var       .353   .331   .316
#> Cum Prop Comm Var   .353   .684  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> χ²(102) = 123.75, p = .070
#> CFI: .99
#> TLI: .98
#> RMSEA [90% CI]: .02 [.00; .03]
#> AIC: -80.25
#> BIC: -510.14
#> ECVI: 0.52
#> CAF: .50
#> SRMR: .03
#> 
#> ── Residual Diagnostics ────────────────────────────────────────────────────────
#> 
#> Residual cutoff: |r| > .100
#> Number of large residuals: 0
#> Largest absolute residual: .069
#> 
#> No absolute residuals > .100 occurred.
#> 
#> Inspect the residual matrix for details (e.g., with residuals()).

# Tuning knobs are supplied through the control objects. Here the SPSS preset is
# used for the estimation and rotation, with the maximum PAF iterations raised.
mod_spss <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
                    rotation = "promax",
                    estimate_control = estimate_control(type = "SPSS", max_iter = 500),
                    rotate_control = rotate_control(type = "SPSS"))
#> Warning: An argument was set together with `type` = "SPSS"; the supplied value is used
#> and may differ from the "SPSS" preset:
#> • max_iter = 500
mod_spss
#> 
#> EFA performed with estimator = 'PAF' and rotation = 'promax'.
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>        F1     F2     F3    h2    u2
#> V1   -.048   .035   .613  .367  .633
#> V2   -.001   .067   .482  .277  .723
#> V3    .060   .056   .453  .283  .717
#> V4    .101  -.009   .551  .378  .622
#> V5    .157  -.018   .438  .293  .707
#> V6   -.072  -.049   .704  .399  .601
#> V7    .001   .533   .093  .357  .643
#> V8   -.016   .581   .030  .349  .651
#> V9    .038   .550  -.001  .330  .670
#> V10  -.021   .674  -.071  .383  .617
#> V11   .015   .356   .232  .297  .703
#> V12   .020   .651  -.010  .432  .568
#> V13   .614   .086  -.067  .394  .606
#> V14   .548  -.068   .088  .322  .678
#> V15   .561   .128  -.070  .363  .637
#> V16   .555  -.050   .091  .344  .656
#> V17   .664  -.037  -.027  .390  .610
#> V18   .555   .004   .050  .350  .650
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1     F2     F3
#> F1  1.000
#> F2   .617  1.000
#> F3   .648   .632  1.000
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F2     F3
#> SS loadings        2.199  2.074  2.034
#> Prop Tot Var        .122   .115   .113
#> Cum Prop Tot Var    .122   .237   .350
#> Prop Comm Var       .349   .329   .323
#> Cum Prop Comm Var   .349   .677  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> CAF: .50
#> SRMR: .02
#> df: 102

# Analytic (expected-information) standard errors for the above
ML_info <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
                   estimator = "ML", rotation = "oblimin", se = "information")
ML_info
#> 
#> EFA performed with estimator = 'ML' and rotation = 'oblimin'.
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>        F1     F2     F3    h2    u2
#> V1   -.036   .043   .607  .373  .627
#> V2    .013   .087   .458  .274  .726
#> V3    .074   .074   .430  .280  .720
#> V4    .111   .007   .536  .379  .621
#> V5    .164   .005   .418  .290  .710
#> V6   -.055  -.036   .687  .402  .598
#> V7    .017   .524   .095  .355  .645
#> V8   -.003   .562   .044  .345  .655
#> V9    .044   .535   .017  .328  .672
#> V10  -.019   .661  -.051  .385  .615
#> V11   .030   .352   .230  .296  .704
#> V12   .034   .649  -.015  .437  .563
#> V13   .612   .095  -.068  .397  .603
#> V14   .540  -.053   .086  .320  .680
#> V15   .552   .137  -.065  .363  .637
#> V16   .550  -.039   .092  .345  .655
#> V17   .652  -.035  -.013  .390  .610
#> V18   .549   .012   .052  .349  .651
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1     F2     F3
#> F1  1.000
#> F2   .591  1.000
#> F3   .621   .596  1.000
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F2     F3
#> SS loadings        2.225  2.088  1.994
#> Prop Tot Var        .124   .116   .111
#> Cum Prop Tot Var    .124   .240   .350
#> Prop Comm Var       .353   .331   .316
#> Cum Prop Comm Var   .353   .684  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> χ²(102) = 123.75, p = .070
#> CFI: .99
#> TLI: .98
#> RMSEA [90% CI]: .02 [.00; .03]
#> AIC: -80.25
#> BIC: -510.14
#> ECVI: 0.52
#> CAF: .50
#> SRMR: .03
summary(ML_info)
#> 
#> EFA performed with estimator = 'ML' and rotation = 'oblimin'.
#> 
#> ── Model Diagnostics ───────────────────────────────────────────────────────────
#> 
#> Factors: 3
#> Variables: 18
#> N: 500
#> Rotation local optima: 1 distinct from 6 of 101 starts
#> Heywood cases: 0
#> Cross-loading items (|loading| >= .300): 0
#> Items without salient loading (|loading| >= .300): 0
#> Factors with fewer than 3 salient indicators: 0
#> Items with primary-loading gap < .200: 1
#> Largest |residual|: .069
#> Factor intercorrelations > .85: none
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>        F1     F2     F3    h2    u2
#> V1   -.036   .043   .607  .373  .627
#> V2    .013   .087   .458  .274  .726
#> V3    .074   .074   .430  .280  .720
#> V4    .111   .007   .536  .379  .621
#> V5    .164   .005   .418  .290  .710
#> V6   -.055  -.036   .687  .402  .598
#> V7    .017   .524   .095  .355  .645
#> V8   -.003   .562   .044  .345  .655
#> V9    .044   .535   .017  .328  .672
#> V10  -.019   .661  -.051  .385  .615
#> V11   .030   .352   .230  .296  .704
#> V12   .034   .649  -.015  .437  .563
#> V13   .612   .095  -.068  .397  .603
#> V14   .540  -.053   .086  .320  .680
#> V15   .552   .137  -.065  .363  .637
#> V16   .550  -.039   .092  .345  .655
#> V17   .652  -.035  -.013  .390  .610
#> V18   .549   .012   .052  .349  .651
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── 95% Wald CIs for salient rotated loadings ───────────────────────────────────
#> 
#> Variable  Factor  est    lower  upper
#> V13       F1       .612   .488   .736
#> V14       F1       .540   .411   .669
#> V15       F1       .552   .423   .681
#> V16       F1       .550   .421   .679
#> V17       F1       .652   .535   .769
#> V18       F1       .549   .420   .679
#> V7        F2       .524   .394   .653
#> V8        F2       .562   .437   .687
#> V9        F2       .535   .407   .662
#> V10       F2       .661   .548   .773
#> V11       F2       .352   .211   .493
#> V12       F2       .649   .529   .770
#> V1        F3       .607   .474   .739
#> V2        F3       .458   .313   .604
#> V3        F3       .430   .282   .578
#> V4        F3       .536   .395   .677
#> V5        F3       .418   .269   .567
#> V6        F3       .687   .570   .805
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1     F2     F3
#> F1  1.000
#> F2   .591  1.000
#> F3   .621   .596  1.000
#> 
#> ── 95% Wald CIs for factor intercorrelations ───────────────────────────────────
#> 
#> Factors   est    lower  upper
#> F1 ~~ F2   .591   .499   .683
#> F1 ~~ F3   .621   .531   .712
#> F2 ~~ F3   .596   .503   .690
#> 
#> ── Structure Matrix ────────────────────────────────────────────────────────────
#> 
#>       F1    F2    F3
#> V1   .366  .384  .610
#> V2   .349  .367  .518
#> V3   .385  .374  .520
#> V4   .448  .392  .609
#> V5   .427  .351  .523
#> V6   .350  .341  .631
#> V7   .385  .590  .418
#> V8   .357  .587  .377
#> V9   .371  .571  .363
#> V10  .339  .619  .331
#> V11  .381  .507  .459
#> V12  .409  .661  .393
#> V13  .626  .416  .369
#> V14  .562  .317  .390
#> V15  .593  .425  .360
#> V16  .584  .341  .410
#> V17  .623  .343  .371
#> V18  .589  .368  .401
#> 
#> ── Simple Structure Diagnostics ────────────────────────────────────────────────
#> 
#> Items with primary-loading gap < .200:
#> • V11: F2 = .352, F3 = .230
#> 
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F2     F3
#> SS loadings        2.225  2.088  1.994
#> Prop Tot Var        .124   .116   .111
#> Cum Prop Tot Var    .124   .240   .350
#> Prop Comm Var       .353   .331   .316
#> Cum Prop Comm Var   .353   .684  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> χ²(102) = 123.75, p = .070
#> CFI: .99
#> TLI: .98
#> RMSEA [90% CI]: .02 [.00; .03]
#> AIC: -80.25
#> BIC: -510.14
#> ECVI: 0.52
#> CAF: .50
#> SRMR: .03
#> 
#> ── Residual Diagnostics ────────────────────────────────────────────────────────
#> 
#> Residual cutoff: |r| > .100
#> Number of large residuals: 0
#> Largest absolute residual: .069
#> 
#> No absolute residuals > .100 occurred.
#> 
#> Inspect the residual matrix for details (e.g., with residuals()).

# \donttest{
# Robust (sandwich) standard errors and a scaled chi-square for ordinal raw data.
# These need a polychoric/tetrachoric correlation method and estimator ML, ULS, or DWLS.
DWLS_rob <- efa_fit(DOSPERT_raw, n_factors = 6, cor_method = "poly",
                    estimator = "DWLS", rotation = "oblimin", se = "sandwich")
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
#> Warning: Some response-category combinations are empty despite a non-negligible expected
#> count.
#> ℹ The polychoric asymptotic covariance (and any DWLS weights or robust standard
#>   errors derived from it) can be unreliable for such structurally sparse cells;
#>   interpret them with caution.
#> Warning: Analytic standard errors could not be computed for all parameters.
#> ℹ The factor solution's rotational orientation is only weakly determined (two
#>   canonical variances nearly coincide), so the unrotated loadings have no
#>   well-defined standard error. The rotated loadings and communalities are
#>   gauge-invariant and unaffected.
DWLS_rob
#> 
#> EFA performed with estimator = 'DWLS' and rotation = 'oblimin'.
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>           F1     F2     F3     F4     F5     F6    h2    u2
#> ethR_1   .022   .508   .005   .076   .071   .121  .388  .612
#> ethR_2   .019   .565   .082   .131   .046   .014  .435  .565
#> ethR_3   .074   .739  -.227   .017   .026   .064  .697  .303
#> ethR_4  -.072   .617  -.052  -.062   .034   .041  .371  .629
#> ethR_5   .143   .544  -.118   .024   .039  -.006  .408  .592
#> ethR_6   .000   .685   .039   .037  -.016  -.066  .467  .533
#> finR_1   .043   .053  -.021   .804   .051   .134  .828  .172
#> finR_2   .011  -.025   .045  -.011  -.099   .699  .486  .514
#> finR_3   .055   .065   .019   .820   .033   .115  .840  .160
#> finR_4  -.017   .060  -.074   .122   .012   .779  .681  .319
#> finR_5   .038   .078   .015   .834   .013   .131  .866  .134
#> finR_6  -.004   .006   .055   .132   .054   .741  .663  .337
#> heaR_1   .091   .373   .118   .127   .179  -.034  .344  .656
#> heaR_2   .027   .337   .165   .087   .267  -.060  .336  .664
#> heaR_3  -.040   .082  -.076   .072   .676  -.001  .530  .470
#> heaR_4   .021   .053  -.119   .136   .705  -.026  .626  .374
#> heaR_5  -.030   .141   .075  -.040   .463  -.004  .281  .719
#> heaR_6   .122   .234   .122   .038   .397   .038  .427  .573
#> recR_1   .364  -.152   .198  -.094   .308   .069  .368  .632
#> recR_2   .489  -.054  -.178   .102   .354   .118  .596  .404
#> recR_3   .549  -.098  -.041   .019   .328   .113  .577  .423
#> recR_4   .997   .088   .016   .020  -.180  -.036  .904  .096
#> recR_5   .914   .133  -.030   .056  -.105  -.043  .827  .173
#> recR_6   .621  -.007   .052   .051   .074   .120  .534  .466
#> socR_1  -.067  -.058   .693  -.073  -.003  -.020  .507  .493
#> socR_2   .002   .023   .805   .173   .000  -.040  .623  .377
#> socR_3  -.042  -.110   .648  -.070  -.022   .046  .475  .525
#> socR_4   .040  -.078   .729   .159  -.009  -.068  .515  .485
#> socR_5   .121   .115   .436  -.185  -.001   .121  .299  .701
#> socR_6  -.006   .051   .592  -.147   .010   .145  .434  .566
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1     F2     F3     F4     F5     F6
#> F1  1.000
#> F2   .231  1.000
#> F3   .134  -.074  1.000
#> F4   .236   .397  -.147  1.000
#> F5   .404   .469   .112   .312  1.000
#> F6   .363   .151   .195   .290   .181  1.000
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F2     F3     F4     F5     F6
#> SS loadings        3.280  3.153  2.887  2.708  2.250  2.053
#> Prop Tot Var        .109   .105   .096   .090   .075   .068
#> Cum Prop Tot Var    .109   .214   .311   .401   .476   .544
#> Prop Comm Var       .201   .193   .177   .166   .138   .126
#> Cum Prop Comm Var   .201   .394   .571   .737   .874  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> scaled χ²(270) = 3521.43, p < .001
#> CFI: .96
#> TLI: .93
#> RMSEA [90% CI]: .06 [.06; .06]
#> AIC: NA
#> BIC: NA
#> CAF: .41
#> SRMR: .03
summary(DWLS_rob)
#> 
#> EFA performed with estimator = 'DWLS' and rotation = 'oblimin'.
#> 
#> ── Model Diagnostics ───────────────────────────────────────────────────────────
#> 
#> Factors: 6
#> Variables: 30
#> N: 3123
#> Rotation local optima: 1 distinct from 6 of 101 starts
#> Heywood cases: 0
#> Cross-loading items (|loading| >= .300): 3
#> Items without salient loading (|loading| >= .300): 0
#> Factors with fewer than 3 salient indicators: 0
#> Items with primary-loading gap < .200: 5
#> Largest |residual|: .153
#> Factor intercorrelations > .85: none
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>           F1     F2     F3     F4     F5     F6    h2    u2
#> ethR_1   .022   .508   .005   .076   .071   .121  .388  .612
#> ethR_2   .019   .565   .082   .131   .046   .014  .435  .565
#> ethR_3   .074   .739  -.227   .017   .026   .064  .697  .303
#> ethR_4  -.072   .617  -.052  -.062   .034   .041  .371  .629
#> ethR_5   .143   .544  -.118   .024   .039  -.006  .408  .592
#> ethR_6   .000   .685   .039   .037  -.016  -.066  .467  .533
#> finR_1   .043   .053  -.021   .804   .051   .134  .828  .172
#> finR_2   .011  -.025   .045  -.011  -.099   .699  .486  .514
#> finR_3   .055   .065   .019   .820   .033   .115  .840  .160
#> finR_4  -.017   .060  -.074   .122   .012   .779  .681  .319
#> finR_5   .038   .078   .015   .834   .013   .131  .866  .134
#> finR_6  -.004   .006   .055   .132   .054   .741  .663  .337
#> heaR_1   .091   .373   .118   .127   .179  -.034  .344  .656
#> heaR_2   .027   .337   .165   .087   .267  -.060  .336  .664
#> heaR_3  -.040   .082  -.076   .072   .676  -.001  .530  .470
#> heaR_4   .021   .053  -.119   .136   .705  -.026  .626  .374
#> heaR_5  -.030   .141   .075  -.040   .463  -.004  .281  .719
#> heaR_6   .122   .234   .122   .038   .397   .038  .427  .573
#> recR_1   .364  -.152   .198  -.094   .308   .069  .368  .632
#> recR_2   .489  -.054  -.178   .102   .354   .118  .596  .404
#> recR_3   .549  -.098  -.041   .019   .328   .113  .577  .423
#> recR_4   .997   .088   .016   .020  -.180  -.036  .904  .096
#> recR_5   .914   .133  -.030   .056  -.105  -.043  .827  .173
#> recR_6   .621  -.007   .052   .051   .074   .120  .534  .466
#> socR_1  -.067  -.058   .693  -.073  -.003  -.020  .507  .493
#> socR_2   .002   .023   .805   .173   .000  -.040  .623  .377
#> socR_3  -.042  -.110   .648  -.070  -.022   .046  .475  .525
#> socR_4   .040  -.078   .729   .159  -.009  -.068  .515  .485
#> socR_5   .121   .115   .436  -.185  -.001   .121  .299  .701
#> socR_6  -.006   .051   .592  -.147   .010   .145  .434  .566
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── 95% Wald CIs for salient rotated loadings ───────────────────────────────────
#> 
#> Variable  Factor  est    lower  upper
#> recR_1    F1       .364   .315   .414
#> recR_2    F1       .489   .445   .534
#> recR_3    F1       .549   .506   .592
#> recR_4    F1       .997   .977  1.017
#> recR_5    F1       .914   .894   .934
#> recR_6    F1       .621   .591   .650
#> ethR_1    F2       .508   .467   .550
#> ethR_2    F2       .565   .519   .610
#> ethR_3    F2       .739   .700   .778
#> ethR_4    F2       .617   .580   .654
#> ethR_5    F2       .544   .500   .588
#> ethR_6    F2       .685   .649   .721
#> heaR_1    F2       .373   .323   .423
#> heaR_2    F2       .337   .284   .389
#> socR_1    F3       .693   .667   .719
#> socR_2    F3       .805   .783   .827
#> socR_3    F3       .648   .621   .675
#> socR_4    F3       .729   .704   .754
#> socR_5    F3       .436   .400   .473
#> socR_6    F3       .592   .561   .622
#> finR_1    F4       .804   .773   .835
#> finR_3    F4       .820   .790   .850
#> finR_5    F4       .834   .804   .864
#> heaR_3    F5       .676   .633   .719
#> heaR_4    F5       .705   .658   .752
#> heaR_5    F5       .463   .414   .513
#> heaR_6    F5       .397   .347   .448
#> recR_1    F5       .308   .254   .362
#> recR_2    F5       .354   .306   .402
#> recR_3    F5       .328   .278   .378
#> finR_2    F6       .699   .671   .727
#> finR_4    F6       .779   .751   .806
#> finR_6    F6       .741   .714   .768
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1     F2     F3     F4     F5     F6
#> F1  1.000
#> F2   .231  1.000
#> F3   .134  -.074  1.000
#> F4   .236   .397  -.147  1.000
#> F5   .404   .469   .112   .312  1.000
#> F6   .363   .151   .195   .290   .181  1.000
#> 
#> ── 95% Wald CIs for factor intercorrelations ───────────────────────────────────
#> 
#> Factors   est    lower  upper
#> F1 ~~ F2   .231   .195   .267
#> F1 ~~ F3   .134   .106   .162
#> F1 ~~ F4   .236   .196   .277
#> F1 ~~ F5   .404   .372   .436
#> F1 ~~ F6   .363   .329   .396
#> F2 ~~ F3  -.074  -.104  -.043
#> F2 ~~ F4   .397   .358   .435
#> F2 ~~ F5   .469   .440   .498
#> F2 ~~ F6   .151   .114   .188
#> F3 ~~ F4  -.147  -.191  -.103
#> F3 ~~ F5   .112   .082   .142
#> F3 ~~ F6   .195   .162   .229
#> F4 ~~ F5   .312   .269   .355
#> F4 ~~ F6   .290   .247   .334
#> F5 ~~ F6   .181   .142   .220
#> 
#> ── Structure Matrix ────────────────────────────────────────────────────────────
#> 
#>           F1     F2     F3     F4     F5    F6
#> ethR_1   .231   .595  -.009   .339   .364  .242
#> ethR_2   .216   .639   .032   .366   .371  .169
#> ethR_3   .252   .801  -.258   .388   .394  .168
#> ethR_4   .077   .601  -.086   .196   .276  .086
#> ethR_5   .272   .612  -.139   .301   .345  .119
#> ethR_6   .142   .680  -.032   .279   .310  .053
#> finR_1   .312   .427  -.106   .893   .366  .396
#> finR_2   .222   .030   .175   .148   .022  .687
#> finR_3   .322   .434  -.073   .900   .365  .393
#> finR_4   .303   .233   .055   .383   .204  .805
#> finR_5   .308   .443  -.081   .914   .351  .404
#> finR_6   .326   .191   .185   .358   .237  .799
#> heaR_1   .283   .515   .097   .325   .438  .148
#> heaR_2   .234   .482   .149   .269   .471  .107
#> heaR_3   .258   .424  -.022   .317   .712  .125
#> heaR_4   .325   .447  -.066   .392   .763  .133
#> heaR_5   .189   .330   .118   .141   .513  .094
#> heaR_6   .376   .460   .167   .277   .589  .224
#> recR_1   .483   .035   .320   .019   .389  .245
#> recR_2   .663   .297  -.060   .368   .560  .346
#> recR_3   .699   .211   .096   .252   .526  .355
#> recR_4   .939   .235   .113   .222   .266  .316
#> recR_5   .896   .313   .054   .284   .333  .300
#> recR_6   .711   .205   .160   .245   .365  .383
#> socR_1  -.013  -.158   .695  -.221  -.006  .061
#> socR_2   .141   .027   .770   .052   .149  .171
#> socR_3   .011  -.199   .667  -.213  -.032  .117
#> socR_4   .129  -.073   .702   .008   .090  .122
#> socR_5   .206   .055   .495  -.140   .115  .213
#> socR_6   .108  -.025   .638  -.169   .079  .226
#> 
#> ── Simple Structure Diagnostics ────────────────────────────────────────────────
#> 
#> Items with cross-loadings, |loading| >= .300 on multiple factors:
#> • recR_1: F1 = .364, F5 = .308
#> • recR_2: F1 = .489, F5 = .354
#> • recR_3: F1 = .549, F5 = .328
#> 
#> Items with primary-loading gap < .200:
#> • heaR_1: F2 = .373, F5 = .179
#> • heaR_2: F2 = .337, F5 = .267
#> • heaR_6: F5 = .397, F2 = .234
#> • recR_1: F1 = .364, F5 = .308
#> • recR_2: F1 = .489, F5 = .354
#> 
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F2     F3     F4     F5     F6
#> SS loadings        3.280  3.153  2.887  2.708  2.250  2.053
#> Prop Tot Var        .109   .105   .096   .090   .075   .068
#> Cum Prop Tot Var    .109   .214   .311   .401   .476   .544
#> Prop Comm Var       .201   .193   .177   .166   .138   .126
#> Cum Prop Comm Var   .201   .394   .571   .737   .874  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> scaled χ²(270) = 3521.43, p < .001
#> CFI: .96
#> TLI: .93
#> RMSEA [90% CI]: .06 [.06; .06]
#> AIC: NA
#> BIC: NA
#> CAF: .41
#> SRMR: .03
#> 
#> ── Residual Diagnostics ────────────────────────────────────────────────────────
#> 
#> Residual cutoff: |r| > .100
#> Number of large residuals: 5
#> Largest absolute residual: .153
#> 
#> Largest residuals:
#> • socR_5 ~~ socR_6: .153
#> • heaR_1 ~~ heaR_2: .130
#> • heaR_1 ~~ heaR_4: -.112
#> • heaR_3 ~~ recR_1: -.110
#> • heaR_4 ~~ recR_1: -.100
#> 
#> Inspect the residual matrix for details (e.g., with residuals()).

# The same robust SEs and scaled chi-square for continuous data: a Pearson
# correlation with estimator ML or ULS (the fourth-moment ADF covariance).
ML_rob <- efa_fit(GRiPS_raw, n_factors = 1, cor_method = "pearson",
                  estimator = "ML", rotation = "none", se = "sandwich")
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
ML_rob
#> 
#> EFA performed with estimator = 'ML' and rotation = 'none'.
#> 
#> ── Unrotated Loadings ──────────────────────────────────────────────────────────
#> 
#>             F1    h2    u2
#> fun        .796  .634  .366
#> friends    .851  .725  .275
#> enjoy      .872  .760  .240
#> hurt       .767  .588  .412
#> part       .817  .667  .333
#> commonly   .832  .692  .308
#> chances    .788  .621  .379
#> attracted  .845  .715  .285
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                 F1
#> SS loadings   5.401
#> Prop Tot Var   .675
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> scaled χ²(20) = 33.64, p = .029
#> CFI: 1.00
#> TLI: 1.00
#> RMSEA [90% CI]: .03 [.01; .05]
#> AIC: NA
#> BIC: NA
#> CAF: .50
#> SRMR: .01
summary(ML_rob)
#> 
#> EFA performed with estimator = 'ML' and rotation = 'none'.
#> 
#> ── Model Diagnostics ───────────────────────────────────────────────────────────
#> 
#> Factors: 1
#> Variables: 8
#> N: 810
#> Heywood cases: 0
#> Cross-loading items (|loading| >= .300): 0
#> Items without salient loading (|loading| >= .300): 0
#> Factors with fewer than 3 salient indicators: 0
#> Items with primary-loading gap < .200: 0
#> Largest |residual|: .030
#> 
#> ── Unrotated Loadings ──────────────────────────────────────────────────────────
#> 
#>             F1    h2    u2
#> fun        .796  .634  .366
#> friends    .851  .725  .275
#> enjoy      .872  .760  .240
#> hurt       .767  .588  .412
#> part       .817  .667  .333
#> commonly   .832  .692  .308
#> chances    .788  .621  .379
#> attracted  .845  .715  .285
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── 95% Wald CIs for salient unrotated loadings ─────────────────────────────────
#> 
#> Variable   Factor  est    lower  upper
#> fun        F1       .796   .761   .831
#> friends    F1       .851   .825   .878
#> enjoy      F1       .872   .851   .893
#> hurt       F1       .767   .726   .808
#> part       F1       .817   .782   .851
#> commonly   F1       .832   .800   .864
#> chances    F1       .788   .745   .831
#> attracted  F1       .845   .814   .877
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                 F1
#> SS loadings   5.401
#> Prop Tot Var   .675
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> scaled χ²(20) = 33.64, p = .029
#> CFI: 1.00
#> TLI: 1.00
#> RMSEA [90% CI]: .03 [.01; .05]
#> AIC: NA
#> BIC: NA
#> CAF: .50
#> SRMR: .01
#> 
#> Note: Wald CIs from the robust (Godambe) sandwich covariance.
#> 
#> ── Residual Diagnostics ────────────────────────────────────────────────────────
#> 
#> Residual cutoff: |r| > .100
#> Number of large residuals: 0
#> Largest absolute residual: .030
#> 
#> No absolute residuals > .100 occurred.
#> 
#> Inspect the residual matrix for details (e.g., with residuals()).
# }

# \donttest{
# Two-stage FIML correlations from raw data with missing values: the saturated
# multivariate-normal moments are EM-estimated (assuming the data are missing at
# random) and the standardized covariance is analysed.
x_miss <- GRiPS_raw
x_miss[cbind(1:20, 1)] <- NA
efa_fiml <- efa_fit(x_miss, n_factors = 1, estimator = "ML", cor_method = "fiml")
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
efa_fiml
#> 
#> EFA performed with estimator = 'ML' and rotation = 'none'.
#> Correlations: FIML (two-stage, missing data)
#> 
#> ── Unrotated Loadings ──────────────────────────────────────────────────────────
#> 
#>             F1    h2    u2
#> fun        .795  .632  .368
#> friends    .852  .725  .275
#> enjoy      .871  .759  .241
#> hurt       .767  .588  .412
#> part       .817  .667  .333
#> commonly   .832  .692  .308
#> chances    .788  .620  .380
#> attracted  .845  .715  .285
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                 F1
#> SS loadings   5.400
#> Prop Tot Var   .675
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> scaled χ²(20) = 56.85, p < .001
#> CFI: 1.00
#> TLI: 1.00
#> RMSEA [90% CI]: .05 [.03; .06]
#> AIC: NA
#> BIC: NA
#> CAF: .50
#> SRMR: .01
# }

if (FALSE) { # \dontrun{
# Bootstrap standard errors from raw data, reproducible via a fixed seed and run
# in parallel across replicates.
future::plan(future::multisession, workers = 2)
efa_boot <- efa_fit(GRiPS_raw, n_factors = 1, estimator = "PAF", rotation = "none",
                    se = "np-boot", b_boot = 1000, seed = 42)
future::plan(future::sequential)
} # }
```
