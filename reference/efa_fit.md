# Exploratory factor analysis (EFA)

This function does an EFA with `PAF`, `ML`, `ULS`/`MINRES`, or `DWLS`
estimation, by default followed by an oblimin rotation. Estimation and
rotation are controlled through the control objects built by
[`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
and
[`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md);
each accepts a `type` ("EFAtools", "SPSS", "psych", or "none") that
fills in its remaining settings.

## Usage

``` r
efa_fit(
  x,
  n_factors,
  N = NA,
  estimator = c("PAF", "ML", "ULS", "MINRES", "DWLS"),
  rotation = c("oblimin", "none", "varimax", "equamax", "quartimax", "geominT",
    "bentlerT", "bifactorT", "promax", "quartimin", "simplimax", "bentlerQ", "geominQ",
    "bifactorQ"),
  se = c("none", "information", "sandwich", "np-boot"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra", "fiml"),
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  estimate_control = NULL,
  rotate_control = NULL,
  b_boot = 1000,
  ci = 0.95,
  seed = NULL
)
```

## Source

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
Clogg (Eds.), Latent variables analysis: Applications for developmental
research (pp. 399–419). Sage.

Asparouhov, T., & Muthén, B. (2010). Simple second order chi-square
correction. Mplus Technical Appendix.

Muthén, B., du Toit, S. H. C., & Spisic, D. (1997). Robust inference
using weighted least squares and quadratic estimating equations in
latent variable modeling with categorical and continuous outcomes.
Unpublished manuscript.

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

Bartlett, M. S. (1951). The effect of standardization on a Chi-square
approximation in factor analysis. Biometrika, 38, 337–344.

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

Jöreskog, K. G., & van Thillo, M. (1971). New rapid algorithms for
factor analysis by unweighted least squares, generalized least squares
and maximum likelihood (Research Memorandum RM-71-5). Educational
Testing Service.

Marsh, H. W. (1998). Pairwise deletion for missing data in structural
equation models: Nonpositive definite matrices, parameter estimates,
goodness of fit, and adjusted sample sizes. Structural Equation
Modeling, 5, 22–36. doi: 10.1080/10705519809540087

Yuan, K.-H., Hayashi, K., & Yanagihara, H. (2007). A class of population
covariance matrices in the bootstrap approach to covariance structure
analysis. Multivariate Behavioral Research, 42, 261–281. doi:
10.1080/00273170701360662

Zhang, X., & Savalei, V. (2016). Bootstrapping confidence intervals for
fit indexes in structural equation modeling. Structural Equation
Modeling, 23, 392–408. doi: 10.1080/10705511.2015.1118692

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
  correlation matrix is used. Without `N`, the chi-square-based fit
  indices are not available, with a warning; the same applies to a
  positive `N` that is very small relative to the number of variables.
  With raw data, `N` is found from the data:

  - the number of rows under the default pairwise deletion (`use`), with
    a warning if values are missing, because each correlation then rests
    on fewer cases, so the chi-square is too large and the standard
    errors are too small (Marsh, 1998);

  - the number of complete cases under listwise deletion, and whenever
    an asymptotic covariance is needed (`estimator = "DWLS"` or
    `se = "sandwich"`, except under `cor_method = "fiml"`);

  - the number of rows with at least one observed value under
    `cor_method = "fiml"`.

- estimator:

  character. The estimator used to fit the EFA: "PAF" (principal axis
  factoring), "ML" (maximum likelihood), "ULS" (unweighted least
  squares; "MINRES" is an accepted alias returning identical results),
  or "DWLS" (diagonally weighted least squares, for ordinal data). See
  the *Estimators* section in Details for their properties and data
  requirements. Lower-case versions (e.g., "paf") are also accepted.

- rotation:

  character. An oblique rotation ("oblimin", the default, "promax",
  "quartimin", "simplimax", "bentlerQ", "geominQ", or "bifactorQ"), an
  orthogonal rotation ("varimax", "equamax", "quartimax", "geominT",
  "bentlerT", or "bifactorT"), or no rotation ("none"). A single factor
  cannot be rotated, so with `n_factors = 1` the default is "none". See
  the *Rotations* section in Details for their properties and known
  issues.

- se:

  character. Whether and how to compute standard errors (and matching
  confidence intervals): "none" (default), "information" (analytic
  standard errors from the expected Fisher information of the ML
  solution), "sandwich" (robust standard errors from raw data, which
  stay reliable under non-normality or a misspecified estimator weight),
  or "np-boot" (non-parametric bootstrap). The methods differ in their
  assumptions, their data requirements, and which estimator, rotation,
  and `cor_method` combinations they support; see the *Standard errors*
  section in Details.

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
  given as input. Default is "pairwise.complete.obs". It is ignored when
  `cor_method = "fiml"` (which handles the missingness itself, so every
  case contributes), and it is overridden to listwise deletion whenever
  an asymptotic covariance is required (the "DWLS" estimator, or
  `se = "sandwich"`), because the covariance must describe the same
  cases as the correlation matrix.

- estimate_control:

  a control object from
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  bundling the estimation control arguments: the `type` preset; the
  principal-axis-factoring iteration settings `init_comm`, `criterion`,
  `criterion_type`, `max_iter`, and `abs_eigen`; and the
  maximum-likelihood `start_method`. Defaults to
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md),
  which uses the `"EFAtools"` type. See
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md).

- rotate_control:

  a control object from
  [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  bundling the rotation control arguments: the `type` preset; Kaiser
  `normalize`; the convergence `precision`; the factor `order_type`; the
  varimax/promax settings `varimax_type` and `p_type`; the
  simplimax/promax `k`; `random_starts`; and the rotation-engine
  arguments `gam` (oblimin), `delta` (geomin), and `maxit`. Defaults to
  [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md),
  which uses the `"EFAtools"` type. The estimation and rotation presets
  are independent. See
  [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  for details.

- b_boot:

  numeric. The number of bootstrap samples to draw. Default is 1000.
  Must be at least 2, the smallest number from which a standard error is
  defined. Under `cor_method = "fiml"` each bootstrap sample re-runs the
  EM moment estimation, so a smaller value may be advisable.

- ci:

  numeric. The level of the confidence intervals: the percentile
  intervals from the bootstrap samples under `se = "np-boot"`, and the
  analytic Wald intervals under `se = "information"` and
  `se = "sandwich"`, the corrected two-stage intervals of
  `cor_method = "fiml"` included. Must be greater than 0 and smaller
  than 1. Default is .95 for 95% CIs.

- seed:

  numeric. An optional seed for the random-number generator.

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
  the unrotated loadings. Its rows are `"SS loadings"` and
  `"Prop Tot Var"`, followed by `"Cum Prop Tot Var"`, `"Prop Comm Var"`,
  and `"Cum Prop Comm Var"`; those last three are omitted for a
  single-factor solution, where they would only repeat the two above
  them. So code that indexes a row by name (for example
  `["Prop Comm Var", ]`) must allow for the two-row form at
  `n_factors = 1`.

- fit_indices:

  A named list of fit indices computed from the unrotated loadings (see
  *Fit indices* in Details): the model Chi Square (`chi`, with `p_chi`
  and `df`), `CFI`, `TLI`, `RMSEA` with its 90% confidence interval
  (`RMSEA_LB`, `RMSEA_UB`), `AIC`, `BIC`, `ECVI`, `RMSR`, `SRMR`, `CAF`,
  and the independence-model baseline `chi_null`, `df_null`, and
  `p_null` that CFI and TLI are computed from. Whenever the Chi Square
  is scaled, a few extra scaled-statistic fields are appended for
  advanced diagnostics. The list also carries `Fm`, the estimator's own
  objective value at the solution – on its own scale, not the
  discrepancy the Chi Square is built from.

- model_implied_R:

  The model implied correlation matrix.

- residuals:

  Residual correlations, i.e., orig_R - model_implied_R

- standardized_residuals:

  Residual correlations standardized by their bootstrap standard errors.
  Only returned, if `se = "np-boot"`.

- rot_loadings:

  Loading matrix containing the final rotated loadings. For an oblique
  rotation this is the pattern matrix – each variable's unique
  contribution from each factor, with the other factors partialled out –
  and is the matrix normally interpreted (see *Rotations*).

- Phi:

  The factor intercorrelations (only for oblique rotations).

- Structure:

  The structure matrix `rot_loadings %*% Phi`, holding the plain
  variable-factor correlations, which are inflated by the factor
  intercorrelations (only for oblique rotations).

- rotmat:

  The rotation matrix. The rotated loadings are recovered from the
  unrotated loadings as `unrot_loadings %*% rotmat` for orthogonal
  rotations and for promax, and as `unrot_loadings %*% t(solve(rotmat))`
  for the other oblique rotations.

- vars_accounted_rot:

  Matrix of explained variances and sums of squared loadings. Based on
  rotated loadings and, for oblique rotations, the factor
  intercorrelations. Same rows as `vars_accounted`; it is returned only
  when a rotation was actually applied, which never happens for a
  single-factor solution, so it always has the five-row form.

- settings:

  A list of the settings used, including `seed` (`NULL` when none was
  supplied), `input_type` (`"raw"` or `"correlation"`, what `x` was),
  and `cor_method_used` (the correlation method that actually ran,
  `NA_character_` for a correlation-matrix input, which consumes none).
  `cor_method` keeps the requested value whether or not it was used. For
  the criterion rotations fitted by gradient projection it additionally
  carries `rotation_diagnostics`, a list summarising the multi-start
  run:

  - `converged`: whether the start whose solution is returned reached
    the convergence tolerance (reported separately from the counts below
    because the returned solution is the one with the lowest criterion
    value, which need not be a converged one).

  - `n_starts_total`: the `random_starts` random starts plus the
    rational start.

  - `n_optimized`: how many of those starts were actually optimized –
    fewer than `n_starts_total` whenever the solver screens the random
    starts and optimizes only the most promising ones.

  - `n_converged`: how many optimized starts reached the convergence
    tolerance.

  - `n_distinct_minima`: how many distinct local optima those converged
    starts found; more than one means the criterion is multimodal on
    these data.

  - `criterion_spread`: the range of the criterion values they attained.

  - `criterion_best`: the criterion value of the returned solution.

  When `normalize = TRUE`, `criterion_best` and `criterion_spread` are
  evaluated on the Kaiser-normalized loadings the criterion is optimized
  on, not on the returned `rot_loadings`, so they are not directly
  comparable to a criterion recomputed from the returned loadings.

- fiml:

  Diagnostics of the FIML correlation's EM estimation, present only for
  `cor_method = "fiml"`: whether it `converged` before `fiml_max_iter`,
  how many iterations it took (`iter`), the number of distinct
  missingness patterns (`n_patterns`), and the number of cases used
  (`n`, the reported `N`). A `converged` of `FALSE` means the analysed
  correlation is the EM's last iterate, not the true FIML estimate; this
  is flagged in the printed output. Set `fiml_max_iter` and `fiml_tol`
  through
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md).

- SE:

  A named list of standard error matrices, returned only when `se` is
  not `"none"` (see *Standard errors* in Details for what each method
  assumes and when a component comes back `NA`). For `se = "np-boot"`:
  bootstrap SDs of the loadings (rotated too, if a rotation was
  applied), the residuals, and, where available, RMSEA, CFI, and TLI
  (`fit_indices`, from the transformed cases), plus – for oblique
  rotations – `Phi` and the structure coefficients; `valid_replicates`
  (and, when rotated, `valid_target_rotations`) records how many
  bootstrap replicates each of those is based on (for `fit_indices`,
  count the rows of `replicates$fit_indices` without `NA`). For
  `se = "information"` and `se = "sandwich"`: Wald-type SEs for the
  loadings, the uniquenesses, and the communalities – identical to each
  other, since communality is just `1 - uniqueness` – plus `Phi` and the
  structure coefficients for oblique rotations. `"sandwich"` stays valid
  under non-normality; `"information"` assumes the model is correctly
  specified.

- CI:

  A named list of confidence intervals of width `ci`, matching the
  components of `SE`: percentile intervals for `se = "np-boot"`, Wald
  intervals for `se = "information"` and `se = "sandwich"`. Only
  returned if `se` is not `"none"`. See *Standard errors* in Details for
  when the bootstrap gives `fit_indices`.

- replicates:

  A named list of raw bootstrap replicate arrays – the aligned loadings,
  `Phi`, structure coefficients, residuals, and RMSEA, CFI, and TLI
  behind `SE` and `CI` – with one replicate per array's last dimension
  (first dimension for `fit_indices`). Failed replicates are left `NA`.
  Populated only for `se = "np-boot"`; `NULL` otherwise.

- vcov_unrot_loadings:

  The full unrotated-loading covariance matrix behind
  `SE$unrot_loadings`: a `p * n_factors` by `p * n_factors` matrix in
  column-major `vec(Lambda)` order, with rows and columns labelled
  `"<variable>_<factor>"`. Populated for `se = "information"` and
  `se = "sandwich"` – always the unrotated block, even when a rotation
  is applied – and `NA`-filled if the analytic covariance is unreliable
  (a Heywood case or a singular information matrix); `NULL` for
  `se = "np-boot"` and `se = "none"`. It can be populated even when
  `SE$unrot_loadings` is `NA`: a weakly determined rotational
  orientation invalidates only the marginal SEs, not the underlying
  covariance.

## Details

### Estimators

The estimator is chosen with `estimator`.

- **PAF** (principal axis factoring) iteratively estimates the
  communalities and makes no distributional assumptions, which makes it
  robust and a good general-purpose default. A converged PAF solution is
  a stationary point of the ULS criterion (Jöreskog & van Thillo, 1971),
  so PAF and ULS give practically the same solution and the same fit
  indices (see *Fit indices*).

- **ML** (maximum likelihood) maximises the normal-theory likelihood. It
  is the only estimator with analytic expected-information standard
  errors (`se = "information"`), but it assumes multivariate normality
  and is the most prone to Heywood (improper) cases. Its starting values
  are set by `start_method`.

- **ULS** (unweighted least squares) minimises the sum of squared
  correlation residuals. "MINRES" (minimum residual) is the same
  estimator under a different name and returns identical results. It
  makes no normality assumption and is robust to mild non-normality.

- **DWLS** (diagonally weighted least squares) is the recommended
  estimator for ordinal data. It weights each off-diagonal correlation
  residual by the inverse asymptotic variance of the corresponding
  polychoric correlation (Muthén, du Toit, & Spisic, 1997), reproducing
  the loadings of a diagonally weighted least squares fit. It therefore
  requires raw ordinal data with `cor_method = "poly"` or `"tetra"`.
  Because the weighting follows the polychoric asymptotic covariance,
  the matrix and the weights are estimated on the listwise-complete
  cases. Its fit-index behaviour is described under *Fit indices*.

### Correlation methods

When raw data are supplied, `cor_method` selects how the correlation
matrix is computed (it is ignored when a correlation matrix is entered
directly).

- **"pearson"** (default), **"spearman"**, and **"kendall"** are passed
  to [`stats::cor()`](https://rdrr.io/r/stats/cor.html). The factor
  model assumes a Pearson correlation, but a rank correlation is
  analysed on its own scale, not converted to a Pearson-equivalent
  value; Kendall's tau in particular gives the most attenuated loadings
  of the three. For ordinal items prefer `"poly"` / `"tetra"` below,
  which estimate the correlation of the underlying continuous variables.

- **"poly"** / **"tetra"** compute polychoric / tetrachoric correlations
  for ordinal / binary data, assuming an underlying bivariate-normal
  latent variable. They use a two-step estimator. The polychoric
  asymptotic covariance that underlies both the DWLS weights and the
  scaled (sandwich) statistic relies on large-sample theory that
  degrades for empty or near-empty response-category combinations; with
  very sparse cells the resulting weights and standard errors can be
  unreliable (a warning is issued when empty cells are present), so
  interpret them with caution and consider collapsing rare categories.
  Each of the `p(p - 1)/2` variable pairs is a separate numerical
  optimisation, so a polychoric matrix takes much longer to compute than
  a Pearson one, and the difference grows quadratically in the number of
  variables; with `se = "np-boot"` the whole matrix is re-estimated for
  every bootstrap replicate.

- **"fiml"** estimates a two-stage full-information maximum-likelihood
  correlation. The saturated multivariate-normal mean and covariance are
  estimated from raw data with missing values by an EM algorithm
  assuming the data are missing at random (Yuan, Marshall, & Bentler,
  2002; Little & Rubin, 2002), and the standardized covariance is then
  analysed. The model fit indices are corrected two-stage statistics
  wherever the correction can be formed (see *Fit indices*). `"fiml"`
  uses every case and handles the missingness itself, so `use` is
  ignored; it supplies a continuous (Pearson-type) correlation only and
  is therefore not compatible with `estimator = "DWLS"`. Standard errors
  are available analytically for `estimator = "ML"` or `"ULS"` and, for
  any estimator, by the non-parametric bootstrap (see *Standard
  errors*). For multiply imputed data,
  [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)
  is the alternative route to handling missingness. Both routes assume
  the values are missing at random (MAR), and which one to prefer is
  largely practical: FIML is a single, efficient fit and is the simpler
  default when the analysis model is the whole story, whereas multiple
  imputation is more flexible when the imputation model should draw on
  auxiliary variables not in the factor model, or when the same
  imputations feed several downstream analyses.

### Rotations

A rotation transforms the unrotated loadings toward a simpler, more
interpretable pattern; all rotations are performed by rotation engines
built into the package. Orthogonal rotations keep the factors
uncorrelated, whereas oblique rotations let them correlate (returning a
pattern matrix, a structure matrix, and the factor intercorrelations
`Phi`) and are usually more realistic for psychological constructs. The
default, oblimin, is therefore oblique: if `Phi` comes back near zero
the oblique solution is essentially the orthogonal one anyway, whereas
imposing orthogonality on genuinely correlated factors distorts the
pattern.

For an oblique solution the **pattern** matrix (`rot_loadings`) holds
each variable's unique contribution from each factor, with the other
factors partialled out; it is what is normally interpreted and reported.
The **structure** matrix (`Structure = rot_loadings %*% Phi`) holds the
plain variable-factor correlations, which are inflated by the factor
intercorrelations. The two coincide only when `Phi` is the identity,
which is why an orthogonal rotation returns `rot_loadings` alone.

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
  general factor plus group factors (bifactor simple structure). It
  needs at least two group factors (`n_factors >= 3`): with two factors
  the criterion is identically zero, so no rotation is performed and the
  unrotated loadings are returned with a warning.

Oblique rotations:

- **oblimin** (the default) is quartimin with a tunable argument, `gam`
  (default 0, i.e. quartimin itself); turning it up trades some row
  simplicity for more strongly correlated factors, and can drive the
  solution toward factor collapse, so inspect `Phi` before interpreting
  a fit with `gam > 0`.

- **promax** is a fast two-step rotation: a varimax solution is raised
  to a power (controlled by `k` and `p_type`) to form a target that is
  then fitted obliquely.

- **quartimin** simplifies the variables (rows), like quartimax, so each
  loads mainly on one factor, but factors are allowed to correlate.

- **simplimax** drives the `k` smallest loadings toward zero. Its
  criterion is only piecewise smooth and strongly multimodal, which
  makes it by far the most start-dependent rotation offered here:
  different random seeds can reach noticeably different solutions.
  Because every start is fully optimised (there is no screening stage),
  raising `random_starts` to several hundred costs proportionally more
  time but buys a better optimum.

- **bentlerQ** is the oblique Bentler invariant pattern simplicity
  criterion.

- **geominQ** is the oblique geomin criterion; it handles complex
  (cross-loading) structure well but is multimodal, so it benefits from
  more `random_starts` (and uses a more thorough multi-start search
  internally).

- **bifactorQ** is the oblique (correlated) Jennrich-Bentler bifactor
  criterion, with the same two-group-factor requirement as
  **bifactorT**.

The criterion-based rotations (all except varimax and promax) are fitted
by gradient projection with `random_starts` random starts to guard
against local minima; the complexity criteria (simplimax and geominQ in
particular) are the most multimodal. The starts are drawn from the
random-number generator, so different starts can reach genuinely
different optima and such a fit is reproducible only when the generator
is controlled: pass `seed`, or call
[`base::set.seed()`](https://rdrr.io/r/base/Random.html) beforehand.

### Standard errors

`se` selects whether and how standard errors (and matching confidence
intervals) are computed. The analytic methods (`"information"` and
`"sandwich"`) cover the unrotated loadings, the uniquenesses and the
communalities and, when a rotation is applied, the rotated loadings and
– for oblique rotations – the factor correlations and the structure
coefficients. The bootstrap (`"np-boot"`) covers the unrotated loadings,
the residuals, where available RMSEA, CFI, and TLI, and, when a rotation
is applied, the rotated loadings and – for oblique rotations – the
factor correlations and the structure coefficients; it reports no
uniqueness or communality standard errors.

- **"none"** (default) computes no standard errors.

- **"information"** returns analytic standard errors from the expected
  (Fisher) information matrix of the maximum-likelihood solution, and
  therefore requires `estimator = "ML"` and `cor_method = "pearson"` (or
  `"fiml"`, see below). The rotated standard errors are obtained by
  propagating the unrotated-loading covariance through the rotation
  (Jennrich, 1973). Unlike the bootstrap it also works from a
  correlation matrix as long as `N` is supplied. These standard errors
  assume multivariate normality and a correctly specified model; under
  heavy-tailed data or model misfit they can understate the sampling
  variability, where `"sandwich"` or `"np-boot"` are more robust.

  The rotated loadings, `Phi`, the structure coefficients, the
  uniquenesses, and the communalities do not depend on how the unrotated
  solution happens to be oriented, and so are comparable across
  programs. The **unrotated** loading standard errors are not: a program
  using a different orientation will report different unrotated loading
  standard errors for the same fit. That orientation can also fail to be
  well defined – for example when two factors are of near-equal strength
  – in which case the unrotated loadings have no well-defined sampling
  distribution and their standard errors diverge. `efa_fit()` detects
  this and returns `NA` for the unrotated loading standard errors with a
  warning; the other quantities are unaffected. Use those, or
  `se = "np-boot"`, when the unrotated loadings themselves are the
  quantity of interest. This detection applies to every analytic route.
  A Heywood case (a uniqueness at its lower boundary) is more severe:
  there the whole `SE`/`CI` block comes back `NA` with a warning,
  because the standard-error approximation fails for every parameter.
  The scaled chi-square (from `"sandwich"`, or from
  `cor_method = "fiml"`) does not rely on that approximation and is
  still reported.

- **"sandwich"** returns robust ("sandwich") standard errors estimated
  from raw data, combining the estimator weight with a distribution-free
  covariance of the correlations, so it stays valid under non-normality
  and weight misspecification (Browne, 1984; Satorra & Bentler, 1994).
  It is available either for ordinal data with `cor_method = "poly"` or
  `"tetra"` and `estimator` one of `"ML"`, `"ULS"`, or `"DWLS"`, or for
  continuous data with `cor_method = "pearson"` and `estimator = "ML"`
  or `"ULS"`. It reports the same coverage as `"information"` and
  additionally fills the model fit's chi-square block with a scaled
  chi-square (see *Fit indices*); the statistic reported as `chi` is
  always the scaled-and-shifted one (flagged by `chi_scaled_type`), and
  a mean-adjusted alternative is returned alongside it as
  `chi_mean_adjusted`. Because the covariance must describe the same
  cases as the correlation matrix, the sandwich (like
  `estimator = "DWLS"`) is computed on the listwise-complete cases; on
  data with missing values the reported `N`, the correlation matrix, and
  the point estimate therefore reflect the complete cases regardless of
  `use`.

- **"np-boot"** draws a non-parametric (case-resampling) bootstrap and
  needs raw data. It is the most general method – available for any
  `estimator`, `rotation`, and `cor_method` – and the most robust to
  non-normality and misfit, at the cost of speed; its intervals are
  bootstrap percentile intervals. The replicates run with the `future`
  framework; by default they run sequentially, but registering a plan
  with
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  runs them in parallel instead, as the examples show. With a fixed
  `seed` the bootstrap is reproducible and yields the same result
  regardless of the number of workers. Under `cor_method = "fiml"` each
  resample also re-runs the EM moment estimation and is therefore slow,
  so a smaller `b_boot` may be advisable.

  The intervals of RMSEA, CFI, and TLI come from a second set of fits. A
  resample of the raw cases adds sampling misfit to the misfit of the
  sample, which moves these intervals toward worse fit. Thus,
  `efa_fit()` first transforms the cases so that their misfit is equal
  to the estimated population misfit. Then it resamples the same rows
  (Yuan, Hayashi, & Yanagihara, 2007). If the RMSEA is 0, the
  transformed cases fit the model exactly. In simulations with normal
  data, the intervals of RMSEA and CFI had good coverage. Use them with
  caution if the CFI is below approximately .76 or the RMSEA is above
  approximately .10 (Zhang & Savalei, 2016). TLI is CFI rescaled by a
  constant, so its interval behaves in the same way. The transformation
  needs complete cases and Pearson correlations. Thus, these intervals
  need `cor_method = "pearson"`, no missing values in the analysed cases
  (complete data, or listwise deletion), and a chi-square with df \> 0.
  Otherwise, CFI and TLI are printed without an interval. The other fit
  indices get no bootstrap interval: SRMR had poor coverage with either
  bootstrap, and Zhang and Savalei (2016) did not study the others. The
  second set of fits is not rotated, so it usually takes less time than
  the first.

The analytic methods (`"information"` and `"sandwich"`) are not
available with the `"promax"` or `"simplimax"` rotations, which have no
supported analytic route for the rotated standard errors; use
`"np-boot"` there. They are not available for PAF either; use ULS, which
gives the same solution at convergence, or `"np-boot"`. Under
`cor_method = "fiml"`, `"information"` and `"sandwich"` both return, for
`estimator = "ML"` or `"ULS"`, corrected two-stage sandwich standard
errors (Yuan & Bentler, 2000; Savalei & Bentler, 2009).

### Fit indices

For ML, ULS, and PAF, `efa_fit()` returns the model chi-square (with its
p-value and degrees of freedom), the Comparative Fit Index (CFI;
Bentler, 1990), the Tucker-Lewis Index (TLI, also called the non-normed
fit index; Tucker & Lewis, 1973), the Root Mean Square Error of
Approximation (RMSEA) with its 90% confidence interval (Browne & Cudeck,
1992), the Akaike and Bayesian Information Criteria (AIC, BIC), the
Expected Cross-Validation Index (ECVI; Browne & Cudeck, 1989), the Root
Mean Squared Residual (RMSR), the Standardized Root Mean Squared
Residual (SRMR; Bentler, 1995), and the common-part-accounted-for (CAF)
index (Lorenzo-Seva, Timmerman, & Kiers, 2011). They come with the
independence-baseline statistics `chi_null`, `df_null`, and `p_null`. On
the unscaled paths `chi_null` is Bartlett's test of sphericity; a scaled
chi-square, and the two-stage statistic of a `cor_method = "fiml"` fit,
each carry their own baseline instead. The baseline degrees of freedom
`df_null` are \\p(p - 1)/2\\ for \\p\\ variables.

The print method shows the chi-square, CFI, TLI, RMSEA, and SRMR;
[`summary()`](https://rdrr.io/r/base/summary.html) adds AIC, BIC, and
ECVI, which are only useful to compare models. With `se = "np-boot"`,
CFI and TLI also show their bootstrap interval (see *Standard errors*).
RMSEA always shows its analytic 90% interval; its bootstrap interval is
in `CI$fit_indices`. RMSR is the root mean square of the off-diagonal
residuals; SRMR rescales it by a fixed factor that depends only on the
number of variables, and both are `NA` if a residual is unavailable.
(psych's `rms` uses a different divisor and equals RMSR\\/\sqrt{2}\\.)
CAF is returned for
[`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)
but not printed: although Lorenzo-Seva et al. describe it as ranging
from 0 to 1, a correctly specified model gives a CAF near 0.5, because
the KMO of residuals that hold only sampling noise is near 0.5. Read it
as a relative measure. It is `NA` where that KMO is undefined.

The model chi-square is the Bartlett-corrected discrepancy (matching
[`stats::factanal()`](https://rdrr.io/r/stats/factanal.html) for ML).
For ULS and PAF it is the same maximum-likelihood discrepancy, evaluated
at the fitted solution, as in
[`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html); its chi-square
distribution is then an approximation. It is not the least-squares
criterion lavaan reports as its standard ULS test statistic, so the two
are not comparable. AIC and BIC are built on this chi-square and can
therefore be negative; ECVI (built on the same chi-square plus a
non-negative penalty) cannot. Because of the Bartlett correction, ECVI
differs slightly from the uncorrected Browne-Cudeck form reported by
lavaan and Mplus (their AIC/BIC use an unrelated log-likelihood-based
formula, so they are not comparable). On the unscaled paths, CFI, TLI,
and RMSEA are computed on the uncorrected \\(N - 1)\\ scale rather than
on the Bartlett-corrected chi-square, so you cannot recompute one from
the other by hand there; on the scaled (sandwich) and FIML paths, the
reported chi and chi_null are exactly the pair the indices use.

Two cases differ from the full set:

- **DWLS** by default returns only RMSR, SRMR, CAF, and df, because the
  ordinary maximum-likelihood discrepancy is not its fit function. When
  `se = "sandwich"`, a scaled chi-square and the CFI, TLI, and RMSEA
  derived from it are reported; that statistic is a two-stage correction
  applied to the polychoric correlation residuals (Browne, 1984), not
  identical to the full WLSMV test of lavaan or Mplus, which also
  corrects for the response thresholds.

- **`cor_method = "fiml"`** reports two-stage-corrected statistics
  (Yuan, Marshall, & Bentler, 2002). The correction itself can be
  degenerate – typically with a small sample, a high proportion of
  missing values, or near-collinear variables – in which case an
  uncorrected likelihood-ratio statistic is reported in its place with a
  warning, and the print methods label that line `uncorrected` rather
  than `scaled`; read its p-value and the CFI, TLI, and RMSEA derived
  from it as indicative only.

Whenever the chi-square is a scaled one (`se = "sandwich"`, or a
`cor_method = "fiml"` fit), AIC, BIC, and ECVI are `NA`; see the
`fit_indices` entry in Value for the additional components then
returned. Beyond that, the chi-square and everything derived from it are
`NA` whenever the statistic itself is undefined: when `N` is not
supplied, when the model is underidentified (a negative df), and when a
positive `N` is too small relative to the number of variables and
factors for the small-sample correction to remain valid. Each case
raises its own warning. The residual summaries (RMSR, SRMR, CAF) and the
degrees of freedom are still returned there, but residual size does not
establish that a model is identified: below zero degrees of freedom a
near-zero residual is an artefact of over-parameterisation, not close
fit. A just-identified model (df = 0) has no test of fit, so the print
method shows only SRMR and df.

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

# Principal axis factoring with oblimin rotation (the defaults)
mod_paf <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500)
mod_paf
#> 
#> EFA performed with estimator = 'PAF' and rotation = 'oblimin' (N = 500).
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
#> χ²(102) = 123.92, p = .069
#> CFI: .99
#> TLI: .98
#> RMSEA [90% CI]: .02 [.00; .03]
#> SRMR: .02
summary(mod_paf)
#> 
#> EFA performed with estimator = 'PAF' and rotation = 'oblimin' (N = 500).
#> 
#> ── Model Diagnostics ───────────────────────────────────────────────────────────
#> 
#> N: 500
#> Rotation local optima: 1 distinct from 6 of 101 starts
#> Heywood cases: 0
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
#> χ²(102) = 123.92, p = .069
#> CFI: .99
#> TLI: .98
#> RMSEA [90% CI]: .02 [.00; .03]
#> SRMR: .02
#> AIC: -80.08
#> BIC: -509.97
#> ECVI: 0.52
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

# ML estimation with promax rotation
mod_ml <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
                  estimator = "ML", rotation = "promax")
mod_ml
#> 
#> EFA performed with estimator = 'ML' and rotation = 'promax' (N = 500).
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>        F1     F2     F3    h2    u2
#> V1   -.051   .030   .624  .373  .627
#> V2    .001   .077   .470  .274  .726
#> V3    .064   .064   .441  .280  .720
#> V4    .101  -.008   .551  .379  .621
#> V5    .157  -.008   .429  .290  .710
#> V6   -.070  -.053   .709  .402  .598
#> V7    .005   .532   .089  .355  .645
#> V8   -.015   .573   .037  .345  .655
#> V9    .035   .545   .008  .328  .672
#> V10  -.031   .676  -.063  .385  .615
#> V11   .019   .353   .231  .296  .704
#> V12   .023   .663  -.026  .437  .563
#> V13   .620   .088  -.076  .397  .603
#> V14   .547  -.066   .086  .320  .680
#> V15   .559   .133  -.073  .363  .637
#> V16   .557  -.052   .091  .345  .655
#> V17   .663  -.046  -.018  .390  .610
#> V18   .556   .002   .049  .349  .651
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
#> F2   .619  1.000
#> F3   .647   .631  1.000
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F2     F3
#> SS loadings        2.204  2.079  2.025
#> Prop Tot Var        .122   .115   .112
#> Cum Prop Tot Var    .122   .238   .350
#> Prop Comm Var       .349   .330   .321
#> Cum Prop Comm Var   .349   .679  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> χ²(102) = 123.75, p = .070
#> CFI: .99
#> TLI: .98
#> RMSEA [90% CI]: .02 [.00; .03]
#> SRMR: .03

# Tuning knobs are supplied through the control objects. Here the SPSS preset is
# used for the estimation and rotation, with the maximum PAF iterations raised.
mod_spss <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
                    rotation = "promax",
                    estimate_control = estimate_control(type = "SPSS", max_iter = 500),
                    rotate_control = rotate_control(type = "SPSS"))
mod_spss
#> 
#> EFA performed with estimator = 'PAF' and rotation = 'promax' (N = 500).
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
#> χ²(102) = 123.92, p = .069
#> CFI: .99
#> TLI: .98
#> RMSEA [90% CI]: .02 [.00; .03]
#> SRMR: .02

# Analytic (expected-information) standard errors of an ML fit
ML_info <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
                   estimator = "ML", se = "information")
ML_info
#> 
#> EFA performed with estimator = 'ML' and rotation = 'oblimin' (N = 500).
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
#> SRMR: .03
summary(ML_info)
#> 
#> EFA performed with estimator = 'ML' and rotation = 'oblimin' (N = 500).
#> 
#> ── Model Diagnostics ───────────────────────────────────────────────────────────
#> 
#> N: 500
#> Rotation local optima: 1 distinct from 6 of 101 starts
#> Heywood cases: 0
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
#> ── 95% Wald CIs (expected information) for salient rotated loadings ────────────
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
#> ── 95% Wald CIs (expected information) for factor intercorrelations ────────────
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
#> SRMR: .03
#> AIC: -80.25
#> BIC: -510.14
#> ECVI: 0.52
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
                    estimator = "DWLS", se = "sandwich")
#> Warning: 18 variable pairs have an empty response-category combination despite a
#> non-negligible expected count.
#> ✖ Affected pairs: "ethR_3-recR_3", "ethR_5-heaR_2", "finR_3-finR_6",
#>   "finR_3-heaR_1", "finR_3-socR_3", and 13 more.
#> ℹ The polychoric asymptotic covariance (and any DWLS weights or robust standard
#>   errors derived from it) can be unreliable for such structurally sparse cells;
#>   interpret them with caution and consider collapsing rare response categories
#>   in these variables.
#> Warning: Analytic standard errors could not be computed for all parameters.
#> ℹ The factor solution's rotational orientation is only weakly determined (two
#>   canonical variances nearly coincide), so the unrotated loadings have no
#>   well-defined standard error. The rotated loadings and communalities are
#>   gauge-invariant and unaffected.
DWLS_rob
#> 
#> EFA performed with estimator = 'DWLS' and rotation = 'oblimin'
#>   (polychoric correlations, N = 3123).
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
#> SRMR: .03
summary(DWLS_rob)
#> 
#> EFA performed with estimator = 'DWLS' and rotation = 'oblimin'
#>   (polychoric correlations, N = 3123).
#> 
#> ── Model Diagnostics ───────────────────────────────────────────────────────────
#> 
#> N: 3123
#> Rotation local optima: 1 distinct from 6 of 101 starts
#> Heywood cases: 0
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
#> ── 95% Wald CIs (robust sandwich) for salient rotated loadings ─────────────────
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
#> ── 95% Wald CIs (robust sandwich) for factor intercorrelations ─────────────────
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
                  estimator = "ML", se = "sandwich")
ML_rob
#> 
#> EFA performed with estimator = 'ML' and rotation = 'none' (N = 810).
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
#> SRMR: .01
summary(ML_rob)
#> 
#> EFA performed with estimator = 'ML' and rotation = 'none' (N = 810).
#> 
#> ── Model Diagnostics ───────────────────────────────────────────────────────────
#> 
#> N: 810
#> Heywood cases: 0
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
#> ── 95% Wald CIs (robust sandwich) for salient unrotated loadings ───────────────
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
#> ── Simple Structure Diagnostics ────────────────────────────────────────────────
#> 
#> All items have a salient loading (|loading| >= .300) and no cross-loading
#> (|loading| >= .300), all primary-loading gaps are at least .200, and all
#> factors have at least 3 salient indicators.
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
#> SRMR: .01
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
efa_fiml
#> 
#> EFA performed with estimator = 'ML' and rotation = 'none'
#>   (two-stage FIML correlations, N = 810).
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
#> SRMR: .01
# }

if (FALSE) { # \dontrun{
# Bootstrap standard errors from raw data, reproducible via a fixed seed and run
# in parallel across replicates. future::plan() returns the plan it replaces, so
# on.exit() puts the session back as it was -- also if the fit fails.
efa_boot <- local({
  old_plan <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old_plan), add = TRUE)
  efa_fit(GRiPS_raw, n_factors = 1, estimator = "PAF",
          se = "np-boot", b_boot = 1000, seed = 42)
})
} # }
```
