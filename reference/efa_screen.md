# Screen data for exploratory factor analysis

Computes a set of diagnostics describing how suitable a correlation
matrix (or raw data) is for exploratory factor analysis: the
Kaiser-Meyer-Olkin (KMO) measure of sampling adequacy overall and per
variable, Bartlett's test of sphericity, the determinant and condition
number of the correlation matrix, and the squared multiple correlation
(SMC) of each variable with all the others.

## Usage

``` r
efa_screen(
  x,
  N = NA,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
  mcd_alpha = 0.5,
  outlier_cutoff = 0.975,
  seed = NULL
)
```

## Source

Bartlett, M. S. (1951). The effect of standardization on a Chi-square
approximation in factor analysis. Biometrika, 38, 337-344.

Belsley, D. A., Kuh, E. & Welsch, R. E. (1980). Regression diagnostics:
Identifying influential data and sources of collinearity. Wiley.

Croux, C. & Haesbroeck, G. (1999). Influence function and efficiency of
the minimum covariance determinant scatter matrix estimator. Journal of
Multivariate Analysis, 71, 161-190.

Field, A. (2018). Discovering statistics using IBM SPSS statistics (5th
ed.). Sage.

Henze, N. & Zirkler, B. (1990). A class of invariant consistent tests
for multivariate normality. Communications in Statistics - Theory and
Methods, 19, 3595-3617.

Kaiser, H. F. (1970). A second generation little jiffy. Psychometrika,
35, 401-415.

Kaiser, H. F. & Rice, J. (1974). Little jiffy, mark IV. Educational and
Psychological Measurement, 34, 111-117.

Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis
with applications. Biometrika, 57, 519-530.

Mardia, K. V. (1974). Applications of some measures of multivariate
skewness and kurtosis in testing normality and robustness studies.
Sankhya B, 36, 115-128.

Pison, G., Van Aelst, S. & Willems, G. (2002). Small sample corrections
for LTS and MCD. Metrika, 55, 111-123.

Rhemtulla, M., Brosseau-Liard, P. E. & Savalei, V. (2012). When can
categorical variables be treated as continuous? A comparison of robust
continuous and categorical SEM estimation methods under suboptimal
conditions. Psychological Methods, 17, 354-373.

Rousseeuw, P. J. & Van Driessen, K. (1999). A fast algorithm for the
minimum covariance determinant estimator. Technometrics, 41, 212-223.

## Arguments

- x:

  data.frame or matrix. Data frame or matrix of raw data, or a matrix of
  correlations.

- N:

  numeric. The number of observations. Only needs to be specified when a
  correlation matrix is supplied; it is required for Bartlett's test of
  sphericity and is taken from the data when raw data are supplied.
  Default is `NA`.

- use:

  character. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) if raw data are
  supplied. Default is `"pairwise.complete.obs"`.

- cor_method:

  character. Correlation computed from raw data: `"pearson"`,
  `"spearman"`, or `"kendall"` (passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)), or `"poly"` /
  `"tetra"` for polychoric / tetrachoric correlations of ordinal /
  binary data (a two-step estimator with no empty-cell continuity
  correction). Default is `"pearson"`.

- mcd_alpha:

  numeric. The proportion of observations covered by the minimum
  covariance determinant (MCD) subset used for the robust outlier
  diagnostics, in \[0.5, 1\]. The default, `0.5`, gives the most robust
  (highest-breakdown) estimate; larger values trade robustness for
  efficiency. Only used when raw data are supplied.

- outlier_cutoff:

  numeric. The probability defining the chi-square cutoff for flagging a
  multivariate outlier: an observation is flagged when its squared
  robust distance exceeds `qchisq(outlier_cutoff, p)` for `p` variables.
  Default is `0.975`. Only used when raw data are supplied.

- seed:

  integer. Optional seed for the random subsets drawn by the MCD
  algorithm, making the outlier diagnostics reproducible; the caller's
  random-number-generator state is left unchanged. Default is `NULL`.
  Only used when raw data are supplied.

## Value

An object of class `efa_screen`, a list containing:

- kmo:

  A list with the overall KMO (`KMO`) and the per-variable KMO
  (`KMO_i`).

- bartlett:

  A list with Bartlett's chi-square statistic (`chisq`), its `p_value`,
  and its degrees of freedom (`df`); `chisq` and `p_value` are `NA` when
  `N` is too small for the Bartlett correction. `NULL` when `N` is
  unavailable.

- determinant:

  The determinant of the correlation matrix.

- condition:

  The condition number of the correlation matrix (largest over smallest
  eigenvalue).

- smc:

  The per-variable squared multiple correlations.

- per_item:

  A data frame with one row per variable (row names are the variable
  names) holding the per-item diagnostics: `variance` (over the
  available values), `missing` (percentage of missing values), `smc`,
  `kmo_i`, and `flags` (a comma-separated list of any
  sparse/empty-category issues, the empty string `""` for a categorical
  variable with none, and `NA` for a variable treated as continuous).
  `NULL` when a correlation matrix is supplied.

- normality:

  Multivariate-normality diagnostics computed from the complete cases of
  the raw data: a list with `mardia` (Mardia's multivariate skewness
  statistic `skewness`, its degrees of freedom `skewness_df` and p-value
  `skewness_p`, a `small_sample` flag recording whether the small-sample
  skewness correction was applied, the standardised `kurtosis` statistic
  and its p-value `kurtosis_p`, and the underlying coefficients `b1p`
  and `b2p`), `hz` (the Henze-Zirkler `statistic` and its `p_value`),
  and `n_complete` (the number of complete cases used). When the
  complete-case covariance is singular the tests are skipped and a
  classed note (of class `efa_screen_no_mvn`) is returned instead,
  alongside a warning. `NULL` when a correlation matrix is supplied.

- outliers:

  Multivariate-outlier diagnostics from the complete cases of the raw
  data: a list with `distances` (each complete case's robust Mahalanobis
  distance, named by its row number in the supplied data), `cutoff` (the
  flagging threshold on the distance scale,
  `sqrt(qchisq(outlier_cutoff, p))`, directly comparable to
  `distances`), `flagged` (the row numbers, in the supplied data, whose
  robust distance exceeds `cutoff`), `center` and `cov` (the robust
  location and scatter underlying the distances), `method` (`"mcd"` for
  the robust estimate, or `"classical"` when the fallback was used), and
  `n_complete` (the number of complete cases). When neither a robust nor
  a classical covariance can be formed a classed note (of class
  `efa_screen_no_outliers`) is returned instead. `NULL` when a
  correlation matrix is supplied.

- categories:

  A named list with, for each variable treated as categorical, the
  response-category counts (in category order); `NA` for a variable
  treated as continuous. `NULL` when a correlation matrix is supplied.

- note:

  A classed note explaining that the raw-data diagnostics need raw data;
  `NULL` when raw data are supplied.

- settings:

  A list of the settings used, including `n_obs`, the number of rows in
  the supplied raw data (`NA` for a correlation-matrix input).

## Details

The diagnostics are computed from the analysis correlation matrix \\R\\:

- KMO:

  The Kaiser-Meyer-Olkin measure of sampling adequacy (Kaiser, 1970;
  Kaiser & Rice, 1974), overall and for each variable. Larger values (a
  rough floor of .50) indicate greater suitability for factor analysis;
  see
  [`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md).

- Bartlett:

  Bartlett's (1951) test of sphericity, the likelihood-ratio test of the
  hypothesis that \\R\\ is an identity matrix, with \\df = p(p - 1)/2\\
  for \\p\\ variables. It requires the sample size `N`; when `N` is
  unavailable (a correlation matrix supplied without `N`) the test is
  skipped with a warning, `$bartlett` is `NULL`, and the remaining
  diagnostics are still returned. See
  [`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md).

- Determinant:

  The determinant of \\R\\. A value near zero signals extreme
  multicollinearity or a (near-)singular matrix; as a rough guide, a
  determinant below about 0.00001 is commonly taken as a sign of
  multicollinearity (Field, 2018).

- Condition number:

  The ratio of the largest to the smallest eigenvalue of \\R\\. Large
  values indicate an ill-conditioned matrix with near-linear
  dependencies among the variables. Its square root is the condition
  index: an index above 30 is a conventional sign of strong
  multicollinearity, and 10 to 30 of moderate multicollinearity
  (Belsley, Kuh & Welsch, 1980).

- SMC:

  The squared multiple correlation of each variable with all the others,
  \\1 - 1/(R^{-1})\_{ii}\\. A low SMC flags a variable that shares
  little variance with the rest of the set.

- Variance and missing data:

  The variance of each variable (over its available values) and the
  percentage of missing values. These, and the category tabulation
  below, are computed column by column from the supplied data using
  every non-missing value, and so do not depend on `use`, which governs
  only the correlation matrix. Under a listwise `use` (`"complete.obs"`
  / `"na.or.complete"`) the correlation matrix and `N` are based on the
  complete cases, while the missingness is reported over every supplied
  row so that it explains why `N` is smaller; the number of supplied
  rows is recorded in `settings$n_obs`. Available only from raw data.

- Categories:

  For each variable with fewer than ten distinct values (treated as
  categorical), the response-category counts, flagging a *sparse*
  category (fewer than five responses) and, for integer-coded variables,
  an *empty* interior category (an unused category between the smallest
  and largest observed value). A variable with ten or more distinct
  values is treated as continuous and is not tabulated. As a rough
  guide, items with fewer than five response categories are better
  analysed with an ordinal estimator (polychoric correlations with
  categorical least squares) than with normal-theory maximum likelihood
  (Rhemtulla et al., 2012). Available only from raw data.

- Multivariate normality:

  Two tests of multivariate normality computed from the complete cases
  of the raw data: Mardia's (1970) multivariate skewness and kurtosis,
  and the Henze-Zirkler (1990) omnibus test. A small p-value indicates a
  departure from multivariate normality, a reason to prefer robust or
  ordinal estimation over normal-theory maximum likelihood. Available
  only from raw data, and skipped with a note if the complete-case
  covariance is singular.

- Outliers:

  Multivariate outliers flagged by their robust Mahalanobis distance. A
  high-breakdown robust location and scatter are estimated from the
  complete cases with the fast minimum covariance determinant (MCD)
  algorithm (Rousseeuw & Van Driessen, 1999), using a subset covering a
  proportion `mcd_alpha` of the observations; an observation whose
  squared robust distance exceeds `qchisq(outlier_cutoff, p)` is
  flagged. With too few complete cases (\\n \le 2p\\) or collinear
  variables the robust covariance is undefined, so the classical
  Mahalanobis distance is used instead with a warning; if even that
  covariance is singular the diagnostic is skipped with a note.
  Available only from raw data.

## See also

[`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md)
and
[`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md)
for the individual suitability measures, and
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
for factor retention criteria.

Other factor analysis suitability:
[`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md),
[`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md),
[`print.efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_screen.md)

## Examples

``` r
# From a correlation matrix (supply N for Bartlett's test of sphericity)
efa_screen(test_models$baseline$cormat, N = 500)
#> 
#> ── Sampling adequacy and sphericity ────────────────────────────────────────────
#> 
#> ✔ The overall KMO value for your data is marvellous (Overall KMO = 0.916).
#> These data are probably suitable for factor analysis.
#> 
#> ✔ The Bartlett's test of sphericity was significant at an alpha level of .05.
#> These data are probably suitable for factor analysis.
#> 𝜒²(153) = 2173.28, p < .001
#> 
#> ── Multicollinearity ───────────────────────────────────────────────────────────
#> 
#> ✔ Determinant: 0.0121. No concern (a value near 0 signals multicollinearity).
#> ✔ Condition number: 11.680. No concern (large values signal near-collinear variables).
#> 
#> ── Per-variable diagnostics ────────────────────────────────────────────────────
#> 
#>       MSA   SMC
#> V1  0.900 0.309
#> V2  0.914 0.250
#> V3  0.924 0.260
#> V4  0.932 0.315
#> V5  0.923 0.269
#> V6  0.891 0.305
#> V7  0.928 0.304
#> V8  0.919 0.288
#> V9  0.916 0.282
#> V10 0.892 0.303
#> V11 0.928 0.278
#> V12 0.908 0.352
#> V13 0.922 0.322
#> V14 0.905 0.283
#> V15 0.924 0.308
#> V16 0.934 0.283
#> V17 0.907 0.313
#> V18 0.923 0.299
#> 
#> ── Recommendations ─────────────────────────────────────────────────────────────
#> 
#> ✔ The data appear suitable for factor analysis.
#> ℹ Per-item variance, missing-data, category, normality, and outlier diagnostics
#>   require raw data; only a correlation matrix was supplied.

# From raw data (N is taken from the data)
efa_screen(GRiPS_raw)
#> 
#> ── Sampling adequacy and sphericity ────────────────────────────────────────────
#> 
#> ✔ The overall KMO value for your data is marvellous (Overall KMO = 0.955).
#> These data are probably suitable for factor analysis.
#> 
#> ✔ The Bartlett's test of sphericity was significant at an alpha level of .05.
#> These data are probably suitable for factor analysis.
#> 𝜒²(28) = 5054.06, p < .001
#> 
#> ── Multicollinearity ───────────────────────────────────────────────────────────
#> 
#> ✔ Determinant: 0.00188. No concern (a value near 0 signals multicollinearity).
#> ✔ Condition number: 24.340. No concern (large values signal near-collinear variables).
#> 
#> ── Per-variable diagnostics ────────────────────────────────────────────────────
#> 
#>           variance missing   SMC   MSA  flags
#> fun          1.427       0 0.609 0.955       
#> friends      1.309       0 0.685 0.950       
#> enjoy        1.154       0 0.715 0.947 sparse
#> hurt         1.136       0 0.557 0.968 sparse
#> part         1.271       0 0.636 0.955       
#> commonly     1.099       0 0.653 0.958       
#> chances      1.389       0 0.592 0.961       
#> attracted    1.259       0 0.679 0.950 sparse
#> 
#> ── Multivariate normality ──────────────────────────────────────────────────────
#> 
#> ✖ Mardia's skewness: 𝜒²(120) = 910.1, p < .001.
#> ✖ Mardia's kurtosis: z = 49.32, p < .001.
#> ✖ Henze-Zirkler: HZ = 11.65, p < .001.
#> These data depart from multivariate normality.
#> 
#> ── Outliers ────────────────────────────────────────────────────────────────────
#> 
#> ℹ 214 of 810 observations were flagged as multivariate outliers (robust distance > 4.19).
#> 
#> ── Recommendations ─────────────────────────────────────────────────────────────
#> 
#> ! These data depart from multivariate normality; normal-theory standard errors
#>   and fit statistics may be biased - prefer robust (sandwich) or bootstrapped
#>   standard errors.
#> ! Bartlett's test is significant, but it assumes multivariate normality and
#>   grows more sensitive as N increases; because these data are non-normal, treat
#>   it as uninformative here and rely on the KMO.
#> ! 3 variables have a sparse response category (< 5 responses): enjoy, hurt, and
#>   attracted; a low-frequency category can destabilise polychoric estimates -
#>   consider collapsing it into an adjacent category.
#> ! 214 observations were flagged as potential multivariate outliers; inspect
#>   them (see `$outliers$flagged`) before down-weighting or excluding.
```
