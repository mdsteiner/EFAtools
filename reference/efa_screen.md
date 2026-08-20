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

Cochran, W. G. (1954). Some methods for strengthening the common
\\\chi^2\\ tests. Biometrics, 10, 417-451.

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
  correlations. At least three variables are required, and no variable
  may be a perfect linear combination of the others.

- N:

  numeric. The number of observations. Only needs to be specified when a
  correlation matrix is supplied; it is required for Bartlett's test of
  sphericity and is taken from the data when raw data are supplied.
  Default is `NA`.

- use:

  character. The missing-data policy for raw data. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) for `"pearson"`,
  `"spearman"`, and `"kendall"`; for `"poly"` / `"tetra"` the same
  policies are applied to the raw data before the polychoric estimation,
  where `"all.obs"` and `"everything"` abort on a missing value instead
  of returning `NA` correlations. Default is `"pairwise.complete.obs"`.

- cor_method:

  character. Correlation computed from raw data: `"pearson"`,
  `"spearman"`, or `"kendall"` (passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)), or `"poly"` /
  `"tetra"` for polychoric / tetrachoric correlations of ordinal /
  binary data (a two-step estimator). A Spearman or Kendall matrix is
  screened on its own metric and is not transformed to the Pearson scale
  the factor model assumes; Kendall's tau in particular is not a Pearson
  correlation. Only `"poly"` and `"tetra"` accept ordered factors as
  well as numeric response codes; the other three need numeric data. An
  unordered factor or a character column is refused, because its
  categories carry no response order and the alphabetical order of its
  levels is not a safe substitute. Default is `"pearson"`.

- mcd_alpha:

  numeric. The proportion of observations covered by the minimum
  covariance determinant (MCD) subset used for the robust outlier
  diagnostics, in \[0.5, 1\]. The default, `0.5`, gives the most robust
  (highest-breakdown) estimate; larger values trade robustness for
  efficiency. Only used when raw data are supplied.

- outlier_cutoff:

  numeric. The probability defining the chi-square cutoff for flagging a
  multivariate outlier: an observation is flagged when its squared
  robust distance exceeds `qchisq(outlier_cutoff, p)` for `p` variables,
  in \[0.5, 0.9999\]. Default is `0.975`. Only used when raw data are
  supplied.

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
  and its degrees of freedom (`df`); `chisq` and `p_value` are `NA`,
  with a warning, when `N` is too small for the Bartlett correction.
  `NULL` when `N` is unavailable.

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
  `skewness_p`, the standardised `kurtosis` statistic and its p-value
  `kurtosis_p` (both `NA` when the exact null variance is zero, as at
  \\n = p + 1\\, where the kurtosis coefficient is constant and gives no
  information), and the underlying coefficients `b1p` and `b2p`), `hz`
  (the Henze-Zirkler `statistic` and its `p_value`), and `n_complete`
  (the number of complete cases used). When the complete-case covariance
  is singular the tests are skipped and a classed note (of class
  `efa_screen_no_mvn`) is returned instead, alongside a warning. When
  the Henze-Zirkler null approximation degenerates at many variables,
  `hz` keeps its `statistic`, has an `NA` `p_value`, and carries a
  `message` and the class `efa_screen_no_hz`, also alongside a warning;
  Mardia's tests are unaffected. `NULL` when a correlation matrix is
  supplied.

- outliers:

  Multivariate-outlier diagnostics from the complete cases of the raw
  data: a list with `distances` (each complete case's robust Mahalanobis
  distance, named by its row number in the supplied data), `cutoff` (the
  flagging threshold on the distance scale,
  `sqrt(qchisq(outlier_cutoff, p))`, directly comparable to
  `distances`), `flagged` (the row numbers, in the supplied data, whose
  robust distance exceeds `cutoff`), `center` and `cov` (the robust
  location and scatter underlying the distances), `method` (`"mcd"` for
  the robust estimate, or `"classical"` when the fallback was used),
  `fallback_reason` (why the robust estimate was unavailable, or `NULL`
  when `method` is `"mcd"`), and `n_complete` (the number of complete
  cases). When neither a robust nor a classical covariance can be formed
  a classed note (of class `efa_screen_no_outliers`) is returned
  instead. `NULL` when a correlation matrix is supplied.

- categories:

  A named list with, for each variable treated as categorical, the
  response-category counts (in category order, labelled by the original
  levels for a factor or character column and by the value itself
  otherwise); `NA` for a variable treated as continuous. `NULL` when a
  correlation matrix is supplied.

- note:

  A classed note explaining that the raw-data diagnostics need raw data;
  `NULL` when raw data are supplied.

- settings:

  A list of the settings used, including `n_obs`, the number of rows in
  the supplied raw data (`NA` for a correlation-matrix input),
  `outlier_cutoff`, the probability behind the outlier flagging
  threshold, `mcd_alpha`, the coverage of the minimum covariance
  determinant subset, and `seed` (`NULL` when none was supplied).

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
  diagnostics are still returned. A supplied `N` that is too small
  relative to \\p\\ for the correction \\N - 1 - (2p + 5)/6\\ to be
  positive leaves the statistic itself `NA`, also with a warning. See
  [`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md).

- Determinant:

  The determinant of \\R\\, reported as a number. It is the product of
  the \\p\\ eigenvalues of \\R\\, so it falls geometrically as variables
  are added even when every eigenvalue stays far from zero: on a clean
  one-factor pool with all loadings 0.5 it goes from 0.24 at \\p = 10\\
  to 6.7e-07 at \\p = 60\\ while the condition index moves only from 2.1
  to 4.6. A fixed cut-off on it, such as the 0.00001 that is commonly
  quoted (Field, 2018), is therefore a statement about the number of
  variables as much as about the data, and this package does not use
  one. The multicollinearity verdict rests on the condition index below.

- Condition number:

  The ratio of the largest to the smallest eigenvalue of \\R\\. Large
  values indicate an ill-conditioned matrix with near-linear
  dependencies among the variables. Its square root is the condition
  index: an index above 30 is a conventional sign of strong
  multicollinearity, and 10 to 30 of moderate multicollinearity
  (Belsley, Kuh & Welsch, 1980). The index is a ratio of eigenvalues, so
  unlike the determinant it does not move with \\p\\ on its own; it
  carries the multicollinearity verdict in the printed report and the
  multicollinearity recommendation, which fires above 30.

- SMC:

  The squared multiple correlation of each variable with all the others,
  \\1 - 1/(R^{-1})\_{ii}\\. A low SMC flags a variable that shares
  little variance with the rest of the set.

- Variance and missing data:

  The variance of each variable (over its available values) and the
  percentage of missing values. Ordered-factor columns are recoded to
  their integer level codes, so `variance` and the empty-category check
  below refer to those codes (the category counts themselves are
  labelled by the original levels). Such columns reach the correlation
  matrix only with `cor_method = "poly"` or `"tetra"`; a Pearson,
  Spearman, or Kendall correlation of a factor frame is refused, as is
  any correlation of an unordered factor or a character column. These,
  and the category tabulation below, are computed column by column from
  the supplied data using every non-missing value, and so do not depend
  on `use`, which governs only the correlation matrix. Under a listwise
  `use` (`"complete.obs"` / `"na.or.complete"`) the correlation matrix
  and `N` are based on the complete cases, while the missingness is
  reported over every supplied row so that it explains why `N` is
  smaller; the number of supplied rows is recorded in `settings$n_obs`.
  Available only from raw data.

- Categories:

  For each variable with fewer than ten distinct values (treated as
  categorical), the response-category counts, flagging a *sparse*
  category (fewer than five responses; a heuristic of this package,
  borrowing the conventional minimum cell count of five from
  contingency-table analysis (Cochran, 1954), because a polychoric
  correlation is estimated from the bivariate response table, in which a
  category with very few responses leaves its threshold poorly
  determined) and, for integer-coded variables, an *empty* interior
  category (an unused category between the smallest and largest observed
  value). A variable with ten or more distinct values is treated as
  continuous and is not tabulated. As a rough guide, items with fewer
  than five response categories are better analysed with an ordinal
  estimator (polychoric correlations with categorical least squares)
  than with normal-theory maximum likelihood (Rhemtulla et al., 2012).
  Available only from raw data.

- Multivariate normality:

  Two tests of multivariate normality computed from the complete cases
  of the raw data: Mardia's (1970) multivariate skewness and kurtosis,
  and the Henze-Zirkler (1990) omnibus test. Both use the
  maximum-likelihood (divisor-\\n\\) covariance, so Mardia's
  coefficients differ from implementations using the unbiased divisor by
  a factor of \\(n / (n - 1))^3\\ for the skewness coefficient and \\(n
  / (n - 1))^2\\ for the kurtosis coefficient. The kurtosis coefficient
  is standardised with its exact null moments, \\p(p + 2)(n - 1) / (n +
  1)\\ for the mean (Mardia, 1970) and \\8p(p + 2)(n - 3)(n - p - 1)(n -
  p + 1) / ((n + 1)^2 (n + 3)(n + 5))\\ for the variance (Mardia, 1974).
  Both become the asymptotic pair \\p(p + 2)\\ and \\8p(p + 2) / n\\ as
  the sample size increases. Most other implementations use that
  asymptotic pair. It overstates both moments, by an amount that
  increases with the ratio of variables to observations. With 40
  variables and 200 observations it rejects data from an exact
  multivariate normal more than half of the time at the .05 level, where
  the exact moments hold the nominal rate. The skewness statistic always
  carries Mardia's (1974) correction, which makes its expectation equal
  its degrees of freedom at every sample size. Implementations that
  apply the correction only below 20 observations leave the statistic
  biased low, and the bias increases with the ratio of variables to
  observations until the test has no power. The correction sets the mean
  exactly, but the chi-square approximation keeps an inexact variance,
  so the skewness test holds a rejection rate near, but not exactly at,
  its nominal level. A small p-value indicates a departure from
  multivariate normality, a reason to prefer robust or ordinal
  estimation over normal-theory maximum likelihood. The Henze-Zirkler
  p-value comes from a lognormal approximation whose null variance falls
  geometrically with the number of variables. Beyond roughly 50 to 60
  variables (the exact point depends on the sample size) that variance
  is no longer resolvable in double precision, and the p-value would be
  exactly 0 or exactly 1 as a rounding-level residual decides, even for
  data from an exact multivariate normal. It is then withheld, with a
  warning, and the statistic alone is reported; the printed report says
  how many of the available tests reject multivariate normality.
  Available only from raw data, and skipped with a note if the
  complete-case covariance is singular.

- Outliers:

  Multivariate outliers flagged by their robust Mahalanobis distance. A
  high-breakdown robust location and scatter are estimated from the
  complete cases with the fast minimum covariance determinant (MCD)
  algorithm (Rousseeuw & Van Driessen, 1999), using a subset covering a
  proportion `mcd_alpha` of the observations; an observation whose
  squared robust distance exceeds `qchisq(outlier_cutoff, p)` is
  flagged. The search runs on columns divided by a robust scale and the
  estimate is returned to the supplied units, so the diagnostic does not
  depend on how the variables happen to be measured. The robust
  covariance is undefined with too few complete cases (\\n \le 2p\\),
  when the variables are so nearly linearly dependent that no covering
  subset has a usable covariance, and when a whole covering subset lies
  exactly on a lower-dimensional hyperplane (an *exact fit*, common with
  coarse discrete items on which many respondents answer an item pair
  identically). In all three cases the classical Mahalanobis distance is
  used instead, with a warning naming which of the three applies; those
  distances are computed from a covariance the outliers themselves
  inflate, so the diagnostic is no longer high-breakdown and tends to
  under-flag. If even the classical covariance is singular the
  diagnostic is skipped with a note. Available only from raw data.

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
#> These data are probably suitable for factor analysis (verbal bands: Kaiser &
#> Rice, 1974).
#> 
#> ✔ The Bartlett's test of sphericity was significant at an alpha level of .05.
#> These data are probably suitable for factor analysis.
#> χ²(153) = 2173.28, p < .001
#> 
#> ── Multicollinearity ───────────────────────────────────────────────────────────
#> 
#> ℹ Determinant: 0.0121. It falls as variables are added, so the condition index
#> below carries the verdict.
#> ✔ Condition number: 11.680 (condition index 3.418). No concern (index below 10;
#> Belsley, Kuh & Welsch, 1980).
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

# From raw data (N is taken from the data; the seed makes the outlier
# diagnostics reproducible)
efa_screen(GRiPS_raw, seed = 1)
#> Warning: A robust (MCD) covariance could not be computed; classical Mahalanobis
#> distances were used.
#> ℹ At least half the complete cases lie exactly on a lower-dimensional
#>   hyperplane (an "exact fit"). This is common with coarse discrete items, where
#>   many respondents give identical answers on an item pair; it does not mean the
#>   data are collinear at the correlation level.
#> 
#> ── Sampling adequacy and sphericity ────────────────────────────────────────────
#> 
#> ✔ The overall KMO value for your data is marvellous (Overall KMO = 0.955).
#> These data are probably suitable for factor analysis (verbal bands: Kaiser &
#> Rice, 1974).
#> 
#> ✔ The Bartlett's test of sphericity was significant at an alpha level of .05.
#> These data are probably suitable for factor analysis.
#> χ²(28) = 5054.06, p < .001
#> 
#> ── Multicollinearity ───────────────────────────────────────────────────────────
#> 
#> ℹ Determinant: 0.00188. It falls as variables are added, so the condition index
#> below carries the verdict.
#> ✔ Condition number: 24.340 (condition index 4.934). No concern (index below 10;
#> Belsley, Kuh & Welsch, 1980).
#> 
#> ── Per-variable diagnostics ────────────────────────────────────────────────────
#> 
#>           variance missing%   SMC   MSA  flags
#> fun          1.427        0 0.609 0.955       
#> friends      1.309        0 0.685 0.950       
#> enjoy        1.154        0 0.715 0.947 sparse
#> hurt         1.136        0 0.557 0.968 sparse
#> part         1.271        0 0.636 0.955       
#> commonly     1.099        0 0.653 0.958       
#> chances      1.389        0 0.592 0.961       
#> attracted    1.259        0 0.679 0.950 sparse
#> 
#> ── Multivariate normality ──────────────────────────────────────────────────────
#> 
#> ✖ Mardia's skewness: χ²(120) = 914.23, p < .001.
#> ✖ Mardia's kurtosis: z = 50.44, p < .001.
#> ✖ Henze-Zirkler: HZ = 11.65, p < .001.
#> These data depart from multivariate normality: 3 of the 3 tests reject it.
#> 
#> ── Outliers ────────────────────────────────────────────────────────────────────
#> 
#> ! A robust (MCD) covariance could not be computed; classical Mahalanobis
#> distances were used.
#> At least half the complete cases lie exactly on a lower-dimensional hyperplane
#> (an "exact fit"). This is common with coarse discrete items, where many
#> respondents give identical answers on an item pair; it does not mean the data
#> are collinear at the correlation level.
#> These distances come from a covariance the outliers themselves inflate, so the
#> diagnostic is no longer high-breakdown and tends to under-flag.
#> ℹ 71 of 810 observations were flagged as multivariate outliers (Mahalanobis
#> distance > 4.19).
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
#> ! 71 observations were flagged as potential multivariate outliers; inspect them
#>   (see `$outliers$flagged`) before down-weighting or excluding.
```
