# Screen data for exploratory factor analysis

Checks whether your data are suitable for exploratory factor analysis
(EFA). From a correlation matrix or raw data, it reports the
Kaiser-Meyer-Olkin (KMO) measure of sampling adequacy, Bartlett's test
of sphericity, the condition index of the correlation matrix, and each
variable's squared multiple correlation (SMC). When you supply raw data,
it also reports each variable's variance and percentage of missing
values, category counts for categorical variables, Mardia's tests of
multivariate normality, and multivariate outliers.

## Usage

``` r
efa_screen(
  x,
  N = NA,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
  mcd_alpha = 0.75,
  outlier_cutoff = 0.975,
  seed = lifecycle::deprecated()
)
```

## Source

Bartlett, M. S. (1951). The effect of standardization on a Chi-square
approximation in factor analysis. Biometrika, 38, 337-344.

Belsley, D. A. (1991). A guide to using the collinearity diagnostics.
Computer Science in Economics and Management, 4, 33-50.

Belsley, D. A., Kuh, E. & Welsch, R. E. (1980). Regression diagnostics:
Identifying influential data and sources of collinearity. Wiley.

Croux, C. & Haesbroeck, G. (1999). Influence function and efficiency of
the minimum covariance determinant scatter matrix estimator. Journal of
Multivariate Analysis, 71, 161-190.

Field, A. (2018). Discovering statistics using IBM SPSS statistics (5th
ed.). Sage.

Hubert, M., Debruyne, M. & Rousseeuw, P. J. (2018). Minimum covariance
determinant and extensions. WIREs Computational Statistics, 10, e1421.

Hubert, M., Rousseeuw, P. J. & Verdonck, T. (2012). A deterministic
algorithm for robust location and scatter. Journal of Computational and
Graphical Statistics, 21, 618-637.

Kaiser, H. F. (1970). A second generation little jiffy. Psychometrika,
35, 401-415.

Kaiser, H. F. & Rice, J. (1974). Little jiffy, mark IV. Educational and
Psychological Measurement, 34, 111-117.

Leys, C., Delacre, M., Mora, Y. L., Lakens, D. & Ley, C. (2019). How to
classify, detect, and manage univariate and multivariate outliers, with
emphasis on pre-registration. International Review of Social Psychology,
32, 5.

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

Zijlstra, W. P., van der Ark, L. A. & Sijtsma, K. (2007). Outlier
detection in test and questionnaire data. Multivariate Behavioral
Research, 42, 531-555.

## Arguments

- x:

  data.frame or matrix. Raw data, or a correlation matrix. Needs at
  least three variables, none of which is a perfect linear combination
  of the others.

- N:

  numeric. The number of observations. Set this only when you supply a
  correlation matrix; it is needed for Bartlett's test of sphericity and
  is taken from the data automatically when you supply raw data. Default
  is `NA`.

- use:

  character. How to handle missing values in raw data. For
  `cor_method = "pearson"`, `"spearman"`, or `"kendall"` this is passed
  to [`stats::cor()`](https://rdrr.io/r/stats/cor.html). For `"poly"` or
  `"tetra"` the same rule is applied to the raw data before the
  correlations are estimated; `"all.obs"` and `"everything"` then stop
  with an error on any missing value, instead of returning `NA`
  correlations. Default is `"pairwise.complete.obs"`.

- cor_method:

  character. How to compute correlations from raw data: `"pearson"`,
  `"spearman"`, or `"kendall"` (via
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)), or `"poly"` /
  `"tetra"` for polychoric / tetrachoric correlations of ordinal or
  binary data. Spearman and Kendall matrices are screened as they are,
  without conversion to the Pearson scale. Default is `"pearson"`.

- mcd_alpha:

  numeric. The proportion of cases used to build the robust outlier
  estimate, between 0.5 and 1. `0.5` resists the most outliers (up to
  half the cases). The default, `0.75`, resists up to a quarter but is
  more efficient (Hubert et al., 2018), and Leys et al. (2019) recommend
  it. Used only with raw data in which all variables are continuous.

- outlier_cutoff:

  numeric. The probability used to set the cutoff for flagging a
  multivariate outlier, between 0.5 and 0.9999. Default is `0.975`. Used
  only with raw data in which all variables are continuous.

- seed:

  **\[deprecated\]** Accepted and ignored. The outlier detection is
  deterministic, so it needs no seed.

## Value

An object of class `efa_screen`, a list containing:

- kmo:

  A list with the overall KMO (`KMO`) and the per-variable KMO
  (`KMO_i`).

- bartlett:

  A list with Bartlett's chi-square statistic (`chisq`), its `p_value`,
  and its degrees of freedom (`df`); `chisq` and `p_value` are `NA` when
  `N` was too small for the correction. `NULL` when `N` is unavailable.

- determinant:

  The determinant of the correlation matrix.

- condition:

  The condition number of the correlation matrix (largest eigenvalue
  over smallest). Its square root is the condition index.

- collinear:

  The names of the variables involved in a near linear dependency (see
  Details). `NULL` when no condition index is above 30.

- smc:

  The per-variable squared multiple correlations.

- per_item:

  A data frame with one row per variable (row names are the variable
  names): `variance`, `missing` (percentage), `smc`, `kmo_i`, and
  `flags` (any sparse/empty-category issues). `NULL` when a correlation
  matrix is supplied instead of raw data.

- normality:

  A list with `mardia` (skewness `skewness`, `skewness_df`,
  `skewness_p`, kurtosis `kurtosis` and `kurtosis_p`, and the underlying
  `b1p`/`b2p`) and `n_complete` (the number of complete cases used).
  `NULL` without raw data, or a note (`$message`) that tells why the
  tests were skipped.

- outliers:

  A list with `distances` (each complete case's distance, named by its
  row number), `cutoff` (the flagging threshold, on the same scale as
  `distances`), `flagged` (the row numbers above `cutoff`), `center` and
  `cov` (the location and scatter used), `method`, `fallback_reason`,
  and `n_complete`. `method` is `"mcd"`, `"classical"` (the fallback for
  continuous data; `fallback_reason` tells why), or `"ranking"`
  (categorical data: classical distances, `cutoff` is `NA`, and
  `flagged` is empty). `NULL` without raw data, or a note (`$message`)
  that tells why no distances could be computed.

- categories:

  A named list with the response-category counts for each categorical
  variable (in category order); `NA` for a variable treated as
  continuous. `NULL` without raw data.

- note:

  A note (`$message`) that tells why the raw-data diagnostics
  (`per_item`, `normality`, `outliers`, `categories`) are missing, when
  a correlation matrix is supplied instead of raw data. `NULL` when raw
  data are supplied.

- settings:

  The settings used: `N`, `n_obs` (rows in the raw data supplied, `NA`
  for a correlation-matrix input), `use`, `cor_method`, `mcd_alpha`, and
  `outlier_cutoff`.

## Details

The diagnostics are computed from the analysis correlation matrix \\R\\:

- KMO:

  The Kaiser-Meyer-Olkin measure of sampling adequacy (Kaiser, 1970;
  Kaiser & Rice, 1974), overall and for each variable; see
  [`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md).
  It shows how much common variance your variables share. Higher values
  are better; a common rule of thumb treats values below .50 as
  unacceptable.

- Bartlett:

  Bartlett's (1951) test of sphericity: the likelihood-ratio test of
  whether the correlation matrix is an identity matrix, i.e., whether
  your variables correlate with each other at all; see
  [`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md).
  At usual sample sizes the test is almost always significant, so a
  significant result is weak evidence of suitability, but a
  non-significant result is a strong warning. The test needs the sample
  size `N`; without it, `$bartlett` is `NULL`. If `N` is too small
  relative to the number of variables, the statistic is `NA`.

- Determinant:

  The determinant of \\R\\, returned but not used for a verdict. It
  falls as you add variables, also when they are not collinear, so a
  fixed cut-off on it (such as the 0.00001 often quoted from
  Field, 2018) mostly reflects the number of variables.

- Condition index:

  The square root of the condition number (the ratio of the largest to
  the smallest eigenvalue of \\R\\): the collinearity diagnostic of
  Belsley, Kuh & Welsch (1980). An index of 10 or less is rarely of
  interest; an index above 30 flags a near linear dependency (Belsley,
  1991). The report grades such an index as moderate (30 to 100), strong
  (100 to 300), or very strong (above 300), and names the variables
  involved: those whose variance-decomposition proportions, summed over
  all indices above 30, are above .5 (Belsley, Kuh & Welsch, 1980).
  These values come from regression on data that are not centred, but a
  correlation matrix is centred, so use them as a guide and not as a
  test.

- SMC:

  The squared multiple correlation of each variable with all the others.
  A low value flags a variable that has little in common with the rest
  of your set.

- Variance and missing data:

  For raw data: each variable's variance (over its available values) and
  percentage of missing values over all rows you supplied. Factor
  columns are recoded to integer levels first, so `variance` reflects
  those codes. With `use = "complete.obs"` or `"na.or.complete"`, `N` is
  the number of complete rows; otherwise it is the number of rows.

- Categories:

  For raw data: for each variable with fewer than ten distinct values
  (treated as categorical), the count of responses in each category. A
  category with fewer than five responses is flagged as sparse, and an
  unused category between the smallest and largest response is flagged
  as empty. As a rule of thumb, items with fewer than five response
  categories are better analysed with `cor_method = "poly"` or `"tetra"`
  than with Pearson correlations (Rhemtulla et al., 2012).

- Multivariate normality:

  For raw data, using only complete cases: Mardia's (1970) tests of
  multivariate skewness and kurtosis, with the small-sample corrections
  of Mardia (1974). A small p-value suggests your data depart from a
  multivariate normal distribution, a reason to prefer robust standard
  errors over normal-theory ones. In a very small sample the kurtosis
  statistic is `NA`. With `cor_method = "poly"` or `"tetra"` the tests
  are skipped: they test the observed codes, but a polychoric analysis
  assumes normality only of the underlying continuous variables.

- Outliers:

  For raw data, using only complete cases: each case's distance from the
  centre of your data.

  - If all variables are continuous, the distance is a robust one: it
    uses the reweighted minimum covariance determinant (MCD; Rousseeuw &
    Van Driessen, 1999), found with the deterministic DetMCD algorithm
    (Hubert, Rousseeuw & Verdonck, 2012) and corrected for consistency
    (Croux & Haesbroeck, 1999) and small samples (Pison et al., 2002). A
    case is flagged when its squared distance is above the
    `outlier_cutoff` quantile of the chi-square distribution. When there
    are too few complete cases or the variables are nearly collinear,
    the classical (non-robust) distance is used instead, and
    `$outliers$fallback_reason` tells why.

  - If any variable is categorical (fewer than ten distinct values, e.g.
    a rating-scale item), no case is flagged, and the report names the
    cases with the largest classical distances. The MCD theory needs
    continuous data (Hubert, Debruyne & Rousseeuw, 2018): with tied
    answers, many cases can share a hyperplane, which makes the MCD
    degenerate. Also, the distances of rating-scale data do not follow
    the chi-square distribution, so a chi-square cutoff flags many clean
    cases when items are skewed (Zijlstra, van der Ark & Sijtsma, 2007).
    To screen for careless responding, use indices made for that
    purpose.

  Distance-based screening does not detect straight-lining (the same
  answer to all items).

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
#> ✔ Overall KMO = 0.916: marvellous (Kaiser & Rice, 1974).
#> ✔ Bartlett's test of sphericity: χ²(153) = 2173.28, p < .001 (significant).
#> 
#> ── Multicollinearity ───────────────────────────────────────────────────────────
#> 
#> ✔ Condition index: 3.418. An index of 10 or less is rarely of interest
#> (Belsley, 1991).
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
#> ✔ Overall KMO = 0.955: marvellous (Kaiser & Rice, 1974).
#> ✔ Bartlett's test of sphericity: χ²(28) = 5054.06, p < .001 (significant).
#> 
#> ── Multicollinearity ───────────────────────────────────────────────────────────
#> 
#> ✔ Condition index: 4.934. An index of 10 or less is rarely of interest
#> (Belsley, 1991).
#> 
#> ── Per-variable diagnostics ────────────────────────────────────────────────────
#> 
#>           variance   SMC   MSA  flags
#> fun          1.427 0.609 0.955       
#> friends      1.309 0.685 0.950       
#> enjoy        1.154 0.715 0.947 sparse
#> hurt         1.136 0.557 0.968 sparse
#> part         1.271 0.636 0.955       
#> commonly     1.099 0.653 0.958       
#> chances      1.389 0.592 0.961       
#> attracted    1.259 0.679 0.950 sparse
#> 
#> ── Multivariate normality ──────────────────────────────────────────────────────
#> 
#> ✖ Mardia's skewness: χ²(120) = 914.23, p < .001.
#> ✖ Mardia's kurtosis: z = 50.44, p < .001.
#> These data depart from multivariate normality: both tests reject it.
#> 
#> ── Outliers ────────────────────────────────────────────────────────────────────
#> 
#> ℹ No case is flagged: some variables have fewer than ten distinct values
#> (categorical), and a chi-square cutoff is not valid for their Mahalanobis
#> distances. The largest distances are in rows 561, 238, 622, 549, and 186.
#> 
#> ── Recommendations ─────────────────────────────────────────────────────────────
#> 
#> ! These data depart from multivariate normality; normal-theory standard errors
#>   and fit statistics may be biased - prefer robust (sandwich) or bootstrapped
#>   standard errors.
#> ! 3 variables have a sparse response category (< 5 responses): enjoy, hurt, and
#>   attracted; a low-frequency category can destabilise polychoric estimates -
#>   consider collapsing it into an adjacent category.
```
