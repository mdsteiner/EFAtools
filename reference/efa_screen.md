# Screen data for exploratory factor analysis

Checks whether your data are suitable for exploratory factor analysis
(EFA). From a correlation matrix or raw data, it reports the
Kaiser-Meyer-Olkin (KMO) measure of sampling adequacy, Bartlett's test
of sphericity, the determinant and condition number of the correlation
matrix, and each variable's squared multiple correlation (SMC). When you
supply raw data, it also reports each variable's variance and percentage
of missing values, category counts for categorical variables, tests of
multivariate normality, and multivariate outliers.

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

Belsley, D. A. (1991). A guide to using the collinearity diagnostics.
Computer Science in Economics and Management, 4, 33-50.

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
  binary data. A Spearman or Kendall correlation matrix is screened on
  its own scale, not converted to look like a Pearson correlation
  matrix; Kendall's tau in particular measures something different from
  a Pearson correlation, not just a rescaled version of it. Default is
  `"pearson"`.

- mcd_alpha:

  numeric. The proportion of cases used to build the robust outlier
  estimate, between 0.5 and 1. The default, `0.5`, is the most robust
  choice; a larger value uses more of the data but resists outliers less
  well. Used only with raw data.

- outlier_cutoff:

  numeric. The probability used to set the cutoff for flagging a
  multivariate outlier, between 0.5 and 0.9999. Default is `0.975`. Used
  only with raw data.

- seed:

  integer. A seed for the random subsets used by the outlier detection,
  so the result is reproducible. Does not affect your random-number
  generator elsewhere. Default is `NULL`. Used only with raw data.

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
  over smallest).

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
  `b1p`/`b2p`), `hz` (the Henze-Zirkler `statistic` and its `p_value`),
  and `n_complete` (the number of complete cases used). `NULL` without
  raw data, or a note explaining why when the complete-case data cannot
  support the tests.

- outliers:

  A list with `distances` (each complete case's robust distance, named
  by its row number), `cutoff` (the flagging threshold, on the same
  scale as `distances`), `flagged` (the row numbers exceeding `cutoff`),
  `center` and `cov` (the robust location and scatter), `method`
  (`"mcd"` or the `"classical"` fallback), `fallback_reason` (why the
  robust estimate was unavailable, when `method` is `"classical"`), and
  `n_complete`. `NULL` without raw data, or a note explaining why when
  no covariance can be formed.

- categories:

  A named list with the response-category counts for each categorical
  variable (in category order); `NA` for a variable treated as
  continuous. `NULL` without raw data.

- note:

  Explains why the raw-data diagnostics (`per_item`, `normality`,
  `outliers`, `categories`) are missing, when a correlation matrix is
  supplied instead of raw data. `NULL` when raw data are supplied.

- settings:

  The settings used: `N`, `n_obs` (rows in the raw data supplied, `NA`
  for a correlation-matrix input), `use`, `cor_method`, `mcd_alpha`,
  `outlier_cutoff`, and `seed`.

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
  A significant result supports doing a factor analysis. The test needs
  the sample size `N`; without it, this diagnostic is skipped with a
  warning and `$bartlett` is `NULL`. If `N` is too small relative to the
  number of variables, the statistic is `NA`, also with a warning.

- Determinant:

  The determinant of \\R\\, reported as a number only. It falls as you
  add variables even when the variables are not collinear, so a fixed
  cut-off on it (such as the 0.00001 often quoted from Field, 2018) says
  more about how many variables you have than about your data. Use the
  condition number instead.

- Condition number:

  The ratio of the largest to the smallest eigenvalue of \\R\\. Its
  square root, the condition index, is the collinearity diagnostic of
  Belsley, Kuh & Welsch (1980); it drives the printed report and its
  recommendation. An index of 10 or less is rarely of interest. An index
  above 30 flags a near linear dependency: two or more variables that
  together carry much the same information. An index between the two is
  not negligible, but it stays below the value that flags a dependency.
  Belsley (1991) gives 30 as one example value and calls the choice of a
  cut-off "somewhat of an art form", so the report grades an index above
  30 by its position on the scale 1, 3, 10, 30, 100, 300, 1000: moderate
  (30 to 100), strong (100 to 300), or very strong (above 300). These
  values come from regression diagnostics on data that are not centred,
  but a correlation matrix is centred, so use them as a guide and not as
  a test.

- SMC:

  The squared multiple correlation of each variable with all the others.
  A low value flags a variable that has little in common with the rest
  of your set.

- Variance and missing data:

  For raw data: each variable's variance (over its available values) and
  percentage of missing values, computed from every row you supplied.
  These missing-value percentages explain why the correlation matrix's
  sample size (`N`) can be smaller than the number of rows in your data.
  Ordered-factor columns are recoded to integer levels first, so
  `variance` reflects those codes.

- Categories:

  For raw data: for each variable with fewer than ten distinct values
  (treated as categorical), the count of responses in each category. A
  category with fewer than five responses is flagged as sparse, and an
  unused category between the smallest and largest response is flagged
  as empty. As a rule of thumb, items with fewer than five response
  categories are better analysed with `cor_method = "poly"` or `"tetra"`
  than with Pearson correlations (Rhemtulla et al., 2012).

- Multivariate normality:

  For raw data, using only complete cases: two tests of multivariate
  normality, Mardia's (1970) test of skewness and kurtosis and the
  Henze-Zirkler (1990) test. A small p-value on either test suggests
  your data depart from a multivariate normal distribution, a reason to
  prefer a robust or ordinal method over normal-theory maximum
  likelihood. In a very small sample the kurtosis statistic is `NA`. The
  Henze-Zirkler p-value is not available with more than about 50 to 60
  variables; its test statistic is still reported.

- Outliers:

  For raw data, using only complete cases: multivariate outliers, found
  from a robust estimate of each case's distance from the centre of your
  data (the minimum covariance determinant method; Rousseeuw & Van
  Driessen, 1999). A flagged case is unusually far from the rest of your
  sample. When there are too few complete cases, the variables are too
  collinear, or too many cases share identical answers, a plain
  (non-robust) distance is used instead, with a warning explaining why.

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
#> ✔ Condition number: 11.680 (condition index 3.418). An index of 10 or less is
#> rarely of interest (Belsley, 1991).
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
#> ✔ Condition number: 24.340 (condition index 4.934). An index of 10 or less is
#> rarely of interest (Belsley, 1991).
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
