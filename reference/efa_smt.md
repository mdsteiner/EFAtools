# Sequential chi square model tests, RMSEA lower bound, and AIC

Sequential chi square model tests (SMT) are a factor retention method
where multiple EFAs with increasing numbers of factors are fitted and
the number of factors for which the Chi Square value first becomes
non-significant is taken as the suggested number of factors. Preacher,
Zhang, Kim, & Mels (2013) suggested a similar approach with the lower
bound of the 90% confidence interval of the Root Mean Square Error of
Approximation (RMSEA; Browne & Cudeck, 1992; Steiger & Lind, 1980), and
with the Akaike Information Criterion (AIC). For the RMSEA, the number
of factors for which this lower bound first falls below .05 is the
suggested number of factors to retain. For the AIC, it is the number of
factors where the AIC is lowest.

## Usage

``` r
efa_smt(
  x,
  N = NA,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
  estimate_control = NULL
)
```

## Source

Auerswald, M., & Moshagen, M. (2019). How to determine the number of
factors to retain in exploratory factor analysis: A comparison of
extraction methods under realistic conditions. Psychological Methods,
24(4), 468–491. https://doi.org/10.1037/met0000200

Browne, M.W., & Cudeck, R. (1992). Alternative ways of assessing model
fit. Sociological Methods and Research, 21, 230–258.

Preacher, K. J., Zhang G., Kim, C., & Mels, G. (2013). Choosing the
Optimal Number of Factors in Exploratory Factor Analysis: A Model
Selection Perspective, Multivariate Behavioral Research, 48(1), 28-56,
doi:10.1080/00273171.2012.710386

Steiger, J. H., & Lind, J. C. (1980, May). Statistically based tests for
the number of common factors. Paper presented at the annual meeting of
the Psychometric Society, Iowa City, IA.

## Arguments

- x:

  data.frame or matrix. Dataframe or matrix of raw data or matrix with
  correlations.

- N:

  numeric. The number of observations. Needs only be specified if a
  correlation matrix is used.

- use:

  character. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) if raw data is
  given as input. Default is "pairwise.complete.obs".

- cor_method:

  character. One of `"pearson"`, `"spearman"`, or `"kendall"`, passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html). `"poly"` and
  `"tetra"` are not supported because `SMT` rests on a normal-theory
  chi-square test that is not valid for polychoric / tetrachoric
  correlations. Default is `"pearson"`.

- estimate_control:

  an
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  object with the estimation settings for the sequential
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  fits. `NULL` (default) uses the
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  defaults. The sequential models are fitted with maximum likelihood
  (the chi-square, RMSEA, and AIC the SMT is built on are defined for
  it), so of the estimation knobs only `start_method` takes effect; the
  ones governing principal axis factoring do not apply. The models are
  unrotated, so no rotation settings apply either.

## Value

An object of class `efa_retention` (see
[`print.efa_retention()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_retention.md)
for the print method). SMT has no plot;
[`plot.efa_retention()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_retention.md)
returns `NULL` with a message for it. Its main fields are:

- n_factors:

  A named numeric vector (`"chi"`, `"RMSEA"`, `"AIC"`) with the
  suggested number of factors from the sequential chi-square model
  tests, the RMSEA lower bound, and the AIC.

- results:

  A list with one record per criterion, each holding the criterion
  values for the null model (zero factors) through the maximum number of
  factors.

- settings:

  A list of the settings used.

## Details

As a first step in the procedure, a maximum number of factors to extract
is determined for which the model is still over-identified (df \> 0).

Then, EFAs with increasing numbers of factors from 1 to the maximum
number are fitted with maximum likelihood estimation.

For the SMT, first the significance of the chi square value for a model
with 0 factors is determined. If this value is not significant, 0
factors are suggested to retain. If it is significant, a model with 1
factor is estimated and the significance of its chi square value is
determined, and so on, until a non-significant result is obtained. The
suggested number of factors is the number of factors for the model where
the chi square value first becomes non-significant.

Regarding the RMSEA, the suggested number of factors is the number of
factors for the model where the lower bound of the 90% confidence
interval of the RMSEA first falls below the .05 threshold.

Regarding the AIC, the suggested number of factors is the number of
factors for the model with the lowest AIC.

The sequential models are fitted without inequality constraints, so a
solution can be inadmissible (a Heywood case, or a fit that did not
converge). Only the models the three rules actually select are checked
for this; if one of them is inadmissible a warning is raised and the
corresponding suggestion should be interpreted with caution.

In comparison with other prominent factor retention criteria, SMT
performed well at determining the number of factors to extract in EFA
(Auerswald & Moshagen, 2019). The RMSEA lower bound also performed well
at determining the true number of factors, while the AIC performed well
at determining the most generalizable model (Preacher, Zhang, Kim, &
Mels, 2013).

## See also

[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
as a wrapper function for this and the other factor retention criteria.

Other factor retention criteria:
[`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md),
[`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md),
[`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md),
[`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
[`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md),
[`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md),
[`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md),
[`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md)

## Examples

``` r
SMT_base <- efa_smt(test_models$baseline$cormat, N = 500)
SMT_base
#> ── Sequential model tests ──────────────────────────────────────────────────────
#> 
#> • Sequential chi-square model tests: 3
#> • Lower bound of RMSEA 90% CI: 2
#> • Akaike Information Criterion: 3
```
