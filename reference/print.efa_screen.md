# Print and format an efa_screen object

[`print()`](https://rdrr.io/r/base/print.html) turns the factor-analysis
screening diagnostics computed by
[`efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/efa_screen.md)
into a sectioned report with banded, colour-coded verdicts: sampling
adequacy and sphericity (the Kaiser-Meyer-Olkin measure and Bartlett's
test of sphericity), multicollinearity (the determinant and condition
number of the correlation matrix), the per-variable diagnostics, and,
when raw data were supplied, multivariate normality and multivariate
outliers. It closes with a consolidated list of actionable
recommendations (for example, which items to consider dropping, whether
to prefer an ordinal or a robust estimator, and a caveat that keeps an
over-powered Bartlett's test from being over-trusted).
[`format()`](https://rdrr.io/r/base/format.html) assembles the same
report and returns it as a character vector;
[`print()`](https://rdrr.io/r/base/print.html) is
`cat(format(x), sep = "\n")`. The lines follow the active console theme,
so they are plain when colours are disabled (for example when captured
into a file or stripped with
[`cli::ansi_strip()`](https://cli.r-lib.org/reference/ansi_strip.html)).
[`print()`](https://rdrr.io/r/base/print.html) does not draw a plot.

## Usage

``` r
# S3 method for class 'efa_screen'
print(x, digits = 3, ...)

# S3 method for class 'efa_screen'
format(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `efa_screen` (output from
  [`efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/efa_screen.md)).

- digits:

  Integer. The number of decimal places the reported values are rounded
  to. Default is 3.

- ...:

  Not used; for consistency with the generic.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its argument `x`
invisibly. [`format()`](https://rdrr.io/r/base/format.html) returns a
character vector with the report lines.

## See also

Other factor analysis suitability:
[`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md),
[`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md),
[`efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/efa_screen.md)

## Examples

``` r
# From raw data
efa_screen(iris[, 1:4])
#> 
#> ── Sampling adequacy and sphericity ────────────────────────────────────────────
#> 
#> ✖ The overall KMO value for your data is miserable (Overall KMO = 0.54).
#> These data are hardly suitable for factor analysis.
#> 
#> ✔ The Bartlett's test of sphericity was significant at an alpha level of .05.
#> These data are probably suitable for factor analysis.
#> χ²(6) = 706.96, p < .001
#> 
#> ── Multicollinearity ───────────────────────────────────────────────────────────
#> 
#> ✔ Determinant: 0.00811. No concern (a value near 0 signals multicollinearity).
#> ! Condition number: 140.889 (condition index 11.870). Moderate
#> multicollinearity (index 10 to 30; Belsley, Kuh & Welsch, 1980).
#> 
#> ── Per-variable diagnostics ────────────────────────────────────────────────────
#> 
#>              variance missing%   SMC   MSA
#> Sepal.Length    0.686        0 0.859 0.584
#> Sepal.Width     0.190        0 0.524 0.270
#> Petal.Length    3.116        0 0.968 0.531
#> Petal.Width     0.581        0 0.938 0.634
#> 
#> ── Multivariate normality ──────────────────────────────────────────────────────
#> 
#> ✖ Mardia's skewness: χ²(20) = 67.43, p < .001.
#> ✔ Mardia's kurtosis: z = -0.23, p = 0.818.
#> ✖ Henze-Zirkler: HZ = 2.34, p < .001.
#> These data depart from multivariate normality.
#> 
#> ── Outliers ────────────────────────────────────────────────────────────────────
#> 
#> ℹ 55 of 150 observations were flagged as multivariate outliers (robust distance
#> > 3.34).
#> 
#> ── Recommendations ─────────────────────────────────────────────────────────────
#> 
#> ! Overall sampling adequacy is miserable (KMO = 0.54); the variables may not
#>   form a factorable set.
#> ! 1 variable has a low individual MSA (< .5): Sepal.Width; consider removing it
#>   (little shared variance).
#> ! These data depart from multivariate normality; normal-theory standard errors
#>   and fit statistics may be biased - prefer robust (sandwich) or bootstrapped
#>   standard errors.
#> ! Bartlett's test is significant, but it assumes multivariate normality and
#>   grows more sensitive as N increases; because these data are non-normal, treat
#>   it as uninformative here and rely on the KMO.
#> ! 55 of 150 observations (37%) exceed the outlier cutoff, far above the 2.5%
#>   expected under multivariate normality; this usually means the data are not
#>   elliptically distributed (subgroups or a mixture) rather than that this many
#>   cases are contaminated.

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
#> χ²(153) = 2173.28, p < .001
#> 
#> ── Multicollinearity ───────────────────────────────────────────────────────────
#> 
#> ✔ Determinant: 0.0121. No concern (a value near 0 signals multicollinearity).
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

# format() returns the same lines as a character vector:
writeLines(format(efa_screen(test_models$baseline$cormat, N = 500)))
#> 
#> ── Sampling adequacy and sphericity ────────────────────────────────────────────
#> 
#> ✔ The overall KMO value for your data is marvellous (Overall KMO = 0.916).
#> These data are probably suitable for factor analysis.
#> 
#> ✔ The Bartlett's test of sphericity was significant at an alpha level of .05.
#> These data are probably suitable for factor analysis.
#> χ²(153) = 2173.28, p < .001
#> 
#> ── Multicollinearity ───────────────────────────────────────────────────────────
#> 
#> ✔ Determinant: 0.0121. No concern (a value near 0 signals multicollinearity).
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
```
