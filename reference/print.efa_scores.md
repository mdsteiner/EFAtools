# Print and format an efa_scores object

[`print()`](https://rdrr.io/r/base/print.html) shows a concise overview
of an
[`efa_scores()`](https://mdsteiner.github.io/EFAtools/reference/efa_scores.md)
result: a header naming the method and whether factor scores were
computed, and the per-factor determinacy table (determinacy, squared
determinacy, and Guttman index).
[`summary()`](https://rdrr.io/r/base/summary.html) returns a
`summary.efa_scores` object whose print method adds the full
factor-weight matrix, the score validity/univocality matrix, and the
score intercorrelations.
[`format()`](https://rdrr.io/r/base/format.html) assembles the same
report and returns it as a character vector;
[`print()`](https://rdrr.io/r/base/print.html) is
`cat(format(x), sep = "\n")`. The lines follow the active console theme,
so they are plain when colours are disabled (for example when captured
into a file or stripped with
[`cli::ansi_strip()`](https://cli.r-lib.org/reference/ansi_strip.html)).

## Usage

``` r
# S3 method for class 'efa_scores'
print(x, digits = 3, ...)

# S3 method for class 'efa_scores'
format(x, digits = 3, ...)

# S3 method for class 'efa_scores'
summary(object, digits = 3, ...)

# S3 method for class 'summary.efa_scores'
print(x, ...)

# S3 method for class 'summary.efa_scores'
format(x, digits = x$opts$digits, ...)
```

## Arguments

- x, object:

  An object of class `efa_scores`; for the `summary.efa_scores` methods,
  the object returned by
  [`summary()`](https://rdrr.io/r/base/summary.html).

- digits:

  numeric. Number of decimal places for the printed tables. Default is
  3.

- ...:

  Not used; for consistency with the generics.

## Value

[`print()`](https://rdrr.io/r/base/print.html) and the print method for
`summary.efa_scores` objects return their argument invisibly.
[`format()`](https://rdrr.io/r/base/format.html) returns a character
vector with the report lines (styled to the active console theme; plain
when colours are disabled).
[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.efa_scores`.

## See also

Other factor scoring:
[`efa_scores()`](https://mdsteiner.github.io/EFAtools/reference/efa_scores.md)

## Examples

``` r
efa <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
               estimator = "PAF", rotation = "oblimin")
fs <- efa_scores(test_models$baseline$cormat, f = efa)
#> ℹ `x` is a correlation matrix; factor scores cannot be computed. Only factor
#>   weights and score diagnostics are returned. Enter raw data to get factor
#>   scores.
fs
#> 
#> ── Factor scores (regression) ──────────────────────────────────────────────────
#> 
#> Weights and diagnostics only (correlation-matrix input; no scores).
#> 
#> ── Score determinacy ───────────────────────────────────────────────────────────
#> 
#>      rho  rho2 guttman
#> F1 0.894 0.798   0.597
#> F2 0.888 0.788   0.576
#> F3 0.883 0.780   0.561
summary(fs)
#> 
#> ── Factor scores (regression) ──────────────────────────────────────────────────
#> 
#> Weights and diagnostics only (correlation-matrix input; no scores).
#> 
#> ── Score determinacy ───────────────────────────────────────────────────────────
#> 
#>      rho  rho2 guttman
#> F1 0.894 0.798   0.597
#> F2 0.888 0.788   0.576
#> F3 0.883 0.780   0.561
#> 
#> ── Factor weights ──────────────────────────────────────────────────────────────
#> 
#>        F1    F2     F3
#> V1  0.016 0.037  0.206
#> V2  0.023 0.036  0.146
#> V3  0.038 0.033  0.140
#> V4  0.060 0.025  0.194
#> V5  0.062 0.014  0.138
#> V6  0.009 0.013  0.247
#> V7  0.024 0.177  0.053
#> V8  0.017 0.187  0.031
#> V9  0.031 0.173  0.020
#> V10 0.016 0.222 -0.002
#> V11 0.026 0.115  0.084
#> V12 0.035 0.240  0.030
#> V13 0.202 0.051  0.010
#> V14 0.163 0.002  0.050
#> V15 0.177 0.059  0.007
#> V16 0.170 0.006  0.050
#> V17 0.214 0.013  0.016
#> V18 0.173 0.025  0.041
#> 
#> ── Score validity and univocality ──────────────────────────────────────────────
#> 
#> Diagonal: validity (score-factor correlation). Off-diagonal: univocality.
#> 
#>       F1    F2    F3
#> F1 0.894 0.638 0.668
#> F2 0.643 0.888 0.650
#> F3 0.676 0.653 0.883
#> 
#> ── Score intercorrelations ─────────────────────────────────────────────────────
#> 
#>       F1    F2    F3
#> F1 1.000 0.719 0.756
#> F2 0.719 1.000 0.735
#> F3 0.756 0.735 1.000

# format() returns the same lines as plain text:
writeLines(format(fs))
#> 
#> ── Factor scores (regression) ──────────────────────────────────────────────────
#> 
#> Weights and diagnostics only (correlation-matrix input; no scores).
#> 
#> ── Score determinacy ───────────────────────────────────────────────────────────
#> 
#>      rho  rho2 guttman
#> F1 0.894 0.798   0.597
#> F2 0.888 0.788   0.576
#> F3 0.883 0.780   0.561
```
