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
vector with the report lines.
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
#>      rho  rho2  guttman
#> F1  .894  .798     .597
#> F2  .888  .788     .576
#> F3  .883  .780     .561
summary(fs)
#> 
#> ── Factor scores (regression) ──────────────────────────────────────────────────
#> 
#> Weights and diagnostics only (correlation-matrix input; no scores).
#> 
#> ── Score determinacy ───────────────────────────────────────────────────────────
#> 
#>      rho  rho2  guttman
#> F1  .894  .798     .597
#> F2  .888  .788     .576
#> F3  .883  .780     .561
#> 
#> ── Factor weights ──────────────────────────────────────────────────────────────
#> 
#>       F1    F2     F3
#> V1   .016  .037   .206
#> V2   .023  .036   .146
#> V3   .038  .033   .140
#> V4   .060  .025   .194
#> V5   .062  .014   .138
#> V6   .009  .013   .247
#> V7   .024  .177   .053
#> V8   .017  .187   .031
#> V9   .031  .173   .020
#> V10  .016  .222  -.002
#> V11  .026  .115   .084
#> V12  .035  .240   .030
#> V13  .202  .051   .010
#> V14  .163  .002   .050
#> V15  .177  .059   .007
#> V16  .170  .006   .050
#> V17  .214  .013   .016
#> V18  .173  .025   .041
#> 
#> ── Score validity and univocality ──────────────────────────────────────────────
#> 
#> Diagonal: validity (score-factor correlation). Off-diagonal: univocality.
#> 
#>      F1    F2    F3
#> F1  .894  .638  .668
#> F2  .643  .888  .650
#> F3  .676  .653  .883
#> 
#> ── Score intercorrelations ─────────────────────────────────────────────────────
#> 
#>       F1     F2     F3
#> F1  1.000   .719   .756
#> F2   .719  1.000   .735
#> F3   .756   .735  1.000

# format() returns the same lines as a character vector:
writeLines(format(fs))
#> 
#> ── Factor scores (regression) ──────────────────────────────────────────────────
#> 
#> Weights and diagnostics only (correlation-matrix input; no scores).
#> 
#> ── Score determinacy ───────────────────────────────────────────────────────────
#> 
#>      rho  rho2  guttman
#> F1  .894  .798     .597
#> F2  .888  .788     .576
#> F3  .883  .780     .561
```
