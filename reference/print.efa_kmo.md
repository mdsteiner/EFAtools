# Print and format an efa_kmo object

[`print()`](https://rdrr.io/r/base/print.html) shows the
Kaiser-Meyer-Olkin (KMO) criterion computed by
[`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md):
a titled section with a verdict on the overall KMO value (and what it
implies for the suitability of the data for factor analysis), the
overall value, and the per-variable KMO values.
[`format()`](https://rdrr.io/r/base/format.html) assembles the same
report and returns it as a character vector;
[`print()`](https://rdrr.io/r/base/print.html) is
`cat(format(x), sep = "\n")`. The lines follow the active console theme,
so they are plain when colours are disabled (for example when captured
into a file or stripped with
[`cli::ansi_strip()`](https://cli.r-lib.org/reference/ansi_strip.html)).

## Usage

``` r
# S3 method for class 'efa_kmo'
print(x, ...)

# S3 method for class 'efa_kmo'
format(x, ...)
```

## Arguments

- x:

  An object of class `efa_kmo` (output from
  [`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md)).

- ...:

  Not used; for consistency with the generic.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its argument `x`
invisibly. [`format()`](https://rdrr.io/r/base/format.html) returns a
character vector with the report lines.

## Examples

``` r
KMO_base <- efa_kmo(test_models$baseline$cormat)
KMO_base
#> 
#> ── Kaiser-Meyer-Olkin criterion (KMO) ──────────────────────────────────────────
#> 
#> ✔ The overall KMO value for your data is marvellous.
#> These data are probably suitable for factor analysis.
#> 
#> Overall: 0.916
#> 
#> For each variable:
#>    V1    V2    V3    V4    V5    V6    V7    V8    V9   V10   V11   V12   V13 
#> 0.900 0.914 0.924 0.932 0.923 0.891 0.928 0.919 0.916 0.892 0.928 0.908 0.922 
#>   V14   V15   V16   V17   V18 
#> 0.905 0.924 0.934 0.907 0.923 

# format() returns the same lines as a character vector:
writeLines(format(KMO_base))
#> 
#> ── Kaiser-Meyer-Olkin criterion (KMO) ──────────────────────────────────────────
#> 
#> ✔ The overall KMO value for your data is marvellous.
#> These data are probably suitable for factor analysis.
#> 
#> Overall: 0.916
#> 
#> For each variable:
#>    V1    V2    V3    V4    V5    V6    V7    V8    V9   V10   V11   V12   V13 
#> 0.900 0.914 0.924 0.932 0.923 0.891 0.928 0.919 0.916 0.892 0.928 0.908 0.922 
#>   V14   V15   V16   V17   V18 
#> 0.905 0.924 0.934 0.907 0.923 
```
