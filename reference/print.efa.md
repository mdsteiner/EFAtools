# Print and summarise an efa object

[`print()`](https://rdrr.io/r/base/print.html) shows a concise overview
of an
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
or
[`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)
solution: a model header, the loading matrix (with the factor
intercorrelations for oblique solutions), the variances accounted for,
and the model fit. [`summary()`](https://rdrr.io/r/base/summary.html)
returns a `summary.efa` object whose print method adds the full
diagnostics: model and simple-structure diagnostics, confidence-interval
tables, the structure matrix, multiple-imputation uncertainty (for
pooled objects), and residual diagnostics.
[`format()`](https://rdrr.io/r/base/format.html) assembles the same
report and returns it as a character vector;
[`print()`](https://rdrr.io/r/base/print.html) is
`cat(format(x), sep = "\n")`. The lines follow the active console theme,
so they are plain when colours are disabled (for example when captured
into a file or stripped with
[`cli::ansi_strip()`](https://cli.r-lib.org/reference/ansi_strip.html)).

## Usage

``` r
# S3 method for class 'efa'
print(x, ...)

# S3 method for class 'efa_mi'
print(x, ...)

# S3 method for class 'efa'
format(
  x,
  cutoff = 0.3,
  digits = 3,
  max_name_length = 10,
  sort_loadings = c("none", "primary", "clustered"),
  show_loading_legend = TRUE,
  max_factors_per_block = NULL,
  ...
)

# S3 method for class 'efa_mi'
format(x, ...)

# S3 method for class 'efa'
summary(
  object,
  cutoff = 0.3,
  digits = 3,
  max_name_length = 10,
  ci = c("auto", "none", "separate"),
  ci_filter = c("salient", "all", "nonzero"),
  diagnostics_top_n = 10,
  residual_cutoff = 0.1,
  residual_top_n = 10,
  show_structure = TRUE,
  sort_loadings = c("none", "primary", "clustered"),
  show_loading_legend = TRUE,
  cross_loading_cutoff = cutoff,
  min_primary_gap = 0.2,
  min_salient_per_factor = 3,
  max_factors_per_block = NULL,
  show_mi_diagnostics = NULL,
  ...
)

# S3 method for class 'efa_mi'
summary(object, ...)

# S3 method for class 'summary.efa'
print(x, ...)

# S3 method for class 'summary.efa'
format(x, ...)
```

## Arguments

- x, object:

  An object of class `efa` (from
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md))
  or `efa_mi` (from
  [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md));
  for the `summary.efa` methods, the object returned by
  [`summary()`](https://rdrr.io/r/base/summary.html).

- ...:

  Further arguments passed to
  [`print.efa_loadings()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_loadings.md).

- cutoff:

  numeric. The absolute value at or above which loadings are emphasised
  in the loading table. Default is .3.

- digits:

  numeric. Number of decimal places for the printed tables. Default is
  3.

- max_name_length:

  numeric. Maximum length of the variable names to display; longer names
  are cut from the right.

- sort_loadings:

  character. Optional row sorting for the loading table. See
  [`print.efa_loadings()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_loadings.md).

- show_loading_legend:

  logical. Whether to print a short legend for the loading-table
  styling. Default is `TRUE`.

- max_factors_per_block:

  numeric or `NULL`. Maximum number of factor columns per loading-table
  block. If `NULL`, chosen from the console width.

- ci:

  character. Which confidence intervals
  [`summary()`](https://rdrr.io/r/base/summary.html) shows, if
  available. `"auto"` and `"separate"` print CI sections when CIs were
  computed; `"none"` suppresses them. Default is `"auto"`.

- ci_filter:

  character. Which loading CIs
  [`summary()`](https://rdrr.io/r/base/summary.html) prints: `"salient"`
  (default), `"all"`, or `"nonzero"`; see Details.

- diagnostics_top_n:

  numeric. Maximum number of item-level entries
  [`summary()`](https://rdrr.io/r/base/summary.html) prints per
  simple-structure diagnostic.

- residual_cutoff:

  numeric. Absolute residual cutoff for the residual diagnostics in
  [`summary()`](https://rdrr.io/r/base/summary.html). Default is .1.

- residual_top_n:

  numeric. Maximum number of residuals
  [`summary()`](https://rdrr.io/r/base/summary.html) prints. Use `Inf`
  to print all residuals above `residual_cutoff`.

- show_structure:

  logical. Whether [`summary()`](https://rdrr.io/r/base/summary.html)
  prints the structure matrix for oblique solutions when available.
  Default is `TRUE`.

- cross_loading_cutoff:

  numeric. Cutoff for counting cross-loadings in the
  [`summary()`](https://rdrr.io/r/base/summary.html) diagnostics.
  Defaults to `cutoff`.

- min_primary_gap:

  numeric. Minimum desired absolute difference between the largest and
  second-largest absolute loading of an item, used in the
  [`summary()`](https://rdrr.io/r/base/summary.html) diagnostics.

- min_salient_per_factor:

  numeric. Minimum number of salient indicators per factor used in the
  [`summary()`](https://rdrr.io/r/base/summary.html) diagnostics.
  Default is 3.

- show_mi_diagnostics:

  logical or `NULL`. Whether
  [`summary()`](https://rdrr.io/r/base/summary.html) prints a
  multiple-imputation uncertainty summary for pooled EFAs. `NULL` shows
  it for pooled objects.

## Value

[`print()`](https://rdrr.io/r/base/print.html) and the print method for
`summary.efa` objects return their argument invisibly.
[`format()`](https://rdrr.io/r/base/format.html) returns a character
vector with the report lines (styled to the active console theme; plain
when colours are disabled).
[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.efa`.

## Details

The methods are shared by single-imputation `efa` objects and pooled
`efa_mi` objects. For `efa_mi` objects the header reports the number of
imputations and the alignment/pooling settings; confidence intervals and
a multiple-imputation uncertainty summary are shown by
[`summary()`](https://rdrr.io/r/base/summary.html) when the pooled
object carries bootstrap/MI quantities.

In [`summary()`](https://rdrr.io/r/base/summary.html), `ci_filter`
controls which loading intervals are shown: `"salient"` reports
intervals for loadings whose absolute point estimate is at least
`cutoff`, `"nonzero"` reports intervals excluding zero, and `"all"`
reports every finite interval.

## Examples

``` r
mod <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
               method = "PAF", rotation = "promax")
mod
#> 
#> EFA performed with type = 'EFAtools', method = 'PAF', and rotation = 'promax'.
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
#> CAF: .50
#> SRMR: .02
#> df: 102

# The full diagnostics, CI tables, and residual diagnostics:
summary(mod)
#> 
#> EFA performed with type = 'EFAtools', method = 'PAF', and rotation = 'promax'.
#> 
#> ── Model Diagnostics ───────────────────────────────────────────────────────────
#> 
#> Factors: 3
#> Variables: 18
#> N: 500
#> Heywood cases: 0
#> Cross-loading items (|loading| >= .300): 0
#> Items without salient loading (|loading| >= .300): 0
#> Factors with fewer than 3 salient indicators: 0
#> Items with primary-loading gap < .200: 1
#> Largest |residual|: .069
#> Factor intercorrelations > .85: none
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
#> ── Structure Matrix ────────────────────────────────────────────────────────────
#> 
#>       F1    F2    F3
#> V1   .371  .394  .605
#> V2   .353  .371  .524
#> V3   .388  .379  .527
#> V4   .453  .402  .611
#> V5   .430  .356  .529
#> V6   .354  .352  .627
#> V7   .391  .593  .431
#> V8   .362  .590  .387
#> V9   .378  .574  .372
#> V10  .349  .616  .341
#> V11  .385  .512  .467
#> V12  .416  .657  .415
#> V13  .624  .424  .386
#> V14  .563  .326  .400
#> V15  .595  .430  .375
#> V16  .583  .350  .419
#> V17  .623  .355  .380
#> V18  .590  .378  .412
#> 
#> ── Simple Structure Diagnostics ────────────────────────────────────────────────
#> 
#> Items with primary-loading gap < .200:
#> • V11: F2 = .356, F3 = .232
#> 
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
#> CAF: .50
#> SRMR: .02
#> df: 102
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

# format() returns plain text, e.g. for embedding in a report:
writeLines(format(mod))
#> 
#> EFA performed with type = 'EFAtools', method = 'PAF', and rotation = 'promax'.
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
#> CAF: .50
#> SRMR: .02
#> df: 102
```
