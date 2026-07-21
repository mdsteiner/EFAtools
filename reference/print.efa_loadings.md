# Print a loading matrix

Print a loading matrix

## Usage

``` r
# S3 method for class 'efa_loadings'
print(
  x,
  cutoff = 0.3,
  digits = 3,
  max_name_length = 10,
  h2 = NULL,
  color = TRUE,
  name_style = c("truncate", "abbreviate", "full"),
  max_factor_name_length = NULL,
  max_factors_per_block = NULL,
  sort_loadings = c("none", "primary", "clustered"),
  legend = FALSE,
  ...
)

# S3 method for class 'efa_loadings'
format(x, ...)
```

## Arguments

- x:

  a loading matrix of class `efa_loadings` (or the legacy `LOADINGS`).

- cutoff:

  numeric. The value at or above which loadings are emphasized; default
  is .3.

- digits:

  numeric. Passed to `round`. Number of digits to round the loadings to
  (default is 3).

- max_name_length:

  numeric. The maximum length of the variable names to display.
  Everything beyond this will be cut from the right unless
  `name_style = "abbreviate"` or `name_style = "full"` is used.

- h2:

  numeric. Vector of communalities to print. If named and `x` has row
  names, names are used to align communalities to rows.

- color:

  logical. Whether to apply console styling using cli. Default is
  `TRUE`.

- name_style:

  character. How to shorten variable names longer than
  `max_name_length`. `"truncate"` cuts names from the right,
  `"abbreviate"` uses
  [`base::abbreviate()`](https://rdrr.io/r/base/abbreviate.html), and
  `"full"` prints full names.

- max_factor_name_length:

  numeric or `NULL`. Optional maximum length of factor names. If `NULL`,
  factor names are not shortened.

- max_factors_per_block:

  numeric or `NULL`. Maximum number of factor columns to print per
  block. If `NULL`, the number is chosen from the console width.

- sort_loadings:

  character. Optional row sorting. `"none"` preserves the input order,
  `"primary"` groups rows by the factor with the largest absolute
  loading, and `"clustered"` additionally sorts within each factor by
  the size of the primary loading.

- legend:

  logical. Whether to append a short explanation of the styling. Default
  is `FALSE` for standalone loading matrices. The legend is only printed
  when the styling it describes is actually rendered (`color = TRUE` and
  a colour-capable console); in plain output it is omitted and this
  argument has no effect.

- ...:

  additional arguments passed to print or format

## Details

The method prints a loading matrix in a compact, console-oriented table.
Loadings with absolute value greater than or equal to `cutoff` are
emphasized, smaller loadings are de-emphasized, and Heywood-relevant
communality/ uniqueness values are marked when `h2` is supplied. Long
variable names can be truncated, abbreviated, or printed in full. If the
matrix has many factor columns, the table is split into column blocks so
that the output remains readable in narrower consoles.

If `h2` is named and `x` has row names, `h2` is matched to the row names
of `x` before any optional row sorting is applied. If `x` has no row
names, a named `h2` vector is used in the supplied order.

## Examples

``` r
EFAtools_PAF <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
                        estimator = "PAF", rotation = "promax")
EFAtools_PAF
#> 
#> EFA performed with estimator = 'PAF' and rotation = 'promax'.
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
