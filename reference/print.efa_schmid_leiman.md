# Print and format an efa_schmid_leiman object

[`print()`](https://rdrr.io/r/base/print.html) shows a summarised output
of the
[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
function: a model header (when the settings are available), the
Schmid-Leiman loading matrix, and the variances accounted for.
[`format()`](https://rdrr.io/r/base/format.html) assembles the same
report and returns it as a character vector;
[`print()`](https://rdrr.io/r/base/print.html) is
`cat(format(x), sep = "\n")`. The lines follow the active console theme,
so they are plain when colours are disabled (for example when captured
into a file or stripped with
[`cli::ansi_strip()`](https://cli.r-lib.org/reference/ansi_strip.html)).

## Usage

``` r
# S3 method for class 'efa_schmid_leiman'
print(x, ...)

# S3 method for class 'efa_schmid_leiman'
format(x, ...)
```

## Arguments

- x:

  An object of class `efa_schmid_leiman` (output from
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)).

- ...:

  Not used; for consistency with the generic.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its argument `x`
invisibly. [`format()`](https://rdrr.io/r/base/format.html) returns a
character vector with the report lines (styled to the active console
theme; plain when colours are disabled).

## Examples

``` r
EFA_mod <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
                   method = "PAF", rotation = "promax")
sl_mod <- efa_schmid_leiman(EFA_mod, method = "PAF")
sl_mod
#> 
#> EFA for second-order loadings performed with type = 'EFAtools' and method = 'PAF'
#> 
#> ── Schmid-Leiman Solution ──────────────────────────────────────────────────────
#> 
#>        g     F1     F2     F3    h2    u2
#> V1   .489  -.029   .022   .356  .367  .633
#> V2   .444  -.001   .042   .280  .277  .723
#> V3   .459   .036   .035   .263  .283  .717
#> V4   .522   .061  -.005   .320  .378  .622
#> V5   .468   .095  -.011   .254  .293  .707
#> V6   .478  -.044  -.031   .409  .399  .601
#> V7   .491   .001   .336   .054  .357  .643
#> V8   .463  -.010   .366   .018  .349  .651
#> V9   .457   .023   .347   .000  .330  .670
#> V10  .449  -.013   .425  -.041  .383  .617
#> V11  .477   .009   .224   .135  .297  .703
#> V12  .513   .012   .410  -.006  .432  .568
#> V13  .502   .372   .054  -.039  .395  .605
#> V14  .455   .332  -.043   .051  .322  .678
#> V15  .489   .340   .081  -.041  .363  .637
#> V16  .476   .336  -.032   .053  .343  .657
#> V17  .477   .402  -.023  -.016  .390  .610
#> V18  .485   .336   .003   .029  .350  .650
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      g     F1    F2     F3
#> SS loadings        4.111  .770  .783   .642
#> Prop Tot Var        .228  .043  .044   .036
#> Cum Prop Tot Var    .228  .271  .315   .350
#> Prop Comm Var       .652  .122  .124   .102
#> Cum Prop Comm Var   .652  .774  .898  1.000

# format() returns the same lines as plain text:
writeLines(format(sl_mod))
#> 
#> EFA for second-order loadings performed with type = 'EFAtools' and method = 'PAF'
#> 
#> ── Schmid-Leiman Solution ──────────────────────────────────────────────────────
#> 
#>        g     F1     F2     F3    h2    u2
#> V1   .489  -.029   .022   .356  .367  .633
#> V2   .444  -.001   .042   .280  .277  .723
#> V3   .459   .036   .035   .263  .283  .717
#> V4   .522   .061  -.005   .320  .378  .622
#> V5   .468   .095  -.011   .254  .293  .707
#> V6   .478  -.044  -.031   .409  .399  .601
#> V7   .491   .001   .336   .054  .357  .643
#> V8   .463  -.010   .366   .018  .349  .651
#> V9   .457   .023   .347   .000  .330  .670
#> V10  .449  -.013   .425  -.041  .383  .617
#> V11  .477   .009   .224   .135  .297  .703
#> V12  .513   .012   .410  -.006  .432  .568
#> V13  .502   .372   .054  -.039  .395  .605
#> V14  .455   .332  -.043   .051  .322  .678
#> V15  .489   .340   .081  -.041  .363  .637
#> V16  .476   .336  -.032   .053  .343  .657
#> V17  .477   .402  -.023  -.016  .390  .610
#> V18  .485   .336   .003   .029  .350  .650
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      g     F1    F2     F3
#> SS loadings        4.111  .770  .783   .642
#> Prop Tot Var        .228  .043  .044   .036
#> Cum Prop Tot Var    .228  .271  .315   .350
#> Prop Comm Var       .652  .122  .124   .102
#> Cum Prop Comm Var   .652  .774  .898  1.000
```
