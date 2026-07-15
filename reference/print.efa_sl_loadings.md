# Print an efa_sl_loadings object

Print an efa_sl_loadings object

## Usage

``` r
# S3 method for class 'efa_sl_loadings'
print(x, cutoff = 0.2, digits = 3, color = TRUE, ...)

# S3 method for class 'efa_sl_loadings'
format(x, ...)
```

## Arguments

- x:

  class efa_sl_loadings matrix.

- cutoff:

  numeric. The number above which to print loadings in bold (default is
  .2).

- digits:

  numeric. Passed to `round`. Number of digits to round the loadings to
  (default is 3).

- color:

  logical. Whether to apply console styling using cli. Default is
  `TRUE`.

- ...:

  additional arguments passed to print or format.

## Details

Prints a Schmid-Leiman loading matrix (general factor, group factors,
and the communality/uniqueness columns) as a styled, decimal-aligned
table. Loadings with absolute value greater than or equal to `cutoff`
are emphasised, smaller loadings are de-emphasised, and Heywood-relevant
cells (a loading or communality above 1, or a negative uniqueness) are
highlighted. If the matrix has many columns or the console is narrow,
the table is split into stacked column blocks so the output stays
readable.

## Examples

``` r
EFA_mod <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
                   method = "PAF", rotation = "promax")
efa_schmid_leiman(EFA_mod, method = "PAF")
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
