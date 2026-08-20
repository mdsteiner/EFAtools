# Print and format an efa_average object

[`print()`](https://rdrr.io/r/base/print.html) shows a summarised output
of the
[`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md)
function: the averaging settings, the
error/convergence/Heywood/admissibility rates, the indicator-to-factor
correspondences, the averaged loadings (and, for oblique solutions, the
factor intercorrelations), the variances accounted for, and the model
fit. [`format()`](https://rdrr.io/r/base/format.html) assembles the same
report and returns it as a character vector; by default (`plot = FALSE`)
[`print()`](https://rdrr.io/r/base/print.html) is
`cat(format(x), sep = "\n")`. With `plot = TRUE` it additionally draws
the loading plot, which is the one thing
[`format()`](https://rdrr.io/r/base/format.html) cannot return: the
printed lines are the same, but the call has a side effect beyond them.
The lines follow the active console theme, so they are plain when
colours are disabled (for example when captured into a file or stripped
with
[`cli::ansi_strip()`](https://cli.r-lib.org/reference/ansi_strip.html)).

## Usage

``` r
# S3 method for class 'efa_average'
print(x, stat = c("average", "range"), plot = FALSE, ...)

# S3 method for class 'efa_average'
format(x, stat = c("average", "range"), ...)
```

## Arguments

- x:

  An object of class `efa_average` (output from
  [`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md)).

- stat:

  character. A vector with the statistics to print. Possible inputs are
  "average", "sd", "range", "min", and "max". Default is "average" and
  "range".

- plot:

  logical. Whether a plot of the average and min- max loadings should be
  created. Default is FALSE. If more than 10 factors are extracted, no
  plot is created. Only used by
  [`print()`](https://rdrr.io/r/base/print.html);
  [`plot.efa_average()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_average.md)
  draws the same plot and returns it, so it is the route to take when
  the plot object itself is wanted.

- ...:

  Not used; for consistency with the generic.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its argument `x`
invisibly. [`format()`](https://rdrr.io/r/base/format.html) returns a
character vector with the report lines.

## Examples

``` r
# \donttest{
EFA_aver <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500)
#> ℹ Extracting data
#> ✔ Extracting data [16ms]
#> 
#> ℹ Reordering factors
#> ✔ Reordering factors [29ms]
#> 
#> ℹ Averaging data
#> ✔ Averaging data [23ms]
#> 
EFA_aver
#> 
#> Averaging performed with averaging method mean (trim = 0) across 72 EFAs,
#> varying the following settings: init_comm, criterion_type, k_promax, p_type,
#> and varimax_type.
#> 
#> The error rate is at 0%. Of the solutions that did not result in an error, 100%
#> converged. Of the solutions that converged, 0% contained Heywood cases and 100%
#> were admissible.
#> 
#> ══ Indicator-to-Factor Correspondences ═════════════════════════════════════════
#> 
#> For each cell, the proportion of solutions including the respective
#> indicator-to-factor correspondence. A salience threshold of 0.3 was used to
#> determine indicator-to-factor correspondences.
#> 
#>       F1    F2    F3
#> V1    .00   .00  1.00
#> V2    .00   .00  1.00
#> V3    .00   .00  1.00
#> V4    .00   .00  1.00
#> V5    .00   .00  1.00
#> V6    .00   .00  1.00
#> V7    .00  1.00   .00
#> V8    .00  1.00   .00
#> V9    .00  1.00   .00
#> V10   .00  1.00   .00
#> V11   .00  1.00   .00
#> V12   .00  1.00   .00
#> V13  1.00   .00   .00
#> V14  1.00   .00   .00
#> V15  1.00   .00   .00
#> V16  1.00   .00   .00
#> V17  1.00   .00   .00
#> V18  1.00   .00   .00
#> 
#> ══ Loadings ════════════════════════════════════════════════════════════════════
#> 
#> ── Mean ────────────────────────────────────────────────────────────────────────
#> 
#>        F1     F2     F3
#> V1   -.019   .056   .588
#> V2    .022   .083   .466
#> V3    .079   .073   .440
#> V4    .119   .017   .532
#> V5    .169   .006   .426
#> V6   -.042  -.023   .672
#> V7    .027   .515   .110
#> V8    .010   .558   .051
#> V9    .060   .530   .022
#> V10   .004   .644  -.043
#> V11   .039   .351   .237
#> V12   .046   .625   .016
#> V13   .590   .105  -.035
#> V14   .527  -.039   .104
#> V15   .540   .143  -.039
#> V16   .534  -.022   .108
#> V17   .634  -.010   .000
#> V18   .535   .029   .071
#> 
#> ── Range (max − min) ───────────────────────────────────────────────────────────
#> 
#>       F1    F2    F3
#> V1   .098  .077  .092
#> V2   .076  .061  .063
#> V3   .064  .064  .052
#> V4   .069  .087  .069
#> V5   .046  .079  .047
#> V6   .108  .097  .121
#> V7   .059  .043  .052
#> V8   .059  .057  .064
#> V9   .048  .052  .072
#> V10  .059  .079  .090
#> V11  .061  .011  .015
#> V12  .058  .067  .084
#> V13  .072  .044  .099
#> V14  .056  .071  .049
#> V15  .062  .035  .098
#> V16  .055  .071  .052
#> V17  .081  .064  .084
#> V18  .055  .060  .064
#> 
#> ══ Factor Intercorrelations from Oblique Solutions ═════════════════════════════
#> 
#> ── Mean ────────────────────────────────────────────────────────────────────────
#> 
#>       F1     F2     F3
#> F1  1.000
#> F2   .536  1.000
#> F3   .568   .553  1.000
#> 
#> ── Range (max − min) ───────────────────────────────────────────────────────────
#> 
#>      F1    F2    F3
#> F1  .000
#> F2  .204  .000
#> F3  .255  .248  .000
#> 
#> ══ Variances Accounted for ═════════════════════════════════════════════════════
#> 
#> ── Mean ────────────────────────────────────────────────────────────────────────
#> 
#>                  F1     F2     F3
#> SS loadings    2.195  2.074  2.038
#> Prop Tot Var    .122   .115   .113
#> Prop Comm Var   .348   .329   .323
#> 
#> ── Range (max − min) ───────────────────────────────────────────────────────────
#> 
#>                 F1    F2    F3
#> SS loadings    .036  .052  .063
#> Prop Tot Var   .002  .003  .003
#> Prop Comm Var  .006  .008  .010
#> 
#> ══ Model Fit ═══════════════════════════════════════════════════════════════════
#> 
#> The fit indices are the mean of the per-solution fit indices, not the fit of
#> the averaged loadings printed above, which are a cell-wise summary rather than
#> a fitted solution.
#> 
#> CAF, RMSR, and SRMR averaged over 72 of 72 solutions.
#> 
#>        M (SD) [Min; Max]
#> CAF:  .50 (.00) [.50; .50]
#> RMSR: .03 (.00) [.03; .03]
#> SRMR: .02 (.00) [.02; .02]
#> df: 102

# format() returns the same lines as a character vector:
writeLines(format(EFA_aver))
#> 
#> Averaging performed with averaging method mean (trim = 0) across 72 EFAs,
#> varying the following settings: init_comm, criterion_type, k_promax, p_type,
#> and varimax_type.
#> 
#> The error rate is at 0%. Of the solutions that did not result in an error, 100%
#> converged. Of the solutions that converged, 0% contained Heywood cases and 100%
#> were admissible.
#> 
#> ══ Indicator-to-Factor Correspondences ═════════════════════════════════════════
#> 
#> For each cell, the proportion of solutions including the respective
#> indicator-to-factor correspondence. A salience threshold of 0.3 was used to
#> determine indicator-to-factor correspondences.
#> 
#>       F1    F2    F3
#> V1    .00   .00  1.00
#> V2    .00   .00  1.00
#> V3    .00   .00  1.00
#> V4    .00   .00  1.00
#> V5    .00   .00  1.00
#> V6    .00   .00  1.00
#> V7    .00  1.00   .00
#> V8    .00  1.00   .00
#> V9    .00  1.00   .00
#> V10   .00  1.00   .00
#> V11   .00  1.00   .00
#> V12   .00  1.00   .00
#> V13  1.00   .00   .00
#> V14  1.00   .00   .00
#> V15  1.00   .00   .00
#> V16  1.00   .00   .00
#> V17  1.00   .00   .00
#> V18  1.00   .00   .00
#> 
#> ══ Loadings ════════════════════════════════════════════════════════════════════
#> 
#> ── Mean ────────────────────────────────────────────────────────────────────────
#> 
#>        F1     F2     F3
#> V1   -.019   .056   .588
#> V2    .022   .083   .466
#> V3    .079   .073   .440
#> V4    .119   .017   .532
#> V5    .169   .006   .426
#> V6   -.042  -.023   .672
#> V7    .027   .515   .110
#> V8    .010   .558   .051
#> V9    .060   .530   .022
#> V10   .004   .644  -.043
#> V11   .039   .351   .237
#> V12   .046   .625   .016
#> V13   .590   .105  -.035
#> V14   .527  -.039   .104
#> V15   .540   .143  -.039
#> V16   .534  -.022   .108
#> V17   .634  -.010   .000
#> V18   .535   .029   .071
#> 
#> ── Range (max − min) ───────────────────────────────────────────────────────────
#> 
#>       F1    F2    F3
#> V1   .098  .077  .092
#> V2   .076  .061  .063
#> V3   .064  .064  .052
#> V4   .069  .087  .069
#> V5   .046  .079  .047
#> V6   .108  .097  .121
#> V7   .059  .043  .052
#> V8   .059  .057  .064
#> V9   .048  .052  .072
#> V10  .059  .079  .090
#> V11  .061  .011  .015
#> V12  .058  .067  .084
#> V13  .072  .044  .099
#> V14  .056  .071  .049
#> V15  .062  .035  .098
#> V16  .055  .071  .052
#> V17  .081  .064  .084
#> V18  .055  .060  .064
#> 
#> ══ Factor Intercorrelations from Oblique Solutions ═════════════════════════════
#> 
#> ── Mean ────────────────────────────────────────────────────────────────────────
#> 
#>       F1     F2     F3
#> F1  1.000
#> F2   .536  1.000
#> F3   .568   .553  1.000
#> 
#> ── Range (max − min) ───────────────────────────────────────────────────────────
#> 
#>      F1    F2    F3
#> F1  .000
#> F2  .204  .000
#> F3  .255  .248  .000
#> 
#> ══ Variances Accounted for ═════════════════════════════════════════════════════
#> 
#> ── Mean ────────────────────────────────────────────────────────────────────────
#> 
#>                  F1     F2     F3
#> SS loadings    2.195  2.074  2.038
#> Prop Tot Var    .122   .115   .113
#> Prop Comm Var   .348   .329   .323
#> 
#> ── Range (max − min) ───────────────────────────────────────────────────────────
#> 
#>                 F1    F2    F3
#> SS loadings    .036  .052  .063
#> Prop Tot Var   .002  .003  .003
#> Prop Comm Var  .006  .008  .010
#> 
#> ══ Model Fit ═══════════════════════════════════════════════════════════════════
#> 
#> The fit indices are the mean of the per-solution fit indices, not the fit of
#> the averaged loadings printed above, which are a cell-wise summary rather than
#> a fitted solution.
#> 
#> CAF, RMSR, and SRMR averaged over 72 of 72 solutions.
#> 
#>        M (SD) [Min; Max]
#> CAF:  .50 (.00) [.50; .50]
#> RMSR: .03 (.00) [.03; .03]
#> SRMR: .02 (.00) [.02; .02]
#> df: 102
# }
```
