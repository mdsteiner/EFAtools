# Compare two vectors or matrices (communalities or loadings)

The function takes two objects of the same dimensions containing numeric
information (loadings or communalities) and returns a list of class
`efa_compare` containing summary information of the differences of the
objects.

## Usage

``` r
efa_compare(
  x,
  y,
  reorder = c("congruence", "names", "none"),
  corres = TRUE,
  thresh = 0.3,
  digits = 4,
  m_red = 0.001,
  range_red = 0.001,
  round_red = 3,
  print_diff = TRUE,
  na.rm = FALSE,
  x_labels = c("x", "y"),
  plot = TRUE,
  plot_red = 0.01
)
```

## Arguments

- x:

  matrix, or vector. Loadings or communalities of a factor analysis
  output.

- y:

  matrix, or vector. Loadings or communalities of another factor
  analysis output to compare to x.

- reorder:

  character. Whether and how elements / columns should be reordered. If
  "congruence" (default), the columns of `y` are matched to those of `x`
  by a joint one-to-one assignment that maximizes the total Tucker's
  congruence coefficient (a standard measure of similarity between two
  loading vectors) across all columns at once, and each matched column's
  sign is flipped if needed. This way, mismatched factor order or sign
  between two solutions does not distort the comparison. It applies to
  matrices only, and warns when `x` and `y` are vectors. If "names", the
  columns of a matrix – or the elements of a vector – are put in
  alphabetical order of their names; the rows of a matrix are assumed to
  be aligned already and are left untouched. If "none", no reordering is
  done.

- corres:

  logical. Whether factor correspondences should be compared if a matrix
  is entered. Default is TRUE.

- thresh:

  numeric. The threshold at or above which a loading is classified as
  substantial. Default is .3.

- digits:

  numeric. Number of decimals to print in the output. Default is 4.

- m_red:

  numeric. Number above which the mean and median should be printed in
  red (i.e., if .001 is used, the mean will be in red if it is larger
  than .001, otherwise it will be displayed in green.) Default is .001.

- range_red:

  numeric. Number above which the min and max should be printed in red
  (i.e., if .001 is used, min and max will be in red if the max is
  larger than .001, otherwise it will be displayed in green). Default is
  .001. Note that the color of min also depends on max, that is min will
  be displayed in the same color as max.

- round_red:

  numeric. The number of agreeing decimals below which the report
  highlights the agreement in red (i.e., if 3 is used, the value is
  shown in red when the compared numbers agree to fewer than 3 decimals,
  otherwise in green). Default is 3.

- print_diff:

  logical. Whether the difference vector or matrix should be printed or
  not. Default is TRUE.

- na.rm:

  logical. Whether NAs should be removed from the difference summaries
  and factor-correspondence classifications. With `FALSE`, a missing
  loading makes the correspondence counts undefined (`NA`). Default is
  FALSE.

- x_labels:

  character. A vector of length two containing identifying labels for
  the two objects x and y that will be compared. These will be used as
  labels on the x-axis of the plot, and to name the direction of the
  signed elementwise differences in the printed report (see
  [`print.efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_compare.md)).
  Default is "x" and "y".

- plot:

  **\[superseded\]** Accepted and validated, but without effect;
  retained for backwards compatibility. The difference plot is drawn by
  [`plot.efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_compare.md).
  Default is TRUE.

- plot_red:

  numeric. Threshold above which to plot the absolute differences in
  red. Default is .01.

## Value

A list of class `efa_compare` with the following components:

- diff:

  The vector or matrix containing the differences between x and y.

- mean_abs_diff:

  The mean absolute difference between x and y.

- median_abs_diff:

  The median absolute difference between x and y.

- min_abs_diff:

  The minimum absolute difference between x and y.

- max_abs_diff:

  The maximum absolute difference between x and y.

- max_dec:

  The maximum number of decimals to which a comparison makes sense. For
  example, if x contains only values up to the third decimals, and y is
  a normal double, max_dec will be three.

- are_equal:

  The maximal number of decimals to which all elements of x and y agree
  in absolute value. The comparison is on magnitudes, so two elements
  that are equal in size but opposite in sign count as agreeing; signed
  disagreements are reflected in `diff` and the mean / median / min /
  max absolute differences. `0` means the two agree in their integer
  parts but in no decimal place. `NA` means there is no agreement at
  all: either they already differ in their integer parts, or
  `na.rm = FALSE` and an element is missing.

- diff_corres:

  The number of differing variable-to-factor correspondences between x
  and y, when only the highest loading is considered. `NA` whenever the
  correspondences were not compared: for vector input, for a matrix with
  a single column, with `corres = FALSE`, and when a loading is missing
  under `na.rm = FALSE`.

- diff_corres_cross:

  The number of differing variable-to-factor correspondences between x
  and y when all loadings `>= thresh` are considered. `NA` under the
  same conditions as `diff_corres`.

- g:

  The root mean squared distance (RMSE) between x and y.

- settings:

  List of the settings used.

## Details

`digits`, `m_red`, `range_red`, `round_red`, `print_diff`, and
`plot_red` only control how the result is displayed; each is stored in
the returned object's `settings` and can be overridden later without
recomputing the comparison – `digits`, `m_red`, `range_red`,
`round_red`, and `print_diff` in a call to
[`print.efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_compare.md),
and `plot_red` in a call to
[`plot.efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_compare.md).

## See also

[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
for the solutions being compared, and
[`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md)
to rotate one solution onto another before comparing.

## Examples

``` r
# A type SPSS EFA to mimick the SPSS implementation
EFA_SPSS_6 <- efa_fit(test_models$case_11b$cormat, n_factors = 6,
                      estimate_control = estimate_control(type = "SPSS"),
                      rotate_control = rotate_control(type = "SPSS"))
#> Warning: Reached the maximum number of iterations without convergence; results may not
#> be interpretable.

# A type psych EFA to mimick the psych::fa() implementation
EFA_psych_6 <- efa_fit(test_models$case_11b$cormat, n_factors = 6,
                       estimate_control = estimate_control(type = "psych"),
                       rotate_control = rotate_control(type = "psych"))

# compare the two
efa_compare(EFA_SPSS_6$unrot_loadings, EFA_psych_6$unrot_loadings,
            x_labels = c("SPSS", "psych"))
#> 
#> ── Summary statistics ──────────────────────────────────────────────────────────
#> 
#> Mean [min, max] absolute difference:  .0025 [ .0000,  .0215]
#> Median absolute difference:  .0008
#> Root mean squared distance (RMSE):  .0048
#> Max decimals where all numbers agree in absolute value: 0
#> Differing indicator-to-factor correspondences: 0 (highest loading),
#>   0 (all |loadings| >= 0.3)
#> 
#> ── Elementwise differences ─────────────────────────────────────────────────────
#> 
#> Differences: SPSS - psych.
#> 
#>        F1      F2      F3      F4      F5      F6
#> V1    .0000  -.0001  -.0003  -.0053  -.0015  -.0012
#> V2    .0001   .0005   .0008  -.0080  -.0048  -.0019
#> V3   -.0005  -.0002  -.0015  -.0097  -.0131   .0000
#> V4    .0001  -.0001  -.0011   .0049   .0128   .0020
#> V5    .0001  -.0001  -.0007   .0027   .0131   .0007
#> V6   -.0008  -.0030  -.0053  -.0178   .0215   .0008
#> V7    .0000  -.0005  -.0009  -.0047   .0041  -.0001
#> V8    .0000   .0002   .0000  -.0085  -.0027  -.0018
#> V9    .0002   .0005   .0006  -.0067  -.0063  -.0002
#> V10   .0000  -.0002  -.0007   .0005   .0008   .0006
#> V11   .0003   .0007  -.0012   .0025  -.0024  -.0006
#> V12  -.0004  -.0018   .0020  -.0014  -.0021  -.0001
#> V13   .0001   .0018   .0018   .0137  -.0031   .0002
#> V14   .0000   .0006   .0004   .0093   .0035   .0001
#> V15  -.0001   .0011   .0006   .0171   .0025   .0004
#> V16   .0001   .0002   .0011   .0042  -.0027   .0001
#> V17   .0001   .0000   .0008   .0021   .0008   .0006
#> V18   .0001   .0001   .0007   .0003   .0016   .0006
```
