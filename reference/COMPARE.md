# Compare two vectors or matrices (communalities or loadings)

**\[superseded\]**

`COMPARE()` has been superseded by
[`efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/efa_compare.md),
which is the recommended interface going forward. It remains available
and unchanged so existing code keeps working.

## Usage

``` r
COMPARE(
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

A list of class `c("efa_compare", "COMPARE")`, identical to the value of
[`efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/efa_compare.md);
see there for the components.

## See also

[`efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/efa_compare.md)
