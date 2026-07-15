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
  by an optimal one-to-one assignment that maximizes Tucker's congruence
  coefficient, so each factor is matched exactly once; if "names",
  objects are reordered according to their names; if "none", no
  reordering is done.

- corres:

  logical. Whether factor correspondences should be compared if a matrix
  is entered.

- thresh:

  numeric. The threshold to classify a pattern coefficient as
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
  larger than .001, otherwise it will be displayed in green. Default is
  .001). Note that the color of min also depends on max, that is min
  will be displayed in the same color as max.

- round_red:

  numeric. Number above which the max decimals to round to where all
  corresponding elements of x and y are still equal are displayed in red
  (i.e., if 3 is used, the number will be in red if it is smaller than
  3, otherwise it will be displayed in green). Default is 3.

- print_diff:

  logical. Whether the difference vector or matrix should be printed or
  not. Default is TRUE.

- na.rm:

  logical. Whether NAs should be removed in the mean, median, min, and
  max functions. Default is FALSE.

- x_labels:

  character. A vector of length two containing identifying labels for
  the two objects x and y that will be compared. These will be used as
  labels on the x-axis of the plot. Default is "x" and "y".

- plot:

  logical. Retained for backwards compatibility; the difference plot is
  now drawn with
  [`plot.efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_compare.md)
  rather than when printing. Default is TRUE.

- plot_red:

  numeric. Threshold above which to plot the absolute differences in
  red. Default is .01.

## Value

A list of class `c("efa_compare", "COMPARE")`, identical to the value of
[`efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/efa_compare.md);
see there for the components.

## See also

[`efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/efa_compare.md)
