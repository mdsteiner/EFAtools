# Control objects for estimation and rotation settings

`estimate_control()` and `rotate_control()` collect the estimation and
rotation *tuning* arguments of a factor analysis into two small,
validated objects. They are a declarative surface over the same settings
resolved internally by the package's estimation and rotation engines, so
that a fit's many tuning knobs can be prepared, inspected, and reused as
a single value instead of being passed one by one.

## Usage

``` r
estimate_control(
  type = c("EFAtools", "psych", "SPSS", "none"),
  init_comm = NA,
  criterion = NA,
  criterion_type = NA,
  max_iter = NA,
  abs_eigen = NA,
  start_method = "psych"
)

rotate_control(
  type = c("EFAtools", "psych", "SPSS", "none"),
  normalize = TRUE,
  precision = 1e-05,
  order_type = NA,
  varimax_type = NA,
  p_type = NA,
  k = NA,
  random_starts = 100,
  ...
)
```

## Arguments

- type:

  character. One of `"EFAtools"` (default), `"psych"`, `"SPSS"`, or
  `"none"`. Selects the preset that fills the `NA`-defaulted knobs below
  when the control is used to fit a model.

- init_comm:

  character. Method for the initial communalities in principal axis
  factoring: `"smc"` (squared multiple correlations), `"mac"` (maximum
  absolute correlations), or `"unity"`. `NA` (default) resolves from
  `type`.

- criterion:

  numeric. The convergence criterion for principal axis factoring:
  iteration stops once the change in communalities falls below it. A
  single number greater than 0 and smaller than 1; `NA` (default)
  resolves from `type`.

- criterion_type:

  character. The convergence criterion type for principal axis
  factoring: `"max_individual"` (the largest change in any communality,
  as in SPSS) or `"sum"` (the change in the summed communalities, as in
  [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html)). `NA`
  (default) resolves from `type`.

- max_iter:

  numeric. The maximum number of principal-axis-factoring iterations
  before the procedure is halted with a warning. A single whole number
  of at least 1; `NA` (default) resolves from `type`.

- abs_eigen:

  logical. Which algorithm the principal-axis-factoring iterations use:
  `FALSE` computes the loadings from the eigenvalues (as in
  [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html)); `TRUE` uses
  the absolute eigenvalues (as in SPSS). `NA` (default) resolves from
  `type`.

- start_method:

  character. Starting values for the maximum-likelihood optimiser:
  `"psych"` (default, the
  [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html) starts) or
  `"factanal"` (the
  [`stats::factanal()`](https://rdrr.io/r/stats/factanal.html) starts);
  abbreviations are matched. Not governed by `type`. Only maximum
  likelihood uses it, so `NA` leaves it unset and is rejected only by a
  fit that is actually run with `method = "ML"`.

- normalize:

  logical. If `TRUE` (default), a Kaiser normalization is performed
  before the rotation. The one knob that is always on unless you turn it
  off with `FALSE`.

- precision:

  numeric. The convergence tolerance of the rotation procedure. A single
  number greater than 0 and at most 1; default `1e-5`.

- order_type:

  character. How the factors are ordered: `"eigen"` (by descending sum
  of squared loadings, as in
  [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html)) or
  `"ss_factors"` (by descending unweighted sum of squared loadings).
  `NA` (default) resolves from `type`.

- varimax_type:

  character. The varimax variant used (for the varimax and promax
  rotations): `"svd"` (as in
  [`stats::varimax()`](https://rdrr.io/r/stats/varimax.html)) or
  `"kaiser"` (the SPSS / Kaiser (1958) procedure). `NA` (default)
  resolves from `type`.

- p_type:

  character. How the promax target matrix is computed: `"unnorm"` (the
  unnormalized target of Hendrickson & White (1964), also used by psych
  and stats) or `"norm"` (the normalized target used by SPSS). `NA`
  (default) resolves from `type`.

- k:

  numeric. The promax power (for the target matrix) or the number of
  near-zero loadings for simplimax. A single number greater than 0; `NA`
  (default) leaves it to the fit (the `type`-dependent promax value, or
  `nrow(loadings)` for simplimax).

- random_starts:

  numeric. The number of random starts used by the criterion-based
  rotations to guard against local minima. A single whole number of at
  least 0, where `0` runs the rotation from its warm start only; default
  `100`.

- ...:

  Additional arguments forwarded to the rotation engine (for example a
  criterion-specific `gam` for oblimin or `delta` for geomin, or
  `maxit`). They are stored in `extra_args` and passed on to the
  rotation engine when the control is used to fit a model. An estimation
  knob (which belongs in `estimate_control()`) or one of the former
  spellings `P_type` and `randomStarts` is rejected here, because the
  fit would silently drop it.

## Value

`estimate_control()` returns a list of class `efa_estimate_control` with
the components `type`, `init_comm`, `criterion`, `criterion_type`,
`max_iter`, `abs_eigen`, and `start_method`. `rotate_control()` returns
a list of class `efa_rotate_control` with the components `type`,
`normalize`, `precision`, `order_type`, `varimax_type`, `p_type`, `k`,
`random_starts`, and `extra_args` (a named list of any additional
arguments forwarded to the rotation engine).

## Details

Each argument that governs a `type` preset defaults to `NA`, meaning
"leave this knob to the preset". Setting `type` to one of `"EFAtools"`,
`"psych"`, or `"SPSS"` fills those knobs from the corresponding preset
when the fit is run; setting `type = "none"` requires the relevant knobs
to be supplied explicitly. The control object only records the chosen
`type` and the knobs you set: the preset is resolved (and any "argument
set alongside `type`" warning issued) when the object is used to fit a
model, exactly as it is today, because which preset applies depends on
the estimator and rotation.

## See also

[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md),
which takes both controls;
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md),
the retention criteria, and
[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md),
which take an `estimate_control` for the fits they run.

Other Control functions:
[`print.efa_control`](https://mdsteiner.github.io/EFAtools/reference/print.efa_control.md)

## Examples

``` r
# Estimation knobs taken entirely from a preset:
estimate_control(type = "SPSS")
#> Estimation control (type: "SPSS")
#> 
#> init_comm: <from type preset>
#> criterion: <from type preset>
#> criterion_type: <from type preset>
#> max_iter: <from type preset>
#> abs_eigen: <from type preset>
#> start_method: psych

# A preset with one knob pinned to a non-preset value:
estimate_control(type = "EFAtools", max_iter = 500)
#> Estimation control (type: "EFAtools")
#> 
#> init_comm: <from type preset>
#> criterion: <from type preset>
#> criterion_type: <from type preset>
#> max_iter: 500
#> abs_eigen: <from type preset>
#> start_method: psych

# Every knob supplied explicitly (type = "none"):
estimate_control(type = "none", init_comm = "smc", criterion = 1e-3,
                 criterion_type = "sum", max_iter = 300, abs_eigen = TRUE)
#> Estimation control (type: "none")
#> 
#> init_comm: smc
#> criterion: 0.001
#> criterion_type: sum
#> max_iter: 300
#> abs_eigen: TRUE
#> start_method: psych


# Rotation knobs taken from a preset:
rotate_control(type = "psych")
#> Rotation control (type: "psych")
#> 
#> normalize: TRUE
#> precision: 1e-05
#> order_type: <from type preset>
#> varimax_type: <from type preset>
#> p_type: <from type preset>
#> k: <from type preset>
#> random_starts: 100

# A criterion-specific extra argument, forwarded to the rotation engine:
rotate_control(type = "EFAtools", k = 3, gam = 0.5)
#> Rotation control (type: "EFAtools")
#> 
#> normalize: TRUE
#> precision: 1e-05
#> order_type: <from type preset>
#> varimax_type: <from type preset>
#> p_type: <from type preset>
#> k: 3
#> random_starts: 100
#> extra_args: gam = 0.5
```
