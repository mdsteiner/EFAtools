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
  start_method = "psych",
  fiml_max_iter = 500,
  fiml_tol = 1e-05
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
  fit that is actually run with `estimator = "ML"`.

- fiml_max_iter:

  numeric. The maximum number of EM iterations used to estimate the
  two-stage full-information maximum-likelihood moments from raw data
  with missing values (`cor_method = "fiml"`); the last iterate is
  returned, with a warning, if the cap is reached. A single whole number
  of at least 1; default `500`. Not governed by `type`, and unused by
  every other correlation method. The EM converges linearly and needs
  more iterations the larger the fraction of missing information, so
  raise it when a fit reports that the moments did not converge.

- fiml_tol:

  numeric. The convergence tolerance of that EM: iteration stops once
  the largest change in the standardized moments (the standardized
  means, log-variances, and correlations) falls below it, so it does not
  depend on the variables' measurement scale. A single number greater
  than 0 and smaller than 1 (at or above 1 the criterion is met
  immediately and the starting moments would be returned as converged);
  default `1e-5`. Not governed by `type`.

- normalize:

  logical. If `TRUE` (default), a Kaiser normalization is performed
  before the rotation. The one knob that is always on unless you turn it
  off with `FALSE`.

- precision:

  numeric. The convergence tolerance of the rotation procedure. A single
  number greater than 0 and at most 1; default `1e-5`. Each rotation
  stage monitors its own quantity, so the same number is not the same
  tolerance everywhere: `varimax_type = "kaiser"` stops on the
  *absolute* change in the varimax simplicity criterion, which is an
  average over variables (and so does not scale with how many there are)
  but rises with the number of factors, roughly toward
  `1 - 1 / n_factors`, so a fixed value is a relatively weaker tolerance
  the more factors are extracted; `varimax_type = "svd"` stops on the
  *relative* change in the singular values (as in
  [`stats::varimax()`](https://rdrr.io/r/stats/varimax.html)); and the
  criterion rotations fitted by gradient projection stop when the
  *projected-gradient norm* falls below it. Promax inherits whichever of
  the two varimax tests its `varimax_type` selects, because it rotates a
  varimax base.

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
  `nrow(loadings)` for simplimax). Simplimax counts loadings, so a fit
  using it additionally requires a whole number no larger than the
  number of loadings in the solution; promax's power has no such
  restriction.

- random_starts:

  numeric. The number of random starts used by the criterion-based
  rotations to guard against local minima. A single whole number of at
  least 0, where `0` runs the rotation from its warm start only; default
  `100`. The default suffices for the smooth criteria; `simplimax`
  remains materially start-dependent at it, so raise it there (see the
  *Rotations* section of
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)).

- ...:

  Additional arguments forwarded to the rotation engine. Only the names
  a rotation engine can consume are accepted: `maxit` (a whole number of
  at least 0 bounding a *single* gradient-projection optimization – the
  multi-start search runs several of them and each is bounded
  separately, so it is not a budget for the run as a whole; varimax and
  promax have no such stage and take only `precision`), and the
  criterion parameters `gam` (oblimin; `gam = 0` is the recommended
  default, and larger values increasingly reward correlated factors and
  can drive the solution toward factor collapse, so inspect `Phi` before
  interpreting a fit with `gam > 0`) and `delta` (geomin; a positive
  number, default `0.01`); anything else is rejected as a misspelling.
  They are stored in `extra_args` and passed on to the rotation engine
  when the control is used to fit a model; an extra a given fit's
  rotation does not consume is ignored by that fit, so one control can
  serve fits with different rotations. An estimation knob (which belongs
  in `estimate_control()`) or one of the former spellings `P_type` and
  `randomStarts` is likewise rejected here, because the fit would
  silently drop it.

## Value

`estimate_control()` returns a list of class `efa_estimate_control` with
the components `type`, `init_comm`, `criterion`, `criterion_type`,
`max_iter`, `abs_eigen`, `start_method`, `fiml_max_iter`, and
`fiml_tol`. `rotate_control()` returns a list of class
`efa_rotate_control` with the components `type`, `normalize`,
`precision`, `order_type`, `varimax_type`, `p_type`, `k`,
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
#> fiml_max_iter: 500
#> fiml_tol: 1e-05

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
#> fiml_max_iter: 500
#> fiml_tol: 1e-05

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
#> fiml_max_iter: 500
#> fiml_tol: 1e-05


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
