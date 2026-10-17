# Schmid-Leiman transformation

**\[superseded\]**

`SL()` has been superseded by
[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md),
which is the recommended interface going forward. It remains available
so existing code keeps working.

## Usage

``` r
SL(
  x,
  Phi = NULL,
  type = c("EFAtools", "psych", "SPSS", "none"),
  method = c("PAF", "ML", "ULS", "MINRES"),
  g_name = "g",
  ...
)
```

## Arguments

- x:

  An oblique
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  or [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html) solution, a
  pattern matrix (give its factor intercorrelations in `Phi`), or a
  second-order CFA solution from `lavaan`.

- Phi:

  matrix. The factor intercorrelations of the pattern matrix in `x`.
  With an
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  or [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html) solution in
  `x`, `Phi` replaces the fitted factor intercorrelations. Not used for
  a `lavaan` solution.

- type:

  character. One of "EFAtools" (default), "psych", "SPSS", or "none".
  This is used to control the procedure of the second-order factor
  analysis. In
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  it is set through the `type` of the
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  object.

- method:

  character. The estimator for the second-order factor analysis; passed
  to
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  as its `estimator` argument. One of `"PAF"`, `"ML"`, `"ULS"`, or
  `"MINRES"`.

- g_name:

  character. The name of the general factor of a `lavaan` second-order
  solution. Default is `"g"`.

- ...:

  The estimation tuning knobs of the second-order
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  (`init_comm`, `criterion`, `criterion_type`, `max_iter`, `abs_eigen`,
  `start_method`). They are repacked, together with `type`, into an
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  object, so that they tune that fit exactly as they always did. Other
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  arguments have no effect, and `se`, `b_boot`, `ci`, and `seed` give an
  error. The estimator is selected with `method`.

## Value

A list of class `c("efa_schmid_leiman", "SL")`, identical to the value
of
[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md);
see there for the components.

## See also

[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
