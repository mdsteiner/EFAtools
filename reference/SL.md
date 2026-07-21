# Schmid-Leiman transformation

**\[superseded\]**

`SL()` has been superseded by
[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md),
which is the recommended interface going forward. It remains available
and unchanged so existing code keeps working.

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

  object of class
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md),
  class [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html), class
  [`lavaan::lavaan()`](https://rdrr.io/pkg/lavaan/man/lavaan.html), a
  matrix, or an `efa_loadings`/`loadings` object. If class
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  or class [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html),
  pattern coefficients and factor intercorrelations are taken from this
  object. If class
  [`lavaan::lavaan()`](https://rdrr.io/pkg/lavaan/man/lavaan.html), it
  must be a second-order CFA solution. In this case first-order and
  second-order factor loadings are taken from this object and the
  `g_name` argument has to be specified. x can also be a pattern matrix
  from an oblique factor solution (see `Phi`).

- Phi:

  matrix. A matrix of factor intercorrelations from an oblique factor
  solution. Only needs to be specified if a pattern matrix is entered
  directly into `x`.

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

  character. The name of the general factor. This needs only be
  specified if `x` is a `lavaan` second-order solution. Default is "g".

- ...:

  Further arguments passed on to the second-order
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md),
  including the estimation tuning knobs (`init_comm`, `criterion`,
  `criterion_type`, `max_iter`, `abs_eigen`, `start_method`), which are
  repacked, together with `type`, into an
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  object so that they tune that fit exactly as they always did. The
  estimator is selected with `method`.

## Value

A list of class `c("efa_schmid_leiman", "SL")`, identical to the value
of
[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md);
see there for the components.

## See also

[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
