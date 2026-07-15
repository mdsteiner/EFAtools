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
  from an oblique factor solution (see `Phi`) or a matrix of first-order
  factor loadings from a higher-order confirmatory factor analysis (see
  `L2`).

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

  character. One of "PAF", "ML", or "ULS" to use principal axis
  factoring, maximum likelihood, or unweighted least squares,
  respectively, used in
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  to find the second-order loadings. "MINRES" is accepted as a synonym
  for "ULS" (the same estimator).

- g_name:

  character. The name of the general factor. This needs only be
  specified if `x` is a `lavaan` second-order solution. Default is "g".

- ...:

  Arguments to be passed to
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md).
  The estimation tuning knobs are not passed here; they live in
  `estimate_control`.

## Value

A list of class `c("efa_schmid_leiman", "SL")`, identical to the value
of
[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md);
see there for the components.

## See also

[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
