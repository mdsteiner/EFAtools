# Schmid-Leiman transformation

This function implements the Schmid-Leiman (SL) transformation (Schmid &
Leiman, 1957). It takes the pattern coefficients and factor
intercorrelations from an oblique factor solution as input and can
reproduce the results from
[`psych::schmid()`](https://rdrr.io/pkg/psych/man/schmid.html) and from
the SPSS implementation from Wolff & Preising (2005). Other arguments
from
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
can be used to control the procedure to find the second-order loadings
more flexibly. The function can also be used on a second-order
confirmatory factor analysis (CFA) solution from lavaan.

## Usage

``` r
efa_schmid_leiman(
  x,
  Phi = NULL,
  method = c("PAF", "ML", "ULS", "MINRES"),
  g_name = "g",
  estimate_control = NULL,
  ...
)
```

## Source

Schmid, J. & Leiman, J. M. (1957). The development of hierarchical
factor solutions. Psychometrika, 22(1), 53–61. doi:10.1007/BF02289209

Wolff, H.-G., & Preising, K. (2005). Exploring item and higher order
factor structure with the Schmid-Leiman solution: Syntax codes for SPSS
and SAS. Behavior Research Methods, 37 , 48–58. doi:10.3758/BF03206397

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

- estimate_control:

  an
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  object with the estimation settings for the second-order
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  fit, including the `type` preset. `NULL` (default) uses the
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  defaults. The second-order fit is unrotated, so no rotation settings
  apply.

- ...:

  Arguments to be passed to
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md).
  The estimation tuning knobs are not passed here; they live in
  `estimate_control`.

## Value

A list of class `c("efa_schmid_leiman", "SL")` containing the following

- orig_R:

  Original correlation matrix.

- sl:

  A matrix with general factor loadings, group factor loadings,
  communalities, and uniquenesses.

- L2:

  Second-order factor loadings.

- vars_accounted:

  A matrix of explained variances and sums of squared loadings.

- iter:

  The number of iterations needed for convergence in EFA.

- convergence:

  Integer convergence code of the second-order EFA (0 = converged); `NA`
  for a lavaan input. See
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md).

- settings:

  list. The settings (arguments) used in EFA to get the second-order
  loadings.

## Details

The SL transformation (also called SL orthogonalization) is a procedure
with which an oblique factor solution is transformed into a
hierarchical, orthogonalized solution. As a first step, the factor
intercorrelations are factor analyzed to extract a single second-order
(general) factor, yielding a two-level hierarchical structure. The
first-order factor and the second-order factor are then orthogonalized,
resulting in an orthogonalized factor solution with proportionality
constraints. The procedure thus makes a suggested hierarchical data
structure based on factor intercorrelations explicit. One major
advantage of SL transformation is that it enables variance partitioning
between higher-order and first-order factors, including the calculation
of McDonald's omegas (see
[`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md)).

## See also

Other factor rotation:
[`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md)

## Examples

``` r
## Use with an output from the EFAtools::efa_fit function, both with type EFAtools
EFA_mod <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
                   method = "PAF", rotation = "promax")
SL_EFAtools <- efa_schmid_leiman(EFA_mod, method = "PAF",
                                 estimate_control = estimate_control(type = "EFAtools"))

# \donttest{
## Use with an output from the psych::fa function with type psych
fa_mod <- psych::fa(test_models$baseline$cormat, nfactors = 3, n.obs = 500,
                    fm = "pa", rotate = "Promax")
SL_psych <- efa_schmid_leiman(fa_mod, method = "PAF",
                              estimate_control = estimate_control(type = "psych"))
# }

## Use more flexibly by entering a pattern matrix and phi directly (useful if
## a factor solution found with another program should be subjected to SL
## transformation)

## For demonstration, take pattern matrix and phi from an EFA output
## This gives the same solution as the first example
EFA_mod <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
                   method = "PAF", rotation = "promax")
SL_flex <- efa_schmid_leiman(EFA_mod$rot_loadings, Phi = EFA_mod$Phi, method = "PAF",
                             estimate_control = estimate_control(type = "EFAtools"))

# \donttest{
## Use with a lavaan second-order CFA output
if (requireNamespace("lavaan", quietly = TRUE)) {

# Create and fit model in lavaan (assume all variables have SDs of 1)
mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
        F2 =~ V7 + V8 + V9 + V10 + V11 + V12
        F3 =~ V13 + V14 + V15 + V16 + V17 + V18
        g =~ F1 + F2 + F3'
fit <- lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
                   sample.nobs = 500, estimator = "ml")

SL_lav <- efa_schmid_leiman(fit, g_name = "g")

}
# }
```
