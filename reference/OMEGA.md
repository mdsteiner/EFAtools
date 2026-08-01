# McDonald's omega

**\[superseded\]**

`OMEGA()` has been superseded by
[`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md),
which is the recommended interface going forward. It remains available
and unchanged so existing code keeps working.

This function finds omega total, hierarchical, and subscale, as well as
additional model-based indices of interpretive relevance (H index, ECV,
PUC) from a Schmid-Leiman (SL) solution or lavaan single factor,
second-order (see below), or bifactor solution. The SL-based omegas can
either be found from a
[`psych::schmid()`](https://rdrr.io/pkg/psych/man/schmid.html),
[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md),
or, in a more flexible way, by leaving `model = NULL` and specifying
additional arguments. By setting the `type` argument, results from
[`psych::omega()`](https://rdrr.io/pkg/psych/man/omega.html) can be
reproduced.

## Usage

``` r
OMEGA(
  model = NULL,
  type = c("EFAtools", "psych"),
  g_name = "g",
  group_names = NULL,
  add_ind = TRUE,
  factor_corres = NULL,
  var_names = NULL,
  fac_names = NULL,
  g_load = NULL,
  s_load = NULL,
  u2 = NULL,
  cormat = NULL,
  pattern = NULL,
  Phi = NULL,
  variance = c("correlation", "sums_load")
)
```

## Source

McDonald, R. P. (1978). Generalizability in factorable domains: ‘‘Domain
validity and generalizability’’. Educational and Psychological
Measurement, 38, 75–79.

McDonald, R. P. (1985). Factor analysis and related methods. Hillsdale,
NJ: Erlbaum.

McDonald, R. P. (1999). Test theory: A unified treatment. Mahwah, NJ:
Erlbaum.

Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016a). Applying
bifactor statistical indices in the evaluation of psychological
measures. Journal of Personality Assessment, 98, 223-237.

Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016b). Evaluating
bifactor models: Calculating and interpreting statistical indices.
Psychological Methods, 21, 137-150.

Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
reliability within latent variable systems. In R. Cudeck, S. du Toit, &
D. Sörbom (Eds.), Structural equation modeling: Present and future—A
Festschrift in honor of Karl Jöreskog (pp. 195–216). Lincolnwood, IL:
Scientific Software International.

Reise, S. P., Scheines, R., Widaman, K. F., & Haviland, M. G. (2013).
Multidimensionality and structural coefficient bias in structural
equation modeling: A bifactor perspective. Educational and Psychological
Measurement, 73, 5–26.

Bonifay, W. E., Reise, S. P., Scheines, R., & Meijer, R. R. (2015). When
are multidimensional data unidimensional enough for structural equation
modeling?: An evaluation of the DETECT multidimensionality index.
Structural Equation Modeling, 22, 504—516.

Gignac, G. E. (2014). On the Inappropriateness of Using Items to
Calculate Total Scale Score Reliability via Coefficient Alpha for
Multidimensional Scales. European Journal of Psychological Assessment,
30, 130-139.

## Arguments

- model:

  class
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md),
  class `schmid`, or class `lavaan` object. That is, an output object
  from
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  or [`psych::schmid()`](https://rdrr.io/pkg/psych/man/schmid.html), or
  a `lavaan` fit object with a single factor, second-order, or bifactor
  solution. If of class `lavaan`, only `g_name` needs to be specified
  additionally. If of class
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  or `schmid`, only the arguments `factor_corres` and `cormat` need to
  be specified additionally.

- type:

  character. Either `"EFAtools"` (default) or `"psych"` (see details)

- g_name:

  character. The name of the general factor from the lavaan solution.
  This needs only be specified if `model` is a `lavaan` second-order or
  bifactor solution. Default is "g".

- group_names:

  character. An optional vector of group names. The length must
  correspond to the number of groups for which the `lavaan` model was
  fitted.

- add_ind:

  logical. Whether additional indices (H index, ECV, PUC) should be
  calculated or not (see details for these indices). If FALSE, only
  omegas are returned. Default is `TRUE`.

- factor_corres:

  matrix. A logical matrix or a numeric matrix containing 0's and 1's
  that indicates which variable corresponds to which group factor. Must
  have the same dimensions as the matrix of group factor loadings from
  the SL solution. Cross-loadings are allowed here. See examples for
  use.

- var_names:

  character. A vector with subtest names in the order of the rows from
  the SL solution. This needs only be specified if `model` is left
  `NULL`.

- fac_names:

  character. An optional vector of group factor names in the order of
  the columns of the SL solution. If left `NULL`, names of the group
  factors from the entered solution are taken.

- g_load:

  numeric. A vector of general factor loadings from an SL solution. This
  needs only be specified if `model` is left `NULL`.

- s_load:

  matrix. A matrix of group factor loadings from an SL solution. This
  needs only be specified if `model` is left `NULL`.

- u2:

  numeric. A vector of uniquenesses from an SL solution. This needs only
  be specified if `model` is left `NULL`.

- cormat:

  matrix. A correlation matrix to be used when
  `variance = "correlation"`. If left `NULL` and an
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  output is entered in `model`, the correlation matrix is taken from the
  output. If left `NULL` and a
  [`psych::schmid()`](https://rdrr.io/pkg/psych/man/schmid.html) output
  is entered, the correlation matrix will be found based on the pattern
  matrix and Phi from the
  [`psych::schmid()`](https://rdrr.io/pkg/psych/man/schmid.html) output
  using
  [`psych::factor.model()`](https://rdrr.io/pkg/psych/man/factor.model.html).
  If left `NULL` and model is also left `NULL`, the correlation matrix
  is found based on the pattern matrix and Phi entered. However, if the
  correlation matrix is available, `cormat` should be specified instead
  of `Phi` and `pattern`.

- pattern:

  matrix. Pattern coefficients from an oblique factor solution. This
  needs only be specified if `model` is left `NULL`,
  `variance = "correlation"` and `cormat` is also left `NULL`.

- Phi:

  matrix. Factor intercorrelations from an oblique factor solution. This
  needs only be specified if `model` is left `NULL`,
  `variance = "correlation"` and `cormat` is also left `NULL`.

- variance:

  character. If `"correlation"` (default), then total variances for the
  whole scale as well as for the subscale composites are calculated
  based on the correlation matrix. If `"sums_load"`, then total
  variances are calculated using the squared sums of general factor
  loadings and group factor loadings and the sum of uniquenesses (see
  details).

## Value

If found for an SL or `lavaan` second-order or bifactor solution without
multiple groups: A matrix with omegas for the whole scale and for the
subscales and (only if `add_ind = TRUE`) with the H index, ECV, and PUC.

- tot:

  Omega total.

- hier:

  Omega hierarchical.

- sub:

  Omega subscale.

- H:

  H index.

- ECV:

  Explained common variance.

- PUC:

  Percent of uncontaminated correlations.

If found for a `lavaan` single factor solution without multiple groups:
A (named) vector with omega total and (if `add_ind = TRUE`) the H index
for the single factor.

If found for a `lavaan` output from a multiple group analysis: A list
containing the output described above for each group.

## Details

### What this function does

All types of McDonald's omegas (total, hierarchical, and subscale;
McDonald, 1978, 1985, 1999) are calculated for the general factor as
well as for the subscales / group factors (see, e.g., Gignac, 2014;
Rodriguez et al., 2016a, 2016b). Omegas refer to the correlation between
a factor and a unit-weighted composite score and thus the true score
variance in a unit-weighted composite based on the respective
indicators. Omega total is the total true score variance in a composite.
Omega hierarchical is the true score variance in a composite that is
attributable to the general factor, and omega subscale is the true score
variance in a composite attributable to all subscales / group factors
(for the whole scale) or to the specific subscale / group factor (for
subscale composites).

Accordingly, on a subscale row the `hier` column reports the share of
that subscale's composite variance due to the general factor and the
`sub` column the share due to the subscale-specific factor; the latter
corresponds to the omega hierarchical subscale of Rodriguez et al.
(2016a, 2016b).

The H index (also construct reliability or replicability index) is the
correlation between an optimally-weighted composite score and a factor
(Hancock & Mueller, 2001; Rodriguez et al., 2016a, 2016b). It, too, can
be calculated for the whole scale / general factor as well as for the
subscales / group factors. Low values indicate that a latent variable is
not well defined by its indicators.

The ECV (Rodriguez et al., 2016a, 2016b) is the ratio of the variance
explained by the general factor and the variance explained by the
general factor and the group factors.

The PUC (Bonifay et al., 2015; Reise et al., 2013, Rodriguez et al.,
2016a, 2016b) refers to the proportion of correlations in the underlying
correlation matrix that is not contaminated by variance of both the
general factor and the group factors (i.e., correlations between
indicators from different group factors, which reflect only general
factor variance). The higher the PUC, the more similar a general factor
from a multidimensional model will be to the single factor from a
unidimensional model.

### How to use this function

If `model` is a `lavaan` second-order or bifactor solution, only the
name of the general factor from the lavaan model needs to be specified
additionally with the `g_name` argument. It is then determined whether
this general factor is a second-order factor (second-order model with
one second-order factor assumed) or a breadth factor (bifactor model
assumed). Please note that this function only works for second-order
models if they contain no more than one second-order factor. In case of
a second-order solution, a Schmid-Leiman transformation is performed on
the first- and second-order loadings and omega coefficients are obtained
from the transformed (orthogonalized) solution (see
[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
for more information on Schmid-Leiman transformation). There is also the
possibility to enter a `lavaan` single factor solution. In this case,
`g_name` is not needed. Finally, if a solution from a `lavaan` multiple
group analysis is entered, the indices are computed for each group. For
`lavaan` input the composite variances entering the omegas are
model-implied (computed from the fitted loadings and residual
variances), so the coefficients coincide with the observed-score
versions when the model fits perfectly. The type argument is not
evaluated if `model` is of class `lavaan`.

If `model` is of class
[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
or [`psych::schmid()`](https://rdrr.io/pkg/psych/man/schmid.html) only
the `type` and, depending on the type (see below), the `factor_corres`
arguments need to be specified additionally. If model is of class
[`psych::schmid()`](https://rdrr.io/pkg/psych/man/schmid.html) and
`variance = "correlation"` (default), it is recommended to also provide
the original correlation matrix in `cormat` to get more accurate
results. Otherwise, the correlation matrix will be found based on the
pattern matrix and Phi from the
[`psych::schmid()`](https://rdrr.io/pkg/psych/man/schmid.html) output
using the
[`psych::factor.model()`](https://rdrr.io/pkg/psych/man/factor.model.html)
function.

If `model = NULL`, the arguments `type`, `factor_corres` (depending on
the type, see below), `var_names`, `g_load`, `s_load`, and `u2` and
either `cormat` (recommended) or `Phi` and `pattern` need to be
specified. If `Phi` and `pattern` are specified instead of `cormat`, the
correlation matrix is found using the
[`psych::factor.model()`](https://rdrr.io/pkg/psych/man/factor.model.html)
function.

The only difference between `type = "EFAtools"` and `type = "psych"` is
the determination of variable-to-factor correspondences.
`type = "psych"` reproduces the
[`psych::omega()`](https://rdrr.io/pkg/psych/man/omega.html) results,
where variable-to-factor correspondences are found by taking the highest
group factor loading for each variable as the relevant group factor
loading. To do this, `factor_corres` must be left `NULL`.

`variance = "correlation"` gives the observed-variance form of omega,
which reproduces
[`psych::omega()`](https://rdrr.io/pkg/psych/man/omega.html);
`"sums_load"` gives McDonald's model-implied total, which partitions
exactly into omega hierarchical plus omega subscale. The two settings
agree on the whole-scale omega total and omega hierarchical up to model
misfit, and differ mainly in the whole-scale omega subscale, which
counts all group-factor variance under `"sums_load"` but only the
assigned subscale composites under `"correlation"`. On the subscale rows
the two conventions agree when simple structure is well-achieved.

## See also

[`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md)
for the same coefficients in a tidy, long-format result.

## Examples

``` r
# \donttest{
## Use with lavaan outputs
if (requireNamespace("lavaan", quietly = TRUE)) {

# Create and fit bifactor model in lavaan (assume all variables have SDs of 1)
mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
        F2 =~ V7 + V8 + V9 + V10 + V11 + V12
        F3 =~ V13 + V14 + V15 + V16 + V17 + V18
        g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
             V13 + V14 + V15 + V16 + V17 + V18'
fit_bi <- lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
                      sample.nobs = 500, estimator = "ml", orthogonal = TRUE)

# Compute omegas and additional indices for bifactor solution
OMEGA(fit_bi, g_name = "g")

# Compute only omegas
OMEGA(fit_bi, g_name = "g", add_ind = FALSE)

# Create and fit second-order model in lavaan (assume all variables have SDs of 1)
mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
        F2 =~ V7 + V8 + V9 + V10 + V11 + V12
        F3 =~ V13 + V14 + V15 + V16 + V17 + V18
        g =~ F1 + F2 + F3'
fit_ho <- lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
                      sample.nobs = 500, estimator = "ml")

# Compute omegas and additional indices for second-order solution
OMEGA(fit_ho, g_name = "g")
}
#> ℹ The specified general factor is a second-order factor; omegas are computed on
#>   the Schmid-Leiman transformed second-order solution.
#> Omega total, omega hierarchical, omega subscale, H index, explained common
#> variance (ECV), and percent of uncontaminated correlations (PUC) for the
#> general factor (top row) and omegas and H index for the group factors:
#> 
#>      tot  hier   sub     H   ECV   PUC
#> g  0.882 0.764 0.118 0.848 0.684 0.706
#> F1 0.744 0.549 0.195 0.361            
#> F2 0.764 0.504 0.259 0.448            
#> F3 0.768 0.505 0.263 0.454            
# }

## Use with an output from the SL function, with type EFAtools
efa_mod <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
                   estimator = "PAF", rotation = "promax")
sl_mod <- efa_schmid_leiman(efa_mod, estimator = "PAF")

# Indicator-to-factor correspondences from a salience threshold (here: .20):
factor_corres_1 <- sl_mod$sl[, c("F1", "F2", "F3")] >= .2

OMEGA(sl_mod, type = "EFAtools", factor_corres = factor_corres_1)
#> Omega total, omega hierarchical, omega subscale, H index, explained common
#> variance (ECV), and percent of uncontaminated correlations (PUC) for the
#> general factor (top row) and omegas and H index for the group factors:
#> 
#>      tot  hier   sub     H   ECV   PUC
#> g  0.883 0.740 0.125 0.842 0.652 0.706
#> F1 0.769 0.500 0.269 0.463            
#> F2 0.763 0.494 0.270 0.472            
#> F3 0.744 0.519 0.225 0.408            

## Use with an output from the psych::schmid function, with type psych for
## OMEGA
schmid_mod <- psych::schmid(test_models$baseline$cormat, nfactors = 3,
                            n.obs = 500, fm = "pa", rotate = "Promax")
#> Loading required namespace: GPArotation
# Find correlation matrix from phi and pattern matrix from psych::schmid output
OMEGA(schmid_mod, type = "psych")
#> Omega total, omega hierarchical, omega subscale, H index, explained common
#> variance (ECV), and percent of uncontaminated correlations (PUC) for the
#> general factor (top row) and omegas and H index for the group factors:
#> 
#>       tot  hier   sub     H   ECV   PUC
#> g   0.883 0.750 0.122 0.845 0.661 0.706
#> F1* 0.769 0.498 0.272 0.466            
#> F2* 0.765 0.494 0.271 0.473            
#> F3* 0.745 0.543 0.202 0.379            
# Use specified correlation matrix
OMEGA(schmid_mod, type = "psych", cormat = test_models$baseline$cormat)
#> Omega total, omega hierarchical, omega subscale, H index, explained common
#> variance (ECV), and percent of uncontaminated correlations (PUC) for the
#> general factor (top row) and omegas and H index for the group factors:
#> 
#>       tot  hier   sub     H   ECV   PUC
#> g   0.883 0.750 0.122 0.845 0.661 0.706
#> F1* 0.769 0.498 0.272 0.466            
#> F2* 0.764 0.494 0.270 0.473            
#> F3* 0.745 0.543 0.202 0.379            

## Manually specify components (useful if omegas should be computed for a SL
## or bifactor solution found with another program)
## As an example, we extract the elements from an SL output here. This gives
## the same results as in the second example above.

factor_corres <- matrix(c(rep(0, 12), rep(1, 6), rep(0, 6), rep(1, 6),
                        rep(0, 6), rep(1, 6), rep(0, 12)), ncol = 3,
                        byrow = FALSE)

OMEGA(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
      g_load = sl_mod$sl[, "g"], s_load = sl_mod$sl[, c("F1", "F2", "F3")],
      u2 = sl_mod$sl[, "u2"], cormat = test_models$baseline$cormat,
      factor_corres = factor_corres)
#> Omega total, omega hierarchical, omega subscale, H index, explained common
#> variance (ECV), and percent of uncontaminated correlations (PUC) for the
#> general factor (top row) and omegas and H index for the group factors:
#> 
#>     tot  hier   sub     H   ECV   PUC
#> g 0.883 0.740 0.125 0.842 0.652 0.706
#> 1 0.769 0.500 0.269 0.463            
#> 2 0.763 0.494 0.270 0.472            
#> 3 0.744 0.519 0.225 0.408            
```
