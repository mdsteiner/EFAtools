# Reliability and common-variance coefficients for a factor solution

Computes model-based reliability coefficients (McDonald's omega total,
hierarchical, and subscale; standardized Cronbach's alpha; and the H
index) together with the bifactor common-variance indices (explained
common variance, ECV, and percent of uncontaminated correlations, PUC)
for the general factor and each group factor, and returns them as a
tidy, long-format table. The coefficients can be obtained from a
Schmid-Leiman solution
([`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
or [`psych::schmid()`](https://rdrr.io/pkg/psych/man/schmid.html)), an
oblique
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
(correlated-factors) solution, a `lavaan` single-factor, second-order,
or bifactor fit, a raw bifactor loading matrix, or manually supplied
components.

## Usage

``` r
efa_reliability(
  model = NULL,
  coefficients = NULL,
  g_name = "g",
  group_names = NULL,
  factor_map = NULL,
  variance = c("correlation", "sums_load"),
  var_names = NULL,
  fac_names = NULL,
  g_load = NULL,
  s_load = NULL,
  u2 = NULL,
  cormat = NULL,
  pattern = NULL,
  Phi = NULL
)
```

## Source

McDonald, R. P. (1978). Generalizability in factorable domains: Domain
validity and generalizability. Educational and Psychological
Measurement, 38, 75-79.

McDonald, R. P. (1985). Factor analysis and related methods. Hillsdale,
NJ: Erlbaum.

McDonald, R. P. (1999). Test theory: A unified treatment. Mahwah, NJ:
Erlbaum.

Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
tests. Psychometrika, 16, 297-334.

Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
reliability within latent variable systems. In R. Cudeck, S. du Toit, &
D. Sörbom (Eds.), Structural equation modeling: Present and future (pp.
195-216). Lincolnwood, IL: Scientific Software International.

Bonifay, W. E., Reise, S. P., Scheines, R., & Meijer, R. R. (2015). When
are multidimensional data unidimensional enough for structural equation
modeling? An evaluation of the DETECT multidimensionality index.
Structural Equation Modeling, 22, 504-516.

Reise, S. P., Scheines, R., Widaman, K. F., & Haviland, M. G. (2013).
Multidimensionality and structural coefficient bias in structural
equation modeling: A bifactor perspective. Educational and Psychological
Measurement, 73, 5-26.

Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016a). Applying
bifactor statistical indices in the evaluation of psychological
measures. Journal of Personality Assessment, 98, 223-237.

Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016b). Evaluating
bifactor models: Calculating and interpreting statistical indices.
Psychological Methods, 21, 137-150.

## Arguments

- model:

  an
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md),
  `schmid`
  ([`psych::schmid()`](https://rdrr.io/pkg/psych/man/schmid.html)),
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  (oblique), or `lavaan` object; a raw bifactor loading matrix (general
  factor first); or `NULL` to supply the components manually via
  `g_load`, `s_load`, `u2`, and `var_names`.

- coefficients:

  character. An optional subset of the coefficients to return, any of
  `"omega_total"`, `"omega_hierarchical"`, `"omega_subscale"`,
  `"alpha"`, `"H"`, `"ECV"`, and `"PUC"`. Default `NULL` returns all of
  them.

- g_name:

  character. The name of the general factor in the `lavaan` solution.
  Only needed for a `lavaan` second-order or bifactor fit. Default is
  `"g"`.

- group_names:

  character. An optional vector of group names for a `lavaan`
  multiple-group fit. Its length must match the number of groups.

- factor_map:

  matrix. A logical or 0/1 matrix indicating which variable corresponds
  to which group factor, with the same dimensions as the group loading
  matrix (cross-loadings are allowed). Its columns are matched to the
  group factors by position, so a map in a different factor order than
  the solution yields well-formed but meaningless subscale coefficients;
  a map whose assigned items hardly load on the factor they are assigned
  to while loading on another one is warned about. If `NULL` (default),
  each variable is assigned to the group factor on which it loads most
  strongly. Not used for `lavaan` input.

- variance:

  character. The total-variance denominator for the coefficients:
  `"correlation"` (default) takes it from the correlation matrix (the
  observed-score omega, as in
  [`psych::omega()`](https://rdrr.io/pkg/psych/man/omega.html));
  `"sums_load"` uses the model-implied composite variance from the
  squared loading sums and the uniquenesses, which needs no correlation
  matrix and so is the way to score a bare loading matrix or manual
  components given without one. Some inputs fix the convention: `lavaan`
  is always model-implied and an oblique
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  is always correlation-based.

- var_names:

  character. Subtest names in the row order of the loadings. Only needed
  when `model` is `NULL`.

- fac_names:

  character. An optional vector of group-factor names in the column
  order of the loadings. Taken from the input if `NULL`.

- g_load:

  numeric. General-factor loadings. Only needed when `model` is `NULL`.

- s_load:

  matrix. Group-factor loadings. Only needed when `model` is `NULL`.

- u2:

  numeric. Uniquenesses. Only needed when `model` is `NULL` (or to
  override the communality-based default for a raw bifactor matrix).

- cormat:

  matrix or data.frame. A correlation matrix used when
  `variance = "correlation"`. If `NULL`, it is taken from the input
  object or reconstructed from `pattern` and `Phi` where possible.

- pattern:

  matrix. Pattern coefficients from an oblique solution, used with `Phi`
  to reconstruct `cormat` when `model` is `NULL` and no `cormat` is
  given.

- Phi:

  matrix. Factor intercorrelations from an oblique solution, used with
  `pattern`.

## Value

An object of class `efa_reliability`: a long-format data frame with one
row per computed coefficient, with columns

- coefficient:

  the coefficient name (e.g. `"omega_total"`).

- level:

  `"general"` for the general-factor row, `"group"` otherwise.

- factor:

  the factor label (`"g"` for the general factor).

- group:

  the group label, or `NA` for a single ungrouped solution.

- value:

  the coefficient value.

Structurally undefined cells (for example ECV and PUC on a group factor)
are omitted. The object carries a `settings` attribute (the
total-variance convention used, and whether the solution has a general
factor) and a `kind` attribute tagging each coefficient as a reliability
coefficient or a common-variance index, and has a
[`print()`](https://rdrr.io/r/base/print.html) method.

## Details

### Coefficients

The reliability coefficients are McDonald's omegas (McDonald, 1978,
1985, 1999), standardized Cronbach's alpha (Cronbach, 1951), and the H
index (construct replicability; Hancock & Mueller, 2001). The omegas
give the share of true score variance in a unit-weighted composite:
omega total the share due to all factors together, omega hierarchical
the share due to the general factor only, and omega subscale the share
due to the group factors (for the whole scale) or to the specific group
factor (for a subscale composite). Alpha is computed from the
standardized correlation matrix. The H index is the reliability of an
optimally weighted composite; low values indicate a factor that is not
well defined by its indicators.

The common-variance indices ECV and PUC (Bonifay et al., 2015; Reise et
al., 2013; Rodriguez et al., 2016a, 2016b) describe the general factor
and so are reported for the general factor only. The ECV is the share of
the common variance that is explained by the general factor. The PUC is
the proportion of correlations that reflect general-factor variance
alone (those between indicators of different group factors); the higher
it is, the more similar the general factor is to the single factor of a
unidimensional model.

### Input

An oblique
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
object is scored as the correlated-factors model it is (having no
general factor, it omits the bifactor indices – omega hierarchical, ECV,
and PUC), and a bare loading matrix is read as a raw bifactor solution
(general factor in the first column). For a correlated-factors
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
solution `variance` is always `"correlation"`. The indicator-to-factor
correspondences come from `factor_map` when it is supplied; otherwise
each variable is assigned to the group factor on which it loads most
strongly. For `lavaan` input the composite variances are model-implied
(`variance` is not used), and the coefficients are computed per group.

## See also

[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
for the solution these are computed from, and
[`OMEGA()`](https://mdsteiner.github.io/EFAtools/reference/OMEGA.md),
the superseded function that returns these same coefficients in a wide,
per-factor layout.

Other reliability coefficients:
[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md),
[`print.efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_reliability.md)

## Examples

``` r
## From an oblique EFA (correlated-factors) solution. With no factor_map, each
## item is auto-assigned to its highest-loading factor.
efa_mod <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
                   estimator = "PAF", rotation = "promax")
efa_reliability(efa_mod)
#> 
#> Total variance from the correlation matrix.
#> 
#> Correlated-factors solution: with no general factor, each factor's subscale
#> omega equals its total omega.
#> 
#> ── Reliability coefficients ────────────────────────────────────────────────────
#> 
#>      tot   sub  alpha    H
#> g   .883         .868
#> F1  .734  .734   .768  .760
#> F2  .680  .680   .763  .753
#> F3  .667  .667   .743  .738

## From a Schmid-Leiman solution, with an explicit indicator-to-factor map.
sl_mod <- efa_schmid_leiman(efa_mod, estimator = "PAF")
fc <- sl_mod$sl[, c("F1", "F2", "F3")] >= .2
efa_reliability(sl_mod, factor_map = fc)
#> 
#> Total variance from the correlation matrix.
#> 
#> ── Reliability coefficients ────────────────────────────────────────────────────
#> 
#>      tot  hier   sub  alpha    H
#> g   .883  .740  .125   .868  .842
#> F1  .769  .500  .269   .768  .463
#> F2  .763  .494  .270   .763  .472
#> F3  .744  .519  .225   .743  .408
#> 
#> ── Common-variance indices ─────────────────────────────────────────────────────
#> 
#>     ECV   PUC
#> g  .652  .706

## Request a subset of the coefficients only.
efa_reliability(sl_mod, factor_map = fc,
                coefficients = c("omega_total", "alpha"))
#> 
#> Total variance from the correlation matrix.
#> 
#> ── Reliability coefficients ────────────────────────────────────────────────────
#> 
#>      tot  alpha
#> g   .883   .868
#> F1  .769   .768
#> F2  .763   .763
#> F3  .744   .743

## From a lavaan bifactor solution.
# \donttest{
if (requireNamespace("lavaan", quietly = TRUE)) {
mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
        F2 =~ V7 + V8 + V9 + V10 + V11 + V12
        F3 =~ V13 + V14 + V15 + V16 + V17 + V18
        g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
             V13 + V14 + V15 + V16 + V17 + V18'
fit <- lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
                   sample.nobs = 500, estimator = "ml", orthogonal = TRUE)
efa_reliability(fit, g_name = "g")
}
#> 
#> Total variance from the model-implied composite variance.
#> 
#> ── Reliability coefficients ────────────────────────────────────────────────────
#> 
#>      tot  hier   sub  alpha    H
#> g   .883  .765  .118   .868  .849
#> F1  .748  .570  .177   .742  .379
#> F2  .767  .501  .265   .762  .482
#> F3  .769  .494  .274   .768  .473
#> 
#> ── Common-variance indices ─────────────────────────────────────────────────────
#> 
#>     ECV   PUC
#> g  .672  .706
# }
```
