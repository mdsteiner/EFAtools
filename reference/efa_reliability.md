# Reliability coefficients for a factor solution

Computes model-based reliability coefficients for a factor solution:
McDonald's omega (total, hierarchical, and subscale), standardized
Cronbach's alpha, and the H index. For a solution with a general factor,
it also computes the common-variance indices ECV and PUC. The result is
a long-format table with one row per coefficient.

## Usage

``` r
efa_reliability(
  model = NULL,
  coefficients = NULL,
  g_name = "g",
  group_names = NULL,
  factor_map = NULL,
  variance = c("correlation", "sums_load"),
  var_names = lifecycle::deprecated(),
  fac_names = NULL,
  g_load = lifecycle::deprecated(),
  s_load = lifecycle::deprecated(),
  u2 = NULL,
  cormat = NULL,
  pattern = lifecycle::deprecated(),
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

Gignac, G. E. (2014). On the inappropriateness of using items to
calculate total scale score reliability via coefficient alpha for
multidimensional scales. European Journal of Psychological Assessment,
30, 130-139.

Flora, D. B. (2020). Your coefficient alpha is probably wrong, but which
coefficient omega is right? A tutorial on using R to obtain better
reliability estimates. Advances in Methods and Practices in
Psychological Science, 3, 484-501.

Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
structural equation modeling: An alternative to coefficient alpha.
Psychometrika, 74, 155-167.

Zinbarg, R. E., Revelle, W., Yovel, I., & Li, W. (2005). Cronbach's
alpha, Revelle's beta, and McDonald's omega H: Their relations with each
other and two alternative conceptualizations of reliability.
Psychometrika, 70, 123-133.

Zinbarg, R. E., Yovel, I., Revelle, W., & McDonald, R. P. (2006).
Estimating generalizability to a latent variable common to all of a
scale's indicators: A comparison of estimators for omega H. Applied
Psychological Measurement, 30, 121-144.

Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
reliability within latent variable systems. In R. Cudeck, S. du Toit, &
D. Sörbom (Eds.), Structural equation modeling: Present and future - A
Festschrift in honor of Karl Jöreskog (pp. 195-216). Lincolnwood, IL:
Scientific Software International.

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

  An
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md),
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md),
  [`psych::schmid()`](https://rdrr.io/pkg/psych/man/schmid.html), or
  `lavaan` solution, or a loading matrix (see Input).

- coefficients:

  character. A subset of the coefficients to return: any of
  `"omega_total"`, `"omega_hierarchical"`, `"omega_subscale"`,
  `"alpha"`, `"H"`, `"ECV"`, and `"PUC"`. The default `NULL` returns all
  that the solution defines.

- g_name:

  character. The name of the general factor of a `lavaan` bifactor or
  second-order fit. Default is `"g"`.

- group_names:

  character. Labels for the groups of a `lavaan` multiple-group fit, one
  per group.

- factor_map:

  matrix. A logical or 0/1 matrix, variables by group factors, that
  assigns the variables to the group factors. A variable can belong to
  more than one factor. The columns must be in the column order of the
  solution. If `NULL` (default), see Input. Not used for `lavaan` input.

- variance:

  character. The variance of each composite: `"correlation"` (default)
  takes it from the correlation matrix; `"sums_load"` uses the
  model-implied variance and needs no correlation matrix. Not used for
  `lavaan` input (see Details).

- var_names, g_load, s_load, pattern:

  **\[deprecated\]** Give the loadings as one matrix in `model` instead.

- fac_names:

  character. Labels for the group factors, in the column order of the
  loadings; for a single-factor solution, the label of its factor. Not
  used for `lavaan` input.

- u2:

  numeric. Uniquenesses, one per variable, to replace the ones derived
  from a loading matrix in `model`.

- cormat:

  matrix or data.frame. The correlation matrix of the variables, on the
  metric of the loadings, for `variance = "correlation"`. Named
  variables are matched by name, else by position. If `NULL`, it is
  taken from the solution where available.

- Phi:

  matrix. The factor intercorrelations of a pattern matrix in `model`.
  Not used for a fitted solution, and refused for the loading table of
  an
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  solution.

## Value

An object of class `efa_reliability`: a data frame with one row per
coefficient and the columns

- coefficient:

  the coefficient name.

- level:

  `"general"` for the general factor (or the factor of a single-factor
  solution), `"total"` for the whole scale of a correlated-factors
  solution, and `"group"` for a group factor.

- factor:

  the factor label: `"g"` for the general factor, `"total"` for the
  whole scale of a correlated-factors solution.

- group:

  the group of a `lavaan` multiple-group fit, else `NA`.

- value:

  the coefficient value.

Coefficients that the solution does not define are omitted. The
attributes are `settings` (the variance used, and whether the solution
has no general factor), `kind` (reliability coefficient or
common-variance index), and `factor_map` (the assignment of the
variables to the group factors that was used: a logical matrix, a list
of them for a multiple-group fit, or `NULL` for a single factor).

## Details

### Coefficients

Each row of the result describes one composite, that is, the
unit-weighted sum of a set of variables: the whole scale, or the
variables assigned to one group factor. Each omega is the share of the
composite's variance that the model gives to the factors (McDonald,
1978, 1985, 1999; Zinbarg et al., 2005, 2006):

- **Omega total** (`omega_total`): the share due to all factors. On a
  group-factor row, it is the omega S of Rodriguez et al. (2016a,
  2016b).

- **Omega hierarchical** (`omega_hierarchical`): the share due to the
  general factor. On the whole-scale row, it tells whether the total
  score measures one construct.

- **Omega subscale** (`omega_subscale`): the share due to the group
  factors; on a group-factor row, the share due to its own group factor.
  There, it is the omega HS of Rodriguez et al. (2016a, 2016b): what the
  subscale score measures in addition to the general factor.

On the whole-scale row, omega total is omega hierarchical plus omega
subscale. On a group-factor row, omega total can be larger, because
cross-loadings on other group factors add true-score variance.

The other coefficients are:

- **Alpha** (`alpha`): standardized Cronbach's alpha (Cronbach, 1951).
  It assumes essentially tau-equivalent variables and is a lower bound
  otherwise. For a multidimensional scale, it is rarely the coefficient
  to report (Gignac, 2014).

- **H index** (`H`): the reliability of an optimally weighted composite
  (Hancock & Mueller, 2001). A low value shows that the variables do not
  define the factor well.

- **ECV** and **PUC** (general factor only): the share of the common
  variance due to the general factor, and the proportion of correlations
  between variables of different group factors (Bonifay et al., 2015;
  Reise et al., 2013; Rodriguez et al., 2016a, 2016b). Together, they
  tell whether a unidimensional model is acceptable.

A solution without a general factor (a correlated-factors solution)
gives omega total, alpha, and the H index of each factor, and omega
total and alpha of the whole scale, whose row is labelled `"total"`.
Omega subscale is not given: it is the share in addition to a general
factor, and when the factors correlate, it is not a share of the
variance. A single-factor solution gives omega total, alpha, and the H
index. Use
[`semTools::compRelSEM()`](https://rdrr.io/pkg/semTools/man/compRelSEM.html)
for composite reliability and
[`semTools::AVE()`](https://rdrr.io/pkg/semTools/man/AVE.html) for the
average variance extracted of a `lavaan` fit.

### Total variance

With `variance = "correlation"`, the variance of a composite is the sum
of its block of the correlation matrix, so residual correlations that
the model does not reproduce change omega total. With `"sums_load"`, it
is the model-implied variance. `lavaan` input always uses the
model-implied variance, which includes freed residual covariances.
[`psych::omega()`](https://rdrr.io/pkg/psych/man/omega.html) gives other
values: it computes the whole-scale omega total from the residual
variances, and it leaves cross-loadings out of the whole-scale omega
subscale and of the omega total of a subscale.

### Reverse-keyed and ordinal variables

The coefficients describe the raw sum of the variables as fitted.
Variables keyed in the opposite direction decrease the true-score
variance of that sum, and the function warns when it finds such
variables. Reverse-code them before you fit the solution (Flora, 2020).

With polychoric or tetrachoric correlations, or a `lavaan` fit with
`ordered` variables, the coefficients describe the sum of the latent
continuous responses, not the sum of the ordinal answers. The two can be
very different when the answers have few categories or are skewed. Green
and Yang (2009) give an omega for the ordinal sum score; this function
does not compute it.

### Input

- **An
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  solution.** An oblique solution is a correlated-factors solution. An
  orthogonal solution is scored in the same way, with factor
  correlations of 0. A `"bifactorT"` solution is a bifactor solution,
  whose general factor is the factor that loads on every variable (the
  column whose smallest absolute loading is largest). A `"bifactorQ"`
  solution, whose group factors are correlated, and an unrotated
  solution with two or more factors are refused.

- **An
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  or [`psych::schmid()`](https://rdrr.io/pkg/psych/man/schmid.html)
  solution**, or the loading table of an
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  solution: a bifactor solution.

- **A loading matrix.** Without `Phi`, the first column is the general
  factor and the other columns are the group factors. With `Phi`, the
  matrix is the pattern matrix of a correlated-factors solution (not the
  structure matrix). Give the variable names as row names. A pattern
  matrix needs `cormat` for `variance = "correlation"`.

- **A `lavaan` fit**: single-factor, correlated-factors (cross-loadings
  are permitted), second-order, or bifactor. A fit is a bifactor
  solution if variables load on a factor named `g_name` and on other
  factors. Fit it with `orthogonal = TRUE`, because the coefficients
  need uncorrelated factors. A second-order solution is Schmid-Leiman
  transformed. A multiple-group fit is scored per group.

Without `factor_map`, each variable is assigned to the group factor with
its largest absolute loading, and a variable without a group loading to
none. A loading matrix or a Schmid-Leiman solution in which most
variables have an exact zero group loading (a target matrix, or a
transformed `lavaan` solution) uses its non-zero pattern instead. So
does a `lavaan` fit, but in an ESEM block, each variable is assigned to
the block's factor with its largest loading.

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
## From an oblique EFA (correlated-factors) solution. Each variable is assigned to the
## factor with its largest loading.
efa_mod <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
                   estimator = "PAF", rotation = "promax")
efa_reliability(efa_mod)
#> 
#> ── Reliability coefficients ────────────────────────────────────────────────────
#> 
#>         tot  alpha    H
#> total  .883   .868
#> F1     .769   .768  .760
#> F2     .765   .763  .753
#> F3     .745   .743  .738

## From a Schmid-Leiman solution, with an explicit variable-to-factor map.
sl_mod <- efa_schmid_leiman(efa_mod, estimator = "PAF")
fc <- sl_mod$sl[, c("F1", "F2", "F3")] >= .2
efa_reliability(sl_mod, factor_map = fc)
#> 
#> ── Reliability coefficients ────────────────────────────────────────────────────
#> 
#>      tot  hier   sub  alpha    H
#> g   .883  .740  .143   .868  .842
#> F1  .769  .500  .269   .768  .463
#> F2  .765  .494  .270   .763  .472
#> F3  .745  .519  .225   .743  .408
#> 
#> ── Common-variance indices ─────────────────────────────────────────────────────
#> 
#>     ECV   PUC
#> g  .652  .706

## Request a subset of the coefficients only.
efa_reliability(sl_mod, factor_map = fc,
                coefficients = c("omega_total", "alpha"))
#> 
#> ── Reliability coefficients ────────────────────────────────────────────────────
#> 
#>      tot  alpha
#> g   .883   .868
#> F1  .769   .768
#> F2  .765   .763
#> F3  .745   .743

## From an oblique pattern matrix and its factor intercorrelations. This is
## the same correlated-factors solution, and gives the same coefficients.
efa_reliability(efa_mod$rot_loadings, Phi = efa_mod$Phi,
                cormat = test_models$baseline$cormat)
#> 
#> ── Reliability coefficients ────────────────────────────────────────────────────
#> 
#>         tot  alpha    H
#> total  .883   .868
#> F1     .769   .768  .760
#> F2     .765   .763  .753
#> F3     .745   .743  .738

## From lavaan fits: a bifactor solution, and a correlated-factors one.
# \donttest{
if (requireNamespace("lavaan", quietly = TRUE)) {
mod_cf <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
           F2 =~ V7 + V8 + V9 + V10 + V11 + V12
           F3 =~ V13 + V14 + V15 + V16 + V17 + V18'
mod <- paste(mod_cf, 'g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
                           V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18',
             sep = "\n")
fit <- lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
                   sample.nobs = 500, estimator = "ml", orthogonal = TRUE)
efa_reliability(fit, g_name = "g")

## No general factor: omega hierarchical, omega subscale, ECV, and PUC are omitted.
fit_cf <- lavaan::cfa(mod_cf, sample.cov = test_models$baseline$cormat,
                      sample.nobs = 500, estimator = "ml")
efa_reliability(fit_cf)
}
#> 
#> ── Reliability coefficients ────────────────────────────────────────────────────
#> 
#>         tot  alpha    H
#> total  .882   .868
#> F1     .744   .743  .747
#> F2     .764   .763  .765
#> F3     .768   .768  .769
# }
```
