# Reliability and common-variance coefficients for a factor solution

Computes model-based reliability coefficients for a factor solution:
McDonald's omega (total, hierarchical, and subscale), standardized
Cronbach's alpha, and the H index. For a bifactor solution, it also
computes two common-variance indices, ECV and PUC, for the general
factor. The result is a tidy, long-format table with one row per
coefficient.

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

  an
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md),
  `schmid`
  ([`psych::schmid()`](https://rdrr.io/pkg/psych/man/schmid.html)),
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  (oblique), or `lavaan` object; a raw bifactor loading matrix (general
  factor first), the loading table of an
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  solution, or the pattern matrix of an oblique solution together with
  its `Phi`; or `NULL` to supply the components manually via `g_load`,
  `s_load`, `u2`, and `var_names`.

- coefficients:

  character. An optional subset of the coefficients to return, any of
  `"omega_total"`, `"omega_hierarchical"`, `"omega_subscale"`,
  `"alpha"`, `"H"`, `"ECV"`, and `"PUC"`. Default `NULL` returns all of
  them.

- g_name:

  character. The name of the general factor in the `lavaan` solution.
  Only needed for a `lavaan` second-order or bifactor fit. A fit in
  which every variable loads on a single factor has no general factor,
  and reads none. No other input is affected by it. Default is `"g"`.

- group_names:

  character. An optional vector of group names for a `lavaan`
  multiple-group fit. Its length must match the number of groups. Not
  used for any other input – including a single-group `lavaan` fit,
  since every one of those is scored as a single ungrouped solution with
  no group label.

- factor_map:

  matrix. A logical or 0/1 matrix indicating which variable belongs to
  which group factor, with the same dimensions as the group loading
  matrix (cross-loadings are allowed). Match its columns to the group
  factors by position – a map given in a different factor order than the
  solution still runs, but produces meaningless subscale coefficients.
  The function warns if a mapped item loads weakly on its assigned
  factor and more strongly on another one. If `NULL` (default), each
  variable is assigned to the group factor on which it loads most
  strongly. Not used for `lavaan` input.

- variance:

  character. The total-variance denominator for the coefficients.
  `"correlation"` (default) takes each composite's variance from the
  correlation matrix, giving the observed-score omega. `"sums_load"`
  uses the model-implied composite variance from the loadings and the
  uniquenesses instead; it needs no correlation matrix, so it is the way
  to score a bare loading matrix or manual components given without one.
  `lavaan` input fixes the convention: its composite variances are
  always model-implied, and include any freed residual covariance as
  well as the residual variances. Neither convention changes the metric
  the coefficients are on – with polychoric or tetrachoric correlations,
  or an `ordered` `lavaan` fit, both describe the latent-response
  composite rather than the ordinal sum score (see Details).

- var_names:

  character. Subtest names in the row order of the loadings. Only needed
  when `model` is `NULL`.

- fac_names:

  character. An optional vector of group-factor names in the column
  order of the loadings. Taken from the input if `NULL`. A single-factor
  solution has no group factors; `fac_names` then labels its one factor
  instead. If neither `fac_names` nor the solution names that factor, it
  is labelled `"F1"`. Not used for `lavaan` input, whose factor labels
  come from the model syntax.

- g_load:

  numeric. General-factor loadings. Only needed when `model` is `NULL`.

- s_load:

  matrix. Group-factor loadings. Only needed when `model` is `NULL`.

- u2:

  numeric. Uniquenesses. Only needed when `model` is `NULL`, or to
  override the communality-based default for a loading matrix. Under
  `variance = "correlation"`, the coefficients follow from the loadings
  and the correlation matrix alone, so `u2` only enters the check for an
  improper solution. Under `"sums_load"`, it is part of every
  composite's variance.

- cormat:

  matrix or data.frame. A correlation matrix used when
  `variance = "correlation"`. It must hold the same variables as the
  solution: named variables in a different order are reordered to match;
  a different set of variables, or a different number of them, is an
  error. The matrix must use the solution's own variable names – for
  manually supplied components, that means the row names of `s_load`,
  not `var_names`, which only labels the output. Without row names, the
  variables are matched by position, so supply the matrix in the row
  order of the loadings. If `NULL`, it is taken from the input object,
  or reconstructed from `pattern` and `Phi` where possible. Supply the
  matrix on the same metric as the loadings: a polychoric or tetrachoric
  matrix keeps the coefficients on the latent-response metric, while a
  Pearson matrix given with loadings fitted to a polychoric one mixes
  the two (see Details).

- pattern:

  matrix. Pattern coefficients from a separate oblique solution, used
  with `Phi` to reconstruct a correlation matrix when `model` is `NULL`.
  Supply it for a Schmid-Leiman input, whose `s_load` holds the
  orthogonalized group loadings rather than the oblique ones. It is an
  alternative to `cormat`, not a supplement: giving both is an error.

- Phi:

  matrix. Factor intercorrelations. `NULL` (default) means uncorrelated
  group factors, as a Schmid-Leiman or bifactor solution has. It cannot
  be combined with a general factor: supplied together with a non-zero
  `g_load`, or with the loading table of an
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  solution, it is an error.

  Supply it together with `s_load` (manually supplied components, with
  `g_load` zero throughout), or with a loading matrix of two or more
  factors in `model`, to score the input as a correlated-factors
  solution instead.

  Without `Phi`, a loading matrix in `model` is read as a bifactor
  solution (general factor first). The one exception is a matrix that
  still carries the loading class
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  gives its output – a class that subsetting a matrix or reordering its
  rows drops. Such a matrix is rejected instead; supply `Phi` to score
  it as a correlated-factors solution.

  With `Phi`, a matrix in `model` that carries no loading class is read
  as a correlated-factors solution, and a warning names this reading
  (drop `Phi` to read it as a bifactor matrix instead).

  Given together with `pattern` instead, see `pattern` above. With a
  fitted solution already in `model`, `Phi` is ignored, with a warning,
  since that solution carries its own factor intercorrelations.

## Value

An object of class `efa_reliability`: a long-format data frame with one
row per computed coefficient, with columns

- coefficient:

  the coefficient name (e.g. `"omega_total"`).

- level:

  `"general"` for the general-factor row, `"total"` for the whole-scale
  row of a correlated-factors solution, and `"group"` for every other
  row. A single-factor solution has only the general-factor row.

- factor:

  the factor label: `"g"` for the general factor of a solution with
  group factors. A single-factor solution's one factor takes its own
  name, or `"F1"` if neither the input nor `fac_names` names it. A
  correlated-factors solution has no general factor, so its first row is
  labelled `"total"` instead – it describes the composite of every
  variable, not a factor of the model.

- group:

  the group label, or `NA` for a single ungrouped solution.

- value:

  the coefficient value.

Structurally undefined cells (for example, ECV and PUC on a group
factor) are omitted. The object also carries a `settings` attribute (the
total-variance convention used, and whether the solution has a general
factor) and a `kind` attribute tagging each coefficient as a reliability
coefficient or a common-variance index. It has a
[`print.efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_reliability.md)
method.

## Details

The function reads many kinds of input: a Schmid-Leiman solution
([`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
or [`psych::schmid()`](https://rdrr.io/pkg/psych/man/schmid.html)), an
oblique
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
(correlated-factors) solution, a `lavaan` fit (single-factor,
correlated-factors, second-order, or bifactor), a raw bifactor loading
matrix, an oblique pattern matrix given with its factor
intercorrelations, or manually supplied components.

### Coefficients

The reliability coefficients are McDonald's omegas (McDonald, 1978,
1985, 1999; Zinbarg et al., 2005, 2006, for omega hierarchical
specifically), standardized Cronbach's alpha (Cronbach, 1951), and the H
index (construct replicability; Hancock & Mueller, 2001).

The omegas give the share of true score variance in a unit-weighted
composite. Omega total is the share due to all factors together. Omega
hierarchical is the share due to the general factor only. Omega subscale
is the share due to the group factors: for the whole scale, or for one
specific group factor in a subscale composite.

Alpha is the standardized coefficient, computed from the correlation
matrix of the items. Where no such matrix is available – for `lavaan`
input, and for components supplied without one – alpha is computed from
the model-implied correlation matrix instead, so it then reflects the
fitted model rather than the raw data.

The H index is the reliability of an optimally weighted composite. A low
value means the factor is not well defined by its indicators.

All of these coefficients describe the raw sum of the variables as
fitted, without reverse-coding. If some items are keyed in the opposite
direction – for example, a reverse-worded item that was not
reverse-scored – they lower that sum's true-score variance. The
coefficients then look poor even though the model fits well. The
function warns when it detects this. Reverse-code such items before
fitting the solution (Flora, 2020).

The sum these coefficients describe is not always the sum of the raw
answers. Polychoric and tetrachoric correlations describe the continuous
latent responses assumed to underlie ordinal answers; a `lavaan` fit
that declares its indicators `ordered` does the same. In that case the
loadings, the uniquenesses, the correlation matrix, and the coefficients
are all on that latent-response metric. They give the reliability of the
unit-weighted sum of the latent responses, which is not observed – not
the reliability of the ordinal sum score the user actually computes from
the answers. The two can differ substantially, especially where the
answers use few categories or are strongly skewed. Green and Yang (2009)
give an omega for the ordinal sum score itself, computed from the fitted
model and its thresholds; this package does not compute it. Pearson
correlations raise no such distinction, because their metric is the
answers as scored.

Omega total is lower when a solution reproduces the observed
correlations poorly, because residual covariance does not count as true
score. [`psych::omega()`](https://rdrr.io/pkg/psych/man/omega.html)
computes the whole-scale omega total differently: residually, from the
observed total-score variance. The two agree when the model reproduces
the correlations exactly, and diverge otherwise.

The three coefficients are not generally additive. Omega total need not
equal omega hierarchical plus omega subscale, except on the whole-scale
row under `variance = "sums_load"`.

A single-factor solution is scored as such on every input route, and
returns omega total, alpha, and the H index. Alpha assumes essentially
tau-equivalent items, an assumption nested within a one-factor model –
so a single factor is the case where reporting alpha is defensible, not
merely possible.

The other coefficients are omitted because a single factor does not
define them: omega subscale is the variance due to the group factors,
and there are none; omega hierarchical would equal omega total, since
the one factor accounts for all common variance; and ECV and PUC would
both be 1 by construction, which reflects the number of factors in the
model rather than evidence of unidimensionality.

Each coefficient answers a different question:

- Omega hierarchical: can the total score be read as a measure of one
  construct?

- Omega subscale: does a subscale score add anything beyond the general
  factor?

- The H index: is a factor well defined by its indicators?

- ECV together with PUC: is a unidimensional model defensible?

Alpha assumes essentially tau-equivalent items. Factor analysis instead
yields congeneric solutions, for which alpha is only a lower bound. For
a multidimensional scale, alpha is rarely the coefficient to report
(Gignac, 2014).

Composite reliability and average variance extracted (AVE) are not among
the coefficients this function computes. Use
[`semTools::compRelSEM()`](https://rdrr.io/pkg/semTools/man/compRelSEM.html)
and [`semTools::AVE()`](https://rdrr.io/pkg/semTools/man/AVE.html) to
compute them from a `lavaan` fit.

The common-variance indices ECV and PUC (Bonifay et al., 2015; Reise et
al., 2013; Rodriguez et al., 2016a, 2016b) describe the general factor,
so they are reported for the general factor only. ECV is the share of
the common variance explained by the general factor. PUC is the
proportion of correlations that reflect general-factor variance alone –
correlations between indicators of different group factors. The higher
the PUC, the more the general factor resembles the single factor of a
unidimensional model.

### Input

`efa_reliability()` reads several kinds of input, illustrated in the
examples below.

- **An oblique
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  solution** is scored as the correlated-factors model it is. It has no
  general factor, so it omits the bifactor indices (omega hierarchical,
  ECV, and PUC). Its whole-scale row is labelled `"total"` rather than
  `"g"`.

- **An
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  solution** – or a `schmid` object from
  [`psych::schmid()`](https://rdrr.io/pkg/psych/man/schmid.html), or a
  raw bifactor loading matrix with the general factor in its first
  column – is scored as a bifactor solution, with the whole-scale row
  labelled `"g"`. Indicator-to-factor correspondences come from
  `factor_map` if supplied (see `factor_map` below). Without one, an
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  or [`psych::schmid()`](https://rdrr.io/pkg/psych/man/schmid.html)
  solution is mapped by each variable's strongest group loading. A raw
  bifactor matrix instead defaults to its exact zero pattern – supply
  `factor_map` explicitly if that matrix has no exact zeros (e.g. an
  estimated rather than a target matrix).

- **A `lavaan` fit** – single-factor, correlated-factors, second-order,
  or bifactor – is scored per its structure. The structure is detected
  automatically, not taken from `g_name`. The general-factor
  coefficients need the latent factors to be uncorrelated: fit a
  bifactor model with `orthogonal = TRUE` (`lavaan`'s default leaves the
  factors correlated), and leave a second-order model's first-order
  factor covariances at zero. A fit whose factors correlate is rejected
  rather than scored as though they did not. `variance` is not used for
  `lavaan` input: its composite variances are always model-implied,
  computed separately per group for a multiple-group fit.

- **A pattern matrix with its `Phi`** – for example, the loadings and
  `Phi` of an oblique
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  solution – is scored as the correlated-factors solution it represents.
  Uniquenesses are derived from the loadings under `Phi`, unless `u2` is
  given. Supply the pattern matrix, not the structure matrix – the
  structure matrix would be read as a pattern too, giving the wrong
  result. Unlike a fitted solution, this pair carries no correlation
  matrix. `variance = "correlation"` then needs one, supplied in
  `cormat`, and errors without it. `variance = "sums_load"` needs none.

A one-factor solution – a one-column loading matrix, a single-factor
`lavaan` fit, a single-factor
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
solution, or single-factor components – is scored the same way
regardless of route (see Coefficients above for what it returns). Its
row is labelled with the factor's own name, with `fac_names` if
supplied, or with `"F1"` if neither names it; it is never labelled
`"g"`.

Manually supplied components (`g_load`, `s_load`, `u2`, `var_names`,
`Phi`) follow the same reading rules as the matrix routes above. Do not
pass a correlation matrix as `s_load` (or as `model`, for the matrix
routes) – supply it as `cormat` instead.

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
#> Correlated-factors solution: a factor's total omega counts the true score
#> variance its composite receives from every factor, through its cross-loadings
#> and any factor correlations; its subscale omega counts only that factor's own
#> contribution.
#> 
#> ── Reliability coefficients ────────────────────────────────────────────────────
#> 
#>         tot   sub  alpha    H
#> total  .883         .868
#> F1     .769  .734   .768  .760
#> F2     .765  .680   .763  .753
#> F3     .745  .667   .743  .738

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
#> Total variance from the correlation matrix.
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
#> Total variance from the correlation matrix.
#> 
#> Correlated-factors solution: a factor's total omega counts the true score
#> variance its composite receives from every factor, through its cross-loadings
#> and any factor correlations; its subscale omega counts only that factor's own
#> contribution.
#> 
#> ── Reliability coefficients ────────────────────────────────────────────────────
#> 
#>         tot   sub  alpha    H
#> total  .883         .868
#> F1     .769  .734   .768  .760
#> F2     .765  .680   .763  .753
#> F3     .745  .667   .743  .738

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
print(efa_reliability(fit, g_name = "g"))

## No general factor: omega hierarchical, ECV, and PUC are omitted.
fit_cf <- lavaan::cfa(mod_cf, sample.cov = test_models$baseline$cormat,
                      sample.nobs = 500, estimator = "ml")
efa_reliability(fit_cf)
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
#> ℹ Each variable loads on one factor only, so the fit is scored as a
#>   correlated-factors solution; omega hierarchical, ECV, and PUC are not defined
#>   for it and are omitted.
#> 
#> Total variance from the model-implied composite variance.
#> 
#> Correlated-factors solution: with no general factor, each factor's subscale
#> omega equals its total omega.
#> 
#> ── Reliability coefficients ────────────────────────────────────────────────────
#> 
#>         tot   sub  alpha    H
#> total  .882         .868
#> F1     .744  .744   .743  .747
#> F2     .764  .764   .763  .765
#> F3     .768  .768   .768  .769
# }
```
