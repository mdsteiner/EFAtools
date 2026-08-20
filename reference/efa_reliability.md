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
(correlated-factors) solution, a `lavaan` single-factor,
correlated-factors, second-order, or bifactor fit, a raw bifactor
loading matrix, an oblique pattern matrix given with its factor
intercorrelations, or manually supplied components.

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
  Only needed for a `lavaan` second-order or bifactor fit; a fit in
  which every variable loads on a single factor has no general factor
  and reads none. Every other input is silently unaffected by it, and –
  unlike the other arguments a route does not read, which are reported –
  no note says so: `g_name` carries a default rather than `NULL`, so a
  name the caller supplied cannot be told from the one they did not.
  Default is `"g"`.

- group_names:

  character. An optional vector of group names for a `lavaan`
  multiple-group fit. Its length must match the number of groups. Not
  used for any other input, a single-group `lavaan` fit included: every
  one of them is scored as one ungrouped solution, which needs no group
  label.

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
  `"correlation"` (default) takes each composite's variance from the
  correlation matrix, giving the observed-score omega; `"sums_load"`
  uses the model-implied composite variance from the loadings and the
  uniquenesses, which needs no correlation matrix and so is the way to
  score a bare loading matrix or manual components given without one.
  `lavaan` input fixes the convention: its composite variances are
  always model-implied, and count any freed residual covariance as well
  as the residual variances. Neither convention changes the metric the
  coefficients are on: with polychoric or tetrachoric correlations, or
  an `ordered` `lavaan` fit, both describe the latent-response composite
  rather than the ordinal sum score (see Details).

- var_names:

  character. Subtest names in the row order of the loadings. Only needed
  when `model` is `NULL`.

- fac_names:

  character. An optional vector of group-factor names in the column
  order of the loadings. Taken from the input if `NULL`. A single-factor
  solution has no group factors, and `fac_names` labels its one factor
  instead; where neither `fac_names` nor the solution names that factor,
  it is labelled `"F1"`. Not used for `lavaan` input, whose factor
  labels come from the model syntax.

- g_load:

  numeric. General-factor loadings. Only needed when `model` is `NULL`.

- s_load:

  matrix. Group-factor loadings. Only needed when `model` is `NULL`.

- u2:

  numeric. Uniquenesses. Only needed when `model` is `NULL` (or to
  override the communality-based default for a loading matrix). Under
  `variance = "correlation"` the coefficients follow from the loadings
  and the correlation matrix alone, so `u2` then only enters the check
  for an improper solution; under `"sums_load"` it is part of every
  composite's variance.

- cormat:

  matrix or data.frame. A correlation matrix used when
  `variance = "correlation"`. It must hold the variables of the
  solution: named variables in another order are reordered to it, a
  different set of variables or a different number of them is an error.
  The names it is matched against are the solution's own – for manually
  supplied components, the row names of `s_load` rather than
  `var_names`, which only labels the output; without them the variables
  are matched by position, so supply the matrix in the row order of the
  loadings. If `NULL`, it is taken from the input object or
  reconstructed from `pattern` and `Phi` where possible. Supply the
  matrix on the same metric as the loadings: a polychoric or tetrachoric
  matrix keeps the coefficients on the latent-response metric, and a
  Pearson matrix given with loadings fitted to a polychoric one mixes
  the two (see Details).

- pattern:

  matrix. Pattern coefficients from a separate oblique solution, used
  with `Phi` to reconstruct a correlation matrix when `model` is `NULL`.
  Supply it for a Schmid-Leiman input, whose `s_load` holds the
  orthogonalized group loadings rather than the oblique ones. It is an
  alternative to `cormat`, not a supplement: giving both is an error.

- Phi:

  matrix. Factor intercorrelations, describing whichever loadings it is
  supplied with: given together with `pattern` it belongs to that
  solution and only reconstructs `cormat`, given without one it is the
  correlation matrix of the group factors of `s_load`, which makes the
  manually supplied components a correlated-factors solution and enters
  the coefficients. It must then be a proper correlation matrix of as
  many factors as `s_load` has columns, and the components must carry no
  general factor; for a single group factor that leaves only the 1 by 1
  identity, which is accepted and changes nothing, the components being
  scored as the single factor they are either way. `NULL` (default)
  means uncorrelated group factors, as a Schmid-Leiman or bifactor
  solution has. Supplied with a loading matrix of two or more factors in
  `model`, it is that matrix's factor intercorrelations, and makes the
  pair a correlated-factors solution rather than a rejected pattern
  matrix. It does this whether or not the matrix carries the loading
  class
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  gives its loadings, a class that subsetting a matrix or reordering its
  rows drops. A matrix that carries no such class is read as a bifactor
  one when no `Phi` accompanies it, so one that is read as a pattern
  warns that it was. Supplied with the loading table of an
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  solution, it is an error: such a table states that it is a hierarchy,
  whose group factors are uncorrelated by construction, so there is
  nothing for it to describe, and the same combination given through the
  components is refused in the same terms. With a fitted solution in
  `model` it is ignored, with a warning, that solution carrying its own.

## Value

An object of class `efa_reliability`: a long-format data frame with one
row per computed coefficient, with columns

- coefficient:

  the coefficient name (e.g. `"omega_total"`).

- level:

  `"general"` for the general-factor row, `"total"` for the whole-scale
  row of a correlated-factors solution, `"group"` otherwise. A
  single-factor solution has only the general-factor row.

- factor:

  the factor label (`"g"` for the general factor of a solution with
  group factors; the one factor of a single-factor solution takes its
  own name, or `"F1"` where the input gives it none). A
  correlated-factors solution has no general factor, so its first row is
  labelled `"total"`: it describes the composite of every variable
  rather than a factor of the model.

- group:

  the group label, or `NA` for a single ungrouped solution.

- value:

  the coefficient value.

Structurally undefined cells (for example ECV and PUC on a group factor)
are omitted. The object carries a `settings` attribute (the
total-variance convention used, and whether the solution has a general
factor) and a `kind` attribute tagging each coefficient as a reliability
coefficient or a common-variance index, and has a
[`print.efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_reliability.md)
method.

## Details

### Coefficients

The reliability coefficients are McDonald's omegas (McDonald, 1978,
1985, 1999), standardized Cronbach's alpha (Cronbach, 1951), and the H
index (construct replicability; Hancock & Mueller, 2001). The omegas
give the share of true score variance in a unit-weighted composite:
omega total the share due to all factors together, omega hierarchical
the share due to the general factor only, and omega subscale the share
due to the group factors (for the whole scale) or to the specific group
factor (for a subscale composite). Alpha is the standardized
coefficient, computed from the correlation matrix of the items; where
none is available – for `lavaan` input, and for components supplied
without one – it is computed from the model-implied correlation matrix
and so reflects the fitted model rather than the data. The H index is
the reliability of an optimally weighted composite; low values indicate
a factor that is not well defined by its indicators.

All of them describe the raw unit-weighted sum of the variables the
solution was fitted on. A variable keyed against the rest therefore
subtracts from that sum's true score variance instead of adding to it,
and the coefficients collapse – correctly for the sum that was scored,
but in a way that reads as a poor scale. The variables are not reflected
first, which would describe a different sum than the one scored; a
solution whose variables are not all keyed in the same direction is
warned about instead, and such variables should be reverse-coded before
the solution is fitted (Flora, 2020).

The variables that sum is over are not always the answers the
respondents gave. Polychoric and tetrachoric correlations describe the
continuous latent responses that are assumed to underlie ordinal
answers, and so does a `lavaan` fit that declares its indicators
`ordered`. The loadings, the uniquenesses and the correlation matrix are
then all on that latent-response metric, and so are the coefficients.
They give the reliability of the unit-weighted sum of the latent
responses, which is not observed. They do not give the reliability of
the ordinal sum score the user computes from the answers. The two
quantities differ, and they differ most where the answers use few
categories or are strongly skewed. Green and Yang (2009) give the
alternative: an omega for the ordinal sum score itself, computed from
the fitted model and its thresholds. This package does not compute it.
Pearson correlations raise no such distinction, their metric being the
answers as scored.

Every omega total – for the whole scale and for each subscale composite
– is the true score variance the fitted model attributes to that
composite, counting every factor its variables load on, over that
composite's variance. Residual covariance among a composite's variables
is not true score under the model and is therefore not counted, so a
solution that reproduces the observed correlations poorly returns a
lower omega total than one that reproduces them well.
[`psych::omega()`](https://rdrr.io/pkg/psych/man/omega.html) instead
obtains the whole-scale omega total residually, as the observed
total-score variance less the unique variances, and counts only the
general factor and the composite's own group factor for a subscale.
Applied to the same solution, its whole-scale value agrees with the one
here when the model reproduces the correlations exactly, and its
subscale values only when a composite's variables load on no other group
factor; it also fits its own Schmid-Leiman solution, which moves the
values again.

The three omegas add up only on the whole-scale row under
`variance = "sums_load"`, which is a hierarchical variance partition:
there the whole-scale omega subscale counts all group-factor variance,
and omega total is exactly omega hierarchical plus omega subscale.
Everywhere else omega subscale counts only the assigned subscale
composites (whole-scale row) or the factor's own contribution (subscale
rows), so omega total need not equal omega hierarchical plus omega
subscale.

A solution with a single factor is scored as such on every input route,
and returns omega total, alpha, and the H index. Alpha assumes
essentially tau-equivalent items, which is nested in a one-factor model,
so the unidimensional case is the one for which reporting it is
defensible rather than merely possible. The other coefficients are
omitted because a single factor does not define them: omega subscale is
the variance due to the group factors, of which there are none; omega
hierarchical is the same quantity as omega total, the one factor
accounting for all of the common variance; and the ECV and the PUC are 1
by construction, which describes the number of factors in the model
rather than the evidence for unidimensionality it would be read as.

Which coefficient answers which question: omega hierarchical for whether
a total score can be read as a measure of one construct, omega subscale
for whether a subscale score adds anything beyond the general factor,
the H index for whether a factor is well defined by its indicators, and
ECV with PUC for whether a unidimensional model is defensible. Alpha
assumes essentially tau-equivalent items, whereas factor analysis yields
congeneric solutions, for which it is a lower bound; for a
multidimensional scale it is rarely the coefficient to report (Gignac,
2014). Composite reliability and the average variance extracted are not
returned, being coefficients of a fitted confirmatory measurement model;
[`semTools::compRelSEM()`](https://rdrr.io/pkg/semTools/man/compRelSEM.html)
and [`semTools::AVE()`](https://rdrr.io/pkg/semTools/man/AVE.html)
compute them from a `lavaan` fit.

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
and PUC). It still reports a whole-scale row, which describes the
composite of every variable: that row is labelled `"total"`, not `"g"`,
there being no general factor for it to be. A bare loading matrix is
read as a raw bifactor solution (general factor in the first column),
whose whole-scale row is the general factor and keeps the label `"g"`. A
one-factor
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
solution is scored from its unrotated loadings, having nothing to
rotate, and needs no oblique rotation as a solution with more factors
does; a one-column loading matrix (with or without the loading class an
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
solution's loadings carry), a single-factor `lavaan` fit, and components
carrying one factor are all read as the same single-factor solution and
return the same coefficients for it, whether that factor is given as the
general one or as the only group factor. Its row is labelled with the
name the input gives the factor – `fac_names` where that is supplied,
the solution's own factor label otherwise – and `"F1"` where the input
names it none. It is not labelled `"g"` as the general-factor row of a
solution with group factors is: a single factor is the whole model
rather than the general factor of a bifactor or hierarchical one. The
loading table of an
[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
solution is read as such a bifactor matrix too, without its `h2` and
`u2` display columns and with its uniquenesses; unlike the solution
itself it carries no correlation matrix, so supply `cormat` to score it
against the observed correlations rather than the model-implied ones.
The loading matrix of an ordinary factor solution is not a bifactor
matrix and is rejected – unless it is supplied with `Phi`, which makes
the pair a correlated-factors solution: it is then scored as the
solution it came from is, with the uniquenesses derived from the
loadings under `Phi` unless `u2` is given. Unlike the solution itself
the pair carries no correlation matrix, and none is reconstructed for
it, so `variance = "correlation"` needs one in `cormat` and errors
without it; `variance = "sums_load"` needs none. The matrix to pass is
the pattern, which is what `Phi` accompanies; a structure matrix would
be read as a pattern. `Phi` decides that reading on its own, whether or
not the matrix still carries the loading class
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
gives it: the class does not survive subsetting a matrix or reordering
its rows, while `Phi` stays, and the group factors of a hierarchy are
uncorrelated and have no factor correlations to supply. A matrix given
no `Phi` is a bifactor one, its first column the general factor, unless
it is an ordinary factor solution's loading matrix, which is rejected as
above. Because the two readings report different coefficients, a matrix
that carries no loading class and is read as a pattern warns that it
was. A Schmid-Leiman table states that it is a hierarchy, so `Phi` with
one is an error and not an argument that is dropped. The rejection
applies to a matrix of two or more factors: one of a single factor is
read as that factor, there being no hierarchy to mistake it for and no
factor intercorrelations to supply. A correlation matrix is rejected,
which belongs in `cormat`. Both `variance` conventions apply to a
correlated-factors solution: its factor intercorrelations enter the
composite's common variance either way. The indicator-to-factor
correspondences come from `factor_map` when it is supplied; otherwise
each variable is assigned to the group factor on which it loads most
strongly. For a bare bifactor loading matrix they default instead to the
non-zero group-loading pattern, which assumes the matrix carries
structural zeros – supply `factor_map` for an estimated bifactor matrix,
which has none. A Schmid-Leiman loading table is estimated and has no
such zeros, so it is mapped by strongest loading as the solution itself
is. For `lavaan` input the composite variances are model-implied
(`variance` is not used), and the coefficients are computed per group.
Which model a fit is follows from its structure rather than from
`g_name`: one that routes the covariances of its factors through a
higher-order factor is a second-order model, one in which some variable
loads on two or more factors is a bifactor model, and one in which every
variable loads on a single factor is a correlated-factors model, which
has no general factor and so, like an oblique
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
solution, omits omega hierarchical, ECV, and PUC while its factor
intercorrelations enter the remaining coefficients. The general-factor
coefficients split a composite's variance into a general part and one
part per group factor, which needs the latent variables to be
uncorrelated, so a bifactor model has to be fitted with
`orthogonal = TRUE` (not `lavaan`'s default) and a second-order model
must leave the covariances between its first-order factors at zero; such
a fit whose factors correlate is rejected rather than scored as though
they did not.

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

## From an oblique pattern matrix and its factor intercorrelations, which is the
## same correlated-factors solution and gives the same coefficients.
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
