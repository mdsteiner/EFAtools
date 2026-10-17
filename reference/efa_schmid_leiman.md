# Schmid-Leiman transformation

Transforms an oblique factor solution into a hierarchical solution with
one general factor and group factors that are uncorrelated with it and
with each other (Schmid & Leiman, 1957). With the same settings, the
results agree with
[`psych::schmid()`](https://rdrr.io/pkg/psych/man/schmid.html) and with
the SPSS syntax of Wolff and Preising (2005).

## Usage

``` r
efa_schmid_leiman(
  x,
  Phi = NULL,
  estimator = c("PAF", "ML", "ULS", "MINRES"),
  g_name = NULL,
  estimate_control = NULL
)
```

## Source

Schmid, J. & Leiman, J. M. (1957). The development of hierarchical
factor solutions. Psychometrika, 22(1), 53–61. doi:10.1007/BF02289209

Wolff, H.-G., & Preising, K. (2005). Exploring item and higher order
factor structure with the Schmid-Leiman solution: Syntax codes for SPSS
and SAS. Behavior Research Methods, 37 , 48–58. doi:10.3758/BF03206397

Mansolf, M., & Reise, S. P. (2016). Exploratory bifactor analysis: The
Schmid-Leiman orthogonalization and Jennrich-Bentler analytic rotations.
Multivariate Behavioral Research, 51(5), 698–717.
doi:10.1080/00273171.2016.1215898

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

- estimator:

  character. The estimator of the second-order
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md):
  `"PAF"` (default), `"ML"`, or `"ULS"`. `"MINRES"` is a synonym for
  `"ULS"`. Not used for a `lavaan` solution and for two first-order
  factors, for which no second-order EFA is done.

- g_name:

  character. The name of the general factor in a `lavaan` solution. If
  `NULL` (default), the latent variable that loads on the other latent
  variables is used.

- estimate_control:

  An
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  object with the settings of the second-order
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md).
  `NULL` (default) uses the
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  defaults. Not used for a `lavaan` solution and for two first-order
  factors.

## Value

A list of class `c("efa_schmid_leiman", "SL")`:

- orig_R:

  The correlation matrix of the input solution. `NA` for a pattern
  matrix.

- sl:

  The Schmid-Leiman loadings: the general factor (`g`), the group
  factors, the communalities (`h2`), and the uniquenesses (`u2`).

- L2:

  The second-order loadings, with one row for each first-order factor.

- vars_accounted:

  The variances explained by the general factor and the group factors.

- iter:

  The number of iterations of the second-order EFA.

- convergence:

  The convergence code of the second-order EFA (0 = converged; see
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)).

- settings:

  The estimator and the estimation settings of the second-order EFA.

`iter`, `convergence`, and `settings` are `NA` if no second-order EFA
was done, that is, for a `lavaan` solution and for two first-order
factors.

## Details

The factor intercorrelations are factor analyzed to get one second-order
(general) factor. The general-factor loadings of the items are their
first-order loadings times the second-order loadings. The group-factor
loadings are the first-order loadings times the residual standard
deviations of the first-order factors. The general factor and the group
factors thus divide the common variance, for example for McDonald's
omegas (see
[`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md)).

The transformation has proportionality constraints: in each group
factor, the ratio of the group-factor loading to the general-factor
loading is the same for all items. Thus, it recovers a bifactor
structure exactly only if each item loads on one group factor and these
ratios are equal. If items have cross-loadings, it overestimates their
general-factor loadings, and the result depends on the oblique rotation
of the first-order solution (Mansolf & Reise, 2016). The orthogonal
bifactor rotation of
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
(`rotation = "bifactorT"`) does not have these constraints.

With two first-order factors, the second-order loadings are not
identified. The two loadings then get the same absolute value, and their
product is the factor intercorrelation, as in
[`psych::schmid()`](https://rdrr.io/pkg/psych/man/schmid.html).

A `lavaan` solution must have one second-order factor and no other paths
or covariances between the latent variables. Its general-factor loadings
also include the direct loadings of the general factor on the items, if
the model has them.

The group factors keep the column order and the names of the input
solution. If the columns do not have unique names (or a name is `"g"`,
`"h2"`, or `"u2"`), all group factors are labelled `"F1"` to `"Fk"`.

## See also

Other factor rotation:
[`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md)

Other reliability coefficients:
[`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md),
[`print.efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_reliability.md)

## Examples

``` r
efa_mod <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
                   estimator = "PAF", rotation = "promax")
efa_schmid_leiman(efa_mod)
#> 
#> EFA for second-order loadings performed with estimator = 'PAF'
#> 
#> ── Schmid-Leiman Solution ──────────────────────────────────────────────────────
#> 
#>        g     F1     F2     F3    h2    u2
#> V1   .489  -.029   .022   .356  .367  .633
#> V2   .444  -.001   .042   .280  .277  .723
#> V3   .459   .036   .035   .263  .283  .717
#> V4   .522   .061  -.005   .320  .378  .622
#> V5   .468   .095  -.011   .254  .293  .707
#> V6   .478  -.044  -.031   .409  .399  .601
#> V7   .491   .001   .336   .054  .357  .643
#> V8   .463  -.010   .366   .018  .349  .651
#> V9   .457   .023   .347   .000  .330  .670
#> V10  .449  -.013   .425  -.041  .383  .617
#> V11  .477   .009   .224   .135  .297  .703
#> V12  .513   .012   .410  -.006  .432  .568
#> V13  .502   .372   .054  -.039  .395  .605
#> V14  .455   .332  -.043   .051  .322  .678
#> V15  .489   .340   .081  -.041  .363  .637
#> V16  .476   .336  -.032   .053  .343  .657
#> V17  .477   .402  -.023  -.016  .390  .610
#> V18  .485   .336   .003   .029  .350  .650
#> 
#> ── Second-order Loadings ───────────────────────────────────────────────────────
#> 
#>     F1    F2    F3
#> g  .796  .777  .814
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      g     F1    F2     F3
#> SS loadings        4.111  .770  .783   .642
#> Prop Tot Var        .228  .043  .044   .036
#> Cum Prop Tot Var    .228  .271  .315   .350
#> Prop Comm Var       .652  .122  .124   .102
#> Cum Prop Comm Var   .652  .774  .898  1.000

# A pattern matrix with its factor intercorrelations, for example from another
# program, gives the same solution
efa_schmid_leiman(efa_mod$rot_loadings, Phi = efa_mod$Phi)
#> 
#> EFA for second-order loadings performed with estimator = 'PAF'
#> 
#> ── Schmid-Leiman Solution ──────────────────────────────────────────────────────
#> 
#>        g     F1     F2     F3    h2    u2
#> V1   .489  -.029   .022   .356  .367  .633
#> V2   .444  -.001   .042   .280  .277  .723
#> V3   .459   .036   .035   .263  .283  .717
#> V4   .522   .061  -.005   .320  .378  .622
#> V5   .468   .095  -.011   .254  .293  .707
#> V6   .478  -.044  -.031   .409  .399  .601
#> V7   .491   .001   .336   .054  .357  .643
#> V8   .463  -.010   .366   .018  .349  .651
#> V9   .457   .023   .347   .000  .330  .670
#> V10  .449  -.013   .425  -.041  .383  .617
#> V11  .477   .009   .224   .135  .297  .703
#> V12  .513   .012   .410  -.006  .432  .568
#> V13  .502   .372   .054  -.039  .395  .605
#> V14  .455   .332  -.043   .051  .322  .678
#> V15  .489   .340   .081  -.041  .363  .637
#> V16  .476   .336  -.032   .053  .343  .657
#> V17  .477   .402  -.023  -.016  .390  .610
#> V18  .485   .336   .003   .029  .350  .650
#> 
#> ── Second-order Loadings ───────────────────────────────────────────────────────
#> 
#>     F1    F2    F3
#> g  .796  .777  .814
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      g     F1    F2     F3
#> SS loadings        4.111  .770  .783   .642
#> Prop Tot Var        .228  .043  .044   .036
#> Cum Prop Tot Var    .228  .271  .315   .350
#> Prop Comm Var       .652  .122  .124   .102
#> Cum Prop Comm Var   .652  .774  .898  1.000

# \donttest{
# A psych::fa() solution, with the psych settings for the second-order EFA
fa_mod <- psych::fa(test_models$baseline$cormat, nfactors = 3, n.obs = 500,
                    fm = "pa", rotate = "Promax")
efa_schmid_leiman(fa_mod, estimate_control = estimate_control(type = "psych"))
#> 
#> EFA for second-order loadings performed with estimator = 'PAF'
#> 
#> ── Schmid-Leiman Solution ──────────────────────────────────────────────────────
#> 
#>        g    PA2    PA3    PA1    h2    u2
#> V1   .501  -.037   .015   .338  .367  .633
#> V2   .453  -.006   .036   .265  .277  .723
#> V3   .468   .031   .029   .248  .283  .717
#> V4   .533   .054  -.012   .302  .378  .622
#> V5   .476   .090  -.017   .240  .293  .707
#> V6   .493  -.053  -.040   .390  .399  .601
#> V7   .492   .003   .335   .044  .357  .643
#> V8   .463  -.007   .366   .009  .349  .651
#> V9   .457   .026   .348  -.008  .330  .670
#> V10  .446  -.009   .427  -.048  .383  .617
#> V11  .482   .008   .222   .123  .297  .703
#> V12  .512   .016   .411  -.014  .432  .568
#> V13  .499   .375   .055  -.046  .395  .605
#> V14  .456   .332  -.044   .042  .322  .678
#> V15  .486   .343   .082  -.048  .363  .637
#> V16  .477   .336  -.033   .044  .343  .657
#> V17  .475   .404  -.023  -.023  .390  .610
#> V18  .485   .337   .002   .020  .350  .650
#> 
#> ── Second-order Loadings ───────────────────────────────────────────────────────
#> 
#>     PA2   PA3   PA1
#> g  .801  .780  .844
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      g     PA2   PA3   PA1
#> SS loadings        4.170  .775  .786   .575
#> Prop Tot Var        .232  .043  .044   .032
#> Cum Prop Tot Var    .232  .275  .318   .350
#> Prop Comm Var       .661  .123  .125   .091
#> Cum Prop Comm Var   .661  .784  .909  1.000

# A lavaan second-order CFA solution
if (requireNamespace("lavaan", quietly = TRUE)) {
  mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
          F2 =~ V7 + V8 + V9 + V10 + V11 + V12
          F3 =~ V13 + V14 + V15 + V16 + V17 + V18
          g =~ F1 + F2 + F3'
  fit <- lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
                     sample.nobs = 500)
  efa_schmid_leiman(fit)
}
#> 
#> ── Schmid-Leiman Solution ──────────────────────────────────────────────────────
#> 
#>        g    F1    F2    F3    h2    u2
#> V1   .515  .307  .000  .000  .360  .640
#> V2   .453  .270  .000  .000  .278  .722
#> V3   .459  .274  .000  .000  .286  .714
#> V4   .536  .320  .000  .000  .390  .610
#> V5   .478  .285  .000  .000  .310  .690
#> V6   .500  .299  .000  .000  .340  .660
#> V7   .491  .000  .352  .000  .366  .634
#> V8   .476  .000  .342  .000  .344  .656
#> V9   .468  .000  .336  .000  .332  .668
#> V10  .476  .000  .341  .000  .343  .657
#> V11  .448  .000  .321  .000  .304  .696
#> V12  .524  .000  .375  .000  .415  .585
#> V13  .505  .000  .000  .365  .388  .612
#> V14  .454  .000  .000  .328  .314  .686
#> V15  .485  .000  .000  .350  .358  .642
#> V16  .479  .000  .000  .346  .349  .651
#> V17  .495  .000  .000  .357  .373  .627
#> V18  .481  .000  .000  .347  .352  .648
#> 
#> ── Second-order Loadings ───────────────────────────────────────────────────────
#> 
#>     F1    F2    F3
#> g  .859  .813  .811
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      g     F1    F2     F3
#> SS loadings        4.239  .515  .714   .731
#> Prop Tot Var        .236  .029  .040   .041
#> Cum Prop Tot Var    .236  .264  .304   .344
#> Prop Comm Var       .684  .083  .115   .118
#> Cum Prop Comm Var   .684  .767  .882  1.000
# }
```
