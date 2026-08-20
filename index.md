# EFAtools

The EFAtools package performs exploratory factor analysis (EFA) and
compares EFA solutions. It offers current factor retention methods and
many estimation, rotation, and correlation options. The package
implements core iterative procedures, including principal axis factoring
(PAF), rotation, and polychoric correlation estimation, in C++ to
increase speed.

## Installation

You can install the release version from
[CRAN](https://cran.r-project.org/) with:

``` r

install.packages("EFAtools")
```

You can install the development version from
[GitHub](https://github.com/) with:

``` r

# install.packages("pak")
pak::pak("mdsteiner/EFAtools")
```

## Package Overview

The `efa_*` functions cover the steps of an EFA workflow:

- **Screening and suitability**:
  [`efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/efa_screen.md)
  checks the data for multivariate normality, outliers, and suitability
  for factor analysis.
  [`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md)
  and
  [`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md)
  run the Kaiser-Meyer-Olkin criterion and Bartlett’s test of sphericity
  separately.
- **Factor retention**:
  [`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
  runs several factor retention criteria with a single call. They are
  also available on their own:
  [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md),
  [`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md),
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md),
  [`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
  [`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md),
  [`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md),
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
  [`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
  and
  [`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md).
- **Fitting**:
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  fits the factor model.
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  and
  [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  configure the estimation and rotation settings.
- **Rotation and transformation**:
  [`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md)
  rotates a solution onto a target.
  [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  transforms an oblique solution into a hierarchical one.
- **Reliability**:
  [`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md)
  computes reliability and common-variance coefficients for a factor
  solution.
- **Factor scores**:
  [`efa_scores()`](https://mdsteiner.github.io/EFAtools/reference/efa_scores.md)
  estimates factor scores together with score-quality diagnostics.
- **Comparison**:
  [`efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/efa_compare.md)
  compares two solutions (loadings or communalities).
- **Averaging**:
  [`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md)
  averages solutions across implementations and methods to assess their
  stability.
- **Multiple groups**:
  [`efa_group()`](https://mdsteiner.github.io/EFAtools/reference/efa_group.md)
  fits a solution per group and compares them.
- **Multiple imputation**:
  [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)
  fits and pools solutions across multiply imputed data sets.
- **Simulation**:
  [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md)
  simulates data from a common-factor population model.
- **Power**:
  [`efa_power()`](https://mdsteiner.github.io/EFAtools/reference/efa_power.md)
  performs analytic and simulation-based power analysis.

The uppercase names
([`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md),
[`N_FACTORS()`](https://mdsteiner.github.io/EFAtools/reference/N_FACTORS.md),
…) still work and keep their arguments. Use the `efa_*` names for new
code.

The following vignettes and articles cover these in detail:

- [EFAtools](https://mdsteiner.github.io/EFAtools/articles/EFAtools.html):
  a complete EFA workflow, from screening to reliability and factor
  scores.
- [Migrating to the efa\_\*
  interface](https://mdsteiner.github.io/EFAtools/articles/Migrating_to_efa.html):
  what the old names and arguments map to.
- [EFA with ordinal and missing
  data](https://mdsteiner.github.io/EFAtools/articles/Ordinal_and_missing_data.html):
  polychoric and tetrachoric correlations, DWLS, FIML, and multiply
  imputed data.
- [Multigroup exploratory factor
  analysis](https://mdsteiner.github.io/EFAtools/articles/Multigroup_EFA.html):
  comparing solutions across groups and assessing approximate
  invariance.
- [Simulating data and power analysis for
  EFA](https://mdsteiner.github.io/EFAtools/articles/Simulation_and_power.html):
  simulating data and planning sample sizes.

## Examples

The following examples show some EFAtools functions.

``` r

# load the package
library(EFAtools)
```

### Factor Retention

Use
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
to test if your data are suitable for factor analysis. It runs multiple
factor retention criteria in a single call.

With raw data,
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
runs every criterion:

``` r


# Run multiple factor retention methods
efa_retain(GRiPS_raw)
#> Warning: The suggested maximum number of factors was 2, but the Hull method needs at
#> least 3.
#> ℹ Setting it to 3.
#> ── Tests for the suitability of the data for factor analysis ───────────────────
#> 
#> ✔ The Bartlett's test of sphericity was significant at an alpha level of .05:
#>   χ²(28) = 5054.06, p < .001. These data are probably suitable for factor
#>   analysis.
#> ✔ The Kaiser-Meyer-Olkin criterion is marvellous (KMO = 0.955). These data are
#>   probably suitable for factor analysis.
#> 
#> ── Suggested number of factors ─────────────────────────────────────────────────
#> 
#> 9 suggestions from 6 criteria, all suggesting 1 factor.
#> 
#> Comparison data
#> • Suggested number of factors: 1
#> 
#> Empirical Kaiser Criterion
#> • Braeken & van Assen (2017): 1
#> 
#> Hull method
#> Estimator: ML
#> • CAF: 1
#> • CFI: 1
#> • RMSEA: 1
#> 
#> Minimum average partial
#> • Original implementation (TR2): 1
#> • Revised implementation (TR4): 1
#> 
#> Next Eigenvalue Sufficiency Test
#> • Suggested number of factors: 1
#> 
#> Parallel analysis
#> Eigenvalues found using SMC; 1000 simulated datasets.
#> • SMC: 1
```

With a correlation matrix,
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
skips the criteria that need raw data:

``` r


efa_retain(DOSPERT$cormat, N = DOSPERT$N)
#> Warning: `x` is a correlation matrix, but "CD" needs raw data.
#> ℹ Skipping "CD".
#> ── Tests for the suitability of the data for factor analysis ───────────────────
#> 
#> ✔ The Bartlett's test of sphericity was significant at an alpha level of .05:
#>   χ²(780) = 16071.13, p < .001. These data are probably suitable for factor
#>   analysis.
#> ✔ The Kaiser-Meyer-Olkin criterion is meritorious (KMO = 0.9). These data are
#>   probably suitable for factor analysis.
#> 
#> ── Suggested number of factors ─────────────────────────────────────────────────
#> 
#> 8 suggestions from 5 criteria, ranging from 1 to 12 factors (most common: 10).
#> 
#> Empirical Kaiser Criterion
#> • Braeken & van Assen (2017): 10
#> 
#> Hull method
#> Estimator: ML
#> • CAF: 12
#> • CFI: 1
#> • RMSEA: 1
#> 
#> Minimum average partial
#> • Original implementation (TR2): 5
#> • Revised implementation (TR4): 6
#> 
#> Next Eigenvalue Sufficiency Test
#> • Suggested number of factors: 10
#> 
#> Parallel analysis
#> Eigenvalues found using SMC; 1000 simulated datasets.
#> • SMC: 12
#> 
#> ── Criteria that could not be run ──────────────────────────────────────────────
#> 
#> ! CD: needs raw data, but a correlation matrix was supplied
```

### EFA

#### Raw-Data Input

With raw data, you can use every feature: sandwich and bootstrap
standard errors, DWLS estimation with polychoric correlations, and
two-stage FIML estimation of correlations.

In the first example below, which uses bootstrap SEs (`se = "np-boot"`),
the confidence intervals are percentile intervals over refitted
resamples. For the loadings and factor correlations, the intervals are
centred on the point estimate. For the indices derived from the
chi-square (RMSEA, AIC, BIC, ECVI), the intervals are located above the
point estimate, because each resample carries the sample’s own misfit
plus fresh sampling noise. If a point estimate falls below its own lower
bound, this shift is the reason, not a miscomputed interval. CFI and TLI
are unaffected because they are ratios: their baseline chi-square shifts
together with the model chi-square.

``` r


# ULS / MINRES estimation with oblimin rotation and bootstrap SEs
mod <- efa_fit(DOSPERT_raw, n_factors = 5, estimator = "uls", rotation = "oblimin",
               se = "np-boot", seed = 1)
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
mod
#> 
#> EFA performed with estimator = 'ULS' and rotation = 'oblimin'.
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>           F1     F2     F3     F4     F5    h2    u2
#> ethR_1   .513  -.018   .030  -.016   .130  .309  .691
#> ethR_2   .518  -.044   .078   .019   .045  .304  .696
#> ethR_3   .639  -.001   .024  -.223   .081  .490  .510
#> ethR_4   .586  -.122  -.050  -.060   .046  .295  .705
#> ethR_5   .477   .065  -.010  -.127   .032  .267  .733
#> ethR_6   .621  -.098  -.007  -.016  -.021  .345  .655
#> finR_1  -.004  -.005   .841  -.021   .025  .717  .283
#> finR_2  -.066   .029  -.045   .068   .688  .476  .524
#> finR_3  -.005  -.013   .856   .010   .016  .730  .270
#> finR_4   .072   .041   .090  -.040   .710  .600  .400
#> finR_5  -.005  -.029   .873   .000   .040  .768  .232
#> finR_6   .054   .064   .093   .085   .683  .599  .401
#> heaR_1   .426   .087   .101   .087  -.036  .273  .727
#> heaR_2   .453   .053   .050   .136  -.050  .262  .738
#> heaR_3   .415   .130   .071   .019  -.052  .257  .743
#> heaR_4   .362   .163   .123  -.015  -.066  .254  .746
#> heaR_5   .382   .091  -.019   .123  -.057  .185  .815
#> heaR_6   .430   .206   .026   .138   .003  .338  .662
#> recR_1   .017   .407  -.035   .217   .026  .254  .746
#> recR_2   .117   .531   .111  -.101   .038  .410  .590
#> recR_3   .060   .619   .026   .003   .054  .452  .548
#> recR_4  -.072   .861  -.033  -.059   .027  .682  .318
#> recR_5  -.008   .805   .013  -.091  -.003  .628  .372
#> recR_6  -.020   .637   .031   .025   .102  .467  .533
#> socR_1  -.029  -.085  -.071   .646  -.004  .419  .581
#> socR_2   .093  -.027   .031   .679   .039  .474  .526
#> socR_3  -.133  -.058   .018   .640  -.005  .416  .584
#> socR_4  -.004   .018   .032   .614   .007  .383  .617
#> socR_5   .049   .103  -.045   .379   .051  .185  .815
#> socR_6   .004  -.008   .016   .549   .041  .308  .692
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1     F2     F3     F4     F5
#> F1  1.000
#> F2   .372  1.000
#> F3   .448   .319  1.000
#> F4   .006   .200  -.042  1.000
#> F5   .154   .290   .344   .145  1.000
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F2     F3     F4     F5
#> SS loadings        3.163  2.940  2.457  2.323  1.664
#> Prop Tot Var        .105   .098   .082   .077   .055
#> Cum Prop Tot Var    .105   .203   .285   .363   .418
#> Prop Comm Var       .252   .234   .196   .185   .133
#> Cum Prop Comm Var   .252   .486   .682   .867  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> χ²(295) = 3604.56, p < .001
#> CFI [95% bootstrap-CI]: .90 [.88, .90]
#> TLI [95% bootstrap-CI]: .85 [.82, .85]
#> RMSEA [90% CI] [95% bootstrap-CI]: .06 [.06; .06] [.06, .07]
#> AIC [95% bootstrap-CI]: 3014.56 [3034.98, 3700.48]
#> BIC [95% bootstrap-CI]: 1230.83 [1251.25, 1916.74]
#> ECVI [95% bootstrap-CI]: 1.26 [1.27, 1.48]
#> CAF [95% bootstrap-CI]: .43 [.43, .45]
#> SRMR [95% bootstrap-CI]: .03 [.03, .04]
#> 
#> Note: Bootstrap CIs based on 1000 bootstrap samples.
# detailed output with summary()
summary(mod)
#> 
#> EFA performed with estimator = 'ULS' and rotation = 'oblimin'.
#> 
#> ── Model Diagnostics ───────────────────────────────────────────────────────────
#> 
#> Factors: 5
#> Variables: 30
#> N: 3123
#> Bootstrap samples: 1000
#> Valid target-rotated samples: 1000 out of 1000
#> Rotation local optima: 1 distinct from 6 of 101 starts
#> Heywood cases: 0
#> Cross-loading items (|loading| >= .300): 0
#> Items without salient loading (|loading| >= .300): 0
#> Factors with fewer than 3 salient indicators: 0
#> Items with primary-loading gap < .200: 2
#> Largest |residual|: .246
#> Factor intercorrelations > .85: none
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>           F1     F2     F3     F4     F5    h2    u2
#> ethR_1   .513  -.018   .030  -.016   .130  .309  .691
#> ethR_2   .518  -.044   .078   .019   .045  .304  .696
#> ethR_3   .639  -.001   .024  -.223   .081  .490  .510
#> ethR_4   .586  -.122  -.050  -.060   .046  .295  .705
#> ethR_5   .477   .065  -.010  -.127   .032  .267  .733
#> ethR_6   .621  -.098  -.007  -.016  -.021  .345  .655
#> finR_1  -.004  -.005   .841  -.021   .025  .717  .283
#> finR_2  -.066   .029  -.045   .068   .688  .476  .524
#> finR_3  -.005  -.013   .856   .010   .016  .730  .270
#> finR_4   .072   .041   .090  -.040   .710  .600  .400
#> finR_5  -.005  -.029   .873   .000   .040  .768  .232
#> finR_6   .054   .064   .093   .085   .683  .599  .401
#> heaR_1   .426   .087   .101   .087  -.036  .273  .727
#> heaR_2   .453   .053   .050   .136  -.050  .262  .738
#> heaR_3   .415   .130   .071   .019  -.052  .257  .743
#> heaR_4   .362   .163   .123  -.015  -.066  .254  .746
#> heaR_5   .382   .091  -.019   .123  -.057  .185  .815
#> heaR_6   .430   .206   .026   .138   .003  .338  .662
#> recR_1   .017   .407  -.035   .217   .026  .254  .746
#> recR_2   .117   .531   .111  -.101   .038  .410  .590
#> recR_3   .060   .619   .026   .003   .054  .452  .548
#> recR_4  -.072   .861  -.033  -.059   .027  .682  .318
#> recR_5  -.008   .805   .013  -.091  -.003  .628  .372
#> recR_6  -.020   .637   .031   .025   .102  .467  .533
#> socR_1  -.029  -.085  -.071   .646  -.004  .419  .581
#> socR_2   .093  -.027   .031   .679   .039  .474  .526
#> socR_3  -.133  -.058   .018   .640  -.005  .416  .584
#> socR_4  -.004   .018   .032   .614   .007  .383  .617
#> socR_5   .049   .103  -.045   .379   .051  .185  .815
#> socR_6   .004  -.008   .016   .549   .041  .308  .692
#> 
#> ── 95% bootstrap CIs for salient rotated loadings ──────────────────────────────
#> 
#> Variable  Factor  est    lower  upper
#> ethR_1    F1       .513   .465   .553
#> ethR_2    F1       .518   .472   .558
#> ethR_3    F1       .639   .589   .677
#> ethR_4    F1       .586   .532   .630
#> ethR_5    F1       .477   .431   .523
#> ethR_6    F1       .621   .572   .662
#> heaR_1    F1       .426   .381   .472
#> heaR_2    F1       .453   .407   .497
#> heaR_3    F1       .415   .361   .472
#> heaR_4    F1       .362   .306   .419
#> heaR_5    F1       .382   .326   .435
#> heaR_6    F1       .430   .382   .475
#> recR_1    F2       .407   .366   .447
#> recR_2    F2       .531   .487   .572
#> recR_3    F2       .619   .578   .654
#> recR_4    F2       .861   .825   .891
#> recR_5    F2       .805   .771   .835
#> recR_6    F2       .637   .600   .669
#> finR_1    F3       .841   .801   .868
#> finR_3    F3       .856   .819   .881
#> finR_5    F3       .873   .837   .899
#> socR_1    F4       .646   .612   .680
#> socR_2    F4       .679   .648   .708
#> socR_3    F4       .640   .602   .674
#> socR_4    F4       .614   .581   .648
#> socR_5    F4       .379   .339   .418
#> socR_6    F4       .549   .511   .587
#> finR_2    F5       .688   .646   .723
#> finR_4    F5       .710   .667   .740
#> finR_6    F5       .683   .642   .715
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1     F2     F3     F4     F5
#> F1  1.000
#> F2   .372  1.000
#> F3   .448   .319  1.000
#> F4   .006   .200  -.042  1.000
#> F5   .154   .290   .344   .145  1.000
#> 
#> ── 95% bootstrap CIs for factor intercorrelations ──────────────────────────────
#> 
#> Factors   est    lower  upper
#> F1 ~~ F2   .372   .327   .405
#> F1 ~~ F3   .448   .400   .481
#> F1 ~~ F4   .006  -.037   .049
#> F1 ~~ F5   .154   .108   .196
#> F2 ~~ F3   .319   .270   .359
#> F2 ~~ F4   .200   .160   .236
#> F2 ~~ F5   .290   .241   .326
#> F3 ~~ F4  -.042  -.082   .002
#> F3 ~~ F5   .344   .293   .373
#> F4 ~~ F5   .145   .098   .185
#> 
#> ── Structure Matrix ────────────────────────────────────────────────────────────
#> 
#>           F1    F2     F3     F4    F5
#> ethR_1   .540  .217   .299   .001  .212
#> ethR_2   .543  .190   .310   .017  .142
#> ethR_3   .661  .223   .348  -.209  .155
#> ethR_4   .525  .081   .192  -.072  .075
#> ethR_5   .500  .223   .240  -.106  .102
#> ethR_6   .578  .121   .233  -.035  .042
#> finR_1   .374  .265   .846  -.054  .309
#> finR_2   .031  .203   .168   .175  .681
#> finR_3   .376  .265   .854  -.027  .307
#> finR_4   .237  .294   .382   .068  .759
#> finR_5   .382  .260   .875  -.037  .331
#> finR_6   .225  .329   .369   .193  .754
#> heaR_1   .499  .284   .303   .097  .102
#> heaR_2   .488  .250   .247   .140  .072
#> heaR_3   .487  .296   .279   .037  .077
#> heaR_4   .467  .315   .315   .005  .077
#> heaR_5   .399  .235   .156   .136  .039
#> heaR_6   .519  .402   .279   .181  .158
#> recR_1   .158  .453   .103   .304  .166
#> recR_2   .370  .601   .350   .007  .234
#> recR_3   .310  .666   .269   .134  .252
#> recR_4   .237  .820   .221   .118  .246
#> recR_5   .297  .787   .270   .069  .221
#> recR_6   .247  .674   .260   .165  .298
#> socR_1  -.089  .010  -.140   .631  .036
#> socR_2   .107  .164   .049   .678  .155
#> socR_3  -.143  .025  -.089   .626  .057
#> socR_4   .022  .152   .013   .617  .111
#> socR_5   .078  .197   .012   .409  .128
#> socR_6   .018  .120   .006   .553  .124
#> 
#> ── Simple Structure Diagnostics ────────────────────────────────────────────────
#> 
#> Items with primary-loading gap < .200:
#> • heaR_4: F1 = .362, F2 = .163
#> • recR_1: F2 = .407, F4 = .217
#> 
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F2     F3     F4     F5
#> SS loadings        3.163  2.940  2.457  2.323  1.664
#> Prop Tot Var        .105   .098   .082   .077   .055
#> Cum Prop Tot Var    .105   .203   .285   .363   .418
#> Prop Comm Var       .252   .234   .196   .185   .133
#> Cum Prop Comm Var   .252   .486   .682   .867  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> χ²(295) = 3604.56, p < .001
#> CFI [95% bootstrap-CI]: .90 [.88, .90]
#> TLI [95% bootstrap-CI]: .85 [.82, .85]
#> RMSEA [90% CI] [95% bootstrap-CI]: .06 [.06; .06] [.06, .07]
#> AIC [95% bootstrap-CI]: 3014.56 [3034.98, 3700.48]
#> BIC [95% bootstrap-CI]: 1230.83 [1251.25, 1916.74]
#> ECVI [95% bootstrap-CI]: 1.26 [1.27, 1.48]
#> CAF [95% bootstrap-CI]: .43 [.43, .45]
#> SRMR [95% bootstrap-CI]: .03 [.03, .04]
#> 
#> Note: Bootstrap CIs based on 1000 bootstrap samples.
#> 
#> ── Residual Diagnostics ────────────────────────────────────────────────────────
#> 
#> Residual cutoff: |r| > .100
#> Number of large residuals: 7
#> Largest absolute residual: .246
#> 
#> Largest residuals:
#> • heaR_3 ~~ heaR_4: .246
#> • socR_5 ~~ socR_6: .190
#> • socR_2 ~~ socR_4: .145
#> • recR_4 ~~ recR_5: .138
#> • heaR_1 ~~ heaR_2: .135
#> • recR_2 ~~ recR_3: .126
#> • recR_1 ~~ recR_3: .112
#> 
#> Inspect the residual matrix for details (e.g., with residuals()).

# inspect residuals with residuals()
residuals(mod)
#>               ethR_1       ethR_2       ethR_3       ethR_4        ethR_5
#> ethR_1  0.0000000000  0.001919361  0.040110633 -0.004366305  0.0398686298
#> ethR_2  0.0019193609  0.000000000  0.051351673 -0.020430720 -0.0049429995
#> ethR_3  0.0401106327  0.051351673  0.000000000  0.081672100  0.0288818382
#> ethR_4 -0.0043663054 -0.020430720  0.081672100  0.000000000  0.0333149594
#> ethR_5  0.0398686298 -0.004942999  0.028881838  0.033314959  0.0000000000
#> ethR_6  0.0405841568  0.026282970  0.045450944  0.028149579  0.0187898709
#> finR_1  0.0111270742 -0.006434044  0.014245152  0.009380118  0.0025759402
#> finR_2  0.0075531702  0.003704974 -0.015347877  0.017744475  0.0176028890
#> finR_3  0.0062334676  0.003308966  0.001315495  0.009312462  0.0034945485
#> finR_4 -0.0166237499 -0.014663821 -0.017357093 -0.016698897 -0.0174312455
#> finR_5 -0.0082532566  0.003043982 -0.004687753  0.004285391  0.0092379622
#> finR_6 -0.0092860185 -0.003711480 -0.006966906 -0.013540434 -0.0301857466
#> heaR_1  0.0049738778  0.008359162 -0.055880188  0.007130092  0.0010109939
#> heaR_2 -0.0513834281  0.055164201 -0.070564379 -0.036729296 -0.0331386742
#> heaR_3 -0.0097218604 -0.049226343 -0.048071095 -0.053699894 -0.0525279163
#> heaR_4 -0.0330385584 -0.008960833 -0.015555983 -0.039906944 -0.0221868272
#> heaR_5 -0.0326686399 -0.078645168 -0.046818454  0.042178391 -0.0501154768
#> heaR_6 -0.0324404536 -0.012615694 -0.084461504 -0.064143974 -0.0056839539
#> recR_1  0.0171481594 -0.020115444 -0.023081276 -0.018231529 -0.0215742190
#> recR_2 -0.0078417438 -0.019215500  0.026180414 -0.018336628  0.0005970064
#> recR_3  0.0213627795 -0.023521702  0.007778518  0.023307837 -0.0194798459
#> recR_4 -0.0014847694  0.019319774  0.023063666  0.031610969  0.0222082596
#> recR_5 -0.0059703080  0.014000766  0.020587444  0.031113565  0.0168282430
#> recR_6  0.0209481909  0.027891094  0.015364308 -0.029059408  0.0169439828
#> socR_1  0.0008069591 -0.012454685 -0.018898845  0.029253070  0.0035236268
#> socR_2  0.0117411813  0.004707962  0.012155553 -0.002291726  0.0282702203
#> socR_3  0.0229908383 -0.005414279  0.026567346  0.025314364 -0.0142252492
#> socR_4  0.0224634597  0.006150381  0.042416616 -0.019458383  0.0157349168
#> socR_5 -0.0150925009  0.014061906  0.002214296 -0.012377131  0.0357741533
#> socR_6 -0.0148244794  0.030629585  0.024032270  0.039389222 -0.0112568143
#>              ethR_6        finR_1        finR_2        finR_3        finR_4
#> ethR_1  0.040584157  1.112707e-02  0.0075531702  6.233468e-03 -0.0166237499
#> ethR_2  0.026282970 -6.434044e-03  0.0037049741  3.308966e-03 -0.0146638210
#> ethR_3  0.045450944  1.424515e-02 -0.0153478775  1.315495e-03 -0.0173570934
#> ethR_4  0.028149579  9.380118e-03  0.0177444747  9.312462e-03 -0.0166988973
#> ethR_5  0.018789871  2.575940e-03  0.0176028890  3.494548e-03 -0.0174312455
#> ethR_6  0.000000000 -5.661974e-03 -0.0162182518 -1.020901e-03 -0.0043628684
#> finR_1 -0.005661974  0.000000e+00  0.0087007058  7.563691e-05 -0.0056162714
#> finR_2 -0.016218252  8.700706e-03  0.0000000000  2.498749e-03  0.0008732958
#> finR_3 -0.001020901  7.563691e-05  0.0024987491  0.000000e+00 -0.0069267035
#> finR_4 -0.004362868 -5.616271e-03  0.0008732958 -6.926703e-03  0.0000000000
#> finR_5  0.025056131  3.479851e-03  0.0004840636  3.513885e-03  0.0041245381
#> finR_6 -0.016471269 -1.355585e-02 -0.0082562023  2.963387e-03  0.0214565406
#> heaR_1 -0.018437189  1.792566e-02  0.0036665710  8.973994e-04  0.0137074928
#> heaR_2 -0.016473594 -2.387758e-02  0.0155674333  6.578790e-03  0.0077786798
#> heaR_3 -0.038072331 -1.429475e-02 -0.0123398249 -7.788668e-03  0.0336369458
#> heaR_4 -0.060953169 -3.755198e-03 -0.0399505434 -6.037948e-03  0.0222392295
#> heaR_5 -0.065974441 -6.846035e-03 -0.0034818851 -8.308054e-03  0.0246862895
#> heaR_6 -0.003199384 -7.804389e-03  0.0103731217 -5.703732e-03  0.0123498811
#> recR_1 -0.037810156  6.737341e-03  0.0258258476  5.928859e-03 -0.0112637689
#> recR_2 -0.032038662  1.225974e-02  0.0006637449 -1.310818e-02  0.0020505490
#> recR_3 -0.016985690 -3.360856e-03  0.0081718050 -1.803233e-03 -0.0009077597
#> recR_4  0.041722071 -9.300863e-03 -0.0062793787  9.980525e-03  0.0062919397
#> recR_5  0.035964299  4.252238e-04 -0.0071202488  9.513055e-03  0.0011025022
#> recR_6  0.010032347  6.758971e-03  0.0046227367 -3.828779e-03 -0.0324267616
#> socR_1  0.016903192  2.591856e-03  0.0336326821  7.609337e-03  0.0021698992
#> socR_2 -0.001469239  1.023198e-03  0.0135072705  6.453139e-03 -0.0085527698
#> socR_3  0.012457793 -3.306369e-03 -0.0127770353  9.893623e-03 -0.0147406883
#> socR_4 -0.001953734 -7.945170e-03 -0.0248896127 -1.207071e-02  0.0125014063
#> socR_5  0.043324393  2.112224e-02  0.0038238407 -1.807355e-02 -0.0084603641
#> socR_6  0.011257856  7.598806e-03 -0.0343171349 -1.676254e-03 -0.0055717516
#>               finR_5       finR_6        heaR_1       heaR_2       heaR_3
#> ethR_1 -8.253257e-03 -0.009286018  0.0049738778 -0.051383428 -0.009721860
#> ethR_2  3.043982e-03 -0.003711480  0.0083591619  0.055164201 -0.049226343
#> ethR_3 -4.687753e-03 -0.006966906 -0.0558801885 -0.070564379 -0.048071095
#> ethR_4  4.285391e-03 -0.013540434  0.0071300921 -0.036729296 -0.053699894
#> ethR_5  9.237962e-03 -0.030185747  0.0010109939 -0.033138674 -0.052527916
#> ethR_6  2.505613e-02 -0.016471269 -0.0184371889 -0.016473594 -0.038072331
#> finR_1  3.479851e-03 -0.013555851  0.0179256603 -0.023877580 -0.014294749
#> finR_2  4.840636e-04 -0.008256202  0.0036665710  0.015567433 -0.012339825
#> finR_3  3.513885e-03  0.002963387  0.0008973994  0.006578790 -0.007788668
#> finR_4  4.124538e-03  0.021456541  0.0137074928  0.007778680  0.033636946
#> finR_5  1.110223e-16 -0.004350768 -0.0025407724  0.008688483 -0.019196106
#> finR_6 -4.350768e-03  0.000000000 -0.0086506865  0.010828731  0.032743193
#> heaR_1 -2.540772e-03 -0.008650687  0.0000000000  0.134739202 -0.037920792
#> heaR_2  8.688483e-03  0.010828731  0.1347392019  0.000000000  0.011491920
#> heaR_3 -1.919611e-02  0.032743193 -0.0379207917  0.011491920  0.000000000
#> heaR_4 -1.814601e-02  0.042818853 -0.0879576040 -0.039916181  0.246046503
#> heaR_5  2.483989e-03  0.017696479  0.0069513895  0.053173834  0.097484207
#> heaR_6 -4.840824e-03  0.019470073  0.0400382493  0.048967368  0.042175647
#> recR_1 -1.332485e-02 -0.010110670  0.0385138743  0.042180221 -0.028312126
#> recR_2 -9.949879e-03 -0.006382141  0.0042548181 -0.019268210  0.014912310
#> recR_3 -1.926819e-03 -0.019009888  0.0077585766 -0.004347300 -0.016209976
#> recR_4  1.754674e-02 -0.011312457  0.0020272674 -0.007770434 -0.032917724
#> recR_5  7.495640e-03 -0.001281005  0.0126986626 -0.010735364 -0.019892009
#> recR_6 -4.329788e-04  0.016609493 -0.0481766164 -0.020457652 -0.015608991
#> socR_1  7.722306e-03 -0.026656562  0.0290079981  0.018384301 -0.016464683
#> socR_2 -8.089515e-03 -0.008774178 -0.0175652234 -0.010820069 -0.007281281
#> socR_3  9.901727e-04  0.013649984 -0.0130453313 -0.033654039 -0.007074202
#> socR_4  7.128660e-03 -0.004181621 -0.0559967679 -0.022966696  0.002606542
#> socR_5 -9.945617e-04 -0.010900238  0.0065213961 -0.018231448 -0.021954127
#> socR_6  2.060514e-04  0.019249088 -0.0149993691 -0.047981685 -0.004446855
#>              heaR_4       heaR_5       heaR_6       recR_1        recR_2
#> ethR_1 -0.033038558 -0.032668640 -0.032440454  0.017148159 -0.0078417438
#> ethR_2 -0.008960833 -0.078645168 -0.012615694 -0.020115444 -0.0192155004
#> ethR_3 -0.015555983 -0.046818454 -0.084461504 -0.023081276  0.0261804145
#> ethR_4 -0.039906944  0.042178391 -0.064143974 -0.018231529 -0.0183366281
#> ethR_5 -0.022186827 -0.050115477 -0.005683954 -0.021574219  0.0005970064
#> ethR_6 -0.060953169 -0.065974441 -0.003199384 -0.037810156 -0.0320386620
#> finR_1 -0.003755198 -0.006846035 -0.007804389  0.006737341  0.0122597355
#> finR_2 -0.039950543 -0.003481885  0.010373122  0.025825848  0.0006637449
#> finR_3 -0.006037948 -0.008308054 -0.005703732  0.005928859 -0.0131081803
#> finR_4  0.022239229  0.024686290  0.012349881 -0.011263769  0.0020505490
#> finR_5 -0.018146006  0.002483989 -0.004840824 -0.013324848 -0.0099498787
#> finR_6  0.042818853  0.017696479  0.019470073 -0.010110670 -0.0063821411
#> heaR_1 -0.087957604  0.006951390  0.040038249  0.038513874  0.0042548181
#> heaR_2 -0.039916181  0.053173834  0.048967368  0.042180221 -0.0192682103
#> heaR_3  0.246046503  0.097484207  0.042175647 -0.028312126  0.0149123101
#> heaR_4  0.000000000  0.055888592  0.052126969 -0.017062979  0.0235161399
#> heaR_5  0.055888592  0.000000000  0.088708963  0.064717712  0.0069537566
#> heaR_6  0.052126969  0.088708963  0.000000000  0.027032976  0.0055625248
#> recR_1 -0.017062979  0.064717712  0.027032976  0.000000000  0.0421362847
#> recR_2  0.023516140  0.006953757  0.005562525  0.042136285  0.0000000000
#> recR_3 -0.007881816 -0.013440453 -0.004225714  0.111558311  0.1261720989
#> recR_4 -0.049566716 -0.032688075 -0.039437873 -0.055620434 -0.0780331049
#> recR_5 -0.034114831 -0.030578314 -0.035273311 -0.068600428 -0.0419696778
#> recR_6  0.025719051 -0.016720621  0.011850316  0.001493897 -0.0186904376
#> socR_1 -0.033282944  0.022436783 -0.021808583  0.044539158 -0.0412123383
#> socR_2 -0.001017320 -0.031735813  0.008396804 -0.039516375 -0.0069890310
#> socR_3 -0.010387563 -0.003902217 -0.020379358  0.006295940  0.0025606841
#> socR_4  0.021094781 -0.028690728 -0.036873922 -0.047504851  0.0200858684
#> socR_5 -0.021543449 -0.046677307  0.015586106 -0.010230758  0.0044974004
#> socR_6  0.014141970 -0.015473665 -0.035944127 -0.016774471  0.0080665231
#>               recR_3       recR_4        recR_5        recR_6        socR_1
#> ethR_1  0.0213627795 -0.001484769 -0.0059703080  0.0209481909  0.0008069591
#> ethR_2 -0.0235217023  0.019319774  0.0140007659  0.0278910937 -0.0124546847
#> ethR_3  0.0077785183  0.023063666  0.0205874440  0.0153643084 -0.0188988448
#> ethR_4  0.0233078374  0.031610969  0.0311135645 -0.0290594075  0.0292530699
#> ethR_5 -0.0194798459  0.022208260  0.0168282430  0.0169439828  0.0035236268
#> ethR_6 -0.0169856905  0.041722071  0.0359642985  0.0100323466  0.0169031921
#> finR_1 -0.0033608558 -0.009300863  0.0004252238  0.0067589714  0.0025918557
#> finR_2  0.0081718050 -0.006279379 -0.0071202488  0.0046227367  0.0336326821
#> finR_3 -0.0018032332  0.009980525  0.0095130549 -0.0038287794  0.0076093366
#> finR_4 -0.0009077597  0.006291940  0.0011025022 -0.0324267616  0.0021698992
#> finR_5 -0.0019268186  0.017546741  0.0074956402 -0.0004329788  0.0077223064
#> finR_6 -0.0190098880 -0.011312457 -0.0012810050  0.0166094933 -0.0266565623
#> heaR_1  0.0077585766  0.002027267  0.0126986626 -0.0481766164  0.0290079981
#> heaR_2 -0.0043473004 -0.007770434 -0.0107353643 -0.0204576519  0.0183843007
#> heaR_3 -0.0162099762 -0.032917724 -0.0198920088 -0.0156089909 -0.0164646833
#> heaR_4 -0.0078818158 -0.049566716 -0.0341148310  0.0257190507 -0.0332829443
#> heaR_5 -0.0134404535 -0.032688075 -0.0305783145 -0.0167206211  0.0224367832
#> heaR_6 -0.0042257139 -0.039437873 -0.0352733110  0.0118503157 -0.0218085835
#> recR_1  0.1115583110 -0.055620434 -0.0686004282  0.0014938968  0.0445391583
#> recR_2  0.1261720989 -0.078033105 -0.0419696778 -0.0186904376 -0.0412123383
#> recR_3  0.0000000000 -0.053597834 -0.0558336643 -0.0213076569 -0.0101104824
#> recR_4 -0.0535978342  0.000000000  0.1378809240  0.0330595053  0.0295162153
#> recR_5 -0.0558336643  0.137880924  0.0000000000 -0.0109523874  0.0232207682
#> recR_6 -0.0213076569  0.033059505 -0.0109523874  0.0000000000 -0.0269745409
#> socR_1 -0.0101104824  0.029516215  0.0232207682 -0.0269745409  0.0000000000
#> socR_2 -0.0123826202  0.012570217  0.0140286215  0.0146283687  0.0209785075
#> socR_3  0.0060427671  0.007252872  0.0085777710 -0.0056397962  0.0168906777
#> socR_4  0.0079668791  0.009953680  0.0134165965  0.0030089091 -0.0072011748
#> socR_5 -0.0098538521  0.005302575  0.0053179702  0.0096066819 -0.0351485978
#> socR_6 -0.0087254811 -0.004799714 -0.0051350217  0.0245084501 -0.0419866266
#>              socR_2        socR_3       socR_4        socR_5        socR_6
#> ethR_1  0.011741181  0.0229908383  0.022463460 -0.0150925009 -0.0148244794
#> ethR_2  0.004707962 -0.0054142788  0.006150381  0.0140619063  0.0306295848
#> ethR_3  0.012155553  0.0265673465  0.042416616  0.0022142959  0.0240322702
#> ethR_4 -0.002291726  0.0253143638 -0.019458383 -0.0123771310  0.0393892216
#> ethR_5  0.028270220 -0.0142252492  0.015734917  0.0357741533 -0.0112568143
#> ethR_6 -0.001469239  0.0124577931 -0.001953734  0.0433243930  0.0112578556
#> finR_1  0.001023198 -0.0033063690 -0.007945170  0.0211222411  0.0075988064
#> finR_2  0.013507270 -0.0127770353 -0.024889613  0.0038238407 -0.0343171349
#> finR_3  0.006453139  0.0098936231 -0.012070712 -0.0180735515 -0.0016762542
#> finR_4 -0.008552770 -0.0147406883  0.012501406 -0.0084603641 -0.0055717516
#> finR_5 -0.008089515  0.0009901727  0.007128660 -0.0009945617  0.0002060514
#> finR_6 -0.008774178  0.0136499841 -0.004181621 -0.0109002377  0.0192490883
#> heaR_1 -0.017565223 -0.0130453313 -0.055996768  0.0065213961 -0.0149993691
#> heaR_2 -0.010820069 -0.0336540395 -0.022966696 -0.0182314479 -0.0479816855
#> heaR_3 -0.007281281 -0.0070742019  0.002606542 -0.0219541268 -0.0044468546
#> heaR_4 -0.001017320 -0.0103875632  0.021094781 -0.0215434489  0.0141419698
#> heaR_5 -0.031735813 -0.0039022170 -0.028690728 -0.0466773074 -0.0154736648
#> heaR_6  0.008396804 -0.0203793580 -0.036873922  0.0155861065 -0.0359441270
#> recR_1 -0.039516375  0.0062959397 -0.047504851 -0.0102307579 -0.0167744712
#> recR_2 -0.006989031  0.0025606841  0.020085868  0.0044974004  0.0080665231
#> recR_3 -0.012382620  0.0060427671  0.007966879 -0.0098538521 -0.0087254811
#> recR_4  0.012570217  0.0072528718  0.009953680  0.0053025745 -0.0047997140
#> recR_5  0.014028622  0.0085777710  0.013416597  0.0053179702 -0.0051350217
#> recR_6  0.014628369 -0.0056397962  0.003008909  0.0096066819  0.0245084501
#> socR_1  0.020978507  0.0168906777 -0.007201175 -0.0351485978 -0.0419866266
#> socR_2  0.000000000 -0.0473998150  0.145141225 -0.0500361053 -0.0588975674
#> socR_3 -0.047399815  0.0000000000 -0.007742205 -0.0054723176  0.0775395800
#> socR_4  0.145141225 -0.0077422045  0.000000000 -0.0483058266 -0.0487045862
#> socR_5 -0.050036105 -0.0054723176 -0.048305827  0.0000000000  0.1897362412
#> socR_6 -0.058897567  0.0775395800 -0.048704586  0.1897362412  0.0000000000

# DWLS estimation based on polychoric correlations, with robust sandwich SEs
mod <- efa_fit(GRiPS_raw, n_factors = 1, estimator = "dwls", cor_method = "poly",
               se = "sandwich")
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
#> Warning: 22 variable pairs have an empty response-category combination despite a
#> non-negligible expected count.
#> ✖ Affected pairs: "fun-friends", "fun-attracted", "friends-enjoy",
#>   "friends-hurt", "friends-part", and 17 more.
#> ℹ The polychoric asymptotic covariance (and any DWLS weights or robust standard
#>   errors derived from it) can be unreliable for such structurally sparse cells;
#>   interpret them with caution and consider collapsing rare response categories
#>   in these variables.
mod
#> 
#> EFA performed with estimator = 'DWLS' and rotation = 'none'.
#> 
#> ── Unrotated Loadings ──────────────────────────────────────────────────────────
#> 
#>             F1    h2    u2
#> fun        .818  .669  .331
#> friends    .855  .731  .269
#> enjoy      .893  .797  .203
#> hurt       .775  .601  .399
#> part       .824  .679  .321
#> commonly   .843  .711  .289
#> chances    .817  .668  .332
#> attracted  .859  .738  .262
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                 F1
#> SS loadings   5.594
#> Prop Tot Var   .699
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> scaled χ²(20) = 237.21, p < .001
#> CFI: .99
#> TLI: .99
#> RMSEA [90% CI]: .12 [.10; .13]
#> AIC: NA
#> BIC: NA
#> CAF: .49
#> SRMR: .02
summary(mod)
#> 
#> EFA performed with estimator = 'DWLS' and rotation = 'none'.
#> 
#> ── Model Diagnostics ───────────────────────────────────────────────────────────
#> 
#> Factors: 1
#> Variables: 8
#> N: 810
#> Heywood cases: 0
#> Cross-loading items (|loading| >= .300): 0
#> Items without salient loading (|loading| >= .300): 0
#> Factors with fewer than 3 salient indicators: 0
#> Items with primary-loading gap < .200: 0
#> Largest |residual|: .038
#> 
#> ── Unrotated Loadings ──────────────────────────────────────────────────────────
#> 
#>             F1    h2    u2
#> fun        .818  .669  .331
#> friends    .855  .731  .269
#> enjoy      .893  .797  .203
#> hurt       .775  .601  .399
#> part       .824  .679  .321
#> commonly   .843  .711  .289
#> chances    .817  .668  .332
#> attracted  .859  .738  .262
#> 
#> ── 95% Wald CIs for salient unrotated loadings ─────────────────────────────────
#> 
#> Variable   Factor  est    lower  upper
#> fun        F1       .818   .798   .838
#> friends    F1       .855   .835   .875
#> enjoy      F1       .893   .875   .910
#> hurt       F1       .775   .750   .800
#> part       F1       .824   .802   .846
#> commonly   F1       .843   .825   .861
#> chances    F1       .817   .795   .840
#> attracted  F1       .859   .842   .877
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                 F1
#> SS loadings   5.594
#> Prop Tot Var   .699
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> scaled χ²(20) = 237.21, p < .001
#> CFI: .99
#> TLI: .99
#> RMSEA [90% CI]: .12 [.10; .13]
#> AIC: NA
#> BIC: NA
#> CAF: .49
#> SRMR: .02
#> 
#> Note: Wald CIs from the robust (Godambe) sandwich covariance.
#> 
#> ── Residual Diagnostics ────────────────────────────────────────────────────────
#> 
#> Residual cutoff: |r| > .100
#> Number of large residuals: 0
#> Largest absolute residual: .038
#> 
#> No absolute residuals > .100 occurred.
#> 
#> Inspect the residual matrix for details (e.g., with residuals()).
```

#### Correlation Input

When you do not have raw data, you can enter a correlation matrix and
sample size instead. With ML estimation, you can still get
information-based SEs (from the expected information matrix), but these
assume multivariate normality.

``` r


# ML estimation with oblimin rotation and information SEs, based on correlation
# matrix and N
mod <- efa_fit(test_models$baseline$cormat, N = 500,  n_factors = 3, estimator = "ml",
           rotation = "oblimin", se = "information")
mod
#> 
#> EFA performed with estimator = 'ML' and rotation = 'oblimin'.
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>        F1     F2     F3    h2    u2
#> V1   -.036   .043   .607  .373  .627
#> V2    .013   .087   .458  .274  .726
#> V3    .074   .074   .430  .280  .720
#> V4    .111   .007   .536  .379  .621
#> V5    .164   .005   .418  .290  .710
#> V6   -.055  -.036   .687  .402  .598
#> V7    .017   .524   .095  .355  .645
#> V8   -.003   .562   .044  .345  .655
#> V9    .044   .535   .017  .328  .672
#> V10  -.019   .661  -.051  .385  .615
#> V11   .030   .352   .230  .296  .704
#> V12   .034   .649  -.015  .437  .563
#> V13   .612   .095  -.068  .397  .603
#> V14   .540  -.053   .086  .320  .680
#> V15   .552   .137  -.065  .363  .637
#> V16   .550  -.039   .092  .345  .655
#> V17   .652  -.035  -.013  .390  .610
#> V18   .549   .012   .052  .349  .651
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1     F2     F3
#> F1  1.000
#> F2   .591  1.000
#> F3   .621   .596  1.000
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F2     F3
#> SS loadings        2.225  2.088  1.994
#> Prop Tot Var        .124   .116   .111
#> Cum Prop Tot Var    .124   .240   .350
#> Prop Comm Var       .353   .331   .316
#> Cum Prop Comm Var   .353   .684  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> χ²(102) = 123.75, p = .070
#> CFI: .99
#> TLI: .98
#> RMSEA [90% CI]: .02 [.00; .03]
#> AIC: -80.25
#> BIC: -510.14
#> ECVI: 0.52
#> CAF: .50
#> SRMR: .03
summary(mod)
#> 
#> EFA performed with estimator = 'ML' and rotation = 'oblimin'.
#> 
#> ── Model Diagnostics ───────────────────────────────────────────────────────────
#> 
#> Factors: 3
#> Variables: 18
#> N: 500
#> Rotation local optima: 1 distinct from 6 of 101 starts
#> Heywood cases: 0
#> Cross-loading items (|loading| >= .300): 0
#> Items without salient loading (|loading| >= .300): 0
#> Factors with fewer than 3 salient indicators: 0
#> Items with primary-loading gap < .200: 1
#> Largest |residual|: .069
#> Factor intercorrelations > .85: none
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>        F1     F2     F3    h2    u2
#> V1   -.036   .043   .607  .373  .627
#> V2    .013   .087   .458  .274  .726
#> V3    .074   .074   .430  .280  .720
#> V4    .111   .007   .536  .379  .621
#> V5    .164   .005   .418  .290  .710
#> V6   -.055  -.036   .687  .402  .598
#> V7    .017   .524   .095  .355  .645
#> V8   -.003   .562   .044  .345  .655
#> V9    .044   .535   .017  .328  .672
#> V10  -.019   .661  -.051  .385  .615
#> V11   .030   .352   .230  .296  .704
#> V12   .034   .649  -.015  .437  .563
#> V13   .612   .095  -.068  .397  .603
#> V14   .540  -.053   .086  .320  .680
#> V15   .552   .137  -.065  .363  .637
#> V16   .550  -.039   .092  .345  .655
#> V17   .652  -.035  -.013  .390  .610
#> V18   .549   .012   .052  .349  .651
#> 
#> ── 95% Wald CIs for salient rotated loadings ───────────────────────────────────
#> 
#> Variable  Factor  est    lower  upper
#> V13       F1       .612   .488   .736
#> V14       F1       .540   .411   .669
#> V15       F1       .552   .423   .681
#> V16       F1       .550   .421   .679
#> V17       F1       .652   .535   .769
#> V18       F1       .549   .420   .679
#> V7        F2       .524   .394   .653
#> V8        F2       .562   .437   .687
#> V9        F2       .535   .407   .662
#> V10       F2       .661   .548   .773
#> V11       F2       .352   .211   .493
#> V12       F2       .649   .529   .770
#> V1        F3       .607   .474   .739
#> V2        F3       .458   .313   .604
#> V3        F3       .430   .282   .578
#> V4        F3       .536   .395   .677
#> V5        F3       .418   .269   .567
#> V6        F3       .687   .570   .805
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1     F2     F3
#> F1  1.000
#> F2   .591  1.000
#> F3   .621   .596  1.000
#> 
#> ── 95% Wald CIs for factor intercorrelations ───────────────────────────────────
#> 
#> Factors   est    lower  upper
#> F1 ~~ F2   .591   .499   .683
#> F1 ~~ F3   .621   .531   .712
#> F2 ~~ F3   .596   .503   .690
#> 
#> ── Structure Matrix ────────────────────────────────────────────────────────────
#> 
#>       F1    F2    F3
#> V1   .366  .384  .610
#> V2   .349  .367  .518
#> V3   .385  .374  .520
#> V4   .448  .392  .609
#> V5   .427  .351  .523
#> V6   .350  .341  .631
#> V7   .385  .590  .418
#> V8   .357  .587  .377
#> V9   .371  .571  .363
#> V10  .339  .619  .331
#> V11  .381  .507  .459
#> V12  .409  .661  .393
#> V13  .626  .416  .369
#> V14  .562  .317  .390
#> V15  .593  .425  .360
#> V16  .584  .341  .410
#> V17  .623  .343  .371
#> V18  .589  .368  .401
#> 
#> ── Simple Structure Diagnostics ────────────────────────────────────────────────
#> 
#> Items with primary-loading gap < .200:
#> • V11: F2 = .352, F3 = .230
#> 
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F2     F3
#> SS loadings        2.225  2.088  1.994
#> Prop Tot Var        .124   .116   .111
#> Cum Prop Tot Var    .124   .240   .350
#> Prop Comm Var       .353   .331   .316
#> Cum Prop Comm Var   .353   .684  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> χ²(102) = 123.75, p = .070
#> CFI: .99
#> TLI: .98
#> RMSEA [90% CI]: .02 [.00; .03]
#> AIC: -80.25
#> BIC: -510.14
#> ECVI: 0.52
#> CAF: .50
#> SRMR: .03
#> 
#> ── Residual Diagnostics ────────────────────────────────────────────────────────
#> 
#> Residual cutoff: |r| > .100
#> Number of large residuals: 0
#> Largest absolute residual: .069
#> 
#> No absolute residuals > .100 occurred.
#> 
#> Inspect the residual matrix for details (e.g., with residuals()).
```

## Citation

If you use this package in your research, please acknowledge it by
citing:

Steiner, M.D., & Grieder, S.G. (2020). EFAtools: An R package with fast
and flexible implementations of exploratory factor analysis tools.
*Journal of Open Source Software*, *5*(53), 2521.
<https://doi.org/10.21105/joss.02521>

## Contribute or Report Bugs

If you want to contribute or report bugs, please open an issue on GitHub
or email us at <markus.d.steiner@gmail.com> or
<silvia.steiner.grieder@gmail.com>.
