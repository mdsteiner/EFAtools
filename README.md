
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EFAtools

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/EFAtools)](https://CRAN.R-project.org/package=EFAtools)
[![R-CMD-check](https://github.com/mdsteiner/EFAtools/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mdsteiner/EFAtools/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/mdsteiner/EFAtools/graph/badge.svg)](https://app.codecov.io/gh/mdsteiner/EFAtools)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02521/status.svg)](https://doi.org/10.21105/joss.02521)

<!-- badges: end -->

The EFAtools package provides functions to perform exploratory factor
analysis (EFA) procedures and compare their solutions. The goal is to
provide state-of-the-art factor retention methods and a high degree of
flexibility in the EFA procedures. This way, implementations from R
psych and SPSS can be compared. Moreover, functions for Schmid-Leiman
transformation, and computation of omegas are provided. To speed up the
analyses, some of the iterative procedures like principal axis factoring
(PAF) are implemented in C++.

## Installation

You can install the release version from
[CRAN](https://cran.r-project.org/) with:

``` r
install.packages("EFAtools")
```

You can install the development version from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("mdsteiner/EFAtools")
```

To also build the vignette when installing the development version, use:

``` r
install.packages("devtools")
devtools::install_github("mdsteiner/EFAtools", build_vignettes = TRUE)
```

## Example

Here are a few examples on how to perform the analyses with the
different types and how to compare the results using the `COMPARE`
function. For more details, see the vignette by running
`vignette("EFAtools", package = "EFAtools")`. The vignette provides a
high-level introduction into the functionalities of the package.

``` r
# load the package
library(EFAtools)

# Run multiple factor retention methods
N_FACTORS(test_models$baseline$cormat, N = 500)
#> Warning: `x` is a correlation matrix, but "CD" needs raw data.
#> ℹ Skipping "CD".
#> ── Tests for the suitability of the data for factor analysis ───────────────────
#> 
#> ✔ The Bartlett's test of sphericity was significant at an alpha level of .05:
#>   χ²(153) = 2173.28, p < .001. These data are probably suitable for factor
#>   analysis.
#> ✔ The Kaiser-Meyer-Olkin criterion is marvellous (KMO = 0.916). These data are
#>   probably suitable for factor analysis.
#> 
#> ── Number of factors suggested by the factor retention criteria ────────────────
#> 
#> Empirical Kaiser Criterion
#> • Original implementation (Braeken & van Assen, 2017): 3
#> 
#> Hull method
#> • CAF: 3
#> • CFI: 1
#> • RMSEA: 1
#> 
#> Minimum average partial
#> • Original implementation (TR2): 1
#> • Revised implementation (TR4): 3
#> 
#> Next Eigenvalue Sufficiency Test
#> • Suggested number of factors: 3
#> 
#> Parallel analysis
#> • SMC eigenvalues: 3
#> 
#> ── Criteria that could not be run ──────────────────────────────────────────────
#> 
#> ! CD: needs raw data, but a correlation matrix was supplied

# A type SPSS EFA to mimick the SPSS implementation with
# promax rotation
EFA_SPSS <- EFA(test_models$baseline$cormat, n_factors = 3, type = "SPSS",
                  rotation = "promax")

# look at solution
EFA_SPSS
#> 
#> EFA performed with type = 'SPSS', method = 'PAF', and rotation = 'promax'.
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>        F1     F2     F3    h2    u2
#> V1   -.048   .035   .613  .367  .633
#> V2   -.001   .067   .482  .277  .723
#> V3    .060   .056   .453  .283  .717
#> V4    .101  -.009   .551  .378  .622
#> V5    .157  -.018   .438  .293  .707
#> V6   -.072  -.049   .704  .399  .601
#> V7    .001   .533   .093  .357  .643
#> V8   -.016   .581   .030  .349  .651
#> V9    .038   .550  -.001  .330  .670
#> V10  -.022   .674  -.071  .383  .617
#> V11   .015   .356   .232  .297  .703
#> V12   .020   .651  -.010  .432  .568
#> V13   .614   .086  -.067  .394  .606
#> V14   .548  -.068   .088  .322  .678
#> V15   .561   .128  -.070  .363  .637
#> V16   .555  -.050   .091  .344  .656
#> V17   .664  -.037  -.027  .390  .610
#> V18   .555   .004   .050  .350  .650
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1     F2     F3
#> F1  1.000
#> F2   .617  1.000
#> F3   .648   .632  1.000
#> 
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F2     F3
#> SS loadings        2.198  2.074  2.034
#> Prop Tot Var        .122   .115   .113
#> Cum Prop Tot Var    .122   .237   .350
#> Prop Comm Var       .349   .329   .323
#> Cum Prop Comm Var   .349   .677  1.000
#> 
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> CAF  : .50
#> RMSR  : .03
#> SRMR  : .02
#> df: 102

# A type psych EFA to mimick the psych::fa() implementation with
# promax rotation
EFA_psych <- EFA(test_models$baseline$cormat, n_factors = 3, type = "psych",
                  rotation = "promax")

# compare the type psych and type SPSS implementations
COMPARE(EFA_SPSS$rot_loadings, EFA_psych$rot_loadings,
        x_labels = c("SPSS", "psych"))
#> Mean [min, max] absolute difference:  0.0090 [ 0.0001,  0.0245]
#> Median absolute difference:  0.0095
#> Max decimals where all numbers are equal: 0
#> Minimum number of decimals provided: 17
#> 
#>        F1      F2      F3
#> V1    .0150   .0142  -.0195
#> V2    .0109   .0109  -.0138
#> V3    .0095   .0103  -.0119
#> V4    .0118   .0131  -.0154
#> V5    .0084   .0105  -.0109
#> V6    .0183   .0169  -.0245
#> V7   -.0026  -.0017   .0076
#> V8   -.0043  -.0035   .0102
#> V9   -.0055  -.0040   .0117
#> V10  -.0075  -.0066   .0151
#> V11   .0021   .0029   .0001
#> V12  -.0064  -.0050   .0136
#> V13  -.0109  -.0019   .0163
#> V14  -.0049   .0028   .0070
#> V15  -.0107  -.0023   .0161
#> V16  -.0051   .0028   .0074
#> V17  -.0096  -.0001   .0136
#> V18  -.0066   .0014   .0098

# Run EFA and compute standard errors and confidence intervals based
# on non-parametric bootstrap:
EFA(GRiPS_raw, n_factors = 1, method = "ML", se = "np-boot")
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
#> 
#> EFA performed with type = 'EFAtools', method = 'ML', and rotation = 'none'.
#> 
#> ── Unrotated Loadings ──────────────────────────────────────────────────────────
#> 
#>             F1    h2    u2
#> fun        .796  .634  .366
#> friends    .851  .725  .275
#> enjoy      .872  .760  .240
#> hurt       .767  .588  .412
#> part       .817  .667  .333
#> commonly   .832  .692  .308
#> chances    .788  .621  .379
#> attracted  .845  .715  .285
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                 F1
#> SS loadings   5.401
#> Prop Tot Var   .675
#> 
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> χ²(20) = 59.61, p < .001
#> CFI  [95% bootstrap-CI] : .99 [.98, .99]
#> TLI  : .99
#> RMSEA [90% CI]  [95% bootstrap-CI] : .05 [.04; .06] [.05, .09]
#> AIC  [95% bootstrap-CI] : 19.61 [13.43, 106.96]
#> BIC  [95% bootstrap-CI] : -74.33 [-80.51, 13.02]
#> ECVI  : 0.11
#> CAF  [95% bootstrap-CI] : .50 [.48, .51]
#> RMSR  [95% bootstrap-CI] : .02 [.02, .03]
#> SRMR  : .01
#> 
#> Note: Bootstrap CIs based on 1000 bootstrap samples.


# Average solution across many different EFAs with oblique rotations
EFA_AV <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                      method = c("PAF", "ML", "ULS"), rotation = "oblique",
                      show_progress = FALSE)

# look at solution
EFA_AV
#> 
#> Averaging performed with averaging method mean (trim = 0) across 162 EFAs, varying the following settings: method, init_comm, criterion_type, start_method, rotation, k_promax, P_type, and varimax_type.
#> 
#> The error rate is at 0%. Of the solutions that did not result in an error, 100% converged, 0% contained Heywood cases, and 100% were admissible.
#> 
#> 
#> ══ Indicator-to-Factor Correspondences ═════════════════════════════════════════
#> 
#> For each cell, the proportion of solutions including the respective indicator-to-factor correspondence. A salience threshold of 0.3 was used to determine indicator-to-factor correspondences.
#> 
#>       F1    F2    F3
#> V1    .00   .01  1.00
#> V2    .00   .01  1.00
#> V3    .00   .01  1.00
#> V4    .00   .01  1.00
#> V5    .00   .01  1.00
#> V6    .00   .01  1.00
#> V7    .00  1.00   .09
#> V8    .00  1.00   .09
#> V9    .00  1.00   .09
#> V10   .00  1.00   .09
#> V11   .00   .91   .09
#> V12   .00  1.00   .09
#> V13  1.00   .00   .06
#> V14   .94   .00   .06
#> V15   .99   .00   .06
#> V16   .94   .00   .06
#> V17  1.00   .00   .06
#> V18   .94   .00   .06
#> 
#> 
#> ══ Loadings ════════════════════════════════════════════════════════════════════
#> 
#> ── Mean ────────────────────────────────────────────────────────────────────────
#> 
#>        F1     F2     F3
#> V1   -.036   .052   .592
#> V2    .008   .082   .468
#> V3    .065   .072   .444
#> V4    .102   .015   .537
#> V5    .154   .006   .431
#> V6   -.059  -.025   .669
#> V7    .019   .507   .141
#> V8    .002   .548   .086
#> V9    .051   .521   .058
#> V10  -.004   .635  -.003
#> V11   .028   .345   .259
#> V12   .038   .618   .052
#> V13   .580   .103  -.008
#> V14   .514  -.040   .122
#> V15   .529   .141  -.010
#> V16   .521  -.024   .127
#> V17   .622  -.013   .025
#> V18   .522   .026   .093
#> 
#> 
#> ── Range ───────────────────────────────────────────────────────────────────────
#> 
#>       F1    F2    F3
#> V1   .300  .630  .159
#> V2   .245  .492  .142
#> V3   .263  .436  .181
#> V4   .343  .506  .193
#> V5   .309  .390  .212
#> V6   .336  .703  .174
#> V7   .092  .254  .460
#> V8   .078  .218  .494
#> V9   .069  .177  .474
#> V10  .098  .207  .596
#> V11  .147  .327  .332
#> V12  .082  .210  .574
#> V13  .308  .283  .655
#> V14  .346  .211  .519
#> V15  .283  .253  .633
#> V16  .350  .210  .526
#> V17  .366  .300  .626
#> V18  .331  .206  .556
#> 
#> 
#> 
#> ══ Factor Intercorrelations from Oblique Solutions ═════════════════════════════
#> 
#> ── Mean ────────────────────────────────────────────────────────────────────────
#> 
#>       F1     F2     F3
#> F1  1.000
#> F2   .503  1.000
#> F3   .544   .500  1.000
#> 
#> 
#> ── Range ───────────────────────────────────────────────────────────────────────
#> 
#>      F1    F2    F3
#> F1  .000
#> F2  .827  .000
#> F3  .712  .995  .000
#> 
#> 
#> 
#> ══ Variances Accounted for ═════════════════════════════════════════════════════
#> 
#> ── Mean ────────────────────────────────────────────────────────────────────────
#> 
#>                  F1     F2     F3
#> SS loadings    2.105  1.993  2.209
#> Prop Tot Var    .117   .111   .123
#> Prop Comm Var   .334   .316   .350
#> 
#> 
#> ── Range ───────────────────────────────────────────────────────────────────────
#> 
#>                  F1     F2     F3
#> SS loadings    1.761  2.272  3.887
#> Prop Tot Var    .098   .126   .216
#> Prop Comm Var   .279   .360   .616
#> 
#> 
#> 
#> ══ Model Fit ═══════════════════════════════════════════════════════════════════
#> 
#>        M (SD) [Min; Max]
#> 𝜒²: 123.80 ( 0.08) [123.75; 123.92]
#> df: 102
#> p: .070 (.001) [.069; .070]
#> CFI: .99 (.00) [.99; .99]
#> TLI: .98 (.00) [.98; .98]
#> RMSEA: .02 (.00) [.02; .02]
#> AIC: -80.20 ( 0.08) [-80.25; -80.08]
#> BIC: -510.09 ( 0.08) [-510.14; -509.97]
#> ECVI:  0.52 ( 0.00) [ 0.52;  0.52]
#> CAF: .50 (.00) [.50; .50]
#> RMSR: .03 (.00) [.03; .03]
#> SRMR: .02 (.00) [.02; .03]

# Perform a Schmid-Leiman transformation
SL <- SL(EFA_psych)

# Based on a specific salience threshold for the loadings (here: .20):
factor_corres <- SL$sl[, c("F1", "F2", "F3")] >= .2

# Compute omegas from the Schmid-Leiman solution
OMEGA(SL, factor_corres = factor_corres)
#> Omega total, omega hierarchical, omega subscale, H index, explained common variance (ECV), and percent of uncontaminated correlations (PUC) for the general factor (top row) and omegas and H index for the group factors:
#> 
#>      tot  hier   sub     H   ECV   PUC
#> g  0.883 0.750 0.122 0.845 0.661 0.706
#> F1 0.769 0.498 0.272 0.470            
#> F2 0.764 0.494 0.270 0.477            
#> F3 0.745 0.543 0.202 0.391
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
<silvia.grieder@gmail.com>.
