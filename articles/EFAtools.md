# EFAtools

This vignette provides an overview for the functionalities of the
EFAtools package. The general aim of the package is to provide flexible
implementations of different algorithms for an exploratory factor
analyses (EFA) procedure, including factor retention methods, factor
extraction and rotation methods, as well as the computation of a
Schmid-Leiman solution and McDonald’s omega coefficients.

The package was first designed to enable a comparison of EFA
(specifically, principal axis factoring with subsequent promax rotation)
performed in R using the
[**psych**](https://CRAN.R-project.org/package=psych) package and EFA
performed in SPSS. That is why some functions allow the specification of
a type, including `"psych"` and `"SPSS"`, such that the respective
procedure will be executed to match the output of these implementations
(which do not always lead to the same results; see separate vignette
[**Replicate_SPSS_psych**](https://mdsteiner.github.io/EFAtools/articles/Replicate_SPSS_psych.md "Replicate SPSS and R psych results with EFAtools")
for a demonstration of the replication of original results). This
vignette will go through a complete example, that is, we will first show
how to determine the number of factors to retain, then perform different
factor extraction methods, run a Schmid-Leiman transformation and
compute omegas.

The package can be installed from CRAN using
`install.packages("EFAtools")`, or from GitHub using
`devtools::install_github("mdsteiner/EFAtools")`, and then loaded using:

``` r

library(EFAtools)
```

In this vignette, we will use the `DOSPERT_raw` dataset, which contains
responses to the Domain Specific Risk Taking Scale (DOSPERT) of 3123
participants. The dataset is contained in the `EFAtools` package, for
details, see
[`?DOSPERT_raw`](https://mdsteiner.github.io/EFAtools/reference/DOSPERT_raw.md).
Note that this vignette is to provide a general overview and it is
beyond its scope to explain all methods and functions in detail. If you
want to learn more on the details and methods, please see the respective
help functions for explanations and literature references. However, the
dataset is rather large, so, just to save time when building the
vignette, we will only use the first 500 observations. When you normally
do your analyses, you use the full dataset.

``` r

# only use a subset to make analyses faster
DOSPERT_sub <- DOSPERT_raw[1:500,]
```

## Test Suitability of Data

The first step in an EFA procedure is to test whether your data is
suitable for factor analysis. To this end, the `EFAtools` package
provides the
[`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md)
and the
[`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md)
functions. The Bartlett’s test of sphericity tests whether a correlation
matrix is significantly different from an identity matrix (a correlation
matrix with zero correlations between all variables). This test should
thus be significant. The Kaiser-Meyer-Olkin criterion (KMO) represents
the degree to which each observed variable is predicted by the other
variables in the dataset and thus is another indicator for how
correlated the different variables are.

We can test whether our `DOSPERT_sub` dataset is suitable for factor
analysis as follows.

``` r

# Bartlett's test of sphericity
efa_bartlett(DOSPERT_sub)
#> ℹ `x` is not a correlation matrix; computing correlations from
#> the raw data.
#> 
#> ✔ The Bartlett's test of sphericity was significant at an alpha level of .05.
#> These data are probably suitable for factor analysis.
#> 
#> 𝜒²(435) = 5843.05, p < .001

# KMO criterion
efa_kmo(DOSPERT_sub)
#> ℹ `x` is not a correlation matrix; computing correlations from
#> the raw data.
#> 
#> ── Kaiser-Meyer-Olkin criterion (KMO) ──────────────────────────────────────────
#> 
#> ✔ The overall KMO value for your data is meritorious.
#> These data are probably suitable for factor analysis.
#> 
#> Overall: 0.87
#> 
#> For each variable:
#> ethR_1 ethR_2 ethR_3 ethR_4 ethR_5 ethR_6 finR_1 finR_2 finR_3 finR_4 finR_5 
#>  0.908  0.927  0.906  0.818  0.926  0.837  0.880  0.816  0.853  0.857  0.853 
#> finR_6 heaR_1 heaR_2 heaR_3 heaR_4 heaR_5 heaR_6 recR_1 recR_2 recR_3 recR_4 
#>  0.861  0.912  0.898  0.892  0.889  0.882  0.935  0.888  0.927  0.905  0.849 
#> recR_5 recR_6 socR_1 socR_2 socR_3 socR_4 socR_5 socR_6 
#>  0.840  0.920  0.714  0.784  0.735  0.791  0.840  0.786
```

Note that these tests can also be run in the
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
function.

## Factor Retention Methods

As the goal of EFA is to determine the underlying factors from a set of
multiple variables, one of the most important decisions is how many
factors can or should be extracted. There exists a plethora of factor
retention methods to use for this decision. The problem is that there is
no method that consistently outperforms all other methods. Rather, which
factor retention method to use depends on the structure of the data: are
there few or many indicators, are factors strong or weak, are the factor
intercorrelations weak or strong. For rules on which methods to use,
see, for example, [Auerswald and Moshagen,
(2019)](https://doi.org/10.1037/met0000200).

There are multiple factor retention methods implemented in the
`EFAtools` package. They can either be called with separate functions,
or all (or a selection) of them using the
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
function.

### Calling Separate Functions

Let’s first look at how to determine the number of factors to retain by
calling separate functions. For example, if you would like to perform a
parallel analysis based on squared multiple correlations (SMC; sometimes
also called a parallel analysis with principal factors), you can do the
following:

``` r

# determine the number of factors to retain using parallel analysis
efa_parallel(DOSPERT_sub, eigen_type = "SMC")
#> ℹ `x` is not a correlation matrix; computing correlations from
#> the raw data.
#> ── Parallel analysis ───────────────────────────────────────────────────────────
#> Eigenvalues found using SMC; 1000 simulated datasets.
#> 
#> • SMC eigenvalues: 10
#> 
#> ℹ Number of factors retained using the "means" decision rule.
```

Generating the plot can also be suppressed if the output is printed
explicitly:

``` r

# determine the number of factors to retain using parallel analysis
print(efa_parallel(DOSPERT_sub, eigen_type = "SMC"), plot = FALSE)
#> ℹ `x` is not a correlation matrix; computing correlations from
#> the raw data.
#> ── Parallel analysis ───────────────────────────────────────────────────────────
#> Eigenvalues found using SMC; 1000 simulated datasets.
#> 
#> • SMC eigenvalues: 10
#> 
#> ℹ Number of factors retained using the "means" decision rule.
```

Other factor retention methods can be used accordingly. For example, to
use the empirical Kaiser criterion, use the `efa_ekc` function:

``` r

# determine the number of factors to retain using the empirical Kaiser criterion
print(efa_ekc(DOSPERT_sub), plot = FALSE)
#> ℹ `x` is not a correlation matrix; computing correlations from
#> the raw data.
#> ── Empirical Kaiser Criterion ──────────────────────────────────────────────────
#> 
#> • Original implementation (Braeken & van Assen, 2017): 8
#> 
#> ℹ Multiple implementations of EKC exist; make sure to report which one you used (see the efa_ekc help page for details).
```

The following factor retention methods are currently implemented:
comparison data
([`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md)),
empirical Kaiser criterion
([`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md)),
the hull method
([`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)),
the Kaiser-Guttman criterion
([`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md)),
the minimum average partial test
([`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md)),
the next eigenvalue sufficiency test
([`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md)),
parallel analysis
([`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)),
scree test
([`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md)),
and sequential model tests
([`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md)).
Many of these functions have multiple versions of the respective factor
retention method implemented, for example, the parallel analysis can be
done based on eigenvalues found using unity (principal components) or
SMCs, or on an EFA procedure. Another example is the hull method, which
can be used with different fitting methods (principal axis factoring
\[PAF\], maximum likelihood \[ML\], or unweighted least squares
\[ULS\]), and different goodness of fit indices. Please see the
respective function documentations for details.

### Run Multiple Factor Retention Methods With `efa_retain()`

If you want to use multiple factor retention methods, for example, to
compare whether different methods suggest the same number of factors, it
is easier to use the
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
function. This is a wrapper around all the implemented factor retention
methods. Moreover, it also enables to run the Bartlett’s test of
sphericity and compute the KMO criterion.

For example, to test the suitability of the data for factor analysis and
to determine the number of factors to retain based on parallel analysis
(but only using eigen values based on SMCs and PCA), the EKC, and the
sequential model test, we can run the following code:

``` r

efa_retain(DOSPERT_sub, criteria = c("PARALLEL", "EKC", "SMT"),
           eigen_type_other = c("SMC", "PCA"))
#> ── Tests for the suitability of the data for factor analysis ───────────────────
#> 
#> ✔ The Bartlett's test of sphericity was significant at an alpha level of .05:
#>   χ²(435) = 5843.05, p < .001. These data are probably suitable for factor
#>   analysis.
#> ✔ The Kaiser-Meyer-Olkin criterion is meritorious (KMO = 0.87). These data are
#>   probably suitable for factor analysis.
#> 
#> ── Number of factors suggested by the factor retention criteria ────────────────
#> 
#> Empirical Kaiser Criterion
#> • Original implementation (Braeken & van Assen, 2017): 8
#> 
#> Parallel analysis
#> • PCA eigenvalues: 5
#> • SMC eigenvalues: 10
#> 
#> Sequential model tests
#> • Sequential chi-square model tests: 13
#> • Lower bound of RMSEA 90% CI: 6
#> • Akaike Information Criterion: 13
```

If all possible factor retention methods should be used, it is
sufficient to provide the data object (note that this takes a while, as
the comparison data is computationally expensive and therefore
relatively slow method, especially if larger datasets are used). We
additionally specify the method argument to use unweighted least squares
(ULS) estimation. This is a bit faster than using principle axis
factoring (PAF) and it enables the computation of more goodness of fit
indices:

``` r

efa_retain(DOSPERT_sub, method = "ULS")
#> ── Tests for the suitability of the data for factor analysis ───────────────────
#> 
#> ✔ The Bartlett's test of sphericity was significant at an alpha level of .05:
#>   χ²(435) = 5843.05, p < .001. These data are probably suitable for factor
#>   analysis.
#> ✔ The Kaiser-Meyer-Olkin criterion is meritorious (KMO = 0.87). These data are
#>   probably suitable for factor analysis.
#> 
#> ── Number of factors suggested by the factor retention criteria ────────────────
#> 
#> Comparison data
#> • Suggested number of factors: 10
#> 
#> Empirical Kaiser Criterion
#> • Original implementation (Braeken & van Assen, 2017): 8
#> 
#> Hull method
#> • CAF: 1
#> • CFI: 1
#> • RMSEA: 1
#> 
#> Minimum average partial
#> • Original implementation (TR2): 6
#> • Revised implementation (TR4): 6
#> 
#> Next Eigenvalue Sufficiency Test
#> • Suggested number of factors: 10
#> 
#> Parallel analysis
#> • SMC eigenvalues: 10
```

Now, this is not the scenario one is happy about, but it still does
happen: There is no obvious convergence between the methods and thus the
choice of the number of factors to retain becomes rather difficult (and
to some extent arbitrary). We will proceed with 6 factors, as it is what
is typically used with DOSPERT data, but this does not mean that other
number of factors are not just as plausible.

Note that all factor retention methods, except comparison data (CD), can
also be used with correlation matrices. We use `method = "ULS"` and
`eigen_type_other = c("SMC", "PCA")` to skip the slower criteria. In
this case, the sample size has to be specified:

``` r

efa_retain(test_models$baseline$cormat, N = 500,
           method = "ULS", eigen_type_other = c("SMC", "PCA"))
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
#> • PCA eigenvalues: 3
#> • SMC eigenvalues: 3
#> 
#> ── Criteria that could not be run ──────────────────────────────────────────────
#> 
#> ! CD: needs raw data, but a correlation matrix was supplied
```

## Exploratory Factor Analysis: Factor Extraction

Multiple algorithms to perform an EFA and to rotate the found solutions
are implemented in the `EFAtools` package. All of them can be used using
the
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
function. To perform the EFA, you can use one of principal axis
factoring (PAF), maximum likelihood estimation (ML), and unweighted
least squares (ULS; also sometimes referred to as MINRES). To rotate the
solutions, the `EFAtools` package offers varimax and promax rotations,
as well as a range of other orthogonal and oblique rotations (e.g.,
quartimax, equamax, oblimin, quartimin, geomin, bentler, bifactor, and
simplimax), all computed by rotation engines built into the package.

You can run an EFA with PAF and no rotation like this:

``` r

efa_fit(DOSPERT_sub, n_factors = 6)
#> ℹ `x` is not a correlation matrix; computing correlations from
#> the raw data.
#> 
#> EFA performed with type = 'EFAtools', method = 'PAF', and rotation = 'none'.
#> 
#> ── Unrotated Loadings ──────────────────────────────────────────────────────────
#> 
#>          F1     F2     F3     F4     F5     F6    h2    u2
#> ethR_1  .558  -.193   .141   .232  -.129   .049  .442  .558
#> ethR_2  .472  -.157   .073   .242   .037   .059  .316  .684
#> ethR_3  .507  -.402   .164   .112  -.087   .062  .470  .530
#> ethR_4  .308  -.278   .189   .294  -.103   .224  .356  .644
#> ethR_5  .459  -.173   .202   .042  -.082   .045  .292  .708
#> ethR_6  .412  -.308   .230   .201  -.066   .207  .405  .595
#> finR_1  .602  -.304  -.461  -.045   .244   .028  .730  .270
#> finR_2  .294   .281  -.255  -.066  -.341   .086  .359  .641
#> finR_3  .626  -.285  -.471  -.031   .264   .114  .778  .222
#> finR_4  .557  -.012  -.357  -.023  -.437  -.063  .633  .367
#> finR_5  .657  -.342  -.491  -.064   .218   .028  .842  .158
#> finR_6  .551   .157  -.341  -.007  -.429  -.082  .635  .365
#> heaR_1  .464  -.079   .047   .133   .010   .046  .244  .756
#> heaR_2  .409  -.020   .163   .229   .040   .068  .252  .748
#> heaR_3  .497  -.101   .172   .127   .007  -.325  .409  .591
#> heaR_4  .547  -.110   .183   .038   .070  -.322  .455  .545
#> heaR_5  .382   .002   .175   .250   .028   .010  .240  .760
#> heaR_6  .504   .063   .102   .015   .006  -.196  .307  .693
#> recR_1  .443   .337   .057  -.064   .007  -.063  .322  .678
#> recR_2  .617  -.012   .147  -.141   .126  -.232  .493  .507
#> recR_3  .612   .178   .103  -.247   .080  -.153  .507  .493
#> recR_4  .634   .228   .254  -.447   .033   .235  .774  .226
#> recR_5  .640   .135   .267  -.463   .051   .261  .784  .216
#> recR_6  .629   .238   .135  -.244  -.060  -.014  .534  .466
#> socR_1  .075   .528  -.024   .233   .213   .124  .400  .600
#> socR_2  .358   .428  -.083   .247   .147   .013  .401  .599
#> socR_3  .118   .528  -.160   .122   .044   .116  .348  .652
#> socR_4  .235   .481  -.057   .214   .166  -.093  .372  .628
#> socR_5  .317   .330   .014   .154  -.030   .002  .234  .766
#> socR_6  .260   .462  -.058   .267  -.023   .082  .363  .637
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F2     F3     F4    F5     F6
#> SS loadings        7.013  2.417  1.529  1.246  .850   .641
#> Prop Tot Var        .234   .081   .051   .042  .028   .021
#> Cum Prop Tot Var    .234   .314   .365   .407  .435   .457
#> Prop Comm Var       .512   .177   .112   .091  .062   .047
#> Cum Prop Comm Var   .512   .689   .800   .891  .953  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> CAF: .48
#> SRMR: .03
#> df: 270
```

To rotate the loadings (e.g., using a promax rotation) adapt the
`rotation` argument:

``` r

efa_fit(DOSPERT_sub, n_factors = 6, rotation = "promax")
#> ℹ `x` is not a correlation matrix; computing correlations from
#> the raw data.
#> 
#> EFA performed with type = 'EFAtools', method = 'PAF', and rotation = 'promax'.
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>           F1     F4     F3     F6     F2     F5    h2    u2
#> ethR_1   .547  -.025  -.020   .130   .009   .122  .442  .558
#> ethR_2   .447  -.055   .132   .098   .115  -.043  .316  .684
#> ethR_3   .532   .045   .068   .117  -.212   .030  .470  .530
#> ethR_4   .695  -.035  -.039  -.162   .000  -.002  .356  .644
#> ethR_5   .382   .159  -.051   .130  -.100   .036  .292  .708
#> ethR_6   .664   .086  -.003  -.089  -.063  -.044  .405  .595
#> finR_1  -.015  -.019   .853   .022   .003   .008  .730  .270
#> finR_2  -.050   .121  -.044  -.211   .120   .599  .359  .641
#> finR_3   .046   .035   .895  -.094   .061  -.011  .778  .222
#> finR_4   .039  -.071   .098   .035  -.102   .768  .633  .367
#> finR_5  -.006  -.010   .887   .023  -.038   .057  .842  .158
#> finR_6  -.022  -.051   .023   .061   .031   .774  .635  .365
#> heaR_1   .318   .042   .111   .091   .092   .019  .244  .756
#> heaR_2   .416   .020  -.003   .093   .191  -.081  .252  .748
#> heaR_3   .133  -.140  -.071   .674  -.036  -.009  .409  .591
#> heaR_4   .074  -.031  -.006   .697  -.067  -.074  .455  .545
#> heaR_5   .385  -.037  -.053   .170   .197  -.072  .240  .760
#> heaR_6   .044   .062  -.026   .455   .059   .048  .307  .693
#> recR_1  -.067   .253  -.060   .219   .247   .094  .322  .678
#> recR_2  -.046   .231   .089   .574  -.046  -.090  .493  .507
#> recR_3  -.155   .402   .053   .423   .035   .011  .507  .493
#> recR_4   .042   .926  -.010  -.089   .002  -.017  .774  .226
#> recR_5   .088   .945   .034  -.114  -.065  -.058  .784  .216
#> recR_6  -.020   .510  -.067   .212   .053   .154  .534  .466
#> socR_1   .016   .029   .025  -.137   .679  -.162  .400  .600
#> socR_2   .057  -.035   .103   .072   .602  -.020  .401  .599
#> socR_3  -.084   .061   .029  -.199   .559   .110  .348  .652
#> socR_4  -.090  -.082   .023   .197   .583  -.057  .372  .628
#> socR_5   .097   .032  -.086   .077   .370   .117  .234  .766
#> socR_6   .134  -.035  -.072  -.079   .568   .142  .363  .637
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1     F2     F3     F4     F5     F6
#> F1  1.000
#> F2   .357  1.000
#> F3   .455   .349  1.000
#> F4   .567   .594   .502  1.000
#> F5   .031   .316   .046   .287  1.000
#> F6   .329   .441   .463   .486   .324  1.000
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F4     F3     F6     F2     F5
#> SS loadings        2.565  2.518  2.457  2.225  2.213  1.717
#> Prop Tot Var        .086   .084   .082   .074   .074   .057
#> Cum Prop Tot Var    .086   .169   .251   .326   .399   .457
#> Prop Comm Var       .187   .184   .179   .162   .162   .125
#> Cum Prop Comm Var   .187   .371   .551   .713   .875  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> CAF: .48
#> SRMR: .03
#> df: 270
```

This now performed PAF with promax rotation with the specification, on
average, we found to produce the most accurate results in a simulation
analysis (see function documentation). If you want to replicate the
implementation of the *psych* R package, you can supply the `"psych"`
preset through
[`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
and
[`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md):

``` r

efa_fit(DOSPERT_sub, n_factors = 6, rotation = "promax",
        estimate_control = estimate_control(type = "psych"),
        rotate_control = rotate_control(type = "psych"))
#> ℹ `x` is not a correlation matrix; computing correlations from
#> the raw data.
#> 
#> EFA performed with type = 'psych', method = 'PAF', and rotation = 'promax'.
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>           F1     F3     F4     F6     F2     F5    h2    u2
#> ethR_1   .554  -.024   .116  -.029  -.002   .136  .442  .558
#> ethR_2   .441   .124   .088  -.054   .119  -.035  .316  .684
#> ethR_3   .556   .063   .113   .052  -.227   .048  .470  .530
#> ethR_4   .704  -.042  -.195   .003  -.005  -.006  .356  .644
#> ethR_5   .395  -.052   .131   .156  -.111   .044  .292  .708
#> ethR_6   .677  -.007  -.109   .117  -.070  -.047  .405  .595
#> finR_1  -.002   .827   .040  -.018   .017   .029  .730  .270
#> finR_2  -.034  -.041  -.229   .116   .103   .593  .359  .641
#> finR_3   .057   .868  -.082   .050   .079   .000  .778  .222
#> finR_4   .078   .096   .026  -.100  -.134   .801  .633  .367
#> finR_5   .013   .860   .042  -.010  -.027   .081  .842  .158
#> finR_6   .005   .023   .051  -.088   .003   .801  .635  .365
#> heaR_1   .315   .105   .087   .038   .093   .024  .244  .756
#> heaR_2   .401  -.007   .082   .018   .197  -.083  .252  .748
#> heaR_3   .125  -.072   .700  -.210  -.044   .036  .409  .591
#> heaR_4   .068  -.009   .732  -.103  -.074  -.032  .455  .545
#> heaR_5   .367  -.056   .162  -.048   .201  -.068  .240  .760
#> heaR_6   .036  -.028   .477   .008   .054   .071  .307  .693
#> recR_1  -.086  -.059   .233   .215   .250   .089  .322  .678
#> recR_2  -.051   .085   .618   .170  -.049  -.065  .493  .507
#> recR_3  -.160   .051   .464   .348   .032   .019  .507  .493
#> recR_4   .051  -.010  -.063   .933  -.002  -.059  .774  .226
#> recR_5   .102   .033  -.087   .961  -.069  -.100  .784  .216
#> recR_6  -.018  -.065   .236   .476   .044   .147  .534  .466
#> socR_1  -.043   .022  -.155   .032   .713  -.200  .400  .600
#> socR_2   .007   .097   .062  -.057   .627  -.034  .401  .599
#> socR_3  -.123   .028  -.217   .061   .580   .079  .348  .652
#> socR_4  -.144   .020   .197  -.121   .608  -.065  .372  .628
#> socR_5   .070  -.085   .068   .010   .377   .109  .234  .766
#> socR_6   .094  -.073  -.104  -.044   .583   .122  .363  .637
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1     F2     F3     F4     F5     F6
#> F1  1.000
#> F2   .435  1.000
#> F3   .585   .463  1.000
#> F4   .333   .316   .637  1.000
#> F5   .149   .083   .387   .398  1.000
#> F6   .306   .422   .483   .501   .436  1.000
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F3     F4     F6     F2     F5
#> SS loadings        2.612  2.364  2.353  2.332  2.256  1.778
#> Prop Tot Var        .087   .079   .078   .078   .075   .059
#> Cum Prop Tot Var    .087   .166   .244   .322   .397   .457
#> Prop Comm Var       .191   .173   .172   .170   .165   .130
#> Cum Prop Comm Var   .191   .363   .535   .705   .870  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> CAF: .48
#> SRMR: .03
#> df: 270
```

If you want to use the *SPSS* implementation, you can supply the
`"SPSS"` preset through
[`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
and
[`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md):

``` r

efa_fit(DOSPERT_sub, n_factors = 6, rotation = "promax",
        estimate_control = estimate_control(type = "SPSS"),
        rotate_control = rotate_control(type = "SPSS"))
#> ℹ `x` is not a correlation matrix; computing correlations from
#> the raw data.
#> 
#> EFA performed with type = 'SPSS', method = 'PAF', and rotation = 'promax'.
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>           F1     F4     F3     F2     F6     F5    h2    u2
#> ethR_1   .547  -.025  -.020   .009   .130   .122  .442  .558
#> ethR_2   .447  -.055   .132   .115   .098  -.043  .316  .684
#> ethR_3   .532   .045   .068  -.212   .117   .030  .470  .530
#> ethR_4   .695  -.035  -.039   .000  -.162  -.002  .356  .644
#> ethR_5   .382   .159  -.051  -.100   .130   .036  .292  .708
#> ethR_6   .664   .086  -.003  -.063  -.089  -.044  .405  .595
#> finR_1  -.015  -.019   .853   .003   .022   .008  .730  .270
#> finR_2  -.050   .121  -.044   .120  -.211   .599  .359  .641
#> finR_3   .046   .035   .895   .061  -.094  -.011  .778  .222
#> finR_4   .039  -.071   .098  -.102   .035   .768  .633  .367
#> finR_5  -.006  -.010   .887  -.038   .023   .057  .842  .158
#> finR_6  -.022  -.051   .023   .031   .061   .774  .635  .365
#> heaR_1   .318   .042   .111   .092   .091   .019  .244  .756
#> heaR_2   .416   .020  -.003   .191   .093  -.081  .252  .748
#> heaR_3   .133  -.140  -.071  -.036   .674  -.009  .409  .591
#> heaR_4   .074  -.031  -.006  -.067   .697  -.074  .455  .545
#> heaR_5   .385  -.037  -.053   .197   .170  -.072  .240  .760
#> heaR_6   .044   .062  -.026   .059   .455   .048  .307  .693
#> recR_1  -.067   .253  -.060   .247   .219   .094  .322  .678
#> recR_2  -.046   .231   .089  -.046   .574  -.090  .493  .507
#> recR_3  -.155   .402   .053   .035   .423   .011  .507  .493
#> recR_4   .042   .926  -.010   .002  -.089  -.017  .774  .226
#> recR_5   .088   .945   .034  -.065  -.114  -.058  .784  .216
#> recR_6  -.020   .510  -.067   .053   .212   .154  .534  .466
#> socR_1   .016   .029   .025   .679  -.137  -.162  .400  .600
#> socR_2   .057  -.035   .103   .602   .072  -.020  .401  .599
#> socR_3  -.084   .061   .029   .559  -.199   .110  .348  .652
#> socR_4  -.090  -.082   .023   .583   .197  -.057  .372  .628
#> socR_5   .097   .032  -.086   .370   .077   .117  .234  .766
#> socR_6   .134  -.035  -.072   .568  -.079   .142  .363  .637
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1     F2     F3     F4     F5     F6
#> F1  1.000
#> F2   .357  1.000
#> F3   .455   .349  1.000
#> F4   .031   .316   .046  1.000
#> F5   .567   .594   .502   .287  1.000
#> F6   .329   .441   .463   .324   .486  1.000
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F4     F3     F2     F6     F5
#> SS loadings        2.565  2.518  2.457  2.213  2.225  1.717
#> Prop Tot Var        .086   .084   .082   .074   .074   .057
#> Cum Prop Tot Var    .086   .169   .251   .325   .399   .457
#> Prop Comm Var       .187   .184   .179   .162   .162   .125
#> Cum Prop Comm Var   .187   .371   .551   .712   .875  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> CAF: .48
#> SRMR: .03
#> df: 270
```

This enables comparisons of different implementations. The
[`efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/efa_compare.md)
function provides an easy way to compare how similar two loading
(pattern) matrices are:

``` r

efa_compare(
  efa_fit(DOSPERT_sub, n_factors = 6, rotation = "promax",
          estimate_control = estimate_control(type = "psych"),
          rotate_control = rotate_control(type = "psych"))$rot_loadings,
  efa_fit(DOSPERT_sub, n_factors = 6, rotation = "promax",
          estimate_control = estimate_control(type = "SPSS"),
          rotate_control = rotate_control(type = "SPSS"))$rot_loadings
)
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
#> Mean [min, max] absolute difference:  0.0159 [ 0.0001,  0.0712]
#> Median absolute difference:  0.0102
#> Max decimals where all numbers are equal: 0
#> Minimum number of decimals provided: 18
#> 
#>           F1      F3      F4      F6      F2      F5
#> ethR_1   .0078  -.0035  -.0138  -.0042  -.0106   .0148
#> ethR_2  -.0058  -.0078  -.0100   .0008   .0038   .0082
#> ethR_3   .0245  -.0052  -.0037   .0066  -.0146   .0186
#> ethR_4   .0090  -.0028  -.0325   .0375  -.0055  -.0045
#> ethR_5   .0128  -.0011   .0007  -.0035  -.0116   .0073
#> ethR_6   .0129  -.0038  -.0208   .0315  -.0067  -.0033
#> finR_1   .0131  -.0257   .0183   .0012   .0143   .0209
#> finR_2   .0162   .0025  -.0177  -.0048  -.0171  -.0055
#> finR_3   .0110  -.0272   .0121   .0150   .0182   .0107
#> finR_4   .0392  -.0023  -.0093  -.0286  -.0319   .0331
#> finR_5   .0191  -.0266   .0192   .0007   .0113   .0243
#> finR_6   .0265  -.0002  -.0097  -.0376  -.0279   .0271
#> heaR_1  -.0024  -.0062  -.0038  -.0042   .0009   .0053
#> heaR_2  -.0153  -.0039  -.0105  -.0011   .0059  -.0020
#> heaR_3  -.0075  -.0012   .0258  -.0703  -.0086   .0450
#> heaR_4  -.0063  -.0028   .0356  -.0712  -.0065   .0424
#> heaR_5  -.0182  -.0024  -.0087  -.0112   .0048   .0034
#> heaR_6  -.0085  -.0015   .0222  -.0537  -.0045   .0230
#> recR_1  -.0191   .0004   .0146  -.0382   .0028  -.0053
#> recR_2  -.0045  -.0046   .0439  -.0615  -.0025   .0251
#> recR_3  -.0054  -.0024   .0413  -.0546  -.0023   .0080
#> recR_4   .0091   .0001   .0262   .0078  -.0036  -.0421
#> recR_5   .0148  -.0012   .0270   .0154  -.0042  -.0418
#> recR_6   .0021   .0011   .0247  -.0342  -.0090  -.0068
#> socR_1  -.0595  -.0031  -.0180   .0032   .0346  -.0379
#> socR_2  -.0493  -.0061  -.0095  -.0217   .0255  -.0144
#> socR_3  -.0390  -.0017  -.0185   .0001   .0208  -.0311
#> socR_4  -.0543  -.0032  -.0005  -.0381   .0254  -.0083
#> socR_5  -.0270   .0003  -.0092  -.0213   .0071  -.0082
#> socR_6  -.0401  -.0005  -.0251  -.0089   .0156  -.0205
```

*Why would you want to do this?* One of us has had the experience that a
reviewer asked whether the results can be reproduced in another
statistical program than R. We therefore implemented this possibility in
the package for an easy application of large scale, systematic
comparisons.

Note that the `type` preset only affects the implementations of
principal axis factoring (PAF), varimax and promax rotations. The other
procedures are not affected (except the order of the rotated factors for
the other rotation methods).

As indicated previously, it is also possible to use different estimation
and rotation methods. For example, to perform an EFA with ULS and an
oblimin rotation, you can use the following code:

``` r

efa_fit(DOSPERT_sub, n_factors = 6, rotation = "oblimin", method = "ULS")
#> ℹ `x` is not a correlation matrix; computing correlations from
#> the raw data.
#> 
#> EFA performed with type = 'EFAtools', method = 'ULS', and rotation = 'oblimin'.
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>           F1     F2     F3     F4     F5     F6    h2    u2
#> ethR_1   .039   .009   .522  -.007   .137   .146  .441  .559
#> ethR_2   .166  -.025   .419   .100   .097  -.008  .316  .684
#> ethR_3   .120   .064   .511  -.220   .115   .052  .470  .530
#> ethR_4  -.013  -.038   .644  -.031  -.104   .008  .356  .644
#> ethR_5  -.003   .177   .371  -.103   .126   .059  .292  .708
#> ethR_6   .029   .084   .620  -.086  -.048  -.025  .405  .595
#> finR_1   .850  -.013  -.012   .003   .004   .028  .729  .271
#> finR_2  -.022   .102  -.039   .106  -.137   .551  .358  .642
#> finR_3   .883   .028   .039   .054  -.088   .007  .778  .222
#> finR_4   .160  -.054   .063  -.109   .068   .725  .634  .366
#> finR_5   .890  -.005  -.002  -.037   .006   .073  .842  .158
#> finR_6   .085  -.028   .005   .023   .092   .737  .636  .364
#> heaR_1   .144   .064   .302   .083   .091   .046  .244  .756
#> heaR_2   .029   .047   .389   .176   .095  -.043  .252  .748
#> heaR_3   .011  -.044   .154  -.016   .558   .051  .409  .591
#> heaR_4   .073   .060   .101  -.041   .571  -.008  .455  .545
#> heaR_5  -.016   .002   .362   .185   .158  -.031  .240  .760
#> heaR_6   .035   .123   .064   .074   .382   .092  .307  .693
#> recR_1  -.024   .282  -.050   .254   .193   .122  .322  .678
#> recR_2   .153   .295  -.014  -.016   .467  -.028  .493  .507
#> recR_3   .105   .441  -.118   .062   .348   .056  .507  .493
#> recR_4   .010   .881   .053   .014  -.060   .001  .774  .226
#> recR_5   .052   .897   .095  -.053  -.084  -.039  .786  .214
#> recR_6  -.016   .524   .004   .068   .191   .178  .534  .466
#> socR_1  -.004   .033  -.010   .656  -.104  -.134  .400  .600
#> socR_2   .109   .000   .043   .585   .074   .018  .401  .599
#> socR_3   .009   .054  -.095   .539  -.147   .110  .348  .652
#> socR_4   .031  -.034  -.091   .576   .169  -.015  .372  .628
#> socR_5  -.061   .059   .091   .359   .086   .137  .234  .766
#> socR_6  -.064  -.017   .114   .542  -.036   .155  .363  .637
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1     F2     F3     F4     F5     F6
#> F1  1.000
#> F2   .303  1.000
#> F3   .379   .272  1.000
#> F4   .017   .245   .030  1.000
#> F5   .367   .412   .436   .173  1.000
#> F6   .353   .353   .186   .282   .238  1.000
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F2     F3     F4     F5     F6
#> SS loadings        2.749  2.655  2.465  2.162  1.869  1.797
#> Prop Tot Var        .092   .088   .082   .072   .062   .060
#> Cum Prop Tot Var    .092   .180   .262   .334   .397   .457
#> Prop Comm Var       .201   .194   .180   .158   .136   .131
#> Cum Prop Comm Var   .201   .395   .574   .732   .869  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> χ²(270) = 618.06, p < .001
#> CFI: .93
#> TLI: .89
#> RMSEA [90% CI]: .05 [.05; .06]
#> AIC: 78.06
#> BIC: -1059.88
#> ECVI: 2.02
#> CAF: .48
#> SRMR: .03
```

Of course,
[`efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/efa_compare.md)
can also be used to compare results from different estimation or
rotation methods (in fact, to compare any two matrices), not just from
different implementations:

``` r

efa_compare(
  efa_fit(DOSPERT_sub, n_factors = 6, rotation = "promax")$rot_loadings,
  efa_fit(DOSPERT_sub, n_factors = 6, rotation = "oblimin", method = "ULS")$rot_loadings,
  x_labels = c("PAF and promax", "ULS and oblimin")
)
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
#> Mean [min, max] absolute difference:  0.0278 [ 0.0001,  0.1262]
#> Median absolute difference:  0.0226
#> Max decimals where all numbers are equal: 0
#> Minimum number of decimals provided: 18
#> 
#>           F1      F4      F3      F6      F2      F5
#> ethR_1   .0244  -.0335  -.0590  -.0073   .0157  -.0248
#> ethR_2   .0274  -.0299  -.0345   .0006   .0153  -.0348
#> ethR_3   .0208  -.0194  -.0522   .0025   .0073  -.0224
#> ethR_4   .0509   .0033  -.0253  -.0578   .0306  -.0099
#> ethR_5   .0110  -.0174  -.0482   .0045   .0028  -.0227
#> ethR_6   .0434   .0018  -.0325  -.0408   .0224  -.0190
#> finR_1  -.0027  -.0061   .0021   .0175  -.0001  -.0201
#> finR_2  -.0116   .0192  -.0214  -.0742   .0142   .0478
#> finR_3   .0072   .0077   .0114  -.0060   .0069  -.0179
#> finR_4  -.0240  -.0170  -.0625  -.0334   .0071   .0427
#> finR_5  -.0046  -.0054  -.0032   .0162  -.0006  -.0168
#> finR_6  -.0267  -.0230  -.0624  -.0316   .0075   .0378
#> heaR_1   .0156  -.0223  -.0331   .0003   .0095  -.0269
#> heaR_2   .0272  -.0276  -.0316  -.0019   .0148  -.0376
#> heaR_3  -.0213  -.0955  -.0820   .1154  -.0202  -.0595
#> heaR_4  -.0269  -.0912  -.0791   .1262  -.0266  -.0660
#> heaR_5   .0226  -.0393  -.0376   .0123   .0119  -.0404
#> heaR_6  -.0202  -.0618  -.0610   .0728  -.0150  -.0437
#> recR_1  -.0171  -.0288  -.0354   .0253  -.0074  -.0280
#> recR_2  -.0327  -.0634  -.0637   .1076  -.0301  -.0612
#> recR_3  -.0370  -.0393  -.0520   .0750  -.0277  -.0450
#> recR_4  -.0107   .0446  -.0201  -.0293  -.0118  -.0187
#> recR_5  -.0075   .0487  -.0177  -.0300  -.0117  -.0186
#> recR_6  -.0249  -.0134  -.0510   .0203  -.0149  -.0238
#> socR_1   .0259  -.0045   .0287  -.0333   .0225  -.0286
#> socR_2   .0137  -.0351  -.0068  -.0020   .0166  -.0377
#> socR_3   .0111   .0072   .0203  -.0524   .0206   .0005
#> socR_4   .0005  -.0483  -.0079   .0283   .0065  -.0417
#> socR_5   .0058  -.0272  -.0242  -.0084   .0114  -.0191
#> socR_6   .0199  -.0179  -.0076  -.0427   .0257  -.0124
```

Finally, if you are interested in factor scores from the EFA solution,
these can be obtained with
[`efa_scores()`](https://mdsteiner.github.io/EFAtools/reference/efa_scores.md),
which estimates the scores and their scoring weights directly from an
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
output, together with score-quality diagnostics (determinacy,
univocality, and the Guttman indeterminacy index):

``` r

efa_mod <- efa_fit(DOSPERT_sub, n_factors = 6, rotation = "promax")
#> ℹ `x` is not a correlation matrix; computing correlations from
#> the raw data.
fac_scores <- efa_scores(DOSPERT_sub, f = efa_mod)
fac_scores
#> 
#> ── Factor scores (regression) ──────────────────────────────────────────────────
#> 
#> Scored 500 observations on 6 factors.
#> 
#> ── Score determinacy ───────────────────────────────────────────────────────────
#> 
#>      rho  rho2 guttman
#> F1 0.903 0.815   0.629
#> F4 0.951 0.905   0.810
#> F3 0.961 0.923   0.846
#> F6 0.910 0.829   0.658
#> F2 0.887 0.787   0.574
#> F5 0.907 0.822   0.644
```

### Performance

To improve performance of the iterative procedures (currently the
parallel analysis, and the PAF, ML, and ULS methods) we implemented some
of them in C++. For example, the following code compares the EFAtools
parallel analysis with the corresponding one implemented in the psych
package (the default of
[`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
is to use 1000 datasets, but 25 is enough to show the difference):

``` r

microbenchmark::microbenchmark(
  efa_parallel(DOSPERT_sub, eigen_type = "SMC", n_datasets = 25),
  psych::fa.parallel(DOSPERT_sub, SMC = TRUE, plot = FALSE, n.iter = 25)
)
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  5 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4 
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Parallel analysis suggests that the number of factors =  10  and the number of components =  4
#> Unit: milliseconds
#>                                                                    expr
#>          efa_parallel(DOSPERT_sub, eigen_type = "SMC", n_datasets = 25)
#>  psych::fa.parallel(DOSPERT_sub, SMC = TRUE, plot = FALSE, n.iter = 25)
#>       min       lq     mean   median       uq      max neval
#>  101.7306 124.0934 133.1380 126.6526 142.9834 226.3888   100
#>  113.0272 119.6717 126.9424 125.2425 131.5727 180.6947   100
```

Moreover, the following code compares the PAF implementation (of type
“psych”) of the EFAtools package with the one from the psych package:

``` r

microbenchmark::microbenchmark(
  efa_fit(DOSPERT_raw, 6),
  psych::fa(DOSPERT_raw, 6, rotate = "none", fm = "pa")
)
#> Unit: milliseconds
#>                                                   expr      min       lq
#>                                efa_fit(DOSPERT_raw, 6) 15.29864 15.54146
#>  psych::fa(DOSPERT_raw, 6, rotate = "none", fm = "pa") 22.00119 22.44153
#>      mean   median       uq       max neval
#>  15.82454 15.68703 15.89799  19.89316   100
#>  25.15466 23.33760 26.12161 111.08977   100
```

While these differences are not large, they grow larger the more
iterations the procedures need, which is usually the case if solutions
are more tricky to find. Especially for simulations this might come in
handy. For example, in one simulation analysis we ran over 10,000,000
EFAs, thus a difference of about 25 milliseconds per EFA leads to a
difference in runtime of almost three days.

### Model Averaging

Instead of relying on one of the many possible implementations of, for
example, PAF, and of using just one rotation (e.g., promax), it may be
desirable to average different solutions to potentially arrive at a more
robust, average solution. The
[`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md)
function provides this possibility. In addition to the average solution
it provides the variation across solutions, a matrix indicating the
robustness of indicator-to-factor correspondences, and a visualisation
of the average solution and the variability across solutions. For
example, to average across all available factor extraction methods and
across all available oblique rotations, the following code can be run:

``` r

# Average solution across many different EFAs with oblique rotations
EFA_AV <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                      method = c("PAF", "ML", "ULS"), rotation = "oblique",
                      show_progress = FALSE)

# look at solution
EFA_AV
#> 
#> Averaging performed with averaging method mean (trim = 0) across 162 EFAs,
#> varying the following settings: method, init_comm, criterion_type,
#> start_method, rotation, k_promax, P_type, and varimax_type.
#> 
#> The error rate is at 0%. Of the solutions that did not result in an error, 100%
#> converged, 0% contained Heywood cases, and 100% were admissible.
#> 
#> ══ Indicator-to-Factor Correspondences ═════════════════════════════════════════
#> 
#> For each cell, the proportion of solutions including the respective
#> indicator-to-factor correspondence. A salience threshold of 0.3 was used to
#> determine indicator-to-factor correspondences.
#> 
#>       F1    F2    F3
#> V1    .00   .02  1.00
#> V2    .00   .02  1.00
#> V3    .00   .02  1.00
#> V4    .00   .02  1.00
#> V5    .00   .02  1.00
#> V6    .00   .02  1.00
#> V7    .00  1.00   .09
#> V8    .00  1.00   .09
#> V9    .00  1.00   .09
#> V10   .00  1.00   .09
#> V11   .00   .92   .09
#> V12   .00  1.00   .09
#> V13  1.00   .00   .06
#> V14   .94   .00   .06
#> V15   .99   .00   .06
#> V16   .94   .00   .06
#> V17  1.00   .00   .06
#> V18   .94   .00   .06
#> 
#> ══ Loadings ════════════════════════════════════════════════════════════════════
#> 
#> ── Mean ────────────────────────────────────────────────────────────────────────
#> 
#>        F1     F2     F3
#> V1   -.039   .062   .592
#> V2    .006   .090   .468
#> V3    .063   .079   .444
#> V4    .100   .023   .537
#> V5    .152   .012   .430
#> V6   -.063  -.014   .669
#> V7    .018   .510   .139
#> V8    .002   .550   .084
#> V9    .051   .522   .056
#> V10  -.004   .636  -.005
#> V11   .027   .349   .258
#> V12   .038   .620   .050
#> V13   .582   .098  -.010
#> V14   .515  -.043   .120
#> V15   .531   .136  -.013
#> V16   .523  -.027   .126
#> V17   .624  -.018   .023
#> V18   .524   .023   .091
#> 
#> ── Range ───────────────────────────────────────────────────────────────────────
#> 
#>       F1    F2    F3
#> V1   .280  .654  .131
#> V2   .230  .517  .126
#> V3   .242  .476  .127
#> V4   .304  .557  .123
#> V5   .278  .435  .166
#> V6   .303  .736  .173
#> V7   .092  .243  .495
#> V8   .087  .199  .536
#> V9   .074  .170  .511
#> V10  .099  .205  .637
#> V11  .147  .335  .353
#> V12  .093  .210  .615
#> V13  .343  .154  .669
#> V14  .341  .080  .521
#> V15  .317  .147  .646
#> V16  .345  .079  .529
#> V17  .389  .134  .635
#> V18  .336  .077  .558
#> 
#> ══ Factor Intercorrelations from Oblique Solutions ═════════════════════════════
#> 
#> ── Mean ────────────────────────────────────────────────────────────────────────
#> 
#>       F1     F2     F3
#> F1  1.000
#> F2   .512  1.000
#> F3   .543   .491  1.000
#> 
#> ── Range ───────────────────────────────────────────────────────────────────────
#> 
#>      F1    F2    F3
#> F1  .000
#> F2  .690  .000
#> F3  .702  .963  .000
#> 
#> ══ Variances Accounted for ═════════════════════════════════════════════════════
#> 
#> ── Mean ────────────────────────────────────────────────────────────────────────
#> 
#>                  F1     F2     F3
#> SS loadings    2.110  2.014  2.183
#> Prop Tot Var    .117   .112   .121
#> Prop Comm Var   .335   .319   .346
#> 
#> ── Range ───────────────────────────────────────────────────────────────────────
#> 
#>                  F1     F2     F3
#> SS loadings    1.652  2.461  3.762
#> Prop Tot Var    .092   .137   .209
#> Prop Comm Var   .262   .390   .596
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
```

The first matrix of the output tells us that the indicators are mostly
allocated to the same factors. However, that some rowsums are larger
than one also tells as that there likely are some cross loadings present
in some solutions. Moreover, the relatively high percentages of salient
pattern coefficients all loading on the first factor may indicate that
some rotation methods failed to achieve simple structure and it might be
desirable to exclude these from the model averaging procedure. The rest
of the output is similar to the normal
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
outputs shown above, only that in addition to the average coefficients
their range is also shown. Finally, the plot shows the average pattern
coefficients and their ranges.

**Important disclaimer:** While it is possible that this approach
provides more robust results, we are unaware of simulation studies that
have investigated and shown this. Therefore, it might make sense to for
now use this approach mainly to test the robustness of the results
obtained with one single EFA implementation.

## Exploratory Factor Analysis: Schmid-Leiman transformation and McDonald’s Omegas

For the Schmid-Leiman transformation and computation of omegas, we will
use PAF and promax rotation:

``` r

efa_dospert <- efa_fit(DOSPERT_sub, n_factors = 6, rotation = "promax")
#> ℹ `x` is not a correlation matrix; computing correlations from
#> the raw data.
efa_dospert
#> 
#> EFA performed with type = 'EFAtools', method = 'PAF', and rotation = 'promax'.
#> 
#> ── Rotated Loadings ────────────────────────────────────────────────────────────
#> 
#>           F1     F4     F3     F6     F2     F5    h2    u2
#> ethR_1   .547  -.025  -.020   .130   .009   .122  .442  .558
#> ethR_2   .447  -.055   .132   .098   .115  -.043  .316  .684
#> ethR_3   .532   .045   .068   .117  -.212   .030  .470  .530
#> ethR_4   .695  -.035  -.039  -.162   .000  -.002  .356  .644
#> ethR_5   .382   .159  -.051   .130  -.100   .036  .292  .708
#> ethR_6   .664   .086  -.003  -.089  -.063  -.044  .405  .595
#> finR_1  -.015  -.019   .853   .022   .003   .008  .730  .270
#> finR_2  -.050   .121  -.044  -.211   .120   .599  .359  .641
#> finR_3   .046   .035   .895  -.094   .061  -.011  .778  .222
#> finR_4   .039  -.071   .098   .035  -.102   .768  .633  .367
#> finR_5  -.006  -.010   .887   .023  -.038   .057  .842  .158
#> finR_6  -.022  -.051   .023   .061   .031   .774  .635  .365
#> heaR_1   .318   .042   .111   .091   .092   .019  .244  .756
#> heaR_2   .416   .020  -.003   .093   .191  -.081  .252  .748
#> heaR_3   .133  -.140  -.071   .674  -.036  -.009  .409  .591
#> heaR_4   .074  -.031  -.006   .697  -.067  -.074  .455  .545
#> heaR_5   .385  -.037  -.053   .170   .197  -.072  .240  .760
#> heaR_6   .044   .062  -.026   .455   .059   .048  .307  .693
#> recR_1  -.067   .253  -.060   .219   .247   .094  .322  .678
#> recR_2  -.046   .231   .089   .574  -.046  -.090  .493  .507
#> recR_3  -.155   .402   .053   .423   .035   .011  .507  .493
#> recR_4   .042   .926  -.010  -.089   .002  -.017  .774  .226
#> recR_5   .088   .945   .034  -.114  -.065  -.058  .784  .216
#> recR_6  -.020   .510  -.067   .212   .053   .154  .534  .466
#> socR_1   .016   .029   .025  -.137   .679  -.162  .400  .600
#> socR_2   .057  -.035   .103   .072   .602  -.020  .401  .599
#> socR_3  -.084   .061   .029  -.199   .559   .110  .348  .652
#> socR_4  -.090  -.082   .023   .197   .583  -.057  .372  .628
#> socR_5   .097   .032  -.086   .077   .370   .117  .234  .766
#> socR_6   .134  -.035  -.072  -.079   .568   .142  .363  .637
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── Factor Intercorrelations ────────────────────────────────────────────────────
#> 
#>       F1     F2     F3     F4     F5     F6
#> F1  1.000
#> F2   .357  1.000
#> F3   .455   .349  1.000
#> F4   .567   .594   .502  1.000
#> F5   .031   .316   .046   .287  1.000
#> F6   .329   .441   .463   .486   .324  1.000
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      F1     F4     F3     F6     F2     F5
#> SS loadings        2.565  2.518  2.457  2.225  2.213  1.717
#> Prop Tot Var        .086   .084   .082   .074   .074   .057
#> Cum Prop Tot Var    .086   .169   .251   .326   .399   .457
#> Prop Comm Var       .187   .184   .179   .162   .162   .125
#> Cum Prop Comm Var   .187   .371   .551   .713   .875  1.000
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> CAF: .48
#> SRMR: .03
#> df: 270
```

The indicator names in the output (i.e., the rownames of the rotated
loadings section) tell us which domain (out of ethical, financial,
health, recreational, and social risks) an indicator stems from. From
the pattern coefficients it can be seen that these theoretical domains
are recovered relatively well in the six factor solution, that is,
usually, the indicators from the same domain load onto the same factor.
When we take a look at the factor intercorrelations, we can see that
there are some strong and some weak correlations. It might be worthwhile
to explore whether a general factor can be obtained, and which factors
load more strongly on it. To this end, we will use a Schmid-Leiman (SL)
transformation.

## Schmid-Leiman Transformation

The SL transformation or orthogonalization transforms an oblique
solution into a hierarchical, orthogonalized solution. To do this, the
`EFAtools` package provides the
[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
function.

``` r

sl_dospert <- efa_schmid_leiman(efa_dospert)
sl_dospert
#> 
#> EFA for second-order loadings performed with type = 'EFAtools' and method = 'PAF'
#> 
#> ── Schmid-Leiman Solution ──────────────────────────────────────────────────────
#> 
#>           g     F1     F2     F3     F4     F5     F6    h2    u2
#> ethR_1  .487   .440   .009  -.016  -.018   .093   .067  .445  .555
#> ethR_2  .400   .360   .110   .104  -.041  -.033   .050  .318  .682
#> ethR_3  .441   .428  -.202   .054   .033   .023   .060  .427  .573
#> ethR_4  .225   .559   .000  -.031  -.026  -.001  -.084  .372  .628
#> ethR_5  .406   .308  -.095  -.041   .119   .028   .067  .290  .710
#> ethR_6  .325   .534  -.060  -.003   .064  -.034  -.046  .402  .598
#> finR_1  .523  -.012   .003   .676  -.014   .006   .011  .731  .269
#> finR_2  .266  -.041   .114  -.035   .090   .458  -.109  .316  .684
#> finR_3  .528   .037   .058   .709   .026  -.008  -.049  .789  .211
#> finR_4  .529   .032  -.097   .078  -.053   .587   .018  .643  .357
#> finR_5  .574  -.005  -.036   .703  -.008   .043   .012  .827  .173
#> finR_6  .528  -.018   .029   .018  -.038   .592   .031  .633  .367
#> heaR_1  .403   .256   .088   .088   .031   .014   .047  .246  .754
#> heaR_2  .344   .335   .182  -.002   .015  -.062   .048  .270  .730
#> heaR_3  .503   .107  -.034  -.056  -.104  -.007   .347  .400  .600
#> heaR_4  .548   .060  -.064  -.005  -.023  -.057   .359  .441  .559
#> heaR_5  .331   .310   .187  -.042  -.028  -.055   .088  .254  .746
#> heaR_6  .491   .036   .056  -.021   .046   .037   .235  .304  .696
#> recR_1  .417  -.054   .235  -.047   .189   .072   .113  .288  .712
#> recR_2  .601  -.037  -.044   .071   .172  -.069   .296  .492  .508
#> recR_3  .589  -.125   .033   .042   .299   .008   .218  .503  .497
#> recR_4  .550   .034   .002  -.008   .689  -.013  -.046  .780  .220
#> recR_5  .549   .071  -.062   .027   .704  -.044  -.059  .812  .188
#> recR_6  .585  -.016   .050  -.053   .380   .118   .109  .518  .482
#> socR_1  .030   .013   .646   .020   .021  -.124  -.071  .439  .561
#> socR_2  .306   .046   .573   .081  -.026  -.015   .037  .433  .567
#> socR_3  .081  -.068   .532   .023   .045   .084  -.103  .315  .685
#> socR_4  .217  -.073   .555   .018  -.061  -.044   .102  .376  .624
#> socR_5  .282   .078   .352  -.068   .023   .090   .040  .225  .775
#> socR_6  .211   .108   .540  -.057  -.026   .109  -.040  .366  .634
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                      g      F1     F2     F3     F4     F5     F6
#> SS loadings        5.710  1.550  1.994  1.518  1.327  1.000   .553
#> Prop Tot Var        .190   .052   .066   .051   .044   .033   .018
#> Cum Prop Tot Var    .190   .242   .308   .359   .403   .437   .455
#> Prop Comm Var       .418   .113   .146   .111   .097   .073   .040
#> Cum Prop Comm Var   .418   .532   .678   .789   .886   .960  1.000
```

From the output, it can be seen that all, except the social domain
indicators substantially load on the general factor. That is, the other
domains covary substantially.

## McDonald’s Omegas

Finally, we can compute omega estimates and additional indices of
interpretive relevance based on the SL solution, using
[`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md).
We can either specify the indicator-to-factor correspondences through
the `factor_map` argument, or let them be determined automatically (in
which case each indicator is assigned to its highest-loading factor,
which might lead to a different solution than what is desired, in the
presence of cross-loadings). Given that no cross-loadings are present
here, it is easiest to let the function determine the correspondences
automatically by leaving `factor_map` unspecified.

``` r

efa_reliability(sl_dospert)
#> 
#> Total variance from the correlation matrix.
#> 
#> ── Reliability coefficients ────────────────────────────────────────────────────
#> 
#>      tot  hier   sub  alpha    H
#> g   .922  .718  .177   .886  .886
#> F1  .791  .377  .415   .788  .652
#> F2  .739  .124  .615   .741  .726
#> F3  .912  .344  .568   .914  .739
#> F4  .867  .474  .393   .845  .683
#> F5  .745  .295  .451   .742  .571
#> F6  .709  .531  .177   .716  .305
#> 
#> ── Common-variance indices ─────────────────────────────────────────────────────
#> 
#>     ECV   PUC
#> g  .418  .828
```

If we wanted to specify the indicator-to-factor correspondences
explicitly (for example, according to theoretical expectations), we
could pass them as a `factor_map` matrix:

``` r

efa_reliability(sl_dospert, factor_map = matrix(c(rep(0, 18), rep(1, 6), rep(0, 30), 
                                         rep(1, 6), rep(0, 6), 1, 0, 1, 0, 1,
                                         rep(0, 19), rep(1, 6), rep(0, 31), 1, 0,
                                         1, 0, 1, rep(0, 30), rep(1, 6), 
                                         rep(0, 12)), ncol = 6, byrow = FALSE))
#> 
#> Total variance from the correlation matrix.
#> 
#> ── Reliability coefficients ────────────────────────────────────────────────────
#> 
#>      tot  hier   sub  alpha    H
#> g   .922  .718  .089   .886  .886
#> F1  .536  .535  .001   .844  .026
#> F2  .742  .082  .660   .735  .722
#> F3  .912  .344  .568   .914  .739
#> F4  .318  .317  .001   .763  .022
#> F5  .745  .295  .451   .742  .571
#> F6  .567  .479  .088   .697  .262
#> 
#> ── Common-variance indices ─────────────────────────────────────────────────────
#> 
#>     ECV   PUC
#> g  .418  .848
```
