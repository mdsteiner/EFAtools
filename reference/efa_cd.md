# Comparison data

Factor retention method introduced by Ruscio and Roche (2012). The code
was adapted from the CD code published by Auerswald and Moshagen (2019),
available at <https://osf.io/x5cz2/>.

## Usage

``` r
efa_cd(
  x,
  n_factors_max = NA,
  N_pop = 10000,
  N_samples = 500,
  alpha = 0.3,
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
  max_iter = 50
)
```

## Source

Auerswald, M., & Moshagen, M. (2019). How to determine the number of
factors to retain in exploratory factor analysis: A comparison of
extraction methods under realistic conditions. Psychological Methods,
24(4), 468–491. https://doi.org/10.1037/met0000200

Ruscio, J., & Roche, B. (2012). Determining the number of factors to
retain in an exploratory factor analysis using comparison data of known
factorial structure. Psychological Assessment, 24, 282–292. doi:
10.1037/a0025697

## Arguments

- x:

  data.frame or matrix. Dataframe or matrix of raw data.

- n_factors_max:

  numeric. The maximum number of factors to test against. Larger numbers
  will increase the duration the procedure takes, but test more possible
  solutions. If left NA (default) the maximum number of factors for
  which the model is still over-identified (df \> 0) is used.

- N_pop:

  numeric. Size of finite populations of comparison data. Default is
  10000.

- N_samples:

  numeric. Number of samples drawn from each population. Default is 500.

- alpha:

  numeric. The alpha level used to test the significance of the
  improvement added by an additional factor. Default is .30.

- cor_method:

  character. One of `"pearson"`, `"spearman"`, or `"kendall"`, passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html). `"poly"` and
  `"tetra"` are not supported because `CD` compares the data against
  simulated continuous reference data. Default is "pearson".

- max_iter:

  numeric. The maximum number of iterations after which the iterative
  PAF procedure inside the comparison-data generation is halted; it does
  not cap an EFA of `x`. Default is 50.

## Value

An object of class `efa_retention` (see
[`print.efa_retention()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_retention.md)
and
[`plot.efa_retention()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_retention.md)
for the print and plot methods). Its main fields are:

- n_factors:

  A named numeric vector (`"CD"`) with the suggested number of factors
  according to comparison data results.

- results:

  A list with a single record holding the mean RMSE between the
  eigenvalues of the generated and the entered data per number of
  factors (used for the plot) and, in `rmse_eigenvalues`, the per-sample
  RMSE matrix (rows are samples, columns are factor counts; columns
  beyond the last tested factor count are left as zero).

- settings:

  A list of the settings used.

## Details

Comparison data (CD) extends parallel analysis by reproducing the
observed correlation matrix rather than generating random data: datasets
with a known factor structure are generated with an increasing number of
factors, and the smallest number for which adding a further factor no
longer significantly improves the reproduction of the observed
eigenvalues is retained (Ruscio & Roche, 2012).

Because it reproduces the observed correlation structure instead of a
null model, CD was among the more accurate criteria across a broad range
of conditions in Ruscio and Roche (2012). It is, however, the only
criterion in this family that requires raw data, and by some margin the
most computationally intensive one: a finite population of `N_pop` cases
is generated and `N_samples` samples are drawn from it at every
candidate factor count. It is therefore a good choice when the raw data
are at hand and the runtime is acceptable, and a poor one for a quick
look at a correlation matrix.

The comparison data are obtained by simulation, so the suggested number
of factors varies slightly from run to run. Call
[`base::set.seed()`](https://rdrr.io/r/base/Random.html) beforehand to
make a run reproducible.

Note that if the data contains missing values, these will be removed for
the comparison data procedure using
[`stats::na.omit()`](https://rdrr.io/r/stats/na.fail.html). If missing
data should be treated differently, e.g., by imputation, do this outside
`efa_cd()` and then pass the complete data.

## See also

[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
as a wrapper function for this and the other factor retention criteria.

Other factor retention criteria:
[`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md),
[`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md),
[`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
[`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md),
[`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md),
[`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md),
[`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
[`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md)

## Examples

``` r
# \donttest{
# determine n factors of the GRiPS
efa_cd(GRiPS_raw, N_pop = 500, N_samples = 20)
#> ── Comparison data ─────────────────────────────────────────────────────────────
#> 
#> • Suggested number of factors: 1

# determine n factors of the DOSPERT risk subscale
efa_cd(DOSPERT_raw, N_pop = 500, N_samples = 20)
#> ── Comparison data ─────────────────────────────────────────────────────────────
#> 
#> • Suggested number of factors: 5
# }
```
