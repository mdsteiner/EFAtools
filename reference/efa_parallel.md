# Parallel analysis

Various methods for performing parallel analysis. This function uses
[future_lapply()](https://future.apply.futureverse.org/reference/future_lapply.html)
for which a parallel processing plan can be selected. To do so, register
a plan with
[`future::plan()`](https://future.futureverse.org/reference/plan.html),
for example `future::plan(future::multisession, workers = 2)`; see
examples.

## Usage

``` r
efa_parallel(
  x = NULL,
  N = NA,
  n_vars = NA,
  n_datasets = 1000,
  percent = 95,
  eigen_type = c("PCA", "SMC", "EFA"),
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
  decision_rule = c("means", "percentile", "crawford"),
  n_factors = 1,
  estimate_control = NULL,
  ...
)
```

## Source

Braeken, J., & van Assen, M. A. (2017). An empirical Kaiser criterion.
Psychological Methods, 22, 450–466. https://doi.org/10.1037/met0000074

Crawford, A. V., Green, S. B., Levy, R., Lo, W. J., Scott, L., Svetina,
D., & Thompson, M. S. (2010). Evaluation of parallel analysis methods
for determining the number of factors. Educational and Psychological
Measurement, 70(6), 885-901.

Glorfeld, L. W. (1995). An improvement on Horn's parallel analysis
methodology for selecting the correct number of factors to retain.
Educational and Psychological Measurement, 55(3), 377-393.

Horn, J. L. (1965). A rationale and test for the number of factors in
factor analysis. Psychometrika, 30(2), 179–185.
https://doi.org/10.1007/BF02289447

## Arguments

- x:

  matrix or data.frame. The real data to compare the simulated
  eigenvalues against. Must not contain variables of classes other than
  numeric. Can be a correlation matrix or raw data.

- N:

  numeric. The number of cases / observations to simulate. Only has to
  be specified if `x` is either a correlation matrix or `NULL`. If x
  contains raw data, `N` is found from the dimensions of `x`. Must be
  larger than the number of variables.

- n_vars:

  numeric. The number of variables / indicators to simulate. Only has to
  be specified if `x` is left as `NULL` as otherwise the dimensions are
  taken from `x`.

- n_datasets:

  numeric. The number of datasets to simulate. Must be at least 1.
  Default is 1000.

- percent:

  numeric. The percentile to take from the simulated eigenvalues.
  Default is 95.

- eigen_type:

  character. On what the eigenvalues should be found. Can be either
  "SMC", "PCA", or "EFA". If using "SMC", the diagonal of the
  correlation matrix is replaced by the squared multiple correlations
  (SMCs) of the indicators. If using "PCA", the diagonal values of the
  correlation matrices are left to be 1. If using "EFA", eigenvalues are
  found on the correlation matrices with the final communalities of an
  EFA solution as diagonal. Default is `c("PCA", "SMC", "EFA")`, i.e.
  all three, which costs roughly six times a single non-EFA type:
  `"EFA"` fits an EFA to every simulated dataset and dominates that
  total. Pass a single type if the run is time-critical.

- use:

  character. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) if raw data is
  given as input. Default is "pairwise.complete.obs".

- cor_method:

  character. One of `"pearson"`, `"spearman"`, or `"kendall"`, passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html). `"poly"` and
  `"tetra"` are not supported because `PARALLEL` compares the data
  against simulated continuous reference data. Default is "pearson".

- decision_rule:

  character. Which rule to use to determine the number of factors to
  retain. Default is `"means"`, which will use the average simulated
  eigenvalues. `"percentile"`, uses the percentiles specified in
  percent. `"crawford"` uses the 95th percentile for the first factor
  and the mean afterwards (based on Crawford et al, 2010). All three
  rules retain the factors up to the first observed eigenvalue that
  fails to exceed its reference value; an eigenvalue further down the
  series that rises above its own reference again therefore adds no
  factor. Because the average simulated eigenvalue is a lower reference
  than the percentile, `"means"` tends to retain more factors than the
  more conservative `"percentile"` rule (Glorfeld, 1995).

- n_factors:

  numeric. Number of factors to extract if "EFA" is included in
  `eigen_type`. Default is 1.

- estimate_control:

  an
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  object with the estimation settings for the
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  fits (of both the real and the simulated data) when `"EFA"` is
  included in `eigen_type`. `NULL` (default) uses the
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  defaults. The fits are unrotated, so no rotation settings apply.

- ...:

  Additional arguments passed to
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md).
  For example, `estimator`, to change the estimator (default is "PAF").
  PAF is more robust, but it will take longer compared to the other
  estimators available ("ML" and "ULS"). The estimation tuning knobs are
  not passed here; they live in `estimate_control`, and the
  standard-error arguments (`se`, `b_boot`, `ci`, `seed`) are not
  accepted because the fits are internal steps that keep only their
  eigenvalues.

## Value

An object of class `efa_retention` (see
[`print.efa_retention()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_retention.md)
and
[`plot.efa_retention()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_retention.md)
for the print and plot methods). Its main fields are:

- n_factors:

  A named numeric vector with the suggested number of factors for each
  requested eigenvalue type (`"PCA"`, `"SMC"`, and/or `"EFA"`). These
  are `NA` when no real data are supplied (i.e. only `N` and `n_vars`
  are given). When every observed eigenvalue exceeds its reference value
  (no crossing is found), all `n_vars` components are retained and a
  warning is issued.

- results:

  A list with one record per eigenvalue type, each holding the observed
  eigenvalues (when real data were supplied) and the simulated reference
  values (means and percentiles) used for printing and plotting.

- settings:

  A list of the settings used.

## Details

Parallel analysis (Horn, 1965) compares the eigenvalues obtained from
the sample correlation matrix against those of null model correlation
matrices (i.e., with uncorrelated variables) of the same sample size.
This way, it accounts for the variation in eigenvalues introduced by
sampling error and thus eliminates the main problem inherent in the
Kaiser-Guttman criterion
([`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md)).

Parallel analysis is often argued to be one of the most accurate factor
retention criteria. However, for highly correlated factor structures it
has been shown to underestimate the correct number of factors. The
reason for this is that a null model (uncorrelated variables) is used as
reference. However, when factors are highly correlated, the first
eigenvalue will be much larger compared to the following ones, as later
eigenvalues are conditional on the earlier ones in the sequence and thus
the shared variance is already accounted in the first eigenvalue (e.g.,
Braeken & van Assen, 2017).

The reference eigenvalues are obtained from simulated data, so the
suggested number of factors varies slightly from run to run. Call
[`base::set.seed()`](https://rdrr.io/r/base/Random.html) beforehand to
make a run reproducible; the result is then also independent of the
parallel plan set via
[`future::plan()`](https://future.futureverse.org/reference/plan.html),
so it can be reproduced on a machine with a different number of cores.
For `"PCA"` and `"SMC"` the simulation is drawn in independently seeded
blocks; a block that fails – which happens when a simulated correlation
matrix is singular, so that no eigenvalues can be taken from it – is
redrawn on its own, leaving the blocks that succeeded with the draws
they already made. The `"EFA"` series instead redraws the single dataset
that could not be fitted; if that dataset still cannot be fitted, the
call stops with an error.

When both `"PCA"` and `"SMC"` are requested, the two are read off the
*same* simulated datasets rather than from two independent simulations:
they differ only in the diagonal substituted into the simulated
correlation matrix, so one set of draws serves both and the two
reference series are paired dataset by dataset. A draw that cannot be
used for the SMC series – a simulated matrix with no inverse, and hence
no squared multiple correlations – is discarded for the `"PCA"` series
as well, so that the pairing stays exact. `"EFA"` fits a model to each
simulated dataset and draws its own.

The `efa_parallel` function can also be called together with other
factor retention criteria in the
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
function.

## See also

[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
as a wrapper function for this and the other factor retention criteria.

Other factor retention criteria:
[`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md),
[`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md),
[`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md),
[`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
[`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md),
[`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md),
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md),
[`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
[`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md)

## Examples

``` r
# \donttest{
# example without real data
pa_unreal <- efa_parallel(N = 500, n_vars = 10, n_datasets = 100)

# example with correlation matrix with all eigen_types and PAF estimation
pa_paf <- efa_parallel(test_models$case_11b$cormat, N = 500, n_datasets = 100)

# example with correlation matrix with all eigen_types and ML estimation
# this will be faster than the above with PAF)
pa_ml <- efa_parallel(test_models$case_11b$cormat, N = 500, estimator = "ML",
                      n_datasets = 100)
# }

if (FALSE) { # \dontrun{
# for parallel computation. future::plan() returns the plan it replaces, so
# on.exit() puts the session back as it was -- also if the call fails.
pa_faster <- local({
  old_plan <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old_plan), add = TRUE)
  efa_parallel(test_models$case_11b$cormat, N = 500)
})
} # }
```
