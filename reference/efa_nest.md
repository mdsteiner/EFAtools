# Next eigenvalue sufficiency test (NEST)

NEST uses many synthetic datasets to generate reference eigenvalues
against which to compare the empirical eigenvalues. This is similar to
parallel analysis, but other than parallel analysis, NEST does not just
rely on synthetic eigenvalues based on an identity matrix as null model.
It was introduced by Achim (2017), see also Brandenburg and Papenberg
(2024) and Caron (2025) for further simulation studies including NEST.

## Usage

``` r
efa_nest(
  x,
  N = NA,
  alpha = 0.05,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
  n_datasets = 1000,
  estimate_control = NULL,
  ...
)
```

## Source

Achim, A. (2017). Testing the number of required dimensions in
exploratory factor analysis. The Quantitative Methods for Psychology,
13(1), 64–74. https://doi.org/10.20982/tqmp.13.1.p064

Brandenburg, N., & Papenberg, M. (2024). Reassessment of innovative
methods to determine the number of factors: A simulation-based
comparison of exploratory graph analysis and Next Eigenvalue Sufficiency
Test. Psychological Methods, 29(1), 21–47.
https://doi.org/10.1037/met0000527

Caron, P.-O. (2025). A Comparison of the Next Eigenvalue Sufficiency
Test to Other Stopping Rules for the Number of Factors in Factor
Analysis. Educational and Psychological Measurement, Online-first
publication. https://doi.org/10.1177/00131644241308528

## Arguments

- x:

  data.frame or matrix. Dataframe or matrix of raw data or matrix with
  correlations.

- N:

  numeric. The number of observations. Only needed if x is a correlation
  matrix.

- alpha:

  numeric. The alpha level to use (i.e., 1-alpha percentile of
  eigenvalues is used for reference values).

- use:

  character. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) if raw data is
  given as input. Default is "pairwise.complete.obs".

- cor_method:

  character. One of `"pearson"`, `"spearman"`, or `"kendall"`, passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html). `"poly"` and
  `"tetra"` are not supported because `NEST` compares the data against
  simulated continuous reference data. Default is `"pearson"`.

- n_datasets:

  numeric. The number of datasets to simulate. Default is 1000.

- estimate_control:

  an
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  object with the estimation settings for the
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  reference-model fits. `NULL` (default) uses the
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  defaults. The reference models are unrotated, so no rotation settings
  apply.

- ...:

  Additional arguments passed to
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md).
  For example, `estimator`, to change the estimator (default is "PAF").
  PAF is more robust, but it will take longer compared to the other
  estimators available ("ML" and "ULS"). The estimation tuning knobs are
  not passed here; they live in `estimate_control`, and the
  standard-error arguments (`se`, `b_boot`, `ci`, `seed`) are not
  accepted because the reference-model fits are internal steps.

## Value

An object of class `efa_retention` (see
[`print.efa_retention()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_retention.md)
and
[`plot.efa_retention()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_retention.md)
for the print and plot methods). Its main fields are:

- n_factors:

  A named numeric vector (`"NEST"`) with the suggested number of factors
  according to the NEST procedure.

- results:

  A list with a single record holding the empirical eigenvalues and the
  reference eigenvalues. Only the positions the search actually tested
  carry a reference value; beyond the position at which it stopped the
  series is `NA`.

- settings:

  A list of control settings used.

## Details

NEST compares the first empirical eigenvalue against the first
eigenvalues of `n_dataset` synthetic datasets based on a null model
(i.e., with uncorrelated variables; same as in parallel analysis, see
[`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)).
The following eigenvalues are compared against synthetic datasets based
on an EFA-model with one fewer factors than the position of the
respective empirical eigenvalue. E.g, the second empirical eigenvalue is
compared against synthetic data based on a one-factor model. In each
comparison the \\k\\-th empirical eigenvalue is tested against the
\\k\\-th largest eigenvalue of the synthetic datasets. The `alpha`-level
defines against which percentile of the synthetic eigenvalue
distribution to compare the empirical eigenvalues against, i.e., an
alpha of .05 (the default) uses the 95th percentile as reference value.

The number of factors tested is capped at \\\lfloor 0.8 \times p
\rfloor\\ (with \\p\\ the number of variables; Achim, 2017) and
additionally limited so that the \\(k - 1)\\-factor reference model used
at each step stays over-identified. If no empirical eigenvalue falls at
or below its reference within this range, every tested factor is
accepted and this capped number is returned.

Because each reference model carries the factors already retained, NEST
does not lose accuracy for the strongly correlated factor structures
where parallel analysis tends to under-extract, and it was among the
more accurate criteria in the simulation studies of Brandenburg and
Papenberg (2024) and Caron (2025). The price is runtime: a fresh set of
`n_datasets` reference datasets is drawn and eigen-decomposed at every
candidate factor count, which makes NEST one of the slowest criteria
available here.

The reference models are fitted without inequality constraints. A
Heywood case in one of them leaves no unique variance to simulate the
reference data from, so NEST aborts rather than continuing from an
inadmissible reference.

The reference eigenvalues are obtained from simulated data, so the
suggested number of factors varies slightly from run to run. Call
[`set.seed()`](https://rdrr.io/r/base/Random.html) beforehand to make a
run reproducible.

For details on the method, including simulation studies, see Achim
(2017), Brandenburg and Papenberg (2024), and Caron (2025).

The `efa_nest` function can also be called together with other factor
retention criteria in the
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
[`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md),
[`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
[`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md)

## Examples

``` r
# \donttest{
# with correlation matrix
efa_nest(test_models$baseline$cormat, N = 500)
#> ── Next Eigenvalue Sufficiency Test ────────────────────────────────────────────
#> 
#> • Suggested number of factors: 3

# with raw data
efa_nest(GRiPS_raw)
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
#> ── Next Eigenvalue Sufficiency Test ────────────────────────────────────────────
#> 
#> • Suggested number of factors: 1
# }
```
