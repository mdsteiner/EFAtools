# Estimate factor scores and score-quality diagnostics for an EFA model

Computes factor-score weights, and (from raw data) the factor scores
themselves, for an
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
solution or a directly supplied loading matrix. It also returns
score-quality diagnostics: the score intercorrelations, the determinacy
(validity) and univocality of each score, and Guttman's indeterminacy
index. Factor scores are returned only when raw data are supplied; a
correlation matrix yields the weights and diagnostics alone.

## Usage

``` r
efa_scores(
  x,
  f,
  Phi = NULL,
  rho = NULL,
  method = c("regression", "Bartlett", "Anderson", "tenBerge", "Harman", "components")
)
```

## Source

Thurstone, L. L. (1935). *The vectors of mind*. University of Chicago
Press.

Bartlett, M. S. (1937). The statistical conception of mental factors.
*British Journal of Psychology, 28*, 97-104.

Anderson, T. W., & Rubin, H. (1956). Statistical inference in factor
analysis. In *Proceedings of the Third Berkeley Symposium on
Mathematical Statistics and Probability* (Vol. 5, pp. 111-150).
University of California Press.

Guttman, L. (1955). The determinacy of factor score matrices with
implications for five other basic problems of common-factor theory.
*British Journal of Statistical Psychology, 8*, 65-81.

ten Berge, J. M. F., Krijnen, W. P., Wansbeek, T., & Shapiro, A. (1999).
Some new results on correlation-preserving factor scores prediction
methods. *Linear Algebra and its Applications, 289*, 311-318.

Grice, J. W. (2001). Computing and evaluating factor scores.
*Psychological Methods, 6*, 430-450.

## Arguments

- x:

  data.frame or matrix. Raw data (needed to obtain factor scores) or a
  correlation matrix (yields weights and diagnostics only). When `f` is
  a directly supplied loading matrix, a correlation-matrix `x` also
  supplies the correlations the weights are derived from; when `f` is an
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  object, its own fitted correlations are used instead (supply `rho` to
  derive the weights from another matrix). `x` describes the model
  variables either way. When raw data carry column names, they are
  matched to the model variables by name (any extra columns are ignored,
  and a model variable missing from `x` is an error). A named
  correlation matrix is likewise matched to the loading rows by name;
  its row and column names must use the same order, and it must carry
  one row and column for each model variable. Unnamed input is matched
  by position.

  Raw data are scored as supplied: no imputation is performed, so a case
  with a missing value on any model variable receives `NA` scores (and
  is reported as not scored). A model variable that carries no usable
  spread in `x` – constant, infinite, or observed fewer than twice – is
  an error.

- f:

  object of class
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md),
  an `efa_loadings` object, or a matrix of factor loadings.

- Phi:

  matrix. Factor intercorrelations. Only used when a loading matrix is
  supplied directly in `f`; taken from the `efa` object otherwise, in
  which case a supplied `Phi` is ignored with a warning. Named rows and
  columns are matched to the loading columns and must use the same
  order. Default is `NULL`, in which case the factors are assumed
  uncorrelated.

- rho:

  matrix. Correlation matrix used to derive the scoring weights.
  Defaults to `NULL`, in which case `f$orig_R` is used for an `efa`
  object; for a directly supplied loading matrix, `x` itself when it is
  a correlation matrix, and `cor(x, use = "pairwise")` otherwise. Pass a
  matrix here to score against a correlation other than the one implied
  by `f`/`x`. Named rows and columns are matched to the loading rows;
  row and column names must use the same order.

- method:

  character. The factor-score method: one of `"regression"` (default),
  `"Bartlett"`, `"Anderson"`, `"tenBerge"`, `"Harman"`, or
  `"components"`.

## Value

An object of class `efa_scores`, a list containing:

- weights:

  The `p` by `m` factor-score weight matrix.

- scores:

  The factor scores (`n` by `m`), or `NULL` when a correlation matrix
  was supplied. A case with a missing value on any model variable is not
  scored and keeps `NA` in every column.

- r.scores:

  The `m` by `m` correlations of the factor-score estimates (see Details
  for the `"components"`-method scale caveat).

- score_cor:

  The `m` by `m` score-factor correlation matrix; its diagonal is the
  determinacy (validity) of each score and its off-diagonals the
  univocality.

- determinacy:

  A data frame with, per factor, the determinacy `rho`, the squared
  determinacy `rho2`, and Guttman's indeterminacy index `guttman`.

- settings:

  A list of the settings used, including the number of supplied
  observations `n_obs` and the number of them that could be scored
  `n_scored`.

## Details

Each method combines the loadings with some or all of the factor
correlations and the scoring correlation matrix into weights in a
different way:

- `"regression"`:

  Thurstone's (1935) regression scores.

- `"Bartlett"`:

  Bartlett's (1937) scores.

- `"Anderson"`:

  Anderson & Rubin's (1956) scores.

- `"tenBerge"`:

  ten Berge, Krijnen, Wansbeek & Shapiro's (1999) scores.

- `"Harman"`:

  Harman's (1976) scores, based on an idealized variable (a hypothetical
  variable that would correlate perfectly with the factor).

- `"components"`:

  Component scores. These are formed from the raw, uncentered data
  (`X %*% W`) rather than the standardized data, so unlike the other
  methods they are on the scale of the input variables. The diagnostics
  below describe the standardized combination `scale(X) %*% W`, and
  therefore differ from the realized correlations of the returned scores
  whenever the variables have unequal variances.

The determinacy (validity) of a score is its correlation with the factor
it estimates, computed from the returned weights; for regression scores
it is the multiple correlation between the factor and the observed
variables (Guttman, 1955; Grice, 2001). The off-diagonal score-factor
correlations give the univocality (the correlation of a score with the
*other* factors). Guttman's (1955) indeterminacy index, `2 rho^2 - 1`,
is the minimum correlation between two equally valid sets of scores. For
a `method` other than `"regression"` both quantities are specific to
those scores: the determinacy is the method's own score-factor
correlation (never larger than the regression value), and the reported
`guttman` follows from it.

Determinacies close to 1 mean the scores stand in for the factor with
little loss; Grice (2001) regards values of about .90 and above as the
level required before scores are interpreted for individual cases, and
treats lower values as usable only for group-level research. The Guttman
index makes the same point more sharply, because a factor score is never
the factor: at `rho = .90` two equally valid sets of scores can still
correlate as low as .62, and at `rho = .80` as low as .28, so the rank
order of cases is not unique.

Which `method` to prefer follows from what the scores are for.
Regression scores correlate most highly with the factor, but they are
biased towards it and correlate across factors even when the model is
orthogonal. Bartlett scores are conditionally unbiased, which makes them
the choice when the scores stand in for the factor in a later model.
`"tenBerge"` reproduces the factor intercorrelations `Phi`, so it is the
choice when the scores will be correlated with each other or with
external variables. `"Anderson"` forces the scores to be uncorrelated
with unit variance and is appropriate only when the factors themselves
are orthogonal. `"components"` is a weighted sum of the observed
variables rather than an estimate of a common factor.

## See also

[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
for the solution these are computed from.

Other factor scoring:
[`print.efa_scores()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_scores.md)

## Examples

``` r
# Weights and score diagnostics from an EFA on a correlation matrix
efa <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
               estimator = "PAF", rotation = "oblimin")
fs <- efa_scores(test_models$baseline$cormat, f = efa)
#> ℹ `x` is a correlation matrix; factor scores cannot be computed. Only factor
#>   weights and score diagnostics are returned. Enter raw data to get factor
#>   scores.
fs
#> 
#> ── Factor scores (regression) ──────────────────────────────────────────────────
#> 
#> Weights and diagnostics only (correlation-matrix input; no scores).
#> 
#> ── Score determinacy ───────────────────────────────────────────────────────────
#> 
#>      rho  rho2  guttman
#> F1  .894  .798     .597
#> F2  .888  .788     .576
#> F3  .883  .780     .561
summary(fs)
#> 
#> ── Factor scores (regression) ──────────────────────────────────────────────────
#> 
#> Weights and diagnostics only (correlation-matrix input; no scores).
#> 
#> ── Score determinacy ───────────────────────────────────────────────────────────
#> 
#>      rho  rho2  guttman
#> F1  .894  .798     .597
#> F2  .888  .788     .576
#> F3  .883  .780     .561
#> 
#> ── Factor weights ──────────────────────────────────────────────────────────────
#> 
#>       F1    F2     F3
#> V1   .016  .037   .206
#> V2   .023  .036   .146
#> V3   .038  .033   .140
#> V4   .060  .025   .194
#> V5   .062  .014   .138
#> V6   .009  .013   .247
#> V7   .024  .177   .053
#> V8   .017  .187   .031
#> V9   .031  .173   .020
#> V10  .016  .222  -.002
#> V11  .026  .115   .084
#> V12  .035  .240   .030
#> V13  .202  .051   .010
#> V14  .163  .002   .050
#> V15  .177  .059   .007
#> V16  .170  .006   .050
#> V17  .214  .013   .016
#> V18  .173  .025   .041
#> 
#> ── Score validity and univocality ──────────────────────────────────────────────
#> 
#> Diagonal: validity (score-factor correlation). Off-diagonal: univocality.
#> 
#>      F1    F2    F3
#> F1  .894  .638  .668
#> F2  .643  .888  .650
#> F3  .676  .653  .883
#> 
#> ── Score intercorrelations ─────────────────────────────────────────────────────
#> 
#>       F1     F2     F3
#> F1  1.000   .719   .756
#> F2   .719  1.000   .735
#> F3   .756   .735  1.000

# Factor scores from raw data (Bartlett method)
# \donttest{
efa_raw <- efa_fit(GRiPS_raw, n_factors = 1, estimator = "PAF")
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
efa_scores(GRiPS_raw, f = efa_raw, method = "Bartlett")
#> 
#> ── Factor scores (Bartlett) ────────────────────────────────────────────────────
#> 
#> Scored 810 of 810 observations on 1 factor (see `$scores`).
#> 
#> ── Score determinacy ───────────────────────────────────────────────────────────
#> 
#>      rho  rho2  guttman
#> F1  .972  .946     .891
# }

# Loadings supplied directly, with the factor intercorrelations
efa_scores(test_models$baseline$cormat, f = efa$rot_loadings, Phi = efa$Phi)
#> ℹ `x` is a correlation matrix; factor scores cannot be computed. Only factor
#>   weights and score diagnostics are returned. Enter raw data to get factor
#>   scores.
#> 
#> ── Factor scores (regression) ──────────────────────────────────────────────────
#> 
#> Weights and diagnostics only (correlation-matrix input; no scores).
#> 
#> ── Score determinacy ───────────────────────────────────────────────────────────
#> 
#>      rho  rho2  guttman
#> F1  .894  .798     .597
#> F2  .888  .788     .576
#> F3  .883  .780     .561
```
