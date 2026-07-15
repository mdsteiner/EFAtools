# Estimate factor scores and score-quality diagnostics for an EFA model

Computes factor-score weights (and, from raw data, the factor scores
themselves) natively for an
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
solution or a directly supplied loading matrix, together with the
score-quality diagnostics that describe how well the estimated scores
represent the factors: the score intercorrelations, the determinacy
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
  correlation matrix (yields weights and diagnostics only). When raw
  data carry column names, they are matched to the model variables by
  name (any extra columns are ignored, and a model variable missing from
  `x` is an error); unnamed data are matched by position.

- f:

  object of class
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md),
  an `efa_loadings` object, or a matrix of factor loadings.

- Phi:

  matrix. Factor intercorrelations. Only used when a loading matrix is
  supplied directly in `f`; taken from the `efa` object otherwise.
  Default is `NULL`, in which case the factors are assumed uncorrelated.

- rho:

  matrix. Correlation matrix used to derive the scoring weights.
  Defaults to `NULL`, in which case `f$orig_R` is used for an `efa`
  object and `cor(x, use = "pairwise")` otherwise. Pass a matrix here to
  score against a correlation other than the one implied by `f`/`x`.

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
  was supplied.

- r.scores:

  The `m` by `m` correlations of the factor-score estimates.

- score_cor:

  The `m` by `m` score-factor correlation matrix; its diagonal is the
  determinacy (validity) of each score and its off-diagonals the
  univocality.

- determinacy:

  A data frame with, per factor, the determinacy `rho`, the squared
  determinacy `rho2`, and Guttman's indeterminacy index `guttman`.

- settings:

  A list of the settings used.

## Details

The `p` by `m` weight matrix `W` (standardized scores are
`scale(X) %*% W`) is computed from the structure matrix
`S = Lambda %*% Phi`, the model uniquenesses `Psi = diag(1 - h2)`, and
the scoring correlation matrix `R` according to `method`:

- `"regression"`:

  Thurstone's (1935) regression scores, `W = R^-1 S`.

- `"Bartlett"`:

  Bartlett's (1937) conditionally unbiased scores.

- `"Anderson"`:

  Anderson & Rubin's (1956) uncorrelated, unit-variance scores; defined
  for orthogonal factors only.

- `"tenBerge"`:

  ten Berge, Krijnen, Wansbeek & Shapiro's (1999) scores, which preserve
  the factor intercorrelations.

- `"Harman"`:

  Harman's (1976) idealized-variable scores.

- `"components"`:

  component scores, `W = Lambda`.

The determinacy (validity) of a score is its correlation with the factor
it estimates, computed from the returned weights; for regression scores
it is the multiple correlation between the factor and the observed
variables (Guttman, 1955; Grice, 2001). The off-diagonal score-factor
correlations give the univocality (the correlation of a score with the
*other* factors), and `2 rho^2 - 1` is Guttman's (1955) indeterminacy
index, the minimum correlation between two equally valid sets of scores.
For a `method` other than `"regression"` both quantities are specific to
those scores: the determinacy is the method's own score-factor
correlation (never larger than the regression value), and the reported
`guttman` follows from it.

## See also

Other factor scoring:
[`print.efa_scores()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_scores.md)

## Examples

``` r
# Weights and score diagnostics from an EFA on a correlation matrix
efa <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
               method = "PAF", rotation = "oblimin")
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
#>      rho  rho2 guttman
#> F1 0.894 0.798   0.597
#> F2 0.888 0.788   0.576
#> F3 0.883 0.780   0.561
summary(fs)
#> 
#> ── Factor scores (regression) ──────────────────────────────────────────────────
#> 
#> Weights and diagnostics only (correlation-matrix input; no scores).
#> 
#> ── Score determinacy ───────────────────────────────────────────────────────────
#> 
#>      rho  rho2 guttman
#> F1 0.894 0.798   0.597
#> F2 0.888 0.788   0.576
#> F3 0.883 0.780   0.561
#> 
#> ── Factor weights ──────────────────────────────────────────────────────────────
#> 
#>        F1    F2     F3
#> V1  0.016 0.037  0.206
#> V2  0.023 0.036  0.146
#> V3  0.038 0.033  0.140
#> V4  0.060 0.025  0.194
#> V5  0.062 0.014  0.138
#> V6  0.009 0.013  0.247
#> V7  0.024 0.177  0.053
#> V8  0.017 0.187  0.031
#> V9  0.031 0.173  0.020
#> V10 0.016 0.222 -0.002
#> V11 0.026 0.115  0.084
#> V12 0.035 0.240  0.030
#> V13 0.202 0.051  0.010
#> V14 0.163 0.002  0.050
#> V15 0.177 0.059  0.007
#> V16 0.170 0.006  0.050
#> V17 0.214 0.013  0.016
#> V18 0.173 0.025  0.041
#> 
#> ── Score validity and univocality ──────────────────────────────────────────────
#> 
#> Diagonal: validity (score-factor correlation). Off-diagonal: univocality.
#> 
#>       F1    F2    F3
#> F1 0.894 0.638 0.668
#> F2 0.643 0.888 0.650
#> F3 0.676 0.653 0.883
#> 
#> ── Score intercorrelations ─────────────────────────────────────────────────────
#> 
#>       F1    F2    F3
#> F1 1.000 0.719 0.756
#> F2 0.719 1.000 0.735
#> F3 0.756 0.735 1.000

# Factor scores from raw data (Bartlett method)
# \donttest{
efa_raw <- efa_fit(GRiPS_raw, n_factors = 1, method = "PAF")
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
efa_scores(GRiPS_raw, f = efa_raw, method = "Bartlett")
#> 
#> ── Factor scores (Bartlett) ────────────────────────────────────────────────────
#> 
#> Scored 810 observations on 1 factor.
#> 
#> ── Score determinacy ───────────────────────────────────────────────────────────
#> 
#>      rho  rho2 guttman
#> F1 0.972 0.946   0.891
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
#>      rho  rho2 guttman
#> F1 0.894 0.798   0.597
#> F2 0.888 0.788   0.576
#> F3 0.883 0.780   0.561
```
