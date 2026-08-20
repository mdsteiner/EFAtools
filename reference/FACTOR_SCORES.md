# Estimate factor scores for an EFA model

**\[superseded\]**

`FACTOR_SCORES()` has been superseded by
[`efa_scores()`](https://mdsteiner.github.io/EFAtools/reference/efa_scores.md),
which is the recommended interface going forward. It remains available
so existing code keeps working. Note that `R2` is now the squared
factor-score determinacy of the *requested* `method`: the squared
correlation between a factor and the scores that method produces. For
`method = "Thurstone"` this is each factor's squared multiple
correlation with the observed variables (the value
[`psych::factor.scores()`](https://rdrr.io/pkg/psych/man/factor.scores.html)
returns with `Grice = TRUE`); for every other method it is smaller,
because no estimator correlates more highly with the factor than the
regression estimator does. Earlier versions returned psych's default
`Grice = FALSE` validity coefficient, so the slot is not comparable
across versions.

A convenience wrapper around
[`efa_scores()`](https://mdsteiner.github.io/EFAtools/reference/efa_scores.md)
that returns factor scores and weights in a compact list. Factor scores
are calculated according to the specified method if raw data are
provided, and only factor weights if a correlation matrix is provided.

## Usage

``` r
FACTOR_SCORES(
  x,
  f,
  Phi = NULL,
  rho = NULL,
  method = c("Thurstone", "tenBerge", "Anderson", "Bartlett", "Harman", "components")
)
```

## Arguments

- x:

  data.frame or matrix. Dataframe or matrix of raw data (needed to get
  factor scores) or matrix with correlations.

- f:

  object of class
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  or matrix.

- Phi:

  matrix. A matrix of factor intercorrelations. Only needs to be
  specified if a factor loadings matrix is entered directly into `f`;
  for an
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  object the intercorrelations are taken from the object, and a supplied
  `Phi` is ignored with a warning. Default is `NULL`, in which case the
  intercorrelations of a directly supplied loading matrix are assumed to
  be zero.

- rho:

  matrix. Correlation matrix used to derive the scoring weights.
  Defaults to `NULL`, in which case the matrix the EFA in `f` was fit on
  (`f$orig_R`) is used, so the weights stay consistent with the loadings
  even for a non-Pearson correlation (e.g. polychoric); for a directly
  supplied loading matrix, `x` itself is used when it is a correlation
  matrix, otherwise the Pearson correlation of `x`. Pass a matrix here
  to score against a different correlation.

- method:

  character. The method used to calculate factor scores. One of
  "Thurstone" (regression-based; default), "tenBerge", "Anderson",
  "Bartlett", "Harman", or "components".

## Value

A list of class FACTOR_SCORES containing the following:

- scores:

  The factor scores (only if raw data are provided.)

- weights:

  The factor weights.

- r.scores:

  The correlations of the factor score estimates.

- missing:

  Whether the raw data contained missing values (only if raw data are
  provided).

- R2:

  The squared factor-score determinacy for each factor: the squared
  correlation between a factor and the score the requested `method`
  produces. For `method = "Thurstone"` this equals the squared multiple
  correlation between the factor and the observed variables; for every
  other method it is specific to those scores and smaller. See
  [`efa_scores()`](https://mdsteiner.github.io/EFAtools/reference/efa_scores.md)
  for the underlying score-quality diagnostics.

- settings:

  A list of the settings used.

## See also

[`efa_scores()`](https://mdsteiner.github.io/EFAtools/reference/efa_scores.md)
for the factor-score weights together with the full set of score-quality
diagnostics (determinacy, univocality, and Guttman indeterminacy index)
and a print/summary method.

## Examples

``` r
# Example with raw data with method "Bartlett"
EFA_raw <- efa_fit(DOSPERT_raw, n_factors = 10, estimator = "PAF",
                   rotation = "oblimin",
                   rotate_control = rotate_control(random_starts = 1))
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
fac_scores_raw <- FACTOR_SCORES(DOSPERT_raw, f = EFA_raw, method = "Bartlett")

# Same as above, but with raw data AND a correlation matrix
cor_pearson <- cor(DOSPERT_raw)
EFA_cor_pearson <- efa_fit(cor_pearson, n_factors = 10, N = nrow(DOSPERT_raw),
                           estimator = "PAF", rotation = "oblimin",
                           rotate_control = rotate_control(random_starts = 1))
fac_scores_cor_pearson <- FACTOR_SCORES(DOSPERT_raw, f = EFA_cor_pearson,
                                        rho = cor_pearson,
                                        method = "Bartlett")

# Scores between two alternatives above are identical
isTRUE(all.equal(fac_scores_raw$scores, fac_scores_cor_pearson$scores,
                 check.attributes = FALSE))
#> [1] TRUE

# Example with a correlation matrix only (does not return factor scores)
EFA_cor <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
                   estimator = "PAF", rotation = "oblimin")
fac_scores_cor <- FACTOR_SCORES(test_models$baseline$cormat, f = EFA_cor)
#> ℹ `x` is a correlation matrix; factor scores cannot be computed. Only factor
#>   weights and score diagnostics are returned. Enter raw data to get factor
#>   scores.
```
