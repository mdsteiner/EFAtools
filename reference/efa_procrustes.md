# Rotate a loading matrix to a target using Procrustes alignment

`efa_procrustes()` aligns one loading matrix to a target loading matrix
with the same dimensions. It is used internally by
[`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md),
but can also be used directly when factor columns must be brought into a
common orientation before averaging or comparing solutions.

## Usage

``` r
efa_procrustes(
  A,
  Target,
  rotation = c("orthogonal", "oblique"),
  S = NULL,
  T_init = NULL,
  oblique_eps = 1e-05,
  oblique_maxit = 1000,
  oblique_max_line_search = 10,
  oblique_step0 = 1,
  oblique_normalize = FALSE,
  oblique_random_starts = 0,
  oblique_screen_keep = 2,
  oblique_triage_maxit = 25,
  oblique_triage_improve_tol = 0
)
```

## Arguments

- A:

  Numeric loading matrix to be aligned.

- Target:

  Numeric target matrix with the same dimensions as `A`.

- rotation:

  Character string, either `"orthogonal"` or `"oblique"`.

- S:

  Optional `k x k` cross-product matrix `crossprod(A)`. Supplying this
  is useful when the same `A` is rotated repeatedly. `S` is used only
  when `oblique_normalize = FALSE`; if Kaiser normalization is
  requested, the cross-product must be recomputed on the normalized
  matrix.

- T_init:

  Optional `k x k` starting transformation matrix for the oblique
  solver. Its columns are normalized internally, and the normalized
  matrix must be well enough conditioned to define a proper factor
  correlation matrix: its smallest singular value must be at least
  `1e-4`, the same floor the solver applies to every candidate it
  evaluates. If `NULL` (the default), the oblique solver is warm-started
  from the closed-form orthogonal Procrustes solution.

- oblique_eps:

  Positive convergence tolerance for the projected-gradient norm in the
  oblique solver.

- oblique_maxit:

  Non-negative integer. Maximum number of projected-gradient updates in
  the full oblique solver.

- oblique_max_line_search:

  Non-negative integer. Maximum number of step-halving attempts after
  the initial line-search step.

- oblique_step0:

  Positive initial step size for the oblique solver.

- oblique_normalize:

  Logical; if `TRUE`, apply Kaiser row normalization to the loadings
  (only) in the oblique solver and back-transform the aligned loadings
  afterwards, leaving `Target` unnormalized (as in
  `GPArotation::targetQ(normalize = TRUE)`).

- oblique_random_starts:

  Non-negative integer. Number of additional random starts used by the
  oblique solver.

- oblique_screen_keep:

  Non-negative integer. Number of random starts retained after cheap
  objective screening and sent to triage optimization.

- oblique_triage_maxit:

  Non-negative integer. Number of short optimization iterations used in
  the triage stage.

- oblique_triage_improve_tol:

  Non-negative scalar. Relative improvement required for a triaged start
  to be promoted to full optimization.

## Value

A list containing aligned `loadings`, transformation matrix `T`, factor
intercorrelation matrix `Phi`, target criterion `value`, convergence
diagnostics, line-search diagnostics, and multi-start summaries. Row and
column names are preserved where possible. When
`oblique_normalize = TRUE` the returned `loadings` are back-transformed
to the original scale, but `value` is the criterion on the
Kaiser-normalized loadings, so it is not
`0.5 * sum((loadings - Target)^2)`.

## Details

For `rotation = "orthogonal"`, the function solves the closed-form
orthogonal Procrustes problem

\$\$\min_T \frac{1}{2}\\A T - B\\\_F^2 \quad \textrm{subject to}\quad
T'T = I,\$\$

where `A` is the loading matrix and `B` is `Target`.

For `rotation = "oblique"`, the function calls the compiled
[`.oblique_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/dot-oblique_procrustes.md)
optimizer. The oblique convention is the same as in
[`GPArotation::targetQ()`](https://rdrr.io/pkg/GPArotation/man/rotations.html):

\$\$L = A T^{-T}, \qquad \Phi = T'T, \qquad diag(\Phi) = 1.\$\$

By default the oblique solver is warm-started from the closed-form
orthogonal Procrustes solution, which resolves the factor permutation
and sign indeterminacy and avoids the poor local minima an identity
start can fall into. Supply `T_init` to override this start. Random
starts are only used for oblique alignment. For one-factor models,
oblique and orthogonal alignment are equivalent, so the function uses
the stable one-factor orthogonal solution instead of calling the oblique
optimizer.

## See also

Other factor rotation:
[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)

## Examples

``` r
## Align an estimated loading matrix to a known target pattern: fit an
## unrotated three-factor model, then rotate its loadings toward the true
## population pattern.
efa_mod <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
                   estimator = "PAF", rotation = "none")
target <- population_models$loadings$baseline

## Orthogonal target rotation (rigid rotation/reflection):
efa_procrustes(efa_mod$unrot_loadings, target, rotation = "orthogonal")
#> $loadings
#>            F1        F2        F3
#> V1  0.5448962 0.2170358 0.1513942
#> V2  0.4507236 0.2171357 0.1636118
#> V3  0.4405048 0.2148779 0.2065093
#> V4  0.5233425 0.1997234 0.2542341
#> V5  0.4375023 0.1708336 0.2686071
#> V6  0.5957496 0.1625323 0.1322423
#> V7  0.2312446 0.5184359 0.1852262
#> V8  0.1828320 0.5369116 0.1648273
#> V9  0.1636529 0.5147575 0.1955207
#> V10 0.1158044 0.5877527 0.1564727
#> V11 0.3114658 0.4055499 0.1874302
#> V12 0.1771391 0.5986127 0.2057337
#> V13 0.1575023 0.2359492 0.5603540
#> V14 0.2326733 0.1257311 0.5018375
#> V15 0.1493127 0.2578064 0.5233723
#> V16 0.2423249 0.1444833 0.5137624
#> V17 0.1742742 0.1498759 0.5806916
#> V18 0.2205114 0.1813649 0.5177058
#> 
#> $T
#>            F1         F2         F3
#> F1  0.5659121  0.5803081  0.5856501
#> F2 -0.0559094  0.7357153 -0.6749794
#> F3  0.8225677 -0.3492357 -0.4487949
#> 
#> $Phi
#>    F1 F2 F3
#> F1  1  0  0
#> F2  0  1  0
#> F3  0  0  1
#> 
#> $value
#> [1] 0.7816832
#> 
#> $convergence
#> [1] TRUE
#> 
#> $valid
#> [1] TRUE
#> 
#> $iterations
#> [1] 0
#> 
#> $kappa_T
#> [1] 1
#> 
#> $Table
#>      iter         f log10_s step
#> [1,]    0 0.7816832      NA   NA
#> 
#> $method
#> [1] "orthogonal_procrustes"
#> 
#> $best_start_index
#> [1] 1
#> 
#> $all_start_indices
#> [1] 1
#> 
#> $all_values
#>   start_1 
#> 0.7816832 
#> 
#> $all_converged
#> start_1 
#>    TRUE 
#> 
#> $all_iterations
#> start_1 
#>       0 
#> 
#> $line_search_failed
#> [1] FALSE
#> 

## Oblique target rotation (lets the aligned factors correlate):
efa_procrustes(efa_mod$unrot_loadings, target, rotation = "oblique")
#> $loadings
#>              F1           F2           F3
#> V1   0.61910495  0.051825776 -0.073905972
#> V2   0.48598086  0.080089984 -0.021446271
#> V3   0.45721811  0.066707991  0.041017771
#> V4   0.55804160  0.003034360  0.078530704
#> V5   0.44494565 -0.010503571  0.140003527
#> V6   0.71323921 -0.031245260 -0.101702122
#> V7   0.08033421  0.545581338 -0.004957233
#> V8   0.01556974  0.593065667 -0.019863392
#> V9  -0.01454789  0.560424060  0.035954147
#> V10 -0.08951131  0.685715730 -0.021788946
#> V11  0.22583866  0.368529478  0.003522170
#> V12 -0.02708029  0.663209734  0.017679319
#> V13 -0.06443691  0.076426030  0.619539738
#> V14  0.09528966 -0.075731827  0.547236063
#> V15 -0.06933901  0.120036105  0.565808245
#> V16  0.09806205 -0.058127457  0.553608593
#> V17 -0.02090578 -0.049484472  0.668138357
#> V18  0.05484079 -0.003288976  0.555443773
#> 
#> $T
#>             F1         F2         F3
#> F1  0.87686430  0.8689951  0.8785807
#> F2 -0.04217862  0.4487928 -0.4042665
#> F3  0.47888408 -0.2084050 -0.2542923
#> 
#> $Phi
#>           F1        F2        F3
#> F1 1.0000000 0.6432594 0.6656709
#> F2 0.6432594 1.0000000 0.6350462
#> F3 0.6656709 0.6350462 1.0000000
#> 
#> $value
#> [1] 0.1727218
#> 
#> $convergence
#> [1] TRUE
#> 
#> $valid
#> [1] TRUE
#> 
#> $iterations
#> [1] 10
#> 
#> $kappa_T
#> [1] 2.623131
#> 
#> $Table
#>       iter         f    log10_s      step
#>  [1,]    0 0.7816832  0.2411910 1.0000000
#>  [2,]    1 0.3275075 -0.1609113 0.2500000
#>  [3,]    2 0.1846550 -0.5678217 0.4074261
#>  [4,]    3 0.1987034 -0.1883165 0.6394707
#>  [5,]    4 0.1755404 -0.8040980 0.1900866
#>  [6,]    5 0.1729979 -1.3026785 0.1575771
#>  [7,]    6 0.1727250 -2.2207728 0.2204627
#>  [8,]    7 0.1727220 -2.7769860 0.2160014
#>  [9,]    8 0.1727218 -3.9808669 0.1736397
#> [10,]    9 0.1727218 -4.9404047 0.1653971
#> [11,]   10 0.1727218 -5.5246516 0.1702069
#> 
#> $line_search_failed
#> [1] FALSE
#> 
#> $method
#> [1] "oblique_procrustes"
#> 
#> $best_start_index
#> [1] 1
#> 
#> $all_start_indices
#> [1] 1
#> 
#> $all_values
#>   start_1 
#> 0.1727218 
#> 
#> $all_converged
#> start_1 
#>       1 
#> 
#> $all_iterations
#> start_1 
#>      10 
#> 
#> $screen_start_indices
#> integer(0)
#> 
#> $screen_values
#> numeric(0)
#> 
#> $n_random_starts
#> [1] 0
#> 
#> $n_screened
#> [1] 0
#> 
#> $n_triaged
#> [1] 0
#> 
#> $n_fully_optimized
#> [1] 1
#> 
```
