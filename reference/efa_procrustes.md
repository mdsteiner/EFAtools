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

  Optional `k x k` nonsingular starting transformation matrix for the
  oblique solver. Its columns are normalized internally. If `NULL` (the
  default), the oblique solver is warm-started from the closed-form
  orthogonal Procrustes solution.

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
