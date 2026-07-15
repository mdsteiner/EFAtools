# Rotate a loading matrix to a target using Procrustes alignment

**\[superseded\]**

`PROCRUSTES()` has been superseded by
[`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md),
which is the recommended interface going forward. It remains available
and unchanged so existing code keeps working.

## Usage

``` r
PROCRUSTES(
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

A list identical to the value of
[`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md);
see there for the components.

## See also

[`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md)
