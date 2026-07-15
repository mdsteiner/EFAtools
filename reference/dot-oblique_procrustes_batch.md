# Batched oblique Procrustes target rotation over a cube of loading matrices

Align each slice of a loading-matrix cube to a single shared target
using the same oblique target rotation as
[`.oblique_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/dot-oblique_procrustes.md),
in one call. This removes the per-replicate marshalling overhead of
looping
[`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md)
in R over bootstrap or multiple-imputation arrays.

## Usage

``` r
.oblique_procrustes_batch(
  A,
  B,
  eps = 1e-05,
  maxit = 1000L,
  max_line_search = 10L,
  step0 = 1,
  normalize = FALSE,
  random_starts = 0L,
  screen_keep = 2L,
  triage_maxit = 25L,
  triage_improve_tol = 0
)
```

## Arguments

- A:

  Numeric array of dimension `n x m x b`: the `b` loading matrices to
  align.

- B:

  Numeric `n x m` target loading matrix shared across all slices.

- eps:

  Numeric scalar. Convergence tolerance for the projected-gradient norm.

- maxit:

  Integer scalar. Maximum number of full projected-gradient updates.

- max_line_search:

  Integer scalar. Maximum number of step-halving attempts after the
  initial trial step in each line-search phase.

- step0:

  Numeric scalar. Initial step size used in the projected-gradient
  update.

- normalize:

  Logical scalar. If `TRUE`, apply Kaiser normalization to the loadings
  (only) before rotation and reverse it afterwards, leaving the target
  unnormalized (ignored for single-factor slices).

- random_starts:

  Integer scalar. Number of additional random starts per slice.

- screen_keep:

  Integer scalar. Number of screened random starts retained for triage
  optimization.

- triage_maxit:

  Integer scalar. Number of short optimization iterations used in the
  triage stage.

- triage_improve_tol:

  Numeric scalar. Relative improvement required for a triaged start to
  be promoted to full optimization.

## Value

A named list with the aligned-loadings array `loadings` (`n x m x b`),
the factor-correlation array `Phi` (`m x m x b`), and the per-slice
diagnostics `valid`, `convergence`, `value`, `iterations`, and
`line_search_failed`.

## Details

Each slice `A[, , i]` is aligned to `B`. For a single-factor cube the
alignment reduces to the closed-form sign match
`T = sign(crossprod(A_i, B))` with factor correlation `1`, matching the
one-factor short-circuit in
[`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md).
For two or more factors the slice is warm-started from the closed-form
orthogonal Procrustes solution (mirroring
[`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md))
and optimized with the same multi-start oblique solver as
[`.oblique_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/dot-oblique_procrustes.md).
Random starts are drawn serially with `R::rnorm` in the calling process.

Slices are aligned independently. A slice that cannot be aligned (a
non-finite loading matrix, a failed warm-start decomposition, an invalid
fit, or any linear-algebra exception) is reported with `valid = FALSE`
and `NA` for the loadings, factor correlations, and all other per-slice
diagnostics, rather than aborting the whole call, so one degenerate
replicate does not discard the rest.
