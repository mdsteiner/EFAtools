# Oblique Procrustes target rotation using a k x k inner objective

Compute an oblique target rotation for a loading matrix using a
`targetQ`-compatible parameterization and a `k x k` objective.

## Usage

``` r
.oblique_procrustes(
  A,
  B,
  S_r = NULL,
  T_init_r = NULL,
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

  Numeric matrix. Loading matrix to be rotated.

- B:

  Numeric matrix. Target loading matrix with the same dimensions as `A`.

- S_r:

  Optional numeric `k x k` matrix containing `crossprod(A)`. Supplying
  this is useful when the same `A` is rotated repeatedly. Ignored when
  `normalize = TRUE` because the normalized cross-product is different.

- T_init_r:

  Optional numeric `k x k` starting transformation matrix. If `NULL`,
  the identity matrix is used for the primary start.

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
  (only) before rotation and reverse it afterwards; the target is left
  unnormalized, matching `GPArotation::targetQ(normalize = TRUE)`.

- random_starts:

  Integer scalar. Number of additional random starts.

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

A named list containing the rotated loadings, transformation matrix,
factor correlation matrix, target criterion value, convergence
diagnostics, line-search diagnostics, and multi-start summaries.

## Details

The rotated loading matrix is defined as `L = A %*% solve(t(T))`, and
the corresponding factor correlation matrix is `Phi = t(T) %*% T`. The
optimization is carried out over the transformation matrix `T` under the
oblique normalization constraint `diag(t(T) %*% T) = 1`.

Non-invertible candidate transformations are rejected rather than
evaluated through a pseudo-inverse.

Additional random starts may be requested. To reduce runtime, the solver
uses a two-stage strategy for extra starts: cheap objective screening,
followed by short triage optimization, followed by full optimization
only for starts that improve on the current incumbent by at least
`triage_improve_tol`.

The routine is intended for repeated oblique target rotations in
workflows such as bootstrap alignment or consensus alignment of
exploratory factor solutions across multiply imputed datasets. It
follows the same oblique transformation convention as
[`GPArotation::targetQ()`](https://rdrr.io/pkg/GPArotation/man/rotations.html).

## References

Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection
algorithms and software for arbitrary rotation criteria in factor
analysis. *Educational and Psychological Measurement*, 65, 676-696.

Browne, M. W. (2001). An overview of analytic rotation in exploratory
factor analysis. *Multivariate Behavioral Research*, 36, 111-150.
