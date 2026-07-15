# Orthogonal geomin factor rotation

Rotate a loading matrix orthogonally under the geomin criterion using a
gradient-projection optimizer along the orthogonal (Stiefel) manifold.

## Usage

``` r
.rotate_geomin_orth(
  L,
  delta = 0.01,
  eps = 1e-05,
  normalize = TRUE,
  random_starts = 0L,
  maxit = 1000L,
  max_line_search = 10L,
  step0 = 1,
  screen_keep = 5L,
  triage_maxit = 25L,
  triage_improve_tol = 0
)
```

## Arguments

- L:

  Numeric matrix. The unrotated loading matrix (variables by factors).

- delta:

  Numeric scalar. The geomin offset added to the squared loadings; must
  be a positive finite scalar. `delta = 0.01` is the usual default.

- eps:

  Numeric scalar. Convergence tolerance for the projected-gradient norm.

- normalize:

  Logical scalar. If `TRUE`, apply Kaiser normalization before rotation
  and reverse it afterwards.

- random_starts:

  Integer scalar. Number of additional random orthogonal starts.

- maxit:

  Integer scalar. Maximum number of projected-gradient updates.

- max_line_search:

  Integer scalar. Maximum number of step-halving attempts after the
  initial trial step in each line-search phase.

- step0:

  Numeric scalar. Initial step size used in the projected-gradient
  update.

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

A named list with the rotated loadings, the orthogonal rotation matrix
`Th` (with `L %*% Th` reproducing the rotated loadings), the attained
criterion value, and the convergence and validity flags. The list
additionally reports the criterion value reached at each optimized start
in `all_values`, with a per-start convergence flag in `all_converged`.

## Details

The criterion value `f` and its gradient `dQ/dL` at the rotated loadings
`L = A %*% T` define the search; the engine maps the gradient to the
orthogonal transformation `T`, projects it onto the tangent space,
performs a non-monotone line search, and retracts back onto the
orthogonal group via a polar (singular value) projection. The geomin
criterion sums the per-variable geometric mean of the squared loadings
offset by `delta`; it is prone to local minima, so additional random
starts are recommended.

Additional random orthogonal starts may be requested. To bound runtime
the solver screens each random start by its objective, runs a short
triage optimization on the best-screened starts, and fully optimizes
only those that improve on the current incumbent by at least
`triage_improve_tol`.

## References

Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection
algorithms and software for arbitrary rotation criteria in factor
analysis. *Educational and Psychological Measurement*, 65, 676-696.

Browne, M. W. (2001). An overview of analytic rotation in exploratory
factor analysis. *Multivariate Behavioral Research*, 36, 111-150.
