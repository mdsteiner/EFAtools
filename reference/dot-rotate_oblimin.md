# Oblique oblimin factor rotation

Rotate a loading matrix obliquely under the oblimin criterion using a
gradient-projection optimizer along the oblique (column-normalized)
manifold.

## Usage

``` r
.rotate_oblimin(
  L,
  gam = 0,
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

- gam:

  Numeric scalar. The oblimin parameter; `gam = 0` is the quartimin
  criterion.

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

A named list with the rotated loadings, the transformation matrix `Th`
(with `L %*% t(solve(Th))` reproducing the rotated loadings), the factor
correlation matrix `Phi` (`t(Th) %*% Th`), the attained criterion value,
and the convergence and validity flags. The list additionally reports
the criterion value reached at each optimized start in `all_values`,
with a per-start convergence flag in `all_converged`.

## Details

The criterion value `f` and its gradient `dQ/dL` at the rotated loadings
`L = A %*% solve(t(T))` define the search; the engine maps the gradient
to the transformation `T` on the manifold `diag(t(T) %*% T) = 1`,
projects it onto the tangent space, performs a non-monotone line search,
and retracts back onto the manifold by column normalization. `gam = 0`
is the quartimin criterion.

Additional random starts may be requested. To bound runtime the solver
screens each random start by its objective, runs a short triage
optimization on the best-screened starts, and fully optimizes only those
that improve on the current incumbent by at least `triage_improve_tol`.

## References

Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection
algorithms and software for arbitrary rotation criteria in factor
analysis. *Educational and Psychological Measurement*, 65, 676-696.

Jennrich, R. I., & Sampson, P. F. (1966). Rotation for simple loadings.
*Psychometrika*, 31, 313-323.
