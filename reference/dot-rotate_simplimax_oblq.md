# Oblique simplimax factor rotation

Rotate a loading matrix obliquely under the simplimax criterion using a
gradient-projection optimizer along the oblique (column-normalized)
manifold.

## Usage

``` r
.rotate_simplimax_oblq(
  L,
  k,
  eps = 1e-05,
  normalize = TRUE,
  random_starts = 0L,
  maxit = 1000L,
  max_line_search = 10L,
  step0 = 1
)
```

## Arguments

- L:

  Numeric matrix. The unrotated loading matrix (variables by factors).

- k:

  Integer scalar. The number of "close-to-zero" loadings the criterion
  targets; must be in `[1, nrow(L) * ncol(L)]`. `k = nrow(L)` is the
  usual default.

- eps:

  Numeric scalar. Convergence tolerance for the projected-gradient norm.
  Because the simplimax criterion is only piecewise smooth, the
  projected gradient need not reach this tolerance at the optimum;
  convergence is then reported when the criterion value stalls (the
  non-monotone search described above), so `eps` mainly governs the
  smooth phases of the search.

- normalize:

  Logical scalar. If `TRUE`, apply Kaiser normalization before rotation
  and reverse it afterwards.

- random_starts:

  Integer scalar. Number of random orthogonal starts fully optimized in
  addition to the identity start.

- maxit:

  Integer scalar. Maximum number of projected-gradient updates per
  start.

- max_line_search:

  Integer scalar. Maximum number of step-halving attempts after the
  initial trial step in each line-search phase.

- step0:

  Numeric scalar. Initial step size used in the projected-gradient
  update.

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
and retracts back onto the manifold by column normalization. The
simplimax criterion sums the `k` smallest squared loadings, so it is
minimized when the `k` "close-to-zero" loadings are driven toward zero;
the count `k` is a tuning parameter. Because the set of `k` smallest
loadings is reselected at every evaluation, the criterion is only
piecewise smooth: its gradient jumps as loadings cross the kth-smallest
threshold, so the line search accepts a step whenever it decreases the
largest objective over a short window of recent iterations (a
non-monotone test; Grippo, Lampariello, & Lucidi, 1986), letting the
optimizer step across the kinks where a strictly monotone descent would
stall.

The criterion is strongly prone to local minima, so the solver fully
optimizes the identity start together with `random_starts` random
orthogonal starts and keeps the solution with the lowest criterion
value. Fully optimizing every start – rather than the screen-and-triage
strategy used for the smooth criteria, which assumes the rational start
lies in the global basin – is the standard remedy for the local minima
of complexity-based rotation criteria (Kiers, 1994; Browne, 2001).

## References

Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection
algorithms and software for arbitrary rotation criteria in factor
analysis. *Educational and Psychological Measurement*, 65, 676-696.

Browne, M. W. (2001). An overview of analytic rotation in exploratory
factor analysis. *Multivariate Behavioral Research*, 36, 111-150.

Grippo, L., Lampariello, F., & Lucidi, S. (1986). A nonmonotone line
search technique for Newton's method. *SIAM Journal on Numerical
Analysis*, 23, 707-716.

Kiers, H. A. L. (1994). Simplimax: Oblique rotation to an optimal target
with simple structure. *Psychometrika*, 59, 567-579.
