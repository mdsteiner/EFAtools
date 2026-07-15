# Rotation Jacobians for analytic rotation standard errors

Forward-difference the warm-started re-rotation map
`A -> (rotated loadings, Phi)` over the unrotated loadings `A` to obtain
the rotation Jacobians used by the analytic standard errors for rotated
loadings (`se = "information"` in
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)).
The full `nrow(A) * ncol(A)` finite- difference loop runs in compiled
code, re-solving the rotation from the converged transformation `T_init`
at each perturbation; the caller forms `J V J'` in R.

## Usage

``` r
.rotation_se_jacobian(
  A,
  T_init,
  method,
  param,
  normalize,
  oblique,
  eps,
  general_col = 0L
)
```

## Arguments

- A:

  Numeric matrix. The unrotated loading matrix at the solution.

- T_init:

  Numeric matrix. The converged transformation that warm-starts each
  re-rotation.

- method:

  Character scalar. The criterion family: one of `"cf"`, `"oblimin"`,
  `"geomin"`, `"bentler"`, `"bifactor"`.

- param:

  Numeric scalar. The criterion's tuning argument (`kappa` for `"cf"`,
  `gam` for `"oblimin"`, `delta` for `"geomin"`); ignored for
  `"bentler"` and `"bifactor"`.

- normalize:

  Logical scalar. Apply Kaiser normalization before rotation and reverse
  it after.

- oblique:

  Logical scalar. Use the oblique (column-normalized) manifold;
  otherwise orthogonal.

- eps:

  Numeric scalar. The forward-difference step on the loadings.

- general_col:

  Integer scalar. For `method = "bifactor"`, the zero-based column
  holding the general factor in the (factor-reordered) reported
  solution; ignored by the other criteria.

## Value

A named list with the Jacobian `J_L` (`pk x pk`), the re-rotated
`base_loadings`, a validity flag, and – when oblique – the Jacobian
`J_Phi` (`k^2 x pk`) and `base_Phi`.

## References

Jennrich, R. I. (1973). Standard errors for obliquely rotated factor
loadings. *Psychometrika*, 38, 593-604.
