# Generalized Procrustes Analysis consensus target across loading matrices

Internal helper that constructs a Generalized Procrustes Analysis (GPA)
consensus target across a list of loading matrices and returns the
aligned loadings, the centroid target, and convergence diagnostics. Used
by
[`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)
under `target_method = "consensus"` to build a common rotation target
across imputations. Oblique rotations are not supported here: the
iteration is degenerate for oblique transforms with more than one factor
(cf. Lorenzo-Seva & Van Ginkel 2016, who use a Promin step on top of the
centroid rather than iterated oblique Procrustes); callers should pass
the unrotated solutions of an orthogonal rotation, or use
`target_method = "first_target"`.

## Usage

``` r
.gpa_consensus_target(
  unrotated_list,
  init_targets = NULL,
  rotation = c("orthogonal", "oblique"),
  start = 1,
  multi_start = FALSE,
  starts = NULL,
  tol = 0.001,
  loss_tol = 1e-06,
  loss_patience = 5,
  convergence = c("either", "target", "loss", "both"),
  min_iter = 2,
  max_iter = 200,
  alpha = 1,
  match_target = TRUE,
  hyper_cutoff = 0.15,
  verbose = FALSE
)
```

## Arguments

- unrotated_list:

  List of unrotated loading matrices to be aligned. All matrices must be
  numeric, finite, and have identical dimensions.

- init_targets:

  Optional list of starting target matrices. These are typically rotated
  loading matrices from the corresponding analyses. If `NULL`,
  `unrotated_list` is used.

- rotation:

  Character string, either `"orthogonal"` or `"oblique"`.

- start:

  Either a single integer selecting an element of `init_targets`, or an
  explicit target matrix. Used when `multi_start = FALSE`.

- multi_start:

  Logical. If `FALSE`, perform one consensus-target run. If `TRUE`,
  repeat the single-start algorithm for each element of `starts`.

- starts:

  Integer vector selecting elements of `init_targets` used as starting
  targets when `multi_start = TRUE`. If `NULL`, all elements of
  `init_targets` are used. Duplicate entries are removed.

- tol:

  Positive relative Frobenius-norm convergence tolerance for the outer
  target update.

- loss_tol:

  Positive tolerance for the relative change in the outer consensus
  loss. If `NULL`, loss-based convergence is disabled. It cannot be
  `NULL` when `convergence` is `"loss"` or `"both"`.

- loss_patience:

  Positive integer. Number of consecutive iterations with relative loss
  change below `loss_tol` required for loss-based convergence.

- convergence:

  Character string controlling the stopping rule. `"either"` stops when
  either target or loss convergence is satisfied; `"target"` uses only
  target change; `"loss"` uses only loss change; `"both"` requires both.

- min_iter:

  Non-negative integer. Minimum number of outer iterations before
  convergence can be declared.

- max_iter:

  Positive integer. Maximum number of outer consensus iterations.

- alpha:

  Damping factor for the target update. `alpha = 1` uses the full
  centroid update. Smaller values, such as `0.5`, can reduce
  oscillation.

- match_target:

  Logical. If `TRUE`, the updated centroid is signed and column-matched
  to the previous target before convergence is evaluated.

- hyper_cutoff:

  Non-negative cutoff used by `.hyperplane_count()` for summary output.

- verbose:

  Logical; if `TRUE`, print convergence messages for the outer loop.

## Value

A list with the converged target, aligned matrices, pooled loadings,
pooled `Phi`, convergence history, inner-alignment diagnostics, and
hyperplane summaries. If `multi_start = TRUE`, the `multi_start` element
also contains the per-start losses, convergence indicators, run
summaries, all run objects, and between-run Tucker congruence matrices.

## Details

The iteration alternates two steps:

1.  each loading matrix is aligned to the current target with
    [`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md);

2.  the target is updated to the elementwise centroid of the aligned
    matrices.

The outer loop stops when the target stabilises, when the consensus loss
stabilises, or when both criteria are satisfied.

If `multi_start = FALSE`, one consensus run is performed. If
`multi_start = TRUE`, the same engine is repeated for the selected
starting targets and the run with the smallest final mean loss is
returned as the main result; all runs and a between-run congruence
summary are retained in the `multi_start` component.

## References

Gower, J. C. (1975). Generalized Procrustes analysis. *Psychometrika*,
40, 33-51.

Van Ginkel, J. R., & Kroonenberg, P. M. (2014). Using Generalized
Procrustes Analysis for Multiple Imputation in Principal Component
Analysis. *Journal of Classification*, 31, 242-269.

Lorenzo-Seva, U., & Van Ginkel, J. R. (2016). Multiple Imputation of
missing values in exploratory factor analysis of multidimensional
scales: estimating latent trait scores. *Anales de Psicologia*, 32,
596-608.
