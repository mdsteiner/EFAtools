# Internal single-start GPA-consensus engine

Performs a single GPA-consensus run from one starting target. The
multi-start wrapper
[`.gpa_consensus_target()`](https://mdsteiner.github.io/EFAtools/reference/dot-gpa_consensus_target.md)
dispatches here.

## Usage

``` r
.consensus_target_procrustes_single(
  unrotated_list,
  init_targets = NULL,
  rotation = c("orthogonal", "oblique"),
  start = 1,
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
