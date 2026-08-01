# Multigroup exploratory factor analysis

Fit an exploratory factor analysis in each of several groups at a common
number of factors and bring the per-group solutions into one shared
orientation so their loadings can be compared. Each group is fitted with
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md);
the solutions are then aligned either to a symmetric consensus target or
to a chosen reference group (see *Alignment*).

## Usage

``` r
efa_group(
  x,
  groups = NULL,
  n_factors,
  N = NA,
  reference_group = NULL,
  b_boot = 0L,
  ci = 0.95,
  seed = NULL,
  delta = 0.1,
  invariance = FALSE,
  ...
)
```

## Arguments

- x:

  A data frame or matrix of raw data (with `groups`), or a named list of
  per-group data sets – either raw data frames/matrices or correlation
  matrices (all of one kind).

- groups:

  A vector with one value per row of `x`, giving each row's group. Only
  used when `x` is a single raw data set; leave `NULL` when `x` is a
  list. Rows with a missing group value are dropped with a warning.

- n_factors:

  numeric. The common number of factors extracted in every group.

- N:

  numeric. The number of observations per group, used only for
  correlation-matrix input: either a single value applied to all groups
  or one value per group. Ignored for raw data, where `N` is taken from
  each group's data. Default is `NA`.

- reference_group:

  The group to align the others to (a group name or an integer index).
  If `NULL` (default), orthogonal and unrotated solutions use the
  symmetric consensus target; oblique solutions fall back to the first
  group as reference. Supplying a value forces the reference alignment.

- b_boot:

  numeric. The number of non-parametric bootstrap replicates used to
  form percentile confidence intervals for the between-group Tucker
  congruences. `0` (the default) skips the bootstrap and returns the
  congruence point estimates only. Bootstrapping requires raw data; it
  is skipped with a warning for correlation-matrix input.

- ci:

  numeric. The confidence level for the bootstrap congruence intervals,
  a single value in `(0, 1)`. Default is `0.95`.

- seed:

  numeric or `NULL`. An optional seed making the analysis reproducible.
  It covers the per-group fits – whose rotation may draw random starts –
  whether or not a bootstrap is run, and additionally makes the
  bootstrap independent of the number of parallel workers (the replicate
  fits run with
  [future_lapply()](https://future.apply.futureverse.org/reference/future_lapply.html),
  for which a parallel plan can be set via
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)).
  When supplied, the caller's random-number stream is restored
  afterwards, leaving no side effect. Default is `NULL`.

- delta:

  numeric. The salience threshold for the per-item loading-difference
  flag table: an item's loading on a factor is flagged for a group pair
  when the groups' aligned loadings differ by at least `delta` in
  absolute value. This is a descriptive salience heuristic, not a
  significance test; common alternatives are `0.15` and `0.20`. The
  threshold applies to whatever loading metric the chosen rotation
  produces (pattern coefficients for an oblique rotation). `0` flags
  every cell. Default is `0.1`.

- invariance:

  logical. Whether to add an approximate-invariance verdict per factor
  and group pair from the Lorenzo-Seva and ten Berge (2006) congruence
  bands (see *Value*). Default is `FALSE`.

- ...:

  Additional arguments passed to
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  for every group (for example `estimator`, `rotation`, `cor_method`, or
  an
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  /
  [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  to select a preset).

## Value

An object of class `efa_group`, a list containing:

- loadings:

  A named list of the aligned per-group loading matrices. Their columns
  match the columns of `target` in order and sign.

- target:

  The alignment target: the symmetric consensus target, or the reference
  group's own loadings.

- Phi:

  A named list of the aligned per-group factor intercorrelations for an
  oblique rotation; `NULL` otherwise.

- congruence:

  Tucker congruence between the aligned group loadings, a list with:
  `matrices`, a nested list whose `[[g]][[h]]` element is the
  factor-by-factor congruence matrix between the aligned loadings of
  groups `g` and `h`; `matched`, a groups-by-groups-by-factors array of
  the matched-factor congruences (the diagonal of each pairwise matrix);
  and `degenerate`, a groups-by-groups logical matrix flagging pairs
  whose congruence is undefined (for example, a near-zero factor), for
  which the corresponding entries are `NA`. When `b_boot > 0` (raw
  data), three further elements are added: `matched_se`, the bootstrap
  standard error of each matched congruence; `matched_ci`, a list of
  `lower` and `upper` percentile confidence limits (each a
  groups-by-groups-by-factors array); and `n_boot`, the number of
  bootstrap replicates that aligned in every group and so contributed to
  the intervals (a replicate whose fit did not converge is retained, as
  in
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)).

- diffs:

  A data frame with one row per group pair summarising the differences
  between their aligned loadings: the mean, median, minimum, and maximum
  absolute difference, the root-mean-square difference (`rmse`), and
  `n_flagged`, the number of loading cells whose absolute difference
  reaches `delta`.

- flags:

  A data frame with one row per group pair, item, and factor giving the
  signed loading difference (`diff`), its absolute value (`abs_diff`),
  and `flagged`, whether it reaches `delta`. When a bootstrap was run
  (`b_boot > 0`, raw data), `ci_lower`, `ci_upper`, and `ci_excludes_0`
  add the percentile confidence interval for the difference and whether
  it excludes zero; otherwise these are `NA`.

- invariance:

  When `invariance = TRUE`, a data frame with one row per group pair and
  factor giving the matched Tucker congruence (`phi`), its bootstrap CI
  lower bound (`phi_lower`, `NA` without a bootstrap), and an
  approximate-invariance `verdict` based on the Lorenzo-Seva and ten
  Berge (2006) similarity bands: `phi >= 0.95` is "equal" and
  `[0.85, 0.95)` is "fair"; congruences `< 0.85`, below their bands, are
  labelled "incongruent". The verdict is read from `phi_lower` when a
  bootstrap is available (conservative) and from `phi` otherwise.
  Tucker's congruence is invariant to a proportional rescaling of a
  factor's loadings, so a factor can be graded "equal" even when one
  group's loadings on it are uniformly stronger; read the verdict
  alongside `diffs`. `NULL` when `invariance = FALSE`.

- efa:

  The named list of per-group
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  objects (each retains its own diagnostics, e.g. `heywood`).

- alignment:

  The alignment result: the consensus object (see
  [`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md)),
  or a list with the reference group, the target, and the per-group
  Procrustes results. On the consensus path this is the raw record of
  the Generalized Procrustes iteration, so its `target` and
  `aligned_loadings` are in the orientation that iteration converged to,
  before the gauge described under *Alignment* is applied; the gauged
  matrices are `target` and `loadings` above.

- settings:

  A list of the settings used, including the per-group `N`, the
  alignment method, the group that seeded the consensus frame
  (`alignment_start`, `NULL` on the reference path), the orientation the
  shared frame was put in (`gauge`: the rotation's own name,
  `"principal_axes"`, or `"identity"` for a single factor; `NULL` on the
  reference path), the rotation, the estimator, the input type, and
  whether a bootstrap is available (`can_bootstrap`, `FALSE` for
  correlation-matrix input).

## Details

### Input

Groups can be supplied in two ways: raw data together with a grouping
vector (`x` a data frame or matrix, `groups` one value per row), or a
named list of per-group data sets in `x` (with `groups` left `NULL`).
The list may hold raw data frames or correlation matrices (supply `N`),
but not a mix of the two. All groups must contain the same items in the
same order; a different item set or order is an error rather than being
silently reordered.

Every group is fitted at the same `n_factors` (a common-\\k\\ multigroup
model). Extra arguments in `...` (for example `estimator`, `rotation`,
`cor_method`, or an
[`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
/
[`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
carrying the tuning knobs) are forwarded unchanged to each
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
call, so the estimator and rotation are common to all groups.

The \\k\\-factor model must be identified for the shared item set: its
degrees of freedom \\((p - k)^2 - (p + k)) / 2\\ must be non-negative.
Unlike a single
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
fit – which only warns on an under-identified model – a multigroup fit
aborts, because a shared alignment target across an under-identified
group is not interpretable.

### Alignment

A factor solution is identified only up to a rotation of its factors, so
the per-group solutions must be brought into a common orientation before
their loadings can be compared. Two strategies are available and are
chosen automatically:

- **Consensus** (the default for orthogonal rotations and for unrotated
  solutions): a symmetric Generalized Procrustes Analysis target is
  built across all groups (Gower, 1975), and every group's loadings are
  rotated to it. The consensus frame is only identified up to a global
  rotation of its factors, so it is then rotated into a fixed gauge and
  the same transform is applied to every group. The gauge is the simple
  structure of the rotation that was requested, evaluated on the
  consensus target itself, so the shared loadings are in the same kind
  of frame as the per-group solutions they summarise. Where no criterion
  identifies a frame – an unrotated solution, and a `"bifactorT"`
  request with two factors, whose criterion is identically zero with a
  single group factor – the principal-axes orientation of the target is
  used instead. Either way the columns are ordered by decreasing sum of
  squares and signed by their column sums, and the shared orientation,
  and hence every reported congruence, difference, and flag, is
  independent of the order in which the groups are supplied.

- **Reference**: every group's loadings are aligned by Procrustes
  rotation to one reference group's loadings, which are kept fixed. This
  path is used when `reference_group` is given, and is used
  automatically for oblique rotations because the consensus iteration is
  not defined for oblique transforms with more than one factor. When an
  oblique rotation triggers the reference path without an explicit
  `reference_group`, the first group is used and a message reports this;
  the requested rotation is never silently changed.

In both cases the returned per-group loadings share the column order and
sign of the returned `target`.

### Comparing the aligned loadings

Because the per-group loadings share one orientation, they can be
compared cell by cell. `efa_group()` reports a per-pair summary of their
differences (`diffs`) and a per-item, per-factor flag table (`flags`)
marking cells whose absolute difference reaches `delta`. `delta` is a
*descriptive salience heuristic*, not a significance test; a bootstrap
(`b_boot > 0`) additionally reports whether each difference's confidence
interval excludes zero. With `invariance = TRUE`, a per-factor verdict
grades the matched Tucker congruence by the Lorenzo-Seva and ten Berge
(2006) similarity bands – "equal" (`>= .95`) and "fair" (`[.85, .95)`) –
labelling weaker congruences "incongruent"; when a bootstrap is
available the verdict uses the congruence CI lower bound, so a factor is
judged "equal" only if even the lower bound clears the band.

## References

Efron, B., & Tibshirani, R. J. (1993). *An Introduction to the
Bootstrap*. Chapman & Hall.

Gower, J. C. (1975). Generalized Procrustes analysis. *Psychometrika*,
40, 33-51. doi: 10.1007/BF02291478

Lorenzo-Seva, U., and ten Berge, J. M. F. (2006). Tucker's congruence
coefficient as a meaningful index of factor similarity. *Methodology*,
2, 57-64. doi: 10.1027/1614-2241.2.2.57

## See also

Other factor analysis:
[`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md),
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md),
[`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md),
[`plot.efa_group()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_group.md),
[`print.efa_group()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_group.md)

## Examples

``` r
# Raw data split by a grouping vector (unrotated, consensus alignment)
g <- rep(c("g1", "g2"), length.out = nrow(GRiPS_raw))
mg <- efa_group(GRiPS_raw, groups = g, n_factors = 1)
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
mg$loadings
#> $g1
#>                  F1
#> fun       0.7742597
#> friends   0.8561209
#> enjoy     0.8635188
#> hurt      0.7286359
#> part      0.8183418
#> commonly  0.8379167
#> chances   0.7772425
#> attracted 0.8410818
#> 
#> $g2
#>                  F1
#> fun       0.8156701
#> friends   0.8437139
#> enjoy     0.8842662
#> hurt      0.8035149
#> part      0.8161558
#> commonly  0.8247201
#> chances   0.8003353
#> attracted 0.8531488
#> 

# Per-pair difference summary and the per-item salience-flag table
mg$diffs
#>   group_1 group_2 mean_abs_diff median_abs_diff min_abs_diff max_abs_diff
#> 1      g1      g2    0.02499828      0.01697205  0.002185992   0.07487899
#>         rmse n_flagged
#> 1 0.03309813         0
mg$flags
#>   group_1 group_2 indicator factor         diff    abs_diff flagged ci_lower
#> 1      g1      g2       fun     F1 -0.041410308 0.041410308   FALSE       NA
#> 2      g1      g2   friends     F1  0.012407057 0.012407057   FALSE       NA
#> 3      g1      g2     enjoy     F1 -0.020747446 0.020747446   FALSE       NA
#> 4      g1      g2      hurt     F1 -0.074878986 0.074878986   FALSE       NA
#> 5      g1      g2      part     F1  0.002185992 0.002185992   FALSE       NA
#> 6      g1      g2  commonly     F1  0.013196661 0.013196661   FALSE       NA
#> 7      g1      g2   chances     F1 -0.023092808 0.023092808   FALSE       NA
#> 8      g1      g2 attracted     F1 -0.012067010 0.012067010   FALSE       NA
#>   ci_upper ci_excludes_0
#> 1       NA            NA
#> 2       NA            NA
#> 3       NA            NA
#> 4       NA            NA
#> 5       NA            NA
#> 6       NA            NA
#> 7       NA            NA
#> 8       NA            NA

# \donttest{
# Percentile bootstrap confidence intervals for the between-group congruences, with an
# approximate-invariance verdict read conservatively off the congruence CI lower bound
mg_ci <- efa_group(GRiPS_raw, groups = g, n_factors = 1, b_boot = 100, seed = 42,
                   invariance = TRUE)
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
mg_ci$congruence$matched_ci
#> $lower
#> , , F1
#> 
#>           g1        g2
#> g1 1.0000000 0.9973991
#> g2 0.9973991 1.0000000
#> 
#> 
#> $upper
#> , , F1
#> 
#>           g1        g2
#> g1 1.0000000 0.9996716
#> g2 0.9996716 1.0000000
#> 
#> 
mg_ci$invariance
#>   group_1 group_2 factor       phi phi_lower verdict
#> 1      g1      g2     F1 0.9994099 0.9973991   equal

# A named list of correlation matrices sharing the same items, common
# three-factor model, orthogonal rotation -> symmetric consensus target
bands <- list(age_6_8 = WJIV_ages_6_8$cormat, age_14_19 = WJIV_ages_14_19$cormat)
Ns <- c(WJIV_ages_6_8$N, WJIV_ages_14_19$N)
efa_group(bands, n_factors = 3, N = Ns, rotation = "varimax")
#> ── Multigroup exploratory factor analysis ──────────────────────────────────────
#> 
#> 2 groups (age_6_8, age_14_19) · 3 factors · PAF extraction · varimax rotation
#> Aligned to a symmetric consensus target.
#> N: age_6_8 = 825, age_14_19 = 1685
#> 
#> ── Factor congruence ───────────────────────────────────────────────────────────
#> 
#> Tucker's congruence between the aligned group loadings (matched factors):
#> 
#>                        F1    F2    F3
#> age_6_8 vs age_14_19  .981  .989  .995
#> 
#> ── Loading differences ─────────────────────────────────────────────────────────
#> 
#> Absolute loading differences by group pair (salience threshold |Δ| ≥ .1):
#> 
#>                       mean  median   min   max  rmse  flagged
#> age_6_8 vs age_14_19  .048    .033  .001  .206  .065   20/141

# An oblique rotation aligns to a reference group (reported via a message)
efa_group(bands, n_factors = 3, N = Ns, rotation = "promax")
#> Oblique rotations are aligned to a reference group, not a symmetric consensus
#> target.
#> ℹ The consensus target is undefined for oblique rotations with more than one
#>   factor.
#> ℹ Group "age_6_8" is used as the reference; set `reference_group` to choose
#>   another.
#> ── Multigroup exploratory factor analysis ──────────────────────────────────────
#> 
#> 2 groups (age_6_8, age_14_19) · 3 factors · PAF extraction · promax rotation
#> Aligned to reference group age_6_8.
#> N: age_6_8 = 825, age_14_19 = 1685
#> 
#> ── Factor congruence ───────────────────────────────────────────────────────────
#> 
#> Tucker's congruence between the aligned group loadings (matched factors):
#> 
#>                        F1    F2    F3
#> age_6_8 vs age_14_19  .976  .978  .981
#> 
#> ── Loading differences ─────────────────────────────────────────────────────────
#> 
#> Absolute loading differences by group pair (salience threshold |Δ| ≥ .1):
#> 
#>                       mean  median   min   max  rmse  flagged
#> age_6_8 vs age_14_19  .066    .058  .000  .197  .082   32/141
# }
```
