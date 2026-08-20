# Migrating to the efa\_\* interface

`EFAtools` exposes its functionality through a set of lowercase `efa_*`
functions:
[`efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/efa_screen.md),
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md),
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md),
and so on. Earlier releases used uppercase names instead
([`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md),
[`N_FACTORS()`](https://mdsteiner.github.io/EFAtools/reference/N_FACTORS.md),
[`OMEGA()`](https://mdsteiner.github.io/EFAtools/reference/OMEGA.md),
…). The uppercase names have **not** been removed — they still work —
but the `efa_*` functions are now the recommended interface. This
vignette explains the change, gives the full old-to-new mapping, and
walks through the migration with the largest change to the argument
list: [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md)
to
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md).

``` r

library(EFAtools)
```

## Why the efa\_\* Interface

The `efa_*` names are the interface we recommend for new code. They read
consistently, share a common prefix that groups them in tab-completion
and the documentation index, and they are where the package’s
development now happens.

The uppercase names are **superseded**, not deprecated. In practice that
means:

- They keep working with their **original arguments**, and calling them
  emits **no warning or message**. Existing scripts and saved analyses
  continue to run byte-for-byte unchanged.
- Their argument lists are **frozen**. New capabilities are added only
  to the `efa_*` functions, so an uppercase name will never grow a new
  argument.

So there is no urgency to rewrite working code. Migrate a script when
you want the new functions’ additional features (or simply for
consistency); until then the old calls remain valid.

## The Old-to-New Mapping

Every uppercase function has a lowercase counterpart. For most of them
the migration is a pure rename: the arguments are identical and only the
name changes.

| Old name | Recommended name |
|----|----|
| [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md) | [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md) |
| [`N_FACTORS()`](https://mdsteiner.github.io/EFAtools/reference/N_FACTORS.md) | [`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md) |
| [`EFA_AVERAGE()`](https://mdsteiner.github.io/EFAtools/reference/EFA_AVERAGE-superseded.md) | [`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md) |
| [`EFA_POOLED()`](https://mdsteiner.github.io/EFAtools/reference/EFA_POOLED.md) | [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md) |
| [`COMPARE()`](https://mdsteiner.github.io/EFAtools/reference/COMPARE.md) | [`efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/efa_compare.md) |
| [`SL()`](https://mdsteiner.github.io/EFAtools/reference/SL.md) | [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md) |
| [`OMEGA()`](https://mdsteiner.github.io/EFAtools/reference/OMEGA.md) | [`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md) |
| [`FACTOR_SCORES()`](https://mdsteiner.github.io/EFAtools/reference/FACTOR_SCORES.md) | [`efa_scores()`](https://mdsteiner.github.io/EFAtools/reference/efa_scores.md) |
| [`PROCRUSTES()`](https://mdsteiner.github.io/EFAtools/reference/PROCRUSTES.md) | [`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md) |
| [`BARTLETT()`](https://mdsteiner.github.io/EFAtools/reference/BARTLETT.md) | [`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md) |
| [`KMO()`](https://mdsteiner.github.io/EFAtools/reference/KMO.md) | [`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md) |
| [`PARALLEL()`](https://mdsteiner.github.io/EFAtools/reference/PARALLEL.md) | [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md) |
| [`EKC()`](https://mdsteiner.github.io/EFAtools/reference/EKC.md) | [`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md) |
| [`KGC()`](https://mdsteiner.github.io/EFAtools/reference/KGC.md) | [`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md) |
| [`HULL()`](https://mdsteiner.github.io/EFAtools/reference/HULL.md) | [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md) |
| [`SCREE()`](https://mdsteiner.github.io/EFAtools/reference/SCREE.md) | [`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md) |
| [`MAP()`](https://mdsteiner.github.io/EFAtools/reference/MAP.md) | [`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md) |
| [`NEST()`](https://mdsteiner.github.io/EFAtools/reference/NEST.md) | [`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md) |
| [`SMT()`](https://mdsteiner.github.io/EFAtools/reference/SMT.md) | [`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md) |
| [`CD()`](https://mdsteiner.github.io/EFAtools/reference/CD.md) | [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md) |

[`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md)
and
[`efa_scores()`](https://mdsteiner.github.io/EFAtools/reference/efa_scores.md)
are broader than the
[`OMEGA()`](https://mdsteiner.github.io/EFAtools/reference/OMEGA.md) and
[`FACTOR_SCORES()`](https://mdsteiner.github.io/EFAtools/reference/FACTOR_SCORES.md)
they replace, but they cover the same use cases and are the recommended
way to obtain those quantities going forward.

Alongside the renamed functions, the package ships tools that have **no
uppercase predecessor** — they are new, and available only under the
`efa_*` interface (and the two control constructors):

| New function | Purpose |
|----|----|
| [`efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/efa_screen.md) | Data screening and factorability diagnostics in one report |
| [`efa_group()`](https://mdsteiner.github.io/EFAtools/reference/efa_group.md) | Multigroup EFA with factor congruence |
| [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md) | Simulate data from a factor model |
| [`efa_power()`](https://mdsteiner.github.io/EFAtools/reference/efa_power.md) | Power analysis for EFA |
| [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md) | Bundle the estimation tuning knobs (see below) |
| [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md) | Bundle the rotation tuning knobs (see below) |

For a plain rename, the call is unchanged apart from the name. For
example, the Kaiser-Meyer-Olkin criterion:

``` r

cor_mat <- test_models$baseline$cormat

# Recommended name -- exactly the same arguments as KMO():
efa_kmo(cor_mat)
#> 
#> ── Kaiser-Meyer-Olkin criterion (KMO) ──────────────────────────────────────────
#> 
#> ✔ The overall KMO value for your data is marvellous.
#> These data are probably suitable for factor analysis.
#> 
#> Overall: 0.916
#> 
#> For each variable:
#>    V1    V2    V3    V4    V5    V6    V7    V8    V9   V10   V11   V12   V13 
#> 0.900 0.914 0.924 0.932 0.923 0.891 0.928 0.919 0.916 0.892 0.928 0.908 0.922 
#>   V14   V15   V16   V17   V18 
#> 0.905 0.924 0.934 0.907 0.923

# The old name still works and returns the same result:
identical(KMO(cor_mat)$KMO, efa_kmo(cor_mat)$KMO)
#> [1] TRUE
```

## Migrating from `EFA()` to `efa_fit()`

[`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md) exposed
every estimation and rotation setting as a flat argument, which made for
a long and somewhat unwieldy signature.
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
keeps the **primary** choices — the data, factor count, sample size,
estimator, rotation, and a few others — as top-level arguments (see
[`?efa_fit`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
for the full list), and collects the **tuning** knobs into two small
control objects built by
[`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
and
[`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md).

Each flat
[`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md) tuning
argument now lives in one of the two controls:

| Control object | Arguments it holds |
|----|----|
| [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md) | `type`, `init_comm`, `criterion`, `criterion_type`, `max_iter`, `abs_eigen`, `start_method` |
| [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md) | `type`, `normalize`, `precision`, `order_type`, `varimax_type`, `p_type`, `k`, `random_starts` |

A few points make the translation mechanical:

- **The `type` preset feeds both controls.**
  [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md) had a
  single `type` (`"EFAtools"`, `"psych"`, `"SPSS"`, or `"none"`) that
  governed both estimation and rotation. Pass the same `type` to each
  constructor to reproduce it. Because the two controls are independent,
  you *can* now give them different presets, but you do not have to.
- **Three arguments were renamed.** The estimator is selected with
  `estimator` instead of `method` — passing `method` to
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  is an error that points to the new name, while
  [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md) keeps
  its `method` argument. `P_type` is now `p_type`, and `randomStarts` is
  now `random_starts`;
  [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md) still
  accepts those old spellings silently.
- **Rotation-engine extras** — the criterion-specific arguments a
  rotation may take, such as `maxit`, `gam` (oblimin), or `delta`
  (geomin) — are passed through
  [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)’s
  `...` (or, equivalently,
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)’s
  own `...`). Both validate the names: an extra the selected rotation’s
  engine cannot consume (for example `gamma`, a misspelling of oblimin’s
  `gam`) is rejected with an error, where
  [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md)
  silently ignores it.

The side-by-side below reproduces an SPSS-style analysis. The old flat
call and the new control-based call give identical results:

``` r

# Old flat interface
efa_old <- EFA(cor_mat, n_factors = 3, N = 500, type = "SPSS")

# New interface: the SPSS preset travels through the two control objects
efa_new <- efa_fit(cor_mat, n_factors = 3, N = 500,
                   estimate_control = estimate_control(type = "SPSS"),
                   rotate_control = rotate_control(type = "SPSS"))

# Identical numerical result
all.equal(efa_old$rot_loadings, efa_new$rot_loadings)
#> [1] TRUE
```

An individual knob is set on the control it belongs to. A flat
`EFA(..., type = "SPSS", max_iter = 500, k = 3)` becomes:

``` r

efa_fit(cor_mat, n_factors = 3, N = 500, rotation = "promax",
        estimate_control = estimate_control(type = "SPSS", max_iter = 500),
        rotate_control   = rotate_control(type = "SPSS", k = 3))
```

One behavioural difference is worth knowing about.
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
**rejects** a tuning knob passed directly, rather than silently ignoring
it. This turns a common and previously invisible mistake into an
immediate, informative error that names the constructor the knob belongs
to:

``` r

efa_fit(cor_mat, n_factors = 3, N = 500, max_iter = 500)
#> Error: `max_iter` cannot be passed to `efa_fit()` directly.
#> i The estimation and rotation tuning knobs live in `estimate_control()` and `rotate_control()`.
#> i For example: `efa_fit(x, ..., estimate_control = estimate_control(max_iter = 500))`.
```

The same move — collecting flat estimation knobs into
[`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
— applies to the other functions that fit a model internally:
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md),
[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md),
and the retention criteria that fit a model or use EFA-based eigenvalues
—
[`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
[`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
[`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
[`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md),
[`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md),
and
[`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md).
[`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md)
is the exception: it keeps its flat estimation and rotation arguments,
and only renames `P_type` to `p_type`.

## What the `type` Presets Replicate

The `type` presets are more than shorthand for a bundle of defaults.
`"SPSS"` and `"psych"` follow how SPSS’s FACTOR procedure and
[`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html) implement
principal axis factoring and promax rotation, and `"EFAtools"`, the
default, combines the settings that performed best in a systematic
comparison of the three ([Grieder and Steiner,
2022](https://doi.org/10.3758/s13428-021-01581-x)). There is one place a
preset alone is not enough: because every preset keeps EFAtools’ Kaiser
normalization, `"psych"` also needs `normalize = FALSE` to reproduce
[`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html)’s promax — and
even then a small residual difference remains, from ordinary convergence
slack between the two implementations.

## Returned Objects Keep Their Legacy Classes

Migrating a call does not break code that inspects or dispatches on the
result. The renamed functions attach the **legacy class alongside the
new one**, so [`inherits()`](https://rdrr.io/r/base/class.html) checks
and S3 methods written against the old classes keep working — including
on objects saved from earlier sessions.

``` r

class(efa_new)
#> [1] "efa" "EFA"
inherits(efa_new, "EFA")
#> [1] TRUE
```

The same holds for the other direct renames:
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
objects still inherit `"N_FACTORS"`,
[`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md)
still inherits `"BARTLETT"`, and likewise for
[`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md),
[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md),
[`efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/efa_compare.md),
[`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md),
and
[`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md).
The retention criteria
([`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
[`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md),
and the rest) all share the single `efa_retention` class they have
returned since EFAtools 0.8.0. And the uppercase functions themselves
keep working unchanged: an
[`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md) call
returns the same result as before, now additionally classed `efa` so
that it picks up the shared methods —
[`inherits()`](https://rdrr.io/r/base/class.html) checks against `"EFA"`
are unaffected.

## Where to Next

That is the whole migration: rename the function, and — for
[`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md) only —
move the tuning knobs into
[`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
and
[`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md).
For an overview of the full analysis workflow under the `efa_*`
interface, see the
[EFAtools](https://mdsteiner.github.io/EFAtools/articles/EFAtools.md)
vignette; the individual help pages document each function’s arguments
in full.
