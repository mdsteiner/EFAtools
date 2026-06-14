# EFAtools — Rewrite & Expansion Plan

> Status: **in progress** — Phase 0 complete; Phase 1 substantially complete (see the
> checklists and the unit queue in §5a). Originally prepared from a full review of the
> package (v0.7.1.9000) plus a structured best-practice research pass; last synced
> against the code at commit 288d3d9 (2026-06). This document is the single source of
> truth for the rewrite; it is meant to be edited as decisions are made and phases
> completed.
>
> **2026-06-13 pre-Phase-2 review:** a full correctness/robustness/test audit landed in
> `PRE_PHASE2_REVIEW.md`. It re-opened **B1** (quartimax still routes to bentlerT) and
> down-graded **B9** to partial (rotmat half unfixed → **B22**), and added **B20–B54** to
> §6 plus a review-spawned unit queue (U-rev-a … U-rev-h) to run **before Phase 3**. The
> maintainer's three flagged areas: EFA_POOLED SE aggregation **verified correct**; HULL
> has one critical elimination off-by-one (**B20**); the oblique Procrustes **core math is
> correct** but the default identity start (**B29**) and multi-start selection (**B34**)
> need fixing.
>
> This file is a development artefact, **not** part of the shipped package (it is
> listed in `.Rbuildignore`).

---

## 0. How to read this

- **Section 1** records the decisions the maintainer has already locked.
- **Section 2** lists the decisions still needed before / during specific phases, each
  with a recommendation. *These are deliberately not pre-decided.*
- **Section 3** states the guiding principles that make a rewrite of this size safe.
- **Section 4** describes the target end-state architecture.
- **Section 5** is the phased, step-by-step roadmap (the strategic map); **§5a** is the
  operational unit queue for the current release — commit-sized units of work, each
  with its own definition of done.
- **Sections 6–11** are cross-cutting registers: bugs, tests, CI, naming/deprecation,
  dependencies, risks.

Throughout, code anchors use `file: function()` form — function names are greppable
and do not rot the way line numbers do.

**Plan maintenance:** when a unit of work lands, tick its box (in §5a and, where it
completes a phase bullet, in §5) and record the commit hash in the same step. When a
new phase opens, its **first unit** is to expand that phase into a §5a-style unit
queue — later phases stay at their current altitude until then.

---

## 1. Decisions locked

| # | Decision | Choice |
|---|----------|--------|
| D1 | Release strategy | **Phased**: a sequence of shippable dev releases. `0.8.x` = infra + safety net + internal refactor + bug fixes (no user-visible behaviour change beyond documented fixes); `0.9.x` = C++ engine, new estimators, polychoric, SEs; `1.0.0` = `efa_*` rename + new user functions. |
| D2 | Renaming scope | **All exports** get lowercase `efa_*` names; the old UPPERCASE names are retained as **throttled deprecated wrappers** (`lifecycle::deprecate_warn`). |
| D3 | Argument API | **Split control objects**: keep ~5 primary args top-level; move tuning knobs into `estimate_control()` / `rotate_control()` (mirroring `estimate_model()` / `rotate_model()`). |
| D4 | New-feature scope for the rewrite | **All four** heavy areas are in scope: (a) polychoric + ordinal/DWLS, (b) analytic/robust SEs, (c) `efa_group` multigroup, (d) `efa_simulate` / `efa_power` / `efa_screen`. |
| D5 | R version floor | **`Depends: R (>= 4.1.0)`** — enables the native pipe `\|>` and `\(x)` lambdas, allowing `magrittr` (and the tidyverse cluster) to be dropped. May bump to `4.2.0` later *only* if the `_` pipe placeholder is wanted. Do not use features newer than 4.1 (e.g. base `%\|\|%` is 4.4 — define it internally). |

Implication of D4: this is an ambitious 1.0. Phasing (D1) is what keeps it tractable —
the heavy features land across `0.9.x` and `1.0.0`, each behind its own milestone with
its own validation gate.

---

## 2. Open decisions still needed (recommendations, not commitments)

Please confirm/adjust these. Each is flagged at the phase where it first bites.

| # | Decision | Recommendation | Bites at |
|---|----------|----------------|----------|
| ~~O1~~ | **RESOLVED → see D5.** R floor set to `R (>= 4.1.0)`. May bump to 4.2.0 later if the `_` pipe placeholder is wanted. | Locked | Phase 0/1 |
| O2 | Accept new compiled deps? `LinkingTo: roptim` (phase-1 optimizer), `pbv` (bivariate-normal CDF), `dqrng` (+`sitmo`,`BH`) (parallel RNG); `Imports: lifecycle`, `robustbase` (MCD). | **Yes** to all; revisit `RcppEnsmallen` later (O3). | Phases 3,5,7,8 |
| O3 | Optimizer backend: `roptim` (byte-identical to current `optim`) first, then evaluate `RcppEnsmallen` once a regression suite exists? | **Yes** — `roptim` for phase 3, reassess `RcppEnsmallen` in a later minor. | Phase 3 |
| O4 | GPArotation: fully replace with a native C++ GPF engine, or keep it as a fallback during transition? | **Keep as a `Suggests`/fallback**, port criteria to C++ one-by-one with 1e-6 cross-checks, drop only after full parity. | Phase 4 |
| O5 | Move `lavaan` and (eventually) `psych` from `Imports` to `Suggests` (gated with `requireNamespace`)? | **Yes** for `lavaan` now; `psych` once `cor.smooth`/`factor.scores` are reimplemented. | Phase 1/3 |
| O6 | Behaviour-changing bug fixes (quartimax→quartimax, Bartlett χ² correction, ULS objective). Fix outright (NEWS entry), or keep current numbers behind a legacy/`type` flag for continuity? | **Fix outright** with prominent NEWS entries; only the SPSS/psych *compatibility presets* need to preserve their documented behaviour. | Phase 2 |
| O7 | χ² multiplier convention: standardise on lavaan/factanal (Bartlett-corrected ML; `N` vs `N-1`)? | **Yes**, with explicit parity tests vs lavaan and a documented choice. | Phase 2/6 |
| O8 | EFA object schema: introduce `new_efa()`/`validate_efa()` and standardise fields (keeping current names where possible, NA-filling inapplicable ones)? | **Yes** — stable, documented contract; additive where possible to avoid breaking `$`-indexing users. | Phase 1 |
| O9 | Deprecation timeline for the UPPERCASE wrappers. | `1.0.0` warn (once/session) → ~`1.x` (12–18 mo / 2–3 minors) `always=TRUE` → `2.0.0` `deprecate_stop`. | Phase 8 |
| O10 | FIML scope. Research rates full casewise FIML *very-high* effort. | **Two-stage EM-Σ** (`missing = "two.stage"`) in scope for `0.9.x`; **direct casewise FIML** deferred to post-1.0. | Phase 6 |
| O11 | WLS: full WLS vs DWLS only for v1. | **DWLS** for v1 (diagonal weights from the polychoric ACOV); full WLS post-1.0. | Phase 5/6 |
| O12 | EGA / network dimensionality as a retention method. | **Defer**; optionally add an `EGAnet` bridge in `Suggests` post-1.0 (avoids a heavy dep). | Post-1.0 |
| O13 | Reproducibility: adopt `dqrng` per-replicate streams (and `future.seed=TRUE` at R level). This **changes** exact seed→result reproducibility vs current versions. | **Yes** — document the break; guarantee thread-count-independent reproducibility going forward. | Phase 3/7 |
| O14 | Simulation model-error default method (`CB` exact-RMSEA vs `TKL` RMSEA+CFI). Depend on / vendor `noisemaker`? | Default **CB**; offer `TKL`,`WB`; validate against `noisemaker`/`fungible` (cross-check tests only, not a hard dep). | Phase 7 |
| O15 | `efa_group` defaults: common-`k` vs per-group-`k`; consensus target vs reference group. | Default **common-`k`** + **consensus** target (symmetric); expose `reference_group`. | Phase 7 |
| O16 | Analytic bifactor rotation (Jennrich-Bentler) alongside Schmid-Leiman; analytic-SE coverage for promax (which has no GPA Jacobian → bootstrap only). | Add bi-geomin/bi-quartimin when the GPF engine lands; **promax SEs via bootstrap only**, documented. | Phase 4/6 |
| O17 | `efa_screen` MVN-test scope: §5 Phase 7 lists four self-implemented tests (Mardia, Doornik-Hansen, Henze-Zirkler, Royston). | **Trim to Mardia + Henze-Zirkler** for 1.0; defer Doornik-Hansen/Royston to post-1.0 (each is its own validation burden). | Phase 7 |

---

## 3. Guiding principles

1. **Safety net before surgery.** The package currently has **no snapshot tests** and
   only a couple of cross-package comparisons. Before any refactor: migrate to
   testthat 3e, capture `expect_snapshot()` baselines for every print/format method,
   and add `lavaan`/`psych` regression tests with explicit tolerances. These become the
   contract every later phase must preserve (or intentionally update via snapshot review).
2. **Separate "no behaviour change" from "behaviour change".** Internal refactors
   (Phase 1) must reproduce current numbers bit-for-bit (snapshots green). Behaviour
   changes (Phase 2 bug fixes, new χ² conventions) are isolated, each with a NEWS entry
   and a deliberate snapshot update.
3. **Cross-validate every numerical port.** New C++ estimators/rotations/polychorics are
   validated against `psych`, `stats::factanal`, `GPArotation`, `lavaan`, `polycor`,
   `pbv`, `EFAutilities` to documented tolerances *before* the old path is removed.
4. **One source of truth per concept.** Type presets, rotation-family classification,
   fit indices, factor reflection/ordering, congruence, number formatting, and input
   preparation each get exactly one implementation, consumed everywhere.
5. **Keep the slow/optional out of the hot path.** Bootstrap and simulation re-use a
   lean, allocation-light, fully-compiled estimate→rotate→align path; expensive
   diagnostics (KMO-based CAF, `uniroot` RMSEA CIs, naming, residual matrices) are
   computed only for the point estimate.
6. **Reproducibility is a feature.** One `seed` argument; results identical across
   sequential/parallel/backends/thread-counts (L'Ecuyer / `dqrng` streams bound to the
   replicate index, never to the thread).
7. **Non-breaking-in-practice rename.** Ship deprecated wrappers so no existing user
   script errors on upgrade; migrate **all** internal usages (examples/vignettes/tests/
   README) to `efa_*` in the same release the wrappers land, so the package never emits
   its own deprecation warnings during `R CMD check`.

---

## 4. Target architecture (end state)

### 4.1 The EFA object (`new_efa()` / `validate_efa()`)
A documented S3 constructor that always populates a stable field set (NA where
inapplicable), used by **all** producers (`efa_fit`, `efa_mi`, `efa_group`,
`efa_average`). Fields: `orig_R`, `R` (analysis matrix, e.g. polychoric), `h2`,
`h2_init`, `uniquenesses`, `*_eigen`, `iter`, `convergence`, `heywood` (flag + which),
`unrot_loadings`, `rot_loadings`, `Phi`, `Structure`, `rotmat`, `vars_accounted[_rot]`,
`fit_indices`, `model_implied_R`, `residuals`, `se`/`ci` blocks, `settings`, `call`.
Keep existing names where possible (additive) to avoid breaking `$`-indexing users (O8).

### 4.2 `estimate_model()` — one estimation core
Resolves settings once from a **declarative preset table** (EFAtools/psych/SPSS/none),
calls a method-specific fitter that returns *only* `L`, `psi`, objective `Fm`, `iter`,
`convergence`, `heywood`, and converged eigenpairs, then runs **one** shared
post-processor (sign reflection, naming with `V`-fallback, `.compute_vars`, fit indices,
`model_implied_R`, residuals, communalities). This removes ~80% of `PAF.R`/`ML.R`/`ULS.R`
and is the single seam for moving estimation into C++. Methods: `PAF`, `ML`,
`ULS`/`MINRES` (same fit_type, documented as identical), `GLS`, `DWLS`, (`WLS` post-1.0).

### 4.3 `rotate_model()` — one rotation dispatcher
Resolves rotation settings once; selects the engine **by name from a lookup table**
(native C++ GPF, plus `stats::varimax`/promax special cases, plus a `GPArotation`
fallback during transition); runs **one** shared post-processor
(`.reflect_and_order(loadings, Phi, rotmat, order_type)` producing consistent
loadings/Phi/Structure/colnames). Selecting the engine by name structurally prevents the
quartimax→bentlerT class of bug (Bug B1).

### 4.4 Unified `efa_retention` result class
Every retention criterion returns one shape:
`list(criterion = c(id, label), n_factors = <named numeric>, results = list(<record>...),
settings, status)` where each record has `{label, n_factors, plot_type ("eigen"|"hull"),
x, y, reference?, threshold?, highlight}`. Then:
- **one** `print.efa_retention()` iterates records generically (cli);
- **two** ggplot helpers (`.gg_eigen_plot`, `.gg_hull_plot`) cover all six current plots;
- `efa_retain()` holds a list of these objects and prints/plots via `lapply`;
- a **criterion registry** (`id → {fun, control_slice, label, sub_variants}`) drives the
  orchestrator, progress (cli), docs (`@family`), and the rename — adding a criterion
  becomes one registry entry.

### 4.5 `efa_reliability()` — coefficients over a normalized model spec
Decouple reliability from Schmid-Leiman. Internal spec:
`(g_load?, s_load, u2, map, Phi?, cormat?)`. Front-end **adapters** normalize SL/`schmid`/
lavaan/EFA-correlated-factors/raw-bifactor/manual inputs to that spec; **one**
`.reliability_core()` computes the menu (fixing the two divergent omega-total formulas,
Bug B5b):
- ω_total, ω_hierarchical, ω_subscale (existing);
- Cronbach α (total + per subscale);
- composite/construct reliability (CR / congeneric ω) — generalise the lavaan
  single-factor formula;
- AVE (per factor);
- H index (guarded `1/(1-L²)`, fixing Bug B5a), ECV, PUC (existing bifactor diagnostics);
- GLB and Revelle β: **confirm scope** (O-note: β model-based g-saturation surrogate vs
  true split-half; GLB needs an SDP/optimizer — propose β-surrogate + GLB optional).
Output: one tidy (long) tibble (`coefficient, level, factor, group, value`) + settings;
**one** `print.efa_reliability` (kills the `ncol` typo class, Bug B4).

### 4.6 C++ layout (`src/`)
- `estimate.cpp` — single bounded optimizer loop (`roptim` phase 1) for ML/ULS/GLS/DWLS;
  value+gradient from **one** eigendecomposition per `psi`; `const arma::mat&` (no
  by-value copies); returns loadings/psi/Fm/iter/convergence/heywood/eigenpairs.
- `paf.cpp` — `paf_iter` collapsed from 8 loops to **one** parameterized loop; returns
  status codes (no `stop()`/`warning()` in the hot path); fixes the off-by-one (Bug B10).
- `rotate.cpp` — one templated GPF engine (GPForth + GPFoblq) parameterized by a
  criterion functor returning `(f, Gq)`; random starts, line search, projection,
  Kaiser normalization written once. Built on the proven scaffolding already in
  `oblique_procrustes.cpp` (manifold projection, backtracking line search, random
  orthonormal starts, screen→triage→optimize, robust inverse guards).
- `procrustes.cpp` — keep + add a **batched** entry point (loop over the 3rd array dim in
  C++) for bootstrap/MI alignment, removing per-replicate R `PROCRUSTES()` overhead.
- `polychoric.cpp` — thresholds-once + two-step ρ via a fresh in-house **Brent 1-D
  minimiser** (not root-finding, not vendored `zeroin2`); bivariate-normal CDF via
  **`pbv` (LinkingTo)**, GL-5 (not GL-48); OpenMP over the `n(n-1)/2` pairs; returns the
  matrix **and** the ACOV (Γ) for DWLS/robust SEs; explicit empty-cell/zero-variance/
  near-singular handling + nearest-PD projection.
- `sim.cpp` — one `simulate_cfm()` kernel (MVN + Vale-Maurelli/Fleishman + IG +
  Ruscio-Kaczetow rank-matching + thresholding); used by `efa_simulate`, and refactored
  `CD`/`PARALLEL`/`NEST` to call it. `dqrng` per-stream RNG.
- `se.cpp` — information matrix, rotation Jacobian, sandwich (Γ-based) SEs.
- Delete the large commented-out dead blocks in `parallel.cpp`.

### 4.7 Split `helper.R` (2111 lines, ~50 funcs) into focused files
`fit-indices.R` (`.gof`, RMSEA-CI, SRMR/TLI/ECVI), `averaging.R` (EFA_AVERAGE helpers),
`procrustes-consensus.R` (the ~400-line Procrustes/consensus block), `alignment.R`
(`.align_solution`, congruence — merge the duplicate `.factor_congruence` /
`.tucker_congruence`), `format-helpers.R` (number formatting — consolidate the ~6
duplicated formatters into one `.efa_num()`), `cor-input.R`
(`.prepare_cor_input()` shared by every public function), `presets.R` (the declarative
type tables), `rotation-family.R` (one `.rotation_family()` + canonical name vectors).

### 4.8 Messaging / printing / plotting
- All `stop()/warning()/message()` → `cli_abort/cli_warn/cli_inform` with classed
  conditions (`class=`), bullet vectors, `call=`/`.envir=`. (~280 base calls; ~430
  crayon calls across 43 files.)
- One `.efa_style(x, style)` shim mapping logical styles → cli primitives; route **all**
  coloring through it, then drop `crayon`.
- `format.<class>()` via `cli::cli_format_method()` returning a styled character vector;
  `print` = `cat(format(...), sep="\n")`. `format()` guaranteed plain (no ANSI) for
  clean embedding (fix `format.EFA` returning colored `capture.output`, Bug-adjacent).
- Colored loading **matrix**: keep the existing pad-then-style renderer (there is no
  native cli per-cell-threshold table), but pad with `cli::ansi_align`/`ansi_nchar` so
  width math is correct once cells carry ANSI; generalise the `format.LOADINGS` engine
  into one `.efa_format_matrix()` used for loadings, SL loadings, Phi, variances, COMPARE
  (kills the `\t`-misalignment in SLLOADINGS/`.get_compare_matrix`).
- Trim `print.EFA()` (already a `compact/standard/full` engine) to a `standard` default;
  add `summary.EFA()` returning a `summary.EFA` object rendered by `print.summary.EFA`
  with the full diagnostics/CIs.
- All six base-graphics plots → ggplot2 returning ggplot objects; `print` no longer
  auto-plots; `viridisLite`/`graphics` dropped.

### 4.9 Dependency target
Drop `crayon`, `progress`, `viridisLite`, `graphics`, `magrittr`, `dplyr`, `tidyr`,
`tibble`, `stringr` (replace with base/cli). Move `lavaan` (and eventually `psych`) to
`Suggests` (gated). Add `Imports: lifecycle`, `robustbase`; `LinkingTo: roptim`, `pbv`,
`dqrng`, `sitmo`, `BH`; `Suggests: vdiffr`, and (optional) `covsim`/`noisemaker`/`MVN`/
`energy`/`EGAnet` only for cross-checks/bridges. Target: ~21 Imports → ~8–10.

---

## 5. Phased roadmap

Each phase lists **goal**, **steps** (checklist), **validation**, **exit criteria**.
Releases: Phases 0–2 → **0.8.0**; Phases 3–7 → **0.9.x** (one minor per major chunk is
fine); Phase 8–9 → **1.0.0**.

### Phase 0 — Safety net & infrastructure → `0.8.0`
**Goal:** a green, modern CI and a regression net, with zero behaviour change.
- [x] Replace `.github/workflows/R-CMD-check.yaml` with r-lib/actions **v2**
      (`checkout@v4`, `setup-r@v2`, `setup-r-dependencies@v2`, `check-r-package@v2`),
      matrix `{ubuntu-latest release/devel/oldrel-1, macos-latest, windows-latest}`.
      Remove the manual depends/cache/sysreq steps and `options(crayon.enabled=TRUE)`.
- [x] Add `test-coverage.yaml` (covr + Codecov) and `pkgdown.yaml` (gh-pages).
- [x] `Config/testthat/edition: 3` in DESCRIPTION; replace `expect_is()` →
      `expect_s3_class()`/`expect_type()`; fix 3e fallout.
- [x] Capture `expect_snapshot()` baselines for **every** `print.*`/`format.*` against
      bundled objects (`test_models`, `GRiPS_raw`, …), wrapped in
      `local_reproducible_output()` (`cli.num_colors=1`). *(Gap found in the 2026-06
      review: `EFA_POOLED` print/format was missed and has no tests at all — closed by
      unit **U1**.)*
- [x] Add `tests/testthat/test-regression-*` comparing current estimators/rotations to
      `psych::fa`, `stats::factanal`, `GPArotation`, with tolerances; gate with
      `skip_if_not_installed`/`skip_on_cran`. **These capture current numbers as the
      pre-refactor contract.**
- [x] Enable `Roxygen: list(markdown = TRUE)`; run `roxygen2md::roxygen2md()`; fix
      `\eqn`/`\deqn` by hand; re-document; verify man-page count and `\link` targets.
- [x] Add `README` badges (CRAN, R-CMD-check, Codecov, lifecycle); ORCIDs in `Authors@R`.
- [x] Add `lifecycle` to Imports; `usethis::use_lifecycle()`.
- [x] Switch to R-version `>=4.1`, plan the `\|>` switch.

**Exit:** CI green on all platforms; full snapshot + regression suite passing on the
**unchanged** code.

### Phase 1 — Internal refactor, **no behaviour change** → `0.8.x`
**Goal:** collapse duplication behind the snapshot/regression net; numbers unchanged.
- [x] `presets.R`: declarative type tables (EFAtools/psych/SPSS) + one
      `.resolve_settings(type, user, preset)` emitting a single consolidated `cli_warn`.
      Replace the ~25 copy-pasted "type and X specified" blocks in
      `PAF/PROMAX/VARIMAX/ROTATE_OBLQ/ROTATE_ORTH`. *(b0ffa61)*
- [x] Move messages, warnings, and errors from current base versions to `cli_inform`,
      `cli_warn` and `cli_abort`. Strip all crayon aspects and use `cli` conform
      messaging. *(5782836, 34d3d35, 05b5bc3, 9cd73c6 — conditions only; styled
      printing is units U4–U6.)*
- [x] `estimate_model()` + shared `.finalize_fit()` (sign/naming/vars/gof/residuals);
      `.PAF/.ML/.ULS` become thin fitters. (Engine still R/optim here.) *(40132f2)*
- [x] `rotate_model()` + `.reflect_and_order()`; rotation engine selected by name table.
      *(7af40bf — also fixed B1 and B9.)*
- [x] `.rotation_family()` + canonical orth/oblique name vectors; replace the 4 duplicated
      classifications (EFA, EFA_POOLED, `.extract_data`, `.type_grid`). *(e7761e2)*
- [x] `.prepare_cor_input()` shared by EFA + all retention criteria + KMO/BARTLETT/
      FACTOR_SCORES (one PD check + one smoothing; fix double-smoothing). *(b12df62 —
      note: smoothing still surfaces only `psych::cor.smooth`'s own warning, not a
      classed `cli_warn`; tracked as B15 → unit U9a.)*
- [x] Unified `efa_retention` class + criterion registry; refit all 9 criteria to the
      shape; **one** `print.efa_retention()`; two ggplot helpers; `N_FACTORS` orchestrator
      shrinks dramatically. *(245967e, 6c3b616, da1c592, 18cd89a, 660ca34 — snapshot +
      vdiffr coverage in place; B8 and B13 resolved along the way.)*
- [ ] `new_efa()`/`validate_efa()` (O8); route all producers through it. → unit **U7**
- [x] Split `helper.R` per §4.7; merge duplicate congruence funcs; consolidate number
      formatters into `.efa_num()`; delete dead code (commented `.error_ml2`, cat-based
      progress bars, trivial `.boot_fun` wrapper). → unit **U3**
- [x] cli migration of styled output: the `.efa_style()` shim + `cli_format_method`
      formats; remove `crayon` from Imports (148 crayon calls across 10 print files
      remain). → units **U4–U6**
- [x] Rewrite `print.EFA()` to cli (it is the largest crayon site): trim it to a
      `standard` default and add `summary.EFA()` for the full view. Done here, with the
      crayon→cli migration, because both touch the same styled output in `print.EFA.R`;
      the existing `compact/standard/full` engine makes this low-risk and `summary.EFA`
      is additive. Likewise migrate `print.LOADINGS`/`print.SLLOADINGS` and the shared
      matrix renderer (`.efa_format_matrix`). → units **U4, U5**
- [x] ggplot2 migration of the six base plots; `viridisLite`/`graphics` removed.
      *(Done via the retention-class refactor — no base-graphics calls remain in `R/`,
      and neither package is in Imports. Remaining related work: extract the plot from
      `print.COMPARE` into a `plot.COMPARE` ggplot method and fix B11 → unit U6.)*
- [x] Dependency slimming (drop tidyverse cluster/`stringr`/`progress`; O1 base pipe);
      `lavaan` → Suggests (gate `OMEGA`/`SL` lavaan paths; `skip_if_not_installed`).
      → unit **U8**

**Validation:** snapshots may change **only** where output formatting intentionally moved
to cli/ggplot (review each diff); the **regression** suite (numeric) must stay green.
**Exit:** `EFA.R`, `PAF/ML/ULS.R`, the rotation files, and `N_FACTORS.R` substantially
smaller; numbers identical to `0.8.0`.

### Phase 2 — Correctness fixes (behaviour-changing) → `0.8.x`
**Goal:** fix the verified bugs (§6), each with a NEWS entry and snapshot update (O6/O7).
- [x] Work through the **open** bugs in §6 via their assigned units (see the register's
      Status column — several were already fixed by the Phase 1 refactors).
      Behaviour-changing fixes get explicit NEWS notes; compatibility presets preserve
      their documented behaviour where that is the point. → units **U9a–U9b** (and the
      bugs folded into U5/U6/U8 per the register)
- [x] Add **SRMR, TLI/NNFI, ECVI** to the fit-index module; reconcile the χ² multiplier
      with lavaan (O7) + parity tests. → units **U9c/U9d**
- [x] Centralise a Heywood detector in `estimate_model()` post-processing; surface in
      `settings`/`summary` via `cli_warn` (consistent across PAF/ML/ULS). → unit **U9d**

**Exit:** regression suite re-baselined against `psych`/`lavaan`/`factanal` with the
corrected conventions; NEWS documents every changed number.

### §5a — Unit queue for the remainder of `0.8.x` (operational)

One unit = one reviewable, committable change (the per-change workflow). Definition of
done for **every** unit: only the listed files change; `devtools::test()` green;
snapshots change only as stated; **numeric results may change only where the unit says
so** (and then with a NEWS entry). Tick + commit hash on completion. The order below is
dependency order; U1/U2 come first because they extend the safety net the later units
rely on.

- [x] **U1 — EFA_POOLED safety net.** `EFA_POOLED` currently has **no tests at all**:
      add print/format snapshots + basic structure/pooling unit tests before anything
      touches it. Files: `tests/testthat/` only. Numbers: n/a (tests only).
- [x] **U2 — Direct unit tests for the new internals.** Boundary-case tests for
      `.estimate_model()` / `rotate_model()` (single factor, `k = p`, invalid preset
      combinations) and `.prepare_cor_input()` (raw vs cormat, `N_policy`, smoothing
      path). Files: tests only. Numbers: n/a.
- [x] **U3 — Split `helper.R`** per §4.7; merge `.factor_congruence`/
      `.tucker_congruence`; consolidate the ~6 number formatters into `.efa_num()`;
      delete dead code (commented `.error_ml2`, cat-based progress bars, trivial
      `.boot_fun` wrapper). Pure file moves + dedup. Numbers: no; snapshots: no.
- [x] **U4 — `.efa_format_matrix()`** + migrate `print.LOADINGS`/`print.SLLOADINGS` to
      cli (pad with `cli::ansi_align`/`ansi_nchar` so width math survives ANSI).
      Numbers: no; snapshots: formatting diffs only (review each).
- [x] **U5 — `print.EFA` → cli + `summary.EFA`.** Trim `print.EFA` to a `standard`
      default; add `summary.EFA`/`print.summary.EFA`; fix **B16** (`format()` returns
      plain text) and **B17** (`residuals.EFA` becomes a pure extractor). Numbers: no;
      snapshots: formatting diffs only.
- [x] **U6 — Remaining prints → cli; drop `crayon`.** OMEGA (incl. **B4**), COMPARE
      (incl. **B11**; move the plot into a `plot.COMPARE` ggplot method — `print` no
      longer auto-plots), KMO, BARTLETT, SL, EFA_AVERAGE, EFA_POOLED; add the
      `.efa_style()` shim; remove `crayon` from Imports. Numbers: no; snapshots:
      formatting diffs only.
- [ ] **U7 — `new_efa()`/`validate_efa()`** (O8); route all producers through it.
      Additive fields only. Numbers: no; snapshots: no. -> **deferred**
- [x] **U8 — Dependency slimming.** Drop `stringr`/`dplyr`/`tidyr`/`tibble`/`magrittr`
      (base + `|>`)/`progress`; `lavaan` → Suggests (gate the OMEGA/SL lavaan paths;
      `skip_if_not_installed`); pin `GPArotation (>= 2022.4-1)` (**B12**). Numbers: no.
- [x] **U9 — Phase 2 behaviour-changing fixes**, one sub-unit per family, each with
      NEWS + deliberate snapshot updates. Numbers: **yes**, as documented per sub-unit.
      - [x] **U9a:** **B14** (PARALLEL partition/percentile) + **B15** (classed
            `cli_warn` on smoothing in `.prepare_cor_input()`).
      - [x] **U9b:** small guards — **B7** (`isTRUE()` on length>1 in EFA_AVERAGE),
            **B18** (guard `eig_sym` in `nest_sym.cpp`), **B19** (EKC `AM2019` NA
            convention).
      - [x] **U9c:** **B2**/O7 — Bartlett-corrected ML χ², proper residual-based ULS
            statistic, lavaan parity tests.
      - [x] **U9d:** SRMR/TLI/ECVI fit indices + centralised Heywood detector in
            `.estimate_model()` post-processing.
- [ ] **U10 - fix line wrapping in all print methods** e.g., print.KMO or print.BARTLETT
      wrap badly. 

### Phase 3 — C++ estimation engine → `0.9.0`
**Goal:** full estimation in C++; fast, allocation-light bootstrap.
- [ ] Consider doing **U7** constructor and validators
- [ ] `src/estimate.cpp`: single bounded L-BFGS-B loop via **`roptim`** (O2/O3),
      `const arma::mat&`, fused value+gradient from one eigendecomposition, returns
      loadings/psi/Fm/iter/convergence/heywood/eigenpairs. Keep box constraints on `psi`
      (no log/logit reparam).
- [ ] Collapse `paf_iter` to one loop; return status codes; fix off-by-one (B10).
- [ ] Lean bootstrap path (loadings+psi+Fm only; skip naming/vars/gof/residuals/uniroot);
      compute full fit indices only for the point estimate.
- [ ] Batched C++ Procrustes alignment for the bootstrap/MI arrays.
- [ ] `dqrng` per-replicate streams; one `seed` arg; thread-count-independent results
      (O13). Threading model: parallelise across replicates, keep linear algebra serial
      per thread, set BLAS threads to 1 in parallel regions (prefer `RcppParallel`-style
      safety; default ≤2 threads in examples/tests).
- [ ] Expose `minres` as an alias of `uls` (documented identical).
- [ ] Fold the psych-inherited start-value heuristics into the new optimizer entry
      (the `.fit_uls` start vector; the exact float-equality SMC fallback in `ML.R`) —
      do **not** patch these ad hoc before this phase.

**Validation:** loadings byte-identical to the R-optim path (`roptim` uses R's lbfgsb);
bootstrap benchmark on `DOSPERT_raw` shows the target speedup.
**Exit:** ML/ULS/PAF run a single R→C++ call per fit; bootstrap is fast.

### Phase 4 — C++ rotation (GPF) engine → `0.9.x`
**Goal:** native, fast rotations; reduce/remove the GPArotation dependency.
- [ ] `src/rotate.cpp`: one templated GPForth+GPFoblq engine + criterion functors.
      Port order: (1) CF family (covers varimax/quartimax/equamax), (2) oblimin/
      quartimin, (3) geominQ/T, (4) bentler/bifactor, (5) simplimax, (6) infomax,
      (7) **target & pst** (new — §4/research gap).
- [ ] Reimplement promax on top of C++ varimax; retire the bespoke `.VARIMAX_SPSS`
      Jacobi loop + the mismatched `.SV` monitor (keep an SPSS-compat shim/fixtures if
      `type="SPSS"` parity must be exact — O6).
- [ ] Report local-minima count across random starts (cheap once C++).
- [ ] **Reassess the default `randomStarts` (B12).** The R/GPArotation engine makes 100
      random starts — the value local-minima research recommends as sufficient (Hattori,
      Zhang & Preacher, 2017, *Multivariate Behavioral Research*, for geomin) — far too slow
      for a default: a single `EFA()` is fine, but `EFA_AVERAGE` (grid × starts) and
      bootstrap-heavy uses become unusably slow (the full test suite went from a few minutes
      to ~30). So `EFA()`'s default is kept at **10** for now. Once rotation runs in fast C++,
      re-evaluate whether 100 (or another higher value) is feasible as the default, and whether
      `EFA()`'s default and `.rotate_model()`'s internal default (currently 10 vs 100) should
      be unified at it.
- [ ] Replace the `R::rnorm` element loops in the random-start generators
      (`oblique_procrustes.cpp` / the new `rotate.cpp`) with thread-safe draws
      (`arma::randn` or `dqrng` streams) **before** any OpenMP region touches them.
- [ ] Keep `GPArotation` as a `Suggests` fallback per criterion until 1e-6 parity is
      proven across random-start/normalization settings (O4); then drop.

**Validation:** each criterion vs `GPArotation` to 1e-6 on `test_models`;
Structure = pattern·Phi invariant to reflection/ordering.
**Exit:** oblique bootstrap re-rotation runs entirely in C++.

### Phase 5 — Polychoric + ordinal/DWLS → `0.9.x`
**Goal:** ordinal EFA with a bootstrap-fast polychoric matrix and a DWLS estimator.
- [ ] `src/polychoric.cpp`: thresholds-once; two-step ρ via fresh Brent **minimiser**;
      BVN CDF via **`pbv`** (GL-5); OpenMP over pairs; returns matrix **+ ACOV (Γ)**
      (O11). Documented empty-cell correction arg + zero-variance detection (`cli_abort`)
      + nearest-PD projection (`cli_warn`).
- [ ] Polyserial/biserial + Pearson fallback so mixed-type matrices are possible
      (confirm scope — research open question).
- [ ] `cor_method` gains `"poly"`/`"tet"`/`"auto"` across `efa_fit`, `efa_kmo`,
      `efa_bartlett`, retention criteria.
- [ ] `DWLS` fit_type in `estimate_model()` using diagonal weights from Γ; (full `WLS`
      post-1.0, O11).

**Validation:** BVN CDF vs `pbv`/`mnormt` ~1e-7; matrix vs `polycor` two-step / `lavaan`
`lavCor(ordered=TRUE)` / `psych::polychoric(correct=FALSE,global=FALSE)` per the
documented reference ladder; determinism + PD under bootstrap resampling.
**Exit:** ordinal EFA + DWLS validated; ordinal bootstrap fast on `DOSPERT_raw`.

### Phase 6 — Standard errors, fit, missing data → `0.9.x`
**Goal:** analytic SEs that make bootstrap optional; principled missing-data handling.
- [ ] `se.cpp` **Phase A**: information-matrix SEs for unrotated ML loadings/uniquenesses.
- [ ] **Phase B**: rotation Jacobian (Jennrich) → SEs/CIs for rotated loadings, Phi,
      structure coefficients, communalities (oblique CF family first).
- [ ] **Phase C**: sandwich / Satorra-Bentler robust SEs (bread = augmented info; meat =
      Γ, reusing the polychoric/4th-moment Γ) + scaled χ² for the categorical path.
- [ ] `se = c("none","information","sandwich","np-boot")`; bootstrap stays the
      always-available ground-truth fallback (and the only path for promax — O16).
- [ ] Two-stage EM-Σ missing-data path `missing = "two.stage"` (O10); SEs under
      missingness use Γ estimated under the missingness pattern. Direct FIML deferred.

**Validation:** analytic SEs vs `np-boot` and vs `EFAutilities::efa` / `lavaan`
(`se="robust.sem"`, `test="satorra.bentler"`) to documented tolerances.
**Exit:** rotated-loading SEs available analytically and validated.

### Phase 7 — New user functions → `0.9.x` → `1.0.0`
**Goal:** reliability expansion + the four new functions, on the now-fast engine.
- [ ] **`efa_reliability`** per §4.5 (normalized spec + adapters + tidy output + one
      print). Confirm β/GLB scope (O-note in §4.5).
- [ ] **`efa_scores`**: implement regression / Bartlett / Anderson-Rubin / ten Berge
      weights natively in Armadillo (drop the `psych::factor.scores` pass-through);
      report **determinacy** coefficients; add a real print/summary.
- [ ] **`efa_group`**: per-group `efa_fit` → consensus/target alignment (shared
      alignment module with `efa_mi`) → Tucker φ with bootstrap CIs → pairwise
      COMPARE-style diffs (refactor `COMPARE` core into a printless `.compare_loadings`)
      → per-item loading-difference flags → optional approximate-invariance summary.
      Defaults O15. De Roover MGFR (`rotation="multigroup"`) and a lavaan-ESEM bridge are
      **post-1.0** stretch.
- [ ] **`efa_screen`**: KMO (overall+item, reuse `.compute_kmo`), Bartlett, determinant/
      condition number, per-item SMC/variance/%missing, multivariate normality
      (self-implement; scope per **O17** — recommended: Mardia + Henze-Zirkler for 1.0,
      Doornik-Hansen/Royston post-1.0), robust-Mahalanobis
      outliers (`robustbase::covMcd`), sparse/empty-category flags; actionable verdicts
      ("non-normal → use ML-robust or polychoric"); one cli print.
- [ ] **`efa_simulate`**: `simulate_cfm()` kernel (input `(Lambda,Phi,Psi)` or population
      R); marginals `normal`/`VM`(Fleishman)/`IG`/`empirical`(Ruscio-Kaczetow);
      thresholding for ordinal (+ optional `match="polychoric"`); missingness MCAR/MAR/
      MNAR; model error `none/CB/TKL/WB` returning achieved RMSEA/CFI (O14). Refactor
      `CD`/`PARALLEL`/`NEST` to call the shared kernel. Return raw data as plain
      matrix/df; allow returning only the population model.
- [ ] **`efa_power`**: `mode="rmsea"` (closed-form non-central-χ² close/not-close fit +
      required-N, MacCallum-Browne-Sugawara; verify vs `semTools`) and `mode="simulation"`
      (retention hit-rate + structure recovery via congruence, reusing `efa_simulate` +
      the C++ engine, parallel + reproducible).

**Validation:** each new function gets unit tests, snapshot prints, and (where an oracle
exists) cross-checks (`semTools`, `fungible`/`noisemaker`, `psych::factor.scores`,
`lavaan` multigroup).
**Exit:** all four feature areas implemented and validated under their current
(UPPERCASE-free internally) names.

### Phase 8 — Rename + control API + deprecation → `1.0.0`
**Goal:** the public `efa_*` API; non-breaking-in-practice.
- [ ] Make `efa_*` the canonical, exported implementations (full mapping in §9).
- [ ] `estimate_control()` / `rotate_control()` (+ any function-specific controls)
      returning validated classed lists; ~5 primary args stay top-level; old top-level
      tuning args become **deprecated arguments** (`lifecycle::deprecated()`), accepted
      during the transition, one home per knob (D3).
- [ ] `R/EFAtools-deprecated.R`: every UPPERCASE name as a one-line
      `lifecycle::deprecate_warn("1.0.0", "OLD()", "new()")` wrapper forwarding to the new
      function; shared `id=` per rename family; documented under one
      `@keywords internal` topic with `@aliases` redirects.
- [ ] Migrate **all** internal usages (examples, vignettes, tests, README) to `efa_*` so
      `R CMD check` emits no self-deprecation warnings.
- [ ] `_pkgdown.yml` reference grouped by `@family` (Fitting, Retention, Rotation,
      Reliability, Scores & comparison, Screening & simulation, Power, Datasets,
      Deprecated).

**Exit:** both old and new names exported; upgrading breaks no user script (only throttled
warnings).

### Phase 9 — Docs, vignettes, site, release → `1.0.0`
- [ ] `@inheritParams` params-dummy to cut doc duplication; `@returns`/`@family` everywhere.
- [ ] Update the two existing vignettes; add new ones: workflow overview, ordinal/
      polychoric EFA, multigroup (`efa_group`), reliability, simulation & power, and a
      **"Migrating to efa_*"** vignette with an old→new mapping table.
- [ ] Build/deploy pkgdown; `codemeta.json` (+ sync GHA); CITATION/ORCIDs.
- [ ] Prominent NEWS "breaking change" + rename section; CRAN submission;
      announce (R-pkg-devel / social).

### Post-1.0 stretch (tracked, not scheduled)
Direct casewise FIML; full WLS; De Roover MGFR; analytic bifactor (bi-geomin/bi-quartimin)
exposure; EGA retention (Suggests bridge) + Unique Variable Analysis; lavaan-ESEM
syntax bridge; `RcppEnsmallen` backend; GLB via SDP; consensus "recommended n_factors".

---

## 6. Bug-fix register

Status legend: **open** (assigned to a unit/phase), **fixed** (commit), **superseded**
(the buggy code no longer exists), **mitigated** (workaround in place, real fix
scheduled).

| ID | Severity | Bug | Anchor | Status | Fix → where |
|----|----------|-----|--------|--------|-------------|
| B1 | high | `quartimax` called `GPArotation::bentlerT` | `R/rotate_model.R: .orth_engines` (line ~18) | **fixed** (pending commit) | `quartimax` now dispatches to `GPArotation::quartimax` (the true quartimax criterion); `test-regression-rotations.R` re-baselined against it and NEWS entry added. Behaviour-changing: quartimax-rotated loadings now match `psych::fa(rotate="quartimax")`/`GPArotation::quartimax` and differ from the previous bentlerT output. → **U-rev-b**. |
| B2 | high | ML/ULS χ² omits Bartlett correction; ULS χ² uses LS residual SS as if Wishart | `R/helper.R: .gof()` | **open** | Bartlett-corrected ML; proper residual-based ULS statistic. → **U9c** (O7; SE side in Ph6). |
| B3 | med | ULS objective sums lower triangle incl. diagonal (`trimatl`) & differs from reported `Fm`; obj/grad eigenvalue floors inconsistent (`0` vs `eps*100`). *(The earlier "two eigendecompositions" claim is outdated — one each now, but `R` is still passed by value.)* | `src/uls_helper.cpp: uls_residuals()` / `grad_uls()` | **open** | Strictly off-diagonal; one shared floor; `const&`. → **Ph3** (also fixes minres parity). |
| B4 | high | `print.OMEGA` `ncol(x[[1]] == 3)` paren typo | `R/print.OMEGA.R: print.OMEGA()` | **open** | Fixed by the cli rewrite of the OMEGA print. → **U6** (long term: tidy `efa_reliability`, Ph7). |
| B5a | med | H index unguarded `1/(1-L²)` → Inf/NaN on Heywood | `R/OMEGA_helper.R` | **open** | Guard in `.reliability_core`. → **Ph7** (snapshot now). |
| B5b | high | Two divergent ω_total formulas (`.OMEGA_FLEX` vs `.OMEGA_LAVAAN`) | `R/OMEGA_helper.R` | **open** | One core. → **Ph7**. |
| B6 | med | `factor_corres` dimension check is a no-op (`nrow()` of a vector) | `R/OMEGA_helper.R` | **open** | `length(g_load)`/`nrow(s_load)`. → **Ph7**. |
| B7 | high | `isTRUE()` on length>1 vectors in the error-aggregation condition. *(The `temp_corres` half of the original report appears addressed in `helper.R` — re-verify when picked up.)* | `R/EFA_AVERAGE.R` (error aggregation, ~`:606`) | **open** | Elementwise NA-safe condition. → **U9b**. |
| B8 | med | `print.N_FACTORS` read `x$output$…` (field is `x$outputs`) → Bartlett print + scree silently no-op | was `R/print.N_FACTORS.R` (deleted) | **fixed** (660ca34) | Unified `format.N_FACTORS` reads `x$outputs`. |
| B9 | med | Oblique reflection/reorder didn't propagate to Phi/structure/rotmat (`ss_factors` worst) | was `R/ROTATE_OBLQ.R` (deleted) | **fixed** (7af40bf + pending commit) | Phi/Structure propagation fixed & verified (7af40bf); the **rotmat** half fixed via **B22** (U-rev-b, pending commit). |
| B10 | med | PAF `stop()` on neg-eigenvalues aborted whole bootstrap; max-iter off-by-one convergence flag | `src/paf_iter.cpp` | **mitigated** (77da514: replicates tryCatch-guarded in `.boot_fun`) | C++ still `stop()`s; status codes + off-by-one fix land with the engine rewrite. → **Ph3**. |
| B11 | med | `print.COMPARE` `aes_string()` + `size=` (deprecated); plot drawn inside `print` | `R/print.COMPARE.R: print.COMPARE()` | **open** | `aes(.data$…)` + `linewidth`; move into a `plot.COMPARE` method. → **U6**. |
| B12 | med | `GPArotation` `randomStarts` used with no min-version floor | `DESCRIPTION` | **open** | Pin `GPArotation (>= 2022.4-1)` (U8c). `randomStarts` default kept at **10** — raising `EFA()` to the research-recommended 100 is too slow under the R engine; revisit (and unify the 10-vs-100 `EFA()`/`.rotate_model()` defaults) once rotation is fast in C++. → **U8c + Phase 4**. |
| B13 | low | NEST `prob` comparator `<` vs decision `<=` mismatch (+ docstring) | was the old `R/NEST.R` (rewritten) | **superseded** (6c3b616) | Refactored NEST uses one quantile-reference rule (`stats::quantile(…, 1 - alpha)` then `<=`), documented in the help page. |
| B14 | low | PARALLEL `size_vec` partition can go negative (verified: 11 datasets / 7 workers → chunk of −1); percentile off-by-one vs `stats::quantile` | `R/PARALLEL.R` (chunking; percentile) | **open** | Exact integer partition; standardise on `stats::quantile` (matches `psych::fa.parallel`). → **U9a**. |
| B15 | low | Non-PD smoothing not surfaced as a classed condition — relies on `psych::cor.smooth`'s own warning | `R/helper.R: .prepare_cor_input()` | **half-done** (b12df62 centralised the PD check + smoothing) | Add a classed `cli_warn` on smoothing. → **U9a**. |
| B16 | low | `format.EFA` returns ANSI-laden `capture.output` | `R/print.EFA.R: format.EFA()` | **open** | `format()` = plain; `print()` = styled. → **U5**. |
| B17 | low | `residuals.EFA` names a logical arg `print` and prints as a side effect | `R/residuals.EFA.R` | **open** | Pure extractor with `type=`; formatting in `summary.EFA`. → **U5**. |
| B18 | med | `arma::eig_sym` called unguarded — the only C++ entry point without the checked-eigen pattern added in 77da514; a degenerate replicate matrix is undefined behaviour, not a caught error | `src/nest_sym.cpp` (`eig_sym` call) | **open** | Apply the same guarded-eigen pattern as the other C++ helpers. → **U9b**. |
| B19 | low | EKC `AM2019`: `which(lambda <= refs)[1] - 1` yields `NA` when no eigenvalue crosses the reference; `NA` flows into `N_FACTORS` aggregation (downstream highlight guard anticipates it) | `R/EKC.R` (`AM2019` branch) | **open** | Decide the convention (`NA` vs `J`) + guard + boundary test. → **U9b**. |

**B20–B54 — from the 2026-06-13 pre-Phase-2 review** (`PRE_PHASE2_REVIEW.md`; all
**verified numerically** unless noted). Severity = effect in realistic use.

| ID | Severity | Bug | Anchor | Status | Fix → where |
|----|----------|-----|--------|--------|-------------|
| B20 | **critical** | Convex-hull elimination off-by-one (`while (i < nr_s-1)`) never tests the last interior triplet → wrong retained factors on real data (verified on `WJIV_ages_14_19`) | `R/HULL.R: .hull_calc` | **fixed** (pending commit) | Loop bound → `i <= nr_s - 1`; step 5/6 moved above the `<3` check so a hull collapsing below 3 routes to the few-solutions fallback (guarded indexing) instead of crashing. Hand-built geometry test added (remove/keep/collapse). Bundled-data counts unchanged; WJIV RMSEA 12→1. NEWS entry added. → **U-rev-a**. |
| B21 | **critical** | Pooled CAF computed on a zero-diagonal residual matrix (`.compute_kmo` never sees `diag←1`) → NaN or garbage near-1 in **every** call | `R/EFA_POOLED.R: .efa_pooled_fit_indices` | **fixed** (pending commit) | Added `diag(delta_hat) <- 1` before `.compute_kmo` (covers point + bootstrap paths); CAF snapshots re-baselined (`:NA`→value). Single-imputation identity test (pooled CAF == single-EFA CAF) added. No NEWS (function unreleased). → **U-rev-a**. |
| B22 | high | `rotmat` sign flips never propagated (rows+cols permuted, signs dropped) → `L_unrot %*% rotmat ≠ rot_loadings` after any reflection/reorder (the unfixed half of B9) | `R/rotate_model.R: .reflect_and_order` (line ~79); `R/PROMAX.R: .rotate_promax` | **fixed** (pending commit) | `rotmat <- (rotmat %*% diag(signs))[, ord]` in `.reflect_and_order` and both `.rotate_promax` branches; reproduction identity restored to ~1e-16 across all 13 rotations. Only `rotmat` changes (no internal consumer, never printed) — loadings/Phi/Structure/variances/fit identical. Invariant tests added (orth, oblique, promax both branches); `test-ROTATE_OBLQ.R` re-baselined off the buggy expectation. → **U-rev-b**. |
| B23 | high | `FACTOR_SCORES(impute=…)` never imputes (psych `missing=` unset); "means"/"median" are no-ops that silently switch to the sapa algorithm | `R/FACTOR_SCORES.R` | **fixed** (pending commit) | `impute` argument dropped entirely (signature, `psych::factor.scores` call, `settings`, roxygen `@param` + both examples); `FACTOR_SCORES` now scores complete cases only, as documented for the former `"none"`. Native imputation deferred to the Phase-7 scoring rework. Test updated (`settings` now `c("method")`; `impute=` removed from calls). → **U-rev-e**. |
| B24 | high | `.OMEGA_FLEX` overwrites a user `cormat` with `model$orig_R` (`NA` for flexible-input SL) → OMEGA on a flexible-SL object crashes with a raw base error | `R/OMEGA_helper.R: .OMEGA_FLEX` (line ~26) | **fixed** (pending commit) | The SL branch now pulls `orig_R` only when the user supplied no `cormat` and the SL object actually carries one (`is.null(cormat) && is.matrix(model$orig_R)`); a user cormat is honoured and the classed `efa_omega_need_cormat` path stays reachable for flexible SL (`orig_R = NA`). Test added (flexible SL + cormat == manual `model=NULL` path). → **U-rev-e**. |
| B25 | high | `variance="sums_load"` g row: `omega_tot_g` uses the correspondence-zeroed sum denominator while `omega_h_g`/`omega_sub_g` use the unzeroed one → `tot ≠ hier+sub` | `R/OMEGA_helper.R: .OMEGA_FLEX` | **fixed** (pending commit) | `omega_tot_g` now shares the unzeroed total-variance denominator with `omega_h_g`/`omega_sub_g`, restoring `tot = hier + sub` for the g row (verified == 0). Identity test added. → **U-rev-e**. |
| B26 | high | PUC formula assumes complete-g + simple structure; wrong in the two supported violation scenarios (cross-loadings, incomplete g) | `R/OMEGA_helper.R` (`.OMEGA_LAVAAN`/`.OMEGA_FLEX`) | **open** | Compute PUC from the actual cross-pair count. → Ph7. |
| B27 | high | `NEST()` builds `u2 <- sqrt(1 - h2)` from a reference-model EFA with no Heywood guard → NaN → unclassed `Rcpp::stop`; crashes the documented direct PAF call | `R/NEST.R` (line ~111) | **open** | Clamp/guard communalities (use `heywood` field); classed condition + test. → **U-rev-f**. |
| B28 | high | `EFA_POOLED(se="np-boot")` aborts whole pooled bootstrap if any one component-EFA replicate failed | `R/EFA_POOLED.R: .efa_pooled_bootstrap_pool` (loop ~1252) | **open** | Skip/tally non-finite replicates; warn (classed). → **U-rev-f**. |
| B29 | high | Oblique `PROCRUSTES()` defaults to identity start, 0 random starts → silent bad local minima (trapped solutions report `convergence=TRUE`); hits `EFA_POOLED` first_target + EFA oblique bootstrap CIs. Consensus path unaffected (warm orthogonal start) | `src/oblique_procrustes.cpp` (T_primary `eye`); `R/PROCRUSTES.R` | **fixed** (pending commit) | `PROCRUSTES()` now warm-starts the oblique solver from the closed-form orthogonal Procrustes solution when no `T_init` is supplied (`R/PROCRUSTES.R`), as the consensus engine already does; this becomes the primary start for the exported solver, `EFA_POOLED` first_target/bootstrap alignment, and the EFA oblique bootstrap. Self-recovery + identity-trap-rescue tests added; the pre-existing `targetQ` oracle test (which pinned a shared local-minimum trap) now compares at the global optimum via `randomStarts`. Warm start computed via a shared `.procrustes_orthogonal_T()` helper (reused by the consensus engine; `.orthogonal_procrustes` wraps it) and in Kaiser-normalized space when `oblique_normalize = TRUE`. → **U-rev-d**. |
| B30 | high | Greedy `which.max` congruence alignment can map two factors onto one (non-permutation) → corrupts the EFA_AVERAGE mean; correct `.align_solution` (LSAP) exists | `R/averaging.R: .array_reorder` (~344) | **fixed** (pending commit) | Use `.align_solution`/reject non-permutations. → **U-rev-c**. |
| B31 | high | Vector `precision` built into the grid but the worker passes the raw user arg → `assert_number` fails every EFA, swallowed by `try()` → all-NA result | `R/EFA_AVERAGE.R` (future_lapply ~582; single-combo ~546) | **fixed** (pending commit) | Pass `arg_grid$precision[[i]]`; recycling test. → **U-rev-c**. |
| B32 | high | D2 χ² ARIV centred at `sqrt(mean(χ²))` not `mean(sqrt(χ²))` → one-sided inflation of ARIV/FMI/pooled χ² vs Li et al. (1991)/`mitml` | `R/EFA_POOLED.R: .efa_pooled_D2` (line ~647) | **open** | `r <- (1+1/m)*var(sqrt(chis))`; cite Li et al.; re-baseline. → **U-rev-g**. |
| B33 | high→med | CFI/TLI lack `max(δ,0)` noncentrality clamping → deflated/zero CFI when the null over-fits (≠ lavaan/Bentler), internally contradictory with `p_chi`/RMSEA | `R/fit-indices.R: .gof` (CFI ~148); `R/EFA_POOLED.R` ~751 | **open** | Truncate deltas at 0 (per cited Bentler/Kline). Behaviour-changing in over-fit regime; NEWS. → **U-rev-g**. |
| B34 | high→med | Procrustes multi-start selection prefers any converged fit over any non-converged one regardless of objective → discards better solutions; returned flag hides it | `src/oblique_procrustes.cpp: fit_is_better_cpp` (~278) | **fixed** (pending commit) | `fit_is_better_cpp` now selects on the attained objective value first (an invalid candidate never wins; exact ties broken toward the converged fit); the selected fit's own convergence flag is surfaced. Multi-start selection test added (a converged suboptimal local minimum is correctly discarded for a better non-converged fit). Because the selected fit's honest `$convergence` can now be `FALSE` for the best (lowest-objective) alignment, the downstream consumers that previously dropped/flagged on `isFALSE($convergence)` (EFA oblique bootstrap, EFA_POOLED first_target and bootstrap re-alignment) were switched to key on `$valid`, so a valid best-but-non-converged alignment is retained rather than discarded. → **U-rev-d**. |
| B35 | med | Heywood detector (`h2 ≥ 1+eps`) unreachable for ML/ULS (bounded `psi ≥ .005`) → improper ML/ULS boundary solutions never flagged; NEWS "consistent across PAF/ML/ULS" is vacuous | `R/estimate_model.R: .finalize_fit`; `R/ML.R`/`R/ULS.R` bounds | **open** | Flag boundary (`psi` at bound) solutions, or soften the NEWS claim. → **U-rev-g**. |
| B36 | med | ECV/group-H use correspondence-zeroed group loadings (≠ `psych::omega`/Rodriguez 2016) despite the psych-reproduction claim | `R/OMEGA_helper.R: .OMEGA_FLEX` | **fixed** (pending commit) | The group-factor H index and ECV are now computed from the unzeroed group loadings (`input[[j+1]]`); the correspondence map still drives PUC only. `type="psych"` ECV now reproduces `psych::omega` to ~1e-8 (oracle test added); the EFAtools-type SL ECV changes 0.659→0.652 (`test-OMEGA_helper.R` ECV literals re-baselined). Pulled forward from Ph7. → **U-rev-e**. |
| B37 | med | SL second-order Heywood → NaN `s*` columns + cryptic base warning; `suppressWarnings()` swallows the classed `efa_heywood` | `R/SL.R` (EFA branch) | **fixed** (pending commit) | The EFA branch now detects a second-order communality `>= 1+eps` (`diag(L2 %*% t(L2))`) and raises a classed `efa_omega_heywood` abort before the `sqrt(1 - comm)` that produced the NaN columns; non-Heywood results are byte-identical. Classed-error test added (PAF, Heywood-inducing Phi). → **U-rev-e**. |
| B38 | med | SL lavaan branch lacks the psi/theta Heywood guards `.OMEGA_LAVAAN` has; elementwise `sqrt(psi)` valid only for diagonal psi | `R/SL.R` (lavaan branch) | **fixed** (pending commit) | Ported the `any(diag(theta) <= 0) || any(diag(psi) <= 0)` guard (classed `efa_omega_heywood`); the residualized loadings now use `diag(sqrt(diag(psi)), nrow=…)` (diagonal explicit, numerically identical for proper diagonal-psi models). Classed-error test added (variance-Heywood second-order fit, all loadings < 1). The `.OMEGA_LAVAAN` single-source dedup is left to the §5 item. → **U-rev-e**. |
| B39 | med | NEST no-stop boundary reports `max_fac−1` (just-rejected model); reference fits run past the `df>0` bound | `R/NEST.R` (loop exit) | **open** | Report last accepted; stop at df bound. → **U-rev-f**. |
| B40 | med | PARALLEL silently returns `NA` when real eigenvalues never cross the reference | `R/PARALLEL.R: .determine_factors` | **open** | Decide boundary convention + document; classed note. → **U-rev-f**. |
| B41 | med | Smoothed matrices re-flagged NPD on every downstream call (NPD threshold `eps^0.6` stricter than `psych::cor.smooth`) → repeated `efa_cor_smoothed` (B15 incomplete across calls) | `R/cor-input.R: .prepare_cor_input` | **open** | Relax re-check / tag already-smoothed. → **U-rev-f**. |
| B42 | med | `N_FACTORS` suitability gate hard-aborts when `N` missing even for N-free criteria | `R/N_FACTORS.R: N_FACTORS` | **open** | Gate abort on whether an N-requiring item was requested. → **U-rev-f**. |
| B43 | med | Raw-data path: NA in the computed correlation matrix → unclassed base crash | `R/cor-input.R: .prepare_cor_input` | **open** | Classed condition pointing at `use=`. → **U-rev-f**. |
| B44 | med | `.hull_calc` crashes with an opaque base error if any gof value is NA | `R/HULL.R: .hull_calc` | **open** | Drop/short-circuit NA rows + classed warning. → **U-rev-f**. |
| B45 | med | CONSENSUS_PROCRUSTES between-run congruence computed on **unaligned** pooled loadings → uninterpretable diagnostic | `R/PROCRUSTES.R: CONSENSUS_PROCRUSTES` (~325) | **open** | Align before computing congruence. → **U-rev-d**. |
| B46 | med | EFA_AVERAGE permanently overwrites the user's global `progressr` handlers | `R/EFA_AVERAGE.R` (~551) | **fixed** (pending commit) | Scope with `withr`/`on.exit`. → **U-rev-c**. |
| B47 | med | EFA_AVERAGE relays classed EFA warnings once per grid row (spam); EFA.R comment claims a non-existent suppression | `R/EFA_AVERAGE.R` (~577) | **fixed** (pending commit) | Summarize to one warning; fix comment. → **U-rev-c**. |
| B48 | med | EFA_AVERAGE does not extract/tally/average the new fit indices (SRMR, TLI, ECVI, RMSR) | `R/averaging.R: .extract_data` | **fixed** (pending commit) | Add to extraction + averaging maps. → **U-rev-c**. |
| B49 | med | `print.EFA_AVERAGE` auto-plots by default (violates "print never plots") | `R/print.EFA_AVERAGE.R` | **fixed** (pending commit) | `plot = FALSE` default. → **U-rev-c**. |
| B50 | med | Rubin df lacks Barnard–Rubin adjustment → can exceed N; reported FMI is uncorrected λ | `R/EFA_POOLED.R: .efa_pooled_rubin_pool` (~979) | **open** | Apply Barnard–Rubin df (needs complete-data df). → Ph6/U-rev-g. |
| B51 | med | NA `chi_null` unguarded: `det()` underflow in `.null_chisq` crashes EFA/`.gof`/SMT, silently NAs BARTLETT | `R/fit-indices.R: .null_chisq/.gof`; `R/SMT.R:150` | **open** | Guard determinant; classed condition. → **U-rev-f**. |
| B52 | med | `factor_corres.cpp` cross-loading index sets collide for ≥10 factors (indices joined without separator) | `src/factor_corres.cpp` (~88) | **open** | Join with a separator. → **U-rev-f**. |
| B53 | med | Bootstrap re-fires per-replicate warnings (`efa_type_override`, C++ max-iter) `b_boot` times | `R/EFA.R: .boot_fun` | **open** | Suppress within loop; surface once. → **U-rev-f**. |
| B54 | med | Pooled `fit_indices$Fm` is a different estimand than the component `Fm` it is named after (≈2× for ULS); pooled set omits TLI/ECVI; `fit_indices_pooled_algorithm` SEs/CIs statistically incoherent | `R/EFA_POOLED.R: .efa_pooled_fit_indices` | **open** | Recompute/rename consistently; add missing indices; drop or relabel the algorithm-SE block. → **U-rev-g**. |

**B55–B60 — from the print-method audit re-run** (2026-06-13; `PRE_PHASE2_REVIEW.md`
§4/§5/§8/§10; all verified-numerically).

| ID | Severity | Bug | Anchor | Status | Fix → where |
|----|----------|-----|--------|--------|-------------|
| B55 | med | Any 1-factor solution with a rotation prints garbage Phi/Structure/Variances tables — `.rotate_single_factor` NA-fills (scalar `NA`, not `NULL`) and print guards `is.null` only; variances breakage hits orthogonal 1-factor too | `R/rotate_model.R: .rotate_single_factor` (~119); `R/print.EFA.R` (~474/484/512) | **open** | Return `NULL` for inapplicable fields, or add `all(is.na())` print guards. → **U-rev-h**. |
| B56 | med | `format()` leaks ANSI for `LOADINGS`/`SLLOADINGS`/`N_FACTORS`/`efa_retention` (cli_format_method, no `ansi_strip`) when colours on → violates plain-text `format()` contract | `R/print.LOADINGS.R`/`print.SLLOADINGS.R`/`N_FACTORS.R`/`efa_retention.R` (`format.*`) | **open** | Wrap bodies in `cli::ansi_strip` like the EFA family. → **U-rev-h**. |
| B57 | low (dedup) | `.efa_loading_row_order` == `.loadings_row_order` byte-for-byte (both live) | `R/print.EFA.R` (~1951) / `R/print.LOADINGS.R` (~256) | **open** | Keep one, call from both; ~17 lines. → **U-rev-h**. |
| B58 | low (dedup) | Small-primary-gap diagnostic computed twice (count vs rows); `count == length(rows)` → desync risk | `R/print.EFA.R: .efa_small_primary_gap_count` (~1640) / `_rows` (~1801) | **open** | Count returns `length(rows(...))`; ~10 lines. → **U-rev-h**. |
| B59 | low | `print.EFA_AVERAGE` lists `admissible` among "varied settings" (removal vector omits it; masked unless a solution is inadmissible) | `R/print.EFA_AVERAGE.R` (~39–40) | **fixed** (pending commit) | Add `"admissible"` to the removal vector. → **U-rev-c**. |
| B60 | low | Print/formatter cluster: negative-zero rendered with `−` in `.efa_num`/`.numformat`; dead `.print_efa_structure_note` (only a commented call); `.efa_num` padding returns a `cli ansi_string` even when unstyled | `R/format-helpers.R`; `R/print.EFA.R` (~1874) | **open** | Normalise `-0`; remove dead fn; `as.character()` the padded result. → **U-rev-h**. |

*(Further low/doc/test findings live in `PRE_PHASE2_REVIEW.md` §5–§8; not all are mirrored
here.)*

**Review-spawned unit queue (U-rev-a … U-rev-i):** see `PRE_PHASE2_REVIEW.md` §9 for the
atomic, dependency-ordered breakdown. Insert ahead of Phase 3 — these harden the existing
features the C++/SE phases will build on. **All fifteen subsystem reviews are now complete**
(the print-method + dedup reviewer was re-run with verification on 2026-06-13).

---

## 7. Testing strategy

- **testthat 3e** + `_snaps/` for every print/summary/format and for deprecation warnings
  (`expect_snapshot`), under `local_reproducible_output()`.
- **vdiffr** doppelganger SVGs for every ggplot (do **not** `expect_snapshot` plots).
- **Regression/oracle tests** (gated `skip_if_not_installed` + `skip_on_cran`,
  `tolerance=`): estimators vs `psych::fa`/`stats::factanal`/`lavaan`; rotations vs
  `GPArotation`; polychoric vs `polycor`/`lavaan`/`pbv`/`psych(correct=FALSE)`; SEs vs
  `EFAutilities`/`lavaan`; simulation vs `fungible`/`noisemaker`; power vs `semTools`.
- **Reproducibility test**: bit-identical output for a fixed seed at 1 vs N workers.
- **Determinism/PD test** on bootstrap polychoric replicates.
- Edge cases: Heywood/ultra-Heywood, non-PD inputs, empty/sparse cells, single factor,
  underidentified/just-identified models, named non-PD cormat (B-list).
- **Known gaps to close first** (→ units U1/U2): `EFA_POOLED` has **no tests at all**
  (no unit tests, no print/format snapshots — the Phase 0 "every print/format" baseline
  missed it); `.estimate_model()`, `rotate_model()` and `.prepare_cor_input()` are
  exercised only indirectly through `EFA()`.

---

## 8. CI/CD

- `check-standard.yaml` (r-lib v2 matrix), `test-coverage.yaml` (Codecov),
  `pkgdown.yaml` (gh-pages), optional `lint.yaml`.
- Examples/tests/vignettes default to ≤2 threads (CRAN).
- OpenMP gated behind `SHLIB_OPENMP_CXXFLAGS` with a serial fallback.
- Never `git add` `src/*.o`/`*.dll` (already gitignored — keep it that way).

---

## 9. Naming & deprecation

**Mechanism:** `lifecycle::deprecate_warn("1.0.0", "OLD()", "new()")` in one-line wrappers
in `R/EFAtools-deprecated.R`, one `@keywords internal` topic, shared `id=` per family.
Timeline O9: `1.0.0` warn → `1.x` `always=TRUE` → `2.0.0` `deprecate_stop`.

**Mapping (D2 = all exports):**

| Old | New |
|-----|-----|
| `EFA` | `efa_fit` |
| `N_FACTORS` | `efa_retain` |
| `OMEGA` | `efa_reliability` |
| `SL` | `efa_schmid_leiman` |
| `FACTOR_SCORES` | `efa_scores` |
| `EFA_POOLED` | `efa_mi` |
| `COMPARE` | `efa_compare` |
| `KMO` | `efa_kmo` |
| `EFA_AVERAGE` | `efa_average` |
| `BARTLETT` | `efa_bartlett` |
| `PARALLEL` | `efa_parallel` |
| `EKC` | `efa_ekc` |
| `KGC` | `efa_kgc` |
| `HULL` | `efa_hull` |
| `SCREE` | `efa_scree` |
| `MAP` | `efa_map` |
| `NEST` | `efa_nest` |
| `SMT` | `efa_smt` |
| `CD` | `efa_cd` |
| `VARIMAX` | `efa_varimax` |
| `PROMAX` | `efa_promax` |
| `PROCRUSTES` | `efa_procrustes` |
| `CONSENSUS_PROCRUSTES` | `efa_consensus_procrustes` |
| — (new) | `efa_group`, `efa_screen`, `efa_simulate`, `efa_power` |
| — (new control) | `estimate_control`, `rotate_control` |

Also rename the result **classes** (and S3 methods) to match (e.g. `EFA`→`efa`,
retention criteria → `efa_retention`), keeping old class strings on objects during the
deprecation window so existing `inherits(x, "EFA")` / method dispatch keeps working.
*(Confirm class-rename appetite — it is a second back-compat surface beyond function
names; the safe default is to keep the old class string alongside the new one.)*

---

## 10. Dependency migration

| Package | Now | Target | Note |
|---------|-----|--------|------|
| crayon | Imports | **drop** | after cli migration (Ph1) |
| stringr | Imports | **drop** | cli `ansi_*` / base |
| graphics, viridisLite | Imports | **drop** | after ggplot migration |
| magrittr, dplyr, tidyr, tibble | Imports | **drop** | base + `\|>` (O1) |
| progress | Imports | **drop** | unused direct dep |
| lavaan | Imports | **Suggests** (gated) | OMEGA/SL lavaan paths |
| psych | Imports | **Suggests** (later) | after `cor.smooth`/`factor.scores` reimplemented |
| GPArotation | Imports | **Suggests/fallback → drop** | after C++ GPF parity (O4); pin `>=2022.4-1` meanwhile (B12) |
| future, future.apply, progressr | Imports | **keep (consolidate)** | drop redundant `progress`; revisit once C++ threading lands |
| Rcpp, RcppArmadillo, checkmate, rlang, clue, stats | Imports/LinkingTo | **keep** | |
| lifecycle, robustbase | — | **add Imports** | deprecation; MCD outliers |
| roptim, pbv, dqrng, sitmo, BH | — | **add LinkingTo** | optimizer, BVN CDF, parallel RNG (O2) |
| vdiffr | — | **add Suggests** | ggplot snapshots |
| covsim, noisemaker, MVN, energy, EGAnet, semTools, EFAutilities, polycor | — | **Suggests (cross-checks/bridges only)** | never hard deps |

Target Imports: **~8–10** (from 21).

---

## 11. Risks & mitigations

1. **Numerical drift from the C++ optimizer** → start with `roptim` (R's lbfgsb,
   byte-identical); switch backends only behind the regression suite.
2. **Rotation-Jacobian / sandwich SEs are subtle** (multiple criteria, oblique vs orth,
   Kaiser/CF normalization; promax has no GPA Jacobian) → validate against
   `EFAutilities`/`lavaan`; keep `np-boot` as the always-available fallback; ship SEs in
   phases A→B→C.
3. **Polychoric ACOV (Γ) correctness** (harder than the point estimate) → bootstrap is
   the primary uncertainty path; Γ validated vs `lavaan` NACOV before WLS/robust SEs use
   it; explicit empty-cell handling + nearest-PD projection.
4. **Thread safety** (Armadillo LA / R RNG inside parallel regions) → parallelise across
   replicates only, serial LA per thread, BLAS threads=1 in regions, `dqrng` per-stream.
5. **Snapshot fragility** (ANSI/width/locale) → `local_reproducible_output()`,
   `cli.num_colors=1` everywhere.
6. **Self-inflicted deprecation warnings** → migrate all internal usages to `efa_*` in
   the rename release.
7. **Scope (all four feature areas + rename in one rewrite)** → phasing (D1) keeps each
   release shippable; must-haves first, stretch goals explicitly deferred (§5 post-1.0).
8. **Reproducibility break** from new RNG → documented; new guarantee is thread-count
   independence (O13).
9. **De Roover MGFR has no canonical R oracle** → deferred to post-1.0; validate against
   the authors' reference if/when implemented.

---

## 12. Appendix — key current-code anchors

- Orchestration & bootstrap: `R/EFA.R` (method/rotation dispatch now via
  `.estimate_model()` / `rotate_model()`; output assembly; `.boot_*` helpers — boot
  replicates tryCatch-guarded since 77da514).
- Estimators: `R/estimate_model.R` (dispatcher + `.finalize_fit()`); thin fitters in
  `R/PAF.R`, `R/ML.R` (`.fit_ml()`, `.FAout()`), `R/ULS.R` (`.fit_uls()`,
  `.FAout_wls()`); C++ `src/paf_iter.cpp` (8× near-duplicate loop),
  `src/ml_helper.cpp`, `src/uls_helper.cpp` (by-value `R`; one eigendecomposition each
  in objective and gradient, but inconsistent floors + trimatl asymmetry — B3).
- Rotations: `R/rotate_model.R` (engine-by-name tables + `.reflect_and_order()`),
  `R/VARIMAX.R`, `R/PROMAX.R`; `src/oblique_procrustes.cpp` (GPF scaffolding),
  `src/factor_corres.cpp`. (`ROTATE_ORTH.R`/`ROTATE_OBLQ.R` deleted in 7af40bf.)
- Retention: `R/efa_retention.R` (class, criterion registry, `print.efa_retention`,
  `.gg_eigen_plot()`/`.gg_hull_plot()`) + `R/N_FACTORS.R` (registry-driven
  orchestrator) + the 9 criterion files + `src/parallel.cpp`/`nest_sym.cpp`
  (single-threaded; dead commented blocks in `parallel.cpp`; unguarded `eig_sym` in
  `nest_sym.cpp` — B18).
- Fit indices: `R/helper.R: .gof()`, `.rmsr()`, `.compute_caf()`.
- Input preparation: `R/helper.R: .prepare_cor_input()` (shared PD check + smoothing).
- Post-EFA: `R/OMEGA.R`/`OMEGA_helper.R` (two ω paths), `R/SL.R`, `R/FACTOR_SCORES.R`
  (psych pass-through), `R/COMPARE.R`, `R/KMO.R`, `R/BARTLETT.R`.
- Printing: `R/print.EFA.R` (factored `compact/standard/full` engine — easy
  `summary.EFA` split), `R/print.LOADINGS.R` (best table engine — generalise). The
  retention prints are already collapsed into `print.efa_retention`.
- `helper.R` (~2,230 lines, ~50 funcs) — split per §4.7 (unit U3).

---

*End of plan. Edit Section 2 as decisions are confirmed; tick §5/§5a items (with commit
hashes) as units land; update the §6 Status column when a bug is fixed or superseded.*
