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
| O2 | Accept new compiled deps? `LinkingTo: roptim` (phase-1 optimizer), `pbv` (bivariate-normal CDF), `dqrng` (+`sitmo`,`BH`) (parallel RNG); `Imports: lifecycle`, `robustbase` (MCD). | **Yes** to all; **but `dqrng`/`sitmo`/`BH` deferred out of Phase 3** (only `roptim` added there — see Phase 3 decisions); revisit `RcppEnsmallen` later (O3). | Phases 3,5,7,8 |
| ~~O3~~ | **RESOLVED → `roptim`.** Both `roptim` and `RcppEnsmallen` evaluated 2026-06-15: `roptim` wins on the byte-identical gate (wraps R's `lbfgsb`) and native box-constrained L-BFGS-B; `RcppEnsmallen` has no L-BFGS-B (only Augmented Lagrangian) and stays a post-1.0 backend option. | Locked | Phase 3 |
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

### Phase 3 — C++ estimation engine → `0.9.0`
**Goal:** full estimation in C++; fast, allocation-light bootstrap.

**Decisions (resolved 2026-06-15).** Optimizer = **`roptim`** (O3 closed). Because `roptim`
calls R's `lbfgsb` C routine it is **not thread-safe**, so the estimation bootstrap
parallelises at the **R/process level (`future.apply`)**, not via OpenMP across replicates.
Consequently **`dqrng`/`sitmo`/`BH` are deferred out of Phase 3** — estimation is RNG-free,
resampling + Procrustes starts are R-driven, and `future.seed` (L'Ecuyer) already gives the
worker-count-independent guarantee; dqrng lands when a C++ kernel draws RNG in threads
(Phases 5/7). Only `roptim` is added to `LinkingTo` here. Methods in scope: **PAF/ML/ULS**
(GLS/DWLS need the polychoric Γ — Phase 5).

- [ ] `src/estimate.cpp`: single bounded L-BFGS-B loop via **`roptim`**, `const arma::mat&`,
      fused value+gradient from one eigendecomposition (one-slot psi cache), returns
      loadings/psi/Fm/iter/convergence/eigenpairs. Box constraints on `psi` (lower
      `.uniqueness_floor`, upper 1; no log/logit reparam). Replicate `stats::optim`'s exact
      control (parscale/fnscale + L-BFGS-B factr/pgtol/lmm/maxit defaults) so loadings are
      byte-identical. ML then ULS as two functors on the shared driver.
- [ ] Collapse `paf_iter` to one parameterised loop; return status codes; fix off-by-one
      (B10 — recover its exact definition from git history first; the §6 register and
      `PRE_PHASE2_REVIEW.md` no longer hold it).
- [ ] Lean bootstrap path: per replicate skip the analytic RMSEA-CI `uniroot` and the unused
      `orig`/`final` eigendecompositions, naming and vars; **keep** the cheap fit-index point
      values (chi/CFI/TLI/RMSEA/SRMR/RMSR/AIC/BIC/ECVI/CAF) via a new `.gof(..., ci = FALSE)`
      switch, so **bootstrapped fit-index CIs are retained**. The `uniroot` analytic RMSEA CI
      is the one piece that is both the per-replicate hotspot and statistically redundant
      under bootstrap. Full `.gof()` (with analytic CI) runs only for the point estimate;
      `.gof` stays in R (its few O(p³) factorizations are a small fraction of estimation and
      not worth re-validating in C++).
- [ ] Batched C++ Procrustes alignment for the bootstrap/MI arrays (cube entry point;
      `R::rnorm` starts stay serial — thread-safe rotation RNG is a Phase-4 concern).
- [ ] One `seed` arg; bootstrap parallelised with `future.apply` + `future.seed = TRUE`;
      worker-count-independent reproducibility (O13) **without** dqrng. Cap BLAS/worker
      threads ≤2 in examples/tests.
- [x] Expose `minres` as an alias of `uls` — **done** (`MINRES`→`ULS` resolved at parse time).
- [ ] Fold the start-value heuristics into the new optimizer entry (the `.fit_uls` start
      vector; the exact float-equality SMC fallback in `ML.R`) — kept in R, passed to C++.

**Validation:** loadings byte-identical to the R-optim path (`roptim` uses R's lbfgsb);
`test-regression-estimators.R` (vs psych/factanal) green; bootstrap benchmark on
`DOSPERT_raw` shows the target speedup; reproducibility bit-identical at 1 vs N workers.
**Exit:** ML/ULS/PAF run a single R→C++ call per fit; bootstrap is fast and reproducible.

**Unit queue (P3.1–P3.6 — each atomic, plan-mode first, byte-identical where noted):**
- [ ] **P3.1** `estimate.cpp` driver + **ML** via roptim (add `roptim` to LinkingTo).
      *DoD:* ML loadings byte-identical; `test-regression-estimators.R` ML green.
- [ ] **P3.2** **ULS** functor on the shared driver; route `.fit_uls`; drop `uls_helper.cpp`.
      *DoD:* ULS + minres byte-identical.
- [ ] **P3.3** `paf.cpp`: 8 loops → 1; status codes; B10 (recover definition first).
      *DoD:* PAF byte-identical (≈1e-15 vs psych); NEWS only if B10 moves a number.
- [ ] **P3.4** Lean bootstrap: `.gof(ci = FALSE)` per replicate; skip uniroot + unused eigens.
      *DoD:* fit-index/loading/residual CIs preserved within MC tolerance; faster.
- [ ] **P3.5** Batched C++ Procrustes (cube alignment in one call).
      *DoD:* SE/CI statistically equivalent; record RNG-order change for P3.6 NEWS.
- [ ] **P3.6** `seed` arg + `future.apply` bootstrap + reproducibility test + NEWS.
      *DoD:* bit-identical at 1 vs 2 workers; `test-EFA-boot.R` green.

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
