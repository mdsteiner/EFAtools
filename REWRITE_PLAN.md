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
| O4 | GPArotation: fully replace with a native C++ GPF engine, or keep it as a fallback during transition? | **Replace** — Phase 4 ports the **full** reachable set (CF, oblimin/quartimin, geomin, bentler, bifactor, simplimax) to C++ and **demotes GPArotation Imports → Suggests** (gated test oracle only; no runtime dep/fallback). Cross-checks use **relaxed parity** (best-`Q` + ~1e-4 aligned), not 1e-6 iterate-matching. | Phase 4 |
| O5 | Move `lavaan` and (eventually) `psych` from `Imports` to `Suggests` (gated with `requireNamespace`)? | **Yes** for `lavaan` now; `psych` once `cor.smooth`/`factor.scores` are reimplemented. | Phase 1/3 |
| O6 | Behaviour-changing bug fixes (quartimax→quartimax, Bartlett χ² correction, ULS objective). Fix outright (NEWS entry), or keep current numbers behind a legacy/`type` flag for continuity? | **Fix outright** with prominent NEWS entries; only the SPSS/psych *compatibility presets* need to preserve their documented behaviour. | Phase 2 |
| O7 | χ² multiplier convention: standardise on lavaan/factanal (Bartlett-corrected ML; `N` vs `N-1`)? | **Yes**, with explicit parity tests vs lavaan and a documented choice. | Phase 2/6 |
| O8 | EFA object schema: introduce `new_efa()`/`validate_efa()` and standardise fields (keeping current names where possible, NA-filling inapplicable ones)? | **Yes** — stable, documented contract; additive where possible to avoid breaking `$`-indexing users. | Phase 1 |
| O9 | Deprecation timeline for the UPPERCASE wrappers. | `1.0.0` warn (once/session) → ~`1.x` (12–18 mo / 2–3 minors) `always=TRUE` → `2.0.0` `deprecate_stop`. | Phase 8 |
| O10 | FIML scope. Research rates full casewise FIML *very-high* effort. | **Two-stage EM-Σ** (`missing = "two.stage"`) in scope for `0.9.x`; **direct casewise FIML** deferred to post-1.0. | Phase 6 |
| ~~O11~~ | WLS: full WLS vs DWLS only for v1. | **RESOLVED → DWLS** for v1 (per-element weights `= 1/diag(Γ)`, lavaan WLSMV convention); the full Γ machinery is built in Phase 5 (gated `acov` level) so Phase-6 robust SEs can reuse it; full `WLS` post-1.0. Method stays orthogonal to `cor_method`. | Phase 5/6 |
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
- [x] **P3.5** Batched C++ Procrustes (cube alignment in one call). DONE (awaiting commit).
      `.oblique_procrustes_batch(A_cube, B, ...)` loops slices in C++ over the shared
      `run_oblique_multistart` core; `.boot_se_ci()` oblique branch now makes one batched
      call over the non-failed replicates instead of `b` `PROCRUSTES()` round trips.
      Single-factor slices use a closed-form sign match in C++ (no R fallback loop).
      The batch is per-slice robust: an unalignable slice (non-finite loadings, failed
      warm-start SVD, invalid fit, or any linear-algebra exception) is reported
      `valid = FALSE` with NA loadings/Phi and the run continues, preserving the
      per-replicate failure isolation of the loop it replaces. Shared helpers
      (`validate_oblique_scalars`/`kaiser_normalize_rows`/`finalize_oblique`) back both
      the single and batched entry points.
      *Result:* output **numerically identical** to the old per-replicate path
      (max |loading diff| ≈ 1.7e-15, SE diff ≈ 1e-17); `Structure = pattern·Φ` exact;
      ≈1.9× faster on the alignment (after hoisting the shared target cross-product).
      **For P3.6:** the anticipated RNG-order shift did NOT materialise — processing
      non-failed slices in order keeps the `R::rnorm` stream aligned with the old loop.
      The one behavioural change is that the oblique bootstrap **warm start now uses
      Armadillo's `svd` inside C++** (was R's `svd()` via `.procrustes_orthogonal_T`);
      `U Vᵀ` is sign/order-invariant so it matches to ≈1e-15, but exact bit-identity to
      the old path is not guaranteed in SVD-degenerate cases. P3.6's reproducibility
      test should pin the NEW path against itself; cross-version equality is ~1e-15.
- [ ] **P3.6** `seed` arg + `future.apply` bootstrap + reproducibility test + NEWS.
      *DoD:* bit-identical at 1 vs 2 workers; `test-EFA-boot.R` green.
      NEWS: note that oblique bootstrap re-rotation is now a single compiled pass and
      that its warm start moved to C++ `svd` (see P3.5 above).

### Phase 4 — C++ rotation (GPF) engine → `0.9.x`
**Goal:** native, fast rotations; **remove the GPArotation runtime dependency** — port the full
reachable criterion set to C++ and demote GPArotation to a Suggests-only test oracle.

**Decisions (resolved 2026-06-16).** Build one **`vgQ`-criterion-templated GPForth+GPFoblq
engine** on the `oblique_procrustes.cpp` scaffolding (manifold projection, Armijo line search,
screen→triage→optimize multistart, robustness guards, random orthonormal starts are all
reusable — the Procrustes objective is just one criterion functor; **add** the orthogonal
/Stiefel projection variant). **Scope = the full reachable criterion set**, ported
family-by-family: CF (quartimax/equamax), oblimin/quartimin, geomin (T/Q), bentler (T/Q),
bifactor (T/Q), simplimax. While the port is in progress GPArotation remains the engine for
not-yet-ported criteria (hybrid dispatch in `rotate_model.R`); once all are native it is
**demoted Imports → Suggests** and kept only as a gated test oracle (no runtime dependency). **`varimax`/`promax` stay
native R** (already not GPArotation; `.VARIMAX_SPSS`/`.SV` retirement deferred). **Serial C++,
no OpenMP this phase** (the C++-vs-R win is the lever; `EFA_AVERAGE` already parallelises its
grid via `future`, so in-engine OpenMP would nest/oversubscribe; `R::rnorm` starts stay serial
and reproducible — thread-safe draws + OpenMP deferred to post-1.0). `infomax` and `target/pst`
are **new** criteria, not parity work — deferred. The engine is **not** on the bootstrap hot
path (oblique bootstrap alignment already uses the P3.5 `.oblique_procrustes_batch`); it
affects the **point estimate** and **`EFA_AVERAGE`** only.

**Validation philosophy (relaxed per maintainer).** Rotation criteria are non-convex, so
"parity" = the C++ engine reaches a criterion value `Q` **as good as or better than**
GPArotation (same `Q` definition) **and** returns a valid rotation; for well-identified
fixtures, aligned loadings match GPArotation within **~1e-4** (`aligned_max_diff`, sign
/permutation removed) and the Phi fingerprint within ~1e-4. Byte-identity to GPArotation's
specific iterates is **not** required (different RNG/line-search may find a different,
equally-valid minimum). `Structure = pattern·Phi` invariant to reflection/ordering.

- [x] `src/rotate.cpp`: generalize the scaffolding into a templated GPForth+GPFoblq engine +
      `vgQ` criterion interface (Bernaards & Jennrich, 2005); add the orthogonal/Stiefel
      projection; reuse line search/multistart/Kaiser/starts; return multistart local-minima
      diagnostics. Returned `Th` must match GPArotation's convention so `.reflect_and_order()`
      reconstructs unchanged (orth: `loadings = L·Th`; obliq: `loadings = L·(Th⁻¹)'`,
      `Phi = Th'Th`).
- [x] Criterion functors: CF (quartimax κ=0, equamax κ=k/(2p)), oblimin/quartimin (γ),
      geomin (T/Q, ε).
- [x] **Reassess `randomStarts` (B12).** Once benchmarked, if 100 starts keep a representative
      `EFA()`/`EFA_AVERAGE` acceptably fast (very likely — 100 C++ starts should cost ≤ today's
      10 R starts, and the screen/triage path means 100 starts ≠ 100 full optimizations), bump
      `EFA()`'s default 10 → 100 and unify with `.rotate_model()`'s 100; else pick a compromise
      (~30). Behaviour change → NEWS + snapshot updates.
- [x] Surface the local-minima count across starts (cheap; the engine already tracks it).
- [x] Demote `GPArotation` Imports → Suggests once the full criterion set is native; gate the
      parity tests with `skip_if_not_installed` so it stays a dev/CI oracle, not a runtime dep
      (no runtime fallback — all criteria are native).

**Exit:** all eleven GPArotation criteria run entirely in native C++; GPArotation is demoted to
Suggests (gated test oracle only) and `R CMD check` is clean without it; varimax/promax stay native.

**Unit queue (P4.1–P4.8 — each atomic, plan-mode first; relaxed-parity validation). Criterion
functors port GPArotation's `vgQ.<crit>` (value+gradient) as the canonical definition so the
best-`Q` comparison is apples-to-apples; replicate the exact control params EFAtools feeds each
GPArotation call:**
- [x] **P4.1** `rotate.cpp` GPF engine (generalize scaffolding + add GPForth) + **CF family**
      (quartimax, equamax). *DoD:* relaxed parity; `test-ROTATE_ORTH.R`/`test-regression-rotations.R`
      updated; P3.5 Procrustes/bootstrap tests stay green.
- [x] **P4.2** **oblimin/quartimin** functor (oblique). *DoD:* relaxed parity; `test-ROTATE_OBLQ.R` (Phi fingerprint).
- [x] **P4.3** **geomin** functor (geominT + geominQ) — *in progress*. *DoD:* relaxed parity, best-`Q` emphasis (multi-minima).
- [x] **P4.4** **bentler** functor (bentlerT + bentlerQ) — port `vgQ.bentler` (invariant criterion).
      *DoD:* relaxed parity; `test-ROTATE_ORTH.R` + `test-ROTATE_OBLQ.R` updated for those two.
- [x] **P4.5** **bifactor** functor (bifactorT + bifactorQ) — port `vgQ.bifactor` (Jennrich-Bentler).
      *DoD:* relaxed parity, best-`Q` emphasis (multi-minima).
- [x] **P4.6** **simplimax** functor (oblique, `k` target near-zeros) — port `vgQ.simplimax`.
      *DoD:* relaxed parity; `test-ROTATE_OBLQ.R` updated.
- [x] **P4.7** Benchmark + `randomStarts` B12 decision (bump 100 + unify if feasible) +
      local-minima reporting + NEWS. *DoD:* benchmark recorded; default set; snapshots/NEWS updated.
- [x] **P4.8** **Drop GPArotation runtime dep**: confirm no `GPArotation::` runtime calls remain;
      DESCRIPTION Imports → Suggests; gate every parity test/example with `skip_if_not_installed`;
      `devtools::check()` clean without it; NEWS. *DoD:* check clean; suite green.

### Phase 5 — Polychoric + ordinal/DWLS → `0.9.x`
**Goal:** ordinal EFA with a bootstrap-fast polychoric matrix and a DWLS estimator.

**Decisions (resolved 2026-06-17).** **Core ordinal only** for v1: `cor_method` gains `"poly"`
and `"tet"` (all columns treated categorical; tetrachoric = the 2-category case of the same
routine). **Polyserial/biserial mixed-type + `"auto"` detection deferred** to a later sub-phase.
Polychoric returns the matrix always; the ACOV is computed at a **gated level**
`acov = c("none","diag","full")` so callers pay only for what they need — `diag` (per-element
asymptotic variances) for DWLS weights + the bootstrap, `full` Γ (cross-pair covariances) for the
Phase-6 robust/sandwich SEs. The **full Γ is implemented and validated now** (per maintainer) but
not wired to a Phase-5 caller until Phase 6 adds `se = "robust"`. **Method and `cor_method` stay
orthogonal** — `method = "DWLS"` and `cor_method = "poly"` are chosen independently (DWLS+poly is
the documented recommended ordinal combo; no auto-switch). DWLS weights = inverse `diag(Γ)`
(lavaan WLSMV convention), threaded as a **per-correlation-element** weight (NOT per-variable).
`pbv` added to `LinkingTo`. Threading: point-estimate polychoric uses OpenMP over pairs;
per-replicate polychoric runs serial within each `future` worker (no nesting with the
process-level bootstrap).

- [ ] `src/polychoric.cpp`: thresholds-once; two-step ρ via fresh Brent **minimiser**;
      BVN CDF via **`pbv`** (GL-5); OpenMP over pairs; returns the matrix, plus the ACOV at the
      requested `acov` level (O11). Documented empty-cell correction arg + zero-variance detection
      (`cli_abort`) + nearest-PD projection (`cli_warn`), integrated with the existing
      `.prepare_cor_input()` smoothing (no double-smoothing).
- [ ] `cor_method` gains `"poly"`/`"tet"` across `EFA`, `KMO`, `BARTLETT`, `N_FACTORS` + the
      retention criteria (branch at the two `stats::cor()` seams: `.prepare_cor_input()` and the
      bootstrap recompute). (`"auto"` + mixed-type deferred.)
- [ ] `DWLS` method in `estimate_model()` — a `DwlsFunctor` subclass of the Phase-3
      `FactorFunctor` (weighted-ULS, per-element weights `= 1/diag(Γ)`); full `WLS` post-1.0 (O11).

**Validation:** BVN CDF vs `pbv`/`mnormt` ~1e-7; matrix vs `polycor` two-step / `lavaan`
`lavCor(ordered=TRUE)` / `psych::polychoric(correct=FALSE,global=FALSE)`; diag ACOV vs `polycor`
std.err; full Γ vs `lavaan` NACOV; DWLS vs `lavaan` (ordered, DWLS); determinism + PD under
bootstrap resampling.
**Exit:** ordinal EFA + DWLS validated; ordinal bootstrap fast on `DOSPERT_raw`; full Γ ready for
Phase 6.

**Unit queue (P5.1–P5.4 — each atomic, plan-mode first):**
- [ ] **P5.1** `polychoric.cpp` matrix (thresholds + Brent ρ + pbv BVN + OpenMP + robustness),
      `acov = "none"`. Add `pbv` to LinkingTo. *DoD:* matrix vs polycor/lavaan/psych + BVN vs pbv.
- [ ] **P5.2** Plumb `cor_method = "poly"/"tet"` through `.prepare_cor_input()` + the bootstrap
      seam + `match.arg` everywhere → ordinal EFA + bootstrap on the existing ML/ULS.
      *DoD:* ordinal EFA end-to-end; bootstrap determinism/PD; speed on `DOSPERT_raw`.
- [ ] **P5.3** Polychoric **ACOV** (`acov = "diag"|"full"`): per-element asymptotic variances +
      the full Γ (Muthén 1984 / Jöreskog 1994). *DoD:* diag vs polycor std.err; full Γ vs lavaan
      NACOV. (Full Γ consumed in Phase 6.)
- [ ] **P5.4** **DWLS** estimator (`DwlsFunctor`; weights `1/diag(Γ)`; `method = "DWLS"`,
      orthogonal to `cor_method`). *DoD:* DWLS loadings vs lavaan (ordered, DWLS); bootstrap path
      recomputes diag weights per replicate (serial within worker).

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
| GPArotation | Imports | **Suggests** (gated test oracle) | demoted in Phase 4 once the full criterion set is native (O4); pin `>=2022.4-1` while in Imports (B12) |
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
