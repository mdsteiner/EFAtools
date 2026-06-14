# EFAtools — pre-Phase-2 correctness & robustness review

> **Status:** review complete, 2026-06-13, against the code at `288d3d9`.
> **Scope:** statistical correctness, code correctness, documentation correctness,
> deduplication/shortening, and test strength/shortening — for every feature implemented
> so far (Phases 0–1 of `REWRITE_PLAN.md`). Vignettes were excluded by request.
> **Nature:** local working document, **not** part of the package (listed in
> `.Rbuildignore`). Companion to `REWRITE_PLAN.md`; the actionable bugs are mirrored into
> that file's §6 register (IDs B20–B54, plus re-opened B1/B9).

This review was produced by fanning out fifteen subsystem reviewers, each of which ran
live R checks against `psych`, `lavaan`, `GPArotation`, `factanal`, `mitml`/`mice`,
`semTools`, and bundled data; every high/critical and most medium findings were then
re-run by an independent adversarial verifier that tried to refute them. Findings below
marked **[verified]** were reproduced numerically (by the reviewer and confirmed by the
verifier); **[read]** were confirmed by code reading; **[convention]** depend on a
maintainer decision rather than being unambiguous bugs. One reviewer (print-method
correctness + package-wide dedup sweep) died on an API socket error before emitting its
report — see *Known coverage gap in this review* at the end.

---

## 0. Headline verdict on the three flagged areas

**EFA_POOLED standard-error aggregation — CORRECT.** This was the maintainer's chief
worry and it holds up. An independent re-implementation of Rubin's rules reproduced the
pooled SEs and Wald CIs for unrotated/rotated loadings and Phi **to machine precision**:
within-imputation `U` is a true variance (`stats::var` of bootstrap replicates, not an
SD), `B` uses the `m−1` denominator, the `(1 + 1/m)` factor is present, and — the hardest
conceptual trap — per-imputation bootstrap replicates are correctly **re-aligned to the
final MI target** before `U` is estimated, rather than reused from the imputation-specific
arrays. The defects in `EFA_POOLED` are *elsewhere* (pooled CAF, D2 ARIV centering,
bootstrap-failure brittleness, Barnard–Rubin df), not in the SE core.

**Hull method — one real bug, otherwise faithful.** The fit computations, the `q=0` null
row, the df axis (provably equivalent to the paper's free-parameter axis), the dominance
elimination, and the st index all match Lorenzo-Seva, Timmerman & Kiers (2011) and are
byte-identical to `EFA()`'s `.gof()`. But the convex-hull elimination loop has a
**critical off-by-one** that changes the retained number of factors on real bundled data
(B20).

**Oblique Procrustes suite — core math correct, two robustness defects.** The C++
objective, gradient (finite-difference verified), manifold projection, retraction, Armijo
line search, and the orthogonal closed form are all exactly right, and with an adequate
iteration budget the solver matches or beats `GPArotation::targetQ`. The two problems are
(a) the default **identity start** with zero random starts drops `PROCRUSTES()` into bad
local minima on well-conditioned problems (B29), and (b) the multi-start **selection rule**
prefers any converged fit over any non-converged one regardless of objective value, so it
can discard a far better solution (B34).

---

## 1. What is solid (do not touch)

These were checked numerically and verified strong — the regression net is real:

- **Estimators.** PAF/ML/ULS reproduce `psych::fa`, `stats::factanal`, SPSS-27 stored
  output, and `psych fm="uls"/"minres"` to between 1e-14 and 1e-6; the `ml_helper.cpp`
  gradient passes finite differences; ULS==minres confirmed.
- **Reworked χ² / fit core.** Bartlett-corrected ML χ² == `factanal` STATISTIC to 3e-9;
  residual-based ULS χ² == `psych fm="minres"` STATISTIC; SRMR == lavaan to 4e-10;
  `.null_chisq` == `psych::cortest.bartlett`; RMSEA point + noncentrality CI match the
  standard convention and lavaan at large N. (CFI/TLI clamping is the exception — B33.)
- **Rotation invariants.** Model-implied R reproduced to ≤3.3e-16 for all 13 rotations;
  SPSS varimax/promax reproduce bundled SPSS-23 output to ≤0.0072; the B9 *Phi/Structure*
  propagation is correct (only the rotmat residual remains — B22).
- **Retention criteria.** EKC == Braeken & van Assen (2017) exactly; MAP TR2 ==
  `psych::VSS` to machine precision (and the "TR4 is wrong" suspicion was **refuted** —
  TR4 is a documented Caron 2019 variant, not the revised MAP); PARALLEL B14/B18 fixes
  verified correct; NEST implements Achim's NESTip; CD faithful to Ruscio & Roche (2012).
- **KMO/BARTLETT** == `psych` exactly. **Preset tables** reproduce psych-type PAF to
  1e-15 with one consolidated override warning. **SL/OMEGA core** (`variance="correlation"`)
  == `psych::schmid`/`psych::omega` to machine precision. **Rubin pooling** (above).
- **Test suite**: 2931 pass / 0 fail / 1 deliberate skip, ~3 min; exemplary oracle-parity
  regression nets; consistent decimal-scrub snapshots; near-uniform class-based condition
  tests; correct Suggests gating.

---

## 2. Critical (wrong numbers or crash in normal use)

### B20 — Hull elimination off-by-one changes retained factors `[verified]`
`R/HULL.R :: .hull_calc` (step-5/6 loop). The loop is `i <- 2; while (i < nr_s - 1)`, so
the last interior triplet `(nr_s-2, nr_s-1, nr_s)` is never tested and with `nr_s == 3`
the body never runs. An interior point can remain strictly below its neighbours' chord and
is still marked `on_hull = TRUE`, which the published algorithm forbids. Both the roxygen
(line ~59) and the in-code comment (line ~392) say "ALL the triplets … are considered",
so this contradicts intent, not an alternative convention.
**Evidence:** `HULL(WJIV_ages_14_19$cormat, N=…, method="ULS")` returns RMSEA `retain = 12`
with hull `{0,1,3,4,5,7,8,12,27,28}`; point 27 (fit 0.88792) lies below `chord(12,28)` =
0.89127 and is the last interior point the loop never tests; a paper-faithful re-run drops
27 and retains 1. Present since the 2020 commit that added HULL, so shipped on CRAN.
**Fix:** extend the loop bound to test the final interior triplet (`while (i <= nr_s - 1)`
with guarded indexing), or restate as a standard upper-convex-hull pass; add a geometry
unit test on a hand-built fit table (current tests only check `on_hull` is logical).

### B21 — Pooled CAF computed on a zero-diagonal residual matrix `[verified]`
`R/EFA_POOLED.R :: .efa_pooled_fit_indices`. `residuals` has `diag(residuals) <- 0`
(line ~310); `.efa_pooled_fit_indices` symmetrizes but never restores the unit diagonal
before calling `.compute_kmo`, whereas the reference `.gof()` sets `diag(delta_hat) <- 1`.
The inline comment at line ~775 even says "with diagonal set to 1" — the code contradicts
its own stated intent. The `tryCatch` only catches hard errors, not the NaN (negative
inverse-diagonal entries raised to `-0.5`) or the garbage near-1 value.
**Evidence:** `bfi[1:800,1:10]`, 10% MCAR, `mice m=3`: ML/promax pooled CAF = **NaN**
(components ≈ 0.503–0.506; diag-1 recompute 0.5052); PAF/promax pooled CAF = **0.9989**
(fixed 0.5015); ML/unrotated pooled CAF = **0.9968** (fixed 0.5053). Wrong in *every*
`EFA_POOLED` call; the bootstrap path shares the helper so it is wrong there too. (This is
the NaN informally flagged in the U1 notes, now pinned and broader than first thought.)
**Fix:** set `diag <- 1` on the residual matrix inside `.efa_pooled_fit_indices` before
`.compute_kmo` (one line; covers both the point and bootstrap paths). Re-baseline the
`EFA_POOLED` CAF snapshots, which currently enshrine the broken value.

---

## 3. High (wrong numbers in realistic cases, or a broken API promise)

### B1 (RE-OPEN) — `rotation = "quartimax"` still runs Bentler's criterion `[verified/read]`
`R/rotate_model.R :: .orth_engines` line 18: `quartimax = function(L, ...) GPArotation::bentlerT(L, ...)`.
The 7af40bf "fix" built the engine-by-name table but the table **preserves the original
misrouting** (comment: "uses the bentlerT criterion, as in the established
implementation"), and `test-regression-rotations.R` pins the wrong engine. So B1 is *not*
fixed — it is restructured. `quartimax` and `bentlerT` are genuinely different criteria
with different optima.
**Evidence:** `EFA(rotation="quartimax")` == `GPArotation::bentlerT` to 1e-5 but differs
from `GPArotation::quartimax` by max 0.84 (mean 0.25, aligned) and from
`psych::fa(rotate="quartimax")` by 0.263; the quartimax criterion value is −3.022 for the
true quartimax solution vs −2.866 for the shipped bentlerT one. Qualitatively different
(general-factor-dominant vs varimax-like column SS).
**Decision needed `[convention]`:** O6 says fix outright with a NEWS entry. Either (a)
point `quartimax` at `GPArotation::quartimax` / `cfT(kappa = 0)` (correct, behaviour-
changing — re-baseline the rotation regression test), or (b) keep the legacy mapping
**but** rename the exposed option honestly (e.g. document `quartimax` as Bentler-invariant
and add a separate true-quartimax name). Current state (named quartimax, silently Bentler,
register says "fixed") is the one option that should not stand.

### B22 — `rotmat` satisfies no reproduction identity after reflect/reorder `[verified/read]`
`R/rotate_model.R :: .reflect_and_order` line 79 (`rotmat <- rotmat[ord, ord]`) and
`R/PROMAX.R :: .rotate_promax` (lines ~32, ~84). Sign reflection (lines 67–69) and column
reordering are applied to the loadings, and Phi is handled consistently (B9 genuinely
fixed there), but the **sign flips are never applied to `rotmat`** and both its rows and
columns are permuted. So `L_unrot %*% rotmat == rot_loadings` fails whenever any factor is
reflected or reordered. `rotmat` has no internal consumer, so nothing compensates; it is an
output-only contract that is silently broken.
**Evidence:** `test_models$baseline` (3-factor PAF): `max|L %*% rotmat − rot_loadings|` =
1.165 (psych varimax), 1.308 (psych promax); all 11 GPArotation rotations under
`type="EFAtools"` break at 0.50–1.10 under both `L%*%rotmat` and `L%*%t(solve(rotmat))`.
**Fix (verified):** `rotmat <- (rotmat %*% diag(signs))[, ord]` restores the identity to
~1e-16; mirror in `.rotate_promax`. Add an `L_unrot %*% rotmat == rot_loadings` invariant
test that forces a reflection+reorder.

### B23 — `FACTOR_SCORES(impute=…)` never imputes; "means"/"median" are no-ops `[verified]`
`R/FACTOR_SCORES.R`. `psych::factor.scores` is called without `missing=`, so it stays at
psych's default `FALSE` and the imputation block is unreachable. Worse, any non-`"none"`
`impute` silently switches psych to the `factorScoresSapa` available-item algorithm; and
psych tests `impute == "mean"` while EFAtools passes `"means"`, so even if reachable it
would hit the median branch. The returned `$missing` is the logical `FALSE`, not the
documented per-subject count.
**Evidence:** `DOSPERT_raw[1:300,1:12]` + 25 NAs: `impute="none"` → 50 NA score cells;
`"means"`/`"median"` → 0, and the two are byte-identical to each other and to psych's
`missing=FALSE` sapa path; differ from true mean-imputation by up to 0.33. The API promise
("Whether and how missing values … should be imputed") is unfulfilled.
**Fix:** pass `missing = (impute != "none")` through to psych and align the option spelling
(`"mean"`/`"median"`), or drop the `impute` argument until `efa_scores` reimplements
scoring natively (Phase 7). Add a numeric imputation test.

### B24 — `.OMEGA_FLEX` discards a user `cormat` for SL inputs; flexible-SL OMEGA crashes `[verified]`
`R/OMEGA_helper.R :: .OMEGA_FLEX` line 26: `cormat <- model$orig_R` runs unconditionally
in the SL branch *before* the `is.null(cormat)` check, but `SL()` sets `orig_R <- NA` for
the documented flexible-input workflow (loadings + Phi, and the lavaan input path). The
OMEGA roxygen says the SL `orig_R` is used "If left NULL", so the overwrite contradicts the
docs and makes the classed `efa_omega_need_cormat` error unreachable for SL objects.
**Evidence:** `OMEGA(sl_flex, cormat = R)` dies with a raw `simpleError` (no `efa_*` class);
the equivalent `model = NULL` manual path with the same loadings + cormat returns the
correct g-row (tot 0.883, hier 0.740, H 0.842).
**Fix:** only take `cormat` from `model$orig_R` when the user did not supply one and it is
non-NA; otherwise honour the user cormat (and keep the classed error when neither exists).

### B25 — `.OMEGA_FLEX variance="sums_load"` mixes two denominators in the g row `[verified]`
`R/OMEGA_helper.R :: .OMEGA_FLEX`. `omega_tot_g` divides by the correspondence-*zeroed*
loading sums while `omega_h_g`/`omega_sub_g` divide by the *unzeroed* sums — two different
total-variance values for the same general-factor composite, so the `tot = hier + sub`
identity (asserted for the other paths) breaks for the g row only.
**Evidence:** g row tot 0.88072, hier 0.73992, sub 0.12495, `tot−(hier+sub)` = 0.01584;
zeroed denominator 98.033 vs unzeroed 99.829; a consistent denominator restores the
identity to 0. Correspondence zeroes 36/54 loadings in the test setup, so it triggers in
every `sums_load` call.
**Fix:** use one denominator (the unzeroed sums, matching the group rows) for all three g
quantities. This folds naturally into the §4.5 `efa_reliability` `.reliability_core` work,
but is a one-line fix now if the `sums_load` numbers are user-facing before then.

### B26 — PUC computed wrongly when its simple-structure/complete-g assumptions fail `[verified]`
`R/OMEGA_helper.R` (`.OMEGA_LAVAAN` / `.OMEGA_FLEX`). The PUC formula assumes a complete
general factor and simple structure; in the two explicitly supported scenarios that
violate this (cross-loadings via `factor_corres`, incomplete g) it returns a value that no
longer equals the proportion of unique covariance explained by g.
**Fix:** compute PUC from the actual count of cross-factor item pairs implied by the
correspondence map rather than the closed-form simple-structure shortcut; document the
definition. (Schedule with the Phase-7 reliability rework; snapshot the corrected value.)

### B27 — `NEST()` crashes with an unclassed C++ error on a Heywood reference model `[verified]`
`R/NEST.R` (line ~111). For `nf ≥ 2` NEST fits `EFA(R, N, n_factors = nf-1, …)` and builds
`u2 <- diag(sqrt(1 - mm$h2))` with **no Heywood handling** — it ignores the centralized
`heywood` field that `EFA()` now exposes. A Heywood `h2 > 1` makes `sqrt(1-h2)` NaN, and
`nest_sym.cpp`'s (correct) eig guard then throws a plain `Rcpp::stop`. `N_FACTORS` happens
to swallow it via `try()` and defaults to ML, but a **direct `NEST()` call — the documented
usage, with the docs recommending PAF as "more robust" — hard-crashes**.
**Evidence:** a 5×5 matrix with `r12=r13=.90` reproduces it: NaN in `u2`, then
`"Eigendecomposition failed during the NEST simulation…"` with classes `Rcpp::exception`,
`C++Error` (no `efa_*` class).
**Fix:** clamp/guard communalities before forming `u2` (consult the `heywood` field) and
raise a classed `efa_*` condition with an actionable message; add a Heywood NEST test.

### B28 — `EFA_POOLED(se="np-boot")` hard-errors if any one bootstrap replicate failed `[verified]`
`R/EFA_POOLED.R :: .efa_pooled_bootstrap_pool` (loop ~1252–1292). The per-replicate loop
assumes every `boot.arrays` slab is fully populated; a single failed replicate in any
component `EFA` (which `EFA` itself tolerates and NA-fills) propagates NA into the Procrustes
re-alignment and aborts the whole pooled bootstrap.
**Fix:** skip/`na.omit` non-finite replicates per imputation (as the point path already
tolerates), tally them, and warn with the existing `efa_pooled_*` class; only abort if too
few valid replicates remain.

### B29 — Oblique `PROCRUSTES()` identity start → silent bad local minima `[verified]`
`src/oblique_procrustes.cpp` (line ~447, `T_primary` = `fill::eye`) + `R/PROCRUSTES.R`
(default `oblique_random_starts = 0`). The exported oblique `PROCRUSTES()` starts from the
identity with no random starts; trapped solutions still report `convergence = TRUE`, so the
only safety net (the convergence flag) does not catch them. This feeds `EFA_POOLED`'s
`first_target` loop (default `procrustes_args = list()`) and the `EFA` oblique bootstrap
(`oblique_random_starts = 5`, often too few). The **consensus** path is unaffected — it
warm-starts every inner alignment from the orthogonal Procrustes solution, which is exactly
the missing fix.
**Evidence:** synthetic self-recovery (global optimum 0): identity start trapped 23/50
(12/15 of them reporting `convergence=TRUE`), 5 random starts 12/50, orthogonal warm start
**0/50**. `EFA_POOLED` `first_target` on `DOSPERT_raw` (6 factors, promax): alignment
values `[0.118, 0.060, 2.136, 2.194]` vs warm-start `[…, 0.241, 0.252]`; max pooled-loading
difference 0.198, all flagged `convergence=TRUE`.
**Fix:** default the oblique `PROCRUSTES` to a warm orthogonal-Procrustes start (and/or a
few random starts) as the consensus engine already does; tighten the bootstrap default.

### B30 — Greedy congruence alignment can collapse two factors onto one in EFA_AVERAGE `[verified]`
`R/averaging.R :: .array_reorder` (lines ~344–363). `factor_order <- apply(abs(congruence), 1, which.max)`
with no permutation check, then `Ln[, factor_order]` and an unguarded
`diag(sign(diag(crossprod(L1, Ln))))`. When two target factors both argmax to the same
source column, that column is duplicated and another dropped — a non-permutation "alignment"
that silently corrupts the average. The correct helper (`.align_solution`, LSAP + zero-sign
guard) already exists and is used by EFA bootstrap, EFA_POOLED, and consensus.
**Evidence:** cross-loaded 3-factor data, `EFA_AVERAGE(…, rotation="oblique")`: on 1 of 8
seeds the simplimax solution aligns with order `(1,3,1)`; the misaligned solution enters the
average (the `admissible` flag is recorded but never used to exclude it);
`max|greedy − LSAP average|` = 0.023.
**Fix:** replace `.array_reorder`'s greedy match with `.align_solution` (LSAP), or reject
non-permutation orders and exclude/flag that solution.

### B31 — Vector `precision` makes every EFA in the grid fail (silent all-NA) `[verified]`
`R/EFA_AVERAGE.R` (future_lapply call ~582; single-combination return ~546). `precision`
is documented and validated as a vector and expanded into the grid by `.type_grid`, but the
worker passes the **raw user argument** `precision = precision` instead of `arg_grid$precision[i]`,
while every other column is recycled per row. `EFA()`'s `assert_number(precision)` then
rejects the length->1 vector; the failure is swallowed by `try(silent=TRUE)` and the whole
grid returns all-NA.
**Evidence:** `precision = c(1e-5, 1e-3)` → both grid rows `errors = TRUE`, error
"Assertion on 'precision' failed: Must have length 1", `all(is.na(loadings))`; scalar
`precision` succeeds. Also a mixed-type grid silently sends every row the same scalar.
**Fix:** pass `precision = arg_grid$precision[[i]]` in both call sites; add a vector-argument
recycling test.

### B32 — D2 χ² ARIV centred at `sqrt(mean(χ²))` instead of `mean(sqrt(χ²))` `[verified]`
`R/EFA_POOLED.R :: .efa_pooled_D2` (line ~647). The average relative increase in variance
centres the `sqrt(χ²)` deviations around `sqrt(mean(χ²))` rather than `mean(sqrt(χ²))`, so
it deviates from Li, Raghunathan & Rubin (1991) and the canonical `mitml:::.D2`
(`r <- (1 + 1/m) * var(sqrt(d))`). The bias is strictly one-sided (Jensen), inflating ARIV
and FMI and shifting the pooled χ²/p, RMSEA, CFI, AIC, BIC for both target and null models.
**Evidence:** `.efa_pooled_D2(c(30,35,50), 20)` → ARIV 0.9064 vs mitml/Li 0.9038, χ²
1.0888 vs 1.1466 (~5%); the gap vanishes as the χ²s converge; 2000/2000 random draws had
implemented ARIV ≥ literature ARIV.
**Fix:** `r <- (1 + 1/m) * stats::var(sqrt(chis))`. Cite Li et al. (1991) in the roxygen
(currently no references). Re-baseline the D2 snapshots.

### B33 — CFI/TLI lack noncentrality clamping → deflated CFI when the null fits well `[verified, convention]`
`R/fit-indices.R :: .gof` (CFI block ~148–162; same pattern in `EFA_POOLED.R` ~751). CFI is
the raw `(δ_null − δ_m)/δ_null` with a post-hoc clamp to [0,1] but **no** `max(δ, 0)`
truncation of the noncentralities, contrary to the cited Bentler (1990)/Kline definition
(which truncates δ at 0). When both model and null over-fit, CFI returns arbitrary low
values including 0, while `p_chi ≈ 0.9` and `RMSEA = 0` for the same fit — internally
contradictory.
**Evidence:** random 7-var data, 1-factor ML, N=100: EFAtools CFI 0.21/0.13/0.00 across
seeds where the Bentler convention and `lavaan::cfa` both give 1.00; `psych::fa` gives NaN.
On well-fitting `test_models$baseline` the corrected formula is byte-identical (0.989235),
so the fix is regime-preserving.
**Fix `[convention]`:** truncate the deltas at 0 before the ratio (`δ_m ← max(δ_m, 0)`,
`δ_null ← max(δ_null, 0)`), per the cited source. Behaviour-changing only in the over-fit
regime; NEWS entry + snapshot update. TLI carries the same unclamped convention — decide
together.

### B34 — Procrustes multi-start selection prefers convergence over objective value `[verified]`
`src/oblique_procrustes.cpp :: fit_is_better_cpp` (lines ~278–293). A converged incumbent
can never be displaced by a non-converged candidate regardless of objective value, and the
returned `convergence` flag is the selected fit's own (TRUE), so the `EFA.R` `isFALSE(convergence)`
filter cannot catch it. Combined with B29 this means the random-start machinery can find
the global optimum and then throw it away.
**Evidence:** exact-fit constructions, 10 starts: 8/60 seeds returned a converged local
minimum (value 1.1–1.5) while fully-optimized starts at `f≈0` were discarded; raising
`oblique_maxit` rescues some cases (so the maxit budget is also a factor).
**Fix:** select on objective value first, breaking ties by convergence (or among converged
fits only when at least one converged); surface the true best-start convergence status.

---

## 4. Medium (robustness, minor numeric convention, or a meaningful crash on edge input)

- **B35 — Heywood detection is unreachable for ML/ULS `[verified, convention]`.**
  `R/estimate_model.R :: .finalize_fit` flags `h2 ≥ 1+eps`, but the ML/ULS L-BFGS-B bound
  `psi ≥ .005` makes `h2 ≥ 1` unreachable, so boundary (improper) ML/ULS solutions are never
  flagged while equivalent PAF ones are. Matches `factanal`/`psych` silence, **but** the
  unreleased NEWS promises detection "consistent across PAF, ML, and ULS", which is vacuous.
  `summary()` prints "Heywood cases: 0" for improper ML solutions and EFA_AVERAGE
  admissibility filtering lets them through. *Fix:* flag boundary solutions (`psi` at the
  lower bound, or `h2 ≥ 1 − tol`) for ML/ULS, or soften the NEWS claim to PAF.
- **B36 — ECV/group-H use correspondence-zeroed group loadings `[verified]`.**
  `R/OMEGA_helper.R :: .OMEGA_FLEX`. Deviates from `psych::omega`/Rodriguez et al. (2016)
  despite the documented psych-reproduction claim. *Fix:* use unzeroed group loadings for
  ECV (fold into the reliability rework).
- **B37 — SL second-order Heywood → NaN columns, swallowed warning `[verified]`.**
  `R/SL.R` EFA branch: a second-order Heywood yields NaN `s*` columns with only a cryptic
  base warning; `suppressWarnings()` also swallows the package's own classed `efa_heywood`.
  *Fix:* detect and raise a classed condition; do not blanket-suppress.
- **B38 — SL lavaan branch lacks the psi/theta Heywood guards `.OMEGA_LAVAAN` has, and uses
  an elementwise `sqrt(psi)` valid only for diagonal psi `[verified]`.** *Fix:* port the
  guard set; use the diagonal explicitly.
- **B39 — NEST no-stop boundary reports `max_fac − 1` (the just-rejected model) `[verified]`;
  reference fits run past the `df>0` bound.** `R/NEST.R`. *Fix:* report the last *accepted*
  model and stop at the df bound.
- **B40 — PARALLEL silently returns NA when real eigenvalues never cross the reference
  `[verified]`.** `R/PARALLEL.R :: .determine_factors`. *Fix:* decide the boundary
  convention (retain max, or 0) and document; classed note.
- **B41 — Smoothed matrices re-flagged NPD on every downstream call `[verified]`** (B15
  incomplete *across* calls). `R/cor-input.R :: .prepare_cor_input`: the NPD threshold
  (`eps^0.6`) is stricter than `psych::cor.smooth`'s, so a smoothed matrix re-trips
  `efa_cor_smoothed` in each chained call (e.g. inside HULL→PARALLEL→EFA). *Fix:* relax the
  re-check threshold or tag already-smoothed matrices.
- **B42 — `N_FACTORS` suitability gate hard-aborts when `N` is missing, even for N-free
  criteria `[verified]`.** *Fix:* gate the abort on whether an N-requiring criterion/index
  was actually requested.
- **B43 — Raw-data path: an NA-containing computed correlation matrix crashes with an
  unclassed base error `[verified]`.** `R/cor-input.R`. *Fix:* detect and raise a classed
  condition pointing at `use=`.
- **B44 — `.hull_calc` crashes with an opaque base error if any gof value is NA `[verified]`.**
  *Fix:* drop/short-circuit NA rows with a classed warning.
- **B45 — CONSENSUS_PROCRUSTES between-run congruence computed on unaligned loadings
  `[verified]`** (lines ~325–334) — the diagnostic is uninterpretable. *Fix:* align before
  computing congruence.
- **B46 — EFA_AVERAGE permanently overwrites the user's global `progressr` handlers
  `[verified]`** (~551–555). *Fix:* set within `withr`/`on.exit` scope.
- **B47 — EFA_AVERAGE relays classed EFA warnings once per grid row (spam); the EFA.R
  comment claims a suppression that does not exist `[verified]`.** *Fix:* collect and emit
  one summarized warning; correct the comment.
- **B48 — EFA_AVERAGE does not extract/tally/average the new fit indices (SRMR, TLI, ECVI,
  RMSR) `[verified]`.** `R/averaging.R :: .extract_data`. *Fix:* add them to the extraction
  and averaging maps.
- **B49 — `print.EFA_AVERAGE` auto-plots by default `[verified]`,** against the package's
  "print never plots" convention. *Fix:* default `plot = FALSE`; keep `plot.EFA_AVERAGE`.
- **B50 — Rubin df without Barnard–Rubin can exceed N; reported FMI is the uncorrected λ
  `[verified]`.** `R/EFA_POOLED.R :: .efa_pooled_rubin_pool` (~979–987). *Fix:* apply the
  Barnard–Rubin small-sample df adjustment (needs a complete-data df per parameter).
- **B51 — NA `chi_null` is unguarded: `det()` underflow in `.null_chisq` crashes EFA mid-
  `.gof` and SMT, and silently NAs BARTLETT `[verified]`.** `R/fit-indices.R`, `R/SMT.R:150`.
  *Fix:* guard the determinant; classed condition.
- **B52 — `factor_corres.cpp` cross-loading index sets collide for ≥10 factors `[verified]`**
  (indices concatenated without a separator, lines ~88–89). *Fix:* join with a separator.
- **B53 — Bootstrap re-fires per-replicate warnings (`efa_type_override`, C++ max-iter)
  `b_boot` times `[verified]`.** `R/EFA.R :: .boot_fun`. *Fix:* suppress within the replicate
  loop; surface once.
- **B54 — Pooled `fit_indices$Fm` is a different estimand than the component `Fm` it is
  named after (≈2× off for ULS); pooled set omits TLI/ECVI `[verified]`.**
  `R/EFA_POOLED.R :: .efa_pooled_fit_indices`. *Fix:* recompute/rename consistently; add the
  missing indices. (Also: `fit_indices_pooled_algorithm` SEs/CIs are statistically
  incoherent — residual indices ~1/√m too narrow, χ²-based inflated — consider dropping that
  block or documenting it as descriptive only.)
- **B55 — Any 1-factor solution with a rotation prints garbage Phi/Structure/Variances
  tables `[verified]`** (from the re-run print audit). `R/rotate_model.R :: .rotate_single_factor`
  returns `Phi`/`Structure`/`vars_accounted_rot` as scalar **logical `NA`** (not `NULL`), and
  `print.EFA`/`summary.EFA` guard those blocks with `is.null` only, so the `NA` scalars render
  as nonsense. *Evidence:* `EFA(test_models$baseline$cormat, 1, N=500, method="PAF", rotation="promax")`
  prints a 1×1 "Factor Intercorrelations" `F1 NA`, a "Variances Accounted for" `V1 NA`, and
  `summary()` dumps the Structure as `[1] NA`. The **variances** breakage hits the orthogonal
  1-factor path too (varimax → `V1 NA`), since both branches NA-fill `vars_accounted_rot`. Not
  a numeric error — loadings/h2/u2/fit print correctly — but visibly broken on a common path.
  *Fix:* return `NULL` (not `NA`) from `.rotate_single_factor` for the inapplicable fields
  (matching the orthogonal Phi/Structure branch), or add `is.null || all(is.na())` guards to
  the three print sections.
- **B56 — `format()` leaks ANSI for `LOADINGS`, `SLLOADINGS`, `N_FACTORS`, `efa_retention`
  `[verified]`.** These four `format.*` methods build their value with `cli::cli_format_method`
  and return it unstripped (`color = TRUE` default, or always-coloured rules/bullets), so when
  colours are on the returned string carries raw ESC bytes — violating the convention that
  `format()` is plain text. `format.EFA`/`format.EFA_POOLED`/`format.summary.EFA` already do it
  right (`cli::ansi_strip(capture.output(print(x)))`). Conditional: in a plain non-interactive
  `Rscript` (`num_ansi_colors()==1`) nothing leaks; it bites when colours are forced on **and**
  the caller captures `format()` to a file/string. *Fix:* wrap the four bodies in
  `cli::ansi_strip`, as the EFA family does.

A handful of **low** items are catalogued in §8.

---

## 5. Deduplication / shortening opportunities

- **`SMT` re-implements `.gof` machinery** (duplicate p-values + a copy-pasted RMSEA
  noncentrality-CI block, `R/SMT.R` ~136–180). Call the shared helper.
- **`KGC` and `SCREE` duplicate the PCA/SMC/EFA eigenvalue block.** Extract one
  `.three_eigen(R, N, …)` helper consumed by both (and reusable by PARALLEL/EKC).
- **SL transformation logic is duplicated** between `SL()` and `.OMEGA_LAVAAN`'s higher-
  order branch, and the two copies have already diverged (B38). Single source.
- **SMC starting values computed in three places** (`.fit_ml`, `.fit_uls`, `.PAF`) with
  inconsistent singular-matrix fallbacks. One `.smc_start()` (also tidies B-list edges).
- **`COMPARE` duplicates alignment** with its own greedy non-permutation reorder (same
  class of bug as B30). Route through `.align_solution`.
- **Dead code:** `.calc_cis` and `.stat_over_list` in `R/helper.R` have **zero usages**
  (confirmed) yet ship generated `.Rd` files; `.calc_cis`' roxygen falsely claims
  `EFA_POOLED` uses it. The `.hull_calc` few-solutions branch copies `st` then
  unconditionally overwrites it with NA. Remove.
- **`lifecycle` is an unused dependency** — `@importFrom lifecycle deprecated` in
  `R/EFAtools-package.R` with no `deprecated()`/`deprecate_*()` call anywhere. Drop from
  Imports until the Phase-8 deprecation work needs it (or leave a TODO; it is harmless but
  flagged by `R CMD check` tooling).
- **B57 — `.efa_loading_row_order` (`print.EFA.R` ~1951) and `.loadings_row_order`
  (`print.LOADINGS.R` ~256) are byte-identical `[verified]`** (confirmed by `identical()` on
  the full deparse; both live — called at `print.EFA.R:832/1679` and `print.LOADINGS.R:220`).
  Keep one, call it from both sites; ~17 lines, no behaviour change.
- **B58 — the small-primary-gap diagnostic is computed twice `[verified]`**:
  `.efa_small_primary_gap_count` (`print.EFA.R` ~1640) and `.efa_small_primary_gap_rows`
  (~1801) share an identical sort/threshold block, differing only in `sum()` vs `which()`.
  `count == length(rows)` by construction, so editing one predicate desyncs the displayed
  count from the listed items. *Fix:* make the count return `length(rows_helper(...))`; ~10 lines.
- **`print.EFA.R` (~2000 lines) and `EFA_POOLED.R` (~1540 lines)** remain the two largest
  files; the broader `.efa_format_matrix` consolidation is a Phase-3 target. The print-method
  audit itself is now complete (see §10).

---

## 6. Test suite — strengths and where to invest

**Strengths:** 2931 pass / 0 fail / 1 deliberate skip in ~3 min; `test-regression-estimators.R`
and `test-regression-rotations.R` are genuine oracle-parity nets (sign/permutation-invariant,
seeded, gated); the decimal-scrub snapshot convention is applied consistently; ~268 condition
tests are essentially all class-based; Suggests gating is correct down to top-level fixtures;
bootstrap/C++-guard failure paths are well covered; EFA_POOLED tests assert real statistical
identities (Rubin identity on identical imputations), not shapes.

**Highest-value coverage gaps (each would have caught a finding above):**
1. **Procrustes suite** — ~1700 source lines (incl. exported `CONSENSUS_PROCRUSTES`,
   `.align_solution`, `.tucker_congruence`, the C++ optimizer) covered by **2 expectations**.
   A self-recovery test (`A = L0 %*% solve(t(T))`, recover `T`/`Phi`) catches B29/B34.
2. **Hull geometry** — no test on the elimination itself; a hand-built fit table catches B20.
3. **EFA_AVERAGE** — tests are structural only: add a numeric averaging check, an alignment
   (non-permutation) test (B30), and a vector-argument recycling test (B31).
4. **Fit-index oracles** — RMSEA (point + CI bounds), CFI, AIC/BIC have no cross-check;
   add lavaan/`semTools` parity (catches B33) and a just-identified `df=0` case.
5. **SL/FACTOR_SCORES/OMEGA** — the SL identity test is tautological; add a `psych::schmid`/
   `psych::omega` numeric oracle, a `sums_load` g-row identity (B25), and an `impute` test
   (B23).
6. **SPSS fixtures** — bundled `SPSS_23`/`SPSS_27` oracle data are **unused**; wire them into
   the varimax/promax parity tests (the SPSS claim is currently untested in the suite).

**Shortening / hygiene:**
- `test-regression-rotations.R` consumes **37% of suite runtime** (64.8s) for 19
  expectations — narrow the start counts / fixtures.
- Expensive stochastic fixtures (CD, PARALLEL) are recomputed identically in three files —
  hoist to a shared seeded helper.
- The ~600 lines of "11-fixture × N-field" assertion walls in `test-EFA_AVERAGE.R` /
  `test-EFA.R` are loopable; `expect_output(str(x), "List of N")` appears ~70× as a weak
  structure assertion; a dead 48-combination `efa_all_np` fixture runs and is never asserted.
- `test-helper.R` (1016 lines) no longer mirrors the post-split source layout — split to
  match `fit-indices.R` / `averaging.R` / `procrustes-consensus.R` / `alignment.R` / etc.
- Copy-pasted singular-data block (13 files) and 11×11 Burt matrix (4 files) → shared
  fixtures. Remove commented-out expectations (`test-EFA.R:315-317`, `test-SMT.R:132`).

---

## 7. Documentation findings (selected; full list verified)

- **NEWS.md omits the largest user-visible changes** of the dev version (retention-object
  redesign, `summary.EFA`, `plot.COMPARE`, lavaan→Suggests, the cli restyle). **README.md**
  shows pre-rewrite χ² values, omits MAP in the `N_FACTORS` output, and its Model-Fit block
  lacks the new indices.
- **`N_FACTORS`/`HULL` document `percent` as "a vector of percentiles"** but `PARALLEL`
  rejects vectors — following the doc breaks the criterion.
- **`EFA_AVERAGE` `@return vars_accounted`** says "Based on the unrotated loadings"; the
  averaged values are rotation-based whenever a rotation is used.
- **EFA `@return` omits the new `heywood` element** (promised in NEWS) and
  `standardized_residuals`; **EFA Details & HULL Details** still say the ML default start is
  `"factanal"` (actual: `"psych"`); `MAP` Details claims an `N` argument it does not have and
  renders a stray `#'`.
- **`KGC` uses `>= 1`** while its own docs and the Kaiser convention say "greater than 1"
  (line 118/131/146 vs the prose). Decide and align.
- Stale cross-references (`@seealso` lists missing MAP/NEST/SCREE; `@family` not used);
  `COMPARE @param plot_red` documents 0.001 vs the actual 0.01 default; `SL @param x`
  references a non-existent `L2`; "Grieder & Steiner … Manuscript in Preparation" and a
  Preacher DOI typo.

---

## 8. Low-severity catalogue (one-liners)

`.is_cormat` lower-bound sign flip misclassifies an exact −1 correlation as raw data;
`print.BARTLETT` prints "χ²(df)=NA, p NA" after announcing nothing was rendered; BARTLETT
accepts absurdly small N (negative multiplier → negative χ², p=1) silently; `.resolve_settings`
override warning loses pluralization after `{.val {type}}`; `standardized_residuals` has a
NaN diagonal (0/0); `just-identified df=0` reports `p_chi=0` (perfect fit looks significant);
`n_factors=0` / `b_boot ∈ {0,1,negative}` reach C++/give opaque errors; `.rotate_model`
returns NULL for an unknown rotation name; `order_type` is a documented no-op everywhere
except promax; `randomStarts` roxygen omits the default and the run-to-run nondeterminism;
`screen_keep = 2` silently caps requested random starts; `normalize=TRUE` also Kaiser-
normalizes the Procrustes target (diverges from `GPArotation::targetQ(normalize=TRUE)`);
`.SV` varimax criterion mixes `|λ|` and `λ⁴` terms vs the SPSS manual it cites; CD compares
comparison-data eigenvalues against a target from different (pairwise vs listwise) data when
NAs are present; HULL's PAF gof-coercion notice uses `cli_alert_info` + manual `col_cyan`
and is pinned by exact message text. Unclassed warnings in the EFA bootstrap path violate
the classed-condition convention.

From the print-method re-run (all `[verified]`): **B59** — `print.EFA_AVERAGE` lists
`admissible` among "varied settings" because the column-removal vector (`print.EFA_AVERAGE.R`
~39–40) strips the other outcome flags but not `admissible` (masked when every solution is
admissible; surfaces when any is inadmissible — *fix:* add `"admissible"` to the vector).
**B60 (cluster):** `.efa_num`/`.numformat` render a value rounding to zero from below as
negative zero (e.g. `-0.00001` → `-.00` / `-0`; base R and psych normalise it — *fix:* drop a
leading `−` when only `.`/`0` remain); `.print_efa_structure_note` (`print.EFA.R` ~1874) is
dead (its only call site is a commented `else` block at ~496); `.efa_num`'s padding branch
returns a `cli ansi_string` object even when unstyled (*fix:* `as.character()` before
returning so it is always plain character). The `print.BARTLETT` malformed "p NA" line the
re-run flagged is already the §8 BARTLETT item above (same root; fix early-returns after the
no-result branch).

---

## 9. Suggested unit ordering (per-change workflow)

Each is atomic and independently committable. Recommended order — correctness first,
behaviour-changing items grouped with NEWS + deliberate snapshot updates:

1. **U-rev-a (critical, no behaviour debate):** B20 Hull off-by-one + geometry test;
   B21 pooled CAF diagonal + re-baseline. — done (commit 20f1315)
2. **U-rev-b (rotation contract):** B22 rotmat sign propagation + invariant test;
   decide & resolve B1 quartimax (NEWS if pointed at true quartimax). — done (commit d7bbee5)
3. **U-rev-c (EFA_AVERAGE):** B31 precision recycling, B30 LSAP alignment, B48 new fit
   indices, B49 print default, B46/B47 progressr/warning hygiene, B59 `admissible` print +
   numeric/alignment tests. — done (pending commit)
4. **U-rev-d (Procrustes robustness):** B29 warm start + B34 selection rule; add the
   self-recovery test file. Re-run after to confirm EFA_POOLED first_target/bootstrap CIs.
   — done (pending commit)
5. **U-rev-e (OMEGA/SL/scores):** B24 cormat, B25 sums_load denominator, B23 impute,
   B36/B37/B38 guards; psych oracle tests. (Coordinate with the Phase-7 reliability rework —
   some fold in.) — done (pending commit)
6. **U-rev-f (robustness guards):** B27 NEST Heywood, B28 EFA_POOLED boot-failure, B42/B43/
   B44/B51 classed guards, B41 smoothing re-flag, B52 factor_corres separator.
7. **U-rev-g (fit conventions, behaviour-changing):** B33 CFI/TLI clamping, B32 D2 ARIV,
   B35 Heywood-flag/NEWS, B54 pooled Fm — each with NEWS + snapshot.
8. **U-rev-h (print methods):** B55 1-factor NA-fill print, B56 `format()` ANSI strip ×4,
   B57/B58 print-helper dedup, B60 cluster (negative-zero, dead `.print_efa_structure_note`,
   `.efa_num` ansi_string) + 1-factor/NA print snapshots.
9. **U-rev-i (dedup + docs):** SMT/.gof, KGC/SCREE eigen helper, SMC start, dead-code
   removal, lifecycle; NEWS/README refresh, percent/start/heywood doc fixes.

---

## 10. Print-method audit — completed (re-run 2026-06-13)

The fifteenth reviewer originally died on an API socket error; it was **re-run with
verification** and is now complete. It constructed every `print.*`/`format.*`/`summary.*`
object and captured output across 1-factor, NULL-Phi, NA-fit-index, all-failed, and
long-name edge cases, plus a fresh dedup sweep excluding the §5 items. Result: **10
findings, 6 verified-and-confirmed (0 refuted), 4 low** — folded above as **B55–B60**
(§4 print correctness, §5 dedup, §8 lows) and mirrored into `REWRITE_PLAN.md` §6.

**Strengths it confirmed:** the shared matrix renderer is robust (wide matrices split into
stacked column blocks with `h2`/`u2` repeated, display-width alignment survives long names
and styled headers, and NA/1×1/single-column/lower-triangle matrices all render cleanly);
the `format.EFA` family enforces the plain-text contract via `ansi_strip` (the pattern the
four leaking `format.*` methods in B56 should adopt); fit/residual/CI sections coerce
non-finite values to NA and guard empty matrices, the just-identified/underidentified cases
print clear warnings and suppress fit indices, and the COMPARE/EFA_AVERAGE plot methods
return ggplot objects without drawing and use `linewidth`.

**No remaining coverage gap.** All fifteen subsystem reviews are complete.
