# EFAtools — working guide for Claude

EFAtools is an R package for exploratory factor analysis (EFA), currently undergoing a
major rewrite and expansion. Work in small, reviewable increments.

## Working agreement (hard rules)

- **Never commit.** Implement only; the maintainer reviews and commits every change.
- **Never bump the package version** (`DESCRIPTION` `Version`, `cran-comments`, etc.).
- **Never create a branch without explicit approval.** Ask first; create one only when
  told to.
- **Implement exactly what was asked** — no extra features, refactors, renames, or
  "drive-by" fixes. If you notice something else worth doing, mention it (or flag a
  separate task); do not act on it in the same change.
- **No bloat.** Prefer the smallest change that does the job. Reuse existing helpers; do
  not add new dependencies, files, or abstractions without asking.
- **Self-contained code & docs.** Never mention a "phase", "step", a roadmap, milestone,
  or any planning document in code, comments, `NEWS.md`, vignettes, documentation, tests,
  or (if you draft them) commit messages. Every shipped artifact must read as if the
  package was always written this way. Other developers and users have no knowledge of the
  rewrite plan.
- **Document clearly.** Exported functions get roxygen2 (markdown) with `@param`,
  `@returns`, `@examples`, and `@family`. Internal helpers get a one-line purpose comment.
  Explain non-obvious math/algorithms in a comment with a literature reference.

## Per-change workflow

1. The maintainer names **one** unit of work (a single fix or feature).
2. Claude enters **plan mode** and proposes a focused plan for that one unit only.
3. On approval, Claude implements it — including tests and any snapshot updates — then
   runs the relevant tests and reports the result.
4. The maintainer runs `/code-review` and reviews personally.
5. The maintainer requests changes (Claude applies them) or edits directly.
6. The maintainer commits. Then move to the next unit.

Keep each unit **atomic** so it can be reviewed and committed on its own. Do not bundle
unrelated changes.

## Conventions

- Messages/warnings/errors: `cli::cli_abort` / `cli_warn` / `cli_inform` (give errors a
  `class=`). No base `stop`/`warning`/`message`, no `crayon`.
- Printing: `format.<class>()` built with `cli::cli_format_method()`; `print` =
  `cat(format(x, ...), sep = "\n")`. `format()` returns plain text (no ANSI).
- Plots: ggplot2 only; return ggplot objects. No base graphics; `print` must not auto-plot.
- Public API: lowercase `efa_*`; internal helpers are dotted (`.helper`). snake_case args.
- Tests: testthat 3rd edition. `expect_snapshot()` for printed output (wrap in
  `testthat::local_reproducible_output()`); `vdiffr` for plots; cross-check numerics
  against reference packages (psych, lavaan, GPArotation, …) with explicit `tolerance=`,
  gated by `skip_if_not_installed()` / `skip_on_cran()`. Condition tests assert classes
  (`expect_error(..., class = "efa_*")`, same for warnings), never message text via
  regexp.
- C++: RcppArmadillo; pass large matrices by `const&`; keep hot/bootstrap paths
  allocation-light.
- Minimum R is **4.1.0**: the native pipe `|>` and `\(x)` lambdas are fine; do **not** use
  newer features (e.g. the `_` pipe placeholder needs 4.2; base `%||%` needs 4.4 — define
  `%||%` internally if you need it).

## Checks

- Local R: **R 4.6.0** is installed and on `PATH` (multiple versions exist); use it for all
  R / `devtools` work (`Rscript`, `R`).
- Tests: `devtools::test()`.
- Re-document after roxygen edits: `devtools::document()`.
- Full check: `devtools::check()`.

## Roadmap (local only)

The detailed roadmap and design decisions live in `REWRITE_PLAN.md`, a **local working
document that is not part of the package**. Do not reference it from any shipped file. It
may be untracked or absent in some checkouts — if you need it and it is missing, ask.
