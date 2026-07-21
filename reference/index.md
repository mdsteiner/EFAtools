# Package index

## Fitting

Fit exploratory factor models, average across solutions, pool across
imputed data sets, compare groups, and configure the estimation and
rotation engines.

- [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  : Exploratory factor analysis (EFA)
- [`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md)
  : Model averaging across different EFA estimators and types
- [`efa_group()`](https://mdsteiner.github.io/EFAtools/reference/efa_group.md)
  : Multigroup exploratory factor analysis
- [`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md)
  : Exploratory factor analysis on multiple data imputations
- [`plot(`*`<efa_group>`*`)`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_group.md)
  : Plot a multigroup factor analysis
- [`print(`*`<efa_group>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_group.md)
  [`format(`*`<efa_group>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_group.md)
  : Print and format a multigroup factor analysis
- [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  : Control objects for estimation and rotation settings
- [`print(`*`<efa_estimate_control>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_control.md)
  [`format(`*`<efa_estimate_control>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_control.md)
  [`print(`*`<efa_rotate_control>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_control.md)
  [`format(`*`<efa_rotate_control>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_control.md)
  : Print and format a control object
- [`print(`*`<efa>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa.md)
  [`print(`*`<efa_mi>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa.md)
  [`format(`*`<efa>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa.md)
  [`format(`*`<efa_mi>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa.md)
  [`summary(`*`<efa>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa.md)
  [`summary(`*`<efa_mi>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa.md)
  [`print(`*`<summary.efa>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa.md)
  [`format(`*`<summary.efa>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa.md)
  : Print and summarise an efa object
- [`print(`*`<efa_average>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_average.md)
  [`format(`*`<efa_average>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_average.md)
  : Print and format an efa_average object
- [`plot(`*`<efa_average>`*`)`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_average.md)
  : Plot efa_average object
- [`print(`*`<efa_loadings>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_loadings.md)
  [`format(`*`<efa_loadings>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_loadings.md)
  : Print a loading matrix
- [`residuals(`*`<efa>`*`)`](https://mdsteiner.github.io/EFAtools/reference/residuals.efa.md)
  : Extract residuals from an efa object

## Factor retention

Criteria for determining the number of factors to retain, and a wrapper
to run several of them at once.

- [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md)
  : Comparison data
- [`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md)
  : Empirical Kaiser criterion
- [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)
  : Hull method for determining the number of factors to retain
- [`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md)
  : Kaiser-Guttman criterion
- [`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md)
  : Velicer's minimum average partial (MAP) criterion
- [`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md)
  : Next eigenvalue sufficiency test (NEST)
- [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  : Parallel analysis
- [`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
  : Various factor retention criteria
- [`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md)
  : Scree plot
- [`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md)
  : Sequential chi square model tests, RMSEA lower bound, and AIC
- [`print(`*`<efa_retain>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_retain.md)
  : Print method for efa_retain objects
- [`format(`*`<efa_retain>`*`)`](https://mdsteiner.github.io/EFAtools/reference/format.efa_retain.md)
  : Format method for efa_retain objects
- [`plot(`*`<efa_retain>`*`)`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_retain.md)
  : Plot method for efa_retain objects
- [`print(`*`<efa_retention>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_retention.md)
  : Print method for efa_retention objects
- [`format(`*`<efa_retention>`*`)`](https://mdsteiner.github.io/EFAtools/reference/format.efa_retention.md)
  : Format method for efa_retention objects
- [`plot(`*`<efa_retention>`*`)`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_retention.md)
  : Plot method for efa_retention objects

## Rotation and transformation

Align a solution with a target and transform an oblique solution into a
hierarchical one.

- [`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md)
  : Rotate a loading matrix to a target using Procrustes alignment
- [`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)
  : Schmid-Leiman transformation
- [`print(`*`<efa_schmid_leiman>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_schmid_leiman.md)
  [`format(`*`<efa_schmid_leiman>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_schmid_leiman.md)
  : Print and format an efa_schmid_leiman object
- [`print(`*`<efa_sl_loadings>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_sl_loadings.md)
  [`format(`*`<efa_sl_loadings>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_sl_loadings.md)
  : Print an efa_sl_loadings object

## Reliability

Reliability coefficients for a fitted or user-specified factor solution.

- [`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md)
  : Reliability and common-variance coefficients for a factor solution
- [`print(`*`<efa_reliability>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_reliability.md)
  [`format(`*`<efa_reliability>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_reliability.md)
  : Print and format a reliability object

## Scores and comparison

Factor scores with score-quality diagnostics, and comparison of two
solutions.

- [`efa_scores()`](https://mdsteiner.github.io/EFAtools/reference/efa_scores.md)
  : Estimate factor scores and score-quality diagnostics for an EFA
  model
- [`print(`*`<efa_scores>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_scores.md)
  [`format(`*`<efa_scores>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_scores.md)
  [`summary(`*`<efa_scores>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_scores.md)
  [`print(`*`<summary.efa_scores>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_scores.md)
  [`format(`*`<summary.efa_scores>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_scores.md)
  : Print and format an efa_scores object
- [`efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/efa_compare.md)
  : Compare two vectors or matrices (communalities or loadings)
- [`print(`*`<efa_compare>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_compare.md)
  [`format(`*`<efa_compare>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_compare.md)
  : Print and format an efa_compare object
- [`plot(`*`<efa_compare>`*`)`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_compare.md)
  : Plot efa_compare object

## Screening and simulation

Check whether data are suitable for factor analysis and simulate data
from a common factor model.

- [`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md)
  : Bartlett's test of sphericity
- [`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md)
  : Kaiser-Meyer-Olkin criterion
- [`efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/efa_screen.md)
  : Screen data for exploratory factor analysis
- [`print(`*`<efa_screen>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_screen.md)
  [`format(`*`<efa_screen>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_screen.md)
  : Print and format an efa_screen object
- [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md)
  : Simulate data from a common-factor population model
- [`print(`*`<efa_simulated>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_simulated.md)
  [`format(`*`<efa_simulated>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_simulated.md)
  : Print and format an efa_simulated object
- [`print(`*`<efa_bartlett>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_bartlett.md)
  [`format(`*`<efa_bartlett>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_bartlett.md)
  : Print and format an efa_bartlett object
- [`print(`*`<efa_kmo>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_kmo.md)
  [`format(`*`<efa_kmo>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_kmo.md)
  : Print and format an efa_kmo object

## Power

Analytic and simulation-based power analysis.

- [`efa_power()`](https://mdsteiner.github.io/EFAtools/reference/efa_power.md)
  : Power analysis for exploratory factor analysis
- [`plot(`*`<efa_power>`*`)`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_power.md)
  : Plot the RMSEA power curve
- [`print(`*`<efa_power>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_power.md)
  [`format(`*`<efa_power>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.efa_power.md)
  : Print and format an efa_power object

## Datasets

Example data sets and models.

- [`DOSPERT`](https://mdsteiner.github.io/EFAtools/reference/DOSPERT.md)
  : DOSPERT
- [`DOSPERT_raw`](https://mdsteiner.github.io/EFAtools/reference/DOSPERT_raw.md)
  : DOSPERT_raw
- [`GRiPS_raw`](https://mdsteiner.github.io/EFAtools/reference/GRiPS_raw.md)
  : GRiPS_raw
- [`IDS2_R`](https://mdsteiner.github.io/EFAtools/reference/IDS2_R.md) :
  Intelligence subtests from the Intelligence and Development Scales–2
- [`RiskDimensions`](https://mdsteiner.github.io/EFAtools/reference/RiskDimensions.md)
  : RiskDimensions
- [`SPSS_23`](https://mdsteiner.github.io/EFAtools/reference/SPSS_23.md)
  : Various outputs from SPSS (version 23) FACTOR
- [`SPSS_27`](https://mdsteiner.github.io/EFAtools/reference/SPSS_27.md)
  : Various outputs from SPSS (version 27) FACTOR
- [`UPPS_raw`](https://mdsteiner.github.io/EFAtools/reference/UPPS_raw.md)
  : UPPS_raw
- [`WJIV_ages_14_19`](https://mdsteiner.github.io/EFAtools/reference/WJIV_ages_14_19.md)
  : Woodcock Johnson IV: ages 14 to 19
- [`WJIV_ages_20_39`](https://mdsteiner.github.io/EFAtools/reference/WJIV_ages_20_39.md)
  : Woodcock Johnson IV: ages 20 to 39
- [`WJIV_ages_3_5`](https://mdsteiner.github.io/EFAtools/reference/WJIV_ages_3_5.md)
  : Woodcock Johnson IV: ages 3 to 5
- [`WJIV_ages_40_90`](https://mdsteiner.github.io/EFAtools/reference/WJIV_ages_40_90.md)
  : Woodcock Johnson IV: ages 40 to 90 plus
- [`WJIV_ages_6_8`](https://mdsteiner.github.io/EFAtools/reference/WJIV_ages_6_8.md)
  : Woodcock Johnson IV: ages 6 to 8
- [`WJIV_ages_9_13`](https://mdsteiner.github.io/EFAtools/reference/WJIV_ages_9_13.md)
  : Woodcock Johnson IV: ages 9 to 13
- [`population_models`](https://mdsteiner.github.io/EFAtools/reference/population_models.md)
  : population_models
- [`test_models`](https://mdsteiner.github.io/EFAtools/reference/test_models.md)
  : Four test models used in Grieder and Steiner (2022)

## Superseded

The original uppercase names. They remain exported and keep their
original argument lists, but the lowercase `efa_*` equivalents are the
recommended interface.

- [`EFA()`](https://mdsteiner.github.io/EFAtools/reference/EFA.md)
  **\[superseded\]** : Exploratory factor analysis (EFA)
- [`EFA_POOLED()`](https://mdsteiner.github.io/EFAtools/reference/EFA_POOLED.md)
  **\[superseded\]** : Exploratory factor analysis on multiple data
  imputations
- [`EFA_AVERAGE()`](https://mdsteiner.github.io/EFAtools/reference/EFA_AVERAGE-superseded.md)
  **\[superseded\]** : Model averaging across different EFA methods and
  types
- [`N_FACTORS()`](https://mdsteiner.github.io/EFAtools/reference/N_FACTORS.md)
  **\[superseded\]** : Various factor retention criteria
- [`BARTLETT()`](https://mdsteiner.github.io/EFAtools/reference/BARTLETT.md)
  **\[superseded\]** : Bartlett's test of sphericity
- [`KMO()`](https://mdsteiner.github.io/EFAtools/reference/KMO.md)
  **\[superseded\]** : Kaiser-Meyer-Olkin criterion
- [`CD()`](https://mdsteiner.github.io/EFAtools/reference/CD.md)
  **\[superseded\]** : Comparison data
- [`EKC()`](https://mdsteiner.github.io/EFAtools/reference/EKC.md)
  **\[superseded\]** : Empirical Kaiser criterion
- [`HULL()`](https://mdsteiner.github.io/EFAtools/reference/HULL.md)
  **\[superseded\]** : Hull method
- [`KGC()`](https://mdsteiner.github.io/EFAtools/reference/KGC.md)
  **\[superseded\]** : Kaiser-Guttman criterion
- [`MAP()`](https://mdsteiner.github.io/EFAtools/reference/MAP.md)
  **\[superseded\]** : Minimum average partial
- [`NEST()`](https://mdsteiner.github.io/EFAtools/reference/NEST.md)
  **\[superseded\]** : Next eigenvalue sufficiency test
- [`PARALLEL()`](https://mdsteiner.github.io/EFAtools/reference/PARALLEL.md)
  **\[superseded\]** : Parallel analysis
- [`SCREE()`](https://mdsteiner.github.io/EFAtools/reference/SCREE.md)
  **\[superseded\]** : Scree plot
- [`SMT()`](https://mdsteiner.github.io/EFAtools/reference/SMT.md)
  **\[superseded\]** : Sequential model tests
- [`SL()`](https://mdsteiner.github.io/EFAtools/reference/SL.md)
  **\[superseded\]** : Schmid-Leiman transformation
- [`OMEGA()`](https://mdsteiner.github.io/EFAtools/reference/OMEGA.md)
  **\[superseded\]** : McDonald's omega
- [`print(`*`<OMEGA>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.OMEGA.md)
  [`format(`*`<OMEGA>`*`)`](https://mdsteiner.github.io/EFAtools/reference/print.OMEGA.md)
  : Print and format an OMEGA object
- [`FACTOR_SCORES()`](https://mdsteiner.github.io/EFAtools/reference/FACTOR_SCORES.md)
  **\[superseded\]** : Estimate factor scores for an EFA model
- [`COMPARE()`](https://mdsteiner.github.io/EFAtools/reference/COMPARE.md)
  **\[superseded\]** : Compare two vectors or matrices (communalities or
  loadings)
- [`PROCRUSTES()`](https://mdsteiner.github.io/EFAtools/reference/PROCRUSTES.md)
  **\[superseded\]** : Rotate a loading matrix to a target using
  Procrustes alignment
