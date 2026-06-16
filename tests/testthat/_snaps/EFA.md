# print.EFA output is stable (PAF, promax)

    Code
      print(efa_psych)
    Output
      
      EFA performed with type = 'psych', method = 'PAF', and rotation = 'promax'.
      
      -- Rotated Loadings ------------------------------------------------------------
      
             F1     F3     F2    h2    u2
      V1 <num> <num> <num> <num> <num>
      V2 <num> <num> <num> <num> <num>
      V3 <num> <num> <num> <num> <num>
      V4 <num> <num> <num> <num> <num>
      V5 <num> <num> <num> <num> <num>
      V6 <num> <num> <num> <num> <num>
      V7 <num> <num> <num> <num> <num>
      V8 <num> <num> <num> <num> <num>
      V9 <num> <num> <num> <num> <num>
      V10 <num> <num> <num> <num> <num>
      V11 <num> <num> <num> <num> <num>
      V12 <num> <num> <num> <num> <num>
      V13 <num> <num> <num> <num> <num>
      V14 <num> <num> <num> <num> <num>
      V15 <num> <num> <num> <num> <num>
      V16 <num> <num> <num> <num> <num>
      V17 <num> <num> <num> <num> <num>
      V18 <num> <num> <num> <num> <num>
      
      Legend:
        bold = |loading| >= <num>
        grey = below cutoff
        red h2/u2 = Heywood-relevant value
      
      
      -- Factor Intercorrelations ----------------------------------------------------
      
            F1 F2 F3
      F1 <num>
      F2 <num> <num>
      F3 <num> <num> <num>
      
      
      -- Variances Accounted for -----------------------------------------------------
      
                           F1 F3 F2
      SS loadings <num> <num> <num>
      Prop Tot Var <num> <num> <num>
      Cum Prop Tot Var <num> <num> <num>
      Prop Comm Var <num> <num> <num>
      Cum Prop Comm Var <num> <num> <num>
      
      
      -- Model Fit -------------------------------------------------------------------
      
      CAF  : <num>
      RMSR  : <num>
      SRMR  : <num>
      df: 102

# print.EFA output is stable (ML, promax)

    Code
      print(efa_ml_moderate)
    Output
      
      EFA performed with type = 'EFAtools', method = 'ML', and rotation = 'promax'.
      
      -- Rotated Loadings ------------------------------------------------------------
      
            F1    F3    F2    h2    u2
      V1 <num> <num> <num> <num> <num>
      V2 <num> <num> <num> <num> <num>
      V3 <num> <num> <num> <num> <num>
      V4 <num> <num> <num> <num> <num>
      V5 <num> <num> <num> <num> <num>
      V6 <num> <num> <num> <num> <num>
      V7 <num> <num> <num> <num> <num>
      V8 <num> <num> <num> <num> <num>
      V9 <num> <num> <num> <num> <num>
      V10 <num> <num> <num> <num> <num>
      V11 <num> <num> <num> <num> <num>
      V12 <num> <num> <num> <num> <num>
      V13 <num> <num> <num> <num> <num>
      V14 <num> <num> <num> <num> <num>
      V15 <num> <num> <num> <num> <num>
      V16 <num> <num> <num> <num> <num>
      V17 <num> <num> <num> <num> <num>
      V18 <num> <num> <num> <num> <num>
      
      Legend:
        bold = |loading| >= <num>
        grey = below cutoff
        red h2/u2 = Heywood-relevant value
      
      
      -- Factor Intercorrelations ----------------------------------------------------
      
            F1 F2 F3
      F1 <num>
      F2 <num> <num>
      F3 <num> <num> <num>
      
      
      -- Variances Accounted for -----------------------------------------------------
      
                           F1 F3 F2
      SS loadings <num> <num> <num>
      Prop Tot Var <num> <num> <num>
      Cum Prop Tot Var <num> <num> <num>
      Prop Comm Var <num> <num> <num>
      Cum Prop Comm Var <num> <num> <num>
      
      
      -- Model Fit -------------------------------------------------------------------
      
      χ²(102) = <num>, p = <num>
      CFI  : <num>
      TLI  : <num>
      RMSEA [90% CI]  : <num> [ <num>; <num>]
      AIC  : <num>
      BIC  : <num>
      ECVI  : <num>
      CAF  : <num>
      RMSR  : <num>
      SRMR  : <num>

# summary.EFA output is stable (PAF, promax)

    Code
      print(summary(efa_psych))
    Output
      
      EFA performed with type = 'psych', method = 'PAF', and rotation = 'promax'.
      
      -- Model Diagnostics -----------------------------------------------------------
      
      Factors: 3
      Variables: 18
      N: 500
      Heywood cases: 0
      Cross-loading items (|loading| >= <num>): 0
      Items without salient loading (|loading| >= <num>): 0
      Factors with fewer than 3 salient indicators: 0
      Items with primary-loading gap < <num>: 1
      Largest |residual|: <num>
      Factor intercorrelations > <num>: none
      
      -- Rotated Loadings ------------------------------------------------------------
      
             F1     F3     F2    h2    u2
      V1 <num> <num> <num> <num> <num>
      V2 <num> <num> <num> <num> <num>
      V3 <num> <num> <num> <num> <num>
      V4 <num> <num> <num> <num> <num>
      V5 <num> <num> <num> <num> <num>
      V6 <num> <num> <num> <num> <num>
      V7 <num> <num> <num> <num> <num>
      V8 <num> <num> <num> <num> <num>
      V9 <num> <num> <num> <num> <num>
      V10 <num> <num> <num> <num> <num>
      V11 <num> <num> <num> <num> <num>
      V12 <num> <num> <num> <num> <num>
      V13 <num> <num> <num> <num> <num>
      V14 <num> <num> <num> <num> <num>
      V15 <num> <num> <num> <num> <num>
      V16 <num> <num> <num> <num> <num>
      V17 <num> <num> <num> <num> <num>
      V18 <num> <num> <num> <num> <num>
      
      Legend:
        bold = |loading| >= <num>
        grey = below cutoff
        red h2/u2 = Heywood-relevant value
      
      
      -- Factor Intercorrelations ----------------------------------------------------
      
            F1 F2 F3
      F1 <num>
      F2 <num> <num>
      F3 <num> <num> <num>
      
      
      -- Structure Matrix ------------------------------------------------------------
      
            F1 F3 F2
      V1 <num> <num> <num>
      V2 <num> <num> <num>
      V3 <num> <num> <num>
      V4 <num> <num> <num>
      V5 <num> <num> <num>
      V6 <num> <num> <num>
      V7 <num> <num> <num>
      V8 <num> <num> <num>
      V9 <num> <num> <num>
      V10 <num> <num> <num>
      V11 <num> <num> <num>
      V12 <num> <num> <num>
      V13 <num> <num> <num>
      V14 <num> <num> <num>
      V15 <num> <num> <num>
      V16 <num> <num> <num>
      V17 <num> <num> <num>
      V18 <num> <num> <num>
      
      
      -- Simple Structure Diagnostics ------------------------------------------------
      
      Items with primary-loading gap < <num>:
      * V11: F2 = <num>, F3 = <num>
      
      
      -- Variances Accounted for -----------------------------------------------------
      
                           F1 F3 F2
      SS loadings <num> <num> <num>
      Prop Tot Var <num> <num> <num>
      Cum Prop Tot Var <num> <num> <num>
      Prop Comm Var <num> <num> <num>
      Cum Prop Comm Var <num> <num> <num>
      
      
      -- Model Fit -------------------------------------------------------------------
      
      CAF  : <num>
      RMSR  : <num>
      SRMR  : <num>
      df: 102
      
      -- Residual Diagnostics --------------------------------------------------------
      
      Residual cutoff: |r| > <num>
      Number of large residuals: 0
      Largest absolute residual: <num>
      
      No absolute residuals > <num> occurred.
      
      Inspect the residual matrix for details (e.g., with residuals()).

# summary.EFA output is stable (ML, promax)

    Code
      print(summary(efa_ml_moderate))
    Output
      
      EFA performed with type = 'EFAtools', method = 'ML', and rotation = 'promax'.
      
      -- Model Diagnostics -----------------------------------------------------------
      
      Factors: 3
      Variables: 18
      N: 500
      Heywood cases: 0
      Cross-loading items (|loading| >= <num>): 0
      Items without salient loading (|loading| >= <num>): 0
      Factors with fewer than 3 salient indicators: 0
      Items with primary-loading gap < <num>: 0
      Largest |residual|: <num>
      Factor intercorrelations > <num>: none
      
      -- Rotated Loadings ------------------------------------------------------------
      
            F1    F3    F2    h2    u2
      V1 <num> <num> <num> <num> <num>
      V2 <num> <num> <num> <num> <num>
      V3 <num> <num> <num> <num> <num>
      V4 <num> <num> <num> <num> <num>
      V5 <num> <num> <num> <num> <num>
      V6 <num> <num> <num> <num> <num>
      V7 <num> <num> <num> <num> <num>
      V8 <num> <num> <num> <num> <num>
      V9 <num> <num> <num> <num> <num>
      V10 <num> <num> <num> <num> <num>
      V11 <num> <num> <num> <num> <num>
      V12 <num> <num> <num> <num> <num>
      V13 <num> <num> <num> <num> <num>
      V14 <num> <num> <num> <num> <num>
      V15 <num> <num> <num> <num> <num>
      V16 <num> <num> <num> <num> <num>
      V17 <num> <num> <num> <num> <num>
      V18 <num> <num> <num> <num> <num>
      
      Legend:
        bold = |loading| >= <num>
        grey = below cutoff
        red h2/u2 = Heywood-relevant value
      
      
      -- Factor Intercorrelations ----------------------------------------------------
      
            F1 F2 F3
      F1 <num>
      F2 <num> <num>
      F3 <num> <num> <num>
      
      
      -- Structure Matrix ------------------------------------------------------------
      
            F1 F3 F2
      V1 <num> <num> <num>
      V2 <num> <num> <num>
      V3 <num> <num> <num>
      V4 <num> <num> <num>
      V5 <num> <num> <num>
      V6 <num> <num> <num>
      V7 <num> <num> <num>
      V8 <num> <num> <num>
      V9 <num> <num> <num>
      V10 <num> <num> <num>
      V11 <num> <num> <num>
      V12 <num> <num> <num>
      V13 <num> <num> <num>
      V14 <num> <num> <num>
      V15 <num> <num> <num>
      V16 <num> <num> <num>
      V17 <num> <num> <num>
      V18 <num> <num> <num>
      
      
      -- Variances Accounted for -----------------------------------------------------
      
                           F1 F3 F2
      SS loadings <num> <num> <num>
      Prop Tot Var <num> <num> <num>
      Cum Prop Tot Var <num> <num> <num>
      Prop Comm Var <num> <num> <num>
      Cum Prop Comm Var <num> <num> <num>
      
      
      -- Model Fit -------------------------------------------------------------------
      
      χ²(102) = <num>, p = <num>
      CFI  : <num>
      TLI  : <num>
      RMSEA [90% CI]  : <num> [ <num>; <num>]
      AIC  : <num>
      BIC  : <num>
      ECVI  : <num>
      CAF  : <num>
      RMSR  : <num>
      SRMR  : <num>
      
      -- Residual Diagnostics --------------------------------------------------------
      
      Residual cutoff: |r| > <num>
      Number of large residuals: 0
      Largest absolute residual: <num>
      
      No absolute residuals > <num> occurred.
      
      Inspect the residual matrix for details (e.g., with residuals()).

# print/summary.EFA omit the inapplicable tables for a rotated single factor

    Code
      print(efa_1fac)
    Output
      
      EFA performed with type = 'EFAtools', method = 'PAF', and rotation = 'promax'.
      
      -- Rotated Loadings ------------------------------------------------------------
      
            F1    h2    u2
      V1 <num> <num> <num>
      V2 <num> <num> <num>
      V3 <num> <num> <num>
      V4 <num> <num> <num>
      V5 <num> <num> <num>
      V6 <num> <num> <num>
      V7 <num> <num> <num>
      V8 <num> <num> <num>
      V9 <num> <num> <num>
      V10 <num> <num> <num>
      V11 <num> <num> <num>
      V12 <num> <num> <num>
      V13 <num> <num> <num>
      V14 <num> <num> <num>
      V15 <num> <num> <num>
      V16 <num> <num> <num>
      V17 <num> <num> <num>
      V18 <num> <num> <num>
      
      Legend:
        bold = |loading| >= <num>
        grey = below cutoff
        red h2/u2 = Heywood-relevant value
      
      
      -- Variances Accounted for -----------------------------------------------------
      
      No variance-accounted table available.
      
      -- Model Fit -------------------------------------------------------------------
      
      CAF  : <num>
      RMSR  : <num>
      SRMR  : <num>
      df: 135

---

    Code
      print(summary(efa_1fac))
    Output
      
      EFA performed with type = 'EFAtools', method = 'PAF', and rotation = 'promax'.
      
      -- Model Diagnostics -----------------------------------------------------------
      
      Factors: 1
      Variables: 18
      N: 500
      Heywood cases: 0
      Cross-loading items (|loading| >= <num>): 0
      Items without salient loading (|loading| >= <num>): 0
      Factors with fewer than 3 salient indicators: 0
      Items with primary-loading gap < <num>: 0
      Largest |residual|: <num>
      
      -- Rotated Loadings ------------------------------------------------------------
      
            F1    h2    u2
      V1 <num> <num> <num>
      V2 <num> <num> <num>
      V3 <num> <num> <num>
      V4 <num> <num> <num>
      V5 <num> <num> <num>
      V6 <num> <num> <num>
      V7 <num> <num> <num>
      V8 <num> <num> <num>
      V9 <num> <num> <num>
      V10 <num> <num> <num>
      V11 <num> <num> <num>
      V12 <num> <num> <num>
      V13 <num> <num> <num>
      V14 <num> <num> <num>
      V15 <num> <num> <num>
      V16 <num> <num> <num>
      V17 <num> <num> <num>
      V18 <num> <num> <num>
      
      Legend:
        bold = |loading| >= <num>
        grey = below cutoff
        red h2/u2 = Heywood-relevant value
      
      
      -- Variances Accounted for -----------------------------------------------------
      
      No variance-accounted table available.
      
      -- Model Fit -------------------------------------------------------------------
      
      CAF  : <num>
      RMSR  : <num>
      SRMR  : <num>
      df: 135
      
      -- Residual Diagnostics --------------------------------------------------------
      
      Residual cutoff: |r| > <num>
      Number of large residuals: 9
      Largest absolute residual: <num>
      
      Largest residuals:
      * V10 ~~ V12: <num>
      * V1 ~~ V6: <num>
      * V9 ~~ V10: <num>
      * V14 ~~ V17: <num>
      * V8 ~~ V10: <num>
      * V1 ~~ V4: <num>
      * V13 ~~ V17: <num>
      * V17 ~~ V18: <num>
      * V4 ~~ V5: <num>
      
      Inspect the residual matrix for details (e.g., with residuals()).

