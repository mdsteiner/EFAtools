# bootstrap arrays are pooled into MI SEs and CIs

    Code
      print(summary(pooled_boot))
    Output
      
      Pooled EFA across 2 imputations performed with type = 'EFAtools', method = 'ML', and rotation = 'none'.
      Pooling settings: align_unrotated = 'signed_tucker_congruence', fit_pool_method = 'D2'.
      
      -- Pooled Model Diagnostics ----------------------------------------------------
      
      Factors: 1
      Variables: 8
      N: 250
      Imputations: 2
      Pooling: align_unrotated = 'signed_tucker_congruence', fit_pool_method = 'D2'
      Bootstrap samples per imputation: 16
      Valid target-rotated samples: 16 out of 16 per imputation
      Heywood cases: 0
      Cross-loading items (|loading| >= <num>): 0
      Items without salient loading (|loading| >= <num>): 0
      Factors with fewer than 3 salient indicators: 0
      Items with primary-loading gap < <num>: 0
      Largest |residual|: <num>
      
      -- Unrotated Loadings ----------------------------------------------------------
      
                  F1    h2    u2
      fun <num> <num> <num>
      friends <num> <num> <num>
      enjoy <num> <num> <num>
      hurt <num> <num> <num>
      part <num> <num> <num>
      commonly <num> <num> <num>
      chances <num> <num> <num>
      attracted <num> <num> <num>
      
      Legend:
        bold = |loading| >= <num>
        grey = below cutoff
        red h2/u2 = Heywood-relevant value
      
      -- 95% bootstrap/MI CIs for salient unrotated loadings -------------------------
      
      Variable   Factor  est    lower  upper
      fun        F1 <num> <num> <num>
      friends    F1 <num> <num> <num>
      enjoy      F1 <num> <num> <num>
      hurt       F1 <num> <num> <num>
      part       F1 <num> <num> <num>
      commonly   F1 <num> <num> <num>
      chances    F1 <num> <num> <num>
      attracted  F1 <num> <num> <num>
      
      -- Variances Accounted for -----------------------------------------------------
      
                      F1
      SS loadings <num>
      Prop Tot Var <num>
      
      -- Model Fit -------------------------------------------------------------------
      
      χ²(20) = <num>, p = <num>
      CFI [95% bootstrap/MI-CI]: <num> [ <num>, <num>]
      TLI: <num>
      RMSEA [90% CI] [95% bootstrap/MI-CI]: <num> [ <num>; <num>] [ <num>, <num>]
      AIC [95% bootstrap/MI-CI]: <num> [ <num>, <num>]
      BIC [95% bootstrap/MI-CI]: <num> [ <num>, <num>]
      ECVI: <num>
      CAF [95% bootstrap/MI-CI]: <num> [ <num>, <num>]
      RMSR [95% bootstrap/MI-CI]: <num> [ <num>, <num>]
      SRMR: <num>
      
      Note: Bootstrap/MI CIs based on 16 bootstrap samples per imputation.
      
      -- MI Uncertainty Summary ------------------------------------------------------
      
      Largest available FMI: <num>
      Median available FMI: <num>
      Largest available RIV: <num>
      Median available RIV: <num>
      
      -- Residual Diagnostics --------------------------------------------------------
      
      Residual cutoff: |r| > <num>
      Number of large residuals: 0
      Largest absolute residual: <num>
      
      No absolute residuals > <num> occurred.
      
      Inspect the residual matrix for details (e.g., with residuals()).

# print.EFA_POOLED output is stable (PAF, promax)

    Code
      print(pooled_obl)
    Output
      
      Pooled EFA across 3 imputations performed with type = 'EFAtools', method = 'PAF', and rotation = 'promax'.
      Pooling settings: target_method = 'first_target', align_unrotated = 'signed_tucker_congruence', fit_pool_method = 'D2'.
      
      -- Rotated Loadings ------------------------------------------------------------
      
             F1     F2     F3    h2    u2
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
      
                           F1 F2 F3
      SS loadings <num> <num> <num>
      Prop Tot Var <num> <num> <num>
      Cum Prop Tot Var <num> <num> <num>
      Prop Comm Var <num> <num> <num>
      Cum Prop Comm Var <num> <num> <num>
      
      -- Model Fit -------------------------------------------------------------------
      
      CAF: <num>
      RMSR: <num>
      SRMR: <num>
      df: 102

# print.EFA_POOLED output is stable (ML, unrotated)

    Code
      print(pooled_none)
    Output
      
      Pooled EFA across 3 imputations performed with type = 'EFAtools', method = 'ML', and rotation = 'none'.
      Pooling settings: align_unrotated = 'signed_tucker_congruence', fit_pool_method = 'D2'.
      
      -- Unrotated Loadings ----------------------------------------------------------
      
                  F1    h2    u2
      fun <num> <num> <num>
      friends <num> <num> <num>
      enjoy <num> <num> <num>
      hurt <num> <num> <num>
      part <num> <num> <num>
      commonly <num> <num> <num>
      chances <num> <num> <num>
      attracted <num> <num> <num>
      
      Legend:
        bold = |loading| >= <num>
        grey = below cutoff
        red h2/u2 = Heywood-relevant value
      
      -- Variances Accounted for -----------------------------------------------------
      
                      F1
      SS loadings <num>
      Prop Tot Var <num>
      
      -- Model Fit -------------------------------------------------------------------
      
      χ²(20) = <num>, p = <num>
      CFI: <num>
      TLI: <num>
      RMSEA [90% CI]: <num> [ <num>; <num>]
      AIC: <num>
      BIC: <num>
      ECVI: <num>
      CAF: <num>
      RMSR: <num>
      SRMR: <num>

# summary.EFA_POOLED output is stable (PAF, promax)

    Code
      print(summary(pooled_obl))
    Output
      
      Pooled EFA across 3 imputations performed with type = 'EFAtools', method = 'PAF', and rotation = 'promax'.
      Pooling settings: target_method = 'first_target', align_unrotated = 'signed_tucker_congruence', fit_pool_method = 'D2'.
      
      -- Pooled Model Diagnostics ----------------------------------------------------
      
      Factors: 3
      Variables: 18
      N: 500
      Imputations: 3
      Pooling: target_method = 'first_target', align_unrotated = 'signed_tucker_congruence', fit_pool_method = 'D2'
      Alignment: method = 'first_target', converged
      Heywood cases: 0
      Cross-loading items (|loading| >= <num>): 0
      Items without salient loading (|loading| >= <num>): 0
      Factors with fewer than 3 salient indicators: 0
      Items with primary-loading gap < <num>: 1
      Largest |residual|: <num>
      Factor intercorrelations > <num>: none
      
      -- Rotated Loadings ------------------------------------------------------------
      
             F1     F2     F3    h2    u2
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
      
            F1 F2 F3
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
      
                           F1 F2 F3
      SS loadings <num> <num> <num>
      Prop Tot Var <num> <num> <num>
      Cum Prop Tot Var <num> <num> <num>
      Prop Comm Var <num> <num> <num>
      Cum Prop Comm Var <num> <num> <num>
      
      -- Model Fit -------------------------------------------------------------------
      
      CAF: <num>
      RMSR: <num>
      SRMR: <num>
      df: 102
      
      -- Residual Diagnostics --------------------------------------------------------
      
      Residual cutoff: |r| > <num>
      Number of large residuals: 0
      Largest absolute residual: <num>
      
      No absolute residuals > <num> occurred.
      
      Inspect the residual matrix for details (e.g., with residuals()).

# summary.EFA_POOLED output is stable (ML, unrotated)

    Code
      print(summary(pooled_none))
    Output
      
      Pooled EFA across 3 imputations performed with type = 'EFAtools', method = 'ML', and rotation = 'none'.
      Pooling settings: align_unrotated = 'signed_tucker_congruence', fit_pool_method = 'D2'.
      
      -- Pooled Model Diagnostics ----------------------------------------------------
      
      Factors: 1
      Variables: 8
      N: 810
      Imputations: 3
      Pooling: align_unrotated = 'signed_tucker_congruence', fit_pool_method = 'D2'
      Heywood cases: 0
      Cross-loading items (|loading| >= <num>): 0
      Items without salient loading (|loading| >= <num>): 0
      Factors with fewer than 3 salient indicators: 0
      Items with primary-loading gap < <num>: 0
      Largest |residual|: <num>
      
      -- Unrotated Loadings ----------------------------------------------------------
      
                  F1    h2    u2
      fun <num> <num> <num>
      friends <num> <num> <num>
      enjoy <num> <num> <num>
      hurt <num> <num> <num>
      part <num> <num> <num>
      commonly <num> <num> <num>
      chances <num> <num> <num>
      attracted <num> <num> <num>
      
      Legend:
        bold = |loading| >= <num>
        grey = below cutoff
        red h2/u2 = Heywood-relevant value
      
      -- Variances Accounted for -----------------------------------------------------
      
                      F1
      SS loadings <num>
      Prop Tot Var <num>
      
      -- Model Fit -------------------------------------------------------------------
      
      χ²(20) = <num>, p = <num>
      CFI: <num>
      TLI: <num>
      RMSEA [90% CI]: <num> [ <num>; <num>]
      AIC: <num>
      BIC: <num>
      ECVI: <num>
      CAF: <num>
      RMSR: <num>
      SRMR: <num>
      
      -- Residual Diagnostics --------------------------------------------------------
      
      Residual cutoff: |r| > <num>
      Number of large residuals: 0
      Largest absolute residual: <num>
      
      No absolute residuals > <num> occurred.
      
      Inspect the residual matrix for details (e.g., with residuals()).

