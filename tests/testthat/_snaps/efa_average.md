# print output is stable

    Code
      print(efa_def, plot = FALSE)
    Output
      
      Averaging performed with averaging method mean (trim = 0) across 72 EFAs,
      varying the following settings: init_comm, criterion_type, k_promax, p_type,
      and varimax_type.
      
      The error rate is at <pct>. Of the solutions that did not result in an error, <pct>
      converged. Of the solutions that converged, <pct> contained Heywood cases and <pct>
      were admissible.
      
      == Indicator-to-Factor Correspondences =========================================
      
      For each cell, the proportion of solutions including the respective
      indicator-to-factor correspondence. A salience threshold of <num> was used to
      determine indicator-to-factor correspondences.
      
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
      
      == Loadings ====================================================================
      
      -- Mean ------------------------------------------------------------------------
      
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
      
      -- Range (max − min) -----------------------------------------------------------
      
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
      
      == Factor Intercorrelations from Oblique Solutions =============================
      
      -- Mean ------------------------------------------------------------------------
      
            F1 F2 F3
      F1 <num>
      F2 <num> <num>
      F3 <num> <num> <num>
      
      -- Range (max − min) -----------------------------------------------------------
      
           F1 F2 F3
      F1 <num>
      F2 <num> <num>
      F3 <num> <num> <num>
      
      == Variances Accounted for =====================================================
      
      -- Mean ------------------------------------------------------------------------
      
                       F1 F2 F3
      SS loadings <num> <num> <num>
      Prop Tot Var <num> <num> <num>
      Prop Comm Var <num> <num> <num>
      
      -- Range (max − min) -----------------------------------------------------------
      
                      F1 F2 F3
      SS loadings <num> <num> <num>
      Prop Tot Var <num> <num> <num>
      Prop Comm Var <num> <num> <num>
      
      == Model Fit ===================================================================
      
      The fit indices are the mean of the per-solution fit indices, not the fit of
      the averaged loadings printed above, which are a cell-wise summary rather than
      a fitted solution.
      
      CAF, RMSR, and SRMR averaged over <n> of 72 solutions.
      
             M (SD) [Min; Max]
      CAF: <num> ( <num>) [ <num>; <num>]
      RMSR: <num> ( <num>) [ <num>; <num>]
      SRMR: <num> ( <num>) [ <num>; <num>]
      df: 102

---

    Code
      print(efa_all_none, stat = c("average", "sd", "min", "max"), plot = FALSE)
    Output
      
      Averaging performed with averaging method mean (trim = 0) across 12 EFAs,
      varying the following settings: estimator, init_comm, criterion_type,
      abs_eigen, max_iter, and start_method.
      
      The error rate is at <pct>. Of the solutions that did not result in an error, <pct>
      converged. Of the solutions that converged, <pct> contained Heywood cases and <pct>
      were admissible.
      
      == Indicator-to-Factor Correspondences =========================================
      
      For each cell, the proportion of solutions including the respective
      indicator-to-factor correspondence. A salience threshold of <num> was used to
      determine indicator-to-factor correspondences.
      
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
      
      == Loadings ====================================================================
      
      -- Mean ------------------------------------------------------------------------
      
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
      
      -- Standard Deviation ----------------------------------------------------------
      
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
      
      -- Minimum ---------------------------------------------------------------------
      
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
      
      -- Maximum ---------------------------------------------------------------------
      
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
      
      == Variances Accounted for =====================================================
      
      -- Mean ------------------------------------------------------------------------
      
                       F1 F2 F3
      SS loadings <num> <num> <num>
      Prop Tot Var <num> <num> <num>
      Prop Comm Var <num> <num> <num>
      
      -- Standard Deviation ----------------------------------------------------------
      
                      F1 F2 F3
      SS loadings <num> <num> <num>
      Prop Tot Var <num> <num> <num>
      Prop Comm Var <num> <num> <num>
      
      -- Minimum ---------------------------------------------------------------------
      
                       F1 F2 F3
      SS loadings <num> <num> <num>
      Prop Tot Var <num> <num> <num>
      Prop Comm Var <num> <num> <num>
      
      -- Maximum ---------------------------------------------------------------------
      
                       F1 F2 F3
      SS loadings <num> <num> <num>
      Prop Tot Var <num> <num> <num>
      Prop Comm Var <num> <num> <num>
      
      == Model Fit ===================================================================
      
      The fit indices are the mean of the per-solution fit indices, not the fit of
      the averaged loadings printed above, which are a cell-wise summary rather than
      a fitted solution.
      
      Chi-square-based indices averaged over <n> of 12 solutions.
      CAF, RMSR, and SRMR averaged over <n> of 12 solutions.
      
             M (SD) [Min; Max]
      𝜒²: <num> ( <num>) [ <num>; <num>]
      df: 102
      p: <num> ( <num>) [ <num>; <num>]
      CFI: <num> ( <num>) [ <num>; <num>]
      TLI: <num> ( <num>) [ <num>; <num>]
      RMSEA: <num> ( <num>) [ <num>; <num>]
      AIC: <num> ( <num>) [ <num>; <num>]
      BIC: <num> ( <num>) [ <num>; <num>]
      ECVI: <num> ( <num>) [ <num>; <num>]
      CAF: <num> ( <num>) [ <num>; <num>]
      RMSR: <num> ( <num>) [ <num>; <num>]
      SRMR: <num> ( <num>) [ <num>; <num>]

---

    Code
      print(efa_all_md, plot = FALSE)
    Output
      
      Averaging performed with averaging method median across 183 EFAs, varying the
      following settings: estimator, init_comm, criterion_type, abs_eigen, max_iter,
      start_method, rotation, k_promax, p_type, and varimax_type.
      
      The error rate is at <pct>. Of the solutions that did not result in an error, <pct>
      converged. Of the solutions that converged, <pct> contained Heywood cases and <pct>
      were admissible.
      
      == Indicator-to-Factor Correspondences =========================================
      
      For each cell, the proportion of solutions including the respective
      indicator-to-factor correspondence. A salience threshold of <num> was used to
      determine indicator-to-factor correspondences.
      
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
      
      == Loadings ====================================================================
      
      -- Median ----------------------------------------------------------------------
      
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
      
      -- Range (max − min) -----------------------------------------------------------
      
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
      
      == Factor Intercorrelations from Oblique Solutions =============================
      
      -- Median ----------------------------------------------------------------------
      
            F1 F2 F3
      F1 <num>
      F2 <num> <num>
      F3 <num> <num> <num>
      
      -- Range (max − min) -----------------------------------------------------------
      
           F1 F2 F3
      F1 <num>
      F2 <num> <num>
      F3 <num> <num> <num>
      
      == Variances Accounted for =====================================================
      
      -- Median ----------------------------------------------------------------------
      
                       F1 F2 F3
      SS loadings <num> <num> <num>
      Prop Tot Var <num> <num> <num>
      Prop Comm Var <num> <num> <num>
      
      -- Range (max − min) -----------------------------------------------------------
      
                       F1 F2 F3
      SS loadings <num> <num> <num>
      Prop Tot Var <num> <num> <num>
      Prop Comm Var <num> <num> <num>
      
      == Model Fit ===================================================================
      
      The fit indices are the median of the per-solution fit indices, not the fit of
      the averaged loadings printed above, which are a cell-wise summary rather than
      a fitted solution.
      
      Chi-square-based indices averaged over <n> of 183 solutions.
      CAF, RMSR, and SRMR averaged over <n> of 183 solutions.
      
             Md (SD) [Min; Max]
      𝜒²: <num> ( <num>) [ <num>; <num>]
      df: 102
      p: <num> ( <num>) [ <num>; <num>]
      CFI: <num> ( <num>) [ <num>; <num>]
      TLI: <num> ( <num>) [ <num>; <num>]
      RMSEA: <num> ( <num>) [ <num>; <num>]
      AIC: <num> ( <num>) [ <num>; <num>]
      BIC: <num> ( <num>) [ <num>; <num>]
      ECVI: <num> ( <num>) [ <num>; <num>]
      CAF: <num> ( <num>) [ <num>; <num>]
      RMSR: <num> ( <num>) [ <num>; <num>]
      SRMR: <num> ( <num>) [ <num>; <num>]

