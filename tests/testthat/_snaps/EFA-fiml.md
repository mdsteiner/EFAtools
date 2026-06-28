# print/summary.EFA label FIML correlations in the header

    Code
      print(efa)
    Output
      
      EFA performed with type = 'EFAtools', method = 'ML', and rotation = 'none'.
      Correlations: FIML (two-stage, missing data)
      
      -- Unrotated Loadings ----------------------------------------------------------
      
            F1    F2    h2    u2
      V1 <num> <num> <num> <num>
      V2 <num> <num> <num> <num>
      V3 <num> <num> <num> <num>
      V4 <num> <num> <num> <num>
      V5 <num> <num> <num> <num>
      V6 <num> <num> <num> <num>
      
      Legend:
        bold = |loading| >= <num>
        grey = below cutoff
        red h2/u2 = Heywood-relevant value
      
      -- Variances Accounted for -----------------------------------------------------
      
                           F1 F2
      SS loadings <num> <num>
      Prop Tot Var <num> <num>
      Cum Prop Tot Var <num> <num>
      Prop Comm Var <num> <num>
      Cum Prop Comm Var <num> <num>
      
      -- Model Fit -------------------------------------------------------------------
      
      scaled χ²(4) = <num>, p = <num>
      CFI: <num>
      TLI: <num>
      RMSEA [90% CI]: <num> [ <num>; <num>]
      AIC: NA
      BIC: NA
      CAF: <num>
      RMSR: <num>
      SRMR: <num>

