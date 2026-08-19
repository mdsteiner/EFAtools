# print output is stable

    Code
      print(SL_EFAtools)
    Output
      
      EFA for second-order loadings performed with estimator = 'PAF'
      
      -- Schmid-Leiman Solution ------------------------------------------------------
      
             g     F1     F2     F3    h2    u2
      V1 <num> <num> <num> <num> <num> <num>
      V2 <num> <num> <num> <num> <num> <num>
      V3 <num> <num> <num> <num> <num> <num>
      V4 <num> <num> <num> <num> <num> <num>
      V5 <num> <num> <num> <num> <num> <num>
      V6 <num> <num> <num> <num> <num> <num>
      V7 <num> <num> <num> <num> <num> <num>
      V8 <num> <num> <num> <num> <num> <num>
      V9 <num> <num> <num> <num> <num> <num>
      V10 <num> <num> <num> <num> <num> <num>
      V11 <num> <num> <num> <num> <num> <num>
      V12 <num> <num> <num> <num> <num> <num>
      V13 <num> <num> <num> <num> <num> <num>
      V14 <num> <num> <num> <num> <num> <num>
      V15 <num> <num> <num> <num> <num> <num>
      V16 <num> <num> <num> <num> <num> <num>
      V17 <num> <num> <num> <num> <num> <num>
      V18 <num> <num> <num> <num> <num> <num>
      
      -- Variances Accounted for -----------------------------------------------------
      
                           g     F1    F2     F3
      SS loadings <num> <num> <num> <num>
      Prop Tot Var <num> <num> <num> <num>
      Cum Prop Tot Var <num> <num> <num> <num>
      Prop Comm Var <num> <num> <num> <num>
      Cum Prop Comm Var <num> <num> <num> <num>

# the input error names the averaged matrix an efa_average object carries

    Code
      cat(conditionMessage(e))
    Output
      `x` must be an <EFA>, <fa>, or <lavaan> object, a matrix, or a <LOADINGS>/<loadings> object.
      i From an <efa_average> object, use the averaged loading matrix in `$loadings$average`, with the averaged factor correlations in `$Phi$average`.

