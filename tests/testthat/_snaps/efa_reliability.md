# print output is stable

    Code
      print(efa_reliability(sl_mod, factor_map = fc))
    Output
      
      Total variance from the correlation matrix.
      
      -- Reliability coefficients ----------------------------------------------------
      
           tot  hier   sub  alpha    H
      g <num> <num> <num> <num> <num>
      F1 <num> <num> <num> <num> <num>
      F2 <num> <num> <num> <num> <num>
      F3 <num> <num> <num> <num> <num>
      
      -- Common-variance indices -----------------------------------------------------
      
          ECV   PUC
      g <num> <num>

---

    Code
      print(efa_reliability(sl_mod, factor_map = fc, coefficients = c("omega_total",
        "alpha")))
    Output
      
      Total variance from the correlation matrix.
      
      -- Reliability coefficients ----------------------------------------------------
      
           tot  alpha
      g <num> <num>
      F1 <num> <num>
      F2 <num> <num>
      F3 <num> <num>

---

    Code
      print(efa_reliability(efa_mod))
    Output
      
      Total variance from the correlation matrix.
      
      Correlated-factors solution: a factor's total omega counts the true score
      variance its composite receives from every factor, through its cross-loadings
      and any factor correlations; its subscale omega counts only that factor's own
      contribution.
      
      -- Reliability coefficients ----------------------------------------------------
      
              tot   sub  alpha    H
      total <num> <num>
      F1 <num> <num> <num> <num>
      F2 <num> <num> <num> <num>
      F3 <num> <num> <num> <num>

