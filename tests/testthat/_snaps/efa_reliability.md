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
      
      -- Reliability coefficients ----------------------------------------------------
      
           tot   sub  alpha    H
      g <num> <num>
      F1 <num> <num> <num> <num>
      F2 <num> <num> <num> <num>
      F3 <num> <num> <num> <num>

