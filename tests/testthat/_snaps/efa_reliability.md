# print output is stable

    Code
      print(efa_reliability(sl_mod, factor_map = fc))
    Output
      
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
      
      -- Reliability coefficients ----------------------------------------------------
      
              tot  alpha    H
      total <num> <num>
      F1 <num> <num> <num>
      F2 <num> <num> <num>
      F3 <num> <num> <num>

# the model error names the averaged matrix an efa_average object carries

    Code
      cat(conditionMessage(e))
    Output
      Invalid `model`.
      i Enter an `efa_fit()`, `efa_schmid_leiman()`, `psych::schmid()`, or <lavaan> solution, or a loading matrix.
      i From an <efa_average> object, use the averaged loading matrix in `$loadings$average`, with `$Phi$average` as `Phi` and `$orig_R` as `cormat`.

