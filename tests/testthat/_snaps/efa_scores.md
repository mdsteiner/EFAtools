# print and summary output are stable

    Code
      print(fs_grips)
    Output
      
      -- Factor scores (Bartlett) ----------------------------------------------------
      
      Scored 810 of 810 observations on 1 factor (see `$scores`).
      
      -- Score determinacy -----------------------------------------------------------
      
           rho  rho2  guttman
      F1 <num> <num> <num>

---

    Code
      print(fs_ob_cor)
    Output
      
      -- Factor scores (regression) --------------------------------------------------
      
      Weights and diagnostics only (correlation-matrix input; no scores).
      
      -- Score determinacy -----------------------------------------------------------
      
           rho  rho2  guttman
      F1 <num> <num> <num>
      F2 <num> <num> <num>
      F3 <num> <num> <num>

---

    Code
      print(summary(fs_ob_cor))
    Output
      
      -- Factor scores (regression) --------------------------------------------------
      
      Weights and diagnostics only (correlation-matrix input; no scores).
      
      -- Score determinacy -----------------------------------------------------------
      
           rho  rho2  guttman
      F1 <num> <num> <num>
      F2 <num> <num> <num>
      F3 <num> <num> <num>
      
      -- Factor weights --------------------------------------------------------------
      
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
      
      -- Score validity and univocality ----------------------------------------------
      
      Diagonal: validity (score-factor correlation). Off-diagonal: univocality.
      
           F1 F2 F3
      F1 <num> <num> <num>
      F2 <num> <num> <num>
      F3 <num> <num> <num>
      
      -- Score intercorrelations -----------------------------------------------------
      
            F1 F2 F3
      F1 <num> <num> <num>
      F2 <num> <num> <num>
      F3 <num> <num> <num>

