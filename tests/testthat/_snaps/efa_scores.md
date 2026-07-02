# print and summary output are stable

    Code
      print(fs_grips)
    Output
      
      -- Factor scores (Bartlett) ----------------------------------------------------
      
      Scored 810 observations on 1 factor.
      
      -- Score determinacy -----------------------------------------------------------
      
           rho  rho2 guttman
      F1 0.972 0.946   0.891

---

    Code
      print(fs_ob_cor)
    Output
      
      -- Factor scores (regression) --------------------------------------------------
      
      Weights and diagnostics only (correlation-matrix input; no scores).
      
      -- Score determinacy -----------------------------------------------------------
      
           rho  rho2 guttman
      F1 0.894 0.798   0.597
      F2 0.888 0.788   0.576
      F3 0.883 0.780   0.561

---

    Code
      print(summary(fs_ob_cor))
    Output
      
      -- Factor scores (regression) --------------------------------------------------
      
      Weights and diagnostics only (correlation-matrix input; no scores).
      
      -- Score determinacy -----------------------------------------------------------
      
           rho  rho2 guttman
      F1 0.894 0.798   0.597
      F2 0.888 0.788   0.576
      F3 0.883 0.780   0.561
      
      -- Factor weights --------------------------------------------------------------
      
             F1    F2     F3
      V1  0.016 0.037  0.206
      V2  0.023 0.036  0.146
      V3  0.038 0.033  0.140
      V4  0.060 0.025  0.194
      V5  0.062 0.014  0.138
      V6  0.009 0.013  0.247
      V7  0.024 0.177  0.053
      V8  0.017 0.187  0.031
      V9  0.031 0.173  0.020
      V10 0.016 0.222 -0.002
      V11 0.026 0.115  0.084
      V12 0.035 0.240  0.030
      V13 0.202 0.051  0.010
      V14 0.163 0.002  0.050
      V15 0.177 0.059  0.007
      V16 0.170 0.006  0.050
      V17 0.214 0.013  0.016
      V18 0.173 0.025  0.041
      
      -- Score validity and univocality ----------------------------------------------
      
      Diagonal: validity (score-factor correlation). Off-diagonal: univocality.
      
            F1    F2    F3
      F1 0.894 0.638 0.668
      F2 0.643 0.888 0.650
      F3 0.676 0.653 0.883
      
      -- Score intercorrelations -----------------------------------------------------
      
            F1    F2    F3
      F1 1.000 0.719 0.756
      F2 0.719 1.000 0.735
      F3 0.756 0.735 1.000

