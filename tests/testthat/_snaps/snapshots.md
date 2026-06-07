# print.SMT output is stable

    Code
      print(smt)
    Output
      
      Sequential 𝜒² Model Tests suggest 3 factors.
      
      Lower bound of RMSEA 90% confidence interval suggests 2 factors.
      
      AIC suggests 3 factors.
      

---

    Code
      print(smt_id)
    Output
      
      Sequential 𝜒² Model Tests suggest 0 factors.
      
      Lower bound of RMSEA 90% confidence interval suggests 0 factors.
      
      AIC suggests 0 factors.
      

# print.SCREE output is stable

    Code
      print(scree, plot = FALSE)
    Output
      
      Eigenvalues were found using PCA, SMC, and EFA.
      
      

---

    Code
      print(scree_smc, plot = FALSE)
    Output
      
      Eigenvalues were found using SMC.
      
      

# print.KGC output is stable

    Code
      print(kgc, plot = FALSE)
    Output
      
      Eigenvalues were found using PCA, SMC, and EFA.
      
      -- Number of factors suggested by Kaiser-Guttmann criterion --------------------
      
      * With PCA-determined eigenvalues:  3
      * With SMC-determined eigenvalues:  1
      * With EFA-determined eigenvalues:  1
      

---

    Code
      print(kgc_smc, plot = FALSE)
    Output
      
      Eigenvalues were found using SMC.
      
      -- Number of factors suggested by Kaiser-Guttmann criterion --------------------
      
      * With SMC-determined eigenvalues:  1
      

# print.EKC output is stable

    Code
      print(ekc, plot = FALSE)
    Output
      -- Number of factors suggested by EKC ------------------------------------------
      
      * Original implementation (Braeken & van Assen, 2017):  3
    Message
      i Different implementations of EKC exist. Make sure to report which one you relied on (see EKC help page for details)

---

    Code
      print(ekc_both, plot = FALSE)
    Output
      -- Number of factors suggested by EKC ------------------------------------------
      
      * Original implementation (Braeken & van Assen, 2017):  3
      * Adapted implementation (Auerswald & Moshagen, 2019):  2
    Message
      i Different implementations of EKC exist. Make sure to report which one you relied on (see EKC help page for details)

# print.MAP output is stable

    Code
      print(map, plot = FALSE)
    Output
      -- Number of factors suggested by MAP ------------------------------------------
      
      * Original implementation (TR2):  1
      * Revised implementation (TR4):  3

