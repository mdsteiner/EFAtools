# print.efa_retention output is stable for SMT

    Code
      print(smt)
    Output
      -- Sequential model tests ------------------------------------------------------
      
      * Sequential chi-square model tests: 3
      * Lower bound of RMSEA 90% CI: 2
      * Akaike Information Criterion: 3

---

    Code
      print(smt_id)
    Output
      -- Sequential model tests ------------------------------------------------------
      
      * Sequential chi-square model tests: 0
      * Lower bound of RMSEA 90% CI: 0
      * Akaike Information Criterion: 0

# print.efa_retention output is stable for SCREE

    Code
      print(scree)
    Output
      -- Scree plot ------------------------------------------------------------------
      Eigenvalues found using PCA, SMC, and EFA.
      
      i Scree plot is a visual criterion; inspect the plot to identify the elbow.

---

    Code
      print(scree_smc)
    Output
      -- Scree plot ------------------------------------------------------------------
      Eigenvalues found using SMC.
      
      i Scree plot is a visual criterion; inspect the plot to identify the elbow.

# print.efa_retention output is stable for CD

    Code
      print(cd)
    Output
      -- Comparison data -------------------------------------------------------------
      
      * Suggested number of factors: 1

# print.efa_retention output is stable for PARALLEL

    Code
      print(pa)
    Output
      -- Parallel analysis -----------------------------------------------------------
      Eigenvalues found using PCA; 1000 simulated datasets.
      
      * PCA eigenvalues: 3
      
      i Number of factors retained using the "means" decision rule.

---

    Code
      print(pa_nodat)
    Output
      -- Parallel analysis -----------------------------------------------------------
      Eigenvalues found using PCA, SMC, and EFA; 1000 simulated datasets.
      
      i No data were entered; showing the simulated eigenvalues only. No number of factors is suggested.

# print.efa_retention output is stable for KGC

    Code
      print(kgc)
    Output
      -- Kaiser-Guttman criterion ----------------------------------------------------
      
      * PCA eigenvalues: 3
      * SMC eigenvalues: 1
      * EFA eigenvalues: 1

---

    Code
      print(kgc_smc)
    Output
      -- Kaiser-Guttman criterion ----------------------------------------------------
      
      * SMC eigenvalues: 1

# print.efa_retention output is stable for NEST

    Code
      print(nest)
    Output
      -- Next Eigenvalue Sufficiency Test --------------------------------------------
      
      * Suggested number of factors: 3

# print.efa_retention output is stable for EKC

    Code
      print(ekc)
    Output
      -- Empirical Kaiser Criterion --------------------------------------------------
      
      * Original implementation (Braeken & van Assen, 2017): 3
      
      i Multiple implementations of EKC exist; make sure to report which one you used (see the EKC help page for details).

---

    Code
      print(ekc_both)
    Output
      -- Empirical Kaiser Criterion --------------------------------------------------
      
      * Original implementation (Braeken & van Assen, 2017): 3
      * Adapted implementation (Auerswald & Moshagen, 2019): 2
      
      i Multiple implementations of EKC exist; make sure to report which one you used (see the EKC help page for details).

# print.efa_retention output is stable for HULL

    Code
      print(hull)
    Output
      -- Hull method -----------------------------------------------------------------
      Estimation method: ML
      
      * CAF: 3
      * CFI: 1
      * RMSEA: 1

---

    Code
      print(hull_paf)
    Output
      -- Hull method -----------------------------------------------------------------
      Estimation method: PAF
      
      * CAF: 3

# print.efa_retention output is stable for MAP

    Code
      print(map)
    Output
      -- Minimum average partial -----------------------------------------------------
      
      * Original implementation (TR2): 1
      * Revised implementation (TR4): 3

