# printed output is stable

    Code
      print(efa_power(df = 100, N = 200))
    Output
      
      -- RMSEA power analysis --------------------------------------------------------
      
      Test of close fit: H0 RMSEA ≤ .050 vs. H1 RMSEA = .080.
      alpha = .050 · df = 100
      
      Power = .955 at N = 200.
      Critical value χ²(100) = 183.967 · noncentrality H0 = 49.750, H1 = 127.360.

---

    Code
      print(efa_power(df = 100, power = 0.8))
    Output
      
      -- RMSEA power analysis --------------------------------------------------------
      
      Test of close fit: H0 RMSEA ≤ .050 vs. H1 RMSEA = .080.
      alpha = .050 · df = 100
      
      Required N = 132 for a power of .800 (achieved .802).
      Critical value χ²(100) = 163.977 · noncentrality H0 = 32.750, H1 = 83.840.

---

    Code
      print(efa_power(df = 100, N = 200, type = "notclose", group = 2))
    Output
      
      -- RMSEA power analysis --------------------------------------------------------
      
      Test of not-close fit: H0 RMSEA ≥ .050 vs. H1 RMSEA = .010.
      alpha = .050 · df = 100 · groups = 2
      
      Power = .429 at N = 200.
      Critical value χ²(100) = 97.795 · noncentrality H0 = 24.875, H1 = .995.

# simulation-mode printed output is stable

    Code
      print(sim)
    Output
      
      -- EFA power simulation --------------------------------------------------------
      
      18 variables · 3 factors · N = 200 · 10 datasets
      Estimation: PAF · rotation: promax
      Model error: none. The population is exact, so the hit-rate and recovery are
      optimistic; set `target_rmsea` for realism.
      
      Retention hit-rate P(k-hat = 3)
      * EKC_BvA2017: 1.000 (n = 10)
      * MAP_TR2: 1.000 (n = 10)
      * MAP_TR4: 1.000 (n = 10)
      
      Structure recovery (Tucker congruence ≥ .950)
      * min congruence: 1.000 (n = 10)
      * mean congruence: 1.000 (n = 10)
      
      Convergence
      * fits completed: 1.000 (10/10)
      * converged (of completed): 1.000
      * Heywood cases (of completed): .000

