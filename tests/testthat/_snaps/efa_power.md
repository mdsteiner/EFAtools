# printed output is stable

    Code
      print(efa_power(df = 100, N = 200))
    Output
      
      -- RMSEA power analysis --------------------------------------------------------
      
      Test of close fit: H0 RMSEA ≤ .050 vs. H1 RMSEA = .080.
      alpha = .050 · df = 100
      
      Power = .955 at N = 200.

---

    Code
      print(efa_power(df = 100, power = 0.8))
    Output
      
      -- RMSEA power analysis --------------------------------------------------------
      
      Test of close fit: H0 RMSEA ≤ .050 vs. H1 RMSEA = .080.
      alpha = .050 · df = 100
      
      Required N = 132 for a power of .800 (achieved .802).

---

    Code
      print(efa_power(df = 100, N = 200, type = "notclose", group = 2))
    Output
      
      -- RMSEA power analysis --------------------------------------------------------
      
      Test of not-close fit: H0 RMSEA ≥ .050 vs. H1 RMSEA = .010.
      alpha = .050 · df = 100 · groups = 2
      
      Power = .426 at N = 200 (total; 100 per group).

# a criterion that fails on every replicate is reported with a hit-rate of 0

    Code
      print(sim)
    Output
      
      -- EFA power simulation --------------------------------------------------------
      
      18 variables · 3 factors · N = 120 · 3 datasets
      Estimation: PAF · rotation: oblimin
      Model error: none (an exact population gives optimistic rates).
      
      Retention hit-rate P(k-hat = 3)
      * MAP: <num> (no suggestion in 3)
      
      Structure recovery (Tucker congruence ≥ <num>)
      * recovery rate (min congruence): <num>
      * recovery rate (mean congruence): <num>
      * median min congruence: <num>
      
      Convergence: all 3 fits converged, without Heywood cases.
      
      Monte Carlo standard error of each rate: <num> or less.

# simulation-mode printed output is stable

    Code
      print(sim)
    Output
      
      -- EFA power simulation --------------------------------------------------------
      
      18 variables · 3 factors · N = 200 · 10 datasets
      Estimation: PAF · rotation: oblimin
      Model error: none (an exact population gives optimistic rates).
      
      Retention hit-rate P(k-hat = 3)
      * EKC_BvA2017: <num>
      * MAP_TR2: <num>
      * MAP_TR4: <num>
      
      Structure recovery (Tucker congruence ≥ <num>)
      * recovery rate (min congruence): <num>
      * recovery rate (mean congruence): <num>
      * median min congruence: <num>
      
      Convergence: all 10 fits converged, without Heywood cases.
      
      Monte Carlo standard error of each rate: <num> or less.

