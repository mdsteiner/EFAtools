# the input error names where a fitted solution keeps its loadings

    Code
      cat(conditionMessage(e))
    Output
      `x` (<efa/EFA>) and `y` (<efa/EFA>) must be numeric vectors or matrices.
      i From an <efa> solution, use `$rot_loadings`, or `$unrot_loadings` if it was fitted without a rotation.

# print output is stable

    Code
      print(matr)
    Output
      
      -- Summary statistics ----------------------------------------------------------
      
      Mean [min, max] absolute difference:  .2500 [ .0000, 1.0000]
      Median absolute difference:  .0000
      Root mean squared distance (RMSE):  .5000
      Max decimals where all numbers agree in absolute value: none
      Minimum number of decimals provided: 0
      Differing indicator-to-factor correspondences: 1 (highest loading),
        0 (all |loadings| >= 0.3)
      
      -- Elementwise differences -----------------------------------------------------
      
      Differences: x - y.
      
            F1     F2
      V1  .0000   .0000
      V2  .0000  1.0000

---

    Code
      print(int)
    Output
      
      -- Summary statistics ----------------------------------------------------------
      
      Mean [min, max] absolute difference:  .0000 [ .0000,  .0000]
      Median absolute difference:  .0000
      Root mean squared distance (RMSE):  .0000
      Max decimals where all numbers agree in absolute value: 0
      Minimum number of decimals provided: 0
      
      -- Elementwise differences -----------------------------------------------------
      
      Differences: x - y.
      
      V1   .0000
      V2   .0000
      V3   .0000
      V4   .0000
      V5   .0000
      V6   .0000
      V7   .0000
      V8   .0000
      V9   .0000
      V10  .0000

