# a failing group names the group and hands the cause on unaltered

    Code
      cat(conditionMessage(e))
    Output
      The `efa_fit()` fit failed for group "a".
      Caused by error:
      ! `estimator = "DWLS"` requires a polychoric asymptotic covariance.
      x You supplied a correlation matrix, so no asymptotic covariance can be estimated.
      i Supply raw ordinal data with `cor_method = "poly"` or `"tetra"`.

# print and format render the efa_group report

    Code
      print(mg)
    Output
      -- Multigroup exploratory factor analysis --------------------------------------
      
      2 groups (age_6_8, age_14_19) · 3 factors · PAF extraction · varimax rotation
      Aligned to a symmetric consensus target.
      N: age_6_8 = 825, age_14_19 = 1685
      
      -- Factor congruence -----------------------------------------------------------
      
      Tucker's congruence between the aligned group loadings (matched factors):
      
                             F1 F2 F3
      age_6_8 vs age_14_19 <num> <num> <num>
      
      -- Loading differences ---------------------------------------------------------
      
      Absolute loading differences by group pair (salience threshold |Δ| ≥ <num>):
      
                            mean  median   min   max  rmse  flagged
      age_6_8 vs age_14_19 <num> <num> <num> <num> <num> <flagged>
      
      -- Approximate invariance ------------------------------------------------------
      
      Congruence bands (Lorenzo-Seva & ten Berge, 2006: ≥ <num> equal, ≥ <num> fair):
      
                              F1 F2 F3
      age_6_8 vs age_14_19  equal  equal  equal

---

    Code
      print(mg_ident)
    Output
      -- Multigroup exploratory factor analysis --------------------------------------
      
      3 groups (a, b, c) · 2 factors · PAF extraction · varimax rotation
      Aligned to a symmetric consensus target.
      N: a = 500, b = 500, c = 500
      
      -- Factor congruence -----------------------------------------------------------
      
      Tucker's congruence between the aligned group loadings (matched factors):
      
                F1 F2
      a vs b <num> <num>
      a vs c <num> <num>
      b vs c <num> <num>
      
      -- Loading differences ---------------------------------------------------------
      
      Absolute loading differences by group pair (salience threshold |Δ| ≥ <num>):
      
              mean  median   min   max  rmse  flagged
      a vs b <num> <num> <num> <num> <num> <flagged>
      a vs c <num> <num> <num> <num> <num> <flagged>
      b vs c <num> <num> <num> <num> <num> <flagged>
      
      -- Approximate invariance ------------------------------------------------------
      
      Congruence bands (Lorenzo-Seva & ten Berge, 2006: ≥ <num> equal, ≥ <num> fair):
      
                F1 F2
      a vs b  equal  equal
      a vs c  equal  equal
      b vs c  equal  equal

---

    Code
      print(mg_ob)
    Output
      -- Multigroup exploratory factor analysis --------------------------------------
      
      2 groups (age_6_8, age_14_19) · 3 factors · PAF extraction · promax rotation
      Aligned to reference group age_6_8.
      N: age_6_8 = 825, age_14_19 = 1685
      
      -- Factor congruence -----------------------------------------------------------
      
      Tucker's congruence between the aligned group loadings (matched factors):
      
                             F1 F2 F3
      age_6_8 vs age_14_19 <num> <num> <num>
      
      -- Loading differences ---------------------------------------------------------
      
      Absolute loading differences by group pair (salience threshold |Δ| ≥ <num>):
      
                            mean  median   min   max  rmse  flagged
      age_6_8 vs age_14_19 <num> <num> <num> <num> <num> <flagged>

# print reports the bootstrap congruence intervals and verdicts

    Code
      print(mg)
    Output
      -- Multigroup exploratory factor analysis --------------------------------------
      
      2 groups (g1, g2) · 1 factor · PAF extraction · unrotated
      Aligned to a symmetric consensus target.
      N: g1 = 405, g2 = 405
      
      -- Factor congruence -----------------------------------------------------------
      
      Tucker's congruence between the aligned group loadings (matched factors):
      
                 F1
      g1 vs g2 <num>
      
      i 95% percentile bootstrap CIs (50 replicates) are in `$congruence$matched_ci`.
      
      -- Loading differences ---------------------------------------------------------
      
      Absolute loading differences by group pair (salience threshold |Δ| ≥ <num>):
      
                mean  median   min   max  rmse  flagged
      g1 vs g2 <num> <num> <num> <num> <num> <flagged>
      
      i 1 loading has a bootstrap CI excluding 0 (see `$flags`).
      
      -- Approximate invariance ------------------------------------------------------
      
      Congruence bands (Lorenzo-Seva & ten Berge, 2006), read off the CI lower bound:
      
                  F1
      g1 vs g2  equal

# print points out an 'equal' verdict that rests on scale-invariance

    Code
      print(mg)
    Output
      -- Multigroup exploratory factor analysis --------------------------------------
      
      2 groups (g1, g2) · 2 factors · PAF extraction · varimax rotation
      Aligned to a symmetric consensus target.
      N: g1 = 500, g2 = 500
      
      -- Factor congruence -----------------------------------------------------------
      
      Tucker's congruence between the aligned group loadings (matched factors):
      
                  F1 F2
      g1 vs g2 <num> <num>
      
      -- Loading differences ---------------------------------------------------------
      
      Absolute loading differences by group pair (salience threshold |Δ| ≥ <num>):
      
                mean  median   min   max  rmse  flagged
      g1 vs g2 <num> <num> <num> <num> <num> <flagged>
      
      -- Approximate invariance ------------------------------------------------------
      
      Congruence bands (Lorenzo-Seva & ten Berge, 2006: ≥ <num> equal, ≥ <num> fair):
      
                  F1 F2
      g1 vs g2  equal  equal
      
      i Every factor is graded "equal" while the aligned loadings differ by up to
       <num> on average: Tucker's congruence is invariant to a proportional rescaling
      of a factor's loadings, so read the verdicts alongside `$diffs`.

# the report wraps its header and splits its tables to the console width

    Code
      print(mg)
    Output
      -- Multigroup exploratory factor analysis ------------------
      
      4 groups (age_6_8, age_9_13, age_14_19, age_20_39) · 3
        factors · PAF extraction · varimax rotation
      Aligned to a symmetric consensus target.
      N: age_6_8 = 825, age_9_13 = 1572, age_14_19 = 1685,
        age_20_39 = 1251
      
      -- Factor congruence ---------------------------------------
      
      Tucker's congruence between the aligned group loadings
      (matched factors):
      
                               F1 F2 F3
      age_6_8 vs age_9_13 <num> <num> <num>
      age_6_8 vs age_14_19 <num> <num> <num>
      age_6_8 vs age_20_39 <num> <num> <num>
      age_9_13 vs age_14_19 <num> <num> <num>
      age_9_13 vs age_20_39 <num> <num> <num>
      age_14_19 vs age_20_39 <num> <num> <num>
      
      -- Loading differences -------------------------------------
      
      Absolute loading differences by group pair (salience
      threshold |Δ| ≥ <num>):
      
                              mean  median   min   max  rmse
      age_6_8 vs age_9_13 <num> <num> <num> <num> <num>
      age_6_8 vs age_14_19 <num> <num> <num> <num> <num>
      age_6_8 vs age_20_39 <num> <num> <num> <num> <num>
      age_9_13 vs age_14_19 <num> <num> <num> <num> <num>
      age_9_13 vs age_20_39 <num> <num> <num> <num> <num>
      age_14_19 vs age_20_39 <num> <num> <num> <num> <num>
      
                              flagged
      age_6_8 vs age_9_13 <flagged>
      age_6_8 vs age_14_19 <flagged>
      age_6_8 vs age_20_39 <flagged>
      age_9_13 vs age_14_19 <flagged>
      age_9_13 vs age_20_39 <flagged>
      age_14_19 vs age_20_39 <flagged>
      
      -- Approximate invariance ----------------------------------
      
      Congruence bands (Lorenzo-Seva & ten Berge, 2006: ≥ <num>
      equal, ≥ <num> fair):
      
                                F1 F2 F3
      age_6_8 vs age_9_13     equal  equal  equal
      age_6_8 vs age_14_19    equal  equal  equal
      age_6_8 vs age_20_39    equal  equal  equal
      age_9_13 vs age_14_19   equal  equal  equal
      age_9_13 vs age_20_39   equal  equal  equal
      age_14_19 vs age_20_39  equal  equal  equal

