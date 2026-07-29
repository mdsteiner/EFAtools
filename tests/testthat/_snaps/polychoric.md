# the sparse-cell warning names the offending pairs and caps the list

    Code
      cat(sparse_msg(one))
    Output
      1 variable pair has an empty response-category combination despite a non-negligible expected count.
      x Affected pair: "risk-gamble".
      i The polychoric asymptotic covariance (and any DWLS weights or robust standard errors derived from it) can be unreliable for such structurally sparse cells; interpret them with caution and consider collapsing rare response categories in these variables.

---

    Code
      cat(sparse_msg(nudged))
    Output
      15 variable pairs have an empty response-category combination despite a non-negligible expected count.
      x Affected pairs: "v1-v2", "v1-v3", "v1-v4", "v1-v5", "v1-v6", and 10 more.
      i The polychoric asymptotic covariance (and any DWLS weights or robust standard errors derived from it) can be unreliable for such structurally sparse cells; interpret them with caution and consider collapsing rare response categories in these variables.

# the boundary warning names the offending pairs and caps the list

    Code
      cat(boundary_msg(one))
    Output
      1 variable pair has a perfectly monotone response table, so its correlation is estimated at the boundary value 0.9999.
      x Affected pair: "risk-gamble".
      i The two variables are never observed out of order -- either they always vary together, or (for a negative estimate) they always vary in opposite directions -- which the model can reproduce only with a correlation of 1 or -1, so the data put no limit on how close to that boundary the correlation lies and no asymptotic variance or standard error is available.
      i Consider collapsing rare response categories in the variables involved, or dropping one item of a redundant pair.

---

    Code
      cat(boundary_msg(rev))
    Output
      1 variable pair has a perfectly monotone response table, so its correlation is estimated at the boundary value -0.9999.
      x Affected pair: "risk-gamble".
      i The two variables are never observed out of order -- either they always vary together, or (for a negative estimate) they always vary in opposite directions -- which the model can reproduce only with a correlation of 1 or -1, so the data put no limit on how close to that boundary the correlation lies and no asymptotic variance or standard error is available.
      i Consider collapsing rare response categories in the variables involved, or dropping one item of a redundant pair.

---

    Code
      cat(boundary_msg(both))
    Output
      3 variable pairs have a perfectly monotone response table, so their correlation is estimated at the boundary values 0.9999 and -0.9999.
      x Affected pairs: "risk-gamble", "risk-mirror", and "gamble-mirror".
      i The two variables are never observed out of order -- either they always vary together, or (for a negative estimate) they always vary in opposite directions -- which the model can reproduce only with a correlation of 1 or -1, so the data put no limit on how close to that boundary the correlation lies and no asymptotic variance or standard error is available.
      i Consider collapsing rare response categories in the variables involved, or dropping one item of a redundant pair.

# the boundary warning quotes the estimate, not the nearest-PD projection of it

    Code
      cat(conditionMessage(w))
    Output
      1 variable pair has a perfectly monotone response table, so its correlation is estimated at the boundary value 0.9999.
      x Affected pair: "v1-v2".
      i The two variables are never observed out of order -- either they always vary together, or (for a negative estimate) they always vary in opposite directions -- which the model can reproduce only with a correlation of 1 or -1, so the data put no limit on how close to that boundary the correlation lies and no asymptotic variance or standard error is available.
      i Consider collapsing rare response categories in the variables involved, or dropping one item of a redundant pair.

# the continuity-correction warning names the offending pairs and caps the list

    Code
      cat(zero_cell_msg(guttman))
    Output
      21 binary variable pairs have a response combination that never occurs, so a continuity correction of 0.5 was applied before estimating their correlation.
      x Affected pairs: "i1-i2", "i1-i3", "i1-i4", "i1-i5", "i1-i6", and 16 more.
      i Every two-by-two table with an empty cell is reproduced exactly by a correlation of 1 (or -1), so without the correction the estimate would be the boundary value whatever the underlying correlation. The correction adds 0.5 to the empty cell while preserving the table margins, as lavaan and psych do by default.
      i The corrected correlation is a point estimate only: no asymptotic variance or standard error is available for these pairs.

