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

