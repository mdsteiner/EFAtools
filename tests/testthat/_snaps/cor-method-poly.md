# the degenerate-variance warning reads correctly for one and for several pairs

    Code
      cat(degenerate_msg(c(`v1-v2` = 0.01, `v1-v3` = -1)))
    Output
      The polychoric asymptotic covariance is unavailable for variable pair "v1-v3".
      x This pair was either estimated with a continuity correction for a response combination that never occurs -- which repairs the correlation but leaves it without a variance -- or is not identified well enough to have a usable one (e.g. a perfectly ordered or near-empty response table).
      i Any robust/sandwich standard error involving it is withheld. A DWLS fit refuses the data outright when the variance is missing or non-positive, and down-weights it out of the solution when it is merely far too large.
      i Fitting with `estimator = "ULS"` avoids the weights entirely; collapsing rare response categories in the variables involved resolves the under-identified cases.

---

    Code
      cat(degenerate_msg(c(`v1-v2` = NA_real_, `v1-v3` = -1)))
    Output
      The polychoric asymptotic covariance is unavailable for variable pairs "v1-v2" and "v1-v3".
      x These pairs were either estimated with a continuity correction for a response combination that never occurs -- which repairs the correlation but leaves it without a variance -- or are not identified well enough to have a usable one (e.g. a perfectly ordered or near-empty response table).
      i Any robust/sandwich standard error involving them is withheld. A DWLS fit refuses the data outright when the variance is missing or non-positive, and down-weights them out of the solution when it is merely far too large.
      i Fitting with `estimator = "ULS"` avoids the weights entirely; collapsing rare response categories in the variables involved resolves the under-identified cases.

