# print output is stable

    Code
      print(bart_cor)
    Output
      
      v Bartlett's test of sphericity was significant at an alpha level of <num>.
      These data are probably suitable for factor analysis.
      
      χ²(153) = <num>, p < <num>

---

    Code
      print(bart_rand)
    Output
      
      x Bartlett's test of sphericity was not significant at an alpha level of <num>.
      These data are probably not suitable for factor analysis.
      
      χ²(6) = <num>, p = <num>

---

    Code
      print(bart_na)
    Output
      
      ! Bartlett's test of sphericity did not render a result.
      The Bartlett multiplier N - 1 - (2p + 5)/6 is not positive for a sample this
      small relative to the number of variables.
      
      χ²(NA) = NA, p = NA

---

    Code
      print(bart_null)
    Output
      
      ! Bartlett's test of sphericity did not render a result.
      The Bartlett multiplier N - 1 - (2p + 5)/6 is not positive for a sample this
      small relative to the number of variables.
      
      χ²() = NA, p = NA

