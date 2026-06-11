# print output is stable

    Code
      print(bart_cor)
    Output
      
      v The Bartlett's test of sphericity was significant at an alpha level of <num>.
        These data are probably suitable for factor analysis.
      
        𝜒²(153) = <num>, p < <num>

---

    Code
      print(bart_rand)
    Output
      
      x The Bartlett's test of sphericity was not significant at an alpha level of <num>.
        These data are probably not suitable for factor analysis.
      
        𝜒²(6) = <num>, p = <num>

---

    Code
      print(bart_na)
    Output
      
      ! The Bartlett's test of sphericity did not render a result.
      
      
        𝜒²(NA) = NA, pNA

