# print output is stable

    Code
      print(kmo_cor)
    Output
      
      -- Kaiser-Meyer-Olkin criterion (KMO) ------------------------------------------
      
      v The overall KMO value for your data is marvellous.
        These data are probably suitable for factor analysis.
      
        Overall: <num>
      
        For each variable:
         V1    V2    V3    V4    V5    V6    V7    V8    V9   V10   V11   V12   V13 
       <num> <num> <num> <num> <num> <num> <num> <num> <num> <num> <num> <num> <num> 
        V14   V15   V16   V17   V18 
       <num> <num> <num> <num> <num> 

---

    Code
      print(kmo_low)
    Output
      
      -- Kaiser-Meyer-Olkin criterion (KMO) ------------------------------------------
      
      x The overall KMO value for your data is unacceptable.
        These data are not suitable for factor analysis.
      
        Overall: <num>
      
        For each variable:
        V1   V2   V3 
       <num> <num> <num> 

---

    Code
      print(kmo_na)
    Output
      
      -- Kaiser-Meyer-Olkin criterion (KMO) ------------------------------------------
      
      ! Sorry, the KMO value for your data is not available.

