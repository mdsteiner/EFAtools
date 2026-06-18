# print output is stable

    Code
      print(om_sl)
    Output
      Omega total, omega hierarchical, omega subscale, H index, explained common
      variance (ECV), and percent of uncontaminated correlations (PUC) for the
      general factor (top row) and omegas and H index for the group factors:
      
           tot  hier   sub     H   ECV   PUC
      g <num> <num> <num> <num> <num> <num>
      F1 <num> <num> <num> <num>            
      F2 <num> <num> <num> <num>            
      F3 <num> <num> <num> <num>            

---

    Code
      print(structure(0.85, class = "OMEGA"))
    Output
      Omega total for the single factor: <num>

---

    Code
      print(structure(c(0.85, 2.1), class = "OMEGA"))
    Output
      Omega total and H index for the single factor:
      
      Omega: <num>
      H index: <num>

---

    Code
      print(om_mg6)
    Output
      Omega total, omega hierarchical, omega subscale, H index, explained common
      variance (ECV), and percent of uncontaminated correlations (PUC) for the
      general factor (top row) and omegas and H index for the group factors for each
      group:
      
      Group GroupA:
           tot  hier   sub     H   ECV   PUC
      g <num> <num> <num> <num> <num> <num>
      F1 <num> <num> <num> <num>            
      F2 <num> <num> <num> <num>            
      F3 <num> <num> <num> <num>            
      
      Group GroupB:
           tot  hier   sub     H   ECV   PUC
      g <num> <num> <num> <num> <num> <num>
      F1 <num> <num> <num> <num>            
      F2 <num> <num> <num> <num>            
      F3 <num> <num> <num> <num>            

---

    Code
      print(om_mg3)
    Output
      Omega total, omega hierarchical, and omega subscale for the general factor (top
      row) and the group factors for each group:
      
      Group GroupA:
           tot  hier   sub
      g <num> <num> <num>
      F1 <num> <num> <num>
      F2 <num> <num> <num>
      F3 <num> <num> <num>
      
      Group GroupB:
           tot  hier   sub
      g <num> <num> <num>
      F1 <num> <num> <num>
      F2 <num> <num> <num>
      F3 <num> <num> <num>

---

    Code
      print(structure(list(GroupA = 0.85, GroupB = 0.8), class = "OMEGA"))
    Output
      Omega total for the single factor for each group:
      
      Group GroupA: <num>
      
      Group GroupB: <num>

---

    Code
      print(structure(list(GroupA = c(0.85, 2.1), GroupB = c(0.8, 1.9)), class = "OMEGA"))
    Output
      Omega total and H index for the single factor for each group:
      
      Group GroupA:
      Omega: <num>
      H index: <num>
      
      Group GroupB:
      Omega: <num>
      H index: <num>

