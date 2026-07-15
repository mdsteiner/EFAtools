# control objects print

    Code
      print(estimate_control(type = "SPSS"))
    Output
      Estimation control (type: "SPSS")
      
      init_comm: <from type preset>
      criterion: <from type preset>
      criterion_type: <from type preset>
      max_iter: <from type preset>
      abs_eigen: <from type preset>
      start_method: psych

---

    Code
      print(estimate_control(type = "none", init_comm = "smc", criterion = 0.001,
        criterion_type = "sum", max_iter = 300, abs_eigen = TRUE))
    Output
      Estimation control (type: "none")
      
      init_comm: smc
      criterion: 0.001
      criterion_type: sum
      max_iter: 300
      abs_eigen: TRUE
      start_method: psych

---

    Code
      print(rotate_control(type = "psych"))
    Output
      Rotation control (type: "psych")
      
      normalize: TRUE
      precision: 1e-05
      order_type: <from type preset>
      varimax_type: <from type preset>
      p_type: <from type preset>
      k: <from type preset>
      random_starts: 100

---

    Code
      print(rotate_control(type = "none", normalize = FALSE, order_type = "eigen",
        varimax_type = "kaiser", p_type = "norm", k = 4, random_starts = 50, gam = 0.5))
    Output
      Rotation control (type: "none")
      
      normalize: FALSE
      precision: 1e-05
      order_type: eigen
      varimax_type: kaiser
      p_type: norm
      k: 4
      random_starts: 50
      extra_args: gam = 0.5

