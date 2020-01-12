#' Print PROMAX object
#'
#' Print method showing a summarized output of the PROMAX procedure.
#'
#' @param x list. An object of class PROMAX to be printed
#' @param ...  Further arguments for print.
#'
#' @export
#' @method print PROMAX
#'
#' @examples
#' EFA(IDS2_R, n_factors = 5, type = "GS", rotation = "promax")
print.PROMAX <- function(x, ...) {


  if (length(x$settings) == 6) {

    # extract the different settings

    #Varimax
    within_EFA <- FALSE
    type <- x$settings$type
    kaiser <- x$settings$kaiser
    precision <- x$settings$precision
    order_type <- x$settings$order_type

    # Promax
    P_type <- x$settings$P_type
    k <- x$settings$k

  } else {
    # extract the different settings
    N <- x$settings$N
    max_iter <- x$settings$max_iter
    type <- x$settings$type
    init_comm <- x$settings$init_comm
    criterion <- x$settings$criterion
    criterion_type <- x$settings$criterion_type
    abs_eigen <- x$settings$abs_eigen
    signed_loadings <- x$settings$signed_loadings
    use <- x$settings$use

    # Varimax
    within_EFA <- TRUE
    type <- x$settings$type
    kaiser <- x$settings$kaiser
    precision <- x$settings$precision
    order_type <- x$settings$order_type

    # Promax
    P_type <- x$settings$P_type
    k <- x$settings$k
  }



  # get the loadings
  ldngs <- x$rot_loadings

  # get the variances accounted
  vrs_acc <- x$vars_accounted

  if (isTRUE(within_EFA)) {

    if (x$iter > max_iter) {
      conv_warn <- crayon::red$bold(cli::symbol$cross,
                                    "Maximum number of iterations reached",
                                    "without convergence")
    }

    # sample size message
    if (!is.na(N)) {
      N_used <- crayon::blue("Sample Size:", crayon::bold(N))
    } else {
      N_used <- crayon::blue("Sample Size:", crayon::bold("Not Provided"))
    }

    # Settings Intro message
    sttngs_intro <- crayon::blue("PAF with type =", crayon::bold(type),
                                 "performed with the",
                                 "following settings:")
    # Varimax settings Intro message
    sttngs_intro_vrmx <- crayon::blue("Varimax with type =", crayon::bold(type),
                                      "performed with the",
                                      "following settings:")
    # Promax settings Intro message
    sttngs_intro_prmx <- crayon::blue("Promax with type =", crayon::bold(type),
                                      "performed with the",
                                      "following settings:")

    if (type == "GS") {
      # test if all type GS arguments are on default

      if (init_comm == "smc") {
        i_c_used <- crayon::blue("    ",cli::symbol$bullet, "Initial Communalities:",
                                 crayon::bold(init_comm),
                                 " (", crayon::green(cli::symbol$tick), "default)")
        def_i_c_used <- TRUE
      } else {
        i_c_used <- crayon::blue("    ",cli::symbol$bullet, "Initial Communalities:",
                                 crayon::bold(init_comm),
                                 " (", crayon::red(cli::symbol$cross), "default)")
        def_i_c_used <- FALSE
      }

      if (criterion == 1e-9) {
        crit_used <- crayon::blue("    ",cli::symbol$bullet, "Convergence Criterion:",
                                  crayon::bold(criterion),
                                  " (", crayon::green(cli::symbol$tick),
                                  "default)")
        def_crit_used <- TRUE
      } else {
        crit_used <- crayon::blue("    ",cli::symbol$bullet, "Initial Communalities:",
                                  crayon::bold(criterion),
                                  " (", crayon::red(cli::symbol$cross),
                                  "default)")
        def_crit_used <- FALSE
      }


      if (criterion_type == "max_individual") {
        crit_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Criterion Type:",
                                     crayon::bold(criterion_type),
                                     " (", crayon::green(cli::symbol$tick),
                                     "default)")
        def_crit_tp_used <- TRUE
      } else {
        crit_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Criterion Type:",
                                     crayon::bold(criterion_type),
                                     " (", crayon::red(cli::symbol$cross),
                                     "default)")
        def_crit_tp_used <- FALSE
      }

      if (max_iter == 1e5) {
        mx_it_used <- crayon::blue("    ",cli::symbol$bullet, "Maximum Iteration:",
                                   crayon::bold(max_iter),
                                   " (", crayon::green(cli::symbol$tick),
                                   "default)")
        def_mx_it_used <- TRUE
      } else {
        mx_it_used <- crayon::blue("    ",cli::symbol$bullet, "Maximum Iteration:",
                                   crayon::bold(max_iter),
                                   " (", crayon::red(cli::symbol$cross),
                                   "default)")
        def_mx_it_used <- FALSE
      }

      if (isFALSE(abs_eigen)) {
        abs_eigen_used <- crayon::blue("    ",cli::symbol$bullet, "Absolute Eigenvalues:",
                                       crayon::bold("No"),
                                       " (", crayon::green(cli::symbol$tick),
                                       "default)")
        def_abs_eigen_used <- TRUE
      } else {
        abs_eigen_used <- crayon::blue("    ",cli::symbol$bullet, "Absolute Eigenvalues:",
                                       crayon::bold("Yes"),
                                       " (", crayon::red(cli::symbol$cross),
                                       "default)")
        def_abs_eigen_used <- FALSE
      }

      # Varimax settings
      if (isTRUE(kaiser)) {
        kaiser_used <- crayon::blue("    ",cli::symbol$bullet, "Kaiser Normalization:",
                                    crayon::bold("Yes"),
                                    " (", crayon::green(cli::symbol$tick),
                                    "default)")
        def_kaiser_used <- TRUE
      } else {
        kaiser_used <- crayon::blue("    ",cli::symbol$bullet, "Kaiser Normalization:",
                                    crayon::bold("No"),
                                    " (", crayon::red(cli::symbol$cross),
                                    "default)")
        def_kaiser_used <- FALSE
      }

      if (precision == 1e-10) {
        prcsn_used <- crayon::blue("    ",cli::symbol$bullet, "Precision:",
                                   crayon::bold(precision),
                                   " (", crayon::green(cli::symbol$tick),
                                   "default)")
        def_prcsn_used <- TRUE
      } else {
        prcsn_used <- crayon::blue("    ",cli::symbol$bullet, "Precision:",
                                   crayon::bold(precision),
                                   " (", crayon::red(cli::symbol$cross),
                                   "default)")
        def_prcsn_used <- FALSE
      }

      if (order_type == "eigen") {
        ordr_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Factors ordered according to:",
                                     crayon::bold("Eigenvalues"),
                                     " (", crayon::green(cli::symbol$tick),
                                     "default)")
        def_ordr_tp_used <- TRUE
      } else {
        ordr_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Factors ordered according to:",
                                     crayon::bold("Sum of Squared Loadings"),
                                     " (", crayon::red(cli::symbol$cross),
                                     "default)")
        def_ordr_tp_used <- FALSE
      }

      # Promax settings

      if (P_type == "unnorm") {
        trgt_mtrx_used <- crayon::blue("    ",cli::symbol$bullet, "Target matrix:",
                                     crayon::bold("Unnormalized"),
                                     " (", crayon::green(cli::symbol$tick),
                                     "default)")
        def_trgt_mtrx_used <- TRUE
      } else {
        trgt_mtrx_used <- crayon::blue("    ",cli::symbol$bullet, "Target matrix:",
                                     crayon::bold("Normalized"),
                                     " (", crayon::red(cli::symbol$cross),
                                     "default)")
        def_trgt_mtrx_used <- FALSE
      }

      if (k == 3) {
        k_used <- crayon::blue("    ",cli::symbol$bullet, "Power used:",
                                       crayon::bold(k),
                                       " (", crayon::green(cli::symbol$tick),
                                       "default)")
        def_k_used <- TRUE
      } else {
        k_used <- crayon::blue("    ",cli::symbol$bullet, "Power used:",
                                       crayon::bold(k),
                                       " (", crayon::red(cli::symbol$cross),
                                       "default)")
        def_k_used <- FALSE
      }


    } else if (type == "psych") {

      # test if all type psych arguments are on default

      if (init_comm == "smc") {
        i_c_used <- crayon::blue("    ",cli::symbol$bullet, "Initial Communalities:",
                                 crayon::bold(init_comm),
                                 " (", crayon::green(cli::symbol$tick),
                                 "default)")
        def_i_c_used <- TRUE
      } else {
        i_c_used <- crayon::blue("    ",cli::symbol$bullet, "Initial Communalities:",
                                 crayon::bold(init_comm),
                                 " (", crayon::red(cli::symbol$cross),
                                 "default)")
        def_i_c_used <- FALSE
      }

      if (criterion == 1e-3) {
        crit_used <- crayon::blue("    ",cli::symbol$bullet, "Convergence Criterion:",
                                  crayon::bold(criterion),
                                  " (", crayon::green(cli::symbol$tick),
                                  "default)")
        def_crit_used <- TRUE
      } else {
        crit_used <- crayon::blue("    ",cli::symbol$bullet, "Initial Communalities:",
                                  crayon::bold(criterion),
                                  " (", crayon::red(cli::symbol$cross),
                                  "default)")
        def_crit_used <- FALSE
      }


      if (criterion_type == "sums") {
        crit_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Criterion Type:",
                                     crayon::bold(criterion_type),
                                     " (", crayon::green(cli::symbol$tick),
                                     "default)")
        def_crit_tp_used <- TRUE
      } else {
        crit_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Criterion Type:",
                                     crayon::bold(criterion_type),
                                     " (", crayon::red(cli::symbol$cross),
                                     "default)")
        def_crit_tp_used <- FALSE
      }

      if (max_iter == 50) {
        mx_it_used <- crayon::blue("    ",cli::symbol$bullet, "Maximum Iteration:",
                                   crayon::bold(max_iter),
                                   " (", crayon::green(cli::symbol$tick),
                                   "default)")
        def_mx_it_used <- TRUE
      } else {
        mx_it_used <- crayon::blue("    ",cli::symbol$bullet, "Maximum Iteration:",
                                   crayon::bold(max_iter),
                                   " (", crayon::red(cli::symbol$cross),
                                   "default)")
        def_mx_it_used <- FALSE
      }

      if (isFALSE(abs_eigen)) {
        abs_eigen_used <- crayon::blue("    ",cli::symbol$bullet, "Absolute Eigenvalues:",
                                       crayon::bold("No"),
                                       " (", crayon::green(cli::symbol$tick),
                                       "default)")
        def_abs_eigen_used <- TRUE
      } else {
        abs_eigen_used <- crayon::blue("    ",cli::symbol$bullet, "Absolute Eigenvalues:",
                                       crayon::bold("Yes"),
                                       " (", crayon::red(cli::symbol$cross),
                                       "default)")
        def_abs_eigen_used <- FALSE
      }

      # Varimax settings

      if (isTRUE(kaiser)) {
        kaiser_used <- crayon::blue("    ",cli::symbol$bullet, "Kaiser Normalization:",
                                    crayon::bold("Yes"),
                                    " (", crayon::green(cli::symbol$tick),
                                    "default)")
        def_kaiser_used <- TRUE
      } else {
        kaiser_used <- crayon::blue("    ",cli::symbol$bullet, "Kaiser Normalization:",
                                    crayon::bold("No"),
                                    " (", crayon::red(cli::symbol$cross),
                                    "default)")
        def_kaiser_used <- FALSE
      }

      if (precision == 1e-5) {
        prcsn_used <- crayon::blue("    ",cli::symbol$bullet, "Precision:",
                                   crayon::bold(precision),
                                   " (", crayon::green(cli::symbol$tick),
                                   "default)")
        def_prcsn_used <- TRUE
      } else {
        prcsn_used <- crayon::blue("    ",cli::symbol$bullet, "Precision:",
                                   crayon::bold(precision),
                                   " (", crayon::red(cli::symbol$cross),
                                   "default)")
        def_prcsn_used <- FALSE
      }

      if (order_type == "eigen") {
        ordr_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Factors ordered according to:",
                                     crayon::bold("Eigenvalues"),
                                     " (", crayon::green(cli::symbol$tick),
                                     "default)")
        def_ordr_tp_used <- TRUE
      } else {
        ordr_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Factors ordered according to:",
                                     crayon::bold("Sum of Squared Loadings"),
                                     " (", crayon::red(cli::symbol$cross),
                                     "default)")
        def_ordr_tp_used <- FALSE
      }

      # Promax settings

      if (P_type == "unnorm") {
        trgt_mtrx_used <- crayon::blue("    ",cli::symbol$bullet, "Target matrix:",
                                       crayon::bold("Unnormalized"),
                                       " (", crayon::green(cli::symbol$tick),
                                       "default)")
        def_trgt_mtrx_used <- TRUE
      } else {
        trgt_mtrx_used <- crayon::blue("    ",cli::symbol$bullet, "Target matrix:",
                                       crayon::bold("Normalized"),
                                       " (", crayon::red(cli::symbol$cross),
                                       "default)")
        def_trgt_mtrx_used <- FALSE
      }

      if (k == 4) {
        k_used <- crayon::blue("    ",cli::symbol$bullet, "Power used:",
                               crayon::bold(k),
                               " (", crayon::green(cli::symbol$tick),
                               "default)")
        def_k_used <- TRUE
      } else {
        k_used <- crayon::blue("    ",cli::symbol$bullet, "Power used:",
                               crayon::bold(k),
                               " (", crayon::red(cli::symbol$cross),
                               "default)")
        def_k_used <- FALSE
      }


    } else if (type == "SPSS") {

      # test if all type SPSS arguments are on default

      if (init_comm == "smc") {
        i_c_used <- crayon::blue("    ",cli::symbol$bullet, "Initial Communalities:",
                                 crayon::bold(init_comm),
                                 " (", crayon::green(cli::symbol$tick),
                                 "default)")
        def_i_c_used <- TRUE
      } else {
        i_c_used <- crayon::blue("    ",cli::symbol$bullet, "Initial Communalities:",
                                 crayon::bold(init_comm),
                                 " (", crayon::red(cli::symbol$cross),
                                 "default)")
        def_i_c_used <- FALSE
      }

      if (criterion == 1e-3) {
        crit_used <- crayon::blue("    ",cli::symbol$bullet, "Convergence Criterion:",
                                  crayon::bold(criterion),
                                  " (", crayon::green(cli::symbol$tick),
                                  "default)")
        def_crit_used <- TRUE
      } else {
        crit_used <- crayon::blue("    ",cli::symbol$bullet, "Initial Communalities:",
                                  crayon::bold(criterion),
                                  " (", crayon::red(cli::symbol$cross),
                                  "default)")
        def_crit_used <- FALSE
      }


      if (criterion_type == "max_individual") {
        crit_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Criterion Type:",
                                     crayon::bold(criterion_type),
                                     " (", crayon::green(cli::symbol$tick),
                                     "default)")
        def_crit_tp_used <- TRUE
      } else {
        crit_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Criterion Type:",
                                     crayon::bold(criterion_type),
                                     " (", crayon::red(cli::symbol$cross),
                                     "default)")
        def_crit_tp_used <- FALSE
      }

      if (max_iter == 25) {
        mx_it_used <- crayon::blue("    ",cli::symbol$bullet, "Maximum Iteration:",
                                   crayon::bold(max_iter),
                                   " (", crayon::green(cli::symbol$tick),
                                   "default)")
        def_mx_it_used <- TRUE
      } else {
        mx_it_used <- crayon::blue("    ",cli::symbol$bullet, "Maximum Iteration:",
                                   crayon::bold(max_iter),
                                   " (", crayon::red(cli::symbol$cross),
                                   "default)")
        def_mx_it_used <- FALSE
      }

      if (isTRUE(abs_eigen)) {
        abs_eigen_used <- crayon::blue("    ",cli::symbol$bullet, "Absolute Eigenvalues:",
                                       crayon::bold("Yes"),
                                       " (", crayon::green(cli::symbol$tick),
                                       "default)")
        def_abs_eigen_used <- TRUE
      } else {
        abs_eigen_used <- crayon::blue("    ",cli::symbol$bullet, "Absolute Eigenvalues:",
                                       crayon::bold("No"),
                                       " (", crayon::red(cli::symbol$cross),
                                       "default)")
        def_abs_eigen_used <- FALSE
      }

      # Varimax settings

      if (isTRUE(kaiser)) {
        kaiser_used <- crayon::blue("    ",cli::symbol$bullet, "Kaiser Normalization:",
                                    crayon::bold("Yes"),
                                    " (", crayon::green(cli::symbol$tick),
                                    "default)")
        def_kaiser_used <- TRUE
      } else {
        kaiser_used <- crayon::blue("    ",cli::symbol$bullet, "Kaiser Normalization:",
                                    crayon::bold("No"),
                                    " (", crayon::red(cli::symbol$cross),
                                    "default)")
        def_kaiser_used <- FALSE
      }

      if (precision == 1e-10) {
        prcsn_used <- crayon::blue("    ",cli::symbol$bullet, "Precision:",
                                   crayon::bold(precision),
                                   " (", crayon::green(cli::symbol$tick),
                                   "default)")
        def_prcsn_used <- TRUE
      } else {
        prcsn_used <- crayon::blue("    ",cli::symbol$bullet, "Precision:",
                                   crayon::bold(precision),
                                   " (", crayon::red(cli::symbol$cross),
                                   "default)")
        def_prcsn_used <- FALSE
      }

      if (order_type == "ss_factors") {
        ordr_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Factors ordered according to:",
                                     crayon::bold("Sum of Squared Loadings"),
                                     " (", crayon::green(cli::symbol$tick),
                                     "default)")
        def_ordr_tp_used <- TRUE
      } else {
        ordr_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Factors ordered according to:",
                                     crayon::bold("Eigenvalues"),
                                     " (", crayon::red(cli::symbol$cross),
                                     "default)")
        def_ordr_tp_used <- FALSE
      }

      # Promax settings

      if (P_type == "norm") {
        trgt_mtrx_used <- crayon::blue("    ",cli::symbol$bullet, "Target matrix:",
                                       crayon::bold("Normalized"),
                                       " (", crayon::green(cli::symbol$tick),
                                       "default)")
        def_trgt_mtrx_used <- TRUE
      } else {
        trgt_mtrx_used <- crayon::blue("    ",cli::symbol$bullet, "Target matrix:",
                                       crayon::bold("Unnormalized"),
                                       " (", crayon::red(cli::symbol$cross),
                                       "default)")
        def_trgt_mtrx_used <- FALSE
      }

      if (k == 4) {
        k_used <- crayon::blue("    ",cli::symbol$bullet, "Power used:",
                               crayon::bold(k),
                               " (", crayon::green(cli::symbol$tick),
                               "default)")
        def_k_used <- TRUE
      } else {
        k_used <- crayon::blue("    ",cli::symbol$bullet, "Power used:",
                               crayon::bold(k),
                               " (", crayon::red(cli::symbol$cross),
                               "default)")
        def_k_used <- FALSE
      }

    } else {

      # if no valid type is specified

      # if no type is specified

      i_c_used <- crayon::blue("    ",cli::symbol$bullet, "Initial Communalities:",
                               crayon::bold(init_comm))
      def_i_c_used <- FALSE

      crit_used <- crayon::blue("    ",cli::symbol$bullet, "Convergence Criterion:",
                                crayon::bold(criterion))
      def_crit_used <- FALSE


      crit_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Criterion Type:",
                                   crayon::bold(criterion_type))
      def_crit_tp_used <- FALSE

      mx_it_used <- crayon::blue("    ",cli::symbol$bullet, "Maximum Iteration:",
                                 crayon::bold(max_iter))
      def_mx_it_used <- FALSE


      if (isTRUE(abs_eigen)) {
        abs_eigen_used <- crayon::blue("    ",cli::symbol$bullet, "Absolute Eigenvalues:",
                                       crayon::bold("Yes"))
        def_abs_eigen_used <- FALSE
      } else {
        abs_eigen_used <- crayon::blue("    ",cli::symbol$bullet, "Absolute Eigenvalues:",
                                       crayon::bold("No"))
        def_abs_eigen_used <- FALSE
      }

      # Varimax settings

      if (isTRUE(kaiser)) {
        kaiser_used <- crayon::blue("    ",cli::symbol$bullet, "Kaiser Normalization:",
                                    crayon::bold("Yes"))
        def_kaiser_used <- FALSE
      } else {
        kaiser_used <- crayon::blue("    ",cli::symbol$bullet, "Kaiser Normalization:",
                                    crayon::bold("No"))
        def_kaiser_used <- FALSE
      }


      prcsn_used <- crayon::blue("    ",cli::symbol$bullet, "Precision:",
                                 crayon::bold(precision))
      def_prcsn_used <- FALSE


      if (order_type == "ss_factors") {
        ordr_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Factors ordered according to:",
                                     crayon::bold("Sum of Squared Loadings"))
        def_ordr_tp_used <- FALSE
      } else {
        ordr_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Factors ordered according to:",
                                     crayon::bold("Eigenvalues"))
        def_ordr_tp_used <- FALSE
      }

      # Promax settings

      if (P_type == "norm") {
        trgt_mtrx_used <- crayon::blue("    ",cli::symbol$bullet, "Target matrix:",
                                       crayon::bold("Normalized"))
        def_trgt_mtrx_used <- FALSE
      } else {
        trgt_mtrx_used <- crayon::blue("    ",cli::symbol$bullet, "Target matrix:",
                                       crayon::bold("Unnormalized"))
        def_trgt_mtrx_used <- FALSE
      }

        k_used <- crayon::blue("    ",cli::symbol$bullet, "Power used:",
                               crayon::bold(k))
        def_k_used <- FALSE

    }

    if (all(def_i_c_used, def_crit_used, def_crit_tp_used, def_mx_it_used,
            def_abs_eigen_used)) {

      def_msg <- crayon::green("    ", cli::symbol$tick, "Type", type,
                               "with default", "arguments used.")

    } else {
      def_msg <- crayon::red("    ", cli::symbol$cross, "Type", type,
                             "with customized", "arguments used.")
    }


    if (all(def_kaiser_used, def_prcsn_used,
            def_ordr_tp_used)) {

      def_msg_vrmx <- crayon::green("    ", cli::symbol$tick, "Type", type,
                                    "with default", "arguments used.")

    } else {
      def_msg_vrmx <- crayon::red("    ", cli::symbol$cross, "Type", type,
                                  "with customized", "arguments used.")
    }

    if (all(def_k_used, def_trgt_mtrx_used)) {

      def_msg_prmx <- crayon::green("    ", cli::symbol$tick, "Type", type,
                                    "with default", "arguments used.")

    } else {
      def_msg_prmx <- crayon::red("    ", cli::symbol$cross, "Type", type,
                                  "with customized", "arguments used.")
    }


    # print the different settings, the loadings, and the variances
    cat(crayon::blue$bold("Settings:"))
    cat("\n")
    cat("\n")
    cat(sttngs_intro)
    cat("\n")
    cat(i_c_used)
    cat("\n")
    cat(crit_used)
    cat("\n")
    cat(crit_tp_used)
    cat("\n")
    cat(mx_it_used)
    cat("\n")
    cat(abs_eigen_used)
    cat("\n")
    cat(def_msg)
    cat("\n")
    cat("\n")
    cat(N_used)
    cat("\n")
    if (x$iter > max_iter) {
      cat(conv_warn)
      cat("\n")
    }
    cat("\n")
    cat(sttngs_intro_vrmx)
    cat("\n")
    cat(kaiser_used)
    cat("\n")
    cat(prcsn_used)
    cat("\n")
    cat(ordr_tp_used)
    cat("\n")
    cat(def_msg_vrmx)
    cat("\n")
    cat("\n")
    cat(sttngs_intro_prmx)
    cat("\n")
    cat(trgt_mtrx_used)
    cat("\n")
    cat(k_used)
    cat("\n")
    cat(def_msg_prmx)
    cat("\n")
    cat("\n")
    cat(crayon::blue$bold("Rotated Loadings:"))
    cat("\n")
    cat("\n")
    print(x$rot_loadings)
    cat("\n")
    cat(crayon::blue$bold("Variances Accounted For:"))
    cat("\n")
    cat("\n")
    cat(.get_compare_matrix(vrs_acc, r_red = Inf, n_char = 17))

  } else {

    # Varimax settings Intro message
    sttngs_intro_vrmx <- crayon::blue("Varimax with type =", crayon::bold(type),
                                      "performed with the",
                                      "following settings:")

    # Promax settings Intro message
    sttngs_intro_prmx <- crayon::blue("Promax with type =", crayon::bold(type),
                                      "performed with the",
                                      "following settings:")

    if (type == "GS") {
      # test if all type GS arguments are on default


      # Varimax settings
      if (isTRUE(kaiser)) {
        kaiser_used <- crayon::blue("    ",cli::symbol$bullet, "Kaiser Normalization:",
                                    crayon::bold("Yes"),
                                    " (", crayon::green(cli::symbol$tick),
                                    "default)")
        def_kaiser_used <- TRUE
      } else {
        kaiser_used <- crayon::blue("    ",cli::symbol$bullet, "Kaiser Normalization:",
                                    crayon::bold("No"),
                                    " (", crayon::red(cli::symbol$cross),
                                    "default)")
        def_kaiser_used <- FALSE
      }

      if (precision == 1e-10) {
        prcsn_used <- crayon::blue("    ",cli::symbol$bullet, "Precision:",
                                   crayon::bold(precision),
                                   " (", crayon::green(cli::symbol$tick),
                                   "default)")
        def_prcsn_used <- TRUE
      } else {
        prcsn_used <- crayon::blue("    ",cli::symbol$bullet, "Precision:",
                                   crayon::bold(precision),
                                   " (", crayon::red(cli::symbol$cross),
                                   "default)")
        def_prcsn_used <- FALSE
      }

      if (order_type == "eigen") {
        ordr_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Factors ordered according to:",
                                     crayon::bold("Eigenvalues"),
                                     " (", crayon::green(cli::symbol$tick),
                                     "default)")
        def_ordr_tp_used <- TRUE
      } else {
        ordr_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Factors ordered according to:",
                                     crayon::bold("Sum of Squared Loadings"),
                                     " (", crayon::red(cli::symbol$cross),
                                     "default)")
        def_ordr_tp_used <- FALSE
      }

      # Promax settings

      if (P_type == "unnorm") {
        trgt_mtrx_used <- crayon::blue("    ",cli::symbol$bullet, "Target matrix:",
                                       crayon::bold("Unnormalized"),
                                       " (", crayon::green(cli::symbol$tick),
                                       "default)")
        def_trgt_mtrx_used <- TRUE
      } else {
        trgt_mtrx_used <- crayon::blue("    ",cli::symbol$bullet, "Target matrix:",
                                       crayon::bold("Normalized"),
                                       " (", crayon::red(cli::symbol$cross),
                                       "default)")
        def_trgt_mtrx_used <- FALSE
      }

      if (k == 3) {
        k_used <- crayon::blue("    ",cli::symbol$bullet, "Power used:",
                               crayon::bold(k),
                               " (", crayon::green(cli::symbol$tick),
                               "default)")
        def_k_used <- TRUE
      } else {
        k_used <- crayon::blue("    ",cli::symbol$bullet, "Power used:",
                               crayon::bold(k),
                               " (", crayon::red(cli::symbol$cross),
                               "default)")
        def_k_used <- FALSE
      }



    } else if (type == "psych") {

      # test if all type psych arguments are on default


      # Varimax settings

      if (isTRUE(kaiser)) {
        kaiser_used <- crayon::blue("    ",cli::symbol$bullet, "Kaiser Normalization:",
                                    crayon::bold("Yes"),
                                    " (", crayon::green(cli::symbol$tick),
                                    "default)")
        def_kaiser_used <- TRUE
      } else {
        kaiser_used <- crayon::blue("    ",cli::symbol$bullet, "Kaiser Normalization:",
                                    crayon::bold("No"),
                                    " (", crayon::red(cli::symbol$cross),
                                    "default)")
        def_kaiser_used <- FALSE
      }

      if (precision == 1e-5) {
        prcsn_used <- crayon::blue("    ",cli::symbol$bullet, "Precision:",
                                   crayon::bold(precision),
                                   " (", crayon::green(cli::symbol$tick),
                                   "default)")
        def_prcsn_used <- TRUE
      } else {
        prcsn_used <- crayon::blue("    ",cli::symbol$bullet, "Precision:",
                                   crayon::bold(precision),
                                   " (", crayon::red(cli::symbol$cross),
                                   "default)")
        def_prcsn_used <- FALSE
      }

      if (order_type == "eigen") {
        ordr_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Factors ordered according to:",
                                     crayon::bold("Eigenvalues"),
                                     " (", crayon::green(cli::symbol$tick),
                                     "default)")
        def_ordr_tp_used <- TRUE
      } else {
        ordr_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Factors ordered according to:",
                                     crayon::bold("Sum of Squared Loadings"),
                                     " (", crayon::red(cli::symbol$cross),
                                     "default)")
        def_ordr_tp_used <- FALSE
      }

      # Promax settings

      if (P_type == "unnorm") {
        trgt_mtrx_used <- crayon::blue("    ",cli::symbol$bullet, "Target matrix:",
                                       crayon::bold("Unnormalized"),
                                       " (", crayon::green(cli::symbol$tick),
                                       "default)")
        def_trgt_mtrx_used <- TRUE
      } else {
        trgt_mtrx_used <- crayon::blue("    ",cli::symbol$bullet, "Target matrix:",
                                       crayon::bold("Normalized"),
                                       " (", crayon::red(cli::symbol$cross),
                                       "default)")
        def_trgt_mtrx_used <- FALSE
      }

      if (k == 4) {
        k_used <- crayon::blue("    ",cli::symbol$bullet, "Power used:",
                               crayon::bold(k),
                               " (", crayon::green(cli::symbol$tick),
                               "default)")
        def_k_used <- TRUE
      } else {
        k_used <- crayon::blue("    ",cli::symbol$bullet, "Power used:",
                               crayon::bold(k),
                               " (", crayon::red(cli::symbol$cross),
                               "default)")
        def_k_used <- FALSE
      }



    } else if (type == "SPSS") {

      # test if all type SPSS arguments are on default

      # Varimax settings

      if (isTRUE(kaiser)) {
        kaiser_used <- crayon::blue("    ",cli::symbol$bullet, "Kaiser Normalization:",
                                    crayon::bold("Yes"),
                                    " (", crayon::green(cli::symbol$tick),
                                    "default)")
        def_kaiser_used <- TRUE
      } else {
        kaiser_used <- crayon::blue("    ",cli::symbol$bullet, "Kaiser Normalization:",
                                    crayon::bold("No"),
                                    " (", crayon::red(cli::symbol$cross),
                                    "default)")
        def_kaiser_used <- FALSE
      }

      if (precision == 1e-10) {
        prcsn_used <- crayon::blue("    ",cli::symbol$bullet, "Precision:",
                                   crayon::bold(precision),
                                   " (", crayon::green(cli::symbol$tick),
                                   "default)")
        def_prcsn_used <- TRUE
      } else {
        prcsn_used <- crayon::blue("    ",cli::symbol$bullet, "Precision:",
                                   crayon::bold(precision),
                                   " (", crayon::red(cli::symbol$cross),
                                   "default)")
        def_prcsn_used <- FALSE
      }

      if (order_type == "ss_factors") {
        ordr_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Factors ordered according to:",
                                     crayon::bold("Sum of Squared Loadings"),
                                     " (", crayon::green(cli::symbol$tick),
                                     "default)")
        def_ordr_tp_used <- TRUE
      } else {
        ordr_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Factors ordered according to:",
                                     crayon::bold("Eigenvalues"),
                                     " (", crayon::red(cli::symbol$cross),
                                     "default)")
        def_ordr_tp_used <- FALSE
      }

      # Promax settings

      if (P_type == "norm") {
        trgt_mtrx_used <- crayon::blue("    ",cli::symbol$bullet, "Target matrix:",
                                       crayon::bold("Normalized"),
                                       " (", crayon::green(cli::symbol$tick),
                                       "default)")
        def_trgt_mtrx_used <- TRUE
      } else {
        trgt_mtrx_used <- crayon::blue("    ",cli::symbol$bullet, "Target matrix:",
                                       crayon::bold("Unnormalized"),
                                       " (", crayon::red(cli::symbol$cross),
                                       "default)")
        def_trgt_mtrx_used <- FALSE
      }

      if (k == 4) {
        k_used <- crayon::blue("    ",cli::symbol$bullet, "Power used:",
                               crayon::bold(k),
                               " (", crayon::green(cli::symbol$tick),
                               "default)")
        def_k_used <- TRUE
      } else {
        k_used <- crayon::blue("    ",cli::symbol$bullet, "Power used:",
                               crayon::bold(k),
                               " (", crayon::red(cli::symbol$cross),
                               "default)")
        def_k_used <- FALSE
      }

    } else {
      # Varimax settings

      if (isTRUE(kaiser)) {
        kaiser_used <- crayon::blue("    ",cli::symbol$bullet, "Kaiser Normalization:",
                                    crayon::bold("Yes"))
        def_kaiser_used <- FALSE
      } else {
        kaiser_used <- crayon::blue("    ",cli::symbol$bullet, "Kaiser Normalization:",
                                    crayon::bold("No"))
        def_kaiser_used <- FALSE
      }


      prcsn_used <- crayon::blue("    ",cli::symbol$bullet, "Precision:",
                                 crayon::bold(precision))
      def_prcsn_used <- FALSE


      if (order_type == "ss_factors") {
        ordr_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Factors ordered according to:",
                                     crayon::bold("Sum of Squared Loadings"))
        def_ordr_tp_used <- FALSE
      } else {
        ordr_tp_used <- crayon::blue("    ",cli::symbol$bullet, "Factors ordered according to:",
                                     crayon::bold("Eigenvalues"))
        def_ordr_tp_used <- FALSE
      }

      # Promax settings

      if (P_type == "norm") {
        trgt_mtrx_used <- crayon::blue("    ",cli::symbol$bullet, "Target matrix:",
                                       crayon::bold("Normalized"))
        def_trgt_mtrx_used <- FALSE
      } else {
        trgt_mtrx_used <- crayon::blue("    ",cli::symbol$bullet, "Target matrix:",
                                       crayon::bold("Unnormalized"))
        def_trgt_mtrx_used <- FALSE
      }

      k_used <- crayon::blue("    ",cli::symbol$bullet, "Power used:",
                             crayon::bold(k))
      def_k_used <- FALSE
    }


    if (all(def_kaiser_used, def_prcsn_used,
            def_ordr_tp_used)) {

      def_msg_vrmx <- crayon::green("    ", cli::symbol$tick, "Type", type,
                                    "with default", "arguments used.")

    } else {
      def_msg_vrmx <- crayon::red("    ", cli::symbol$cross, "Type", type,
                                  "with customized", "arguments used.")
    }

    if (all(def_k_used, def_trgt_mtrx_used)) {

      def_msg_prmx <- crayon::green("    ", cli::symbol$tick, "Type", type,
                                    "with default", "arguments used.")

    } else {
      def_msg_prmx <- crayon::red("    ", cli::symbol$cross, "Type", type,
                                  "with customized", "arguments used.")
    }


    # print the different settings, the loadings, and the variances
    cat(crayon::blue$bold("Settings:"))
    cat("\n")
    cat("\n")
    cat(sttngs_intro_vrmx)
    cat("\n")
    cat(kaiser_used)
    cat("\n")
    cat(prcsn_used)
    cat("\n")
    cat(ordr_tp_used)
    cat("\n")
    cat(def_msg_vrmx)
    cat("\n")
    cat("\n")
    cat(sttngs_intro_prmx)
    cat("\n")
    cat(trgt_mtrx_used)
    cat("\n")
    cat(k_used)
    cat("\n")
    cat(def_msg_prmx)
    cat("\n")
    cat("\n")
    cat(crayon::blue$bold("Rotated Loadings:"))
    cat("\n")
    cat("\n")
    print(x$rot_loadings)
    cat("\n")
    cat(crayon::blue$bold("Variances Accounted For:"))
    cat("\n")
    cat("\n")
    cat(.get_compare_matrix(vrs_acc, r_red = Inf, n_char = 17))

  }
}
