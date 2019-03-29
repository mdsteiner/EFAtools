#' Print SL object
#'
#' Print Method showing a summarized output of the SL procedure.
#'
#' @param x list. An object of class SL to be printed
#' @param ...  Further arguments for print.
#'
#' @export
#' @method print SL
#'
#' @examples
#' EFA_SG_5 <- EFA(IDS2_R, n_factors = 5, type = "SG", rotation = "promax")
#' SL(EFA_SG_5)
print.SL <- function(x, ...) {

  # extract the different settings
  max_iter <- x$settings$max_iter
  type <- x$settings$type
  init_comm <- x$settings$init_comm
  criterion <- x$settings$criterion
  criterion_type <- x$settings$criterion_type
  abs_eigen <- x$settings$abs_eigen
  signed_loadings <- x$settings$signed_loadings
  use <- x$settings$use

  # get the loadings
  ldngs <- x$sl

  # get the variances accounted
  vrs_acc <- x$vars_accounted

  if (x$iter > max_iter) {
    conv_warn <- crayon::red$bold(cli::symbol$cross,
                                  "Maximum number of iterations reached",
                                  "without convergence")
  }


  # Settings Intro message
  sttngs_intro <- crayon::blue("PAF for second order loadings with type =",
                               crayon::bold(type),
                               "performed with the",
                               "following settings:")

  if (type == "SG") {
    # test if all type SG arguments are on default

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

    if (criterion == 1e-6) {
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

  }

  if (all(def_i_c_used, def_crit_used, def_crit_tp_used, def_mx_it_used,
          def_abs_eigen_used)) {

    def_msg <- crayon::green("    ", cli::symbol$tick, "Type", type,
                             "with default", "arguments used.")

  } else {
    def_msg <- crayon::red("    ", cli::symbol$cross, "Type", type,
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
  if (x$iter > max_iter) {
    cat(conv_warn)
    cat("\n")
  }
  cat("\n")
  cat(crayon::blue$bold("Schmid-Leiman Solution:"))
  cat("\n")
  cat("\n")
  print(ldngs)
  cat("\n")
  cat(crayon::blue$bold("Variances Accounted For:"))
  cat("\n")
  cat("\n")
  cat(.get_compare_matrix(vrs_acc, r_red = Inf, n_char = 17))

}