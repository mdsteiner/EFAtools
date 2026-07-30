#' Print and format an efa_simulated object
#'
#' `print()` shows a compact summary of the data simulated by [efa_simulate()]: how many
#' datasets were drawn and their dimensions, the marginal distribution, whether the data were
#' discretized into ordered categories or given missing values, and -- when model error was
#' injected -- the method with the target and achieved RMSEA and CFI. The simulated data
#' themselves live in the `data` element and the population correlation matrix in `population`.
#' `format()` returns the same summary as a character vector; `print()` is
#' `cat(format(x), sep = "\n")`. The lines follow the active console theme, so they are plain
#' when colours are disabled.
#'
#' @param x An object of class `efa_simulated` (output from [efa_simulate()]).
#' @param digits Integer. The number of decimal places the reported fit values are rounded to.
#'   Default is 3.
#' @param ... Not used; for consistency with the generic.
#'
#' @returns `print()` returns its argument `x` invisibly. `format()` returns a character vector
#'   with the summary lines.
#'
#' @family data simulation
#'
#' @export
#'
#' @method print efa_simulated
#'
#' @examples
#' Lambda <- population_models$loadings$baseline
#' Phi <- population_models$phis_3$moderate
#' efa_simulate(N = 500, Lambda = Lambda, Phi = Phi, target_rmsea = 0.05, seed = 42)
#'
print.efa_simulated <- function(x, digits = 3, ...) {
  cat(format(x, digits = digits, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.efa_simulated
#' @export
#' @method format efa_simulated
format.efa_simulated <- function(x, digits = 3, ...) {

  set <- x$settings
  p <- ncol(x$population)

  cli::cli_format_method({

    cli::cli_rule(left = "{.strong Simulated data}")
    cli::cli_text("")

    if (is.null(x$data)) {
      cli::cli_text("{cli::qty(p)}Population correlation matrix only ({p} variable{?s}); no data drawn.")
    } else {
      one <- if (is.list(x$data)) x$data[[1L]] else x$data
      nd <- set$n_datasets
      n <- nrow(one)
      kind <- if (is.integer(one)) "ordered-category" else "continuous"
      if (nd > 1L) {
        cli::cli_text("{cli::qty(nd)}{nd} dataset{?s} of {cli::qty(n)}{n} case{?s} by {cli::qty(p)}{p} variable{?s} ({kind}).")
      } else {
        cli::cli_text("{cli::qty(n)}{n} case{?s} by {cli::qty(p)}{p} variable{?s} ({kind}).")
      }
    }

    marg <- switch(set$marginals,
                   normal = "normal",
                   empirical = "empirical",
                   VM = "Vale-Maurelli non-normal",
                   IG = "independent-generator non-normal")
    cli::cli_text("Marginals: {marg}.")
    if (!is.null(set$categories)) {
      cli::cli_text("Discretized into ordered categories.")
    }
    if (!is.null(set$missing) && set$missing != "none") {
      cli::cli_text("Missing data: {set$missing}.")
    }

    cli::cli_text("")
    me <- x$model_error
    if (is.null(me)) {
      cli::cli_text("Model error: none (exact population).")
    } else {
      cli::cli_text("Model error: {me$method}")
      tr <- if (is.null(me$target_rmsea)) "--" else .efa_num(me$target_rmsea, digits, pad = FALSE)
      tc <- if (is.null(me$target_cfi)) "--" else .efa_num(me$target_cfi, digits, pad = FALSE)
      ar <- .efa_num(me$rmsea, digits, pad = FALSE)
      ac <- .efa_num(me$cfi, digits, pad = FALSE)
      cli::cli_text("  RMSEA: target {tr}, achieved {ar}")
      cli::cli_text("  CFI:   target {tc}, achieved {ac}")
    }

  })
}
