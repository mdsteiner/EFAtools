# Declarative `type` presets and the single resolver that applies them.
#
# Every estimator/rotation helper exposes a `type` preset ("EFAtools", "psych",
# "SPSS") that fills a block of tuning arguments, plus "none" which requires the
# user to set them explicitly. `.resolve_settings()` is the one place that turns
# a chosen `type` and the user-supplied (possibly NA) arguments into the final
# settings, warning once about any argument pinned alongside a preset.

# Per-function argument defaults for each preset. One source of truth, consumed
# by .PAF() (via .estimate_model()) and by .rotate_model() for every rotation.
.efa_presets <- list(
  PAF = list(
    EFAtools = list(init_comm = "smc", criterion = 1e-3,
                    criterion_type = "sum", max_iter = 300, abs_eigen = TRUE),
    psych    = list(init_comm = "smc", criterion = 1e-3,
                    criterion_type = "sum", max_iter = 50, abs_eigen = FALSE),
    SPSS     = list(init_comm = "smc", criterion = 1e-3,
                    criterion_type = "max_individual", max_iter = 25,
                    abs_eigen = TRUE)
  ),
  PROMAX = list(
    EFAtools = list(normalize = TRUE, P_type = "norm", order_type = "eigen",
                    varimax_type = "kaiser", k = 4),
    psych    = list(normalize = TRUE, P_type = "unnorm", order_type = "eigen",
                    varimax_type = "svd", k = 4),
    SPSS     = list(normalize = TRUE, P_type = "norm", order_type = "ss_factors",
                    varimax_type = "kaiser", k = 4)
  ),
  VARIMAX = list(
    EFAtools = list(normalize = TRUE, order_type = "eigen",
                    varimax_type = "kaiser"),
    psych    = list(normalize = TRUE, order_type = "eigen",
                    varimax_type = "svd"),
    SPSS     = list(normalize = TRUE, order_type = "ss_factors",
                    varimax_type = "kaiser")
  ),
  ROTATE_OBLQ = list(
    EFAtools = list(normalize = TRUE, order_type = "eigen"),
    psych    = list(normalize = TRUE, order_type = "eigen"),
    SPSS     = list(normalize = TRUE, order_type = "ss_factors")
  ),
  ROTATE_ORTH = list(
    EFAtools = list(normalize = TRUE, order_type = "eigen"),
    psych    = list(normalize = TRUE, order_type = "eigen"),
    SPSS     = list(normalize = TRUE, order_type = "ss_factors")
  )
)

# Resolve a preset block into the final settings list. `type` is the chosen
# preset (or "none"); `user` holds the user-supplied values (NA where unset);
# `preset` is the `.efa_presets` slice for the calling function. Every setting
# uses NA to signal "unset" except `normalize`, the one on/off toggle whose
# preset value is always TRUE, so the user only overrides it by passing FALSE.
.resolve_settings <- function(type, user, preset) {

  setting_names <- names(preset[[1]])
  required <- setdiff(setting_names, "normalize")

  if (type == "none") {

    # Without a preset every NA-defaulted argument must be supplied by the user.
    missing <- required[vapply(user[required], is.na, logical(1))]

    if (length(missing) > 0) {
      cli::cli_abort(
        c(
          paste("{cli::qty(missing)} {.arg type} is {.val none} but",
                "{?a required argument is/required arguments are} unset:",
                "{.arg {missing}}."),
          "i" = paste("Either set {.arg type} to {.val EFAtools},",
                      "{.val psych}, or {.val SPSS}, or supply",
                      "{cli::qty(missing)} {?it/them} explicitly.")
        ),
        class = "efa_type_none"
      )
    }

    return(user[setting_names])
  }

  defaults <- preset[[type]]
  resolved <- user[setting_names]
  pinned <- character(0)

  for (nm in setting_names) {

    if (identical(nm, "normalize")) {
      if (isFALSE(user[[nm]])) pinned <- c(pinned, nm)
    } else if (is.na(user[[nm]])) {
      resolved[[nm]] <- defaults[[nm]]
    } else {
      pinned <- c(pinned, nm)
    }

  }

  if (length(pinned) > 0) {

    # `pairs` is interpolated as a variable (never baked into the glue template)
    # so user-supplied values cannot be parsed as cli markup.
    pairs <- paste0(pinned, " = ",
                    vapply(pinned, function(nm) format(resolved[[nm]]),
                           character(1)))

    cli::cli_warn(
      c(
        paste("{cli::qty(pinned)} {?An argument was/Arguments were} set",
              "together with {.arg type} = {.val {type}}; the supplied",
              "value{?s} {?is/are} used and may differ from the",
              "{.val {type}} preset:"),
        "*" = "{pairs}"
      ),
      class = "efa_type_override"
    )
  }

  resolved
}
