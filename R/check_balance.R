#' @title Checks the variable balance between the probability and non-probability samples
#'
#' @description
#' Function compares totals for auxiliary variables specified in the `x` argument for an `object` that
#' contains either IPW or DR estimator.
#'
#' @param x formula specifying variables to check
#' @param object object of `nonprob` class
#' @param dig number of digits for rounding (default = 2)
#'
#' @importFrom stats aggregate
#' @importFrom survey svytotal
#' @importFrom stats setNames
#'
#' @return A `list` containing totals for non-probability and probability samples and their differences
#'
#' @examples
#'
#' data(admin)
#' data(jvs)
#'
#' jvs_svy <- svydesign(ids = ~ 1,  weights = ~ weight,
#' strata = ~ size + nace + region, data = jvs)
#'
#' ipw_est1 <- nonprob(selection = ~ region + private + nace + size,
#' target = ~ single_shift,
#' svydesign = jvs_svy,
#' data = admin, method_selection = "logit"
#' )
#'
#' ipw_est2 <- nonprob(
#' selection = ~ region + private + nace + size,
#' target = ~ single_shift,
#' svydesign = jvs_svy,
#' data = admin, method_selection = "logit",
#' control_selection = control_sel(est_method = "gee", gee_h_fun = 1))
#'
#' ## check the balance for the standard IPW
#' check_balance(~size+public, ipw_est1)
#'
#' ## check the balance for the calibrated IPW
#' check_balance(~size+public, ipw_est2)
#'
#' ## check balance for a more complicated example
#' check_balance(~ I(size=="M") + I(nace == "C"), ipw_est1)
#'
#' @export
check_balance <- function(x, object, dig) {
  UseMethod("check_balance", object)
}

#' @method check_balance nonprob
#' @exportS3Method
check_balance.nonprob <- function(x, object, dig = 2) {

  # Input validation
  if (!inherits(x, "formula")) {
    stop("The `x` argument must be a formula.")
  }

  if (missing(object) || is.null(object) || !inherits(object, "nonprob")) {
    stop("The `object` argument of class `nonprob` is required.")
  }

  if (object$estimator == "mi") {
    stop("No estimated weights available. Only the IPW or the DR methods are supported.")
  }

  if (!is.numeric(dig) || dig < 0) {
    stop("The `dig` argument must be a non-negative number")
  }

  if (nrow(object$data) == 0) {
    stop("An empty dataset detected (zero rows).")
  }

  if (sum(object$ipw_weights) == 0) {
    stop("The sum of weights is zero.")
  }

  # Extract variables from formula
  vars <- all.vars(x)
  if (length(vars) == 0) {
    stop("No variables specified in the `formula` argument")
  }

  # Check if all variables exist in the data
  missing_vars <- setdiff(vars, names(object$data))
  if (length(missing_vars) > 0) {
    stop(sprintf(
      "The following variables are not present in the dataset: %s.",
      paste(missing_vars, collapse = ", ")
    ))
  }


  nonprob_totals <- colSums(model.matrix(x, object$data)*weights(object))
  totals_names <- names(nonprob_totals)

  prob_totals <- if (!is.null(object$pop_totals)) object$pop_totals else
    colSums(model.matrix(x, object$svydesign$variables)*weights(object$svydesign))


  result <- list(
    nonprob_totals = nonprob_totals,
    prob_totals = prob_totals[totals_names],
    balance = round(nonprob_totals - prob_totals[totals_names], digits = dig)
  )
  return(result)
}
