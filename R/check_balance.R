#' @title Check the variable balance between the probability and non-probability samples
#'
#' @param x Formula specifying variables to check
#' @param object Object of `nonprobsvy` class
#' @param dig Number of digits for rounding (default = 2)
#'
#' @importFrom stats aggregate
#' @importFrom survey svytotal
#' @importFrom stats setNames
#'
#' @return A `list` containing nonprobability totals, probability totals, and their differences
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
#' check_balance(~size, ipw_est1)
#'
#' ## check the balance for the calibrated IPW
#' check_balance(~size, ipw_est2)
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

  if (missing(object) || is.null(object)) {
    stop("The `object` argument is required.")
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

  # Function to calculate totals for one variable
  calculate_totals <- function(var, data) {
    # Check for NAs
    if (sum(is.na(data[[var]])) > 0) {
      warning(sprintf("NA values found in variable %s.", var))
    }

    # For categorical variables, handle each level
    if (is.factor(data[[var]]) || is.character(data[[var]])) {
      levels <- unique(data[[var]][!is.na(data[[var]])])
      totals <- sapply(levels, function(lvl) {
        data_subset <- data[data[[var]] == lvl, ]
        if (nrow(data_subset) < 5) {
          warning(sprintf("Small group size (< 5) for level %s in variable %s.", lvl, var))
        }
        sum(data_subset$ipw_weights)
      })
      names(totals) <- paste0(var, levels)
      return(totals)
    } else {
      # For numeric variables
      return(setNames(sum(data$ipw_weights * data[[var]], na.rm = TRUE), var))
    }
  }

  # Prepare data
  data <- object$data
  data$ipw_weights <- object$ipw_weights

  # Calculate nonprob totals
  nonprob_totals <- tryCatch(
    {
      unlist(lapply(vars, function(var) calculate_totals(var, data)))
    },
    error = function(e) {
      stop(sprintf("Error calculating nonprobability totals: %s.", e$message))
    }
  )

  # Calculate probability totals
  if (!is.null(object$svydesign)) {
    # If svydesign exists
    prob_totals <- tryCatch(
      {
        svy_total <- svytotal(x, object$svydesign)
        svy_totals <- as.vector(svy_total)
        names(svy_totals) <- names(svy_total)
        svy_totals
      },
      error = function(e) {
        stop(sprintf("Error calculating survey totals: %s.", e$message))
      }
    )
  } else {
    # Use population totals
    prob_totals <- object$selection$pop_totals
    if (is.null(prob_totals)) {
      stop("The `pop_totals` argument is null.")
    }
  }

  # Calculate and round differences
  diff <- tryCatch(
    {
      round(nonprob_totals - prob_totals[names(nonprob_totals)], digits = dig)
    },
    error = function(e) {
      stop(sprintf("Error calculating differences: %s.", e$message))
    }
  )

  # Return results with meaningful names
  result <- list(
    nonprob_totals = nonprob_totals,
    prob_totals = prob_totals,
    balance = diff
  )
  return(result)
}
