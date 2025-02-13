#' @title Total values of covariates in subgroups
#'
#' @param x - formula
#' @param object - object of nonprobsvy class
#' @param interaction - logical, if TRUE calculate for all combinations of grouping variables
#'
#' @return A data frame with estimated totals of the given covariates in subgroups
#' @importFrom formula.tools lhs.vars
#' @importFrom formula.tools rhs.vars
#' @importFrom stats aggregate
#' @importFrom survey svytotal
nonprobsvytotal <- function(x, object, interaction) {
  UseMethod("nonprobsvytotal", object)
}
# Internal function - not exported in CRAN version
# Will be exported in future releases after variance estimation is implemented
#' @keywords internal
nonprobsvymean.nonprobsvy <- function(x, object, interaction = FALSE) {
  # TODO variance - NOT READY FOR CRAN
  groups <- rhs.vars(x)
  var <- lhs.vars(x)

  # Validate inputs
  if (!is.null(var)) {
    stop("no dependend variable needed for this method, please remove it and try again")
  }

  if (nrow(object$data) == 0) {
    stop("An empty dataset detected (zero rows).")
  }

  class_nonprob <- class(object)[2]
  if (!class_nonprob %in% c("nonprobsvy_ipw", "nonprobsvy_dr", "nonprobsvy_mi")) {
    stop("Invalid nonprob object class")
  }

  # Input validation checks
  missing_vars <- setdiff(groups, names(object$data))
  if (length(missing_vars) > 0) {
    stop(sprintf(
      "The following variables are not present in the dataset: %s",
      paste(missing_vars, collapse = ", ")
    ))
  }

  if (class_nonprob %in% c("nonprobsvy_ipw", "nonprobsvy_dr")) {
    if (interaction) {
      # Get the data for grouping variables
      group_data <- object$data[, groups, drop = FALSE]

      # Check for NAs in grouping variables
      for (g in groups) {
        current_mask <- !is.na(group_data[[g]])
        if (sum(!current_mask) > 0) {
          warning(sprintf("NA values found in grouping variable %s", g))
        }
      }

      # Calculate weighted means for each combination of groups
      mean_value <- aggregate(object$weights,
                              by = group_data,
                              FUN = function(w) sum(w) / sum(object$weights)
      )

      # Check for small group sizes
      group_sizes <- aggregate(object$weights,
                               by = group_data,
                               FUN = length
      )
      small_groups <- group_sizes$x < 5
      if (any(small_groups)) {
        warning("Small group sizes (< 5) found in some combinations")
      }

      # Rename columns
      names(mean_value) <- c(groups, "mean")
    } else {
      data <- model.matrix(as.formula(paste(x, "- 1")), data = object$data)
      mean_value <- sapply(as.data.frame(data), function(col) weighted.mean(col, object$weights))
      mean_value <- data.frame(
        variable = names(mean_value),
        mean = unname(mean_value)
      )
    }
  } else {
    if (interaction) {
      form <- as.formula(paste("~interaction(", paste(groups, collapse = ", "), ")"))
      mean_value <- svymean(form, object$svydesign)
    } else {
      mean_value <- svymean(x, object$svydesign)
    }
    mean_value <- as.data.frame(mean_value)
  }
  mean_value
}
#' @title Mean values of covariates in subgroups
#'
#' @param x - formula
#' @param object - object of nonprobsvy class
#' @param interaction - logical, if TRUE calculate for all combinations of grouping variables
#'
#' @return A data frame with estimated means of the given covariates in subgroups
#' @importFrom formula.tools lhs.vars
#' @importFrom formula.tools rhs.vars
#' @importFrom stats aggregate
#' @importFrom survey svymean
nonprobsvymean <- function(x, object, interaction) {
  UseMethod("nonprobsvymean", object)
}
# Internal function - not exported in CRAN version
# Will be exported in future releases after variance estimation is implemented
#' @keywords internal
nonprobsvyby.nonprobsvy <- function(y, by, object, FUN) {
  # TODO DR estimator and variances - not ready for CRAN
  # Validate FUN parameter
  if (!FUN %in% c("total", "mean")) {
    stop("FUN must be either 'total' or 'mean'")
  }

  if (nrow(object$data) == 0) {
    stop("An empty dataset detected (zero rows).")
  }

  class_nonprob <- class(object)[2]
  if (!class_nonprob %in% c("nonprobsvy_ipw", "nonprobsvy_dr", "nonprobsvy_mi")) {
    stop("An invalid `nonprob` object class.")
  }

  variables <- rhs.vars(y)
  groups <- rhs.vars(by)

  # Validate inputs
  missing_vars <- setdiff(groups, names(object$data))
  if (length(missing_vars) > 0) {
    stop(sprintf(
      "The following variables are not present in the dataset: %s.",
      paste(missing_vars, collapse = ", ")
    ))
  }

  if (!is.null(variables) && !all(variables %in% names(object$data))) {
    missing_dep_vars <- setdiff(variables, names(object$y))
    stop(sprintf(
      "The following dependent variables are not present: %s.",
      paste(missing_dep_vars, collapse = ", ")
    ))
  }

  # Check for NAs in grouping variables
  for (g in groups) {
    current_mask <- !is.na(object$data[[g]])
    if (sum(!current_mask) > 0) {
      warning(sprintf("NA values found in grouping variable %s", g))
    }
  }

  # Common data preparation
  data <- object$data[, c(variables, groups)]
  weights <- object$weights

  if (FUN == "total") {
    if (class_nonprob == "nonprobsvy_ipw") {
      valid_data <- which(object$R == 1)
      data <- data[valid_data, ]
      weights <- weights[valid_data]

      res <- aggregate(data[, variables] * weights,
                       by = data[, groups, drop = FALSE],
                       sum
      )
      names(res)[ncol(res)] <- paste0("total.", variables)
    } else if (class_nonprob == "nonprobsvy_mi") {
      res <- svyby(formula = ~y_hat_MI, by = by, design = object$svydesign, svytotal)
    }
  } else { # mean
    if (class_nonprob == "nonprobsvy_ipw") {
      res <- aggregate(data[, variables, drop = FALSE],
                       by = data[, groups, drop = FALSE],
                       FUN = function(y, w) weighted.mean(y, w = w[1:length(y)]),
                       w = weights
      )
      names(res) <- c(groups, paste0("mean.", variables))

      # Check for small group sizes
      group_sizes <- aggregate(rep(1, nrow(data)),
                               by = data[, groups, drop = FALSE],
                               FUN = sum
      )
      small_groups <- group_sizes$x < 5
      if (any(small_groups)) {
        warning("Small group sizes (< 5) found in some combinations")
      }
    } else if (class_nonprob == "nonprobsvy_mi") {
      res <- svyby(formula = ~y_hat_MI, by = by, design = object$svydesign, svymean)
    }
  }
  res
}
#' @title Statistics by groups
#'
#' @param y - formula for variable of interest
#' @param by - formula for grouping variables
#' @param object - object of nonprobsvy class
#' @param FUN - string specifying the function to apply ("mean" or "total")
#'
#' @return A data frame with estimated statistics of the given covariates by groups
#' @importFrom formula.tools lhs.vars
#' @importFrom formula.tools rhs.vars
#' @importFrom stats aggregate
#' @importFrom survey svyby
nonprobsvyby <- function(y, by, object, FUN) {
  UseMethod("nonprobsvyby", object)
}
