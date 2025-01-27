# no need for documenting simple functions

#' @method nobs nonprobsvy
#' @importFrom stats nobs
#' @exportS3Method
nobs.nonprobsvy <- function(object,
                            ...) {
  c("prob" = object$prob_size, "nonprob" = object$nonprob_size)
}
#' @method pop_size nonprobsvy
#' @exportS3Method
pop_size.nonprobsvy <- function(object,
                                ...) {
  object$pop_size
}
#' @title Estimate size of population
#' @description Estimate size of population
#' @param object object returned by `nonprobsvy`.
#' @param ... additional parameters
#' @return Vector returning the value of the estimated population size.
#' @export
pop_size <- function(object, ...) {
  UseMethod("pop_size")
}
#' @method residuals nonprobsvy
#' @importFrom stats residuals
#' @exportS3Method
residuals.nonprobsvy <- function(object,
                                 type = c(
                                   "pearson",
                                   "working",
                                   "deviance",
                                   "response",
                                   "pearsonSTD"
                                 ),
                                 ...) { # TODO for pop_totals

  if (length(type) > 1) {
    type <- "response"
  }
  if (any(c("nonprobsvy_dr", "nonprobsvy_mi") %in% class(object))) {
    if (type %in% c(
      "deviance", "pearson", "working",
      "response", "partial"
    )) {
      if (object$control$control_inference$bias_correction == TRUE) {
        r <- object$outcome[[1]]$residuals
        res_out <- switch(type,
          pearson = r,
          working = r * sqrt(object$selection$prior.weights),
          deviance = r * sqrt(abs(object$selection$prior.weights * object$outcome[[1]]$family$variance)),
          response = r / object$outcome[[1]]$family$mu
        )
      } else {
        res_out <- residuals(object$outcome[[1]], type = type)
      }
    } else if (type == "pearsonSTD") {
      r <- object$outcome[[1]]$residuals
      variance <- as.vector((t(r) %*% r) / length(object$outcome[[1]]$coefficients))
      res_out <- r / sqrt((1 - hatvalues(object)$outcome) * variance)
    }
  }
  if (any(c("nonprobsvy_dr", "nonprobsvy_ipw") %in% class(object))) {
    propensity_scores <- object$prob
    if (!is.null(object$prob_size)) {
      R <- c(rep(1, object$nonprob_size), rep(0, object$prob_size))
      s <- c(rep(1, object$nonprob_size), rep(-1, object$prob_size))
    } else {
      R <- rep(1, object$nonprob_size)
      s <- rep(1, object$nonprob_size)
    }
    r <- object$selection$residuals

    res_sel <- switch(type,
      "response" = r,
      "working" = r / propensity_scores * (1 - propensity_scores),
      "pearson" = r / sqrt(propensity_scores * (1 - propensity_scores)),
      "deviance" = s * sqrt(-2 * (R * log(propensity_scores) + (1 - R) * log(1 - propensity_scores))),
      "pearsonSTD" = r / sqrt((1 - hatvalues(object)$selection) * object$selection$variance)
    ) # TODO add partial
  }
  if (class(object)[2] == "nonprobsvy_mi") res <- list(outcome = res_out)
  if (class(object)[2] == "nonprobsvy_ipw") res <- list(selection = res_sel)
  if (class(object)[2] == "nonprobsvy_dr") res <- list(selection = res_sel, outcome = res_out)
  res
}
#' @method cooks.distance nonprobsvy
#' @importFrom stats cooks.distance
#' @exportS3Method
cooks.distance.nonprobsvy <- function(model,
                                      ...) {
  resids <- residuals(model, type = "pearsonSTD")
  hats <- hatvalues(model)
  if (class(model)[2] == "nonprobsvy_mi") {
    residuals_out <- resids$outcome^2
    res_out <- cooks.distance(model$outcome[[1]])
    res <- list(outcome = res_out)
  } else if (class(model)[2] == "nonprobsvy_ipw") {
    residuals_sel <- resids$selection^2
    hats <- hats$selection
    res_sel <- (residuals_sel * (hats / (length(model$selection$coefficients))))
    res <- list(selection = res_sel)
  } else if (class(model)[2] == "nonprobsvy_dr") {
    residuals_out <- resids$outcome^2
    if (model$control$control_inference$bias_correction == TRUE) {
      res_out <- (residuals_out * (hats$outcome / (length(model$outcome[[1]]$coefficients))))
    } else {
      res_out <- cooks.distance(model$outcome[[1]])
    }
    residuals_sel <- resids$selection^2
    hats <- hats$selection
    res_sel <- (residuals_sel * (hats / (length(model$selection$coefficients))))
    res <- list(selection = res_sel, outcome = res_out)
  }
  res
}
#' @method hatvalues nonprobsvy
#' @importFrom stats hatvalues
#' @importFrom Matrix Diagonal
#' @exportS3Method
hatvalues.nonprobsvy <- function(model,
                                 ...) { # TODO reduce execution time and glm.fit object and customise to variable selection
  if (any(c("nonprobsvy_dr", "nonprobsvy_ipw") %in% class(model))) {
    X <- model$X
    propensity_scores <- model$prob
    W <- Matrix::Diagonal(x = propensity_scores * (1 - propensity_scores))
    # W <- as.vector(propensity_scores * (1 - propensity_scores))
    H <- as.matrix(sqrt(W) %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% sqrt(W))
    hat_values_sel <- diag(H)
  }
  if (any(c("nonprobsvy_dr", "nonprobsvy_mi") %in% class(model))) {
    if (model$control$control_inference$bias_correction == TRUE) {
      # TODO for mm consider switch
      X_out <- model$outcome[[1]]$X
      if (model$outcome[[1]]$family$family == "gaussian") {
        hat_matrix <- X_out %*% solve(t(X_out) %*% X_out) %*% t(X_out)
      } else if (model$outcome[[1]]$family$family == "poisson") {
        W <- Matrix::Diagonal(x = 1 / model$outcome[[1]]$fitted_values)
        # W <- as.vector(1 / model$outcome[[1]]$fitted_values)
        hat_matrix <- X_out %*% solve(t(X_out) %*% W %*% X_out) %*% t(X_out)
      } else if (model$outcome[[1]]$family$family == "binomial") {
        W <- Matrix::Diagonal(x = model$outcome[[1]]$family$mu * (1 - model$outcome[[1]]$family$mu))
        # W <- as.vector(model$outcome[[1]]$family$mu * (1 - model$outcome[[1]]$family$mu))
        hat_matrix <- as.matrix(sqrt(W) %*% X_out %*% solve(t(X_out) %*% W %*% X_out) %*% t(X_out) %*% sqrt(W))
      }
      hat_values_out <- diag(hat_matrix)
    } else {
      hat_values_out <- hatvalues(model$outcome[[1]])
    }
  }

  if (class(model)[2] == "nonprobsvy_mi") res <- list(outcome = hat_values_out)
  if (class(model)[2] == "nonprobsvy_ipw") res <- list(selection = hat_values_sel)
  if (class(model)[2] == "nonprobsvy_dr") res <- list(selection = hat_values_sel, outcome = hat_values_out)
  res
}
# CODE MODIFIED FROM stats:::logLik.glm
#' @method logLik nonprobsvy
#' @importFrom stats logLik
#' @exportS3Method
logLik.nonprobsvy <- function(object, ...) { # TODO consider adding appropriate warning/error in case of running on non mle method
  if (any(c("nonprobsvy_dr", "nonprobsvy_ipw") %in% class(object))) {
    val_sel <- object$selection$log_likelihood
    attr(val_sel, "nobs") <- dim(residuals(object, type = "pearson"))[1]
    attr(val_sel, "df") <- length(object$selection$coefficients)
    class(val_sel) <- "logLik"
  }
  if (any(c("nonprobsvy_dr", "nonprobsvy_mi") %in% class(object))) {
    if (class(object$outcome[[1]])[1] == "glm") {
      val_out <- logLik(object$outcome[[1]])
    } else {
      val_out <- "NULL"
    }
  }
  if (class(object)[2] == "nonprobsvy_mi") val <- c("outcome" = val_out)
  if (class(object)[2] == "nonprobsvy_ipw") val <- c("selection" = val_sel)
  if (class(object)[2] == "nonprobsvy_dr") val <- c("selection" = val_sel, "outcome" = val_out)
  val
}
#' @method AIC nonprobsvy
#' @importFrom stats AIC
#' @exportS3Method
AIC.nonprobsvy <- function(object,
                           ...) {
  if (any(c("nonprobsvy_dr", "nonprobsvy_ipw") %in% class(object))) {
    if (!is.na(object$selection$log_likelihood)) {
      res_sel <- 2 * (length(object$selection$coefficients) - object$selection$log_likelihood)
    } else {
      res_sel <- NA
    }
  }
  if (any(c("nonprobsvy_dr", "nonprobsvy_mi") %in% class(object))) res_out <- object$outcome[[1]]$aic
  if (class(object)[2] == "nonprobsvy_mi") res <- c("outcome" = res_out)
  if (class(object)[2] == "nonprobsvy_ipw") res <- c("selection" = res_sel)
  if (class(object)[2] == "nonprobsvy_dr") res <- c("selection" = res_sel, "outcome" = res_out)
  res
}
#' @method BIC nonprobsvy
#' @importFrom stats BIC
#' @exportS3Method
BIC.nonprobsvy <- function(object,
                           ...) {
  if (any(c("nonprobsvy_dr", "nonprobsvy_ipw") %in% class(object))) {
    if (!is.na(object$selection$log_likelihood)) {
      res_sel <- length(object$selection$coefficients) * log(object$nonprob_size + object$prob_size) - 2 * object$selection$log_likelihood
    } else {
      res_sel <- NA
    }
  }
  if (any(c("nonprobsvy_dr", "nonprobsvy_mi") %in% class(object))) {
    if (class(object$outcome[[1]])[1] == "glm") {
      res_out <- BIC(object$outcome[[1]])
    } else {
      res_out <- NA
    }
  }
  if (class(object)[2] == "nonprobsvy_mi") res <- c("outcome" = res_out)
  if (class(object)[2] == "nonprobsvy_ipw") res <- c("selection" = res_sel)
  if (class(object)[2] == "nonprobsvy_dr") res <- list("selection" = res_sel, "outcome" = res_out)
  res
}
#' @title Confidence Intervals for Model Parameters
#'
#' @description A function that computes confidence intervals
#' for selection model coefficients.
#'
#' @param object object of nonprobsvy class.
#' @param parm names of parameters for which confidence intervals are to be
#' computed, if missing all parameters will be considered.
#' @param level confidence level for intervals.
#' @param ... additional arguments
#'
#' @method confint nonprobsvy
#' @return An object with named columns that include upper and
#' lower limit of confidence intervals.
#' @importFrom stats confint
#' @exportS3Method
confint.nonprobsvy <- function(object,
                               parm,
                               level = 0.95,
                               ...) {
  if (any(c("nonprobsvy_dr", "nonprobsvy_ipw") %in% class(object))) {
    std <- object$selection$std_err
    sc <- qnorm(p = 1 - (1 - level) / 2)
    res_sel <- data.frame(object$selection$coefficients - sc * std, object$selection$coefficients + sc * std)
    colnames(res_sel) <- c(
      paste0(100 * (1 - level) / 2, "%"),
      paste0(100 * (1 - (1 - level) / 2), "%")
    )
  }
  if (any(c("nonprobsvy_dr", "nonprobsvy_mi") %in% class(object))) {
    if (class(object$outcome[[1]])[1] == "glm") {
      res_out <- confint(object$outcome[[1]])
    } else {
      std <- sqrt(diag(vcov(object)$outcome))
      sc <- qnorm(p = 1 - (1 - level) / 2)
      res_out <- data.frame(object$outcome[[1]]$coefficients - sc * std, object$outcome[[1]]$coefficients + sc * std)
      colnames(res_out) <- c(
        paste0(100 * (1 - level) / 2, "%"),
        paste0(100 * (1 - (1 - level) / 2), "%")
      )
    }
  }
  if (class(object)[2] == "nonprobsvy_mi") res <- list(outcome = res_out)
  if (class(object)[2] == "nonprobsvy_ipw") res <- list(selection = res_sel)
  if (class(object)[2] == "nonprobsvy_dr") res <- list(selection = res_sel, outcome = res_out)
  res
}
#' @title Obtain Covariance Matrix estimation.
#'
#' @description A \code{vcov} method for \code{nonprobsvy} class.
#'
#' @param object object of nonprobsvy class.
#' @param ... additional arguments for method functions
#'
#' @details  Returns a estimated covariance matrix for model coefficients
#' calculated from analytic hessian or Fisher information matrix usually
#' utilising asymptotic effectiveness of maximum likelihood estimates.
#'
#' @method vcov nonprobsvy
#' @return A covariance matrix for fitted coefficients
#' @importFrom stats vcov
#' @exportS3Method
vcov.nonprobsvy <- function(object,
                            ...) { # TODO consider different vcov methods for selection and outcome models
  if (any(c("nonprobsvy_dr", "nonprobsvy_mi") %in% class(object))) {
    if (class(object$outcome[[1]])[1] == "glm") {
      res_out <- vcov(object$outcome[[1]])
    } else {
      res_out <- object$outcome[[1]]$variance_covariance
    }
  }
  if (any(c("nonprobsvy_dr", "nonprobsvy_ipw") %in% class(object))) res_sel <- object$selection$variance_covariance
  if (class(object)[2] == "nonprobsvy_mi") res <- list(outcome = res_out)
  if (class(object)[2] == "nonprobsvy_ipw") res <- list(selection = res_sel)
  if (class(object)[2] == "nonprobsvy_dr") res <- list(selection = res_sel, outcome = res_out)
  res
}
#' @method deviance nonprobsvy
#' @importFrom stats deviance
#' @exportS3Method
deviance.nonprobsvy <- function(object,
                                ...) {
  if (any(c("nonprobsvy_dr", "nonprobsvy_mi") %in% class(object))) {
    if ((class(object$outcome[[1]])[1] == "glm")) {
      res_out <- deviance(object$outcome[[1]])
    } else {
      res_out <- "NULL"
    }
  }
  if (any(c("nonprobsvy_dr", "nonprobsvy_ipw") %in% class(object))) {
    message("Deviance for selection model not available yet.")
    if (!is.character(object$selection$log_likelihood)) {
      res_sel <- NA
    } else {
      res_sel <- NA
    }
  }
  if (class(object)[2] == "nonprobsvy_mi") res <- c("outcome" = res_out)
  if (class(object)[2] == "nonprobsvy_ipw") res <- c("selection" = res_sel)
  if (class(object)[2] == "nonprobsvy_dr") res <- c("selection" = res_sel, "outcome" = res_out)
  res
}
#' @method check_balance nonprobsvy
#' @exportS3Method
check_balance.nonprobsvy <- function(x, object, dig = 10) {
  # Input validation
  if (!inherits(x, "formula")) {
    stop("The `x` argument must be a formula.")
  }

  if (missing(object) || is.null(object)) {
    stop("The `object` argument is required.")
  }

  if (!any(c("nonprobsvy_dr", "nonprobsvy_ipw") %in% class(object))) {
    stop("No estimated weights available. Only the IPW or the DR methods are supported.")
  }

  if (!is.numeric(dig) || dig < 0) {
    stop("The `dig` argument must be a non-negative number")
  }

  if (nrow(object$data) == 0) {
    stop("An empty dataset detected (zero rows).")
  }

  if (sum(object$weights) == 0) {
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
    stop(sprintf("The following variables are not present in the dataset: %s.",
                 paste(missing_vars, collapse = ", ")))
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
        sum(data_subset$weights)
      })
      names(totals) <- paste0(var, levels)
      return(totals)
    } else {
      # For numeric variables
      return(setNames(sum(data$weights * data[[var]], na.rm = TRUE), var))
    }
  }

  # Prepare data
  data <- object$data
  data$weights <- object$weights

  # Calculate nonprob totals
  nonprob_totals <- tryCatch({
    unlist(lapply(vars, function(var) calculate_totals(var, data)))
  }, error = function(e) {
    stop(sprintf("Error calculating nonprobability totals: %s.", e$message))
  })

  # Calculate probability totals
  if (!is.null(object$svydesign)) {
    # If svydesign exists
    prob_totals <- tryCatch({
      svy_total <- svytotal(x, object$svydesign)
      svy_totals <- as.vector(svy_total)
      names(svy_totals) <- names(svy_total)
      svy_totals
    }, error = function(e) {
      stop(sprintf("Error calculating survey totals: %s.", e$message))
    })
  } else {
    # Use population totals
    prob_totals <- object$selection$pop_totals
    if (is.null(prob_totals)) {
      stop("The `pop_totals` argument is null.")
    }
  }

  # Calculate and round differences
  diff <- tryCatch({
    round(nonprob_totals - prob_totals[names(nonprob_totals)], digits = dig)
  }, error = function(e) {
    stop(sprintf("Error calculating differences: %s.", e$message))
  })

  # Return results with meaningful names
  result <- list(
    nonprob_totals = nonprob_totals,
    prob_totals = prob_totals,
    balance = diff
  )
  result
}
#' @title Check the balance between probability and non-probability samples
#'
#' @param x Formula specifying variables to check
#' @param object Object of nonprobsvy class
#' @param dig Number of digits for rounding (default = 10)
#'
#' @return List containing nonprobability totals, probability totals, and their differences
#' @importFrom stats aggregate
#' @importFrom survey svytotal
#' @importFrom stats setNames
#' @export
check_balance <- function(x, object, dig) {
  UseMethod("check_balance", object)
}
# Internal function - not exported in CRAN version
# Will be exported in future releases after variance estimation is implemented
#' @keywords internal
nonprobsvytotal.nonprobsvy <- function(x, object, interaction = FALSE) {
  # TODO variance - NOT READY FOR CRAN
  groups <- rhs.vars(x)
  var <- lhs.vars(x)
  # Validate inputs
  if (!is.null(var)) {
    stop("No dependend variable needed for this method, please remove it and try again.")
  }
  if (nrow(object$data) == 0) {
    stop("An empty dataset detected (zero rows).")
  }
  class_nonprob <- class(object)[2]
  if (!class_nonprob %in% c("nonprobsvy_ipw", "nonprobsvy_dr", "nonprobsvy_mi")) {
    stop("Invalid nonprob object class")
  }
  # Check if all group variables exist in the dataset
  missing_vars <- setdiff(groups, names(object$data))
  if (length(missing_vars) > 0) {
    stop(sprintf("The following variables are not present in the dataset: %s",
                 paste(missing_vars, collapse = ", ")))
  }
  if (class_nonprob %in% c("nonprobsvy_ipw", "nonprobsvy_dr")) {
    if (interaction) {
      # For single categorical variable, use the same approach as non-interaction
      if (length(groups) == 1) {
        data <- model.matrix(as.formula(paste(x, "- 1")), data = object$data)
        result <- sapply(as.data.frame(data), function(col) sum(col * object$weights))
        result <- data.frame(
          variable = names(result),
          total = unname(result)
        )
        return(result)
      }

      data <- object$data[which(object$R == 1), groups, drop = FALSE]
      weights <- object$weights[which(object$R == 1)]
      # Check for NAs in grouping variables
      for (g in groups) {
        current_mask <- !is.na(data[[g]])
        if (sum(!current_mask) > 0) {
          warning(sprintf("NA values found in grouping variable %s", g))
        }
      }
      # Create all combinations of group variables
      group_combinations <- do.call(expand.grid, lapply(data, unique))
      # Calculate totals for each combination
      result <- data.frame(group_combinations)
      result$total <- NA
      for (i in 1:nrow(result)) {
        mask <- rep(TRUE, nrow(data))
        valid_combination <- TRUE
        for (g in groups) {
          current_val <- result[[g]][i]
          if (is.na(current_val)) {
            valid_combination <- FALSE
            break
          }
          current_mask <- !is.na(data[[g]]) & (data[[g]] == current_val)
          mask <- mask & current_mask
        }
        # Check if we have any valid observations for this combination
        if (valid_combination && sum(mask, na.rm = TRUE) > 0) {
          result$total[i] <- sum(weights[mask])
          # Warning for small group sizes
          if (sum(mask) < 5) {
            warning(sprintf("Small group size (%d) for combination: %s",
                            sum(mask),
                            paste(result[i, groups], collapse = ", ")))
          }
        }
      }
      # Remove rows where combination doesn't exist
      result <- result[!is.na(result$total), ]
      names(result) <- c(groups, "total")
    } else {
      data <- model.matrix(as.formula(paste(x, "- 1")), data = object$data)
      result <- sapply(as.data.frame(data), function(col) sum(col * object$weights))
      result <- data.frame(
        variable = names(result),
        total = unname(result)
      )
    }
  } else {
    if (interaction) {
      form <- as.formula(paste("~interaction(", paste(groups, collapse=", "), ")"))
      result <- svytotal(form, object$svydesign)
    } else {
      result <- svytotal(x, object$svydesign)
    }
    result <- as.data.frame(result)
  }
  result
}
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
    stop(sprintf("The following variables are not present in the dataset: %s",
                 paste(missing_vars, collapse = ", ")))
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
                              FUN = function(w) sum(w) / sum(object$weights))

      # Check for small group sizes
      group_sizes <- aggregate(object$weights,
                               by = group_data,
                               FUN = length)
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
      form <- as.formula(paste("~interaction(", paste(groups, collapse=", "), ")"))
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
    stop(sprintf("The following variables are not present in the dataset: %s.",
                 paste(missing_vars, collapse = ", ")))
  }

  if (!is.null(variables) && !all(variables %in% names(object$data))) {
    missing_dep_vars <- setdiff(variables, names(object$y))
    stop(sprintf("The following dependent variables are not present: %s.",
                 paste(missing_dep_vars, collapse = ", ")))
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
                       sum)
      names(res)[ncol(res)] <- paste0("total.", variables)

    } else if (class_nonprob == "nonprobsvy_mi") {
      res <- svyby(formula = ~ y_hat_MI, by = by, design = object$svydesign, svytotal)
    }
  } else { # mean
    if (class_nonprob == "nonprobsvy_ipw") {
      res <- aggregate(data[, variables, drop = FALSE],
                       by = data[, groups, drop = FALSE],
                       FUN = function(y, w) weighted.mean(y, w = w[1:length(y)]),
                       w = weights)
      names(res) <- c(groups, paste0("mean.", variables))

      # Check for small group sizes
      group_sizes <- aggregate(rep(1, nrow(data)),
                               by = data[, groups, drop = FALSE],
                               FUN = sum)
      small_groups <- group_sizes$x < 5
      if (any(small_groups)) {
        warning("Small group sizes (< 5) found in some combinations")
      }

    } else if (class_nonprob == "nonprobsvy_mi") {
      res <- svyby(formula = ~ y_hat_MI, by = by, design = object$svydesign, svymean)
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
