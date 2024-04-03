# no need for documenting simple functions

#' @method nobs nonprobsvy
#' @importFrom stats nobs
#' @exportS3Method
nobs.nonprobsvy <- function(object,
                            ...) {
  c("prob" = object$prob_size, "nonprob" = object$nonprob_size)
}
#' @method pop.size nonprobsvy
#' @exportS3Method
pop.size.nonprobsvy <- function(object,
                                ...) {
  object$pop_size
}
#' @title Estimate size of population
#' @description Estimate size of population
#' @param object object returned by `nonprobsvy`.
#' @param ... additional parameters
#' @return Vector returning the value of the estimated population size.
#' @export
pop.size <- function(object, ...) {
  UseMethod("pop.size")
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
