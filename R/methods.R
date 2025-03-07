
# my methods --------------------------------------------------------------

#' @method pop_size nonprob
#' @exportS3Method
pop_size.nonprob <- function(object) {
  object$pop_size
}

#' @title Returns population size (estimated or fixed)
#' @description Returns population size that is assumed to be
#'\itemize{
#'  \item{\code{fixed} -- if it is based on the `pop_size` argument,}
#'  \item{\code{estimated} -- if it is based on the probability survey specified in the `svydesign` or based on the estimated propensity scores for the non-probability sample.}
#'}
#' @param object object returned by the `nonprob` function.
#' @return a scalar returning the value of the population size.
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
#' ## estimated population size based on the non-calibrated IPW (MLE)
#' pop_size(ipw_est1)
#'
#' ## estimated population size based on the calibrated IPW (GEE)
#' pop_size(ipw_est2)
#'
#'
#' @export
pop_size <- function(object) {
  UseMethod("pop_size")
}


# base R methods ----------------------------------------------------------

#' @method nobs nonprob
#' @importFrom stats nobs
#' @exportS3Method
nobs.nonprob <- function(object,
                         ...) {
  c("prob" = object$prob_size, "nonprob" = object$nonprob_size)
}

#' @title Extract IPW weights
#' @description A generic function `weights` that returns inverse probability weights (if present)
#'
#' @param object a `nonprob` class object
#' @param ... other arguments passed to methods (currently not supported)
#'
#' @returns A vector of weights or a `NULL` extracted from the `nonprob` object i.e. element `"ipw_weights"`
#'
#' @method weights nonprob
#' @importFrom stats weights
#' @exportS3Method
weights.nonprob <- function(object,
                            ...) {
  object$ipw_weights
}

#' @method residuals nonprob
#' @importFrom stats residuals
#' @exportS3Method
residuals.nonprob <- function(object,
                                 type = c(
                                   "pearson",
                                   "working",
                                   "deviance",
                                   "response",
                                   "pearsonSTD"
                                 ),
                                 ...) { # TODO for pop_totals (only for non-prob part)

  if (length(type) > 1) {
    type <- "response"
  }
  if (object$estimator %in% c("dr", "mi")) {
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
  if (object$estimator %in% c("dr", "ipw")) {
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
  if (object$estimator == "mi") res <- list(outcome = res_out)
  if (object$estimator == "ipw") res <- list(selection = res_sel)
  if (object$estimator == "dr") res <- list(selection = res_sel, outcome = res_out)
  res
}
#' @method cooks.distance nonprob
#' @importFrom stats cooks.distance
#' @exportS3Method
cooks.distance.nonprob <- function(model,
                                      ...) {
  resids <- residuals(model, type = "pearsonSTD")
  hats <- hatvalues(model)
  if (model$estimator == "mi") {
    residuals_out <- resids$outcome^2
    res_out <- cooks.distance(model$outcome[[1]])
    res <- list(outcome = res_out)
  } else if (model$estimator == "ipw") {
    residuals_sel <- resids$selection^2
    hats <- hats$selection
    res_sel <- (residuals_sel * (hats / (length(model$selection$coefficients))))
    res <- list(selection = res_sel)
  } else if (model$estimator == "dr") {
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
#' @method hatvalues nonprob
#' @importFrom stats hatvalues
#' @importFrom Matrix Diagonal
#' @exportS3Method
hatvalues.nonprob <- function(model,
                                 ...) { # TODO reduce execution time and glm.fit object and customise to variable selection
  if (model$estimator %in% c("dr", "ipw")) {
    X <- model$X
    propensity_scores <- model$prob
    W <- Matrix::Diagonal(x = propensity_scores * (1 - propensity_scores))
    # W <- as.vector(propensity_scores * (1 - propensity_scores))
    H <- as.matrix(sqrt(W) %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% sqrt(W))
    hat_values_sel <- diag(H)
  }
  if (model$estimator %in% c("dr", "mi")) {
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

  if (model$estimator == "mi") res <- list(outcome = hat_values_out)
  if (model$estimator == "ipw") res <- list(selection = hat_values_sel)
  if (model$estimator == "dr") res <- list(selection = hat_values_sel, outcome = hat_values_out)
  res
}
# CODE MODIFIED FROM stats:::logLik.glm
#' @method logLik nonprob
#' @importFrom stats logLik
#' @exportS3Method
logLik.nonprob <- function(object, ...) { # TODO consider adding appropriate warning/error in case of running on non mle method
  if (object$estimator %in% c("dr", "ipw")) {
    val_sel <- object$selection$log_likelihood
    attr(val_sel, "nobs") <- dim(residuals(object, type = "pearson"))[1]
    attr(val_sel, "df") <- length(object$selection$coefficients)
    class(val_sel) <- "logLik"
  }
  if (object$estimator %in% c("dr", "mi")) {
    if (class(object$outcome[[1]])[1] == "glm") {
      val_out <- logLik(object$outcome[[1]])
    } else {
      val_out <- "NULL"
    }
  }
  if (object$estimator == "mi") val <- c("outcome" = val_out)
  if (object$estimator == "ipw") val <- c("selection" = val_sel)
  if (object$estimator == "dr") val <- c("selection" = val_sel, "outcome" = val_out)
  val
}
#' @method AIC nonprob
#' @importFrom stats AIC
#' @exportS3Method
AIC.nonprob <- function(object,
                           ...) {
  if (object$estimator %in% c("dr", "ipw")) {
    if (!is.na(object$selection$log_likelihood)) {
      res_sel <- 2 * (length(object$selection$coefficients) - object$selection$log_likelihood)
    } else {
      res_sel <- NA
    }
  }
  if (object$estimator %in% c("dr", "mi")) res_out <- object$outcome[[1]]$aic
  if (object$estimator == "mi") res <- c("outcome" = res_out)
  if (object$estimator == "ipw") res <- c("selection" = res_sel)
  if (object$estimator == "dr") res <- c("selection" = res_sel, "outcome" = res_out)
  res
}
#' @method BIC nonprob
#' @importFrom stats BIC
#' @exportS3Method
BIC.nonprob <- function(object,
                           ...) {
  if (object$estimator %in% c("dr", "ipw")) {
    if (!is.na(object$selection$log_likelihood)) {
      res_sel <- length(object$selection$coefficients) * log(object$nonprob_size + object$prob_size) - 2 * object$selection$log_likelihood
    } else {
      res_sel <- NA
    }
  }
  if (object$estimator %in% c("dr", "mi")) {
    if (class(object$outcome[[1]])[1] == "glm") {
      res_out <- BIC(object$outcome[[1]])
    } else {
      res_out <- NA
    }
  }
  if (object$estimator == "mi") res <- c("outcome" = res_out)
  if (object$estimator == "ipw") res <- c("selection" = res_sel)
  if (object$estimator == "dr") res <- list("selection" = res_sel, "outcome" = res_out)
  res
}
#' @title Confidence Intervals for Model Parameters
#'
#' @description A function that computes confidence intervals
#' for selection model coefficients.
#'
#' @param object object of `nonprob` class.
#' @param parm names of parameters for which confidence intervals are to be
#' computed, if missing all parameters will be considered.
#' @param level confidence level for intervals.
#' @param ... additional arguments
#'
#' @method confint nonprob
#' @return returns a `data.frame` with confidence intervals for the target variables
#' @importFrom stats confint
#' @exportS3Method
confint.nonprob <- function(object,
                            parm,
                            level = 0.95,
                            ...) {

  if (all(is.na(object$output$SE))) stop("Standard errors were not calculated. Please refit the `nonprob` object.")

  if (missing(parm)) parm <- rownames(object$output)

  if (level == 0.95) {
    CIs <- object$confidence_interval
    CIs$target <- rownames(CIs)
    rownames(CIs) <- NULL
  } else {
    CIs <- object$output
    z <- stats::qnorm(1 - (1-level) / 2)
    # confidence interval based on the normal approximation
    CIs$lower_bound <- CIs$mean - z * CIs$SE
    CIs$upper_bound <- CIs$mean + z * CIs$SE
    CIs$target <- rownames(CIs)
    rownames(CIs) <- NULL
  }
  return(CIs[CIs$target %in% parm, c("target", "lower_bound", "upper_bound")])
}
#' @title Obtain Covariance Matrix estimation.
#'
#' @description A \code{vcov} method for \code{`nonprob`} class.
#'
#' @param object object of `nonprob` class.
#' @param ... additional arguments for method functions
#'
#' @details  Returns a estimated covariance matrix for model coefficients
#' calculated from analytic hessian or Fisher information matrix usually
#' utilising asymptotic effectiveness of maximum likelihood estimates.
#'
#' @method vcov nonprob
#' @return A covariance matrix for fitted coefficients
#' @importFrom stats vcov
#' @exportS3Method
vcov.nonprob <- function(object,
                            ...) { # TODO consider different vcov methods for selection and outcome models
  if (object$estimator %in% c("dr", "mi")) {
    if (class(object$outcome[[1]])[1] == "glm") {
      res_out <- vcov(object$outcome[[1]])
    } else {
      res_out <- object$outcome[[1]]$variance_covariance
    }
  }
  if (object$estimator %in% c("dr", "ipw")) res_sel <- object$selection$variance_covariance
  if (object$estimator == "mi") res <- list(outcome = res_out)
  if (object$estimator == "ipw") res <- list(selection = res_sel)
  if (object$estimator == "dr") res <- list(selection = res_sel, outcome = res_out)
  res
}
#' @method deviance nonprob
#' @importFrom stats deviance
#' @exportS3Method
deviance.nonprob <- function(object,
                                ...) {
  if (object$estimator %in% c("dr", "mi")) {
    if ((class(object$outcome[[1]])[1] == "glm")) {
      res_out <- deviance(object$outcome[[1]])
    } else {
      res_out <- "NULL"
    }
  }
  if (object$estimator %in% c("dr", "ipw")) {
    message("Deviance for selection model not available yet.")
    if (!is.character(object$selection$log_likelihood)) {
      res_sel <- NA
    } else {
      res_sel <- NA
    }
  }
  if (object$estimator == "mi") res <- c("outcome" = res_out)
  if (object$estimator == "ipw") res <- c("selection" = res_sel)
  if (object$estimator == "dr") res <- c("selection" = res_sel, "outcome" = res_out)
  res
}

