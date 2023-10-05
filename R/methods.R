#' @title Summary statistics for model of nonprobsvy class.
#'
#' @param object object of nonprobsvy class
#' @param test Type of test for significance of parameters \code{"t"} for t-test
#' and \code{"z"} for normal approximation of students t distribution, by
#' default \code{"z"} is used if there are more than 30 degrees of freedom
#' and \code{"t"} is used in other cases.
#' @param correlation correlation Logical value indicating whether correlation matrix should
#' be computed from covariance matrix by default \code{FALSE}.
#' @param cov Covariance matrix corresponding to regression parameters
#' @param ... Additional optional arguments
#'
#'
#' @method summary nonprobsvy
#' @importFrom stats pt
#' @importFrom stats coef
#' @importFrom stats sd
#' @exportS3Method
summary.nonprobsvy <- function(object,
                               test = c("t", "z"),
                               correlation = FALSE,
                               # regression_confint = FALSE, confint Logical value indicating whether confidence intervals for
                               #                             regression parameters should be constructed TODO
                               cov = NULL, # in case of adding sandwich methods
                               ...) {

  model_specific_info <- specific_summary_info(
    object,
    correlation = correlation,
    ...
  )
  df_residual <- model_specific_info$df_residual
  if (!is.null(df_residual)) {
    if (missing(test)) {if (df_residual > 30) test <- "z" else test <- "t"}
  } else {
    test <- "z"  # TODO, for now just z-test in case of mi estimation
  }

  cf <- list()
  se <- list()
  wald_test_stat <- list()
  p_values <- list()
  crr <- list()
  confidence_interval_coef <- list()

  for (k in model_specific_info) {
    if (attr(k, "glm")) {
      number <- length(se) + 1
      cf[[number]] <- k[, 1]
      se[[number]] <- k[, 2]
      wald_test_stat[[number]] <- k[, 1] / k[, 2]

      p_values[[number]] <- switch (test,
        "t" = 2 *    stats::pt(q = -abs(k[, 1] / k[, 2]), df = df_residual),
        "z" = 2 * stats::pnorm(q = abs(k[, 1] / k[, 2]), lower.tail = FALSE)
      )

      temp_correlation <- if (isFALSE(correlation)) {NULL} else {cov / outer(k[, 2], k[, 2])}
      if(isTRUE(correlation)) {rownames(temp_correlation) <- colnames(temp_correlation) <- names(rownames(k))}

      crr[[number]] <- temp_correlation

      #confidence_interval_coef <- append(confidence_interval_coef,
      #if(isTRUE(confint)) {confint(object, ...)} else {NULL})
    } else {
      # TODO
    }
  }
  res <- structure(
    list(
      call = object$call,
      pop_total = list(
        mean = object$output$mean,
        se = c(object$output[, 2], object$SE$prob, object$SE$nonprob),
        cnf_int = object$confidence_interval
      ),
      sample_size = nobs(object, ...),
      population_size = pop.size(object, ...),
      test = test,
      control = object$control,
      model = switch(
        class(object)[2],
        "nonprobsvy_dr"  = "Doubly-Robust",
        "nonprobsvy_ipw" = "Inverse probability weighted",
        "nonprobsvy_mi"  = "Mass Imputation"
      ),
      aic = ifelse(class(object)[2] %in% c("nonprobsvy_dr", "nonprobsvy_ipw"), AIC(object), "no value for the selected method"),
      bic = ifelse(class(object)[2] %in% c("nonprobsvy_dr", "nonprobsvy_ipw"), BIC(object), "no value for the selected method"),
      residuals = residuals.nonprobsvy(object, type = "response"),
      likelihood = ifelse(class(object)[2] %in% c("nonprobsvy_dr", "nonprobsvy_ipw"), object$selection$log_likelihood, "no value for the selected method"),
      df_residual = ifelse(class(object)[2] %in% c("nonprobsvy_dr", "nonprobsvy_ipw"), object$selection$df_residual, "no value for the selected method"),
      weights = summary(object$weights),
      coef = cf,
      std_err = se,
      w_val = wald_test_stat,
      p_values = p_values,
      crr = crr,
      confidence_interval_coef = confidence_interval_coef,
      names = attr(model_specific_info, "model")
    ),
    class = c("summary_nonprobsvy")
  )
  res
}

# no need for documenting simple functions

#' @method print nonprobsvy
#' @exportS3Method
print.nonprobsvy <- function(x, digits = 8, ...) {
  if (!is.null(x$call)) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
  }

  cat(
    "Estimated population mean: ", format(x$output$mean, digits = digits),
    "\nWith overall std.err of: ", format(x$output$SE, digits = digits),
    "\nConfidence interval:\n", sep = " "
  )
  print(x$confidence_interval)
  invisible(x)
}
#' @method print summary_nonprobsvy
#' @importFrom stats printCoefmat
#' @exportS3Method
print.summary_nonprobsvy <- function(x,
                                     signif.stars = getOption("show.signif.stars"),
                                     digits = max(3L, getOption("digits") - 3L),
                                     ...) { # TODO regression diagnostics divided into outcome and selection models
  if (!is.null(x$call)) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
  }

  #cat("Residuals:\n")
  #print(summary(c(x$residuals[, 1])))

  cat("-------------------------\n")

  cat(
    sep = "",
    "Estimated population mean: ", format(x$pop_total$mean, digits = digits),
    " with overall std.err of: ", format(x$pop_total$se[1], digits = digits),
    "\nAnd std.err for nonprobability and probability samples being respectively:\n",
    format(x$pop_total$se[3], digits = digits), " and ", format(x$pop_total$se[2], digits = digits),
    "\n\nBased on: ", x$model, " method",
    "\n\n",(1 - x$control$control_inference$alpha)*100, "% Confidence inverval for popualtion mean:\n"
  )
  print(x$pop_total$cnf_int)

  cat(
    sep = "",
    "\nFor a population of estimate size: ", x$population_size,
    "\nObtained on a nonprobability sample of size: ", x$sample_size[2],
    "\nWith an auxiliary probability sample of size: ", x$sample_size[1], #TODO
    "\n"
  )

  cat("-------------------------\n\n")

    k <- length(x$names)
    if (k > 0)  cat("Regression coefficients:")
    while(k > 0) {
      cat("\n-----------------------\nFor ", x$names[k], ":\n", sep = "")
      printCoefmat(
        matrix(# TODO:: add conf intervals
          data = c(x$coef[[k]], x$std_err[[k]], x$w_val[[k]], x$p_values[[k]]), #TODO named coefs for selection model
          ncol = 4,
          dimnames = list(
            names(x$coef[[k]]),
            switch(x$test,
            "t" = c("Estimate", "Std. Error", "t value", "P(>|t|)"),
            "z" = c("Estimate", "Std. Error", "z value", "P(>|z|)"))
          )
        ),
        digits = digits,
        signif.stars = signif.stars,
        signif.legend = if (k == length(x$names)) signif.stars else FALSE,
        P.values = TRUE, has.Pvalue = TRUE,
        na.print = "NA",
        ...
      )
      k <- k - 1
    }

    if (length(x$names) > 0) cat("-------------------------\n\n")

  if (x$model %in% c("Doubly-Robust", "Inverse probability weighted")) {

    cat("Weights:\n")
    print(x$weights)

    cat("-------------------------\n\n")

    cat("Residuals:\n")
    print(summary(x$residuals$selection))

    #cat("\nAIC:")
    #print(x$aic)
    #cat("BIC:")
    #print(x$bic)
    cat("\nAIC: ", x$aic[[1]], "\nBIC: ",x$bic[[1]], sep = "")

    cat("\nLog-Likelihood:", x$likelihood, "on", x$df_residual, "Degrees of freedom\n")

    #cat("-------------------------\n\n")

  }


  #cat("\nRegression diagnostics:") #TODO

  invisible(x)
}

# Internal functions, no need for documenting them

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
#' @description - Estimate size of population
#' @param object - object returned by `nonprobsvy`.
#' @param ... additional parameters
#' @export
pop.size <- function(object, ...) {
  UseMethod("pop.size")
}
#' @method residuals nonprobsvy
#' @importFrom stats residuals
#' @exportS3Method
residuals.nonprobsvy <- function(object,
                                 type = c("pearson",
                                          "deviance",
                                          "response"),
                                 ...) { # TODO for pop_totals and variable selection

  if (any(c("nonprobsvy_dr", "nonprobsvy_mi") %in% class(object))) {
    if (object$control$control_inference$vars_selection == FALSE) {
      res_out <- residuals(object$outcome) # TODO for variable selection
    } else { # TODO for variable selection
      r <- object$outcome$family$residuals
      res_out <- switch(type,
                        "response" = r,
                        "working" = r/object$outcome$family$mu,
                        # TODO "deviance" =
                        "pearson" = r/sqrt(object$outcome$family$variance)
                        )
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

    res_sel <- switch(type,
                      "pearson" = (R - propensity_scores)/sqrt(propensity_scores * (1 - propensity_scores)),
                      "deviance" = s * sqrt(-2 * (R * log(propensity_scores) + (1 - R) * log(1 - propensity_scores))),
                      "response" = R - propensity_scores) # TODO studentized_pearson, studentized_deviance
  }
  if (class(object)[2] == "nonprobsvy_mi") res <- list(outcome = res_out)
  if (class(object)[2] == "nonprobsvy_ipw") res <- list(selection = res_sel)
  if (class(object)[2] == "nonprobsvy_dr") res <- list(selection = res_sel, outcome = res_out)
  res
}
#' @method cooks.distance nonprobsvy
#' @importFrom stats cooks.distance
#' @exportS3Method
cooks.distance.nonprobsvy <- function(model, # TODO for variable selection
                                      ...) { # TODO basing on Hat_matrix,
  resids <- residuals(model, type = "pearsonSTD")^2
  hats <- hatvalues(model)
  res_sel <- (resids * (hats / (length(coef(model))))) # TODO
  res_out <- cooks.distance(model$outcome$glm)

}
#' @method hatvalues nonprobsvy
#' @importFrom stats hatvalues
#' @importFrom Matrix Diagonal
#' @exportS3Method
hatvalues.nonprobsvy <- function(model,
                                 ...) { # TODO reduce execution time and glm.fit object and customise to variable selection
  if (any(c("nonprobsvy_dr", "nonprobsvy_ipw") %in% class(model))) {
    propensity_scores <- model$prop_scores
    W <- Matrix::Diagonal(x = propensity_scores * (1 - propensity_scores))
    XWX_inv <-  solve(t(model$X) %*% W %*% model$X)
    hat_values_sel <- vector(mode = "numeric", length = length(propensity_scores))
    for (i in 1:length(hat_values_sel)) {
      hat_values_sel[i] <- W[i,i] * model$X[i,] %*% XWX_inv %*% model$X[i,]
    }
    #hats <- Matrix::Diagonal(x = W %*% object$X %*% XWX_inv %*% t(object$X))
    #names(hat_values) <- row.names(model$parameters)
  }
  if (any(c("nonprobsvy_dr", "nonprobsvy_mi") %in% class(model))) {
  hat_values_out <- hatvalues(model$outcome$glm) # TODO
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
logLik.nonprobsvy <- function(object, ...) {
  if (any(c("nonprobsvy_dr", "nonprobsvy_ipw") %in% class(object))) {
    val_sel <- object$selection$log_likelihood
    attr(val_sel, "nobs") <- dim(residuals(object, type = "pearson"))[1]
    attr(val_sel, "df") <- nrow(object$parameters) #length(object$coefficients)
    class(val_sel) <- "logLik"
  }
  val_out <- ifelse(any(c("nonprobsvy_dr", "nonprobsvy_mi") %in% class(object)), logLik(object$outcome), 0) # TODO for gee
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
  if (!is.character(object$selection$log_likelihood)) {
    if (any(c("nonprobsvy_dr", "nonprobsvy_ipw") %in% class(object))) res_sel <- 2 * (length(object$parameters) - object$selection$log_likelihood)
    if (any(c("nonprobsvy_dr", "nonprobsvy_mi") %in% class(object))) res_out <- object$outcome$aic
    if (class(object)[2] == "nonprobsvy_mi") res <- c("outcome" = res_out)
    if (class(object)[2] == "nonprobsvy_ipw") res <- c("selection" = res_sel)
    if (class(object)[2] == "nonprobsvy_dr") res <- c("selection" = res_sel, "outcome" = res_out)
  } else{
    res <- "not available for this method"
  }
  res
}
#' @method BIC nonprobsvy
#' @importFrom stats BIC
#' @importFrom HelpersMG ExtractAIC.glm
#' @exportS3Method
BIC.nonprobsvy <- function(object,
                           ...) {
  if (!is.character(object$selection$log_likelihood)) {
    if (any(c("nonprobsvy_dr", "nonprobsvy_ipw") %in% class(object))) res_sel <- length(object$parameters) * log(object$nonprob_size + object$prob_size) - 2 * object$selection$log_likelihood
    if (any(c("nonprobsvy_dr", "nonprobsvy_mi") %in% class(object))) {

      if (!is.null(object$outcome$coefficients)) {
        if (object$control$control_inference$vars_selection == TRUE) {
          options(AIC="BIC")
          res_out <- HelpersMG::ExtractAIC.glm(object$outcome)[2]
        } else {
          res_out <- BIC(object$outcome)
        }
      } else {
        res_out <- "not available for this methos"
      }
    }
    if (class(object)[2] == "nonprobsvy_mi") res <- c("outcome" = res_out)
    if (class(object)[2] == "nonprobsvy_ipw") res <- c("selection" = res_sel)
    if (class(object)[2] == "nonprobsvy_dr") res <- list("selection" = res_sel, "outcome" = res_out)
  } else {
    res <- "not available for this method"
  }
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
    std <- sqrt(diag(vcov(object)[[1]]))
    sc <- qnorm(p = 1 - (1 - level) / 2)
    res_sel <- data.frame(object$parameters[,1] - sc * std, object$parameters[,1] + sc * std)
    colnames(res_sel) <- c(paste0(100 * (1 - level) / 2, "%"),
                           paste0(100 * (1 - (1 - level) / 2), "%"))
  }
  if (any(c("nonprobsvy_dr", "nonprobsvy_mi") %in% class(object))) {
    if (object$control$control_inference$vars_selection == FALSE) {
    res_out <- confint(object$outcome)
    } else {
      std <- sqrt(diag(vcov(object)[[1]]))
      sc <- qnorm(p = 1 - (1 - level) / 2)
      res_out <- data.frame(object$beta[,1] - sc * std, object$beta[,1] + sc * std)
      colnames(res_out) <- c(paste0(100 * (1 - level) / 2, "%"),
                             paste0(100 * (1 - (1 - level) / 2), "%"))
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
  if (object$control$control_inference$vars_selection == TRUE) {
    if (any(c("nonprobsvy_dr", "nonprobsvy_mi") %in% class(object))) res_out <- object$outcome$variance_covariance
    if (any(c("nonprobsvy_dr", "nonprobsvy_ipw") %in% class(object)))  res_sel <- object$selection$variance_covariance
    if (class(object)[2] == "nonprobsvy_mi") res <- list(outcome = res_out)
    if (class(object)[2] == "nonprobsvy_ipw") res <- list(selection = res_sel)
    if (class(object)[2] == "nonprobsvy_dr") res <- list(selection = res_sel, outcome = res_out)

  } else {
    if (any(c("nonprobsvy_dr", "nonprobsvy_mi") %in% class(object))) res_out <- vcov(object$outcome)
    if (any(c("nonprobsvy_dr", "nonprobsvy_ipw") %in% class(object)))  res_sel <- object$selection$variance_covariance
    if (class(object)[2] == "nonprobsvy_mi") res <- list(outcome = res_out)
    if (class(object)[2] == "nonprobsvy_ipw") res <- list(selection = res_sel)
    if (class(object)[2] == "nonprobsvy_dr") res <- list(selection = res_sel, outcome = res_out)
    }
  res
}
#' @method deviance nonprobsvy
#' @importFrom stats deviance
#' @exportS3Method
deviance.nonprobsvy <- function(object,
                                ...) {
  res_out <- object$outcome$deviance
  # TODO for selection model - use selection object
}
