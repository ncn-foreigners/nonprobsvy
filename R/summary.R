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
#' @return An object of \code{summary_nonprobsvy} class containing:
#' \itemize{
#' \item \code{call} -- A call which created \code{object}.
#' \item \code{pop_total} -- A list containing information about the estimated population mean, its standard error and confidence interval.
#' \item \code{sample_size} -- The size of the samples used in the model.
#' \item \code{population_size} -- The estimated size of the population from which the nonoprobability sample was drawn.
#' \item \code{test} -- Type of statistical test performed.
#' \item \code{control} -- A List of control parameters used in fitting the model.
#' \item \code{model} -- A descriptive name of the model used, e.g., "Doubly-Robust", "Inverse probability weighted", or "Mass Imputation".
#' \item \code{aic} -- Akaike's information criterion.
#' \item \code{bic} -- Bayesian (Schwarz's) information criterion.
#' \item \code{residuals} -- Residuals from the model, providing information on the difference between observed and predicted values.
#' \item \code{likelihood} -- Logarithm of likelihood function evaluated at coefficients.
#' \item \code{df_residual} -- Residual degrees of freedom.
#' \item \code{weights} -- Distribution of estimated weights obtained from the model.
#' \item \code{coef} -- Regression coefficients estimated by the model.
#' \item \code{std_err} -- Standard errors of the regression coefficients.
#' \item \code{w_val} -- Wald statistic values for the significance testing of coefficients.
#' \item \code{p_values} -- P-values corresponding to the Wald statistic values, assessing the significance of coefficients.
#' \item \code{crr} -- The correlation matrix of the model coefficients, if requested.
#' \item \code{confidence_interval_coef} -- Confidence intervals for the model coefficients.
#' \item \code{names} -- Names of the fitted models.
#' }
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
    if (missing(test)) {
      if (df_residual > 30) test <- "z" else test <- "t"
    }
  } else {
    test <- "z" # TODO, for now just z-test in case of mi estimation
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

      p_values[[number]] <- switch(test,
        "t" = 2 * stats::pt(q = -abs(k[, 1] / k[, 2]), df = df_residual),
        "z" = 2 * stats::pnorm(q = abs(k[, 1] / k[, 2]), lower.tail = FALSE)
      )

      temp_correlation <- if (isFALSE(correlation)) {
        NULL
      } else {
        cov / outer(k[, 2], k[, 2])
      }
      if (isTRUE(correlation)) {
        rownames(temp_correlation) <- colnames(temp_correlation) <- names(rownames(k))
      }

      crr[[number]] <- temp_correlation

      # confidence_interval_coef <- append(confidence_interval_coef,
      # if(isTRUE(confint)) {confint(object, ...)} else {NULL})
    } else {
      # TODO
    }
  }
  if (!is.null(object$SE)) {
    se_mean <- c(object$output[, 2], object$SE$prob, object$SE$nonprob)
  } else {
    se_mean <- NULL
  }
  res <- structure(
    list(
      call = object$call,
      pop_total = list(
        mean = object$output$mean,
        se = se_mean,
        cnf_int = object$confidence_interval
      ),
      sample_size = nobs(object, ...),
      population_size = pop.size(object, ...),
      test = test,
      control = object$control,
      model = switch(class(object)[2],
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


# summary helper functions
# for now just a rough sketch
specific_summary_info <- function(object, ...) {
  UseMethod("specific_summary_info")
}

specific_summary_info.nonprobsvy_ipw <- function(object,
                                                 ...) {
  coeffs_sel <- matrix(c(object$selection$coefficients, object$selection$std_err),
    ncol = 2,
    dimnames = list(
      names(object$selection$coefficients),
      c("Estimate", "Std. Error")
    )
  )
  res <- list(
    coeffs_sel = coeffs_sel,
    weights = object$weights,
    df_residual = object$selection$df_residual
  )

  attr(res$coeffs_sel, "glm") <- TRUE
  attr(res$weights, "glm") <- FALSE
  attr(res$df_residual, "glm") <- FALSE # TODO
  attr(res, "model") <- c("glm regression on selection variable")
  res
}

specific_summary_info.nonprobsvy_mi <- function(object,
                                                ...) {
  if (object$outcome[[1]]$method == "glm") { # TODO for pmm
    coeffs_out <- matrix(c(object$outcome[[1]]$coefficients, object$outcome[[1]]$std_err),
      ncol = 2,
      dimnames = list(
        names(object$outcome[[1]]$coefficients),
        c("Estimate", "Std. Error")
      )
    )
  } else {
    coeffs_out <- "no coefficients"
  }

  res <- list(
    coeffs_out = coeffs_out
  )
  if (object$outcome[[1]]$method == "glm") {
    attr(res$coeffs_out, "glm") <- TRUE
    attr(res, "model") <- "glm regression on outcome variable"
  } else if (object$outcome[[1]]$method == "nn") {
    attr(res$coeffs_out, "glm") <- FALSE
  } else if (object$outcome[[1]]$method == "pmm") { # TODO
    attr(res$coeffs_out, "glm") <- FALSE
    # attr(res, "model") <- "glm regression on outcome variable"
  }
  res
}

specific_summary_info.nonprobsvy_dr <- function(object,
                                                ...) {
  coeffs_sel <- matrix(c(object$selection$coefficients, object$selection$std_err),
    ncol = 2,
    dimnames = list(
      names(object$selection$coefficients),
      c("Estimate", "Std. Error")
    )
  )


  if (object$outcome[[1]]$method == "glm") {
    coeffs_out <- matrix(c(object$outcome[[1]]$coefficients, object$outcome[[1]]$std_err),
      ncol = 2,
      dimnames = list(
        names(object$outcome[[1]]$coefficients),
        c("Estimate", "Std. Error")
      )
    )
  } else {
    coeffs_out <- "no coefficients"
  }

  res <- list(
    coeffs_sel = coeffs_sel,
    coeffs_out = coeffs_out,
    weights = object$weights,
    df_residual = object$selection$df_residual
  )
  attr(res$coeffs_sel, "glm") <- TRUE
  if (object$outcome[[1]]$method == "glm") {
    attr(res$coeffs_out, "glm") <- TRUE
    attr(res, "model") <- c(
      "glm regression on selection variable",
      "glm regression on outcome variable"
    )
  } else if (object$outcome[[1]]$method == "nn") {
    attr(res$coeffs_out, "glm") <- FALSE
    attr(res, "model") <- c("glm regression on selection variable")
  }
  attr(res$weights, "glm") <- FALSE
  attr(res$df_residual, "glm") <- FALSE

  res
}
