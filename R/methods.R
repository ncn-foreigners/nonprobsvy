#' @title Summary statistics for model of nonprobsvy class.
#'
#' @param object object of nonprobsvy class
#' @param test Type of test for significance of parameters \code{"t"} for t-test
#' and \code{"z"} for normal approximation of students t distribution, by
#' default \code{"z"} is used if there are more than 30 degrees of freedom
#' and \code{"t"} is used in other cases.
#' @param regression_confint confint Logical value indicating whether confidence intervals for
#' regression parameters should be constructed
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
                               # regression_confint = FALSE,
                               # cov = NULL, # in case of adding sandwich methods
                               ...) {

  model_specific_info <- specific_summary_info(
    object,
    correlation = correlation,
    ...
  )
  df_residual <- model_specific_info$df_residuals
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
      #if(isTRUE(confint)) {confint.nonprobsvy(object, ...)} else {NULL})
    } else {
      warning("TODO")
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
      population_size = pop_size.nonprobsvy(object, ...), # change to pop_size
      test = test,
      control = object$control,
      model = switch(
        class(object)[2],
        "nonprobsvy_dr"  = "Doubly-Robust",
        "nonprobsvy_ipw" = "Inverse probability weighted",
        "nonprobsvy_mi"  = "Mass Imputation"
      ),
      aic = ifelse(class(object)[2] %in% c("nonprobsvy_dr", "nonprobsvy_ipw"), AIC(object), ""),
      bic = ifelse(class(object)[2] %in% c("nonprobsvy_dr", "nonprobsvy_ipw"), BIC(object), ""),
      residuals = residuals.nonprobsvy(object, type = "response"),
      likelihood = ifelse(class(object)[2] %in% c("nonprobsvy_dr", "nonprobsvy_ipw"), object$log_likelihood, ""),
      df_residual = ifelse(class(object)[2] %in% c("nonprobsvy_dr", "nonprobsvy_ipw"), object$df_residual, ""),
      weights = summary(object$weights),
      coef = cf,
      std_err = se,
      w_val = wald_test_stat,
      p_values = p_values,
      crr = crr,
      confidence_interval_coef = confidence_interval_coef,
      names = attr(model_specific_info, "TODO")
    ),
    class = c("summary_nonprobsvy")
  )

  res
}

# no need for documenting simple functions

#' @method print summary_nonprobsvy
#' @importFrom stats printCoefmat
#' @exportS3Method
print.summary_nonprobsvy <- function(x,
                                     signif.stars = getOption("show.signif.stars"),
                                     digits = max(3L, getOption("digits") - 3L),
                                     ...) {
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

  cat("Regression coefficients:")

  for (k in 1:length(x$names)) {
    cat("\n-----------------------\nFor", x$names[k], ":\n")
    printCoefmat(
      matrix(# TODO:: add conf intervals
        data = c(x$coef[[k]], x$std_err[[k]], x$w_val[[k]], x$p_values[[k]]),
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
  }

  cat("-------------------------\n\n")

  if (x$model %in% c("Doubly-Robust", "Inverse probability weighted")) {

    cat("Weights:\n")
    print(x$weights)

    cat("-------------------------\n\n")

    cat("Residuals:\n")
    print(summary(x$residuals))

    cat("\nAIC:")
    print(x$aic)
    cat("BIC:")
    print(x$bic)

    cat("Log-Likelihood:", x$likelihood, "on", x$df_residual, "Degrees of freedom\n")

    cat("-------------------------\n\n")

  }


  cat("\nRegression diagnostics:") #TODO

  invisible(x)
}

#' @method nobs nonprobsvy
#' @importFrom stats nobs
#' @exportS3Method
nobs.nonprobsvy <- function(object,
                            ...) {
  c("prob" = object$prob_size, "nonprob" = object$nonprob_size)
}
# Internal functions, no need for documenting them

#' @method pop_size nonprobsvy
#' @exportS3Method
pop_size.nonprobsvy <- function(object,
                                ...) {
  object$pop_size
}

#' @export
pop_size <- function(object, ...) {
  UseMethod("pop_size")
}

#' @method residuals nonprobsvy
#' @exportS3Method
residuals.nonprobsvy <- function(object,
                                 type = c("pearson",
                                          "deviance",
                                          "response")) { # TODO for pop_totals

  propensity_scores <- object$prop_scores
  if (!is.null(object$prob_size)) {
    R <- c(rep(1, object$nonprob_size), rep(0, object$prob_size))
    s <- c(rep(1, object$nonprob_size), rep(-1, object$prob_size))
  } else {
    R <- rep(1, object$nonprob_size)
    s <- rep(1, object$nonprob_size)
  }

  res <- switch(type,
                "pearson" = (R - propensity_scores)/sqrt(propensity_scores * (1 - propensity_scores)),
                "deviance" = s * sqrt(-2 * (R * log(propensity_scores) + (1 - R) * log(1 - propensity_scores))),
                "response" = R - propensity_scores) # TODO studentized_pearson, studentized_deviance
  res
}

#' @method pearson nonprobsvy
#' @exportS3Method
cooks.distance.nonprobsvy <- function(object,
                                      ...) { # TODO basing on Hat_matrix

}

#' @method pearson nonprobsvy
#' @importFrom stats hatvalues
#' @importFrom Matrix Diagonal
#' @exportS3Method
hatvalues.nonprobsvy <- function(object,
                                  ...) { #TODO reduce execution time
  propensity_scores <- object$prop_scores
  W <- Matrix::Diagonal(x = propensity_scores * (1 - propensity_scores))
  XWX_inv <-  solve(t(object$X) %*% W %*% object$X)
  hats <- vector(mode = "numeric", length = length(propensity_scores))
  for (i in 1:length(hats)) {
    hats[i] <- W[i,i] * object$X[i,] %*% XWX_inv %*% object$X[i,]
  }
  #hats <- Matrix::Diagonal(x = W %*% object$X %*% XWX_inv %*% t(object$X))
  hats
}

#' @method pearson nonprobsvy
#' @importFrom stats AIC
#' @exportS3Method
AIC.nonprobsvy <- function(object,
                           ...) {
  if (!is.character(object$log_likelihood)) {
    res <- 2 * (length(object$parameters) - object$log_likelihood)
  } else{
    res <- "AIC not available for this method"
  }
  res
}

#' @method pearson nonprobsvy
#' @importFrom stats BIC
#' @exportS3Method
BIC.nonprobsvy <- function(object,
                           ...) {
  if (!is.character(object$log_likelihood)) {
    res <- length(object$parameters) * log(object$nonprob_size + object$prob_size) - 2 * object$log_likelihood
  } else {
    res <- "BIC not available for this method"
  }
  res
}

#' @method confint nonprobsvy
#' @return An object with named columns that include upper and
#' lower limit of confidence intervals.
#' @exportS3Method
confint.nonprobsvy <- function(object,
                               level = 0.95,
                               ...) {
  res <- object$confidence_interval
  colnames(res) <- c(paste0(100 * (1 - level) / 2, "%"),
                     paste0(100 * (1 - (1 - level) / 2), "%"))
  res
}

#' @method vcov nonprobsvy
#' @return A covariance matrix for fitted coefficients
#' @exportS3Method
vcov.nonprobsvy <- function(object,
                            ...) {

}
