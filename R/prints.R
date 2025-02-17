# no print doccumentation

#' @method print nonprob
#' @exportS3Method
print.nonprob <- function(x, digits = 8, ...) {
  if (!is.null(x$call)) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n",
      sep = ""
    )
  }
  cat(
    "Estimated population mean with overall std.err and confidence interval:\n\n"
  )
  print(cbind(mean = x$output$mean, SE = x$output$SE, x$confidence_interval))
  invisible(x)
}
#' @method print summary_nonprob
#' @importFrom stats printCoefmat
#' @exportS3Method
print.summary_nonprob <- function(x,
                                     signif.stars = getOption("show.signif.stars"),
                                     digits = max(3L, getOption("digits") - 3L),
                                     ...) { # TODO regression diagnostics divided into outcome and selection models
  if (!is.null(x$call)) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n",
      sep = ""
    )
  }
  # TODO add printing the Info only for DR and MI models
  if (length(x$pop_total$mean) > 1) {
    cat("Info:\n", "The summary contains information mainly on the first outcome variable.\n",
      "More details on the estimation of this variable and others can be found in the outcome list of nonprob object.",
      "\n\n",
      sep = ""
    )
  }

  # cat("Residuals:\n")
  # print(summary(c(x$residuals[, 1])))

  cat("-------------------------\n")

  # cat(
  #   sep = "",
  #   "Estimated population mean: ", format(x$pop_total$mean, digits = digits),
  #   " with overall std.err of: ", format(x$pop_total$se[1], digits = digits),
  #   "\nAnd std.err for nonprobability and probability samples being respectively:\n",
  #   ifelse(!is.null(x$pop_total$cnf_int), format(x$pop_total$se[3], digits = digits), ""), " and ", ifelse(!is.null(x$pop_total$cnf_int), format(x$pop_total$se[2], digits = digits), ""), # TODO for se = FALSE
  #   "\n\nBased on: ", x$model, " method",
  #   "\n\n",(1 - x$control$control_inference$alpha)*100, "% Confidence inverval for popualtion mean:\n"
  # )
  # print(x$pop_total$cnf_int)

  cat(
    sep = "",
    "Estimated population mean: ", format(x$pop_total$mean[1], digits = digits)
  )
  if (!is.null(x$pop_total$cnf_int)) {
    cat(
      sep = "",
      " with overall std.err of: ", format(x$pop_total$se[1], digits = digits),
      "\nAnd std.err for nonprobability and probability samples being respectively:\n",
      format(x$pop_total$se[3], digits = digits), " and ", format(x$pop_total$se[2], digits = digits), # TODO for se = FALSE
      "\n\n", (1 - x$control$control_inference$alpha) * 100, "% Confidence inverval for popualtion mean:\n"
    )
    print(x$pop_total$cnf_int)
  }

  cat(sep = "", "\n\nBased on: ", x$model, " method")


  cat(
    sep = "",
    "\nFor a population of estimate size: ", x$population_size,
    "\nObtained on a nonprobability sample of size: ", x$sample_size[2],
    "\nWith an auxiliary probability sample of size: ", x$sample_size[1], # TODO
    "\n"
  )

  cat("-------------------------\n\n")

  k <- length(x$names)
  if (k > 0) cat("Regression coefficients:")
  while (k > 0) {
    cat("\n-----------------------\nFor ", x$names[k], ":\n", sep = "")
    printCoefmat(
      matrix( # TODO:: add conf intervals
        data = c(x$coef[[k]], x$std_err[[k]], x$w_val[[k]], x$p_values[[k]]), # TODO named coefs for selection model
        ncol = 4,
        dimnames = list(
          names(x$coef[[k]]),
          switch(x$test,
            "t" = c("Estimate", "Std. Error", "t value", "P(>|t|)"),
            "z" = c("Estimate", "Std. Error", "z value", "P(>|z|)")
          )
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

    # cat("\nAIC:")
    # print(x$aic)
    # cat("BIC:")
    # print(x$bic)
    cat("\nAIC: ", x$aic[[1]], "\nBIC: ", x$bic[[1]], sep = "")

    cat("\nLog-Likelihood:", x$likelihood, "on", x$df_residual, "Degrees of freedom\n")

    # cat("-------------------------\n\n")
  }


  # cat("\nRegression diagnostics:") #TODO

  invisible(x)
}

