#' @title Summary statistics for model of nonprobsvy class.
#'
#' @param object TODO
#' @param test TODO
#' @param regression_confint TODO
#' @param correlation TODO
#' @param cov TODO
#' @param ... TODO
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
  df_residual <- 1000 # TODO:: for now just a big number
  if (missing(test)) {if (df_residual > 30) test <- "z" else test <- "t"}

  model_specific_info <- specific_summary_info(
    object,
    correlation = correlation,
    ...
  )

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
      stop("TODO")
    }
  }


  res <- structure(
    list(
      call = NULL,
      pop_total = list(
        mean = object$output$mean,
        se = c(object$output[, 2], object$SE$prob, object$SE$nonprob),
        cnf_int = object$confidence_interval
      ),
      #sample_size = nobs(object, ...),
      population_size = 0,# TODO:: pop_size(object, ...)
      test = test,
      model = switch(
        class(object)[2],
        "nonprobsvy_dr"  = "Doubly-Robust",
        "nonprobsvy_ipw" = "Invers probability weighted",
        "nonprobsvy_mi"  = "Mass Imputation"
      ),
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
  if (is.null(x$call)) {
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
    "\n\nBased on: ", x$model, " method",# TODO add control to output and add confidence level
    "\n\nConfidence inverval for popualtion mean:\n"
  )
  print(x$pop_total$cnf_int)

  cat(
    sep = "",
    "\nFor a population of estimate size: ", x$population_size,
    "\nObtained on a nonprobability sample of size: ", 0,
    "\nWith and auxiliary probability sample of size: ", 0, #TODO
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

  cat("\nRegression diagnostics:") #TODO

  invisible(x)
}

#' @method nobs nonprobsvy
#' @importFrom stats nobs
#' @exportS3Method
nobs.nonprobsvy <- function(object,
                            ...) {
  c("prob" = 1, "nonprob" = 2)
}
# Internal functions, no need for documenting them


pearson.nonprobsvy <- function(X_nons,
                               X_rand,
                               ps_nons,
                               est_ps_rand) {
  Rnons <- c(rep(1, nrow(X_nons)), rep(0, nrow(X_rand)))
  ps <- c(ps_nons, est_ps_rand)

  res <- (Rnons - ps)/sqrt(ps*(1-ps))
  res
}


deviance.nonprobsvy <- function(X_nons,
                                X_rand,
                                ps_nons,
                                est_ps_rand) {

  Rnons <- c(rep(1, nrow(X_nons)), rep(0, nrow(X_rand)))
  s <- c(rep(1, nrow(X_nons)), rep(-1, nrow(X_rand)))
  ps <- c(ps_nons, est_ps_rand)

  res <- s * sqrt(-2 * (Rnons * log(ps) + (1 - Rnons) * log(1 - ps)))
  res
}

# example using
#pearson_residuals <- pearson_nonprobsvy(X_nons, X_rand, , est_ps_rand) # pearson residuals for propensity score model
#deviance_residuals <- deviance_nonprobsvy(X_nons, X_rand, ps_nons, est_ps_rand) # deviance residuals for propensity score model

cooks_distance.nonprobsvy <- function() {

}


hatvalues.nonprobsvy <- function() {


}


AIC.nonprobsvy <- function() {


}

BIC.nonprobsvy <- function() {



}
