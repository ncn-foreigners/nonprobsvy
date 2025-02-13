#' @method anova nonprobsvy
#' @exportS3Method
anova.nonprobsvy <- function(object, ...) {
  stop("The `anova` method is not implemented for the `nonprobsvy` class. If you would like to compare models (sets of variables), compare point and interval estimates. If you would like to assess variables, do it before applying the `nonprob` function.")
}

#' @method plot nonprobsvy
#' @exportS3Method
plot.nonprobsvy <- function(x, ...) {
  stop("We do not provide tools for visual assessment of the results. If you are interested in covariate balance plots, we recommend using the `cobalt` package. If you are interested in evaluation of the models (e.g. IPW, MI), we recommend using base R functions or the `modelsummary` package.")
}
