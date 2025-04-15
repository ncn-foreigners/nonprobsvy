#' @method anova nonprob
#' @exportS3Method
anova.nonprob <- function(object, ...) {
  stop("The `anova` method is not implemented for the `nonprob` class. If you would like to compare models (sets of variables), compare point and interval estimates. If you would like to assess variables, do it before applying the `nonprob` function.")
}

