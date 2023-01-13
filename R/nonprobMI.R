#' nonprobMI
#
#' nonprobMI: Function for inference based on nonprobability surveys using mass imputation
#
#' @param selection - `formula`, the selection (propensity) equation.
#' @param outcome - `formula`, the outcome equation.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param family.outcome - a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param subset - an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param weights - an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na.action
#' @param control.outcome
#' @param control.inference
#' @param start
#' @param verbose
#' @param contrasts
#' @param model
#' @param x
#' @param y
#' @param ...
#' @importFrom stats glm.fit
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @export

nonprobMI <- function(outcome,
                      data,
                      svydesign,
                      method.outcome,
                      family.outcome,
                      subset,
                      weights,
                      na.action,
                      control.outcome = controlOut(),
                      control.inference = controlInf(),
                      start,
                      verbose,
                      contrasts,
                      model,
                      x,
                      y,
                      ...) {

  XY_nons <- model.frame(outcome, data)
  X_nons <- model.matrix(XY_nons, data)
  X_rand <- model.matrix(outcome, svydesign$variables)

  ## estimation
  model_nons <- nonprobMI.fit()

  model_nons_coefs <- as.matrix(model_nons$coefficients)

  y_rand_pred <-  as.numeric(X_rand %*% model_nons_coefs)

  svydesign <- update(svydesign,
                      .y_hat_MI = y_rand_pred)

  ## inference based on mi method
  infer_nons <- nonprobMI.inference()

  return()
}


nonprobMI.fit <- function(x,
                          y,
                          weights,
                          svydesign,
                          family.outcome,
                          control.outcome,
                          verbose,
                          model,
                          x,
                          y) {


  model_nons <- stats::glm.fit(x = X_nons,
                               y = XY_nons[, 1],
                               weights = weights,
                               start = start,
                               control = control.outcome,
                               family = family.outcome)

  return(model_nons)

}

nonprobMI.inference <- function() {


}


