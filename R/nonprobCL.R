#' nonprobCL
#
#' nonprobCL: Function for inference based on nonprobability surveys using calibration approach
#
#' @param outcome - `formula`, the outcome equation.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param family_outcome - a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param method_outcome - a `character` with method for response variable estimation
#' @param subset - an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param strata - an optional `vector` specifying strata.
#' @param weights - an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na_action a function which indicates what should happen when the data contain `NAs`.
#' @param control_outcome a list indicating parameters to use in fitting model for outcome variable
#' @param control_inference a list indicating parameters to use in inference based on probablity and nonprobability samples, contains parameters such as estimation method or variance method
#' @param start - an optional `list` with starting values for the parameters of the selection and outcome equation
#' @param verbose - verbose, numeric
#' @param contrasts a
#' @param model a
#' @param x a
#' @param y a
#' @param ... a
#'
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats update
#' @importFrom sampling calib
#' @export

nonprobCL <- function(outcome,
                      data,
                      svydesign,
                      method_outcome,
                      family_outcome = "gaussian",
                      subset,
                      strata,
                      weights,
                      na_action,
                      control_outcome = controlOut(),
                      control_inference = controlInf(var_method = "analytic"),
                      start,
                      verbose,
                      contrasts,
                      model,
                      x,
                      y,
                      ...) {

  weights <- rep.int(1, nrow(data)) # to remove

  XY_nons <- model.frame(outcome, data)
  y_name <- colnames(XY_nons)[1]
  X_nons <- model.matrix(XY_nons, data) # matrix of the  nonprobability sample
  svydesign$variables[,y_name] <- rep(0, nrow(svydesign$variables))
  X_rand <- model.matrix(outcome, svydesign$variables) # matrix of the probability sample
  y_nons <- XY_nons[,1]

  R_nons <- rep(1, nrow(X_nons))
  R_rand <- rep(0, nrow(X_rand))
  R <- c(R_nons, R_rand)

  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)

  n_nons <- nrow(X_nons)
  n_rand <- nrow(X_rand)
  X <- rbind(X_nons, X_rand)

  ps_rand <- svydesign$prob
  weights_rand <- 1/ps_rand
  N_est_rand <- sum(weights_rand)

  nonprobCL_inference <- function(...) {


    total <- t(X_rand) %*% as.matrix(weights_rand)
    weights_nons <- sampling::calib(Xs = X_nons, d = weights, total = total, method = "linear")
    N_est_nons <- sum(weights_nons)
    mu_hat <- 1/N_est_nons * sum(weights_nons * y_nons)

    structure(
      list(mean = mu_hat
            ),
      class = "Calibration Weighting")

  }

  #
  infer_nons <- nonprobCL_inference()

  infer_nons




}
