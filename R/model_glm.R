#' Function for the mass imputation model using glm
#'
#' @importFrom stats glm.fit
#' @importFrom ncvreg cv.ncvreg
#' @importFrom stats update
#' @importFrom survey svymean
#'
#' @description
#' Modle for the outcome for the mass imputation estimator
#'
#'
#' @param y_nons target variable from non-probability sample
#' @param X_nons a `model.matrix` with auxiliary variables from non-probability sample
#' @param X_rand a `model.matrix` with auxiliary variables from non-probability sample
#' @param weights case / frequency weights from non-probability sample
#' @param svydesign a svydesign object
#' @param family_outcome family for the glm model
#' @param start_outcome start parameters
#' @param vars_selection whether variable selection should be conducted
#' @param pop_totals population totals from the `nonprob` function
#' @param pop_size population size from the `nonprob` function
#' @param control_outcome controls passed by the `control_out` function
#' @param verbose parameter passed from the main `nonprob` function
#' @param se whether standard errors should be calculated
#'
#' @returns an `nonprob_model` class which is a `list` with the following entries
#'
#' \describe{
#'   \item{model}{fitted model either an `glm.fit` or `cv.ncvreg` object}
#'   \item{y_nons_pred}{predicted values for the non-probablity sample}
#'   \item{y_rand_pred}{predicted values for the probability sample or population totals}
#'   \item{coefficients}{coefficients for the model (if available)}
#'   \item{svydesign}{an updated `surveydesign2` object (new column `y_hat_MI` is added)}
#'   \item{vars_selection}{whether variable selection was performed}
#'   \item{var_prob}{variance for the probability sample component (if available)}
#'   \item{var_nonprob}{variance for the non-probability sampl component}
#'   \item{model}{model type (character `"glm"`)}
#' }
#' @export
model_glm <- function(y_nons,
                      X_nons,
                      X_rand,
                      weights,
                      svydesign,
                      family_outcome,
                      start_outcome,
                      vars_selection,
                      pop_totals,
                      pop_size,
                      control_outcome,
                      verbose,
                      se) {

  predict.glm.fit <- function(object, newdata) {
    coefficients <- object$coefficients
    family <- object$family
    offset <- if (!is.null(object$offset)) object$offset else 0
    eta <- drop(newdata %*% coefficients) + offset
    return(family$linkinv(eta))
  }

  if (!vars_selection) {
    model_fitted <- stats::glm.fit(
      x = X_nons,
      y = y_nons,
      weights = weights,
      family = family_outcome, ## only gaussian() etc
      start = start_outcome,
      control = list(
        epsilon = control_outcome$epsilon,
        maxit = control_outcome$maxit,
        trace = control_outcome$trace
      ),
      intercept = FALSE)
    y_nons_pred <- model_fitted$fitted.values
  } else {
    model_fitted <- ncvreg::cv.ncvreg(
        X = X_nons[, -1, drop = FALSE],
        y = y_nons,
        penalty = control_outcome$penalty,
        family = family_outcome, ## here should be only character vector
        trace = verbose,
        nfolds = control_outcome$nfolds,
        nlambda = control_outcome$nlambda,
        gamma = switch(control_outcome$penalty,
                       SCAD = control_outcome$a_SCAD,
                       control_outcome$a_MCP
        ),
        lambda_min = control_outcome$lambda_min,
        eps = control_outcome$epsilon
      )

      beta_est <- model_fitted$fit$beta[, model_fitted$min]
      beta_selected <- which(abs(beta_est) != 0) - 1
      beta_est <- beta_est[model_fitted$fit$beta[, model_fitted$min] != 0]
      cve_outcome <- model_fitted$cve
      lambda_outcome <- model_fitted$lambda
      lambda_min_outcome <- model_fitted$lambda.min
    }

    if (is.null(pop_totals)) {
      y_rand_pred <- predict.glm.fit(model_fitted$glm, X_rand)
      residuals <- as.vector(y_nons - y_rand_pred)
      svydesign_updated <- stats::update(svydesign, y_hat_MI = y_rand_pred)
    } else {
      eta <- pop_totals %*% model_fitted$coefficients / pop_totals[1]
      y_rand_pred <- model_fitted$family$linkinv(eta)
    }




  ## variance components
  if (se) {

    if (is.null(pop_totals)) {
      svydesign_mean <- survey::svymean( ~ y_hat_MI, svydesign_updated)
      var_prob <- as.vector(attr(svydesign_mean, "var"))

      beta <- model_fitted$coefficients[, 1]
      eta_nons <- X_nons %*% beta
      eta_rand <- X_rand %*% beta

      mx <- 1 / pop_size * colSums(as.data.frame(X_rand) * (weights(svydesign_updated) * model_fitted$family$mu.eta(eta_rand)))
      c <- solve(1 / nrow(X_nons) * t(as.data.frame(X_nons) * model_fitted$family$mu.eta(eta_nons)) %*% X_nons) %*% mx
      var_nonprob <- as.vector(1 / nrow(X_nons)^2 * t(as.matrix(residuals^2)) %*% (X_nons %*% c)^2)

    } else {
      var_prob <- 0

      beta <- model_fitted$coefficients[, 1]
      eta_nons <- X_nons %*% beta
      eta_rand <- pop_totals %*% beta

      mx <- 1 / pop_size * pop_totals * as.vector(model_fitted$family$mu.eta(eta_rand))
      c <- solve(1 / nrow(X_nons) * t(as.data.frame(X_nons) * model_fitted$family$mu.eta(eta_nons)) %*% X_nons) %*% mx
      var_nonprob <- as.vector(1 / nrow(X_nons)^2 * t(as.matrix(residuals^2)) %*% (X_nons %*% c)^2)
      }
    } else {
      var_prob <- NA
      var_nonprob <- NA
    }

  return(
    structure(
      list(
        model = model_fitted,
        y_nons_pred = y_nons_pred,
        y_rand_pred = y_rand_pred,
        coefficients = model_fitted$coefficients,
        svydesign = if (is.null(svydesign)) svydesign else svydesign_updated,
        vars_selection = vars_selection,
        var_prob = var_prob,
        var_nonprob = var_nonprob,
        model = "glm"
    ), class = "nonprob_model")
  )
}

