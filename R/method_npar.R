#' Function for the mass imputation model using nonparametric method
#'
#' @importFrom stats glm.fit
#' @importFrom ncvreg cv.ncvreg
#' @importFrom stats update
#' @importFrom survey svymean
#'
#' @description
#' Model for the outcome for the mass imputation estimator
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
#' @returns an `nonprob_method` class which is a `list` with the following entries
#'
#' \describe{
#'   \item{model_fitted}{fitted model either an `glm.fit` or `cv.ncvreg` object}
#'   \item{y_nons_pred}{predicted values for the non-probablity sample}
#'   \item{y_rand_pred}{predicted values for the probability sample or population totals}
#'   \item{coefficients}{coefficients for the model (if available)}
#'   \item{svydesign}{an updated `surveydesign2` object (new column `y_hat_MI` is added)}
#'   \item{y_mi_hat}{estimated population mean for the target variable}
#'   \item{vars_selection}{whether variable selection was performed}
#'   \item{var_prob}{variance for the probability sample component (if available)}
#'   \item{var_nonprob}{variance for the non-probability sampl component}
#'   \item{model}{model type (character `"npar"`)}
#' }
#' @export
model_npar <- function(y_nons,
                       X_nons,
                       X_rand,
                       svydesign,
                       weights=NULL,
                       family_outcome="gaussian",
                       start_outcome=NULL,
                       vars_selection=FALSE,
                       pop_totals=NULL,
                       pop_size=NULL,
                       control_outcome=control_out(),
                       control_inference=control_inf(),
                       verbose=FALSE,
                       se=TRUE) {

  ## this does not handle aggregated data

  if (missing(y_nons) | missing(X_nons)) {
    stop("Arguments `y_nons` and `X_nons` are required.")
  }

  if (is.null(svydesign)) svydesign <- NULL
  if (is.null(start_outcome)) start_outcome <- numeric(ncol(X_nons))
  if (is.null(weights)) weights <- rep(1, nrow(X_nons))
  if (is.null(pop_size)) pop_size <- sum(weights(svydesign))

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
      family = get(family_outcome)(),
      start = start_outcome,
      control = list(
        epsilon = control_outcome$epsilon,
        maxit = control_outcome$maxit,
        trace = control_outcome$trace
      ),
      intercept = FALSE)

    if (is.null(pop_totals)) {
      y_nons_pred <- predict.glm.fit(model_fitted, X_nons)
      y_rand_pred <- predict.glm.fit(model_fitted, X_rand)
      residuals <- as.vector(y_nons - y_nons_pred)
      svydesign_updated <- stats::update(svydesign, y_hat_MI = y_rand_pred)
      svydesign_mean <- survey::svymean( ~ y_hat_MI, svydesign_updated)
      y_mi_hat <- as.numeric(svydesign_mean)
    } else {
      y_nons_pred <- predict.glm.fit(model_fitted, X_nons)
      y_rand_pred <- NA
      residuals <- as.vector(y_nons - y_nons_pred)
      eta <- drop(pop_totals %*% model_fitted$coefficients / pop_totals[1])
      y_mi_hat <- model_fitted$family$linkinv(eta)
    }
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

    model_fitted$family <- get(model_fitted$fit$family)()
    #beta_est <- model_fitted$fit$beta[, model_fitted$min]
    #beta_selected <- unname(which(abs(beta_est) > 0))
    #beta_est <- beta_est[beta_selected]

    if (is.null(pop_totals)) {
      y_nons_pred <- predict(model_fitted, X_nons[, -1], type = "response")
      y_rand_pred <- predict(model_fitted, X_rand[, -1], type = "response")
      residuals <- as.vector(y_nons - y_nons_pred)
      svydesign_updated <- stats::update(svydesign, y_hat_MI = y_rand_pred)
      svydesign_mean <- survey::svymean( ~ y_hat_MI, svydesign_updated)
      y_mi_hat <- as.numeric(svydesign_mean)

    } else {
      y_nons_pred <- predict(model_fitted, X_nons[, -1], type = "response")
      y_rand_pred <- NA
      residuals <- as.vector(y_nons - y_nons_pred)
      eta <- drop(pop_totals %*% coef(model_fitted) / pop_totals[1])
      y_mi_hat <- model_fitted$family$linkinv(eta)
    }
  }


  ## variance components
  if (se) {

    if (is.null(pop_totals)) {

      var_prob <- as.vector(attr(svydesign_mean, "var"))

      beta <- coef(model_fitted)
      eta_nons <- drop(X_nons %*% beta)
      eta_rand <- drop(X_rand %*% beta)

      mx <- 1 / pop_size * colSums(X_rand * (weights(svydesign_updated) * model_fitted$family$mu.eta(eta_rand)))
      c <- solve(1 / nrow(X_nons) * t(X_nons * model_fitted$family$mu.eta(eta_nons)) %*% X_nons) %*% mx
      var_nonprob <- drop(1 / nrow(X_nons)^2 * t(as.matrix(residuals^2)) %*% (X_nons %*% c)^2)

    } else {
      var_prob <- 0

      beta <- coef(model_fitted)
      eta_nons <- drop(X_nons %*% beta)
      eta_rand <- drop(pop_totals %*% beta)

      mx <- 1 / pop_size * pop_totals * as.vector(model_fitted$family$mu.eta(eta_rand))
      c <- solve(1 / nrow(X_nons) * t(X_nons * model_fitted$family$mu.eta(eta_nons)) %*% X_nons) %*% mx
      var_nonprob <- as.vector(1 / nrow(X_nons)^2 * t(as.matrix(residuals^2)) %*% (X_nons %*% c)^2)
    }
  } else {
    var_prob <- NA
    var_nonprob <- NA
  }

  return(
    structure(
      list(
        model_fitted = model_fitted,
        y_nons_pred = y_nons_pred,
        y_rand_pred = y_rand_pred,
        coefficients = coef(model_fitted),
        svydesign = if (is.null(svydesign)) svydesign else svydesign_updated,
        y_mi_hat = y_mi_hat,
        vars_selection = vars_selection,
        var_prob = var_prob,
        var_nonprob = var_nonprob,
        var_total = var_prob + var_nonprob,
        model = "glm",
        family = family_outcome
      ),
      class = "nonprob_method")
  )
}

