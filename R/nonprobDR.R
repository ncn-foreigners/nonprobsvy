#' nonprobDR
#
#' nonprobDR: Function for inference based on nonprobability surveys, nonprobabilty big data sample and estimated propensoty scores.
#
#' @param selection - `formula`, the selection (propensity) equation.
#' @param outcome - `formula`, the outcome equation.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param pop.totals - an optional `named vector` with population totals.
#' @param pop.means - an optional `named vector` with population means.
#' @param pop.size - an optional `double` with population size.
#' @param method.selection - a `character` with method for propensity scores estimation
#' @param method.outcome - a `character` with method for response variable estimation
#' @param family.selection - a `character` string describing the error distribution and link function to be used in the model. Default is "binomial". Currently only binomial with logit link is supported.
#' @param family.outcome - a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param subset - an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param weights - an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na.action a
#' @param control.selection a
#' @param control.outcome a
#' @param control.inference a
#' @param start a
#' @param verbose a
#' @param contrasts a
#' @param model a
#' @param x a
#' @param y a
#' @param ... a
#'
#' @importFrom stats glm.fit
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @export


nonprobDR <- function(selection,
                      outcome,
                      data,
                      svydesign,
                      pop.totals,
                      pop.means,
                      pop.size,
                      method.selection,
                      method.outcome,
                      family.selection = "binomial",
                      family.outcome = "gaussian",
                      subset,
                      weights,
                      na.action,
                      control.selection,
                      control.outcome,
                      control.inference,
                      start,
                      verbose,
                      contrasts,
                      model,
                      x,
                      y,
                      ...){


  XY_nons <- model.frame(outcome, data)
  X_nons <- model.matrix(XY_nons, data)
  X_rand <- model.matrix(outcome, svydesign$variables)
  y_nons = XY_nons[,1]
  ps_rand <- svydesign$prob
  d_rand <- 1/ps_rand
  ps_method <- method.selection$PropenScore
  loglike <- method.selection$MakeLogLike
  gradient <- method.selection$MakeGradient
  hessian <- method.selection$MakeHessian

  if(is.null(start)){



  }

  ## estimation
  model_nons <- nonprobMI.fit()

  model_nons_coefs <- as.matrix(model_nons$coefficients)

  y_rand_pred <-  as.numeric(X_rand %*% model_nons_coefs) # y_hat for probability sample

  y_nons_pred <- as.numeric(X_nons %*% model_nons_coefs)

  svydesign <- update(svydesign,
                      .y_hat_MI = y_rand_pred) # updating probability sample by adding y_hat variable

  ## inference based on mi method

  infer <- nonprobDR.inference()


  return(infer)

}

nonprobMI.fit <- function(outcome,
                          data,
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

nonprobDR.inference <- function(...){


  log_like <- log_like(X_rand, X_nons, d_rand)
  gradient <- gradient(X_rand, X_nons, d_rand)
  hessian <- hessian(X_rand, X_nons, d_rand)

  ps_nons <- ps_method(X_nons, log_like, gradient, hessian)$ps
  d_nons <- 1/ps_nons
  N_est_nons <- sum(d_nons)
  N_est_rand <- sum(1/ps_rand)

  hess <- ps_method(X_nons, log_like, gradient, hessian)$hess
  theta_hat <- ps_method(X_nons, log_like, gradient, hessian)$theta_hat

  n_nons <- nrow(X_nons)
  n_rand <- nrow(XB_rand)

  weights_nons <- 1/ps_nons

  if(method.selection == "probit"){

    ps_nons_der <- ps_method(X_nons, log_like, gradient, hessian)$psd
    ps_nons <- ps_method(X_nons, log_like, gradient, hessian)$ps

  }


  mu_hat <- 1/N_est_nons * sum(d_nond*(y_nons - y_nons_pred)) + 1/N_est_rand * sum(ps_rand * y_rand_pred)
  h_n <- 1/N_est_nons * sum(y_nons_pred - y_nons)

  b <- switch(method.selection,
                  "logit" = (((1 - ps_nons)/ps_nons) * (y_nons - y_nons_pred - h_n)) %*% as.matrix(X_nons) %*% solve(hess),
                  "cloglog" = (((1 - ps_nons)/ps_nons^2) * (y_nons - y_nons_pred - h_n) * log(1 - ps_nons)) %*% as.matrix(X_nons) %*% solve(hess),
                  "probit" = - (ps_nons_der/ps_nons^2 * (y_nons - y_nons_pred - h_n)) %*% as.matrix(X_nons) %*% solve(hess)
              )

  est_ps_rand <- ps_method(X_rand, log_like, gradient, hessian)$ps

      # a <- 1/N_estA * sum(1 - psA) * (t(as.matrix(yA - y_estA - h_n))  %*% as.matrix(XA))


  s <- as.vector(est_ps_rand/(1 - est_ps_rand)) * as.matrix(X_rand) %*% t(as.matrix(b)) + y_rand_pred - 1/N_est_nons * sum(y_nons_pred)
  ci <- n_rand/(n_rand-1) * (1 - est_ps_rand)
  B_hat <- sum(ci * (s/ps_rand))/sum(ci)
  ei <- (s/psB) - B_hat
  db_var <- sum(ci * ei^2)

  W = 1/N_est_nons^2 * db_var



  var <- switch(method.selection,
                 "logit" = (1/N_est_nons^2) * sum((1 - ps_nons)*(((y_nons - y_nons_pred - h_n)/ps_nons) - b %*% t(as.matrix(X_nons)))^2) + W,
                 "cloglog" = (1/N_est_nons^2) * sum((1 - ps_nons)*(((y_nons - y_nons_pred - h_n)/ps_nons) - b %*% t(as.matrix(log((1 - ps_nons)/ps_nons)*X_nons)))^2) + W,
                 "probit" = (1/N_est_nons^2) * sum((1 - ps_nons) * (((y_nons - y_nons_pred - h_n)/ps_nons) - b %*% t(as.matrix(ps_nons_der/(ps_nons*(1 - ps_nons))*X_nons)))^2) + W,
  )





  ci <- c(mu_hat - 1.96 * sqrt(var), mu_hat + 1.96 * sqrt(var))

  return(list("Population mean estimator" = mu_hat, "variance" = var, "CI" = ci))
}


