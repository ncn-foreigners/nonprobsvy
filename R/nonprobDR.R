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
                      control.selection = controlSel(),
                      control.outcome = controlOut(),
                      control.inference = controlInf(),
                      start,
                      verbose,
                      contrasts,
                      model,
                      x,
                      y,
                      ...){

  weights <- rep.int(1, nrow(data)) # to remove

  XY_nons <- model.frame(outcome, data)
  X_nons <- model.matrix(XY_nons, data) #matrix for nonprobability sample
  X_rand <- model.matrix(selection, svydesign$variables) #matrix for probability sample
  y_nons = XY_nons[,1]
  ps_rand <- svydesign$prob
  d_rand <- 1/ps_rand


  method <- method.selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }

  ps_method <- method$PropenScore # function for propensity score estimation
  loglike <- method$MakeLogLike
  gradient <- method$MakeGradient
  hessian <- method$MakeHessian

  #if(is.null(start)){

  #}

  ## estimation
  model_nons <- nonprobMI.fit(x = X_nons,
                              y = y_nons,
                              weights = weights,
                              family.outcome = family.outcome)

  start <- start.fit(X_nons, X_rand, weights, d_rand,
                     method.selection) #initial values for propensity score estimation

  model_nons_coefs <- as.matrix(model_nons$coefficients)

  y_rand_pred <-  as.numeric(X_rand %*% model_nons_coefs) # y_hat for probability sample

  y_nons_pred <- as.numeric(X_nons %*% model_nons_coefs)

  svydesign <- update(svydesign,
                      .y_hat_MI = y_rand_pred) # updating probability sample by adding y_hat variable

  ## inference based on mi method

  nonprobDR.inference <- function(...){


    log_like <- loglike(X_rand, X_nons, d_rand)
    gradient <- gradient(X_rand, X_nons, d_rand)
    hessian <- hessian(X_rand, X_nons, d_rand)


    ps_nons <- ps_method(X_nons, log_like, gradient, hessian, start)$ps
    d_nons <- 1/ps_nons
    N_est_nons <- sum(d_nons)
    N_est_rand <- sum(1/ps_rand)

    hess <- ps_method(X_nons, log_like, gradient, hessian, start)$hess
    theta_hat <- ps_method(X_nons, log_like, gradient, hessian, start)$theta_hat

    n_nons <- nrow(X_nons)
    n_rand <- nrow(X_rand)

    weights_nons <- 1/ps_nons

    if(method.selection == "probit"){ # for probit model propensity score derivative is needed

      ps_nons_der <- ps_method(X_nons, log_like, gradient, hessian, start)$psd
      est_ps_rand_der <- ps_method(X_rand, log_like, gradient, hessian, start)$psd

    }



    mu_hat <- mu_hatDR(y_nons,
                       y_nons_pred,
                       y_rand_pred,
                       d_nons,
                       d_rand,
                       N_est_nons,
                       N_est_rand) #DR estimator

    h_n <- 1/N_est_nons * sum(y_nons_pred - y_nons) # errors mean

    b <- switch(method.selection,
                "logit" = (((1 - ps_nons)/ps_nons) * (y_nons - y_nons_pred - h_n)) %*% X_nons %*% solve(hess),
                "cloglog" = (((1 - ps_nons)/ps_nons^2) * (y_nons - y_nons_pred - h_n) * log(1 - ps_nons)) %*% X_nons %*% solve(hess),
                "probit" = - (ps_nons_der/ps_nons^2 * (y_nons - y_nons_pred - h_n)) %*% X_nons %*% solve(hess)
    )

    est_ps_rand <- ps_method(X_rand, log_like, gradient, hessian, start)$ps

    pearson_residuals <- pearson.residPS(X_nons, X_rand, ps_nons, est_ps_rand) # pearson residuals for propensity score model
    deviance_residuals <- deviance.residPS(X_nons, X_rand, ps_nons, est_ps_rand) # deviance residuals for propensity score model

    # a <- 1/N_estA * sum(1 - psA) * (t(as.matrix(yA - y_estA - h_n))  %*% as.matrix(XA))


    ## design based variance estimation based on approximations of the second-order inclusion probabilities

    s <- switch(method.selection,
                "logit" = as.vector(est_ps_rand) * X_rand %*% t(as.matrix(b)) + y_rand_pred - 1/N_est_nons * sum(y_nons_pred),
                "cloglog" = as.vector(log(1 - est_ps_rand)) * X_rand %*% t(as.matrix(b)) + y_rand_pred - 1/N_est_nons * sum(y_nons_pred),
                "probit" = as.vector(est_ps_rand_der/(1 - est_ps_rand)) * X_rand %*% t(as.matrix(b)) + y_rand_pred - 1/N_est_nons * sum(y_nons_pred)
    )


    ci <- n_rand/(n_rand-1) * (1 - ps_rand)
    B_hat <- sum(ci * (s/ps_rand))/sum(ci)
    ei <- (s/ps_rand) - B_hat
    db_var <- sum(ci * ei^2)

    W = 1/N_est_nons^2 * db_var # first component



    var <- switch(method.selection, # asymptotic variance by each propensity score method (first component plus second component)
                  "logit" = (1/N_est_nons^2) * sum((1 - ps_nons)*(((y_nons - y_nons_pred - h_n)/ps_nons) - b %*% t(X_nons))^2) + W,
                  "cloglog" = (1/N_est_nons^2) * sum((1 - ps_nons)*(((y_nons - y_nons_pred - h_n)/ps_nons) - b %*% t(as.matrix(log((1 - ps_nons)/ps_nons) * as.data.frame(X_nons))))^2) + W,
                  "probit" = (1/N_est_nons^2) * sum((1 - ps_nons) * (((y_nons - y_nons_pred - h_n)/ps_nons) - b %*% t(as.matrix(ps_nons_der/(ps_nons*(1 - ps_nons)) * as.data.frame(X_nons))))^2) + W,
    )


    ci <- c(mu_hat - 1.96 * sqrt(var), mu_hat + 1.96 * sqrt(var)) # confidence interval

    return(list(populationMean = mu_hat,
                Variance = var,
                CI = ci,
                theta = theta_hat,
                pearson.residuals = pearson_residuals,
                deviance.residuals = deviance_residuals
               ))
  }


  infer <- nonprobDR.inference()


  return(infer)

}

#' mu_hatDR
#
#' mu_hatDR: Function for outcome variable estimation based on doubly robust estimation
#'
#' @param ynons - a
#' @param ynons_pred - a
#' @param yrand_pred - a
#' @param dnons - a
#' @param drand - a
#' @param Nnons - a
#' @param Nrand - a

mu_hatDR <- function(ynons,
                     ynons_pred,
                     yrand_pred,
                     dnons,
                     drand,
                     Nnons,
                     Nrand){

  mu_hat <- 1/Nnons * sum(dnons*(ynons - ynons_pred)) + 1/Nrand * sum(drand * yrand_pred)
  return(mu_hat)

}





#' start.fit
#'
#' start.fit: Function for obtaining initial values for propensity score estimation
#'
#' @param X_nons - a
#' @param X_rand - a
#' @param weights - a
#' @param d - a
#' @param method.selection - a
#' @param control.selection - a

start.fit <- function(X_nons,
                      X_rand,
                      weights,
                      d,
                      method.selection,
                      control.selection = controlSel()){

  glm_mx <- rbind(X_nons, X_rand)
  weights <- c(weights, d)
  Rnons <- c(rep(1, nrow(X_nons)), rep(0, nrow(X_rand)))

  start_model <- glm.fit(x = glm_mx, #glm model for initial values in propensity score estimation
                         y = Rnons,
                         weights = weights,
                         family = binomial(link = method.selection),
                         control = list(control.selection$epsilon,
                                        control.selection$maxit,
                                        control.selection$trace)
                         )

  start <- start_model$coefficients

  return(start)

}

