#' nonprobDR
#
#' nonprobDR: Function for inference based on nonprobability surveys, nonprobabilty big data sample and estimated propensoty scores.
#
#' @param selection - `formula`, the selection (propensity) equation.
#' @param outcome - `formula`, the outcome equation.
#' @param data - an optional `data.frame` with data from the nonprobability sample.
#' @param svydesign - an optional `svydesign` object (from the survey package) containing probability sample.
#' @param pop_totals - an optional `named vector` with population totals.
#' @param pop_means - an optional `named vector` with population means.
#' @param pop_size - an optional `double` with population size.
#' @param method_selection - a `character` with method for propensity scores estimation
#' @param method_outcome - a `character` with method for response variable estimation
#' @param family_selection - a `character` string describing the error distribution and link function to be used in the model. Default is "binomial". Currently only binomial with logit link is supported.
#' @param family_outcome - a `character` string describing the error distribution and link function to be used in the model. Default is "gaussian". Currently supports: gaussian with identity link, poisson and binomial.
#' @param subset - an optional `vector` specifying a subset of observations to be used in the fitting process.
#' @param weights - an optional `vector` of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector. It is assumed that this vector contains frequency or analytic weights
#' @param na_action a function which indicates what should happen when the data contain `NAs`.
#' @param control_selection a list indicating parameters to use in fitting selection model for propensity scores
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
#' @importFrom stats glm.fit
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats update
#' @importFrom stats qnorm
#' @importFrom stats binomial
#' @export


nonprobDR <- function(selection,
                      outcome,
                      data,
                      svydesign,
                      pop_totals,
                      pop_means,
                      pop_size,
                      method_selection,
                      method_outcome,
                      family_selection = "binomial",
                      family_outcome = "gaussian",
                      subset,
                      weights,
                      na_action,
                      control_selection = controlSel(),
                      control_outcome = controlOut(),
                      control_inference = controlInf(),
                      start,
                      verbose,
                      contrasts,
                      model,
                      x,
                      y,
                      ...) {

  weights <- rep.int(1, nrow(data)) # to remove

  XY_nons <- model.frame(outcome, data)
  X_nons <- model.matrix(XY_nons, data) #matrix for nonprobability sample
  X_rand <- model.matrix(selection, svydesign$variables) #matrix for probability sample
  y_nons <- XY_nons[,1]
  ps_rand <- svydesign$prob
  d_rand <- 1/ps_rand


  method <- method_selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }

  ps_method <- method$make_propen_score # function for propensity score estimation
  loglike <- method$make_log_like
  gradient <- method$make_gradient
  hessian <- method$make_hessian

  optim_method <- control_selection$optim_method

  #if(is.null(start)){

  #}

  # estimation
  model_nons <- nonprobMI_fit(x = X_nons,
                              y = y_nons,
                              weights = weights,
                              family_outcome = family_outcome)

  # initial values for propensity score estimation
  start <- start_fit(X_nons, X_rand, weights, d_rand,
                     method_selection)

  model_nons_coefs <- as.matrix(model_nons$coefficients)

  y_rand_pred <-  as.numeric(X_rand %*% model_nons_coefs) # y_hat for probability sample

  y_nons_pred <- as.numeric(X_nons %*% model_nons_coefs)

  # updating probability sample by adding y_hat variable
  svydesign <- stats::update(svydesign,
                             .y_hat_MI = y_rand_pred)

  # inference based on mi method

  nonprobDR_inference <- function(...) {


    log_like <- loglike(X_rand, X_nons, d_rand)
    gradient <- gradient(X_rand, X_nons, d_rand)
    hessian <- hessian(X_rand, X_nons, d_rand)


    ps_nons <- ps_method(X_nons, log_like, gradient, hessian, start, optim_method)$ps
    d_nons <- 1/ps_nons
    N_est_nons <- sum(d_nons)
    N_est_rand <- sum(d_rand)

    hess <- ps_method(X_nons, log_like, gradient, hessian, start, optim_method)$hess
    theta_hat <- ps_method(X_nons, log_like, gradient, hessian, start, optim_method)$theta_hat

    n_nons <- nrow(X_nons)
    n_rand <- nrow(X_rand)

    weights_nons <- 1/ps_nons

    if (method_selection == "probit") { # for probit model, propensity score derivative is required

      ps_nons_der <- ps_method(X_nons, log_like, gradient, hessian, start, optim_method)$psd
      est_ps_rand_der <- ps_method(X_rand, log_like, gradient, hessian, start, optim_method)$psd

    }

   # if(!is.null(pop_size)){

   #  N_est_rand <- pop_size
   #  N_est_nons <- pop_size
   # }


    mu_hat <- mu_hatDR(y = y_nons,
                       y_nons = y_nons_pred,
                       y_rand = y_rand_pred,
                       d_nons = d_nons,
                       d_rand = d_rand,
                       N_nons = N_est_nons,
                       N_rand = N_est_rand) #DR estimator

    h_n <- 1/N_est_nons * sum(y_nons - y_nons_pred) # errors mean

    b <- switch(method_selection,
                "logit" = (((1 - ps_nons)/ps_nons) * (y_nons - y_nons_pred - h_n)) %*% X_nons %*% solve(hess),
                "cloglog" = (((1 - ps_nons)/ps_nons^2) * (y_nons - y_nons_pred - h_n) * log(1 - ps_nons)) %*% X_nons %*% solve(hess),
                "probit" = - (ps_nons_der/ps_nons^2 * (y_nons - y_nons_pred - h_n)) %*% X_nons %*% solve(hess)
    )

    est_ps_rand <- ps_method(X_rand, log_like, gradient, hessian, start, optim_method)$ps

    pearson_residuals <- pearson_nonprobsvy(X_nons, X_rand, ps_nons, est_ps_rand) # pearson residuals for propensity score model
    deviance_residuals <- deviance_nonprobsvy(X_nons, X_rand, ps_nons, est_ps_rand) # deviance residuals for propensity score model

    # a <- 1/N_estA * sum(1 - psA) * (t(as.matrix(yA - y_estA - h_n))  %*% as.matrix(XA))


    # design based variance estimation based on approximations of the second-order inclusion probabilities

    s <- switch(method_selection,
                "logit" = as.vector(est_ps_rand) * X_rand %*% t(as.matrix(b)) + y_rand_pred - 1/N_est_nons * sum(y_nons_pred),
                "cloglog" = as.vector(log(1 - est_ps_rand)) * X_rand %*% t(as.matrix(b)) + y_rand_pred - 1/N_est_nons * sum(y_nons_pred),
                "probit" = as.vector(est_ps_rand_der/(1 - est_ps_rand)) * X_rand %*% t(as.matrix(b)) + y_rand_pred - 1/N_est_nons * sum(y_nons_pred)
    )


    ci <- n_rand/(n_rand-1) * (1 - ps_rand)
    B_hat <- sum(ci * (s/ps_rand))/sum(ci)
    ei <- (s/ps_rand) - B_hat
    db_var <- sum(ci * ei^2)

    W <- 1/N_est_nons^2 * db_var # second component

    #print(sum(is.infinite(log(1 - est_ps_rand))))

    log_likelihood <- log_like(theta_hat) # maximum of the loglikelihood function



    V <- switch(method_selection, # asymptotic variance by each propensity score method (first component)
                  "logit" = (1/N_est_nons^2) * sum((1 - ps_nons)*(((y_nons - y_nons_pred - h_n)/ps_nons) - b %*% t(X_nons))^2),
                  "cloglog" = (1/N_est_nons^2) * sum((1 - ps_nons)*(((y_nons - y_nons_pred - h_n)/ps_nons) - b %*% t(as.matrix(log((1 - ps_nons)/ps_nons) * as.data.frame(X_nons))))^2),
                  "probit" = (1/N_est_nons^2) * sum((1 - ps_nons) * (((y_nons - y_nons_pred - h_n)/ps_nons) - b %*% t(as.matrix(ps_nons_der/(ps_nons*(1 - ps_nons)) * as.data.frame(X_nons))))^2)
    )

    # typically W contains most of the variance value
    var <- V + W

    se <- sqrt(var)

    alpha <- control_inference$alpha
    z <- stats::qnorm(1-alpha/2)

    # confidence interval based on the normal approximation
    ci <- c(mu_hat - z * se, mu_hat + z * se)

    # case when samples overlap - to finish
    if (control_selection$overlap) {

      weights <- nonprobOv(X_nons,
                           X_rand,
                           d_rand,
                           dependent = control_selection$dependence,
                           method_selection)
      N <- sum(weights)

      mu_hat_Ov <-  mu_hatDR(y = y_nons,
                             y_nons = y_nons_pred,
                             y_rand = y_rand_pred,
                             d_nons = weights,
                             d_rand = d_rand,
                             N_nons = N,
                             N_rand = N_est_rand)

    }

    mu_hat <- ifelse(control_selection$overlap, mu_hat_Ov, mu_hat)

    structure(
      list(population_mean = mu_hat,
           variance = var,
           standard_error = se,
           CI = ci,
           theta = theta_hat,
           pearson_residuals = pearson_residuals,
           deviance_residuals = deviance_residuals,
           log_likelihood = log_likelihood,
           beta = model_nons_coefs),
      class = "Doubly-robust")
  }


  infer <- nonprobDR_inference()

  infer

}

#' mu_hatDR
#
#' mu_hatDR: Function for outcome variable estimation based on doubly robust estimation
#'
#' @param y - a
#' @param y_nons - a
#' @param y_rand - a
#' @param d_nons - a
#' @param d_rand - a
#' @param N_nons - a
#' @param N_rand - a

mu_hatDR <- function(y,
                     y_nons,
                     y_rand,
                     d_nons,
                     d_rand,
                     N_nons,
                     N_rand) {

  mu_hat <- 1/N_nons * sum(d_nons*(y - y_nons)) + 1/N_rand * sum(d_rand * y_rand)

  mu_hat

}





#' start_fit
#'
#' start_fit: Function for obtaining initial values for propensity score estimation
#'
#' @param X_nons - a
#' @param X_rand - a
#' @param weights - a
#' @param d - a
#' @param method_selection - a
#' @param control_selection - a

start_fit <- function(X_nons,
                      X_rand,
                      weights,
                      d,
                      method_selection,
                      control_selection = controlSel()) {

  glm_mx <- rbind(X_nons, X_rand)
  weights <- c(weights, d)
  Rnons <- c(rep(1, nrow(X_nons)), rep(0, nrow(X_rand)))

  start_model <- glm.fit(x = glm_mx, #glm model for initial values in propensity score estimation
                         y = Rnons,
                         weights = weights,
                         family = binomial(link = method_selection),
                         control = list(control_selection$epsilon,
                                        control_selection$maxit,
                                        control_selection$trace)
                         )

  start <- start_model$coefficients

  start

}

