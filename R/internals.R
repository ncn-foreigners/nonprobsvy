# These functions are only used internally in the package, so there is no need for documenting them

internal_selection <- function(X,
                               X_nons,
                               X_rand,
                               weights,
                               weights_rand,
                               R,
                               method_selection,
                               optim_method,
                               varcov = FALSE,
                               ...) {

  method <- get_method(method_selection)
  ps_method <- method$make_propen_score # function for propensity score estimation
  loglike <- method$make_log_like
  gradient <- method$make_gradient
  hessian <- method$make_hessian

  # initial values for propensity score estimation
  start <- start_fit(X,
                     R,
                     weights,
                     weights_rand,
                     method_selection)

  log_like <- loglike(X_nons,
                      X_rand,
                      weights_rand)

  gradient <- gradient(X_nons,
                       X_rand,
                       weights_rand)

  hessian <- hessian(X_nons,
                     X_rand,
                     weights_rand)

  maxLik_nons_obj <- ps_method(X_nons,
                               log_like,
                               gradient,
                               hessian,
                               start,
                               optim_method)

  maxLik_rand_obj <- ps_method(X_rand,
                               log_like,
                               gradient,
                               hessian,
                               start,
                               optim_method)

  theta <- maxLik_nons_obj$theta_hat
  log_likelihood <- log_like(theta)

  list(maxLik_rand_obj = maxLik_rand_obj,
       maxLik_nons_obj = maxLik_nons_obj,
       theta = theta,
       log_likelihood = log_likelihood,
       var_cov1 = ifelse(varcov, method$variance_covariance1, "No variance-covariance matrix"),
       var_cov2 = ifelse(varcov, method$variance_covariance2, "No variance-covariance matrix"))

}

internal_outcome <- function(X_nons,
                             X_rand,
                             y,
                             weights,
                             family_outcome,
                              ...) {

  # estimation
  model_nons <- nonprobMI_fit(x = X_nons,
                              y = y,
                              weights = weights,
                              family_outcome = family_outcome)


  model_nons_coefs <- model_nons$coefficients

  y_rand_pred <-  as.numeric(X_rand %*% model_nons_coefs) # y_hat for probability sample
  y_nons_pred <- as.numeric(X_nons %*% model_nons_coefs)

  list(y_rand_pred = y_rand_pred,
       y_nons_pred = y_nons_pred,
       model_nons_coefs = model_nons_coefs)

}

theta_h_estimation <- function(R,
                               X,
                               weights_rand,
                               weights,
                               h,
                               method_selection,
                               maxit,
                               pop_totals = NULL,
                               pop_means = NULL){

  # theta estimation by unbiased estimating function depending on the h_x function TODO
  u_theta <- u_theta(R = R, X = X,
                     weights = c(weights_rand, weights), h = h,
                     method_selection = method_selection)

  u_theta_der <- u_theta_der(R = R, X = X,
                             weights = c(weights_rand, weights), h = h,
                             method_selection = method_selection)
  p <- ncol(X)
  start0 <- rep(0, p)
  for (i in 1:maxit) {
    start <- start0 + MASS::ginv(u_theta_der(start0)) %*% u_theta(start0)
    if (sum(abs(start - start0)) < 0.001) break;
    if (sum(abs(start - start0)) > 1000) break;
    start0 <- start
  }
  theta_h <- as.vector(start)
  theta_h
}


model_frame <- function(formula, data, svydesign) {

  XY_nons <- model.frame(formula, data)
  X_nons <- model.matrix(XY_nons, data) #matrix for nonprobability sample with intercept
  nons_names <- attr(terms(formula, data = data), "term.labels")
  if (all(nons_names %in% colnames(svydesign$variables))) {
    X_rand <- as.matrix(cbind(1, svydesign$variables[,nons_names])) #matrix of probability sample with intercept
  } else {
    stop("variable names in data and svydesign do not match")
  }
  y_nons <- XY_nons[,1]

  list(X_nons = X_nons,
       X_rand = X_rand,
       nons_names = nons_names,
       y_nons = y_nons)
}


get_method <- function(method_selection) {
  method <- method_selection
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }
}
