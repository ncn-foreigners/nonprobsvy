# These functions are only used internally in the package, so there is no need for documenting them.
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom Matrix Matrix
#' @importFrom stats delete.response

# Selection model object
internal_selection <- function(X,
                               X_nons,
                               X_rand,
                               weights,
                               weights_rand,
                               R,
                               method_selection,
                               optim_method,
                               h = h,
                               est_method,
                               maxit,
                               varcov = FALSE,
                               ...) {

  method <- get_method(method = method_selection)

  if (est_method == "mle") {

    ps_method <- method$make_propen_score # function for propensity score estimation
    loglike <- method$make_log_like
    gradient <- method$make_gradient
    hessian <- method$make_hessian

    # initial values for propensity score estimation
    start <- start_fit(X = X,
                       R = R,
                       weight = weights,
                       weights_rand = weights_rand,
                       method_selection = method_selection)

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
  } else if (est_method == "gee"){

    inv_link <- method$make_link_inv
    h_object <- theta_h_estimation(R = R,
                                   X = X,
                                   weights_rand = weights_rand,
                                   weights = weights,
                                   h = h,
                                   method_selection = method_selection,
                                   maxit = maxit) # theta_h estimation for h_x == 2 is equal to the main method for theta estimation
    theta_hat <- h_object$theta_h
    hess <- h_object$hess
    grad <- h_object$grad
    ps_nons <- inv_link(theta_hat %*% t(as.matrix(X_nons)))
    est_ps_rand <- inv_link(theta_hat %*% t(as.matrix(X_rand)))

    if (method_selection == "probit") { # for probit model, propensity score derivative is required
      dinv_link <- method$make_link_inv_der
      ps_nons_der <- dinv_link(theta_hat %*% t(as.matrix(X_nons)))
      est_ps_rand_der <- dinv_link(theta_hat %*% t(as.matrix(X_rand)))
    }

    list(theta_hat = theta_hat,
         hess = hess,
         grad = grad,
         ps_nons = ps_nons,
         est_ps_rand = est_ps_rand,
         ps_nons_der = ifelse(method_selection == "probit", ps_nons_der, NA),
         est_ps_rand_der = ifelse(method_selection == "probit", est_ps_rand_der, NA),
         var_cov1 = ifelse(varcov, method$variance_covariance1, "No variance-covariance matrix"),
         var_cov2 = ifelse(varcov, method$variance_covariance2, "No variance-covariance matrix"))
  }

}
# Outcome model object
internal_outcome <- function(X_nons,
                             X_rand,
                             y,
                             weights,
                             family_outcome,
                             pop_totals = FALSE) {

  # estimation
  model_nons <- nonprobMI_fit(x = X_nons,
                              y = y,
                              weights = weights,
                              family_outcome = family_outcome)


  model_nons_coefs <- model_nons$coefficients

  if (pop_totals) {
    y_rand_pred <- sum(X_rand * model_nons_coefs)
  } else {
     y_rand_pred <-  as.numeric(X_rand %*% model_nons_coefs) # y_hat for probability sample
  }
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

  p <- ncol(X)
  #start0 <- rep(0, p)
  start0 <- start_fit(X = X,
                      R = R,
                      weights = weights,
                      weights_rand = weights_rand,
                      method_selection = method_selection)
  # theta estimation by unbiased estimating function depending on the h_x function TODO
  u_theta <- u_theta(R = R,
                     X = X,
                     weights = c(weights_rand, weights),
                     h = h,
                     method_selection = method_selection,
                     pop_totals = pop_totals)

  u_theta_der <- u_theta_der(R = R,
                             X = X,
                             weights = c(weights_rand, weights),
                             h = h,
                             method_selection = method_selection,
                             pop_totals = pop_totals)


  for (i in 1:maxit) {
    start <- start0 + MASS::ginv(u_theta_der(start0)) %*% u_theta(start0) # consider solve function
    if (sum(abs(start - start0)) < 0.001) break;
    if (sum(abs(start - start0)) > 1000) break;
    start0 <- start
  }
  theta_h <- as.vector(start)
  # opt <- rootSolve::multiroot(f = u_theta, start = start0) <---- just for tests


  list(theta_h = theta_h,
       hess = u_theta_der(theta_h),
       grad = u_theta(theta_h))
}
# Variance for inverse probability weighted estimator
internal_varIPW <- function(X_nons,
                            X_rand,
                            y_nons,
                            ps_nons,
                            mu_hat,
                            hess,
                            ps_nons_der,
                            N,
                            est_ps_rand,
                            ps_rand,
                            est_ps_rand_der,
                            n_rand,
                            pop_size,
                            method_selection,
                            est_method,
                            theta,
                            h,
                            var_cov1 = var_cov1,
                            var_cov2 = var_cov2) {

  eta <- as.vector(X_nons %*% as.matrix(theta))
  method <- get_method(method_selection)
  b_obj <- method$b_vec_ipw(X = X_nons,
                            ps = ps_nons,
                            psd = ps_nons_der,
                            y = y_nons,
                            mu = mu_hat,
                            hess = hess,
                            eta = eta,
                            pop_size = pop_size)
  b <- b_obj$b
  hess_inv <- b_obj$hess_inv

  # sparse matrix
  b_vec <- cbind(-1, b)
  H_mx <- cbind(0, N * hess_inv)
  sparse_mx <- Matrix::Matrix(rbind(b_vec, H_mx), sparse = TRUE)

  if (method_selection == "probit") { # change this chunk of code - condition is redundant

    V1 <- var_cov1(X = X_nons,
                   y = y_nons,
                   mu = mu_hat,
                   ps = ps_nons,
                   psd = ps_nons_der,
                   pop_size = pop_size,
                   est_method = est_method,
                   h = h) # fixed
    V2 <- var_cov2(X = X_rand,
                   eps = est_ps_rand,
                   ps = ps_rand,
                   psd = est_ps_rand_der,
                   n = n_rand,
                   N = N,
                   est_method = est_method,
                   h = h)

  } else {

    V1 <- var_cov1(X = X_nons,
                   y = y_nons,
                   mu = mu_hat,
                   ps = ps_nons,
                   pop_size = pop_size,
                   est_method = est_method,
                   h = h) # fixed
    V2 <- var_cov2(X = X_rand,
                   eps = est_ps_rand,
                   ps = ps_rand,
                   n = n_rand,
                   N = N,
                   est_method = est_method,
                   h = h)

  }

  # variance-covariance matrix for set of parameters (mu_hat and theta_hat)
  V_mx_nonprob <- sparse_mx %*% V1 %*% t(as.matrix(sparse_mx)) # nonprobability component
  V_mx_prob <- sparse_mx %*% V2 %*% t(as.matrix(sparse_mx)) # probability component
  V_mx <- V_mx_nonprob + V_mx_prob

  var_nonprob <- as.vector(V_mx_nonprob[1,1])
  var_prob <- as.vector(V_mx_prob[1,1])
  var <- as.vector(V_mx[1,1])
  # vector of variances for theta_hat
  theta_hat_var <- diag(as.matrix(V_mx[2:ncol(V_mx), 2:ncol(V_mx)]))

  list(var_nonprob = var_nonprob,
       var_prob = var_prob,
       var = var,
       theta_hat_var = theta_hat_var)
}
# Variance for doubly robust estimator
internal_varDR <- function(OutcomeModel,
                           SelectionModel,
                           y_nons_pred,
                           method_selection,
                           theta,
                           ps_nons,
                           hess,
                           ps_nons_der,
                           est_ps_rand,
                           y_rand_pred,
                           N_nons,
                           est_ps_rand_der,
                           svydesign,
                           est_method,
                           h) {

  eta <- as.vector(SelectionModel$X_nons %*% as.matrix(theta))
  h_n <- 1/N_nons * sum(OutcomeModel$y_nons - y_nons_pred) # errors mean
  method <- get_method(method_selection)
  est_method <- get_method(est_method)
  #psd <- method$make_link_inv_der(eta)

  b <- method$b_vec_dr(X = SelectionModel$X_nons,
                       ps = ps_nons,
                       psd = ps_nons_der,
                       y = OutcomeModel$y_nons,
                       mu = mu_hat,
                       hess = hess,
                       eta = eta,
                       h_n = h_n,
                       y_pred = y_nons_pred)

  t <- est_method$make_t(X = SelectionModel$X_rand,
                         ps = est_ps_rand,
                         psd = est_ps_rand_der,
                         b = b,
                         h = h,
                         y_rand = y_rand_pred,
                         y_nons = y_nons_pred,
                         N = N_nons,
                         method_selection = method_selection)
  # asymptotic variance by each propensity score method (nonprobability component)
  var_nonprob <- est_method$make_var_nonprob(ps = ps_nons,
                                             psd = ps_nons_der,
                                             y = OutcomeModel$y_nons,
                                             y_pred = y_nons_pred,
                                             h_n = h_n,
                                             X = SelectionModel$X_nons,
                                             b = b,
                                             N = N_nons,
                                             h = h,
                                             method_selection = method_selection)



  # design based variance estimation based on approximations of the second-order inclusion probabilities
  svydesign <- stats::update(svydesign,
                             t = t)
  svydesign_mean <- survey::svymean(~t, svydesign) #perhaps using survey package to compute prob variance
  var_prob <- as.vector(attr(svydesign_mean, "var"))

  list(var_prob = var_prob,
       var_nonprob = var_nonprob)
}
# create an object with model frames and matrices to preprocess
model_frame <- function(formula, data, svydesign = NULL, pop_totals = NULL, pop_size = NULL) {

  if (!is.null(svydesign)) {
  XY_nons <- model.frame(formula, data)
  X_nons <- model.matrix(XY_nons, data) #matrix for nonprobability sample with intercept
  nons_names <- attr(terms(formula, data = data), "term.labels")
  if (all(nons_names %in% colnames(svydesign$variables))) {
    X_rand <- model.matrix(delete.response(terms(formula)), svydesign$variables) #X_rand <- as.matrix(cbind(1, svydesign$variables[,nons_names])) #matrix of probability sample with intercept
  } else {
    stop("variable names in data and svydesign do not match")
  }
  y_nons <- XY_nons[,1]

  list(X_nons = X_nons,
       X_rand = X_rand,
       nons_names = nons_names,
       y_nons = y_nons)
  } else if (!is.null(pop_totals)) {
    XY_nons <- model.frame(formula, data)
    X_nons <- model.matrix(XY_nons, data) #matrix for nonprobability sample with intercept
    nons_names <- attr(terms(formula, data = data), "term.labels")
    if(all(nons_names %in% names(pop_totals))) { # pop_totals, pop_means defined such as in `calibrate` function
      pop_totals <-  pop_totals[nons_names]
    } else {
      warning("Selection and population totals have different names.")
    }
    y_nons <- XY_nons[,1]

    list(X_nons = X_nons,
         pop_totals = pop_totals,
         nons_names = nons_names,
         y_nons = y_nons)
  }
}
# Function for getting function from the selected method
get_method <- function(method) {
  if (is.character(method)) {
    method <- get(method, mode = "function", envir = parent.frame())
  }
  if (is.function(method)) {
    method <- method()
  }
  method
}
