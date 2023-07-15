# Estimation methods for integrated sources

# Object with output parameters for Maximum likelihood Estimation for propensity scores
mle <- function(...) {

  estimation_model <- function(model, method_selection, ...) {
    method <- get_method(method_selection)
    dinv_link <- method$make_link_inv_der
    maxLik_nons_obj <- model$maxLik_nons_obj
    log_likelihood <- maxLik_nons_obj$log_l # maximum of the loglikelihood function
    theta_hat <- model$theta

    ps_nons <- model$ps
    ps_nons_der <- model$ps_der
    est_ps_rand <- model$ps_rand
    est_ps_rand_der <- model$ps_rand_der
    hess <- maxLik_nons_obj$hess
    grad <- maxLik_nons_obj$grad
    var_cov1 <- model$var_cov1
    var_cov2 <- model$var_cov2
    df_residual <- model$df_residual
    variance_covariance <- solve(-hess) # variance-covariance matrix of estimated parameters

    list(theta_hat = theta_hat,
         grad = grad,
         hess = hess,
         var_cov1 = var_cov1,
         var_cov2 = var_cov2,
         ps_nons = ps_nons,
         est_ps_rand = est_ps_rand,
         ps_nons_der = ps_nons_der,
         est_ps_rand_der = est_ps_rand_der,
         variance_covariance = variance_covariance,
         log_likelihood = log_likelihood,
         df_residual = df_residual)
  }

  make_t <- function(X, ps, psd, b, y_rand, y_nons, h, N, method_selection, weights, weights_sum) {
    method <- get_method(method_selection)
    t <- method$t_vec(X = X,
                      ps = ps,
                      psd = psd,
                      b = b,
                      y_rand = y_rand,
                      y_nons = y_nons,
                      N = N,
                      weights = weights,
                      weights_sum = weights_sum)
    t
  }

  make_var_nonprob <- function(ps, psd, y, y_pred, h_n, X, b, N, h, method_selection, weights = weights, weights_sum, pop_totals = NULL) {
    method <- get_method(method_selection)
    var_nonprob <-  method$var_nonprob(ps = ps,
                                       psd = psd,
                                       y = y,
                                       y_pred = y_pred,
                                       h_n = h_n,
                                       X = X,
                                       b = b,
                                       N = N,
                                       weights = weights,
                                       weights_sum = weights_sum)
    as.numeric(var_nonprob)

  }

  model_selection <- function(X,
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
                              control_selection,
                              ...) {

    method <- get_method(method = method_selection)
    max_lik <- method$make_max_lik # function for propensity score estimation
    loglike <- method$make_log_like
    gradient <- method$make_gradient
    hessian <- method$make_hessian
    inv_link <- method$make_link_inv
    dinv_link <- method$make_link_inv_der

    # initial values for propensity score estimation
    start <- start_fit(X = X,
                       R = R,
                       weights = weights,
                       weights_rand = weights_rand,
                       method_selection = method_selection)

    df_reduced <- nrow(X) - length(start)

    maxLik_nons_obj <- max_lik(X_nons = X_nons,
                               X_rand = X_rand,
                               weights = weights,
                               weights_rand = weights_rand,
                               start = start,
                               control = control_selection)

    theta <- maxLik_nons_obj$theta_hat
    eta_nons <- theta %*% t(X_nons)
    eta_rand <- theta %*% t(X_rand)

    ps_nons <- inv_link(eta_nons)
    est_ps_rand <- inv_link(eta_rand)

    ps_nons_der <- dinv_link(eta_nons)
    est_ps_rand_der <- dinv_link(eta_rand)

    list(maxLik_nons_obj = maxLik_nons_obj,
         theta = theta,
         ps = ps_nons,
         ps_der = ps_nons_der,
         ps_rand = est_ps_rand,
         ps_rand_der = est_ps_rand_der,
         var_cov1 = ifelse(varcov, method$variance_covariance1, "No variance-covariance matrix"),
         var_cov2 = ifelse(varcov, method$variance_covariance2, "No variance-covariance matrix"),
         df_residual = df_reduced)
  }
  structure(
    list(estimation_model = estimation_model,
         make_t = make_t,
         make_var_nonprob = make_var_nonprob,
         model_selection = model_selection),
    class = "method"
  )
}

# Object with output parameters for estimation by Generalized Estimating Equations for propensity scores
gee <- function(...) {

  estimation_model <- function(model, method_selection) {

    theta_hat <- model$theta_hat
    hess <- model$hess
    grad <- model$grad
    ps_nons <- model$ps_nons
    est_ps_rand <- model$est_ps_rand
    ps_nons_der <- model$ps_nons_der
    est_ps_rand_der <- model$est_ps_rand_der
    var_cov1 <- model$var_cov1
    var_cov2 <- model$var_cov2
    df_residual <- model$df_residual
    variance_covariance <- model$variance_covariance # variance-covariance matrix of estimated parameters

    list(theta_hat = theta_hat,
         grad = grad,
         hess = hess,
         var_cov1 = var_cov1,
         var_cov2 = var_cov2,
         ps_nons = ps_nons,
         est_ps_rand = est_ps_rand,
         ps_nons_der = ps_nons_der,
         est_ps_rand_der = est_ps_rand_der,
         variance_covariance = variance_covariance,
         df_residual = df_residual,
         log_likelihood = "NULL")
  }

  make_t <- function(X, ps, psd, b, y_rand, y_nons, h, N, method_selection, weights, weights_sum) {
    if (h == "1") {
      t <- X %*% t(as.matrix(b)) + y_rand - 1/N * sum(weights * y_nons)
    } else if (h == "2") {
      t <- as.vector(ps) * X %*% t(as.matrix(b)) + y_rand - 1/N * sum(weights * y_nons)
    }
    t
  }

  make_var_nonprob <- function(ps, psd, y, y_pred, h_n, X, b, N, h, method_selection, weights, weights_sum, pop_totals) {
    if (!is.null(pop_totals)) h <- "1" # perhaps to remove, just check if appropriate var is calculated
    if (h == "2") {
      var_nonprob <- 1/N^2 * sum((1 - ps) * ((weights*(y - y_pred - h_n)/ps) - b %*% t(X))^2)
    } else if (h == "1") {
      var_nonprob <- 1/N^2 * sum((1 - ps) * ((weights*(y - y_pred - h_n) - b %*% t(X))/ps)^2)
    }
    as.numeric(var_nonprob)
  }

  model_selection <- function(X,
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
                              ...){

    method <- get_method(method = method_selection)
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
    variance_covariance <- solve(-hess)

    # just for testing
    # ps <- c(ps_nons, est_ps_rand)
    # V <- diag(ps * (1 - ps))
    # vcov_selection <- solve(t(X) %*% V %*% X)
    # print(sqrt(diag(vcov_selection)))


    df_reduced <- nrow(X) - length(theta_hat)

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
         variance_covariance = variance_covariance,
         var_cov1 = ifelse(varcov, method$variance_covariance1, "No variance-covariance matrix"),
         var_cov2 = ifelse(varcov, method$variance_covariance2, "No variance-covariance matrix"),
         df_residual = df_reduced)
  }

  structure(
    list(estimation_model = estimation_model,
         make_t = make_t,
         make_var_nonprob = make_var_nonprob,
         model_selection = model_selection),
    class = "method"
  )

}
