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
    variance_covariance <- solve(-hess) # MASS::ginv # variance-covariance matrix of estimated parameters

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
                      weights = weights)
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
                                       weights = weights)
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

  make_t <- function(X, ps, psd, b, y_rand, y_nons, h, N, method_selection, weights) {
    if (h == "1") {
      t <- X %*% t(as.matrix(b)) + y_rand - 1/N * sum(weights * y_nons)
    } else if (h == "2") {
      t <- as.vector(ps) * X %*% t(as.matrix(b)) + y_rand - 1/N * sum(weights * y_nons)
    }
    t
  }

  make_var_nonprob <- function(ps, psd, y, y_pred, h_n, X, b, N, h, method_selection, weights, pop_totals) {
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

mm <- function(X, y, weights, weights_rand, R, n_nons, n_rand, method_selection, family, boot = FALSE) { # TODO

  method <- get_method(method_selection)
  inv_link <- method$make_link_inv
  dinv_link <- method$make_link_inv_der

  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)

  p <- ncol(X)
  par0 <- rep(0, 2*p)
  prior_weights <- c(weights_rand, weights)

  multiroot <- nleqslv::nleqslv(x = par0, # TODO add user-specified parameters to control functions
                                fn = u_theta_beta_dr,
                                method = "Newton", # TODO consider the method Broyden
                                global = "qline", #c("dbldog", "pwldog", cline", "qline", "gline", "hook", "none")
                                xscalm = "fixed", # c("fixed","auto")
                                jacobian = TRUE,
                                control = list(scalex = rep(1, length(par0))), # TODO algorithm did not converge in maxit iterations for cloglog
                                R = R,
                                X = X,
                                y = y,
                                weights = prior_weights,
                                method_selection = method_selection,
                                family_nonprobsvy = family)
  par_sel <- multiroot$x
  if (multiroot$termcd %in% c(2:7, -10)) {
    switch(as.character(multiroot$termcd),
           "2" = warning("Relatively convergent algorithm when fitting selection model by nleqslv, but user must check if function values are acceptably small."),
           "3" = warning("Algorithm did not find suitable point - has stalled cannot find an acceptable new point when fitting selection model by nleqslv."),
           "4" = warning("Iteration limit exceeded when fitting selection model by nleqslv."),
           "5" = warning("ill-conditioned Jacobian when fitting selection model by nleqslv."),
           "6" = warning("Jacobian is singular when fitting selection model by nleqslv."),
           "7" = warning("Jacobian is unusable when fitting selection model by nleqslv."),
           "-10" = warning("user specified Jacobian is incorrect when fitting selection model by nleqslv."))
  }

  theta_hat <- par_sel[1:(p)]
  beta_hat <- par_sel[(p+1):(2*p)]
  names(theta_hat) <- names(beta_hat) <- colnames(X)
  df_residual <- nrow(X) - length(theta_hat)

  ps <- inv_link(theta_hat %*% t(X)) # inv_link(as.vector(X_design %*% as.matrix(theta_hat)))
  ps_der <- dinv_link(theta_hat %*% t(X))
  ps_nons <- ps[loc_nons]
  est_ps_rand <- ps[loc_rand]
  ps_nons_der <- ps_der[loc_nons]
  weights_nons <- 1/ps_nons

  if (!boot) {
    N_nons <- sum(weights * weights_nons)
    # variance-covariance matrix for selection model toFix
    V <- Matrix::Diagonal(n = length(ps), x = ps * (1 - ps))
    vcov_selection <- solve(t(X) %*% V %*% X)
    # vcov_selection <- matrix(0, nrow = nrow(X_design), ncol = ncol(X_design))
    theta_errors <- sqrt(diag(vcov_selection))
  }

  eta <- as.vector(beta_hat %*% t(X))
  y_hat <- family$mu(eta)
  y_rand_pred <- y_hat[loc_rand]
  y_nons_pred <- y_hat[loc_nons]

  if (!boot) {
    sigma <- family$variance(mu = y_hat, y = y[loc_rand])
    residuals <- family$residuals(mu = y_rand_pred, y = y[loc_rand])
  }

  if (!boot) {
    # variance-covariance matrix for outcome model
    # vcov_outcome <- solve(t(X_design) %*% diag(sigma) %*% X_design)
    vcov_outcome <- solve(t(X) %*% (sigma * X))
    beta_errors <- sqrt(diag(vcov_outcome))
  }

  if (!boot) {
    hess <- NULL
    selection <- list(theta_hat = theta_hat,
                      grad = multiroot$fvec[1:(p)], # TODO
                      hess = hess, # TODO
                      ps_nons = ps_nons,
                      variance_covariance = vcov_selection,
                      df_residual = df_residual,
                      log_likelihood = "NULL")

    outcome <- list(beta_hat = beta_hat,
                    grad = multiroot$f.root[(p+1):(2*p)],
                    hess = hess, # TODO
                    variance_covariance = vcov_outcome,
                    df_residual = df_residual,
                    family = list(mu = y_hat,
                                  var = sigma[loc_rand],
                                  residuals = residuals),
                    y_rand_pred = y_rand_pred,
                    y_nons_pred = y_nons_pred,
                    log_likelihood = "NULL")
  } else {
    selection <- list(theta_hat = theta_hat,
                      ps_nons = ps_nons)
    outcome <- list(beta_hat = beta_hat,
                    y_rand_pred = y_rand_pred,
                    y_nons_pred = y_nons_pred)
  }

  list(selection = selection,
       outcome = outcome)
}

##### helpers ########

# joint score equation for theta and beta, used in estimation when variable selections
u_theta_beta_dr <- function(par,
                            R,
                            X,
                            y,
                            weights,
                            method_selection,
                            family_nonprobsvy) {

  method <- get_method(method_selection)

  inv_link <- method$make_link_inv
  inv_link_rev <- method$make_link_inv_rev

  p <- ncol(X)
  theta <- par[1:(p)]
  beta <- par[(p+1):(2*p)]
  eta_pi <- X %*% theta
  ps <- inv_link(eta_pi)
  y[which(is.na(y))] <- 0
  ps <- as.vector(ps)

  eta <- X %*% beta
  mu <- family_nonprobsvy$mu(eta)
  mu_der <- as.vector(family_nonprobsvy$mu_der(eta))
  res <- family_nonprobsvy$residuals(mu = mu, y = y)
  mu_der <- 1

  n <- length(R)
  R_rand <- 1 - R

  utb <- c(apply(X * R/ps * mu_der * weights - X * R_rand * weights * mu_der, 2, sum),
           apply(X * R * weights * as.vector(-inv_link_rev(eta_pi)) * res, 2, sum))/n

  utb
}


u_theta_ipw <- function(par,
                        R,
                        X,
                        y,
                        weights,
                        method_selection) { # TODO

  method <- get_method(method_selection)
  inv_link_rev <- method$make_link_inv_rev
  inv_link <- method$make_link_inv

  p <- ncol(X)
  theta <- par
  X0 <- cbind(1, X)
  eta_pi <- X0 %*% theta
  y[which(is.na(y))] <- 0

  R_rand <- 1- R
  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)
  n <- length(R)
  y_mean <- mean(y[loc_nons])

  #UTB <- apply(X0 * (R * as.vector(inv_link(eta_pi)) - y), 2, sum)/n # TODO
  UTB <- apply(X0 * (R / as.vector(inv_link(eta_pi)) * y - R * y) * as.vector(inv_link_rev(eta_pi)), 2, sum) # TODO

  UTB

}

u_beta_mi <- function(par,
                      R,
                      X,
                      y,
                      weights,
                      family_nonprobsvy) { # TODO

  if(is.character(family_nonprobsvy)) {
    family_nonprobsvy <- paste(family_nonprobsvy, "_nonprobsvy", sep = "")
    family_nonprobsvy <- get(family_nonprobsvy, mode = "function", envir = parent.frame())
    family_nonprobsvy <- family_nonprobsvy()
  }

  p <- ncol(X)
  beta <- par
  eta <- X %*% beta
  mu <- family_nonprobsvy$mu(eta)
  mu_der <- family_nonprobsvy$mu_der(mu)

  n <- length(R)
  R_rand <- 1 - R
  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)
  y_mean <- mean(y[loc_nons])

  UTB <- apply(X * y - X * R_rand * weights * as.vector(mu), 2, sum)
  UTB
}

# TODO Jacobian of the estimating equations for dr method
u_theta_beta_dr_jacob <- function(par,
                                  R,
                                  X,
                                  y,
                                  weights,
                                  method_selection,
                                  family_nonprobsvy){
  method <- get_method(method_selection)

  inv_link <- method$make_link_inv
  inv_link_rev <- method$make_link_inv_rev
  dinv_link_rev <- method$make_link_inv_rev_de

  p <- ncol(X)
  theta <- par[1:(p+1)]
  beta <- par[(p+2):(2*p+2)]
  X0 <- cbind(1, X)
  eta_pi <- X0 %*% theta
  ps <- inv_link(eta_pi)
  y[which(is.na(y))] <- 0
  ps <- as.vector(ps)

  eta <- X0 %*% beta
  mu <- family_nonprobsvy$mu(eta)
  mu_der <- family_nonprobsvy$mu_der(mu)
  mu_der2 <- family_nonprobsvy$mu_der2(mu)
  res <- family_nonprobsvy$residuals(mu = mu, y = y)
  n <- length(R)
  R_rand <- 1 - R

  jac <- c(apply(- X0 * R * weights * as.vector(inv_link_rev(eta_pi)) * mu_der, 2, sum),
           apply(X0 * R/ps * mu_der2 * weights - X0 * R_rand * weights * mu_der2, 2, sum),
           apply(X0 * R * weights * as.vector(dinv_link_rev(eta_pi)) * res * X0, 2, sum),
           apply(X0 * R * weights * as.vector(inv_link_rev(eta_pi)) * mu_der, 2, sum))/n
  jac

}

