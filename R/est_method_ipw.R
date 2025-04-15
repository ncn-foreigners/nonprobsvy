#' ipw maxlik function
#' @noRd
ipw_maxlik <- function(method, X_nons, X_rand, weights, weights_rand, start, control, ...) {

  log_like <- method$make_log_like(
    X_nons,
    X_rand,
    weights,
    weights_rand
  )

  gradient <- method$make_gradient(
    X_nons,
    X_rand,
    weights,
    weights_rand
  )

  hessian <- method$make_hessian(
    X_nons,
    X_rand,
    weights,
    weights_rand
  )

  opt_results <- switch(control$optimizer,
                          "maxLik" = maxLik::maxLik(logLik = log_like,
                                                    grad = gradient,
                                                    hess = hessian,
                                                    method = control$maxlik_method,
                                                    start = start,
                                                    printLevel = control$print_level),
                          "optim" = stats::optim(fn = log_like,
                                                 gr = gradient,
                                                 method = control$optim_method,
                                                 par = start,
                                                 control = list(fnscale = -1,
                                                                trace = control$trace,
                                                                maxit = control$maxit))
                          )


  if (control$optimizer == "maxLik") {
    if (opt_results$code %in% c(3:7, 100)) {
      switch(as.character(opt_results$code),
             "3" = warning("Warning in fitting selection model with the `maxLik` package: probably not converged."),
             "4" = warning("Maxiteration limit reached in fitting selection model by the `maxLik` package."),
             "5" = stop("Infinite value of log_like in fitting selection model by the `maxLik` package, error code 5."),
             "6" = stop("Infinite value of gradient in fitting selection model by the `maxLik` package, error code 6."),
             "7" = stop("Infinite value of hessian in fitting selection model by the `maxLik` package, error code 7."),
             "100" = stop("Error in fitting selection model with the `maxLik` package, error code 100: Bad start."),
      )
    }
    theta <- opt_results$estimate
    grad <- opt_results$gradient
    hess <- opt_results$hessian
    log_likelihood <- log_like(theta)
  }
  if (control$optimizer == "optim") {
    if (as.character(opt_results$convergence) %in% c(1, 10, 51, 52)) {
      switch(as.character(opt_results$convergence),
             "1" = warning("Warning in fitting selection model with the `optim` function:
                           the iteration limit maxit had been reached."),
             "10" = warning("Degeneracy of the Nelder Mead simplex in fitting selection model by the `optim` function."), # TODO -
             "51" = warning("Warning from the L-BFGS-B when fitting by the `optim` function."), # TODO -
             "52" = stop("Indicates an error from the L-BFGS-B method when fitting by the `optim` function.")
      )
    }
    theta <- opt_results$par
    log_likelihood <- log_like(theta)
    grad <- gradient(theta)
    hess <- hessian(theta)
  }


  list(
    log_l = log_likelihood,
    grad = grad,
    hess = hess,
    theta_hat = theta
  )
}

#' gee functions
#' score equation for theta, used in variable selection
#' @noRd
u_theta <- function(R,
                    X,
                    weights,
                    method_selection,
                    gee_h_fun,
                    pop_totals = NULL) {

  method <- switch(method_selection,
                   "logit" = method_ps("logit"),
                   "probit" = method_ps("probit"),
                   "cloglog" = method_ps("cloglog"))

  inv_link <- method$make_link_inv
  function(par) {
    # loc_nons = which(R == 1)
    # loc_rand = which(R == 0)
    theta <- as.matrix(par)
    n <- length(R)
    X0 <- as.matrix(X)
    eta_pi <- X0 %*% theta
    ps <- inv_link(eta_pi)
    R_rand <- 1 - R
    ps <- as.vector(ps)
    N_nons <- sum(1 / ps)

    # "1" = t(X0[loc_nons,]) %*% (1/ps[loc_nons]) - t(X0[loc_rand,]) %*% weights[loc_rand],
    # "2" = c(apply(X0 * R * weights - X0 * R_rand * ps * weights, 2, sum))
    if (is.null(pop_totals)) {
      eq <- switch(gee_h_fun,
                   "1" = c(apply(X0 * R / ps * weights - X0 * R_rand * weights, 2, sum)), # consider division by N_nons
                   "2" = c(apply(X0 * R * weights - X0 * R_rand * ps * weights, 2, sum))
      )
    } else {
      eq <- c(apply(X0 * R / ps * weights, 2, sum)) - pop_totals
    }
    eq
  }
}

# derivative of score equation for theta, used in variable selection
u_theta_der <- function(R,
                        X,
                        weights,
                        method_selection,
                        gee_h_fun,
                        pop_totals = NULL) {

  method <- switch(method_selection,
                   "logit" = method_ps("logit"),
                   "probit" = method_ps("probit"),
                   "cloglog" = method_ps("cloglog"))

  inv_link <- method$make_link_inv
  dinv_link <- method$make_link_inv_der
  inv_link_rev <- method$make_link_inv_rev

  function(par) {
    # loc_nons = which(R == 1)
    # loc_rand = which(R == 0)
    theta <- as.matrix(par)
    X0 <- as.matrix(X)
    p <- ncol(X0)
    eta <- as.numeric(X0 %*% theta)
    ps <- inv_link(eta)
    ps <- as.vector(ps)
    N_nons <- sum(1 / ps)
    R_rand <- 1 - R

    # "1" = t(X0[loc_nons, ]) %*% weights[loc_nons] %*% t(inv_link_rev(eta)[loc_nons]) %*% X0[loc_nons, ],
    # "2" =
    if (!is.null(pop_totals)) {
      mxDer <- t(R * X0 * weights * inv_link_rev(eta)) %*% X0
    } else {
      mxDer <- switch(gee_h_fun,
                      "1" = t(R * X0 * weights * inv_link_rev(eta)) %*% X0, # TODO bug here when solve for some data - probably because of inv_link_rev
                      "2" = -t(R_rand * X0 * weights * dinv_link(eta)) %*% X0
      )
    }
    as.matrix(mxDer, nrow = p) # consider division by N_nons
  }
}

theta_h_estimation <- function(R,
                               X,
                               weights_rand,
                               weights,
                               gee_h_fun,
                               method_selection,
                               maxit,
                               nleqslv_method,
                               nleqslv_global,
                               nleqslv_xscalm,
                               start = NULL,
                               pop_totals = NULL) {

  if (is.null(start)) start <- rep(0, ncol(X))

  p <- ncol(X)
  u_theta <- u_theta(
    R = R,
    X = X,
    weights = c(weights_rand, weights),
    gee_h_fun = gee_h_fun,
    method_selection = method_selection,
    pop_totals = pop_totals
  )

  u_theta_der <- u_theta_der(
    R = R,
    X = X,
    weights = c(weights_rand, weights),
    gee_h_fun = gee_h_fun,
    method_selection = method_selection,
    pop_totals = pop_totals
  )

  root <- nleqslv::nleqslv(
    x = start,
    fn = u_theta,
    method = nleqslv_method,
    global = nleqslv_global,
    xscalm = nleqslv_xscalm,
    jacobian = TRUE,
    jac = if (method_selection == "cloglog") NULL else u_theta_der,
    control = list(maxit = maxit)
  )

  theta_root <- root$x

  if (root$termcd %in% c(2:7, -10)) {
    switch(as.character(root$termcd),
           "2" = warning("Relatively convergent algorithm when fitting selection model by nleqslv,
                         but user must check if function values are acceptably small."),
           "3" = warning("Algorithm did not find suitable point - has stalled cannot find an acceptable
                         new point when fitting selection model by nleqslv."),
           "4" = warning("Iteration limit exceeded when fitting selection model by nleqslv."),
           "5" = warning("Ill-conditioned Jacobian when fitting selection model by nleqslv."),
           "6" = warning("Jacobian is singular when fitting selection model by nleqslv."),
           "7" = warning("Jacobian is unusable when fitting selection model by nleqslv."),
           "-10" = warning("User specified Jacobian is incorrect when fitting selection model by nleqslv.")
    )
  }

  theta_h <- as.vector(theta_root)
  grad <- u_theta(theta_h)
  hess <- if (method_selection == "cloglog") root$jac else u_theta_der(theta_h)

  list(
    theta_h = theta_h,
    hess = hess,
    grad = grad
  )
}




#' estimation method for the IPW estimator
#' @noRd
est_method_ipw <- function(est_method = c("gee", "mle"), ...) {

  est_method <- match.arg(est_method)

  ### generalized estimation equation
  gee <- function(...) {

    estimation_model <- function(model, method_selection, ...) {
      method <- model$method
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
      eta <- c(model$eta_rand, model$eta_nons)
      residuals <- model$residuals
      variance <- as.vector(model$variance)

      list(
        theta_hat = theta_hat,
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
        log_likelihood = NA,
        eta = eta,
        aic = NA,
        variance = variance,
        residuals = residuals,
        method = method
      )
    }

    make_t_comp <- function(X, ps, psd, b, y_rand, y_nons, gee_h_fun, N,
                            method_selection, weights) {
      # common calculations
      mean_nons <- sum(weights * y_nons) / N
      Xb <- X %*% t(as.matrix(b))

      # choose formula based on gee_h_fun
      if (gee_h_fun == 1) {
        t_comp <- Xb + y_rand - mean_nons
      } else if (gee_h_fun == 2) {
        t_comp <- as.vector(ps) * Xb + y_rand - mean_nons
      }
      t_comp
    }

    make_var_nonprob <- function(ps, psd, y, y_pred, h_n, X, b, N, gee_h_fun,
                                 method_selection, weights, pop_totals) {
      # force gee_h_fun=1 if pop_totals available
      if (!is.null(pop_totals)) gee_h_fun <- 1

      # common terms
      resid <- weights * (y - y_pred - h_n)
      model_adj <- b %*% t(X)

      # variance calculation based on h
      if (gee_h_fun == 2) {
        var_nonprob <- 1 / N^2 * sum((1 - ps) * (resid / ps - model_adj)^2)
      } else if (gee_h_fun == 1) {
        var_nonprob <- 1 / N^2 * sum((1 - ps) * ((resid - model_adj) / ps)^2)
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
                                gee_h_fun = gee_h_fun,
                                est_method,
                                maxit,
                                control_selection,
                                start,
                                verbose,
                                varcov = FALSE,
                                ...) {

      method <- switch(method_selection,
                       "logit" = method_ps("logit"),
                       "probit" = method_ps("probit"),
                       "cloglog" = method_ps("cloglog"))

      inv_link <- method$make_link_inv

      if (is.null(start)) {

        start <- switch(control_selection$start_type,
                        "zero" = rep(0, ncol(X)),
                        "mle" = {
                          ipw_maxlik(
                            method,
                            X_nons = X_nons,
                            X_rand = X_rand,
                            weights = weights,
                            weights_rand = weights_rand,
                            start = rep(0, ncol(X)),
                            control = control_selection
                          )$theta_hat
                        })
        }

      h_object <- theta_h_estimation(
        R = R,
        X = X,
        weights_rand = weights_rand,
        weights = weights,
        gee_h_fun = gee_h_fun,
        method_selection = method_selection,
        maxit = maxit,
        nleqslv_method = control_selection$nleqslv_method,
        nleqslv_global = control_selection$nleqslv_global,
        nleqslv_xscalm = control_selection$nleqslv_xscalm,
        start = start
      )

      theta_hat <- h_object$theta_h
      hess <- h_object$hess
      grad <- h_object$grad
      eta_nons <- theta_hat %*% t(as.matrix(X_nons))
      eta_rand <- theta_hat %*% t(as.matrix(X_rand))
      ps_nons <- inv_link(eta_nons)
      est_ps_rand <- inv_link(eta_rand)
      variance_covariance <- try(chol2inv(chol(-hess)), silent = TRUE)

      if (inherits(variance_covariance, "try-error")) {
        if (verbose) message("chol2inv(chol()) failed, using ginv() instead.")
        variance_covariance <- MASS::ginv(-hess)
      }

      resids <- R - c(est_ps_rand, ps_nons)

      df_reduced <- nrow(X) - length(theta_hat)
      variance <- as.vector((t(resids) %*% resids) / df_reduced)

      if (method_selection == "probit") { # for probit model, propensity score derivative is required
        dinv_link <- method$make_link_inv_der
        ps_nons_der <- dinv_link(theta_hat %*% t(as.matrix(X_nons)))
        est_ps_rand_der <- dinv_link(theta_hat %*% t(as.matrix(X_rand)))
      }

      list(
        theta_hat = theta_hat,
        hess = hess,
        grad = grad,
        ps_nons = ps_nons,
        est_ps_rand = est_ps_rand,
        ps_nons_der = ifelse(method_selection == "probit", ps_nons_der, NA),
        est_ps_rand_der = ifelse(method_selection == "probit", est_ps_rand_der, NA),
        variance_covariance = variance_covariance,
        var_cov1 = ifelse(varcov, method$variance_covariance1, "No variance-covariance matrix"),
        var_cov2 = ifelse(varcov, method$variance_covariance2, "No variance-covariance matrix"),
        df_residual = df_reduced,
        eta_nons = eta_nons,
        eta_rand = eta_rand,
        residuals = resids,
        method = method
      )
    }

    return(
      list(
        estimation_model = estimation_model,
        make_t_comp = make_t_comp,
        make_var_nonprob = make_var_nonprob,
        model_selection = model_selection)
      )
  }

  ### maximum likelihood function
  mle <- function(...) {

    estimation_model <- function(model, method_selection, ...) {

      method <- model$method
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
      variance_covariance <- MASS::ginv(-hess) # MASS::ginv # variance-covariance matrix of estimated parameters
      eta <- c(model$eta_rand, model$eta_nons)
      aic <- 2 * (length(theta_hat) - log_likelihood)
      residuals <- model$residuals
      variance <- as.vector(model$variance)
      deviance <- model$deviance
      deviance_null <- model$deviance_null

      list(
        theta_hat = theta_hat,
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
        df_residual = df_residual,
        eta = eta,
        aic = aic,
        variance = variance,
        residuals = residuals,
        method = method
      )
    }

    make_t_comp <- function(X, ps, psd, b, y_rand, y_nons, gee_h_fun, N, method_selection, weights) {

      method <- switch(method_selection,
                       "logit" = method_ps("logit"),
                       "probit" = method_ps("probit"),
                       "cloglog" = method_ps("cloglog"))

      t_comp <- method$t_vec(
        X = X,
        ps = ps,
        psd = psd,
        b = b,
        y_rand = y_rand,
        y_nons = y_nons,
        N = N,
        weights = weights
      )
      t_comp
    }

    make_var_nonprob <- function(ps, psd, y, y_pred, h_n, X, b, N, gee_h_fun, method_selection,
                                 weights = weights, pop_totals = NULL) {

      method <- switch(method_selection,
                       "logit" = method_ps("logit"),
                       "probit" = method_ps("probit"),
                       "cloglog" = method_ps("cloglog"))

      var_nonprob <- method$var_nonprob(
        ps = ps,
        psd = psd,
        y = y,
        y_pred = y_pred,
        h_n = h_n,
        X = X,
        b = b,
        N = N,
        weights = weights
      )
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
                                gee_h_fun = gee_h_fun,
                                est_method,
                                maxit,
                                control_selection,
                                start,
                                verbose = FALSE,
                                varcov = FALSE,
                                ...) {

      method <- switch(method_selection,
                       "logit" = method_ps("logit"),
                       "probit" = method_ps("probit"),
                       "cloglog" = method_ps("cloglog"))

      inv_link <- method$make_link_inv
      dinv_link <- method$make_link_inv_der

      # initial values for propensity score estimation
      if (is.null(start)) start <- rep(0, ncol(X))

      df_reduced <- nrow(X) - length(start)

      maxLik_nons_obj <- ipw_maxlik(
        method,
        X_nons = X_nons,
        X_rand = X_rand,
        weights = weights,
        weights_rand = weights_rand,
        start = start,
        control = control_selection
      )

      theta <- maxLik_nons_obj$theta_hat
      eta_nons <- theta %*% t(X_nons)
      eta_rand <- theta %*% t(X_rand)

      ps_nons <- inv_link(eta_nons)
      est_ps_rand <- inv_link(eta_rand)

      ps_nons_der <- dinv_link(eta_nons)
      est_ps_rand_der <- dinv_link(eta_rand)

      resids <- R - c(est_ps_rand, ps_nons)

      variance <- (t(resids) %*% resids) / df_reduced

      list(
        maxLik_nons_obj = maxLik_nons_obj,
        theta = theta,
        ps = ps_nons,
        ps_der = ps_nons_der,
        ps_rand = est_ps_rand,
        ps_rand_der = est_ps_rand_der,
        var_cov1 = ifelse(varcov, method$variance_covariance1, "No variance-covariance matrix"),
        var_cov2 = ifelse(varcov, method$variance_covariance2, "No variance-covariance matrix"),
        df_residual = df_reduced,
        eta_nons = eta_nons,
        eta_rand = eta_rand,
        residuals = resids,
        variance = variance,
        method = method
      )
    }
    return(
      list(
        estimation_model = estimation_model,
        make_t_comp = make_t_comp,
        make_var_nonprob = make_var_nonprob,
        model_selection = model_selection)
      )
  }



  est_funs <- switch(est_method,
                     "gee" = gee(),
                     "mle" = mle())

  return(est_funs)
}
