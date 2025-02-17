# Object with output parameters for Maximum likelihood Estimation for propensity scores
## this should be moved to MLE() function - that estimate parameters

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

  if (control$optimizer == "maxLik") {
    ########### maxLik ##########
    maxLik_an <- maxLik::maxLik(
      logLik = log_like,
      grad = gradient,
      hess = hessian,
      method = "BFGS",
      start = start,
      printLevel = control$print_level
    )

    if (maxLik_an$code %in% c(3:7, 100)) {
      switch(as.character(maxLik_an$code),
             "3" = warning("Warning in fitting selection model with the `maxLik` package: probably not converged."),
             "4" = warning("Maxiteration limit reached in fitting selection model by the `maxLik` package."),
             "5" = stop("Infinite value of log_like in fitting selection model by the `maxLik` package, error code 5."),
             "6" = stop("Infinite value of gradient in fitting selection model by the `maxLik` package, error code 6."),
             "7" = stop("Infinite value of hessian in fitting selection model by the `maxLik` package, error code 7."),
             "100" = stop("Error in fitting selection model with the `maxLik` package, error code 100: Bad start."),
      )
    }

    theta <- maxLik_an$estimate
    grad <- maxLik_an$gradient
    hess <- maxLik_an$hessian
    log_likelihood <- log_like(theta)
  } else if (control$optimizer == "optim") { # TODO add optimParallel for high-dimensional data
    ########### optim ##########
    maxLik_an <- stats::optim(
      fn = log_like,
      gr = gradient,
      method = control$optim_method,
      par = start,
      control = list(
        fnscale = -1,
        trace = control$trace,
        maxit = control$maxit
      )
    )
    if (maxLik_an$convergence %in% c(1, 10, 51, 52)) {
      switch(as.character(maxLik_an$convergence),
             "1" = warning("Warning in fitting selection model with the `optim` function: the iteration limit maxit had been reached."),
             "10" = warning("Degeneracy of the Nelder Mead simplex in fitting selection model by the `optim` function."), # TODO -
             "51" = warning("Warning from the L-BFGS-B when fitting by the `optim` function."), # TODO -
             "52" = stop("Indicates an error from the L-BFGS-B method when fitting by the `optim` function.")
      )
    }
    theta <- maxLik_an$par
    log_likelihood <- log_like(theta)
    grad <- gradient(theta)
    hess <- hessian(theta)
  } else {
    stop("Provide valid optimizer (`optim` or `maxLik`).")
  }

  list(
    log_l = log_likelihood,
    grad = grad,
    hess = hess,
    theta_hat = theta
  )
}


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
    variance_covariance <- solve(-hess) # MASS::ginv # variance-covariance matrix of estimated parameters
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

  make_t_comp <- function(X, ps, psd, b, y_rand, y_nons, gee_h_fun, N, method_selection, weights, weights_sum) {
    method <- get_method(method_selection)
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
                               weights = weights, weights_sum, pop_totals = NULL) {
    method <- get_method(method_selection)
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

    method_selection_function <- paste(method_selection, "_model_nonprobsvy", sep = "")

    method <- get_method(method = method_selection_function)

    inv_link <- method$make_link_inv
    dinv_link <- method$make_link_inv_der

    # initial values for propensity score estimation
    if (is.null(start)) {
      if (control_selection$start_type == "glm") {
        start <- start_fit(
          X = X,
          R = R,
          weights = weights,
          weights_rand = weights_rand,
          method_selection = method_selection
        )
      } else if (control_selection$start_type == "naive") {

        intercept_start <- suppressWarnings(ipw_maxlik(
          method,
          X_nons = X_nons[, 1, drop = FALSE],
          X_rand = X_rand[, 1, drop = FALSE],
          weights = weights,
          weights_rand = weights_rand,
          start = 0,
          control = control_selection
        )$theta_hat)
        start <- c(intercept_start, rep(0, ncol(X_nons) - 1))
      } else if (control_selection$start_type == "zero") {
        start <- rep(0, ncol(X))
      }
    }

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
  structure(
    list(
      estimation_model = estimation_model,
      make_t_comp = make_t_comp,
      make_var_nonprob = make_var_nonprob,
      model_selection = model_selection
    ),
    class = "method"
  )
}
