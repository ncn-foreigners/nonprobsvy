# Object with output parameters for estimation by Generalized Estimating Equations for propensity scores
gee <- function(...) {
  estimation_model <- function(model, method_selection) {
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

  make_t_comp <- function(X, ps, psd, b, y_rand, y_nons, h, N, method_selection, weights) {
    # common calculations
    mean_nons <- sum(weights * y_nons) / N
    Xb <- X %*% t(as.matrix(b))

    # choose formula based on h
    if (h == 1) {
      t_comp <- Xb + y_rand - mean_nons
    } else if (h == 2) {
      t_comp <- as.vector(ps) * Xb + y_rand - mean_nons
    }
    t_comp
  }

  make_var_nonprob <- function(ps, psd, y, y_pred, h_n, X, b, N, h, method_selection, weights, pop_totals) {
    # force h=1 if pop_totals available
    if (!is.null(pop_totals)) h <- 1

    # common terms
    resid <- weights * (y - y_pred - h_n)
    model_adj <- b %*% t(X)

    # variance calculation based on h
    if (h == 2) {
      var_nonprob <- 1 / N^2 * sum((1 - ps) * (resid / ps - model_adj)^2)
    } else if (h == 1) {
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
                              h = h,
                              est_method,
                              maxit,
                              control_selection,
                              start,
                              verbose,
                              varcov = FALSE,
                              ...) {
    method_selection_function <- paste(method_selection, "_model_nonprobsvy", sep = "")
    method <- get_method(method = method_selection_function)
    inv_link <- method$make_link_inv

    if (is.null(start)) {
      if (control_selection$start_type == "glm") {
        start_to_gee <- start_fit(
          X = X, # <--- does not work with pop_totals
          R = R,
          weights = weights,
          weights_rand = weights_rand,
          method_selection = method_selection
        )
        start <- method$make_max_lik(
          X_nons = X_nons,
          X_rand = X_rand,
          weights = weights,
          weights_rand = weights_rand,
          start = start_to_gee,
          control = control_selection
        )$theta_hat
        ####
      } else if (control_selection$start_type == "naive") {
        start_h <- suppressWarnings(theta_h_estimation(
          R = R,
          X = X[, 1, drop = FALSE],
          weights_rand = weights_rand,
          weights = weights,
          h = h,
          method_selection = method_selection,
          maxit = maxit,
          nleqslv_method = control_selection$nleqslv_method,
          nleqslv_global = control_selection$nleqslv_global,
          nleqslv_xscalm = control_selection$nleqslv_xscalm,
          start = 0
        )$theta_h)
        start <- c(start_h, rep(0, ncol(X) - 1))
      } else if (control_selection$start_type == "zero") {
        start <- rep(0, ncol(X))
      }
    }

    h_object <- theta_h_estimation(
      R = R,
      X = X,
      weights_rand = weights_rand,
      weights = weights,
      h = h,
      method_selection = method_selection,
      maxit = maxit,
      nleqslv_method = control_selection$nleqslv_method,
      nleqslv_global = control_selection$nleqslv_method,
      nleqslv_xscalm = control_selection$nleqslv_method,
      start = start
    )
    theta_hat <- h_object$theta_h
    hess <- h_object$hess
    grad <- h_object$grad
    eta_nons <- theta_hat %*% t(as.matrix(X_nons))
    eta_rand <- theta_hat %*% t(as.matrix(X_rand))
    ps_nons <- inv_link(eta_nons)
    est_ps_rand <- inv_link(eta_rand)
    variance_covariance <- try(solve(-hess), silent = TRUE)
    if (inherits(variance_covariance, "try-error")) {
      if (verbose) message("solve() failed, using ginv() instead.")
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

  structure(
    list(
      estimation_model = estimation_model,
      make_t_comp = make_t_comp,
      make_var_nonprob = make_var_nonprob,
      model_selection = model_selection
    ),
    class = "method" ## this can be removed
  )
}
