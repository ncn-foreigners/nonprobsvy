# bias correction
mm <- function(X, y, weights, weights_rand, R, n_nons, n_rand, method_selection, family, start_selection, start_outcome, boot = FALSE) {
  method_selection_function <- paste(method_selection, "_model_nonprobsvy", sep = "")
  method <- get_method(method_selection_function)
  inv_link <- method$make_link_inv
  dinv_link <- method$make_link_inv_der

  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)

  start <- c(start_outcome, start_selection) # TODO consider add info/error for end-user if one of starts provided only

  p <- ncol(X)
  if (is.null(start)) { # TODO add default start
    par0 <- rep(0, 2 * p)
  } else {
    par0 <- start
  }
  prior_weights <- c(weights_rand, weights)

  # MI - bias correction #########
  # multiroot <- nleqslv::nleqslv( # TODO to fix "Jacobian is completely unusable (all zero entries?)"
  #   x = rep(0, p), # TODO add user-specified parameters to control functions
  #   fn = u_beta_mi,# TODO algorithm did not converge in maxit iterations for cloglog
  #   R = R,
  #   X = X,
  #   y = y,
  #   weights = prior_weights,
  #   family_nonprobsvy = family
  # )
  # print(multiroot$x)
  ##########

  # IPW - bias correction #########
  # multiroot <- nleqslv::nleqslv(
  #   x = rep(0, p), # TODO add user-specified parameters to control functions
  #   fn = u_theta_ipw,
  #   method = "Newton", # TODO consider the method Broyden
  #   global = "qline", # c("dbldog", "pwldog", cline", "qline", "gline", "hook", "none")
  #   xscalm = "fixed", # c("fixed","auto")
  #   jacobian = TRUE,
  #   control = list(scalex = rep(1, length(rep(0, p)))), # TODO algorithm did not converge in maxit iterations for cloglog
  #   R = R,
  #   X = X,
  #   y = y,
  #   weights = weights,
  #   method_selection = method_selection
  # )
  # print(multiroot$x)
  ##########

  ######### BB
  # multiroot <- nleqslv::nleqslv(
  #   par = par0, # TODO add user-specified parameters to control functions
  #   fn = u_theta_beta_dr,
  #   R = R,
  #   X = X,
  #   y = y,
  #   weights = prior_weights,
  #   method_selection = method_selection,
  #   family_nonprobsvy = family
  # )
  # par_sel <- multiroot$par
  ######### NLESQLV
  multiroot <- nleqslv::nleqslv(
    x = par0, # TODO add user-specified parameters to control functions
    fn = u_theta_beta_dr,
    method = "Newton", # TODO consider the method Broyden
    global = "dbldog", # c("dbldog", "pwldog", cline", "qline", "gline", "hook", "none")
    xscalm = "auto", # c("fixed","auto")
    jacobian = TRUE,
    control = list(scalex = rep(1, length(par0))), # TODO algorithm did not converge in maxit iterations for cloglog
    R = R,
    X = X,
    y = y,
    weights = prior_weights,
    method_selection = method_selection,
    family_nonprobsvy = family
  )
  par_sel <- multiroot$x
  if (multiroot$termcd %in% c(2:7, -10)) {
    switch(as.character(multiroot$termcd),
      "2" = warning("Relatively convergent algorithm when fitting selection model by nleqslv, but user must check if function values are acceptably small."),
      "3" = warning("Algorithm did not find suitable point - has stalled cannot find an acceptable new point when fitting selection model by nleqslv."),
      "4" = warning("Iteration limit exceeded when fitting selection model by nleqslv."),
      "5" = warning("ill-conditioned Jacobian when fitting selection model by nleqslv."),
      "6" = warning("Jacobian is singular when fitting selection model by nleqslv."),
      "7" = warning("Jacobian is unusable when fitting selection model by nleqslv."),
      "-10" = warning("user specified Jacobian is incorrect when fitting selection model by nleqslv.")
    )
  }

  theta_hat <- par_sel[1:(p)]
  beta_hat <- par_sel[(p + 1):(2 * p)]
  names(theta_hat) <- names(beta_hat) <- colnames(X)
  df_residual <- nrow(X) - length(theta_hat)

  # selection parameters
  ps <- inv_link(theta_hat %*% t(X)) # inv_link(as.vector(X_design %*% as.matrix(theta_hat)))
  eta_sel <- theta_hat %*% t(X)
  ps_der <- dinv_link(eta_sel)
  ps_nons <- ps[loc_nons]
  est_ps_rand <- ps[loc_rand]
  ps_nons_der <- ps_der[loc_nons]
  weights_nons <- 1 / ps_nons
  resids <- R - c(est_ps_rand, ps_nons)
  variance <- as.vector((t(resids) %*% resids) / df_residual)

  if (!boot) {
    N_nons <- sum(weights * weights_nons)
    # variance-covariance matrix for selection model toFix
    V <- Matrix::Diagonal(n = length(ps), x = ps * (1 - ps))
    vcov_selection <- solve(t(X) %*% V %*% X)
    # vcov_selection <- matrix(0, nrow = nrow(X_design), ncol = ncol(X_design))
    theta_errors <- sqrt(diag(vcov_selection))
  }

  eta_out <- as.vector(beta_hat %*% t(X))
  y_hat <- family$linkinv(eta_out)
  y_rand_pred <- y_hat[loc_rand]
  y_nons_pred <- y_hat[loc_nons]

  if (!boot) {
    # sigma_nons <- family$variance(mu = y_nons_pred, y = y[loc_nons])
    # sigma_rand <- family$variance(mu = y_rand_pred, y = y[loc_rand])
    sigma_nons <- family$variance(mu = y_nons_pred)
    sigma_rand <- family$variance(mu = y_rand_pred)
    residuals <- family$residuals(mu = y_nons_pred, y = y[loc_nons])
  }

  if (!boot) {
    # variance-covariance matrix for outcome model
    # vcov_outcome <- solve(t(X_design) %*% diag(sigma) %*% X_design)
    vcov_outcome <- solve(t(X[loc_nons, ]) %*% (sigma_nons * X[loc_nons, ]))
    beta_errors <- sqrt(diag(vcov_outcome))
  }

  # grad = multiroot$f.root[(p+1):(2*p)]
  if (!boot) {
    hess <- NA
    selection <- list(
      theta_hat = theta_hat, # TODO list as close as possible to SelecttionList
      grad = multiroot$fvec[1:(p)],
      hess = hess, # TODO
      ps_nons = ps_nons,
      est_ps_rand = est_ps_rand,
      variance_covariance = vcov_selection,
      df_residual = df_residual,
      log_likelihood = NA,
      eta = eta_sel,
      aic = NA,
      residuals = resids,
      variance = variance,
      method = method
    )

    outcome <- list(
      coefficients = beta_hat, # TODO list as close as possible to glm
      std_err = beta_errors,
      variance_covariance = vcov_outcome,
      df_residual = df_residual,
      family = list(
        mu = y_nons_pred,
        variance = sigma_nons,
        family = family$family
      ),
      residuals = residuals,
      fitted.values = y_nons_pred,
      sigma_rand = sigma_rand,
      y_rand_pred = y_rand_pred,
      y_nons_pred = y_nons_pred,
      linear.predictors = eta_out[loc_nons],
      X = X[loc_nons, ]
    )
  } else {
    selection <- list(
      coefficients = theta_hat, # TODO list as close as possible to SelecttionList
      ps_nons = ps_nons
    )
    outcome <- list(
      coefficients = beta_hat,
      y_rand_pred = y_rand_pred, # TODO list as close as possible to SelecttionList
      y_nons_pred = y_nons_pred
    )
  }

  list(
    selection = selection,
    outcome = outcome
  )
}
