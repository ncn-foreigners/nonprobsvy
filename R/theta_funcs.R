# score equation for theta, used in variable selection
u_theta <- function(R,
                    X,
                    weights,
                    method_selection,
                    h,
                    N = NULL,
                    pop_totals = NULL,
                    pop_size = NULL) {
  method_selection <- paste(method_selection, "_model_nonprobsvy", sep = "")
  method <- get_method(method_selection)
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
    weights_sum <- sum(weights)

    # "1" = t(X0[loc_nons,]) %*% (1/ps[loc_nons]) - t(X0[loc_rand,]) %*% weights[loc_rand],
    # "2" = c(apply(X0 * R * weights - X0 * R_rand * ps * weights, 2, sum))
    if (is.null(pop_totals)) {
      eq <- switch(h,
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
                        h,
                        N = NULL,
                        pop_totals = NULL) {
  method_selection <- paste(method_selection, "_model_nonprobsvy", sep = "")
  method <- get_method(method_selection)
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
    weights_sum <- sum(weights)

    # "1" = t(X0[loc_nons, ]) %*% weights[loc_nons] %*% t(inv_link_rev(eta)[loc_nons]) %*% X0[loc_nons, ],
    # "2" =
    if (!is.null(pop_totals)) {
      mxDer <- t(R * X0 * weights * inv_link_rev(eta)) %*% X0
    } else {
      mxDer <- switch(h,
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
                               h,
                               method_selection,
                               maxit,
                               nleqslv_method,
                               nleqslv_global,
                               nleqslv_xscalm,
                               start = NULL,
                               pop_totals = NULL,
                               pop_means = NULL) {
  p <- ncol(X)
  # if (is.null(pop_totals) & is.null(pop_means)) {
  #   if (is.null(start)) {
  #     start0 <- start_fit(X = X, # <--- does not work with pop_totals
  #                         R = R,
  #                         weights = weights,
  #                         weights_rand = weights_rand,
  #                         method_selection = method_selection)
  #   } else {
  #     start0 <- start
  #   }
  # } else { # TODO customize start point for fitting with population totals
  #   # start0 <- rep(.8, ncol(X))
  #   # X_pop <- rbind(X, pop_totals)
  #   # weights_randd <- 1
  #   if (is.null(start)) {
  #     start0 <- start_fit(X = X, # <--- does not work with pop_totals
  #                         R = R,
  #                         weights = weights,
  #                         weights_rand = weights_rand,
  #                         method_selection = method_selection)
  #   } else {
  #     start0 <- start
  #   }
  # }
  u_theta <- u_theta(
    R = R,
    X = X,
    weights = c(weights_rand, weights),
    h = h,
    method_selection = method_selection,
    pop_totals = pop_totals
  )

  u_theta_der <- u_theta_der(
    R = R,
    X = X,
    weights = c(weights_rand, weights),
    h = h,
    method_selection = method_selection,
    pop_totals = pop_totals
  )
  # print(start)
  # ######### BB
  # if (method_selection == "cloglog") {
  #   root <- BB::dfsane(
  #     par = start,
  #     fn = u_theta,
  #   )
  # } else {
  #   root <- BB::dfsane(
  #     par = start,
  #     fn = u_theta
  #     # control = list(sigma = 0.1, trace = 1)
  #   )
  # }
  # theta_root <- root$par
  # print(theta_root)
  ######### NLESQLV
  if (method_selection == "cloglog") {
    root <- nleqslv::nleqslv(
      x = start,
      fn = u_theta,
      method = "Newton", # TODO consider the methods
      global = "dbldog", # qline",
      xscalm = "auto",
      jacobian = TRUE,
      control = list(maxit = maxit)
    )
  } else {
    root <- nleqslv::nleqslv(
      x = start,
      fn = u_theta,
      method = "Newton", # TODO consider the methods
      global = "dbldog", # qline",
      xscalm = "auto",
      jacobian = TRUE,
      jac = u_theta_der,
      control = list(maxit = maxit)
    )
  }
  theta_root <- root$x
  # stop("123")
  if (root$termcd %in% c(2:7, -10)) {
    switch(as.character(root$termcd),
      "2" = warning("Relatively convergent algorithm when fitting selection model by nleqslv, but user must check if function values are acceptably small."),
      "3" = warning("Algorithm did not find suitable point - has stalled cannot find an acceptable new point when fitting selection model by nleqslv."),
      "4" = warning("Iteration limit exceeded when fitting selection model by nleqslv."),
      "5" = warning("ill-conditioned Jacobian when fitting selection model by nleqslv."),
      "6" = warning("Jacobian is singular when fitting selection model by nleqslv."),
      "7" = warning("Jacobian is unusable when fitting selection model by nleqslv."),
      "-10" = warning("user specified Jacobian is incorrect when fitting selection model by nleqslv.")
    )
  }
  theta_h <- as.vector(theta_root)
  grad <- u_theta(theta_h)
  if (method_selection == "cloglog") {
    hess <- root$jac
  } else {
    hess <- u_theta_der(theta_h) # TODO compare with root$jac
  }

  list(
    theta_h = theta_h,
    hess = hess,
    grad = grad
  )
}



# joint score equation for theta and beta, used in estimation when variable selections
u_theta_beta_dr <- function(par,
                            R,
                            X,
                            y,
                            weights,
                            method_selection,
                            family_nonprobsvy) {
  method_selection <- paste(method_selection, "_model_nonprobsvy", sep = "")
  method <- get_method(method_selection)

  inv_link <- method$make_link_inv
  inv_link_rev <- method$make_link_inv_rev

  p <- ncol(X)
  theta <- par[1:(p)]
  beta <- par[(p + 1):(2 * p)]
  eta_pi <- X %*% theta
  ps <- inv_link(eta_pi)
  y[which(is.na(y))] <- 0
  ps <- as.vector(ps)

  eta <- X %*% beta
  mu <- family_nonprobsvy$linkinv(eta)
  mu_der <- as.vector(family_nonprobsvy$mu.eta(eta))
  res <- family_nonprobsvy$residuals(mu = mu, y = y)
  mu_der <- 1

  n <- length(R)
  R_rand <- 1 - R

  utb <- c(
    apply(X * R / ps * mu_der * weights - X * R_rand * weights * mu_der, 2, sum),
    apply(X * R * weights * as.vector(-inv_link_rev(eta_pi)) * res, 2, sum)
  ) / n

  utb
}


u_theta_ipw <- function(par,
                        R,
                        X,
                        y,
                        weights,
                        method_selection) { # TODO

  method_selection <- paste(method_selection, "_model_nonprobsvy", sep = "")
  method <- get_method(method_selection)
  inv_link_rev <- method$make_link_inv_rev
  inv_link <- method$make_link_inv

  p <- ncol(X)
  theta <- par
  eta_pi <- X %*% theta
  y[which(is.na(y))] <- 0

  R_rand <- 1 - R
  loc_nons <- which(R == 1)
  loc_rand <- which(R == 0)
  n <- length(R)
  y_mean <- mean(y[loc_nons])

  # UTB <- apply(X0 * (R * as.vector(inv_link(eta_pi)) - y), 2, sum)/n # TODO
  UTB <- apply(X * (R / as.vector(inv_link(eta_pi)) * y - mean(y)) * as.vector(inv_link_rev(eta_pi)), 2, sum) # TODO

  UTB
}

# TODO Jacobian of the estimating equations for dr method
u_theta_beta_dr_jacob <- function(par,
                                  R,
                                  X,
                                  y,
                                  weights,
                                  method_selection,
                                  family_nonprobsvy) {
  method_selection <- paste(method_selection, "_model_nonprobsvy", sep = "")
  method <- get_method(method_selection)

  inv_link <- method$make_link_inv
  inv_link_rev <- method$make_link_inv_rev
  dinv_link_rev <- method$make_link_inv_rev_de

  p <- ncol(X)
  theta <- par[1:(p + 1)]
  beta <- par[(p + 2):(2 * p + 2)]
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

  jac <- c(
    apply(-X0 * R * weights * as.vector(inv_link_rev(eta_pi)) * mu_der, 2, sum),
    apply(X0 * R / ps * mu_der2 * weights - X0 * R_rand * weights * mu_der2, 2, sum),
    apply(X0 * R * weights * as.vector(dinv_link_rev(eta_pi)) * res * X0, 2, sum),
    apply(X0 * R * weights * as.vector(inv_link_rev(eta_pi)) * mu_der, 2, sum)
  ) / n
  jac
}
